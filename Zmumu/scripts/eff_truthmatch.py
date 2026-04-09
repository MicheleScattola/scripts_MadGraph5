import os
import ROOT
import math
import functions as func

# path for ATLAS root file
DATA = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/Zmumu'
TEST_FILE = os.path.join(DATA, 'Events/low_stat/tag_1_delphes_events.root')
INPUT_FILE = os.path.join(DATA,'Events/full_stat/tag_1_delphes_events.root')
TREE_NAME = 'Delphes'
OUTPUT_FOLDER = '../plots/eff_truthmatch'
# loose muons definitions
LOOSE_RECO_PT_MIN = 10.0
LOOSE_ISO_MAX = 0.20
ETA_MAX = 2.5
loose = func.restrictions(LOOSE_RECO_PT_MIN, ETA_MAX, LOOSE_ISO_MAX)
# tight muons definitions
TIGHT_RECO_PT_MIN = 25.0
TIGHT_ISO_MAX = 0.15
tight = func.restrictions(TIGHT_RECO_PT_MIN, ETA_MAX, TIGHT_ISO_MAX)

MUON_MASS = 0.105
DR_TARGET = 0.1

Z_MASS = 91.18


# MAIN
if __name__ == "__main__":

    FILE = INPUT_FILE
    print(f"Input file: {FILE}")
    print(f"Output folder: {OUTPUT_FOLDER}")
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    ROOT.gSystem.mkdir(OUTPUT_FOLDER, True)

    # load ROOT and Delphes
    ROOT.gROOT.SetBatch(True)
    #ROOT.EnableImplicitMT()
    ROOT.gSystem.Load("libDelphes.so")

    ROOT.gInterpreter.Declare(r'''
    #include "classes/SortableObject.h"
    #include "classes/DelphesClasses.h"
    #include "ExRootAnalysis/ExRootTreeReader.h"
    ''')

    # open file
    f = ROOT.TFile.Open(FILE)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open input file: {FILE}")
    tree = f.Get(TREE_NAME)
    if not tree:
        raise RuntimeError(f"Could not find tree '{TREE_NAME}'")
    
    # read tree
    reader = ROOT.ExRootTreeReader(tree)
    branch_muon = reader.UseBranch("Muon")
    branch_particle = reader.UseBranch("Particle")

    n_entries = reader.GetEntries()
    print(f"Total events = {n_entries}")


    # variables and histograms
    
    h_mu_loose_pt = ROOT.TH1F(
        "h_mu_loose_pt",
        "Loose muon p_{T}; p_{T} [GeV]; Entries",
        100, 0, 100
    )
    h_mu_tight_pt = ROOT.TH1F(
        "h_mu_tight_pt",
        "Tight muon p_{T}; p_{T} [GeV]; Entries",
        100, 0, 100
    )
    h_Z_tight_m = ROOT.TH1F(
        "h_Z_tight_m",
        "Z #rightarrow #mu #mu candidates;M [GeV]; Entries",
        100,40,140
    )
    h_eff_tot_pt = ROOT.TEfficiency(
        "h_eff_tot_pt",
        "Loose/Tight: #mu reco eff VS p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_loose_pt = ROOT.TEfficiency(
        "h_eff_loose_pt",
        "",
        100,0,100
    )
    h_eff_loose_eta = ROOT.TEfficiency(
        "h_eff_loose_eta",
        ";#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_loose_phi = ROOT.TEfficiency(
        "h_eff_loose_phi",
        ";#phi;Efficiency",
        100,-math.pi,math.pi
    )
    h_eff_tight_pt = ROOT.TEfficiency(
        "h_eff_tight_pt",
        "Loose/Tight: #mu reco eff VS p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_tight_eta = ROOT.TEfficiency(
        "h_eff_tight_eta",
        "Loose/Tight: #mu reco eff VS #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_tight_phi = ROOT.TEfficiency(
        "h_eff_tight_phi",
        "Loose/Tight: #mu reco eff VS #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )
    h_p_res_loose = ROOT.TH1F(
        "h_p_res_loose",
        "Muon P resolution;(P_{true}-P_{reco})/(P_{true});Entries",
        100,-0.1,0.1
    )
    h_p_res_tight = ROOT.TH1F(
        "h_p_res_tight",
        "Muon P resolution;(P_{true}-P_{reco})/(P_{true});Entries",
        100,-0.1,0.1
    )

    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)

        # Z from dimuons
        if tight.Zsel(branch_muon):
            #construct P4 and fill invariant mass histo
            m1 = branch_muon[0]
            m2 = branch_muon[1]
            p1 = ROOT.TLorentzVector()
            p1.SetPtEtaPhiM(m1.PT,m1.Eta,m1.Phi,MUON_MASS)
            p2 = ROOT.TLorentzVector()
            p2.SetPtEtaPhiM(m2.PT,m2.Eta,m2.Phi,MUON_MASS)
            p3 = p1+p2
            h_Z_tight_m.Fill(p3.M())


        # muon reco efficiency
        # loop MonteCarlo first
        n_particles = branch_particle.GetEntries()
        for i in range(n_particles):
            p=branch_particle.At(i)
            if p.Status != 1:
                continue
            if abs(p.PID) != 13:
                continue
            true_mu = ROOT.TLorentzVector()
            reco_mu = ROOT.TLorentzVector()
            true_mu.SetPxPyPzE(p.Px,p.Py,p.Pz,p.E)

            # initialize candidates and best variables
            isBestFound = False
            isBestLoose = False
            isBestTight = False
            reco_candidates = []
            dR_candidates = []
            isCandidateLoose = []
            isCandidateTight = []
            
            # loop over Reco muons
            for mu in branch_muon:
                #reset for each loop
                isMuLoose = False
                isMuTight = False

                reco_mu.SetPtEtaPhiM(mu.PT,mu.Eta,mu.Phi,MUON_MASS)
                delta_r = reco_mu.DeltaR(true_mu)

                # check if loose or tight
                if loose.isSel(mu):
                    isMuLoose = True
                    h_mu_loose_pt.Fill(mu.PT)

                if tight.isSel(mu):
                    isMuTight = True
                    h_mu_tight_pt.Fill(mu.PT)
                
                # check and append candidate
                if delta_r < DR_TARGET and mu.Charge * p.Charge > 0:
                    cand = ROOT.TLorentzVector()
                    cand.SetPtEtaPhiM(mu.PT, mu.Eta, mu.Phi, MUON_MASS)
                    reco_candidates.append(cand)
                    dR_candidates.append(delta_r)
                    isCandidateLoose.append(isMuLoose)
                    isCandidateTight.append(isMuTight)
                    
            
            if reco_candidates:
                # sort by dR, keep lowest candidate
                paired = sorted(
                    zip(reco_candidates, dR_candidates, isCandidateLoose, isCandidateTight),
                    key=lambda pair: pair[1]
                )
                best_reco, best_dR, isBestLoose, isBestTight = paired[0]
                isBestFound = True
                if isBestLoose:
                    h_p_res_loose.Fill((true_mu.P()-best_reco.P())/true_mu.P())
                if isBestTight:
                    h_p_res_tight.Fill((true_mu.P()-best_reco.P())/true_mu.P())

            
            h_eff_tot_pt.Fill(isBestFound,true_mu.Pt())
            h_eff_loose_pt.Fill((isBestFound and isBestLoose),true_mu.Pt())
            h_eff_loose_eta.Fill((isBestFound and isBestLoose),true_mu.Eta())
            h_eff_tight_pt.Fill((isBestFound and isBestTight),true_mu.Pt())
            h_eff_tight_eta.Fill((isBestFound and isBestTight),true_mu.Eta())
            h_eff_loose_phi.Fill((isBestFound and isBestLoose),true_mu.Eta())
            h_eff_tight_phi.Fill((isBestFound and isBestTight),true_mu.Eta())


        ### Res vs eta and phi


    ############################

    # plotting
    ROOT.gStyle.SetOptStat(0)

    ###
    c1 = ROOT.TCanvas('c1','',800,600)

    h_mu_loose_pt.SetLineColor(ROOT.kBlue)
    h_mu_loose_pt.Draw("HIST")
    h_mu_tight_pt.SetLineColor(ROOT.kRed)
    h_mu_tight_pt.Draw("SAME")

    legend = ROOT.TLegend()
    legend.AddEntry(h_mu_loose_pt, "Loose Mu", "l")
    legend.AddEntry(h_mu_tight_pt, "Tight Mu", "l")
    legend.Draw()

    c1.SaveAs(os.path.join(OUTPUT_FOLDER, "pt_mu.png"))

    ###
    c2 = ROOT.TCanvas('c2','',800,600)
    h_Z_tight_m.Draw("HIST")
    h_Z_tight_m.SetLineColor(ROOT.kRed)
    h_Z_tight_m.SetLineWidth(2)

    c2.SaveAs(os.path.join(OUTPUT_FOLDER, "Z_mass.png"))

    ###

    c3 = ROOT.TCanvas('c3','',800,600)

    h_eff_tot_pt.SetMarkerStyle(20)
    h_eff_tot_pt.SetMarkerColorAlpha(ROOT.kGreen+2,0.6)
    h_eff_tot_pt.SetLineColor(ROOT.kGreen+2)
    h_eff_tot_pt.Draw("AP")
    h_eff_tight_pt.SetMarkerStyle(20)
    h_eff_tight_pt.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_tight_pt.SetLineColor(ROOT.kRed)
    h_eff_tight_pt.Draw("SAME P")
    h_eff_loose_pt.SetMarkerStyle(20)
    h_eff_loose_pt.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_loose_pt.SetLineColor(ROOT.kBlue)
    h_eff_loose_pt.Draw("SAME P")
    

    leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.30)
    leg.AddEntry(h_eff_loose_pt,"Loose #mu","P")
    leg.AddEntry(h_eff_tight_pt,"Tight #mu","P")
    leg.AddEntry(h_eff_tot_pt,"All reco","P")
    leg.Draw()

    c3.Update()
    h_eff_tot_pt.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c3.SaveAs(os.path.join(OUTPUT_FOLDER, "efficiency_pt.png"))

    ###
    c4 = ROOT.TCanvas('c4','',800,600)
    
    h_eff_tight_eta.SetMarkerStyle(20)
    h_eff_tight_eta.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_tight_eta.SetLineColor(ROOT.kRed)
    h_eff_tight_eta.Draw("AP")
    h_eff_loose_eta.SetMarkerStyle(20)
    h_eff_loose_eta.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_loose_eta.SetLineColor(ROOT.kBlue)
    h_eff_loose_eta.Draw("SAME P")

    leg = ROOT.TLegend(0.45, 0.15, 0.65, 0.30)
    leg.AddEntry(h_eff_loose_eta,"Loose #mu","P")
    leg.AddEntry(h_eff_tight_eta,"Tight #mu","P")
    leg.Draw()

    c4.Update()
    h_eff_tight_eta.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c4.SaveAs(os.path.join(OUTPUT_FOLDER, "efficiency_eta.png"))

    ###
    c5 = ROOT.TCanvas('c5','',800,600)
    
    h_eff_tight_phi.SetMarkerStyle(20)
    h_eff_tight_phi.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_tight_phi.SetLineColor(ROOT.kRed)
    h_eff_tight_phi.Draw("AP")
    h_eff_loose_phi.SetMarkerStyle(20)
    h_eff_loose_phi.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_loose_phi.SetLineColor(ROOT.kBlue)
    h_eff_loose_phi.Draw("SAME P")

    leg = ROOT.TLegend(0.45, 0.15, 0.65, 0.30)
    leg.AddEntry(h_eff_loose_phi,"Loose #mu","P")
    leg.AddEntry(h_eff_tight_phi,"Tight #mu","P")
    leg.Draw()

    c5.Update()
    h_eff_tight_phi.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c5.SaveAs(os.path.join(OUTPUT_FOLDER, "efficiency_phi.png"))

    ###
    ROOT.gStyle.SetOptFit(1111)
    c6 = ROOT.TCanvas('c6','',800,600)

    h_p_res_loose.SetMarkerColorAlpha(ROOT.kBlue,0.9)
    h_p_res_loose.SetMarkerStyle(20)
    h_p_res_loose.Draw("HIST P")
    '''h_p_res_loose.Fit("gaus", "RQ", "",-0.1,0.1)
    fit = h_p_res.GetFunction("gaus")
    fit.SetLineColor(ROOT.kRed)
    fit.SetLineWidth(2)
    fit.Draw("SAME")'''
    h_p_res_tight.SetMarkerStyle(22)
    h_p_res_tight.SetMarkerColorAlpha(ROOT.kRed,0.9)
    h_p_res_tight.Draw("SAME P")
    #c6.Modified()
    #c6.Update()

    leg6 = ROOT.TLegend(0.65,0.70,0.88,0.88)
    leg6.AddEntry(h_p_res_loose,"Loose #mu","p")
    leg6.AddEntry(h_p_res_tight,"Tight #mu","p")
    leg6.Draw()
    c6.Update()
    
    c6.SaveAs(os.path.join(OUTPUT_FOLDER, "P_res.png"))


    # cleanup
    del leg, legend, leg6
    del c1, c2, c3, c4, c5, c6
    del h_mu_loose_pt, h_mu_tight_pt, h_Z_tight_m
    del h_eff_tot_pt, h_eff_loose_pt, h_eff_tight_pt
    del h_eff_loose_eta, h_eff_tight_eta
    del h_eff_loose_phi, h_eff_tight_phi
    del h_p_res_loose,h_p_res_tight
    del branch_muon, branch_particle, reader, tree
    del loose_muons, tight_muons, Z_candidate

    f.Close()
    del f


    #os._exit(0)