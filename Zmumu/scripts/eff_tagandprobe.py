import os
import ROOT
import math
import functions as func

# path for ATLAS root file
DATA = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/Zmumu'
TEST_FILE = os.path.join(DATA, 'Events/low_stat/tag_1_delphes_events.root')
INPUT_FILE = os.path.join(DATA,'Events/full_stat/tag_1_delphes_events.root')
TREE_NAME = 'Delphes'
OUTPUT_FOLDER = '../plots/eff_tagandprobe'
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
        "#mu tag & probe vs p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_tight_eta = ROOT.TEfficiency(
        "h_eff_tight_eta",
        "#mu tag & probe vs #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_tight_phi = ROOT.TEfficiency(
        "h_eff_tight_phi",
        "#mu tag & probe vs #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )
    h_m_tight = ROOT.TH1F(
        "h_m_tight",
        "Z candidates (tight+tight);M[GeV];Entries",
        100,40,140
    )
    h_m_loose = ROOT.TH1F(
        "h_m_loose",
        "Z candidates (tight+tight);M[GeV];Entries",
        100,40,140
    )

    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)

        #Muon Tag & Probe
        # loop over reco muons, if tight loop
        N_MUONS = branch_muon.GetEntries()
        for i in range(N_MUONS):
            mu_tag = branch_muon.At(i)

            if tight.isSel(mu_tag):
                tag_tlv = ROOT.TLorentzVector()
                probe_tlv = ROOT.TLorentzVector()
                tag_tlv.SetPtEtaPhiM(mu_tag.PT,mu_tag.Eta,mu_tag.Phi,MUON_MASS)
                isMuLoose = []
                isMuTight = []
                candidate_mass = []
                candidate_p4 = []

                for j in range(N_MUONS):
                    # catch other probes and compute invariant mass
                    if j==i:
                        continue
                    mu_probe = branch_muon.At(j)
                    probe_tlv.SetPtEtaPhiM(mu_probe.PT,mu_probe.Eta,mu_probe.Phi,MUON_MASS)
                    sum_tlv = tag_tlv + probe_tlv
                    
                    # is the probe tight or loose?
                    isMuTight.append(tight.isSel(mu_probe))
                    isMuLoose.append(loose.isSel(mu_probe))
                    candidate_mass.append(sum_tlv.M())
                    candidate_p4.append(sum_tlv)

                # now sort by invariant mass closest to Z mass
                grouped = sorted( zip(candidate_mass,isMuLoose,isMuTight, candidate_p4),
                                 key=lambda pair: abs(Z_MASS - pair[0]) )
                if not grouped:
                    continue
                best_mass, isBestLoose, isBestTight, best_p4 = grouped[0]
                if(isBestLoose):
                    h_m_loose.Fill(best_mass)
                if(isBestTight):
                    h_m_tight.Fill(best_mass)

                h_eff_loose_pt.Fill(isBestLoose,best_p4.Pt())
                h_eff_tight_pt.Fill(isBestTight,best_p4.Pt())
                h_eff_loose_eta.Fill(isBestLoose,best_p4.Eta())
                h_eff_tight_eta.Fill(isBestTight,best_p4.Eta())
                h_eff_loose_phi.Fill(isBestLoose,best_p4.Phi())
                h_eff_tight_phi.Fill(isBestTight,best_p4.Phi())



        ### Res vs eta and phi


    ############################

    # plotting
    ROOT.gStyle.SetOptStat(0)

    ###

    c1 = ROOT.TCanvas('c1','',800,600)

    h_eff_tight_pt.SetMarkerStyle(20)
    h_eff_tight_pt.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_tight_pt.SetLineColor(ROOT.kRed)
    h_eff_tight_pt.Draw("AP")
    h_eff_loose_pt.SetMarkerStyle(20)
    h_eff_loose_pt.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_loose_pt.SetLineColor(ROOT.kBlue)
    h_eff_loose_pt.Draw("SAME P")
    

    leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.30)
    leg.AddEntry(h_eff_loose_pt,"Loose #mu","P")
    leg.AddEntry(h_eff_tight_pt,"Tight #mu","P")
    leg.Draw()

    c1.Update()
    h_eff_tight_pt.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c1.SaveAs(os.path.join(OUTPUT_FOLDER,'efficiency_pt.png'))

    ###
    c2 = ROOT.TCanvas('c2','',800,600)
    
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

    c2.Update()
    h_eff_tight_eta.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c2.SaveAs(os.path.join(OUTPUT_FOLDER,'efficiency_eta.png'))

    ###
    c3 = ROOT.TCanvas('c3','',800,600)
    
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

    c3.Update()
    h_eff_tight_phi.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c3.SaveAs(os.path.join(OUTPUT_FOLDER,'efficiency_phi.png'))

    ###
    c4 = ROOT.TCanvas('c4','',800,600)
    h_m_loose.Draw("HIST")
    h_m_loose.SetLineColor(ROOT.kBlue)
    h_m_loose.SetLineWidth(2)

    c4.SaveAs(os.path.join(OUTPUT_FOLDER,'Z_mass_loose.png'))

    ###
    c5 = ROOT.TCanvas('c5','',800,600)
    h_m_tight.Draw("HIST")
    h_m_tight.SetLineColor(ROOT.kRed)
    h_m_tight.SetLineWidth(2)

    c5.SaveAs(os.path.join(OUTPUT_FOLDER,'Z_mass_tight.png'))

    # cleanup
    del leg
    del c1, c2, c3, c4, c5
    del h_eff_loose_pt, h_eff_tight_pt
    del h_eff_loose_eta, h_eff_tight_eta
    del h_eff_loose_phi, h_eff_tight_phi
    del h_m_loose, h_m_tight
    del branch_muon, branch_particle, reader, tree

    f.Close()
    del f


    #os._exit(0)