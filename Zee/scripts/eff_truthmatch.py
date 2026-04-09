import os
import ROOT
import math

import functions as func

# path for ATLAS root file
DATA = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/Zee'
TEST_FILE = os.path.join(DATA, 'Events/low_stat/tag_1_delphes_events.root')
INPUT_FILE = os.path.join(DATA, 'Events/full_stat/tag_1_delphes_events.root')
OUTPUT_FOLDER = '../plots/eff_truthmatch'
TREE_NAME = 'Delphes'

# electrons selection
PT_MIN = 15.0
ETA_MAX = 2.47
TRANSITION_LOW = 1.37
TRANSITION_HIGH = 1.52
D0_RES_MAX = 5.0
Z0SINTHETA_MAX = 3 # to be defined
DELTA_R_MAX = 0.2
ISO_MAX = 0.05

selection = func.restrictions(
    PT_MIN,
    ETA_MAX,TRANSITION_LOW,TRANSITION_HIGH,
    D0_RES_MAX,
    Z0SINTHETA_MAX
    )

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
    branch_electron = reader.UseBranch("Electron")
    branch_particle = reader.UseBranch("Particle")
    branch_track = reader.UseBranch("Track")

    n_entries = reader.GetEntries()
    print(f"Total events = {n_entries}")


    # variables and histograms
    h_eff_reco_pt = ROOT.TEfficiency(
        "h_eff_reco_pt",
        "Electrons reco eff. VS p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_reco_eta = ROOT.TEfficiency(
        "h_eff_reco_eta",
        "Electrons reco eff. VS #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_reco_phi = ROOT.TEfficiency(
        "h_eff_reco_phi",
        "Electrons reco eff. VS #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )

    h_eff_sel_pt = ROOT.TEfficiency(
        "h_eff_sel_pt",
        "Selection eff. p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_sel_eta = ROOT.TEfficiency(
        "h_eff_sel_eta",
        "Selection eff. #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_sel_phi = ROOT.TEfficiency(
        "h_eff_sel_phi",
        "Isolation eff. #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )

    h_eff_iso_pt = ROOT.TEfficiency(
        "h_eff_iso_pt",
        "Isolation eff. p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_iso_eta = ROOT.TEfficiency(
        "h_eff_iso_eta",
        "Isolation eff. #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_iso_phi = ROOT.TEfficiency(
        "h_eff_iso_phi",
        "Isolation eff. #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )

    h_dR = ROOT.TH1F(
        "h_dR",
        "#DeltaR in best reco candidates;#Delta R;Entries",
        100,0.0,0.2
    )
    h_isoVar = ROOT.TH1F(
        "h_isoVar",
        "Isolation variable in best reco candidates;Iso Var;Entries",
        100,0.0,0.5
    )

    h_p_res_sel = ROOT.TH1F(
        "h_p_res_sel",
        "Electron P resolution;(P_{true}-P_{reco})/(P_{true});Entries",
        100,-0.1,0.1
    )
    h_p_res_iso = ROOT.TH1F(
        "h_p_res_iso",
        "Electron P resolution;(P_{true}-P_{reco})/(P_{true});Entries",
        100,-0.1,0.1
    )



    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)
        #tot_true_el = 0
        #tot_true_el_Z_mother = 0

        # collect true electrons for truth matching
        for p in branch_particle:
            if p.Status != 1:
                continue
            if abs(p.PID) != 11:
                continue

            true_el = ROOT.TLorentzVector()
            true_el.SetPxPyPzE(p.Px,p.Py,p.Pz,p.E)

            isBestFound = False
            isBestSelected = False
            isBestIsolated = False
            isCandidateSelected = []
            isCandidateIsolated = []
            dR_candidates = []
            isoVar_candidates = []
            p4_candidates = []

            # find reconstructed candidate
            for el in branch_electron:
                # does it pass selection?
                isElSel = selection.isSel(el)
                # is it isolated?
                isElIso = func.isIso(el,branch_track,DELTA_R_MAX,ISO_MAX)

                # is the electron matched to a monte carlo particle?
                delta_r = func.delta_r(el,p)
                if delta_r < DELTA_R_MAX and p.Charge * el.Charge > 0:
                    p4 = ROOT.TLorentzVector()
                    p4.SetPtEtaPhiM(el.PT,el.Eta,el.Phi,func.ELECTRON_MASS)
                    p4_candidates.append(p4)
                    isCandidateSelected.append(isElSel)
                    isCandidateIsolated.append(isElIso)
                    isoVar_candidates.append(el.IsolationVar)
                    dR_candidates.append(delta_r)


                
            if (dR_candidates):
                isBestFound = True
                grouped = sorted(zip(isCandidateSelected,isCandidateIsolated,dR_candidates,p4_candidates),key=lambda pair: pair[2])
                isBestSelected, isBestIsolated, best_dR, reco_el = grouped[0]

                if (isBestSelected):
                    h_p_res_sel.Fill((true_el.P()-reco_el.P())/true_el.P())
                if (isBestIsolated):
                    h_p_res_iso.Fill((true_el.P()-reco_el.P())/true_el.P())

            h_eff_reco_pt.Fill(isBestFound,true_el.Pt())
            h_eff_reco_eta.Fill(isBestFound,true_el.Eta())
            h_eff_reco_phi.Fill(isBestFound,true_el.Phi())

            h_eff_sel_pt.Fill((isBestFound and isBestSelected),true_el.Pt())
            h_eff_sel_eta.Fill((isBestFound and isBestSelected),true_el.Eta())
            h_eff_sel_phi.Fill((isBestFound and isBestSelected),true_el.Phi())

            h_eff_iso_pt.Fill((isBestFound and isBestIsolated and isBestSelected),true_el.Pt())
            h_eff_iso_eta.Fill((isBestFound and isBestIsolated and isBestSelected),true_el.Eta())
            h_eff_iso_phi.Fill((isBestFound and isBestIsolated and isBestSelected),true_el.Phi())
      


    ############################

    # plotting
    ROOT.gStyle.SetOptStat(0)

    ###
    c1 = ROOT.TCanvas('c1','',800,600)
    h_eff_reco_pt.SetMarkerStyle(20)
    h_eff_reco_pt.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_reco_pt.SetLineColor(ROOT.kBlue)
    h_eff_reco_pt.Draw("AP")
    h_eff_sel_pt.SetMarkerStyle(20)
    h_eff_sel_pt.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_pt.SetLineColor(ROOT.kRed)
    h_eff_sel_pt.Draw("SAME P")
    h_eff_iso_pt.SetMarkerStyle(20)
    h_eff_iso_pt.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_pt.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_pt.Draw("SAME P")

    leg1 = ROOT.TLegend(0.35,0.12,0.70,0.30)
    leg1.AddEntry(h_eff_reco_pt,f"Reco #DeltaR<{DELTA_R_MAX}","P")
    leg1.AddEntry(h_eff_sel_pt,"Selection cuts","P")
    leg1.AddEntry(h_eff_iso_pt,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg1.Draw()

    c1.Update()
    h_eff_reco_pt.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c1.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_pt.png'))

    ###
    c2 = ROOT.TCanvas('c2','',800,600)
    h_eff_reco_phi.SetMarkerStyle(20)
    h_eff_reco_phi.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_reco_phi.SetLineColor(ROOT.kBlue)
    h_eff_reco_phi.Draw("AP")
    h_eff_sel_phi.SetMarkerStyle(20)
    h_eff_sel_phi.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_phi.SetLineColor(ROOT.kRed)
    h_eff_sel_phi.Draw("SAME P")
    h_eff_iso_phi.SetMarkerStyle(20)
    h_eff_iso_phi.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_phi.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_phi.Draw("SAME P")

    leg2 = ROOT.TLegend(0.35,0.70,0.70,0.88)
    leg2.AddEntry(h_eff_reco_phi,f"Reco #DeltaR<{DELTA_R_MAX}","P")
    leg2.AddEntry(h_eff_sel_phi,"Selection cuts","P")
    leg2.AddEntry(h_eff_iso_phi,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg2.Draw()

    c2.Update()
    h_eff_reco_phi.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c2.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_phi.png'))

    ###
    c3 = ROOT.TCanvas('c3','',800,600)
    h_eff_reco_eta.SetMarkerStyle(20)
    h_eff_reco_eta.SetMarkerColorAlpha(ROOT.kBlue,0.7)
    h_eff_reco_eta.SetLineColor(ROOT.kBlue)
    h_eff_reco_eta.Draw("AP")
    h_eff_sel_eta.SetMarkerStyle(20)
    h_eff_sel_eta.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_eta.SetLineColor(ROOT.kRed)
    h_eff_sel_eta.Draw("SAME P")    
    h_eff_iso_eta.SetMarkerStyle(20)
    h_eff_iso_eta.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_eta.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_eta.Draw("SAME P")

    leg3 = ROOT.TLegend(0.35,0.12,0.66,0.30)
    leg3.AddEntry(h_eff_reco_eta,f"Reco #DeltaR<{DELTA_R_MAX}","P")
    leg3.AddEntry(h_eff_sel_eta,"Selection cuts","P")
    leg3.AddEntry(h_eff_iso_eta,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg3.Draw()

    c3.Update()
    h_eff_reco_eta.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c3.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_eta.png'))

    ###
    #ROOT.gStyle.SetOptFit(1111)
    c4 = ROOT.TCanvas('c4','',800,600)

    h_p_res_sel.SetMarkerColorAlpha(ROOT.kBlue,0.9)
    h_p_res_sel.SetMarkerStyle(20)
    h_p_res_sel.Draw("HIST P")
    h_p_res_iso.SetMarkerStyle(22)
    h_p_res_iso.SetMarkerColorAlpha(ROOT.kRed,0.9)
    h_p_res_iso.Draw("SAME P")

    leg4 = ROOT.TLegend(0.65,0.70,0.88,0.88)
    leg4.AddEntry(h_p_res_sel,"Selection cuts","p")
    leg4.AddEntry(h_p_res_iso,f"ISO^{{varcone}} < {ISO_MAX}","p")
    leg4.Draw()
    c4.Update()
    
    c4.SaveAs(os.path.join(OUTPUT_FOLDER, "P_res.png"))


    # cleanup
    del c1, c2, c3, c4
    del h_eff_reco_pt, h_eff_reco_eta, h_eff_reco_phi
    del h_eff_sel_pt, h_eff_sel_phi, h_eff_sel_eta
    del h_eff_iso_pt, h_eff_iso_phi, h_eff_iso_eta
    del h_p_res_sel, h_p_res_iso
    del leg1, leg2, leg3, leg4
    del branch_electron, branch_particle, branch_track, reader, tree

    f.Close()
    del f


    #os._exit(0)