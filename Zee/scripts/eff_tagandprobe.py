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
    h_eff_sel_pt = ROOT.TEfficiency(
        "h_eff_sel_pt",
        "e^{#pm} tag & probe vs p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_sel_eta = ROOT.TEfficiency(
        "h_eff_sel_eta",
        "e^{#pm} tag & probe vs #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_sel_phi = ROOT.TEfficiency(
        "h_eff_sel_phi",
        "e^{#pm} tag & probe vs #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )

    h_eff_iso_pt = ROOT.TEfficiency(
        "h_eff_iso_pt",
        "e^{#pm} tag & probe vs p_{T};p_{T} [GeV];Efficiency",
        100,0,100
    )
    h_eff_iso_eta = ROOT.TEfficiency(
        "h_eff_iso_eta",
        "e^{#pm} tag & probe vs #eta;#eta;Efficiency",
        100,-2.5,2.5
    )
    h_eff_iso_phi = ROOT.TEfficiency(
        "h_eff_iso_phi",
        "e^{#pm} tag & probe vs #phi;#phi;Efficiency",
        100,-math.pi,math.pi
    )

    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)

        # Electron Tag & Probe
        # loop over electron, select isolated ones
        n_electrons = branch_electron.GetEntries()
        for i in range(n_electrons):
            el_tag = branch_electron.At(i)
            isTagSel = selection.isSel(el_tag)
            isTagIso = func.isIso(el_tag,branch_track,DELTA_R_MAX,ISO_MAX)

            isCandidateSel = []
            isCandidateIso = []
            candidate_mass = []
            candidate_p4 = []

            # select as tag electrons that are both:
            # kinematically selected and isolated
            if isTagSel and isTagIso:
                tag_tlv = ROOT.TLorentzVector()
                probe_tlv = ROOT.TLorentzVector()
                tag_tlv.SetPtEtaPhiM(el_tag.PT,el_tag.Eta,el_tag.Phi,func.ELECTRON_MASS)

                for j in range(n_electrons):
                    # skip itself
                    if i == j:
                        continue
                    el_probe = branch_electron.At(j)
                    probe_tlv.SetPtEtaPhiM(el_probe.PT,el_probe.Eta,el_probe.Phi,func.ELECTRON_MASS)
                    # classify probe and store
                    isProbeSel = selection.isSel(el_probe)
                    isProbeIso = func.isIso(el_probe,branch_track,DELTA_R_MAX,ISO_MAX)
                    isCandidateSel.append(isProbeSel)
                    isCandidateIso.append(isProbeIso)
                    # get candidate invariant mass
                    sum_tlv = tag_tlv + probe_tlv
                    candidate_mass.append(sum_tlv.M())
                    candidate_p4.append(sum_tlv)

                # sort for best invariant mass match
                if (candidate_mass):
                    grouped = sorted(zip(isCandidateSel,isCandidateIso,
                                        candidate_mass,candidate_p4),
                                        key=lambda pair: pair[2])
                    
                    isBestSel, isBestIso, BestMass, BestP4 = grouped[0]

                    h_eff_sel_pt.Fill(isBestSel, BestP4.Pt())
                    h_eff_sel_eta.Fill(isBestSel, BestP4.Eta())
                    h_eff_sel_phi.Fill(isBestSel, BestP4.Phi())
                    h_eff_iso_pt.Fill((isBestSel and isBestIso), BestP4.Pt())
                    h_eff_iso_eta.Fill((isBestSel and isBestIso), BestP4.Eta())
                    h_eff_iso_phi.Fill((isBestSel and isBestIso), BestP4.Phi())

                    

            
        ### Res vs eta and phi


    ############################

    # plotting
    ROOT.gStyle.SetOptStat(0)

    ###
    c1 = ROOT.TCanvas('c1','',800,600)
    h_eff_sel_pt.SetMarkerStyle(20)
    h_eff_sel_pt.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_pt.SetLineColor(ROOT.kRed)
    h_eff_sel_pt.Draw("AP")
    h_eff_iso_pt.SetMarkerStyle(20)
    h_eff_iso_pt.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_pt.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_pt.Draw("SAME P")

    leg1 = ROOT.TLegend(0.35,0.12,0.66,0.30)
    leg1.AddEntry(h_eff_sel_pt,"Selection cuts","P")
    leg1.AddEntry(h_eff_iso_pt,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg1.Draw()

    c1.Update()
    h_eff_sel_pt.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c1.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_pt.png'))

    ###
    c2 = ROOT.TCanvas('c2','',800,600)
    h_eff_sel_phi.SetMarkerStyle(20)
    h_eff_sel_phi.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_phi.SetLineColor(ROOT.kRed)
    h_eff_sel_phi.Draw("AP")
    h_eff_iso_phi.SetMarkerStyle(20)
    h_eff_iso_phi.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_phi.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_phi.Draw("SAME P")

    leg2 = ROOT.TLegend(0.35,0.12,0.66,0.30)
    leg2.AddEntry(h_eff_sel_phi,"Selection cuts","P")
    leg2.AddEntry(h_eff_iso_phi,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg2.Draw()

    c2.Update()
    h_eff_sel_phi.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c2.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_phi.png'))

    ###
    c3 = ROOT.TCanvas('c3','',800,600)
    h_eff_sel_eta.SetMarkerStyle(20)
    h_eff_sel_eta.SetMarkerColorAlpha(ROOT.kRed,0.7)
    h_eff_sel_eta.SetLineColor(ROOT.kRed)
    h_eff_sel_eta.Draw("AP")    
    h_eff_iso_eta.SetMarkerStyle(20)
    h_eff_iso_eta.SetMarkerColorAlpha(ROOT.kGreen+2,0.7)
    h_eff_iso_eta.SetLineColor(ROOT.kGreen + 2)
    h_eff_iso_eta.Draw("SAME P")

    leg3 = ROOT.TLegend(0.35,0.12,0.66,0.30)
    leg3.AddEntry(h_eff_sel_eta,"Selection cuts","P")
    leg3.AddEntry(h_eff_iso_eta,f"ISO^{{varcone}} < {ISO_MAX}","P")
    leg3.Draw()

    c3.Update()
    h_eff_sel_eta.GetPaintedGraph().GetYaxis().SetRangeUser(0.0,1.0)

    c3.SaveAs(os.path.join(OUTPUT_FOLDER, 'efficiency_eta.png'))


    # cleanup
    del leg1, leg2, leg3
    del c1, c3
    del h_eff_sel_pt, h_eff_sel_eta, h_eff_sel_phi
    del h_eff_iso_pt, h_eff_iso_eta, h_eff_iso_phi
    del branch_electron, branch_particle, branch_track, reader, tree

    f.Close()
    del f


    #os._exit(0)