import os
import ROOT
import math

import functions as func

# path for ATLAS root file
DATA = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/ttbar'
TEST_FILE = os.path.join(DATA, 'Events/low_stat/tag_1_delphes_events.root')
INPUT_FILE = os.path.join(DATA, 'Events/full_stat/tag_1_delphes_events.root')
OUTPUT_FOLDER = '../plots/eff_btagging'
TREE_NAME = 'Delphes'

# 1 lep
# 4 jets with 2 btagged
# combinations of light and b-jets to plot inv mass of anti-top

TOP_MASS = 172
# selection
PT_MIN = 25.0
ETA_MAX = 2.5
ISO_MAX = 0.15

selection = func.lepSelection(
    PT_MIN,
    ETA_MAX,
    ISO_MAX
    )

# MAIN
if __name__ == "__main__":
    
    FILE = TEST_FILE
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
    branch_muon = reader.UseBranch("Muon")
    branch_particle = reader.UseBranch("Particle")
    branch_track = reader.UseBranch("Track")
    branch_jet = reader.UseBranch("Jet")

    n_entries = reader.GetEntries()
    print(f"Total events = {n_entries}")


    # variables and histograms
    h_antitop_m = ROOT.TH1F(
        "h_antitop_m",
        "b-jet + b-jet jet jet combinations in p p#rightarrow t#bar{t};Invariant Mass [GeV];Entries",
        60,100,300
    )

    n_selected = 0
    n_jet_events = 0
    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)

        # look for isolated lepton
        lepton_branch = func.lepton_single(branch_electron,branch_muon)
        # skip event if not 1! lepton found
        if lepton_branch is None:
            continue
        # apply cuts to lepton isolation, if not selected continue
        lep = lepton_branch.At(0)
        isLepSel = selection.isSel(lep)
        if not isLepSel:
            continue
        
        n_selected += 1
        # look for number of jets and b-jets
        jet_btag = []
        jet_soft = []
        for j in branch_jet:
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM(j.PT,j.Eta,j.Phi,j.Mass)
            if(j.BTag):
                jet_btag.append(p4)
            if(not j.BTag):
                jet_soft.append(p4)
        if(len(jet_btag) != 2 or len(jet_soft) != 2):
            continue
        n_jet_events += 1

        # get invariant mass
        antitop_candidates = []
        for j in jet_btag:
            sum1 = j
            # get all possible combinations with pairs of soft jets
            for i in range(len(jet_soft)):
                sum2 = sum1 + jet_soft[i]
                for j in range(i+1,len(jet_soft)):
                    sum3 = sum2 + jet_soft[j]
                    antitop_candidates.append(sum3)


        antitop_candidates.sort(key=lambda x: abs(x.M()-TOP_MASS) )
        h_antitop_m.Fill(antitop_candidates[0].M())

    print(f'Total events with selected leptons = {n_selected}')
    print(f'Total selected events with 2+2 jets = {n_jet_events}')

    ############################

    # plotting
    #ROOT.gStyle.SetOptStat(0)

    ###
    c1 = ROOT.TCanvas('c1','',800,600)
    h_antitop_m.Draw("HIST")
    
    c1.SaveAs(os.path.join(OUTPUT_FOLDER,'antitop_m.png'))


    del branch_electron, branch_muon, branch_particle, branch_track, branch_jet, reader, tree
    del c1
    del h_antitop_m

    f.Close()
    del f


    #os._exit(0)