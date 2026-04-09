import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import functions as func

# path for ATLAS root file
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/MuonCP'
INPUT_FILE = os.path.join(folder,'data/Zmumu_ATLAS_delphes_events.root')
TREE_NAME = 'Delphes'
OUTPUT_PLOT = os.path.join(folder, 'plots', 'muon_pt_selection.png')

# loose muons definitions
LOOSE_RECO_PT_MIN = 10.0
LOOSE_ISO_MAX = 0.25
ETA_MAX = 2.5
loose = func.restrictions(LOOSE_RECO_PT_MIN, ETA_MAX, LOOSE_ISO_MAX)
# tight muons definitions
TIGHT_RECO_PT_MIN = 25.0
TIGHT_ISO_MAX = 0.15
tight = func.restrictions(TIGHT_RECO_PT_MIN, ETA_MAX, TIGHT_ISO_MAX)


if __name__ == "__main__":

    # load ROOT and Delphes
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libDelphes.so")

    ROOT.gInterpreter.Declare(r'''
    #include "classes/SortableObject.h"
    #include "classes/DelphesClasses.h"
    #include "ExRootAnalysis/ExRootTreeReader.h"
    ''')

    # open file
    f = ROOT.TFile.Open(INPUT_FILE)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open input file: {INPUT_FILE}")
    tree = f.Get(TREE_NAME)
    if not tree:
        raise RuntimeError(f"Could not find tree '{TREE_NAME}'")
    
    # read tree
    reader = ROOT.ExRootTreeReader(tree)
    branch_muon = reader.UseBranch("Muon")
    branch_particle = reader.UseBranch("Particle")

    n_entries = reader.GetEntries()
    print(f"Total events = {n_entries}")


    # variables: save only PT values to regular Python lists
    loose_muons = []
    tight_muons = []

    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)
    
        for i in range(branch_muon.GetEntries()):
            mu = branch_muon.At(i)
            #select loose and tight
            if loose.isSel(mu):
                loose_muons.append(mu.PT)
            if tight.isSel(mu):
                tight_muons.append(mu.PT)

    # Convert to numpy arrays to isolate from ROOT memory
    loose_muons = np.array(loose_muons)
    tight_muons = np.array(tight_muons)

    # plot with matplotlib
    print(f"Loose muons: {len(loose_muons)}")
    print(f"Tight muons: {len(tight_muons)}")

    os.makedirs(os.path.dirname(OUTPUT_PLOT), exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.hist(loose_muons, bins=50, label=r'loose $\mu$', histtype='step', linewidth=2, color='blue')
    ax.hist(tight_muons, bins=50, label=r'tight $\mu$', histtype='step', linewidth=2, color='red')
    ax.set_xlabel(r'Muon $p_T$ [GeV]')
    ax.set_ylabel('Counts')
    ax.set_title(r'Muon $p_T$ Selection')
    ax.legend()
    fig.savefig(OUTPUT_PLOT, dpi=150)
    print(f"Saved plot: {OUTPUT_PLOT}")
    
    # Force exit without cleanup to avoid PyROOT segfaults
    os._exit(0)






    
    

