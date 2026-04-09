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
ISO_MAX = 0.10

selection = func.restrictions(
    PT_MIN,
    ETA_MAX,TRANSITION_LOW,TRANSITION_HIGH,
    D0_RES_MAX,
    Z0SINTHETA_MAX
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
    branch_particle = reader.UseBranch("Particle")
    branch_track = reader.UseBranch("Track")

    n_entries = reader.GetEntries()
    n_electrons = branch_electron.GetEntries()


    # variables and histograms
    bad_events = 0


    # loop over data
    for event in range(n_entries):
        reader.ReadEntry(event)

        z_id, el1_id, el2_id = func.getDaughters(branch_particle,62,23)
        # func.check_z(branch_particle)
        # check size of daughters
        #print(f"Z (id={z_id}) --> e (id={el1_id}) e (id={el2_id})")
        # skip if z was not found properly
        if z_id<0 or el1_id<0 or el2_id<0:
            bad_events += 1
            continue


    print(f"Total events = {n_entries} , Total bad events = {bad_events}")

        
    ############################

    # plotting
    #ROOT.gStyle.SetOptStat(0)

    
    f.Close()
    del f