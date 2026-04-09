import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT as r


# path for ATLAS and CMS root files
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/myExercises/Zmumu'
ATLAS_events = os.path.join(folder, 'Events/ATLAS/tag_1_delphes_events.root')
CMS_events = os.path.join(folder, 'Events/CMS/tag_1_delphes_events.root')


def save_and_plot(file_paths, labels, output_dir):

    os.makedirs(output_dir, exist_ok=True)

    for file_path, label in zip(file_paths, labels):

        root_file = r.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            raise FileNotFoundError(f"Could not open file: {file_path}")

        tree = root_file.Get("Delphes")

        if not tree:
            raise FileNotFoundError(f"Tree 'Delphes' not found in: {file_path}")

        # Keep only interesting branches active
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("Muon*", 1)
 
        title = rf'#it{{p p #rightarrow Z #rightarrow #mu #mu}} at {label} ;Muon Eta;Muon Phi'
        hist = r.TH2F('hist', title, 20,0,0,20,0,0)
        hist.SetDirectory(0)

        for event in tree:
            for mu in event.Muon:
                hist.Fill(mu.Eta,mu.Phi)


        canvas = r.TCanvas('c','c',600,600)
        canvas.cd()
        hist.Draw("colz")
        #canvas.Draw()
        canvas.SaveAs(f'plots/pyroot_{label}.png')


        root_file.Close()
        del hist
        del tree
        del root_file
        del canvas




if __name__ == "__main__":

    r.gROOT.SetBatch(True)
    r.gSystem.Load("libDelphes.so")

    r.gInterpreter.Declare('#include "classes/DelphesClasses.h"')

    events = [ATLAS_events, CMS_events]

    labels = ["ATLAS", "CMS"]

    save_and_plot(events, labels, output_dir="plots")

