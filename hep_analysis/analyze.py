#!/usr/bin/env python3
import argparse
import math
import os
import sys

FILE = '/mnt/share/physics_data/4lep/MC/mc_345060.ggH125_ZZ4lep.4lep.root'

try:
    import ROOT
except Exception as e:
    print('ERROR: PyROOT is required to run this script.')
    print(e)
    sys.exit(1)

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)


def find_first_tree(root_dir):
    for key in root_dir.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom('TTree'):
            return obj
        if obj.InheritsFrom('TDirectory'):
            tree = find_first_tree(obj)
            if tree:
                return tree
    return None


def branch_names(tree):
    return {b.GetName() for b in tree.GetListOfBranches()}


def get_event_weight(tree, branches, is_data=False):
    if is_data:
        return 1.0

    w = 1.0

    if 'mcWeight' in branches:
        w *= float(tree.mcWeight)

    for name in (
        'scaleFactor_PILEUP',
        'scaleFactor_ELE',
        'scaleFactor_MUON',
        'scaleFactor_PHOTON',
        'scaleFactor_TAU',
        'scaleFactor_BTAG',
        'scaleFactor_LepTRIGGER',
        'scaleFactor_PhotonTRIGGER',
    ):
        if name in branches:
            w *= float(getattr(tree, name))

    return w


def to_list(x):
    if x is None:
        return []
    try:
        return list(x)
    except TypeError:
        return [x]


def get_collection(tree, name, branches):
    if name not in branches:
        return []
    try:
        return to_list(getattr(tree, name))
    except Exception:
        return []


def safe_at(seq, idx, default=0.0):
    return seq[idx] if idx < len(seq) else default


def build_lepton_p4(pt, eta, phi, mass=0.0):
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM(pt, eta, phi, mass)
    return v


def make_hists():
    h = {}
    h['cutflow'] = ROOT.TH1F('cutflow', 'Cutflow;Selection;Events', 9, 0.5, 9.5)
    labels = [
        'all', 'lep_sel', '==4 lep', 'SFOC', 'p_T', 'isolation', 'm12', 'm34','m1234'
    ]
    for i, lab in enumerate(labels, start=1):
        h['cutflow'].GetXaxis().SetBinLabel(i, lab)

    h['nlep'] = ROOT.TH1F('nlep', 'Lepton multiplicity;N_{lep};Events', 8, -0.5, 7.5)
    h['lead_lep_pt'] = ROOT.TH1F('lead_lep_pt', 'Leading lepton p_{T};p_{T}^{lead} [GeV];Events', 50, 0.0, 300.0)
    h['sublead_lep_pt'] = ROOT.TH1F('sublead_lep_pt', 'Subleading lepton p_{T};p_{T}^{sublead} [GeV];Events', 50, 0.0, 200.0)
    h['lead_mass'] = ROOT.TH1F('lead_mass','Leading ll pair invariant mass;M [GeV];Events',50,20.,120.)
    h['sublead_mass'] = ROOT.TH1F('sublead_mass','Sub-leading ll pair invariant mass;M [GeV];Events',50,20.,120.)
    h['candidate_mass'] = ROOT.TH1F('candidate_mass','4 leptons candidate mass;M [GeV];Events',50,100.,170.)

    return h


def save_hists(hists, outdir):
    os.makedirs(outdir, exist_ok=True)
    fout = ROOT.TFile(os.path.join(outdir, 'analysis_histograms.root'), 'RECREATE')
    for hist in hists.values():
        hist.Write()
    fout.Close()

    canvas = ROOT.TCanvas('c', 'c', 800, 600)
    for name, hist in hists.items():
        canvas.Clear()
        if name in ('cutflow',):
            hist.LabelsOption('v', 'X')
        hist.Draw('hist e')
        canvas.SaveAs(os.path.join(outdir, f'{name}.png'))


def main():
    parser = argparse.ArgumentParser(description='Generic PyROOT analysis for an exactly-2-lepton sample.')
    parser.add_argument('-input', default=FILE, help='Input ROOT file')
    parser.add_argument('-o', '--outdir', default='plots', help='Output directory')
    parser.add_argument('--met-cut', type=float, default=20.0, help='MET cut in GeV for cutflow')
    parser.add_argument('--min-jets', type=int, default=2, help='Minimum jet multiplicity for cutflow')
    args = parser.parse_args()

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f'Could not open {args.input}')

    tree = find_first_tree(f)
    if not tree:
        raise RuntimeError('No TTree found in input file')

    is_data = 'data' in args.input.lower()
    print(is_data)

    print(f'Using tree: {tree.GetName()}')
    branches = branch_names(tree)
    print('Found branches:')
    for name in sorted(branches):
        print(f'  {name}')

    required = ['lep_pt', 'lep_eta', 'lep_phi']
    missing = [x for x in required if x not in branches]
    if missing:
        print('WARNING: missing expected branches:', ', '.join(missing))

    h = make_hists()

    nentries = tree.GetEntries()
    print(f'Processing {nentries} entries')

    for i, event in enumerate(tree):
        if i and i % 100000 == 0:
            print(f'  processed {i} / {nentries}')

        w = get_event_weight(tree, branches, is_data=is_data)
        h['cutflow'].Fill(1, w)

        lep_pt = get_collection(tree, 'lep_pt', branches)
        lep_eta = get_collection(tree, 'lep_eta', branches)
        lep_phi = get_collection(tree, 'lep_phi', branches)
        lep_charge = get_collection(tree, 'lep_charge', branches)
        lep_type = get_collection(tree, 'lep_type', branches)

        nlep = len(lep_pt)
        h['nlep'].Fill(nlep, w)

        h['cutflow'].Fill(2, w)

        order = sorted(range(nlep), key=lambda j: lep_pt[j], reverse=True)
        i1, i2, i3, i4 = order[0], order[1], order[2], order[3]
        pt1, pt2, pt3, pt4 = lep_pt[i1], lep_pt[i2], lep_pt[i3], lep_pt[i4]
        eta1, eta2, eta3, eta4 = safe_at(lep_eta, i1), safe_at(lep_eta, i2), safe_at(lep_eta, i3), safe_at(lep_eta, i4)
        phi1, phi2, phi3, phi4 = safe_at(lep_phi, i1), safe_at(lep_phi, i2), safe_at(lep_phi, i3), safe_at(lep_phi, i4)

        # lepton selection
        for i in range(nlep):
            if safe_at(lep_type,i) == 11:
                isSel = checkEl()

        h['lead_lep_pt'].Fill(pt1/1000., w)
        h['sublead_lep_pt'].Fill(pt2/1000., w)

        q1 = safe_at(lep_charge, i1, 0.0)
        q2 = safe_at(lep_charge, i2, 0.0)
        is_os = (q1 * q2 < 0.0)
        
        if not is_os:
            continue

        h['cutflow'].Fill(3, w)

        l1 = build_lepton_p4(pt1, eta1, phi1)
        l2 = build_lepton_p4(pt2, eta2, phi2)
        l3 = build_lepton_p4(pt3, eta3, phi3)
        l4 = build_lepton_p4(pt4, eta4, phi4)

    save_hists(h, args.outdir)
    print(f'Wrote output to: {args.outdir}')


if __name__ == '__main__':
    main()
