import argparse
import os
import sys

try:
    import ROOT
except ImportError as e:
    print('ERROR: PyROOT is required to run this script.')
    sys.exit(1)

ROOT.ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
_this_dir = os.path.dirname(os.path.abspath(__file__))
_functions_h = os.path.join(_this_dir, 'functions.h')
ROOT.gInterpreter.Declare(f'#include "{_functions_h}"')


def find_first_tree(root_file_path):
    """Open the file temporarily to find the first TTree name."""
    f = ROOT.TFile.Open(root_file_path)
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom('TTree'):
            name = obj.GetName()
            f.Close()
            return name
        if obj.InheritsFrom('TDirectory'):
            # Simplified for this example, assuming tree is top-level
            pass
    f.Close()
    return None


def get_available_branches(tree_name, file_path):
    """Retrieve all branch names to dynamically build the weight string."""
    f = ROOT.TFile.Open(file_path)
    tree = f.Get(tree_name)
    branches = [b.GetName() for b in tree.GetListOfBranches()]
    f.Close()
    return branches

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

def save_hists(hists, outdir):
    os.makedirs(outdir, exist_ok=True)
    fout = ROOT.TFile(os.path.join(outdir, 'analysis.root'), 'RECREATE')
    for hist in hists.values():
        hist.Write()
    fout.Close()

    canvas = ROOT.TCanvas('c', 'c', 800, 600)
    for name, hist in hists.items():
        canvas.Clear()
        if name == 'cutflow':
            hist.LabelsOption('v', 'X')
        elif hist.InheritsFrom('TH2'):
            hist.Draw('colz')
        else:
            hist.Draw('hist e')
        canvas.SaveAs(os.path.join(outdir, f'{name}.png'))

def report_to_string(report):
    result = ""
    for entry in report:
        name = entry.GetName()
        passed = entry.GetPass()
        total = entry.GetAll()
        eff = entry.GetEff()
        cumulative = 100.0 * passed / report.begin().__deref__().GetAll() if report.begin() != report.end() else 0
        result += f"{name:10}: pass={passed:<10} all={total:<10} -- eff={eff:.2f} % cumulative eff={cumulative:.2f} %\n"
    return result


def main():
    FILE = '/mnt/share/physics_data/4lep/MC/mc_345060.ggH125_ZZ4lep.4lep.root'
    parser = argparse.ArgumentParser(description='RDataFrame exactly-2-lepton sample.')
    parser.add_argument('-i', '--input', default=FILE, help='Input ROOT file')
    parser.add_argument('-o', '--outdir', default='plots_rdf', help='Output directory')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise RuntimeError(f'Could not find {args.input}')

    tree_name = find_first_tree(args.input)
    if not tree_name:
        raise RuntimeError('No TTree found in input file')

    print(f'Using tree: {tree_name}')
    is_data = 'data' in args.input.lower()
    
    df = ROOT.RDataFrame(tree_name, args.input)
    branches = get_available_branches(tree_name, args.input)

    # WEIGHTS
    if is_data:
        w = "1.0"
    else:
        sf_branches = [
            'mcWeight', 'scaleFactor_PILEUP', 'scaleFactor_ELE', 'scaleFactor_MUON',
            'scaleFactor_PHOTON', 'scaleFactor_TAU', 'scaleFactor_BTAG', 
            'scaleFactor_LepTRIGGER', 'scaleFactor_PhotonTRIGGER'
        ]
        available_sfs = [b for b in sf_branches if b in branches]
        if available_sfs:
            sf = [b for b in available_sfs]
            w = " * ".join(sf)
        else:
            w = "1.0"

    df = df.Define("weight", w)
    df = df.Define("nlep", "lep_pt.size()")

    # HISTOS
    hists = {}

    df = df.Filter("trigE == 1 || trigM == 1","triggered")


    # lepton basic selection: p_T, eta, d0/sigma_d0
    df = df.Define("all_leptons", "ROOT::VecOps::Construct<Particle>(lep_pt, lep_eta, lep_phi,lep_trackd0pvunbiased, lep_tracksigd0pvunbiased,lep_charge, lep_type)") \
            .Define("selected_leptons","selectGoodLeptons(all_leptons)")

    # four leptons
    df = df.Filter("selected_leptons.size() == 4", "Exactly 4 leptons")

    # charge and flavor selection
    df = df.Filter("isGoodQuadruplet(selected_leptons)","only SFOC quadruplets")

    # final reconstruction
    df = df.Define("H_candidates", "leptonPairing(selected_leptons)")

    df = df.Define("m_lead","diLeptonMass(H_candidates,true)")
    df = df.Define("m_sublead","diLeptonMass(H_candidates,false)")
    df = df.Define("m_4l", "fourLeptonMass(H_candidates)")
    df = df.Define("DeltaR_lead","diLeptonDeltaR(H_candidates,true)")
    df = df.Define("DeltaR_sublead","diLeptonDeltaR(H_candidates,false)")
    df = df.Define("tot_weight","weight * 1000 * 139 * XSection / SumWeights")
    
    df = df.Define("eta_lead","diLeptonEta(H_candidates,true)")
    df = df.Define("eta_sublead","diLeptonEta(H_candidates,false)")
    df = df.Define("phi_lead","diLeptonPhi(H_candidates,true)")
    df = df.Define("phi_sublead","diLeptonPhi(H_candidates,false)")

    df = df.Filter("H_candidates.size() == 2", "lepton pairing and pt")
    df = df.Filter("m_lead > 50. && m_lead < 106.","leading mass window")
    df = df.Filter("m_sublead > 12. && m_sublead < 115.","subleading mass window")
    df = df.Filter("m_4l > 5. && m_4l < 160.","4 leptons mass window")
    df = df.Filter("DeltaR_lead > 0.1 && DeltaR_sublead > 0.1","lepton pair isolation")


    hists['m_4l'] = df.Histo1D(
        ('m_4l', '4-Lepton Invariant Mass;m_{4l} [GeV];Events', 50, 80.0, 160.0), 
        'm_4l', 
        'weight'
    )
    hists['m_lead'] = df.Histo1D(
        ('m_lead','leading lepton pair;m_{ll} [GeV];Events',50,40.,110.),
        'm_lead',
        'weight'
    )
    hists['m_sublead'] = df.Histo1D(
        ('m_sublead','subleading lepton pair;m_{ll} [GeV];Events',50,5.,120.),
        'm_sublead',
        'weight'
    )
    hists['eta_phi_lead'] = df.Histo2D(
        ('eta_phi_lead', 'Leading dileptons;#eta;#phi', 50, -2.8, 2.8, 50, -3.2, 3.2),
        'eta_lead',
        'phi_lead',
        'weight'
    )
    hists['eta_phi_sublead'] = df.Histo2D(
        ('eta_phi_sublead','Sub-leading dileptons;#eta;#phi',50,-2.8, 2.8, 50, -3.2, 3.2),
        'eta_sublead',
        'phi_sublead',
        'weight'
    )

    print("Executing event loop...")
    # Resolve histogram pointers

    print("\n--- RDataFrame Filter Report ---")
    df.Report().Print()

    weights_total = df.Sum("tot_weight").GetValue()
    print(f'\n--- TOTAL WEIGHT = {weights_total}\n')

    save_hists(hists, args.outdir)
    with open(os.path.join(args.outdir,'report.txt'),'w') as f:
        f.write(report_to_string(df.Report()))
        f.write(f'\nTOTAL WEIGHT = {weights_total}\n')
    print(f'\nWrote output to: {args.outdir}')

if __name__ == '__main__':
    main()