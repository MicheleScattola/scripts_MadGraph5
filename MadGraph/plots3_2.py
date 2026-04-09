import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import functions as f

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_3_2'
# import lhe and hempc parsers
mg5_root = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7'
sys.path.insert(0, mg5_root)

#import madgraph.various.hepmc_parser as hepmc_parser
import madgraph.various.lhe_parser as lhe_parser


basicreco_events = os.path.join(folder, 'Events/jet1/tag_1_pythia8_BasicReco.lhe.gz')


class MyData:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.events = lhe_parser.EventFile(path)
        # MA5 BasicReco output can contain empty <event> blocks.
        self.events.allow_empty_event = True
        self.pt_ee = []
        self.w_ee = []
        self.pt_j1 = []
        self.w_j1 = []
        self.pt_j2 = []
        self.w_j2 = []

    def add_ee(self, pt, weight):
        self.pt_ee.append(pt)
        self.w_ee.append(weight)

    def add_j1(self, pt, weight):
        self.pt_j1.append(pt)
        self.w_j1.append(weight)

    def add_j2(self, pt, weight):
        self.pt_j2.append(pt)
        self.w_j2.append(weight)


def select_jets(event):
    # Status 1 are final state.
    # when MadAnalysis performs the reconstructions, jets are labeled as status = 1, pid = 21
    jets = [
        p
        for p in event
        if p.status == 1 and p.pid == 21
    ]
    jets.sort(key=f.calc_p_T, reverse=True)
    return jets


def select_electrons(event):
    ep = [p for p in event if p.status == 1 and p.pid == -11]
    em = [p for p in event if p.status == 1 and p.pid == 11]
    return ep, em


def fill_dataset(dataset):
    n_ee = 0
    n_j1 = 0
    n_j2 = 0
    for event in dataset.events:
        weight = event.wgt

        ep, em = select_electrons(event)
        if ep and em:
            # Not really needed, but this picks the hardest of the e+ e- list.
            ep1 = max(ep, key=f.calc_p_T)
            em1 = max(em, key=f.calc_p_T)
            pt_ee = float(np.hypot(ep1.px + em1.px, ep1.py + em1.py))
            dataset.add_ee(pt_ee, weight)
            n_ee += 1

        jets = select_jets(event)
        if len(jets) >= 1:
            dataset.add_j1(f.calc_p_T(jets[0]), weight)
            n_j1 += 1
        if len(jets) >= 2:
            dataset.add_j2(f.calc_p_T(jets[1]), weight)
            n_j2 += 1

    print(f'[INFO] {dataset.name}: events with e+e-={n_ee}, with j1={n_j1}, with j2={n_j2}')




basicreco = MyData(name='BasicReco', path=basicreco_events)

fill_dataset(basicreco)


fig, ax = plt.subplots(3, 1, figsize=(7, 9))
fig.suptitle(r'$pp\to e^+e^- (+j)$ at $\sqrt{S}=13$ TeV')

nbins = 40
pt_range = (0, 80)

ax[0].hist(
    basicreco.pt_ee,
    bins=nbins,
    range=pt_range,
    weights=basicreco.w_ee,
    histtype='stepfilled',
    alpha=0.7,
    linewidth=1.8,
    label='BasicReco',
)
ax[0].set_ylabel('events')
ax[0].set_xlabel(r'$p_T$ [GeV]')
ax[0].set_title(r'$p_T(e^+e^-)$')
ax[0].legend()

ax[1].hist(
    basicreco.pt_j1,
    bins=nbins,
    range=pt_range,
    weights=basicreco.w_j1,
    histtype='stepfilled',
    alpha=0.7,
    linewidth=1.8,
    label='BasicReco',
)
ax[1].set_ylabel('events')
ax[1].set_xlabel(r'$p_T$ [GeV]')
ax[1].set_title(r'$p_T(j_1)$')
ax[1].legend()

ax[2].hist(
    basicreco.pt_j2,
    bins=nbins,
    range=pt_range,
    weights=basicreco.w_j2,
    histtype='stepfilled',
    alpha=0.7,
    linewidth=1.8,
    label='BasicReco',
)
ax[2].set_ylabel('events')
ax[2].set_xlabel(r'$p_T$ [GeV]')
ax[2].set_title(r'$p_T(j_2)$')
ax[2].legend()

os.makedirs('plots', exist_ok=True)
plt.tight_layout()
plt.savefig('plots/pp_ee_j.png', dpi=150, bbox_inches='tight')

