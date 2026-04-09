import sys
import os
import matplotlib.pyplot as plt
import functions as f

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_3_1'
bin_internal = os.path.join(folder, 'bin')
mg5_root = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7'
sys.path.insert(0, mg5_root)
sys.path.insert(0, bin_internal)

import internal.lhe_parser as lhe_parser
import madgraph.various.hepmc_parser as hepmc_parser

# load events into parser
parton_events = os.path.join(folder,'Events/shower/unweighted_events.lhe.gz')
shower_events = os.path.join(folder, 'Events/shower/tag_1_pythia8_events.hepmc.gz')
print('[INFO] loaded events into parser.')

# need different parsers for lhe and hepmc data
def open_event_file(path):
    if path.endswith('.hepmc') or path.endswith('.hepmc.gz'):
        return hepmc_parser.HEPMC_EventFile(path)
    return lhe_parser.EventFile(path)

# class for data
class MyData:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.events = open_event_file(path)
        self.delta_R = []
        self.pT_1 = []
        self.pT_2 = []
        self.weights = []

    def load(self, delta, pt1, pt2, wgt):
        self.delta_R.append(delta)
        self.pT_1.append(pt1)
        self.pT_2.append(pt2)
        self.weights.append(wgt)

# define datasets
parton = MyData(name='parton', path=parton_events)
print('[DEBUG] Parton dataset created')
shower = MyData(name='shower', path=shower_events)
print('[DEBUG] Shower dataset created')
datasets = [parton, shower]

print('[INFO] defined datasets.')

# fill datasets
for data in datasets:
    sample = data.events

    # cycle through events
    for e in sample:

        # find the two e+ e-
        ep = f.find_particle(e,11)
        em = f.find_particle(e,-11)

        if not ep or not em:
            raise RuntimeError('[ERROR] No electron or positron found in decay!!')
        
        # calculate variables
        delta_r = f.calc_delta_R(ep, em)
        pt_p = f.calc_p_T(ep)
        pt_m = f.calc_p_T(em)
        #fill class
        data.load(delta_r, pt_p, pt_m, e.wgt)

print('[INFO] filled datasets.')

# plotting
fig, ax = plt.subplots(3, figsize=(7,10))

fig.suptitle(r'$pp \rightarrow e^+ e^-$ events at $\sqrt{S}=13$ TeV')
nbins=20
range=(0,80)

ax[0].hist(parton.delta_R, weights=parton.weights, bins=nbins, range=(0,5), label=r'$\Delta R$ (parton)', alpha=0.6)
ax[0].hist(shower.delta_R, weights=shower.weights, bins=nbins, range=(0,5), label=r'$\Delta R$ (shower)', alpha=0.6)
ax[0].set_xlabel(r'$\Delta R(e,e)$')
ax[0].set_ylabel('events')
ax[0].set_title(r'$\Delta R(e,e)$')
ax[0].legend()

ax[1].hist(parton.pT_1, weights=parton.weights, bins=nbins, range=range, label=r'$p_T$ (parton)', alpha=0.6)
ax[1].hist(shower.pT_1, weights=shower.weights, bins=nbins, range=range, label=r'$p_T$ (shower)', alpha=0.6)
ax[1].set_xlabel(r'$p_T(e^+)$ [GeV]')
ax[1].set_ylabel('events')
ax[1].set_title(r'$p_T(e^+)$')
ax[1].legend()

ax[2].hist(parton.pT_2, weights=parton.weights, bins=nbins, range=range, label=r'$p_T$ (parton)', alpha=0.6)
ax[2].hist(shower.pT_2, weights=shower.weights, bins=nbins, range=range, label=r'$p_T$ (shower)', alpha=0.6)
ax[2].set_xlabel(r'$p_T(e^-)$ [GeV]')
ax[2].set_ylabel('events')
ax[2].set_title(r'$p_T(e^-)$')
ax[2].legend()

plt.tight_layout()
plt.savefig('plots/pp_ee.png', dpi=150, bbox_inches='tight')

