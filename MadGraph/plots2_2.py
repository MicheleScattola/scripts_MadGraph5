import sys
import os

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_1_3'
bin_internal = os.path.join(folder, 'bin')
sys.path.insert(0, bin_internal)

import internal.lhe_parser as lhe_parser

# load events into parser
eventfile_path = os.path.join(folder, 'Events/run_01/unweighted_events.lhe.gz')
eventfile = lhe_parser.EventFile(eventfile_path)

# initialize event weights and kinematic variables
wgts = []
e_lep = []

# auxiliary function to find particle given pdgid
def find_particle(e, pdg):
    """
    return the first particle in event e with the
    pdg id given
    """
    for p in e:
        if p.pid == pdg: return p
    return False

###

for p in eventfile:
    e = find_particle(p,11)
    p4 = lhe_parser.FourMomentum(e)

    e_lep.append(p4.E)
    wgts.append(p.wgt)
    
    
# plotting
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.hist(e_lep, bins=40, weights=wgts)
ax.set_title(r"$\mu^- \rightarrow e^- \bar{\nu}_e \nu_\mu$ at $\sqrt{S}=13$ TeV")
ax.set_xlabel('Electron energy [GeV]')
ax.set_ylabel('Weighted events')
plt.savefig('plots/electron_e.png', dpi=150, bbox_inches='tight')


