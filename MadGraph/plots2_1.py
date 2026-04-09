import sys
import os
import numpy as np

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_2_1'
bin_internal = os.path.join(folder, 'bin')
sys.path.insert(0, bin_internal)

import internal.lhe_parser as lhe_parser

# load events into parser
eventfile_path = os.path.join(folder, 'Events/13TeV/unweighted_events.lhe.gz')
eventfile = lhe_parser.EventFile(eventfile_path)
eventfile_path2 = os.path.join(folder, 'Events/2TeV/unweighted_events.lhe.gz')
eventfile2 = lhe_parser.EventFile(eventfile_path2)


wgts = []
y = []

wgts2 = []
y2 = []

# auxiliary function to find particle given pdgid
def find_particle(e, pdg):
    """
    return the first particle in event e with the
    pdg id given
    """
    for p in e:
        if p.pid == pdg: return p
    return False

# auxiliary function to return rapidity given a particle
def calc_rapidity(p):
    return 1./2. * np.log((p.E + p.pz)/(p.E - p.pz))


# auxiliary function for the limit
def limit(S,M):
    return np.log(S/M)

###

for e in eventfile:
    z = find_particle(e,23)
    p4_z = lhe_parser.FourMomentum(z)

    y.append(calc_rapidity(p4_z))
    wgts.append(e.wgt)
    
for e in eventfile2:
    z = find_particle(e,23)
    p4_z = lhe_parser.FourMomentum(z)

    y2.append(calc_rapidity(p4_z))
    #print(p4_z)
    wgts2.append(e.wgt)
    
# plotting
import matplotlib.pyplot as plt

fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

ax_left.hist(y, bins=20, range=(-5, 5), weights=wgts, label='data')
ax_left.set_title(r'$p p \rightarrow Z$ at $\sqrt{S}=13$ TeV')
ax_left.set_xlabel('Z Rapidity: y')
ax_left.set_ylabel('Weighted events')

ax_right.hist(y2, bins=20, range=(-5, 5), weights=wgts2, label='data')
ax_right.set_title(r'$p p \rightarrow Z$ at $\sqrt{S}=2$ TeV')
ax_right.set_xlabel('Z Rapidity: y')
ax_right.set_ylabel('Weighted events')

limit1 = limit(13000,91)
limit2 = limit(2000,91)

ax_left.axvline(limit1, color='red', linestyle=':', lw=2, label=r'$|y(z)| = \log \frac{S}{m^2_Z}$')
ax_left.axvline(-limit1, color='red', linestyle=':', lw=2)
ax_left.legend()

ax_right.axvline(limit2, color='red', linestyle=':', lw=2, label=r'$|y(z)| = \log \frac{S}{m^2_Z}$')
ax_right.axvline(-limit2, color='red', linestyle=':', lw=2)
ax_right.legend()

fig.tight_layout()
plt.savefig('plots/rapidity_comparison.png', dpi=150, bbox_inches='tight')

print(len(y), len(y2))
print(sum(wgts), sum(wgts2))
print(min(y), max(y))
print(min(y2), max(y2))


