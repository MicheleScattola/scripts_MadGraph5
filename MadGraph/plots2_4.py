import sys
import os
import numpy as np

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_2_4'
bin_internal = os.path.join(folder, 'bin')
sys.path.insert(0, bin_internal)

import internal.lhe_parser as lhe_parser

# load events into parser
scan = [140,150,160,170,180,190]
paths = []
events = []

for s in scan:
    file = 'Events/m' + str(s) + '/unweighted_events.lhe.gz'
    paths.append(file)
    events.append(os.path.join(folder, file))

# class for data
class MyData:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.events = lhe_parser.EventFile(path)
        self.delta_R = []
        self.weights = []

    def load(self, delta, wgt):
        self.delta_R.append(delta)
        self.weights.append(wgt)

# auxiliary function to find particle given pdgid
def find_particle(e, pdg):
    """
    return the first particle in event e with the
    pdg id given
    """
    for p in e:
        if p.pid == pdg: return p
    return False

# rapidity
def calc_rapidity(p):
    return 1./2. * np.log((p.E + p.pz)/(p.E - p.pz))
# azimuthal psi
def calc_phi(p):
    return np.atan2(p.py, p.px)

# Delta R
def calc_delta_R(a, b):

    diff_y = calc_rapidity(a) - calc_rapidity(b)
    diff_phi = calc_phi(a) - calc_phi(b)

    # normalize diff_phi between -pi and pi
    if abs(diff_phi) > np.pi :
        if diff_phi > 0:
            n = int(diff_phi/(2*np.pi) + 0.5)
            diff_phi -= 2 * np.pi * n
        if diff_phi < 0:
            n = int(0.5 - diff_phi/(2*np.pi))
            diff_phi += 2 * np.pi * n

    return np.sqrt(np.square(diff_y) + np.square(diff_phi))

###

datasets = {
    mtop: MyData(name=f'm{mtop}', path=os.path.join(folder, f'Events/m{mtop}/unweighted_events.lhe.gz'))
    for mtop in scan
}


for mtop in scan:
    sample = datasets[mtop]
    for ev in sample.events:
        t = find_particle(ev, 6)
        g = find_particle(ev, 21)
        if not t or not g:
            raise RuntimeError("[ERROR] no gluon or top found!!")

        p4_t = lhe_parser.FourMomentum(t)
        p4_g = lhe_parser.FourMomentum(g)
        sample.load(calc_delta_R(p4_t, p4_g), ev.wgt)

# plotting
import matplotlib.pyplot as plt

# color palette
colors = plt.get_cmap('tab10', len(scan))(range(len(scan)))

fig, ax = plt.subplots(
    2,
    3,
    figsize=(7, 5),
    sharex=True,
    sharey=True,
    gridspec_kw={"wspace": 0.0, "hspace": 0.0},
)
fig.suptitle(r'$\Delta R(g,t)$ in $e^+e^-\rightarrow t \bar{t} g$ at $\sqrt{s}=1$ TeV')
index = 0
for i in range(2):
    for j in range(3):
        sample = datasets[scan[index]]
        ax[i,j].hist(sample.delta_R,bins=20, weights=sample.weights, label=fr'$m_t={scan[index]}$ GeV', color=colors[index])
        ax[i,j].legend(loc=4)
        # Keep only outer tick labels in a shared-axis grid.
        ax[i,j].set_xlabel(r'$\Delta R(g,t)$')
        ax[i,j].set_ylabel('events')
        ax[i,j].label_outer()
        #ax[i,j].set_yscale('log')
        index += 1

plt.savefig('plots/deltaR_tg_scan.png', dpi=150, bbox_inches='tight')
