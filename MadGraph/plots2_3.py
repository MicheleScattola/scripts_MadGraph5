import sys
import os
from functions import lhe_parser

# path to exercise folder
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/es_1_5'
#bin_internal = os.path.join(folder, 'bin')
#sys.path.insert(0, bin_internal)

#import internal.lhe_parser as lhe_parser

# load events into parser
signal_path = os.path.join(folder, 'sig/Events/signal/unweighted_events.lhe.gz')
signal_events = lhe_parser.EventFile(signal_path)

background_path = os.path.join(folder, 'bkg/Events/background/unweighted_events.lhe.gz')
background_events = lhe_parser.EventFile(background_path)

# initialize event weights and variables
wgts_sig = []
m_sig = []

wgts_bkg = []
m_bkg = []

my_data = [
    {
        "type": "signal",
        "events": signal_events,
        "mass": m_sig,
        "weights": wgts_sig,
    },
    {
        "type": "background",
        "events": background_events,
        "mass": m_bkg,
        "weights": wgts_bkg,
    },
]


def diphoton_mass(event):
    """Return invariant mass of the first two final-state photons in one event."""
    photons = [p for p in event if p.status == 1 and p.pid == 22]
    if len(photons) != 2:
        return None

    p4_1 = lhe_parser.FourMomentum(photons[0])
    p4_2 = lhe_parser.FourMomentum(photons[1])
    return (p4_1 + p4_2).mass

# Single pass over both samples to fill masses and weights
for sample in my_data:
    for ev in sample["events"]:
        m_gg = diphoton_mass(ev)
        if m_gg is None:
            continue

        sample["mass"].append(m_gg)
        sample["weights"].append(ev.wgt)
    
# plotting
import matplotlib.pyplot as plt

nbins = 50
M_min = 100
M_max = 150

invariant_masses = [m_sig,m_bkg]
ev_weights = [wgts_sig, wgts_bkg]
colors = ['red','blue']
labels = [r'$H \rightarrow \gamma \gamma$','background']

fig, ax = plt.subplots()
ax.hist(
    invariant_masses,
    bins=nbins,
    range=(M_min, M_max),
    weights=ev_weights,
    label=labels,
    color=colors,
    histtype='step',
)
ax.set_title(r'$\gamma \gamma$ events at $\sqrt{S} = 13$ TeV')
ax.set_xlabel(r'$m_{\gamma\gamma}$ [GeV]')
ax.set_ylabel(r'$\sigma$ [pb]')
ax.legend(title=r'cut: $p_T (\gamma) > 20$ GeV')

plt.savefig('plots/invariant_H.png', dpi=150, bbox_inches='tight')



