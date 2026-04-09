import numpy as np
import sys

# import lhe and hempc parsers
mg5_root = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7'
sys.path.insert(0, mg5_root)

import madgraph.various.hepmc_parser as hepmc_parser
import madgraph.various.lhe_parser as lhe_parser

# return particle given PDG id
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

# p_T
def calc_p_T(a):
    return np.sqrt(a.px*a.px + a.py*a.py)

# need different parsers for lhe and hepmc data
def open_event_file(path):
    if path.endswith('.hepmc') or path.endswith('.hepmc.gz'):
        return hepmc_parser.HEPMC_EventFile(path)
    return lhe_parser.EventFile(path)
