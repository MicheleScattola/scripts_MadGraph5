import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import scipy.stats as stat
import time
import numba as nb
from math import sin, cos, acos

EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'

if not os.path.exists(EXP_MAP) or not os.path.exists(DATA):
    raise FileNotFoundError('Input file(s) not found!')

parser = argparse.ArgumentParser(
  prog='parallelScan',
  description='Parallelized scan of give astronomical source.')
parser.add_argument('-ra', '--source_ra', type=float, default=201.3,
  help='Astronomical source Right Ascension (degrees)')
parser.add_argument('-dec', '--source_dec', type=float, default=-43.0,
  help='Astronomical source Declination (degrees)')
parser.add_argument('-n', '--nTrials', type=int, default=10000,
  help='Number of Monte Carlo trials (int)')

args = parser.parse_args()
source_ra_deg = np.float64(args.source_ra)
source_dec_deg = np.float64(args.source_dec)
nTrials = np.int64(args.nTrials)

# parallelized events counter
@nb.njit(parallel=True, fastmath=True)
def count_events_for_sources(event_ra_deg, event_dec_deg, source_ra_deg, source_dec_deg, distance_deg, n_events):

  counts = np.int64(0)
  distance_rad = np.deg2rad(distance_deg)
  source_ra = np.deg2rad(source_ra_deg)
  source_dec = np.deg2rad(source_dec_deg)

  for i in nb.prange(n_events):
    ra = np.deg2rad(event_ra_deg[i])
    dec = np.deg2rad(event_dec_deg[i])
    angle = acos(
      sin(dec) * sin(source_dec) +
      cos(dec) * cos(source_dec) * cos(ra - source_ra)
    )
    if angle <= distance_rad:
      counts += 1

  return counts

# parallelized random trial from RA and DEC arrays
@nb.njit(parallel=True, fastmath=True)
def random_trial(nTrials, nEvents, trial_ra_arrays, trial_dec_arrays,
                 source_ra_deg, source_dec_deg, distance_deg):
  """
  Expects pre-computed RA/Dec arrays for each trial
  """
  counts = np.empty(nTrials, np.float64)
  for j in nb.prange(nTrials):
    trial_counts = count_events_for_sources(
      trial_ra_arrays[j],
      trial_dec_arrays[j],
      source_ra_deg,
      source_dec_deg,
      distance_deg=distance_deg,
      n_events=nEvents
    )
    counts[j] = trial_counts
  return counts


# analyze data from Auger
auger_data = np.genfromtxt(DATA, names=True)
labels = auger_data.dtype.names
auger_ra = np.array([obj["RA"] for obj in auger_data],dtype=np.float64)
auger_dec = np.array([obj["Dec"] for obj in auger_data],dtype=np.float64)
auger_e = np.array([obj["E"] for obj in auger_data],dtype=np.float64)

exp_map = hp.read_map(EXP_MAP).astype(np.float64)  # Ensure native byte order
probabilities = exp_map / exp_map.sum()
nside = hp.get_nside(exp_map)
nPositions = np.int64(len(exp_map))

# Pre-compute trial RA/Dec arrays outside the Numba function
def generate_trial_coords(nTrials, nEvents, probabilities, nside):
  """Generate RA/Dec coordinates for all trials upfront"""
  trial_ra_arrays = np.empty((nTrials, nEvents), dtype=np.float64)
  trial_dec_arrays = np.empty((nTrials, nEvents), dtype=np.float64)
  
  for j in range(nTrials):
    pixels = np.random.choice(nPositions, size=nEvents, p=probabilities)
    theta, phi = hp.pix2ang(nside, pixels)
    trial_ra_arrays[j] = np.degrees(phi)
    trial_dec_arrays[j] = 90.0 - np.degrees(theta)
  
  return trial_ra_arrays, trial_dec_arrays

# Generate trial coordinates ONCE before the loop
print(f"Generating {nTrials} trial coordinate sets (this may take a moment)...")
gen_start = time.time()
trial_ra_arrays, trial_dec_arrays = generate_trial_coords(nTrials, len(auger_ra), probabilities, nside)
gen_end = time.time()
print(f"Trial generation took {gen_end - gen_start:.2f} seconds\n")

# scan angles and analyze
pvalues = []
angle = []
start = time.time()
for i in range(40):
  # count events in data
  auger_counts = count_events_for_sources(
    auger_ra.astype(np.float64),
    auger_dec.astype(np.float64),
    source_ra_deg,
    source_dec_deg,
    distance_deg=i+1,
    n_events=np.int64(len(auger_ra))
  )
  # Use pre-computed trial coordinates
  trial_counts = random_trial(
    nTrials=nTrials,
    nEvents=np.int64(len(auger_ra)),
    trial_ra_arrays=trial_ra_arrays,
    trial_dec_arrays=trial_dec_arrays,
    source_ra_deg=source_ra_deg,
    source_dec_deg=source_dec_deg,
    distance_deg=i+1,
  )

  expected = float(np.mean(trial_counts))
  p_hit = float(expected / nPositions)
  pvalue_local = stat.binomtest(
    k=auger_counts,
    n=len(auger_ra),
    p=p_hit,
    alternative='greater'
  ).pvalue
  print(f'Angle: {i+1}, Found: {auger_counts}, Expected: {expected}, p-value: {pvalue_local:.2E}')

  pvalues.append(pvalue_local)
  angle.append(i+1)

  print(f'Calculated {nTrials} trials')

  if i%5==0 and i!=1:
    end = time.time()
    print(f'[INFO] analyzed radius {i}/40 degrees, elapsed : {end - start:.2f} seconds')

plt.figure()
plt.plot(angle,pvalues,label='pvalues vs angle')
plt.ylabel('p-value local')
plt.xlabel('Angle [deg]')
plt.yscale('log')
plt.grid(True, which='both', alpha=0.3)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('pvalueVSangle.png')