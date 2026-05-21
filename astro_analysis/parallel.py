import argparse
import os
import time
from math import acos, cos, sin

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat

try:
  import numba as nb
except ImportError:  # pragma: no cover
  nb = None

start = time.time()

EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'

if not os.path.exists(EXP_MAP) or not os.path.exists(DATA):
  raise FileNotFoundError('Input file(s) not found!')

parser = argparse.ArgumentParser(
  prog='parallel',
  description='Parallelized random extraction from Auger exposure map')
parser.add_argument('-nTrials', '--trials', type=int, default=10000,
  help='Number of trials of independent extractions')
parser.add_argument('-dist', '--distance', type=float, default=15.,
  help='Angular distance from astronomical source')
parser.add_argument('-pVal', '--pValThreshold', type=float, default=0.05,
  help='Threshold p-value for detection')
parser.add_argument('--no-plot', action='store_true',
  help='Disable saving the exposure map plot')

args = parser.parse_args()
n_trials = args.trials
distance_deg = args.distance
p_threshold = args.pValThreshold


def angular_distance(ra1, dec1, ra2, dec2):
  return np.degrees(acos(
    sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1 - ra2)
  ))


if nb is not None:

  @nb.njit(parallel=True, fastmath=True)
  def count_events_for_sources(event_ra_deg, event_dec_deg, source_ra_deg, source_dec_deg, distance_deg):
    n_sources = source_ra_deg.shape[0]
    n_events = event_ra_deg.shape[0]
    counts = np.zeros(n_sources, dtype=np.int64)
    distance_rad = np.deg2rad(distance_deg)

    for i in nb.prange(n_sources):
      source_ra = np.deg2rad(source_ra_deg[i])
      source_dec = np.deg2rad(source_dec_deg[i])
      count = 0
      for j in range(n_events):
        ra = np.deg2rad(event_ra_deg[j])
        dec = np.deg2rad(event_dec_deg[j])
        angle = acos(
          sin(dec) * sin(source_dec) +
          cos(dec) * cos(source_dec) * cos(ra - source_ra)
        )
        if angle <= distance_rad:
          count += 1
      counts[i] = count

    return counts

else:

  def count_events_for_sources(event_ra_deg, event_dec_deg, source_ra_deg, source_dec_deg, distance_deg):
    counts = np.zeros(len(source_ra_deg), dtype=np.int64)
    for i in range(len(source_ra_deg)):
      count = 0
      for ra, dec in zip(event_ra_deg, event_dec_deg):
        theta = angular_distance(
          np.radians(ra), np.radians(dec),
          np.radians(source_ra_deg[i]), np.radians(source_dec_deg[i])
        )
        if theta <= distance_deg:
          count += 1
      counts[i] = count
    return counts


# analyze data from Auger
auger_data = np.genfromtxt(DATA, names=True)
observed_n_events = len(auger_data)
print(f'[INFO] Extracting {observed_n_events} random positions in {n_trials} independent trials')

sources = [
  {"name": "Sagittarius A*", "RA": 266.4168, "DEC": -29.0078, "data_count": 0, "trial_counts": []},
  {"name": "Centaurus A", "RA": 201.3, "DEC": -43.0, "data_count": 0, "trial_counts": []},
  {"name": "Fornax A", "RA": 50.67, "DEC": -37.2, "data_count": 0, "trial_counts": []},
  {"name": "NGC 253", "RA": 11.89, "DEC": -25.29, "data_count": 0, "trial_counts": []},
  {"name": "M83", "RA": 253.47, "DEC": -24.38, "data_count": 0, "trial_counts": []},
  {"name": "M87 (Virgo A)", "RA": 187.7, "DEC": 12.39, "data_count": 0, "trial_counts": []},
  {"name": "Vela SNR", "RA": 128.4, "DEC": -45.18, "data_count": 0, "trial_counts": []},
  {"name": "LMC", "RA": 80.89, "DEC": -69.76, "data_count": 0, "trial_counts": []},
  {"name": "SMC", "RA": 13.16, "DEC": -72.8, "data_count": 0, "trial_counts": []},
]

for source in sources:
  source_count = 0
  for obj in auger_data:
    theta = angular_distance(
      np.radians(obj['RA']), np.radians(obj['Dec']),
      np.radians(source['RA']), np.radians(source['DEC'])
    )
    if theta <= distance_deg:
      source_count += 1
  source['data_count'] = source_count

exp_map = hp.read_map(EXP_MAP)
probabilities = exp_map / exp_map.sum()

if not args.no_plot:
  hp.mollview(exp_map, title='Exposure map and analyzed sources')
  source_ra = [source['RA'] for source in sources]
  source_dec = [source['DEC'] for source in sources]
  hp.projscatter(source_ra, source_dec, lonlat=True, marker='*', color='red', s=50)
  plt.savefig('MapAndSources.png')

source_ra_deg = np.array([source['RA'] for source in sources], dtype=np.float64)
source_dec_deg = np.array([source['DEC'] for source in sources], dtype=np.float64)

print('[INFO] Starting randomized analysis.')

for trial_index in range(n_trials):
  pixels = np.random.choice(len(exp_map), size=observed_n_events, p=probabilities)
  event_theta, event_phi = hp.pix2ang(hp.get_nside(exp_map), pixels)
  event_ra_deg = np.degrees(event_phi)
  event_dec_deg = 90.0 - np.degrees(event_theta)

  trial_counts = count_events_for_sources(
    event_ra_deg.astype(np.float64),
    event_dec_deg.astype(np.float64),
    source_ra_deg,
    source_dec_deg,
    distance_deg,
  )

  for source_index, source in enumerate(sources):
    source['trial_counts'].append(int(trial_counts[source_index]))

print(f'[INFO] Trials finished after {n_trials} iterations')

found = []
for source in sources:
  trial_counts = np.array(source['trial_counts'])
  expected = float(np.mean(trial_counts))
  p_hit = float(expected / observed_n_events)
  pvalue_pre = stat.binomtest(
    k=source['data_count'],
    n=observed_n_events,
    p=p_hit,
    alternative='greater'
  ).pvalue
  pvalue_penalized = 1.0 - pow(1.0 - pvalue_pre, len(sources))
  print(f"Source: {source['name']}, Found: {source['data_count']}, Expected: {expected:.3f}, p-value={pvalue_pre:.5E}")
  if pvalue_penalized < p_threshold:
    found.append(source)

end = time.time()
print(f'[INFO] Elapsed {end - start:.2f} seconds')
