import argparse
import os
import time

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat

start = time.time()

try:
  import numba as nb
except ImportError:  # pragma: no cover
  nb = None

EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'

if not os.path.exists(EXP_MAP) or not os.path.exists(DATA):
  raise FileNotFoundError('Input file(s) not found!')

parser = argparse.ArgumentParser(
  prog='eff',
  description='Efficient random extraction from Auger exposure map')
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


SOURCE_CATALOG = [
  {"name": "Sagittarius A*", "RA": 266.4168, "DEC": -29.0078},
  {"name": "Centaurus A", "RA": 201.3, "DEC": -43.0},
  {"name": "Fornax A", "RA": 50.67, "DEC": -37.2},
  {"name": "NGC 253", "RA": 11.89, "DEC": -25.29},
  {"name": "M83", "RA": 253.47, "DEC": -24.38},
  {"name": "M87 (Virgo A)", "RA": 187.7, "DEC": 12.39},
  {"name": "Vela SNR", "RA": 128.4, "DEC": -45.18},
  {"name": "LMC", "RA": 80.89, "DEC": -69.76},
  {"name": "SMC", "RA": 13.16, "DEC": -72.8},
]


def spherical_to_unit_vectors(ra_deg, dec_deg):
  ra_rad = np.deg2rad(ra_deg)
  dec_rad = np.deg2rad(dec_deg)
  cos_dec = np.cos(dec_rad)
  x = cos_dec * np.cos(ra_rad)
  y = cos_dec * np.sin(ra_rad)
  z = np.sin(dec_rad)
  return np.column_stack((x, y, z)).astype(np.float64)


if nb is not None:

  @nb.njit(parallel=True, fastmath=True)
  def count_hits(event_vectors, source_vectors, cos_distance):
    n_sources = source_vectors.shape[0]
    n_events = event_vectors.shape[0]
    counts = np.zeros(n_sources, dtype=np.int64)

    for i in nb.prange(n_sources):
      sx = source_vectors[i, 0]
      sy = source_vectors[i, 1]
      sz = source_vectors[i, 2]
      count = 0
      for j in range(n_events):
        dot = (
          event_vectors[j, 0] * sx +
          event_vectors[j, 1] * sy +
          event_vectors[j, 2] * sz
        )
        if dot >= cos_distance:
          count += 1
      counts[i] = count

    return counts

else:

  def count_hits(event_vectors, source_vectors, cos_distance):
    counts = np.zeros(source_vectors.shape[0], dtype=np.int64)
    for i, source_vector in enumerate(source_vectors):
      dots = event_vectors @ source_vector
      counts[i] = int(np.sum(dots >= cos_distance))
    return counts


print('[INFO] Loading Auger data.')

auger_data = np.genfromtxt(DATA, names=True)
observed_n_events = len(auger_data)
print(f'[INFO] Extracting {observed_n_events} random positions in {n_trials} independent trials')

source_ra = np.array([source['RA'] for source in SOURCE_CATALOG], dtype=np.float64)
source_dec = np.array([source['DEC'] for source in SOURCE_CATALOG], dtype=np.float64)
source_vectors = spherical_to_unit_vectors(source_ra, source_dec)

observed_event_vectors = spherical_to_unit_vectors(auger_data['RA'], auger_data['Dec'])

exp_map = hp.read_map(EXP_MAP)
probabilities = exp_map / exp_map.sum()

if not args.no_plot:
  hp.mollview(exp_map, title='Exposure map and analyzed sources')
  hp.projscatter(source_ra, source_dec, lonlat=True, marker='*', color='red', s=50)
  plt.savefig('MapAndSources.png')

cos_distance = np.cos(np.deg2rad(distance_deg))

print('[INFO] Counting observed Auger events.')
observed_counts = count_hits(observed_event_vectors, source_vectors, cos_distance)

trial_counts = np.zeros((len(SOURCE_CATALOG), n_trials), dtype=np.int64)
print('[INFO] Starting randomized analysis.')

for trial_index in range(n_trials):
  pixels = np.random.choice(len(exp_map), size=observed_n_events, p=probabilities)
  event_theta, event_phi = hp.pix2ang(hp.get_nside(exp_map), pixels)
  event_ra_deg = np.degrees(event_phi)
  event_dec_deg = 90.0 - np.degrees(event_theta)
  event_vectors = spherical_to_unit_vectors(event_ra_deg, event_dec_deg)

  trial_counts[:, trial_index] = count_hits(event_vectors, source_vectors, cos_distance)

print(f'[INFO] Trials finished after {n_trials} iterations')

found = []
for source_index, source in enumerate(SOURCE_CATALOG):
  source_trials = trial_counts[source_index]
  expected = float(np.mean(source_trials))
  p_hit = float(expected / observed_n_events)
  pvalue_pre = stat.binomtest(
    k=int(observed_counts[source_index]),
    n=observed_n_events,
    p=p_hit,
    alternative='greater'
  ).pvalue
  pvalue_penalized = 1.0 - pow(1.0 - pvalue_pre, len(SOURCE_CATALOG))
  print(
    f"Source: {source['name']}, Found: {int(observed_counts[source_index])}, "
    f"Expected: {expected:.3f}, p-value={pvalue_pre:.5E}"
  )
  if pvalue_penalized < p_threshold:
    found.append(source)

if found:
  print('[INFO] Significant sources:')
  for source in found:
    print(f"  - {source['name']}")
else:
  print('[INFO] No sources passed the threshold')

end = time.time()
print(f'[INFO] Elapsed {end - start:.2f} seconds')
