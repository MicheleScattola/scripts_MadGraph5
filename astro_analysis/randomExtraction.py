import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import scipy.stats as stat
import time
from math import sin, cos, acos

start = time.time()

EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'

if not os.path.exists(EXP_MAP) or not os.path.exists(DATA):
    raise FileNotFoundError('Input file(s) not found!')

parser = argparse.ArgumentParser(
  prog='randomExtraction',
  description='Random extraction from Auger exposure map')
parser.add_argument('-nTrials', '--trials', type=int, default=10000,
  help='Number of trials of independent extractions')
parser.add_argument('-dist','--distance', type=float, default=15.,
  help='Angular distance from astronomical source')
parser.add_argument('-pVal','--pValThreshold', type=float, default=0.05,
  help='Threshold p-value for detection')

args = parser.parse_args()
nTrials = args.trials
distance = args.distance
P_THRESHOLD = args.pValThreshold


def angular_distance(ra1, dec1, ra2, dec2):
  return np.degrees(acos(
    sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1 - ra2)
  ))

# analyze data from Auger
data = np.genfromtxt(DATA, names=True)
labels = data.dtype.names
nPositions = len(data)

print(f'[INFO] Extracting {nPositions} random positions in {nTrials} independent trials')

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
  for obj in data:
    theta = angular_distance(
      np.radians(obj['RA']), np.radians(obj['Dec']),
      np.radians(source["RA"]), np.radians(source["DEC"])
    )
    if theta <= distance:
      source_count += 1
  source["data_count"] = source_count

# analyze exposure map
map = hp.read_map(EXP_MAP)
probabilities = map/map.sum()

# plot exposure map and sources
hp.mollview(map, title='Exposure map and analyzed sources')
source_ra = [source["RA"] for source in sources]
source_dec = [source["DEC"] for source in sources]
hp.projscatter(source_ra, source_dec, lonlat=True, marker='*', color='red', s=50)
plt.savefig('MapAndSources.png')

print('[INFO] Starting randomized analysis.')

for i in range(nTrials):
  # generate the random normalized data points
  pixels = np.random.choice(len(map), size=nPositions, p=probabilities)
  event_theta, event_phi = hp.pix2ang(hp.get_nside(map), pixels)
  event_ra = np.degrees(event_phi)
  event_dec = 90.0 - np.degrees(event_theta)

  for source in sources:
    count = 0
    for ra, dec in zip(event_ra, event_dec):
      theta = angular_distance(
        np.radians(ra), np.radians(dec),
        np.radians(source["RA"]), np.radians(source["DEC"])
      )
      if theta <= distance:
        count += 1
    source["trial_counts"].append(count)


print(f'[INFO] Trials finished after {nTrials} iterations')

# analyze source and calculate p-value
found = []
for source in sources:

  trial_counts = np.array(source["trial_counts"])
  counts = source["data_count"]
  expected = np.mean(trial_counts)
  p_hit = float(expected / nPositions)
  pvalue_pre = stat.binomtest(k=counts, n=nPositions, p=p_hit, alternative='greater').pvalue
  pvalue_penalized = 1.0-pow(1.0-pvalue_pre,len(sources))
  print(f"Source: {source['name']}, Found: {source['data_count']}, Expected: {expected:.3f}, p-value={pvalue_pre:.5E}")
  if pvalue_penalized < P_THRESHOLD:
    found.append(source)
  
end = time.time()
print(f'[INFO] Elapsed {end - start:.2f} seconds')