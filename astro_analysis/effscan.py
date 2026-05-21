import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import scipy.stats as stat
import time
import numba as nb
from numpy import acos, sin, cos

### FUNCTIONS

# parallelized distance counter
@nb.njit(parallel=True, fastmath=True)
def CalculateDistance(event_ra_deg, event_dec_deg, source_ra_deg, source_dec_deg):
  
  source_ra = np.deg2rad(source_ra_deg)
  source_dec = np.deg2rad(source_dec_deg)
  n = np.int64(len(event_ra_deg))
  theta = np.empty(n,dtype=np.float64)

  for i in nb.prange(len(event_ra_deg)):
    ra = np.deg2rad(event_ra_deg[i])
    dec = np.deg2rad(event_dec_deg[i])
    theta[i] = acos(
      sin(dec) * sin(source_dec) +
      cos(dec) * cos(source_dec) * cos(ra - source_ra)
    )

  return theta


# Pre-compute trial RA/Dec arrays outside the Numba function
# returns already sorted by distance
def generate_trial_coords(nTrials, nEvents, probabilities, nside):
  """Generate RA/Dec coordinates for all trials upfront"""
  trial_ra_arrays = np.empty((nTrials, nEvents), dtype=np.float64)
  trial_dec_arrays = np.empty((nTrials, nEvents), dtype=np.float64)
  trial_dist_arrays = np.empty((nTrials, nEvents), dtype=np.float64)
  
  for j in range(nTrials):
    if j!=0 and j%10000 == 0:
      print(f'[Gen - INFO] generated {j} trials')
    pixels = np.random.choice(nPositions, size=nEvents, p=probabilities)
    theta, phi = hp.pix2ang(nside, pixels)
    trial_ra_arrays[j] = np.degrees(phi)
    trial_dec_arrays[j] = 90.0 - np.degrees(theta)
    
    # Sort each trial's coordinates by distance
    trial_dist = CalculateDistance(trial_ra_arrays[j], trial_dec_arrays[j], source_ra_deg, source_dec_deg)
    trial_sort = trial_dist.argsort()
    trial_ra_arrays[j] = trial_ra_arrays[j][trial_sort]
    trial_dec_arrays[j] = trial_dec_arrays[j][trial_sort]
    trial_dist_arrays[j] = trial_dist[trial_sort]
  
  return trial_ra_arrays, trial_dec_arrays, trial_dist_arrays

# return counts in range of angles as array of size (angle)
@nb.njit(parallel=True)
def countsPerAngle(distance, angle_max):
  """Count events cumulatively for each angle threshold.
  Assumes distance array is sorted and in radians."""
  counts = np.zeros(angle_max, dtype=np.int64)
  for angle_deg in nb.prange(1, angle_max + 1):
    angle_rad = np.deg2rad(angle_deg)
    counts[angle_deg - 1] = np.searchsorted(distance, angle_rad)
  return counts



### MAIN
if __name__ == "__main__":

  EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
  DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'
  ANGLE_MAX = np.int64(40)

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

  # analyze data from Auger
  auger_data = np.genfromtxt(DATA, names=True)
  labels = auger_data.dtype.names
  auger_ra = np.array([obj["RA"] for obj in auger_data],dtype=np.float64)
  auger_dec = np.array([obj["Dec"] for obj in auger_data],dtype=np.float64)
  auger_e = np.array([obj["E"] for obj in auger_data],dtype=np.float64)
  auger_dist = CalculateDistance(auger_ra,auger_dec,source_ra_deg,source_dec_deg)
  # sort by distance
  auger_sort = auger_dist.argsort()
  auger_ra = auger_ra[auger_sort]
  auger_dec = auger_dec[auger_sort]
  auger_dist = auger_dist[auger_sort]
  auger_counts = countsPerAngle(auger_dist,ANGLE_MAX)

  exp_map = hp.read_map(EXP_MAP).astype(np.float64) 
  probabilities = exp_map / exp_map.sum()
  nside = hp.get_nside(exp_map)
  nPositions = np.int64(len(exp_map))
  nEvents = np.int64(len(auger_ra))

  # Generate trial coordinates ONCE before the loop
  print(f"[Gen - INFO] Generating {nTrials} trial coordinate sets (ordered by distance)")
  gen_start = time.time()
  trial_ra, trial_dec, trial_dist = generate_trial_coords(nTrials, nEvents, probabilities, nside)
  print(f"[Gen - INFO] Trial generation took {time.time() - gen_start:.2f} seconds\n")

  trial_counter = np.empty((trial_dist.shape[0],ANGLE_MAX),dtype=np.int64)
  print(f'[Trial - INFO] Analyzing trials from 0 to {ANGLE_MAX} degrees')
  trial_start = time.time()
  for j in range(trial_dist.shape[0]):
    trial_counter[j] = countsPerAngle(trial_dist[j], ANGLE_MAX)
  

  expected = trial_counter.mean(axis=0)
  pvalues = np.empty(ANGLE_MAX, dtype=np.float64)

  for i in range(ANGLE_MAX):
    null_counts = trial_counter[:, i]
    pvalues[i] = (np.count_nonzero(null_counts >= auger_counts[i]) + 1.0) / (nTrials + 1.0)
  print(f'[Trial - INFO] Trial analysis took {time.time() - trial_start:.2f} seconds\n')

  angle = [i+1 for i in range(ANGLE_MAX)]
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


