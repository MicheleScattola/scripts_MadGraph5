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
  '''Pre-computing astronomical distance between objects.
     Expects data in degrees.'''
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
def trialDistances(nEvents, probabilities, nside):
  """return array of distances for each trial"""
  trial_ra_arrays = np.empty(nEvents, dtype=np.float64)
  trial_dec_arrays = np.empty(nEvents, dtype=np.float64)
  
  pixels = np.random.choice(nPositions, size=nEvents, p=probabilities)
  trial_ra_arrays, trial_dec_arrays = hp.pix2ang(nside, pixels, lonlat=True)
    
  trial_dist = CalculateDistance(trial_ra_arrays, trial_dec_arrays, source_ra_deg, source_dec_deg).argsort()

  return trial_dist

# return counts in range of angles as array of size (angle)
@nb.njit(parallel=True)
def countsPerAngle(distance, angle_max):
  '''Count events for a given angle.
  Assumes distance array is sorted and in radians.'''
  counts = np.zeros(angle_max, dtype=np.int64)
  for angle_deg in nb.prange(1, angle_max + 1):
    angle_rad = np.deg2rad(angle_deg)
    counts[angle_deg - 1] = np.searchsorted(distance, angle_rad)
  return counts



### MAIN
if __name__ == "__main__":

  EXP_MAP = '/home/mike/Physics/labParticles/data_files/exposure.fits'
  DATA = '/home/mike/Physics/labParticles/data_files/auger.txt'

  if not os.path.exists(EXP_MAP) or not os.path.exists(DATA):
      raise FileNotFoundError('Input file(s) not found!')

  parser = argparse.ArgumentParser(
    prog='trialScan',
    description='Parallelized scan of give astronomical source.')
  parser.add_argument('-ra', '--source_ra', type=float, default=201.3,
    help='Astronomical source Right Ascension (degrees)')
  parser.add_argument('-dec', '--source_dec', type=float, default=-43.0,
    help='Astronomical source Declination (degrees)')
  parser.add_argument('-n', '--nTrials', type=int, default=1000000,
    help='Number of Monte Carlo trials (int)')
  parser.add_argument('-hat', '--topHat', type=int, default=40,
    help='Max angle of top hat region')
  parser.add_argument('-pval', '--targetPvalue', type=float, default=7.99E-6,
    help='Minimum p-value to search for')

  args = parser.parse_args()
  source_ra_deg = np.float64(args.source_ra)
  source_dec_deg = np.float64(args.source_dec)
  nTrials = np.int64(args.nTrials)
  angleMax = np.int64(args.topHat)
  pvalueMin = np.float64(args.targetPvalue)

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
  auger_counts = countsPerAngle(auger_dist,angleMax)

  exp_map = hp.read_map(EXP_MAP).astype(np.float64)
  rotator = hp.Rotator(coord=['G', 'C'])
  exp_map_eq = rotator.rotate_map_pixel(exp_map)
  probabilities = exp_map_eq / exp_map_eq.sum()
  nSide = hp.get_nside(exp_map_eq)
  nPositions = np.int64(len(exp_map_eq))
  nEvents = np.int64(len(auger_ra)) 

  isPvalLower = np.int64(0)

  start = time.time()
  for i in range(nTrials):
    if i%1E3 == 0 and i!=0:
      print(f'[INFO] simulating trial {i}/{nTrials}, elapsed {time.time() - start:.2f}')
      start = time.time()

    trial_distances = trialDistances(nEvents,probabilities,nSide)
    trial_counts = countsPerAngle(trial_distances,angleMax)
    angles = [i+1 for i in range(angleMax)]
    
    pvalues = np.zeros(angleMax,dtype=np.float64)
    for j in range(angleMax):
      pvalues[j] = stat.binomtest(
                                k=int(auger_counts[j]),
                                n=nEvents,
                                p=trial_counts[j]/nEvents,
                                alternative='greater'
                              ).pvalue
      
    p_local_min = np.min(pvalues)
    if p_local_min < pvalueMin:
      isPvalLower += 1

  print(f'Computed {nTrials} trials.\nFound {isPvalLower} trials with smaller p-value with respect to test {pvalueMin:.2E}')
  