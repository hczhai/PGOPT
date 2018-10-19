
import random, numpy as np
from cluster.base import read_zip_clusters

class Stat(object):
  def __init__(self, input_file, sample_number=-1, energy_cut=None, 
    shuffle_input=False, last_only=False):

    if input_file is None: return
    self.clus = read_zip_clusters(input_file, last_only=last_only)
    if energy_cut is not None:
      self.clus = [c for c in self.clus if c.energy < energy_cut]
    if shuffle_input: random.shuffle(self.clus)
    if sample_number != -1:
      self.clus = self.clus[:sample_number]
    print ('%d structures selected. ' % (len(self.clus), ))
  
  def stat(self):
    cg = [c.get_dist() for c in self.clus]
    keys = cg[0].keys()
    eners = np.array([c.energy for c in self.clus])
    emax, emin = eners.max(), eners.min()
    distx = { "mu": {}, "sigma": {} }
    for k in keys:
      cgr = [ci[k] for ci in cg]
      dist = np.concatenate(cgr, axis=0)[:, 0:2]
      distx["mu"][k] = []
      distx["sigma"][k] = []
      for i in range(dist.shape[1]):
        distx["mu"][k].append(np.average(dist[:, i]))
        distx["sigma"][k].append(np.std(dist[:, i]))
    self.dist = dist
    return distx, emax, emin