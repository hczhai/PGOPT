
import numpy as np
from cluster.base import read_zip_clusters, Cluster

def inertia(c):
  ine = np.zeros((3, 3))
  for i in range(3):
    for j in range(i + 1):
      ine[i, j] = ine[j, i] = np.dot(c.atoms[:, i], c.atoms[:, j])
  return np.linalg.eigh(ine)

def reorient(c):
  _, v = inertia(c)
  v = v.T[::-1]
  xc = np.tensordot(c.atoms, v, axes=[1, 1])
  if np.dot(xc[:, 0], abs(xc[:, 0])).sum() < 0:
    v = -v
  c.atoms = np.tensordot(c.atoms, v, axes=[1, 1])
  return v

def inertia_2d(c):
  ine = np.zeros((2, 2))
  for i in range(2):
    for j in range(i + 1):
      ine[i, j] = ine[j, i] = np.dot(c.atoms[:, i], c.atoms[:, j])
  return np.linalg.eigh(ine)

def reorient_2d(c):
  _, v = inertia_2d(c)
  v = v.T[::-1]
  xc = np.tensordot(c.atoms[:, 0:2], v, axes=[1, 1])
  if np.dot(xc[:, 0], abs(xc[:, 0])).sum() < 0:
    v = -v
  return np.array([[v[0, 0], v[0, 1], 0.0], [v[1, 0], v[1, 1], 0.0], [0.0, 0.0, 1.0]])

class Align(object):

  def __init__(self, input_file=None):
    if input_file is None:
      self.clus = []
    else:
      self.clus = read_zip_clusters(input_file)
  
  def align(self):
    for c in self.clus:
      if c.surfnum == 0:
        c.center()
        reorient(c)
      else:
        dc = c.atoms[c.surfnum:] - c.atoms[c.surfnum:].mean(axis=0).reshape((1, 3))
        v = reorient_2d(Cluster(n=len(dc), atoms=dc))
        c.center()
        c.atoms = np.tensordot(c.atoms, v, axes=[1, 1])
