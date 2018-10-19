
import numpy as np
from cluster.coval import AtomicWeight
from cluster.base import Cluster

def read_frequencies(fn):
    f = open(fn, 'r')
    fs = f.readlines()
    f.close()
    fs = [[g for g in f.replace('\n', '').split(' ') if len(g) != 0] for f in fs]
    fs = fs[3:-1]
    return [float(f[2]) for f in fs]

def read_displacements(fn):
    f = open(fn, 'r')
    fs = f.readlines()
    f.close()
    fs = [[g for g in f.replace('\n', '').split(' ') if len(g) != 0] for f in fs]
    fs = fs[1:-1]
    prev = 0
    r = []
    for f in fs:
        if int(f[0]) != prev: r.append([])
        r[-1] += map(float, f[2:])
        prev = int(f[0])
    return np.array(r).T.tolist() # diff row-diff freqs, diff column-diff atom coords

def bin_next(ki):
    bb = False
    for i in range(len(ki))[::-1]:
        if ki[i] == 0:
            bb = True
            ki[i] = 1
            for j in range(i + 1, len(ki)):
                ki[j] = 0
            break
    return bb

def bin_tostr(ki):
    return ''.join(map(lambda x: '+' if x == 1 else '-', ki))

def freqs_net(freqs):
    ff = sorted(freqs, key=np.abs)
    return sorted(ff[6:])

def imag(freqs, thres=None):
    k = 0
    if thres is None:
        ff = sorted(freqs, key=np.abs)
        tt = np.abs(np.array(ff[:6])).max() * 2
        k = 0
        for i in ff[6:]:
            if i < -tt: k += 1
        return k
    else:
        for i in freqs:
            if i < thres: k += 1
        return k

def imag_surface(freqs, thres=None):
    k = 0
    if thres is None:
        thres = 0
    for i in freqs:
        if i < thres: k += 1
    return k

# step in unit angstrom
def make_displacements(clu, freqs, disps, step=0.5, thres=None):
    k = imag(freqs, thres)
    if k == 0: return []
    ki = np.zeros(k)
    rr = []
    ws = np.array(map(lambda x: AtomicWeight.x[x], clu.elems))
    ws = np.zeros(clu.atoms.shape) + 1.0 / np.sqrt(ws * 1822.88853).reshape(clu.n, 1)
    while True:
        r = Cluster(clu.n, elems=clu.elems)
        r.atoms = np.array(clu.atoms)
        for ii, i in enumerate(ki):
            fac = 1 if i == 1 else -1
            dp = disps[ii].reshape(clu.atoms.shape) * ws
            dp = dp / np.linalg.norm(dp)
            vec = step * dp
            r.atoms += vec * fac
        r.tid = clu.tid
        r.label = "TRA:%d%s" % (r.tid, bin_tostr(ki))
        r.energy = clu.energy
        rr.append(r)
        if not bin_next(ki): break
    return rr
