
# Bond Length Distribution Algorithm
# For generating random structures

if __name__ == "__main__":
    import sys
    sys.path.insert(0, sys.path[0] + '/..')

import numpy as np, re
import zipfile, os
from cluster.coval import Covalent

# a set of atomic coordinates
class Cluster(object):
    eps = np.finfo(dtype=np.float64).eps * 100
    def __init__(self, n, elems=None, atoms=None, energy=0.0, label=''):
        self.n = n
        self.atoms = np.zeros((n, 3), dtype=np.float64) if atoms is None else atoms
        self.elems = np.array([''] * n, dtype='|S20') if elems is None else elems
        self.energy = energy
        self.label = label
        self.mag = None
        self.surfnum = 0
        self.el_indices = None
        self.pgroup = None
        self.forces = None

    def pointgroup(self):
        return 'C1'

    # sample uniformly in the unit sphere
    def _random_direction(self):
        d = np.random.random(size=2)
        phi, theta = 2 * np.pi * d[0], np.arccos(2 * d[1] - 1)
        return np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta),
                         np.cos(theta)])

    def to_cluster(self):
        return self

    # return the position of the center
    def get_center(self):
        return np.average(self.atoms, axis=0)

    # will change the coordinates
    def center(self):
        self.atoms -= self.get_center()

    # do a rotation
    def rotate(self, rot_ang, dir_the, dir_phi):
        theta = rot_ang
        the = dir_the
        phi = dir_phi
        u = [np.sin(the) * np.cos(phi), np.sin(the) * np.sin(phi), np.cos(the)]
        r = np.array([[np.cos(theta) + u[0]**2 * (1 - np.cos(theta)),
                       u[0]*u[1] * (1 - np.cos(theta)) - u[2] * np.sin(theta),
                       u[0]*u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
                      [u[1]*u[0] * (1 - np.cos(theta)) + u[2] * np.sin(theta),
                       np.cos(theta) + u[1]**2 * (1 - np.cos(theta)),
                       u[1]*u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
                      [u[2]*u[0] * (1 - np.cos(theta)) - u[1] * np.sin(theta),
                       u[2]*u[1] * (1 - np.cos(theta)) + u[0] * np.sin(theta),
                       np.cos(theta) + u[2]**2 * (1 - np.cos(theta))]])
        self.atoms = self.atoms.dot(r.T)

    # get original data for measuring the bond length distribution
    def get_dist(self):
        ell = []
        eln = []
        for i in range(self.n):
            if self.elems[i] in eln:
                ell[eln.index(self.elems[i])].append(i)
            else:
                eln.append(self.elems[i])
                ell.append([i])
        idel = np.argsort(np.array([len(l) for l in ell]))[::-1]
        ell = [ell[i] for i in idel]
        eln = [eln[i] for i in idel]
        res = {}
        for i in range(len(ell)):
            for j in range(i + 1):
                x = np.linalg.norm(self.atoms[ell[i]].reshape((len(ell[i]), 1, 3)) -
                    self.atoms[ell[j]].reshape((1, len(ell[j]), 3)), axis=2)
                y = eln[i] + '-' + eln[j]
                yp = eln[j] + '-' + eln[i]
                x = np.sort(x, axis=1)
                if i == j: x = x[:, 1:]
                res[y], res[yp] = x, x
        return res

    # write coordinates to file
    def write_xyz(self, fn, append=True):
        f = open(fn, 'a' if append else 'w')
        if self.label == '':
            f.write('%d\nE = %15.8f\n' % (self.n, self.energy))
        elif self.surfnum != 0:
            f.write('%d\n%s = %15.8f SURF = %d\n' %
                    (self.n, self.label, self.energy, self.surfnum))
        elif self.energy is not None:
            f.write('%d\n%s = %15.8f\n' % (self.n, self.label, self.energy))
        else:
            f.write('%d\n%s\n' % (self.n, self.label))
        if self.forces is not None:
            for x, y, z in zip(self.elems, self.atoms, self.forces):
                f.write('%7s%15.8f%15.8f%15.8f %13.5f %13.5f %13.5f\n' %
                        (x, y[0], y[1], y[2], z[0], z[1], z[2]))
        else:
            for x, y in zip(self.elems, self.atoms):
                f.write('%7s%15.8f%15.8f%15.8f\n' % (x, y[0], y[1], y[2]))
        f.close()

    def shuffle(self):
        for i in self.el_indices:
            np.random.shuffle(self.atoms[i])

def write_clusters(fn, clus):
    append = False
    for c in clus:
        c.write_xyz(fn, append)
        append = True

# read xyz structures from xyz file
def read_clusters(fn, fnum=-1, zf=None, iprint=True, bufsize=0):
    clul = []
    if zf is None:
        f = open(fn, 'r', buffering=bufsize)
        fs = f.readlines()
        f.close()
    else:
        fs = zf.open(fn, 'rU').readlines()
    fsx = [f.strip() for f in fs]
    fs = [[g for g in f.strip().split(' ') if len(g) != 0] for f in fs]
    cc = 0
    while cc < len(fs) and (len(clul) < fnum or fnum == -1):
        if not fs[cc][0].isdigit():
            return [None]
        cn = int(fs[cc][0])
        clu = Cluster(cn)
        i = 0
        if len(fs[cc + 1]) == 1:
            clu.label = fs[cc + 1][0]
        if len(fs[cc + 1]) >= 3 and fs[cc + 1][1] == '=':
            clu.energy = float(fs[cc + 1][2])
            clu.label = fs[cc + 1][0]
            if clu.label.startswith("MAG"):
                clu.mag = float(clu.label.split("~")[0].split(":")[2])
        else:
            clu.label = fsx[cc + 1]
            clu.energy = None
        if "SURF" in fs[cc + 1]:
            clu.surfnum = int(fs[cc + 1][fs[cc + 1].index("SURF") + 2])
        if len(fs[cc + 2]) >= 7:
            clu.forces = np.zeros((cn, 3), dtype=float)
            for ixf, f in enumerate(fs[cc + 2:cc + 2 + cn]):
                if len(f) >= 7:
                    clu.forces[ixf] = [float(g) for g in f[4:7]]
        for f in fs[cc + 2:cc + 2 + cn]:
            clu.elems[i] = f[0]
            clu.atoms[i] = np.array([float(g) for g in f[1:4]])
            i = i + 1
        clu.el_indices = elem_indices(clu.elems)
        if "SURF" not in fs[cc + 1]:
            clu.center()
        clul.append(clu)
        cc += 2 + cn
    if iprint:
        print ('structs loaded: %d' % (len(clul), ))
    return clul

# read series of xyz structures from zip file
def read_zip_clusters(zipfn, last_only=False, iprint=True, bufsize=0):
    if not zipfile.is_zipfile(zipfn):
        return read_clusters(zipfn, iprint=iprint, bufsize=bufsize)
    clul = []
    zf = zipfile.ZipFile(zipfn, 'r')
    namel = zf.namelist()
    # sort filenames
    name_idx = [re.findall(r'[0-9]+', name) for name in namel]
    name_idx = [[int(i) for i in n] for n in name_idx]
    idx = sorted(range(len(name_idx)), key=name_idx.__getitem__)
    namel = [namel[i] for i in idx]
    for name in namel:
        xclul = read_clusters(name, zf=zf, iprint=False)
        if last_only:
            clul += xclul[-1:]
        else: clul += xclul
    ce = np.array([c.label for c in clul])
    if (ce == 'Energy').all():
        for i, c in enumerate(clul):
            c.label = 'ORI:{}'.format(i)
    if iprint:
        print ('structs loaded: %d' % (len(clul), ))
    return clul

# slices for elems of same kind
# used when shuffle
def elem_indices(elems):
    idx = []
    n = len(elems)
    i = 0
    while i < n:
        for j in range(i + 1, n + 1):
            if j == n or elems[j] != elems[i]:
                break
        idx.append(slice(i, j))
        i = j
    return idx

# elems list to number used by at_sort
def elem_num(elems):
    lel = []
    cel = []
    for i in elems:
        if i not in cel:
            cel.append(i)
            lel.append(len(cel))
        else:
            j = cel.index(i)
            lel.append(j + 1)
    return lel, len(cel)

def understand_name(name):
    if "|" in name:
        ex_names = name.split("|")
        elems = []
        fl_elems = []
        et = ""
        for en in ex_names:
            gelems, getx = elem_char(en)
            elems.append(gelems)
            fl_elems += list(gelems)
            et += " " + getx
        fl_elems = np.array(fl_elems)
        if len(et) != 0: et = et[1:]
    else:
        elems, et = elem_char(name)
        fl_elems = elems
    return fl_elems, elems, et

# re for resolving cluster name
rxa = r'^\s*([A-Za-z][a-z]*)\s*([0-9]+)(.*)$'
rxb = r'^\s*\(\s*([^\)]+)\s*\)\s*([0-9]+)(.*)$'
rxc = r'^\s*([A-Za-z][a-z]*)(.*)$'
rxbc = r'^\s*\(\s*([^\)]+)\s*\)(.*)$'
rxd = r'^\s*$'

# cluster name resolve
# return elems array and string
def elem_char(name):
    dc = []
    dct = []
    while len(re.findall(rxd, name)) == 0:
        ra = re.findall(rxa, name)
        if len(ra) == 0:
            ra = re.findall(rxb, name)
        if len(ra) != 0:
            ne = ra[0][0][:1].upper() + ra[0][0][1:]
            dc += [ne] * int(ra[0][1])
            dct += [ne, ra[0][1]]
            name = ra[0][2]
        else:
            ra = re.findall(rxc, name)
            if len(ra) == 0:
                ra = re.findall(rxbc, name)
            ne = ra[0][0][:1].upper() + ra[0][0][1:]
            dc += [ne]
            dct += [ne]
            name = ra[0][1]
    return np.array(dc), ' '.join(dct)

class MoleculePool(object):

    def __init__(self, x):
        self.mm = x
        self.ismulti = isinstance(self.mm, list)
        self.smode = "random"
        self.snum = 1
        self.scur = 0
    
    def copy(self):
        import copy
        return MoleculePool(copy.deepcopy(self.mm))
    
    def reject(self):
        self.scur -= 1
    
    def get_prop(self, prop):
        if isinstance(self.mm, list):
            if self.smode == "random":
                ix = np.random.randint(0, len(self.mm))
            else:
                ix = self.scur / self.snum
                self.scur += 1
            return getattr(self.mm[ix], prop)
        else:
            return getattr(self.mm, prop)

    def get_surf(self):
        return self.get_prop("surf")

    def get_atoms(self):
        return self.get_prop("atoms")
    
    def __getattr__(self, attr):
        if attr in [ "atoms", "surf", "mm", "ismulti", "smode", "snum", "scur" ]:
            return getattr(self, attr)
        else:
            return getattr(self.mm[0] if 
                isinstance(self.mm, list) else self.mm, attr)
    
    atoms = property(get_atoms)
    surf = property(get_surf)

# load molecules
libfn = os.path.dirname(__file__) + "/../library/molecules.zip"
zf = zipfile.ZipFile(libfn, 'r')
namel = zf.namelist()
moles = {}
for n in namel:
    if n.endswith('.xyz'):
        key = n[:-4]
        key = key[:1].upper() + key[1:]
        mm = read_clusters(n, zf=zf, iprint=False)
        if len(mm) == 1:
            moles[key] = MoleculePool(mm[0])
        else:
            moles[key] = MoleculePool(mm)
del libfn, zf, namel

# factor for determining second-order sigma
sigma_factor = 4

# use default covalent bonds if necessary
def update_stat(elems, xmean, xsigma, def_sigma):
    elems_ori = list(set(elems))
    elems = []
    for ee in elems_ori:
        if ee in moles:
            elems += [x if x != 'HX' else 'H' for x in moles[ee].elems]
            # if ee.endswith('ne'):
            #     elems += ['C', 'H']
            # else:
            #     elems += list(elem_char(ee.split('_')[0])[0])
        else:
            elems.append(ee)
    elems = list(set(elems))
    for i in range(len(elems)):
        for j in range(len(elems)):
            nn = elems[i] + '-' + elems[j]
            nx = elems[j] + '-' + elems[i]
            if nn in xmean and nx not in xmean:
                xmean[nx] = xmean[nn]
            if nn in xsigma and nx not in xsigma:
                xsigma[nx] = xsigma[nn]
    for i in range(len(elems)):
        for j in range(len(elems)):
            nn = elems[i] + '-' + elems[j]
            if nn in xmean.keys():
                if len(xmean[nn]) == 1: xmean[nn].append(xmean[nn][0])
            if nn not in xmean.keys() or len(xmean[nn]) == 0:
                if elems[i] in Covalent.x and elems[j] in Covalent.x:
                    l = Covalent.x[elems[i]] + Covalent.x[elems[j]]
                    xmean[nn] = [l, l]
            if nn in xsigma.keys():
                if len(xsigma[nn]) == 1: xsigma[nn].append(xsigma[nn][0] * 10)
            if nn not in xsigma.keys() or len(xsigma[nn]) == 0:
                xsigma[nn] = [def_sigma, def_sigma * sigma_factor]

if __name__ == "__main__":
    elem_char('Pt6(CO)')
