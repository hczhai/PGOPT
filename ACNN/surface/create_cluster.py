
from __future__ import print_function
if __name__ == "__main__":
    import sys
    sys.path.insert(0, sys.path[0] + '/..')

import numpy as np
from cluster.base import moles, Cluster
from surface.base import ClusterAtSurface
from surface.surf_symm import to_cellmat
from math import sin, cos


def rotation(atoms, rot_ang, dir_the, dir_phi):
    theta = rot_ang
    the = dir_the
    phi = dir_phi
    u = [sin(the) * cos(phi), sin(the) * sin(phi), cos(the)]
    r = np.array([[cos(theta) + u[0]**2 * (1 - cos(theta)),
                   u[0] * u[1] * (1 - cos(theta)) - u[2] * sin(theta),
                   u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
                  [u[1] * u[0] * (1 - cos(theta)) + u[2] * sin(theta),
                   cos(theta) + u[1]**2 * (1 - cos(theta)),
                   u[1] * u[2] * (1 - cos(theta)) - u[0] * sin(theta)],
                  [u[2] * u[0] * (1 - cos(theta)) - u[1] * sin(theta),
                   u[2] * u[1] * (1 - cos(theta)) + u[0] * sin(theta),
                   cos(theta) + u[2]**2 * (1 - cos(theta))]])
    return atoms.dot(r.T)


def random_angles():
    d = np.random.random(size=2)
    phi, theta = 2 * np.pi * d[0], np.arccos(2 * d[1] - 1)
    return theta, phi


def dir_to_angles(d):
    theta = np.arccos(d[2])
    phi = np.arccos(d[0] / np.sin(theta))
    if d[1] / np.sin(theta) < 0:
        phi = 2 * np.pi - phi
    return theta, phi

# sample uniformly in the (half) unit sphere


def random_direction(half=False, quarter=False, ext_range=None):
    d = np.random.random(size=2)
    phi = 2 * np.pi * d[0]
    if ext_range is not None:
        sta, stb = np.cos(ext_range[0] / 180.0 *
                          np.pi), np.cos(ext_range[1] / 180.0 * np.pi)
        theta = np.arccos(d[1] * (stb - sta) + sta)
    elif quarter:
        st = np.sqrt(2.0) / 2.0
        theta = np.arccos(d[1] * (1 - st) + st)  # st~1
    elif half:
        theta = np.arccos(d[1])  # 0~1
    else:
        theta = np.arccos(2 * d[1] - 1)  # -1~1
    return np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta),
                     np.cos(theta)])


def rotated_molecule(name):
    the, phi = random_angles()
    ang = np.pi * np.random.random()
    return rotation(moles[name].atoms, ang, the, phi)


class CreateCluster(object):
    eps = np.finfo(dtype=float).eps * 100

    def __init__(self, n=0, atoms=None, elems=None, surf=None):
        self.n = n
        self.atoms = np.zeros((n, 3), dtype=float) if atoms is None else atoms
        self.elems = np.array(
            [''] * n, dtype='|S20') if elems is None else elems
        self.surf = surf
        self.ext_range = None

    def to_cluster(self):
        if self.surf is None:
            return Cluster(self.n, elems=self.elems, atoms=self.atoms)
        else:
            c = ClusterAtSurface(self.n, surf=self.surf)
            c.elems = self.elems
            c.atoms = self.atoms
            return c

    @staticmethod
    def generator(surf, xn=0, xpool=None, **opts):
        # xn is the current isomer atom number, not including creation atoms
        # creation atoms is determined by len(elems) in opts
        if "grid" in opts:
            assert surf is not None and xn == 0 and xpool is None
            ucellmat = to_cellmat(surf.unit_cell)
            gx, gy = opts["grid"]
            for i in range(0, gx):
                for j in range(0, gy):
                    ffp = 1.0 * i / gx * ucellmat[0] + 1.0 * j / gy * ucellmat[1]
                    ffx, ffy = ffp[0], ffp[1]
                    c = CreateCluster()
                    c.surf = surf
                    c.create_grid(ffx=ffx, ffy=ffy, **opts)
                    yield c.to_cluster()
        else:
            while True:
                xat = xpool.atoms if xpool is not None else None
                xel = xpool.elems if xpool is not None else None
                c = CreateCluster(n=xn, atoms=xat, elems=xel)
                c.surf = surf
                while c.create(**opts) is None:
                    pass
                yield c.to_cluster()

    def update_core(self, xs_core, xs_other):
        k = np.random.randint(0, len(xs_core) - 1)
        if k < len(xs_other):
            xs_core = np.append(xs_core, xs_other[k])
            xs_other = np.delete(xs_other, k)
        return xs_core, xs_other

    # first order atom 2d
    def add_atom_f_2d(self, telem, clelems, cura, mean, sigma, suelems, suatoms,
                      soelems=None, soatoms=None, atst=None, quarter=False):
        assert soelems is None or len(soelems) == 0
        assert suelems is None or len(suelems) == 0
        assert not quarter
        elemr = clelems[:cura]
        atomr = self.atoms[:cura]
        assert atst is None
        atst = self.atoms[:cura].mean(axis=0)
        x = 0.0
        kdir = random_direction(ext_range=[90.0, 90.0])
        for j in range(0, len(elemr)):
            m = telem + '-' + elemr[j]
            l = np.random.normal(mean[m][0], sigma[m][0])
            va = atomr[j] - atst
            b = -2 * np.dot(va, kdir)
            c = np.square(va).sum() - l ** 2
            delta = b**2 - 4 * c + CreateCluster.eps
            if delta > 0.0:
                xx = (-b + np.sqrt(delta)) / 2
                if xx > x:
                    x = xx
        xatom = atst + x * kdir
        return xatom

    # first order atom
    def add_atom_f(self, telem, clelems, cura, mean, sigma, suelems, suatoms,
                   soelems=None, soatoms=None, atst=None, quarter=False, nocluselems=None):
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        if atst is None:
            atst = self.atoms[:cura].mean(axis=0)
        tt = 0
        while True:
            tt += 1
            if tt > 500:
                return None
            x = 0.0
            kdir = random_direction(quarter=quarter)
            qelem = ''
            for j in range(0, len(elemr)):
                m = telem + '-' + elemr[j]
                l = np.random.normal(mean[m][0], sigma[m][0])
                va = atomr[j] - atst
                b = -2 * np.dot(va, kdir)
                c = np.square(va).sum() - l ** 2
                delta = b**2 - 4 * c + CreateCluster.eps
                if delta > 0.0:
                    xx = (-b + np.sqrt(delta)) / 2
                    if xx > x:
                        x = xx
                        qelem = elemr[j]
            xatom = atst + x * kdir
            okay = True
            if soelems is not None:
                for j in range(0, len(soelems)):
                    m = telem + '-' + soelems[j]
                    kl = np.random.normal(mean[m][0], sigma[m][0])
                    kla = np.linalg.norm(soatoms[j] - xatom)
                    if kla < kl:
                        okay = False
                        break
            if nocluselems is not None and okay:
                if (isinstance(nocluselems, list) and qelem in nocluselems) or \
                    (isinstance(nocluselems, dict) and telem in nocluselems and
                    qelem in nocluselems[telem]):
                    okay = False
                else:
                    for j in range(0, len(elemr)):
                        if (isinstance(nocluselems, list) and elemr[j] in nocluselems) or \
                            (isinstance(nocluselems, dict) and telem in nocluselems and
                            elemr[j] in nocluselems[telem]):
                            m = telem + '-' + elemr[j]
                            kl = np.random.normal(mean[m][0], sigma[m][0])
                            kla = np.linalg.norm(atomr[j] - xatom)
                            if kla < kl:
                                okay = False
                                break
            if okay:
                break
        return xatom

    # second order atom 2d
    def add_atom_s_2d(self, telem, clelems, cura, mix_rate, mean, sigma, suelems, suatoms,
                      soelems=None, soatoms=None, nocluselems=None):
        assert soelems is None or len(soelems) == 0
        assert suelems is None or len(suelems) == 0
        elemr = clelems[:cura]
        atomr = self.atoms[:cura]
        tt = 0
        rk = np.array([0, 0, 1])
        while True:
            tt += 1
            if tt > 500:
                return None
            shu = range(len(elemr))
            np.random.shuffle(shu)
            ix, iy = shu[0:2]
            mx = telem + '-' + elemr[ix]
            my = telem + '-' + elemr[iy]
            lx = np.random.normal(mean[mx][0], sigma[mx][0])
            if my != mx:
                ly = np.random.normal(mean[my][0], sigma[my][0])
            else:
                ly = 0
                while ly < lx:
                    ly = np.random.normal(mean[my][1], sigma[my][1])
            surfn = atomr[iy] - atomr[ix]
            lz = np.linalg.norm(surfn)
            if lz > lx + ly or lx > lz + ly or ly > lz + lx:
                continue
            lxx = ((lz**2 + lx**2 - ly**2) / (2 * lz * lx)) * lx
            lr = np.sqrt(lx ** 2 - lxx ** 2)
            cent = surfn * lxx / lz + atomr[ix]
            surfn = surfn / np.linalg.norm(surfn)
            ok = False
            for i in [-1, 1]:
                ra = np.cross(surfn, rk * i)
                ra = cent + ra * lr / np.linalg.norm(ra)
                okr = True
                for k in shu[2:]:
                    mg = telem + '-' + elemr[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(atomr[k] - ra) < lg:
                        okr = False
                        break
                if okr:
                    ok = True
                    break
            if ok:
                break
        return ra

    # second order atom
    def add_atom_s(self, telem, clelems, cura, mix_rate, mean, sigma, suelems, suatoms,
                   soelems=None, soatoms=None):
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        tt = 0
        while True:
            tt += 1
            if tt > 500:
                return None
            imix = False
            if np.random.random() < mix_rate or cura == 1:
                imix = True
            shu = range(len(elemr))
            if imix and len(elemr) > cura:
                ix = np.random.randint(0, cura)
                iy = cura + np.random.randint(0, len(elemr) - cura)
                del shu[ix], shu[iy - 1]
                shu = [ix, iy] + shu
            else:
                np.random.shuffle(shu)
                ix, iy = shu[0:2]
            mx = telem + '-' + elemr[ix]
            my = telem + '-' + elemr[iy]
            lx = np.random.normal(mean[mx][0], sigma[mx][0])
            if my != mx:
                ly = np.random.normal(mean[my][0], sigma[my][0])
            else:
                ly = 0
                while ly < lx:
                    ly = np.random.normal(mean[my][1], sigma[my][1])
            surfn = atomr[iy] - atomr[ix]
            lz = np.linalg.norm(surfn)
            if lz > lx + ly or lx > lz + ly or ly > lz + lx:
                continue
            lxx = ((lz**2 + lx**2 - ly**2) / (2 * lz * lx)) * lx
            lr = np.sqrt(lx ** 2 - lxx ** 2)
            cent = surfn * lxx / lz + atomr[ix]
            surfn = surfn / np.linalg.norm(surfn)
            ok = False
            for _ in range(5):
                dd = random_direction()
                ra = dd - np.dot(dd, surfn) * surfn
                ra = cent + ra * lr / np.linalg.norm(ra)
                okr = True
                for k in shu[2:]:
                    mg = telem + '-' + elemr[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(atomr[k] - ra) < lg:
                        okr = False
                        break
                if soelems is not None:
                    for k in range(0, len(soelems)):
                        mg = telem + '-' + soelems[k]
                        lg = np.random.normal(mean[mg][0], sigma[mg][0])
                        if np.linalg.norm(soatoms[k] - ra) < lg:
                            okr = False
                            break
                if okr:
                    ok = True
                    break
            if ok:
                break
        return ra

    # molecule: removing hydrogen
    def add_molecule_r(self, mname, clelems, cura, mean, sigma, suelems, suatoms,
                       selected, soelems=None, soatoms=None, ext_range=[0.0, 30.0]):
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        m_selected = []
        for ie, el in enumerate(moles[mname].elems):
            if el == 'HX':
                m_selected.append(ie)
        while True:
            matoms = rotated_molecule(mname)
            masel = selected[np.random.randint(0, len(selected))]
            mhsel = m_selected[np.random.randint(0, len(m_selected))]
            mdir = random_direction(ext_range=ext_range)
            dists = np.linalg.norm(
                matoms - matoms[mhsel].reshape(1, 3), axis=1)
            mbsel = None
            for i in range(matoms.shape[0]):
                if i == mhsel:
                    continue
                if mbsel is None or dists[i] < dists[mbsel]:
                    mbsel = i
            mlb = matoms[mhsel] - matoms[mbsel]
            mlb = mlb / np.linalg.norm(mlb)
            rota = np.cross(np.array([0.0, 0.0, 1.0]), mlb)
            rotd = np.arccos(np.dot(np.array([0.0, 0.0, 1.0]), mlb))
            mrot_the, mrot_phi = dir_to_angles(rota / np.linalg.norm(rota))
            mlbx = rotation(mdir, rotd, mrot_the, mrot_phi)
            mbsele = moles[mname].elems[mbsel]
            masele = elemr[masel]
            m = mbsele + '-' + masele
            klx = np.random.normal(mean[m][0], sigma[m][0])
            mlbx = mlbx * klx
            mzb = atomr[masel] - mlbx
            xatoms = matoms + (mzb - matoms[mbsel]).reshape(1, 3)
            okay = True
            for xi in range(len(xatoms)):
                xee = moles[mname].elems[xi]
                if xee == 'HX':
                    xee = 'H'
                if xi == mhsel:
                    continue
                for j in range(0, len(elemr)):
                    if j == masel and xi == mbsel:
                        kl = klx - CreateCluster.eps
                    else:
                        m = xee + '-' + elemr[j]
                        kl = np.random.normal(mean[m][0], sigma[m][0])
                    kla = np.linalg.norm(atomr[j] - xatoms[xi])
                    if kla < kl:
                        okay = False
                        break
                if soelems is not None:
                    for j in range(0, len(soelems)):
                        m = xee + '-' + soelems[j]
                        kl = np.random.normal(mean[m][0], sigma[m][0])
                        kla = np.linalg.norm(soatoms[j] - xatoms[xi])
                        if kla < kl:
                            okay = False
                            break
                if not okay:
                    break
            if okay:
                break
        xxatoms = np.zeros((xatoms.shape[0] - 1, 3))
        xxatoms[0:mhsel] = xatoms[0:mhsel]
        xxatoms[mhsel:] = xatoms[mhsel + 1:]
        return xxatoms, mhsel, masel

    # first order molecule
    def add_molecule_f(self, mname, clelems, cura, mean, sigma, suelems, suatoms,
                       soelems=None, soatoms=None, atst=None, quarter=False, 
                       nomoleelems=None, nocluselems=None):
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        if atst is None:
            atst = self.atoms[:cura].mean(axis=0)
        if isinstance(nocluselems, dict):
            if mname in nocluselems:
                nocluselems = nocluselems[mname]
            else:
                nocluselems = None
        if isinstance(nomoleelems, dict):
            if mname in nomoleelems:
                nomoleelems = nomoleelems[mname]
            else:
                nomoleelems = None
        while True:
            matoms = rotated_molecule(mname)
            # select one atom from molecule
            if nomoleelems is not None:
                while True:
                    masel = np.random.randint(0, len(matoms))
                    mesel = moles[mname].elems[masel]
                    if mesel not in nomoleelems:
                        break
            else:
                masel = np.random.randint(0, len(matoms))
                mesel = moles[mname].elems[masel]
            mdir = random_direction()
            mrot = np.random.random() * 2 * np.pi
            mla = np.linalg.norm(matoms[masel])  # length
            if np.abs(mla) < CreateCluster.eps:
                maa = mvla = mnla = 0.0
            else:
                maa = np.arccos(np.dot(matoms[masel], mdir) / mla)  # angle
                mvla = np.sin(maa) * mla  # length, vert
                mnla = np.cos(maa) * mla
            while True:
                x = 0.0
                kmxa = None
                kdir = random_direction(quarter=quarter)
                kls = []
                kmxae = None
                for j in range(0, len(elemr)):
                    m = mesel + '-' + elemr[j]
                    kl = np.random.normal(mean[m][0], sigma[m][0])
                    kls.append(kl)
                    kla = np.linalg.norm(atomr[j] - atst)  # length
                    if np.abs(kla) < CreateCluster.eps:
                        if kl > mvla:
                            xx = np.sqrt(kl**2 - mvla**2) + mnla
                            if xx > x:
                                x = xx
                                kmxa = random_direction()
                                kmxae = elemr[j]
                    else:
                        kaa = np.arccos(
                            np.dot(atomr[j] - atst, kdir) / kla)  # angle
                        kvla = np.sin(kaa) * kla  # length, vert
                        kc = np.sqrt(mvla**2 + kvla**2 - 2 *
                                    mvla * kvla * np.cos(mrot))
                        if kl > kc:
                            xx = np.sqrt(kl**2 - kc**2) + mnla + np.cos(kaa) * kla
                            if xx > x:
                                x = xx
                                kmxa = atomr[j] - atst - np.cos(kaa) * kla * kdir
                                kmxa = kmxa / np.linalg.norm(kmxa)
                                kmxae = elemr[j]
                if nocluselems is None or kmxae is None:
                    break
                else:
                    koe = False
                    for ke in nocluselems:
                        if kmxae == ke:
                            koe = True
                            break
                    if koe:
                        continue
                    else:
                        break
            if kmxa is None:
                continue
            mrotn = np.cross(mdir, -kdir)
            mrot_ang = np.arccos(np.dot(mdir, -kdir))
            mrot_the, mrot_phi = dir_to_angles(mrotn / np.linalg.norm(mrotn))
            mgatoms = rotation(matoms, mrot_ang, mrot_the, mrot_phi)
            mxb = mgatoms[masel] + np.cos(maa) * mla * kdir
            if np.linalg.norm(mxb) < CreateCluster.eps:
                mxb = random_direction()
            mxb = mxb / np.linalg.norm(mxb)
            mxa = kmxa
            mxrotn = np.cross(mxa, mxb)
            mxrot_ang = np.arccos(np.dot(mxa, mxb))
            mxrot_the, mxrot_phi = dir_to_angles(
                mxrotn / np.linalg.norm(mxrotn))
            xatoms = rotation(mgatoms, -mxrot_ang + mrot, mxrot_the, mxrot_phi)
            xatoms += kdir * x + atst
            okay = True
            for xi in range(len(xatoms)):
                xee = moles[mname].elems[xi]
                if xi == masel:
                    for j, kl in enumerate(kls):
                        kla = np.linalg.norm(atomr[j] - xatoms[xi])
                        if kla < kl - CreateCluster.eps:
                            okay = False
                            break
                else:
                    for j in range(0, len(elemr)):
                        m = xee + '-' + elemr[j]
                        kl = np.random.normal(mean[m][0], sigma[m][0])
                        kla = np.linalg.norm(atomr[j] - xatoms[xi])
                        if kla < kl:
                            okay = False
                            break
                if soelems is not None:
                    for j in range(0, len(soelems)):
                        m = xee + '-' + soelems[j]
                        kl = np.random.normal(mean[m][0], sigma[m][0])
                        kla = np.linalg.norm(soatoms[j] - xatoms[xi])
                        if kla < kl:
                            okay = False
                            break
                if not okay:
                    break
            if okay:
                break
        return xatoms, masel

    # second order molecule
    def add_molecule_s(self, mname, clelems, cura, mean, sigma, suelems, suatoms,
                       soelems=None, soatoms=None, atst=None, quarter=False, nomoleelems=None):
        xatoms, masel = self.add_molecule_f(mname, clelems, cura, mean, sigma, suelems,
                                            suatoms, soelems, soatoms, atst, quarter, nomoleelems)
        if len(xatoms) == 1:
            return xatoms
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        if isinstance(nomoleelems, dict):
            if mname in nomoleelems:
                nomoleelems = nomoleelems[mname]
            else:
                nomoleelems = None
        tt = 0
        while True:
            tt += 1
            if tt > 500:
                return None
            if nomoleelems is not None:
                while True:
                    mxsel = np.random.randint(0, len(xatoms) - 1)
                    if mxsel >= masel:
                        mxsel += 1
                    mexsel = moles[mname].elems[mxsel]
                    if mexsel not in nomoleelems:
                        break
            else:
                mxsel = np.random.randint(0, len(xatoms) - 1)
                if mxsel >= masel:
                    mxsel += 1
                mexsel = moles[mname].elems[mxsel]
            mcsel = np.random.randint(0, len(atomr))
            my = mexsel + '-' + elemr[mcsel]
            lx = np.linalg.norm(xatoms[masel] - xatoms[mxsel])
            ly = np.random.normal(mean[my][0], sigma[my][0])
            surfn = atomr[mcsel] - xatoms[masel]
            lz = np.linalg.norm(surfn)
            if lz > lx + ly or lx > lz + ly or ly > lz + lx:
                continue
            lxx = ((lz**2 + lx**2 - ly**2) / (2 * lz * lx)) * lx
            lr = np.sqrt(lx ** 2 - lxx ** 2)
            cent = surfn * lxx / lz + xatoms[masel]
            surfn = surfn / np.linalg.norm(surfn)
            ok = False
            mdir = xatoms[mxsel] - xatoms[masel]
            mdir = mdir / np.linalg.norm(mdir)
            for _ in range(5):
                dd = random_direction()
                ra = dd - np.dot(dd, surfn) * surfn
                ra = cent + ra * lr / np.linalg.norm(ra)
                ndir = ra - xatoms[masel]
                ndir = ndir / np.linalg.norm(ndir)
                xxatoms = xatoms - xatoms[masel].reshape((1, 3))
                mrotn = np.cross(mdir, ndir)
                mrot_ang = np.arccos(np.dot(mdir, ndir))
                mrot_the, mrot_phi = dir_to_angles(
                    mrotn / np.linalg.norm(mrotn))
                mgatoms = rotation(xxatoms, mrot_ang, mrot_the, mrot_phi)
                mgatoms = mgatoms + xatoms[masel].reshape((1, 3))
                # print ly, np.linalg.norm(mgatoms[mxsel] - atomr[mcsel])
                okr = True
                for u in range(0, len(xatoms)):
                    if u == masel:
                        continue
                    telem = moles[mname].elems[u]
                    for k in range(0, len(atomr)):
                        if u == mxsel and k == mcsel:
                            continue
                        mg = telem + '-' + elemr[k]
                        lg = np.random.normal(mean[mg][0], sigma[mg][0])
                        if np.linalg.norm(atomr[k] - mgatoms[u]) < lg:
                            okr = False
                            break
                    if soelems is not None:
                        for k in range(0, len(soelems)):
                            mg = telem + '-' + soelems[k]
                            lg = np.random.normal(mean[mg][0], sigma[mg][0])
                            if np.linalg.norm(soatoms[k] - mgatoms[u]) < lg:
                                okr = False
                                break
                    if not okr:
                        break
                if okr:
                    ok = True
                    break
            if ok:
                break
        return mgatoms

    # get adjacent atoms
    def get_adjcent(self, telem, clelems, cura, mean, sigma, suelems, suatoms,
                    xatom, soelems=None, soatoms=None):
        elemr = np.array(list(
            clelems[:cura]) + list(suelems) + ([] if soelems is None else list(soelems)))
        atomr = np.array(list(
            self.atoms[:cura]) + list(suatoms) + ([] if soatoms is None else list(soatoms)))
        res = []
        for k in range(len(elemr)):
            mg = telem + '-' + elemr[k]
            lg = mean[mg][0] + 2 * sigma[mg][0]
            if np.linalg.norm(atomr[k] - xatom) < lg:
                res.append(elemr[k])
        return res

    def get_def_hydro_num(self, d=0.5):
        max_z = self.surf.atoms.max(axis=0)[2]
        return len([x for x in self.surf.atoms if x[2] >= max_z - d])
    
    def hydro_cluster(self, num, elist, refhats, mean, sigma):
        helem = 'H'
        hatoms = []
        print ('num = ', num)
        os_idx = [ix for ix, x in enumerate(self.elems) if x in elist]
        np.random.shuffle(os_idx)
        if len(os_idx) < num:
            return hatoms
        for kk in os_idx:
            xatom = None
            for _ in range(20):
                mg = helem + '-' + self.elems[kk]
                lg = np.random.normal(mean[mg][0], sigma[mg][0])
                kdir = random_direction()
                xatom = kdir * lg + self.atoms[kk]
                okr = True
                for k in range(0, len(self.surf.atoms)):
                    mg = helem + '-' + self.surf.elems[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(self.surf.atoms[k] - xatom) < lg:
                        okr = False
                        break
                for k in range(0, len(self.atoms)):
                    if not okr:
                        break
                    if k == kk:
                        continue
                    mg = helem + '-' + self.elems[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(self.atoms[k] - xatom) < lg:
                        okr = False
                        break
                for k in range(0, len(refhats + hatoms)):
                    if not okr:
                        break
                    mg = helem + '-' + helem
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm((refhats + hatoms)[k] - xatom) < lg:
                        okr = False
                        break
                if not okr:
                    xatom = None
                    continue
            if xatom is not None:
                hatoms.append(xatom)
        print ('hydroed = ', os_idx, len(hatoms))
        return hatoms[:num]

    def hydro_surf(self, mean, sigma, d=0.5, ang=30.0):
        max_z = self.surf.atoms.max(axis=0)[2]
        hs_idx = [ix for ix, x in enumerate(
            self.surf.atoms) if x[2] >= max_z - d]
        helem = 'H'
        hatoms = []
        for kk in hs_idx:
            xatom = None
            for _ in range(10):
                mg = helem + '-' + self.surf.elems[kk]
                lg = np.random.normal(mean[mg][0], sigma[mg][0])
                kdir = random_direction(ext_range=[0.0, ang])
                xatom = kdir * lg + self.surf.atoms[kk]
                okr = True
                for k in range(0, len(self.surf.atoms)):
                    if k == kk:
                        continue
                    mg = helem + '-' + self.surf.elems[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(self.surf.atoms[k] - xatom) < lg:
                        okr = False
                        break
                for k in range(0, len(self.atoms)):
                    if not okr:
                        break
                    mg = helem + '-' + self.elems[k]
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(self.atoms[k] - xatom) < lg:
                        okr = False
                        break
                for k in range(0, len(hatoms)):
                    if not okr:
                        break
                    mg = helem + '-' + helem
                    lg = np.random.normal(mean[mg][0], sigma[mg][0])
                    if np.linalg.norm(hatoms[k] - xatom) < lg:
                        okr = False
                        break
                if not okr:
                    xatom = None
                    continue
            if xatom is not None:
                hatoms.append(xatom)
        return hatoms
    
    def create_grid(self, ffx, ffy, mean, sigma, elems, **opts):
        assert self.n == 0 and len(elems) == 1
        assert self.surf is not None
        # enlarge surface to avoid drop off
        surfabak = np.array(self.surf.atoms)
        surfebak = np.array(self.surf.elems)
        atomsu = np.zeros((surfabak.shape[0] * 9, 3))
        elemsu = []
        cellmat = to_cellmat(self.surf.cell)
        for ix, x in enumerate([-1, 0, 1]):
            for iy, y in enumerate([-1, 0, 1]):
                ik = ix * 3 + iy
                surfadd = np.array(
                    [x * cellmat[0] + y * cellmat[1]])
                atomsu[ik * surfabak.shape[0]:(ik + 1) * surfabak.shape[0]] = \
                    surfabak + surfadd
                elemsu += list(surfebak)
        elemsu = np.array(elemsu)
        hz = self.surf.atoms[:, 2].min()
        xelem = elems[0]
        for iia, iie in zip(atomsu, elemsu):
            lxy = np.sqrt((iia[0] - ffx)**2 + (iia[1] - ffy)**2)
            mg = xelem + '-' + iie
            lg = np.random.normal(mean[mg][0], sigma[mg][0])
            if lg < lxy:
                continue
            lz = np.sqrt(lg ** 2 - lxy ** 2) + iia[2]
            if lz > hz:
                hz = lz
        if np.abs(hz - self.surf.atoms[:, 2].min()) < CreateCluster.eps:
            hz = self.surf.atoms[:, 2].max()
        self.atoms = np.array([[ffx, ffy, hz]], dtype=float)
        self.elems = np.array([xelem], dtype='|S20')
        self.n = 1
    
    # hydro: add hydrogen atoms on the surface
    # nomoleelems: not bind molecule onto cluster via elems in this list
    # rmhydro: remove howmany H atoms because cluster occupied some space
    def create(self, elems, mean, sigma, order=2, morder=1, ncore=4,
               mix_rate=0.5, core_ratio=0.75, ext_range=[0.0, 30.0], twod=0.0,
               hydro=False, rmhydro=3, nomoleelems=None, nocluselems=None, lowpos=-999.0,
               loworderelems=None):
        # generation group-wise order
        if len(elems) == 0:
            return self
        if isinstance(elems[0], np.ndarray):
            gelems = []
            indx = []
            for ex in elems:
                indxt = range(len(ex))
                np.random.shuffle(indxt)
                indx += [i + len(gelems) for i in indxt]
                gelems += list(ex)
            elems = np.array(gelems)
            elemsx = np.array(elems)
            elemsx[indx] = np.array(elemsx[:])
        else:
            elemsx = np.array(elems)
            indx = range(len(elems))
            np.random.shuffle(indx)
            elemsx[indx] = np.array(elemsx[:])
        if loworderelems is None:
            loworderelems = []
        if nocluselems is not None and len(nocluselems) > 1 and nocluselems[0].startswith('-'):
            nox = {}
            lst = None
            for nx in nocluselems:
                if nx.startswith('-'):
                    lst = nx[1:]
                    nox[lst] = []
                else:
                    nox[lst].append(nx)
            nocluselems = nox
        if nomoleelems is not None and len(nomoleelems) > 1 and nomoleelems[0].startswith('-'):
            nox = {}
            lst = None
            for nx in nomoleelems:
                if nx.startswith('-'):
                    lst = nx[1:]
                    nox[lst] = []
                else:
                    nox[lst].append(nx)
            nomoleelems = nox
        # determine n
        # elemx: mol name, disordered, self.elems: at name, ordered
        sx_elemsbak = np.array(self.elems)
        sx_nbak = self.n
        sx_atomsbak = np.array(self.atoms)
        self.elems = []
        for ex in elems:
            if ex in moles:
                self.elems += [x if x !=
                               "HX" else "H" for x in list(moles[ex].elems)]
            else:
                self.elems.append(ex)
        self.n = len(self.elems)
        self.elems = np.array(self.elems)
        self.atoms = np.zeros((self.n, 3))
        # elemr: atom name, disordered
        elemr = []
        for ex in elemsx:
            if ex in moles:
                elemr += [x if x !=
                          "HX" else "H" for x in list(moles[ex].elems)]
            else:
                elemr.append(ex)
        # indxr: index for the disordered, indexr ordered by the ordered
        # for created atoms
        indxrp = []
        ipid = 0
        for ex in elemsx:
            if ex in moles:
                indxrp.append(range(ipid, ipid + moles[ex].n))
                ipid += moles[ex].n
            else:
                indxrp.append([ipid])
                ipid += 1
        indxr = []
        for i in indx:
            indxr += indxrp[i]
        if sx_nbak != 0:
            self.n += sx_nbak
            self.elems = np.array(list(sx_elemsbak) +
                                  list(self.elems), dtype=sx_elemsbak.dtype)
            self.atoms = np.array(list(sx_atomsbak) + list(self.atoms))
            elemr = list(sx_elemsbak) + elemr
            indxr = range(len(sx_elemsbak)) + \
                [i + len(sx_elemsbak) for i in indxr]
        # surface part
        if self.surf is not None:
            assert twod == 0.0
            # enlarge surface to avoid drop off
            surfabak = np.array(self.surf.atoms)
            surfebak = np.array(self.surf.elems)
            atomsu = np.zeros((surfabak.shape[0] * 18, 3))
            elemsu = []
            cellmat = to_cellmat(self.surf.cell)
            ucellmat = to_cellmat(self.surf.unit_cell)
            cellmat[2, 2] = self.surf.cellz
            for ix, x in enumerate([-1, 0, 1]):
                for iy, y in enumerate([-1, 0, 1]):
                    for iz, z in enumerate([-1, 0]):
                        ik = ix * 6 + iy * 2 + iz
                        surfadd = np.array(
                            [x * cellmat[0] + y * cellmat[1] + z * cellmat[2]])
                        atomsu[ik * surfabak.shape[0]:(ik + 1) * surfabak.shape[0]] = \
                            surfabak + surfadd
                        elemsu += list(surfebak)
            elemsu = np.array(elemsu)
            # start point for core selection
            startp = np.array(list(self.surf.atoms.mean(axis=0)[0:2]) +
                              [self.surf.atoms.max(axis=0)[2]])
            min_z = self.surf.atoms.min(axis=0)[2]
            max_z = self.surf.atoms.max(axis=0)[2]
            # lowpos (abs value) indicates below max_z + lowpos position cluster atoms will be rejected
            # lowpos may be set to minus value
            if max_z + lowpos > min_z:
                min_z = max_z + lowpos
            # strat point randomly selected from unit cell
            startp += (np.random.random() - 0.5) * ucellmat[0]
            startp += (np.random.random() - 0.5) * ucellmat[1]
            lsurf = np.linalg.norm(atomsu - startp.reshape(1, 3), axis=1)
            # indices of core and other
            xs_core = np.argsort(lsurf)[0:ncore]
            xs_other = np.argsort(lsurf)[ncore:]
            st_min = atomsu[xs_core].min(axis=0)[0:2]
            st_max = atomsu[xs_core].max(axis=0)[0:2]
            st_z = atomsu[xs_core].mean(axis=0)[2]
            # start point for first atom
            st = np.array([np.random.uniform(st_min[0], st_max[0]),
                           np.random.uniform(st_min[1], st_max[1]), st_z])
            soelems, soatoms = elemsu[xs_other], atomsu[xs_other]
            suelems, suatoms = elemsu[xs_core], atomsu[xs_core]
            if hydro:
                def_hydro_num = self.get_def_hydro_num()
                print ("def_hydro = ", def_hydro_num)
        else:
            xs_core = np.array([], dtype=int)
            xs_other = np.array([], dtype=int)
            soelems, soatoms = None, None
            suelems, suatoms = np.array(
                [], dtype='|S20'), np.array([], dtype=float)
            if hydro:
                def_hydro_num = 0
        cura = sx_nbak  # current atom
        cure = 0  # current element or molecule
        if self.surf is None and cura != 0:
            self.atoms[0:cura] -= self.atoms[0:cura].mean(axis=0).reshape(1, 3)
        satoms_bak = np.array(self.atoms)
        selected = None
        xorder = order
        ltwod = False
        if twod != 0.0:
            if np.random.random() < twod:
                ltwod = True
        while cura != self.n:
            # xorder == 3 means use half first order, half second order
            if cura < self.n / 2 and xorder == 3:
                order = 1
            elif cura >= self.n / 2 and xorder == 3:
                order = 2
            while self.surf is not None and len(xs_core) < core_ratio * cura:
                xs_core, xs_other = self.update_core(xs_core, xs_other)
                suelems, suatoms = elemsu[xs_core], atomsu[xs_core]
                soelems, soatoms = elemsu[xs_other], atomsu[xs_other]
            if cura == 0 and self.surf is None:
                if elemsx[cure] in moles:
                    xatoms = rotated_molecule(elemsx[cure])
                    self.atoms[cura:cura + xatoms.shape[0], :] = xatoms
                    cura += xatoms.shape[0]
                else:
                    self.atoms[cura] = 0.0
                    cura += 1
            elif (((cura == 0 and self.surf is not None) or
                   (cura == 1 and self.surf is None) or
                    order == 1 or elemsx[cure] in loworderelems
                    ) and elemsx[cure] not in moles) or \
                    (morder == 1 and elemsx[cure] in moles):
                # molecule case
                if elemsx[cure] in moles:
                    if cura == 0:
                        atst, qua = st, True
                    else:
                        atst, qua = None, False
                    for _ in range(0, 100):
                        xatoms, _ = self.add_molecule_f(elemsx[cure], elemr, cura, mean,
                                                        sigma, suelems, suatoms, soelems=soelems, soatoms=soatoms,
                                                        atst=atst, quarter=qua, nomoleelems=nomoleelems,
                                                        nocluselems=nocluselems)
                        # special order-2, ensure no low O
                        if xorder == -2 or xorder == -3 and 'O' in moles[elemsx[cure]].elems:
                            if xatoms[:, 2].min() < max_z:
                                print ('low O detected!')
                                continue
                        if self.surf is None or xatoms[:, 2].min() > min_z:
                            break
                        else:
                            xatoms = None
                    self.atoms[cura:cura + xatoms.shape[0], :] = xatoms
                    cura += xatoms.shape[0]
                    if self.surf is None:
                        self.atoms[0:cura] -= self.atoms[0:
                                                         cura].mean(axis=0).reshape(1, 3)
                else:
                    if cura == 0:
                        atst, qua = st, True
                    else:
                        atst, qua = None, False
                    for _ in range(0, 100):
                        if not ltwod:
                            xatom = self.add_atom_f(elemr[cura], elemr, cura, mean, sigma,
                                                    suelems, suatoms, soelems=soelems, soatoms=soatoms, atst=atst,
                                                    quarter=qua, nocluselems=nocluselems)
                        else:
                            xatom = self.add_atom_f_2d(elemr[cura], elemr, cura, mean, sigma,
                                                       suelems, suatoms, soelems=soelems, soatoms=soatoms, atst=atst,
                                                       quarter=qua)
                        if self.surf is None or xatom[2] > min_z:
                            break
                        else:
                            xatom = None
                    if xatom is None:
                        print ("restarted first!")
                        cura = sx_nbak  # current atom
                        cure = 0  # current element or molecule
                        self.atoms = np.array(satoms_bak)
                        if self.surf is not None:
                            st = np.array([np.random.uniform(st_min[0], st_max[0]),
                                           np.random.uniform(st_min[1], st_max[1]), st_z])
                        continue
                    self.atoms[cura] = xatom
                    cura += 1
                    if self.surf is None:
                        self.atoms[0:cura] -= self.atoms[cura - 1] / cura
            elif morder == -1 and elemsx[cure] in moles:
                # requires removing one hydrogen
                if selected is None:
                    selected = range(cura)
                assert len(selected) != 0
                for _ in range(0, 100):
                    xxatoms, mhsel, masel = self.add_molecule_r(elemsx[cure], elemr, cura,
                                                                mean, sigma, suelems, suatoms, selected, soelems=soelems,
                                                                soatoms=soatoms, ext_range=ext_range)
                    if self.surf is None or xxatoms[:, 2].min() > min_z:
                        break
                    else:
                        xxatoms = None
                self.atoms[cura:cura + xxatoms.shape[0], :] = xxatoms
                self.n -= 1
                assert elemr[cura + mhsel] == "H"
                del elemr[cura + mhsel]
                indxrtmp = []
                for idxr in indxr:
                    if idxr == cura + mhsel:
                        continue
                    elif idxr > cura + mhsel:
                        indxrtmp.append(idxr - 1)
                    else:
                        indxrtmp.append(idxr)
                indxr = indxrtmp
                cura += xxatoms.shape[0]
                if self.surf is None:
                    self.atoms[0:cura] -= self.atoms[0:
                                                     cura].mean(axis=0).reshape(1, 3)
                selected.remove(masel)
            else:  # second order
                # molecule case
                if elemsx[cure] in moles:
                    if morder == 2:
                        if cura == 0:
                            atst, qua = st, True
                        else:
                            atst, qua = None, False
                        for _ in range(0, 100):
                            xatoms = self.add_molecule_s(elemsx[cure], elemr, cura, mean, sigma,
                                                         suelems, suatoms, soelems=soelems, soatoms=soatoms, atst=atst,
                                                         quarter=qua, nomoleelems=nomoleelems)
                            if self.surf is None or (xatoms is not None and xatoms[:, 2].min() > min_z):
                                break
                            else:
                                xatoms = None
                        if xatoms is None:
                            print ("restarted second!")
                            cura = sx_nbak  # current atom
                            cure = 0  # current element or molecule
                            self.atoms = np.array(satoms_bak)
                            if self.surf is not None:
                                st = np.array([np.random.uniform(st_min[0], st_max[0]),
                                               np.random.uniform(st_min[1], st_max[1]), st_z])
                            continue
                        self.atoms[cura:cura + xatoms.shape[0], :] = xatoms
                        cura += xatoms.shape[0]
                        if self.surf is None:
                            self.atoms[0:cura] -= self.atoms[0:
                                                             cura].mean(axis=0).reshape(1, 3)
                    else:
                        raise NotImplementedError(
                            'molecule order > 2 not implemented!')
                else:
                    assert not np.isnan(self.atoms[0, 0])
                    foo, flo = False, False
                    for _ in range(0, 100):
                        if not ltwod:
                            xatom = self.add_atom_s(elemr[cura], elemr, cura, mix_rate,
                                                    mean, sigma, suelems, suatoms, soelems=soelems, soatoms=soatoms)
                        else:
                            xatom = self.add_atom_s_2d(elemr[cura], elemr, cura, mix_rate,
                                                       mean, sigma, suelems, suatoms, soelems=soelems, soatoms=soatoms)
                        # special order-2, ensure no - O-O
                        if (xorder == -2 or xorder == -3) and elemr[cura] == 'O':
                            lkk = self.get_adjcent(
                                elemr[cura], elemr, cura, mean, sigma, suelems, suatoms, xatom, soelems, soatoms)
                            if 'O' in lkk:
                                if not foo:
                                    print ('O-O detected!')
                                    foo = True
                                xatom = None
                            elif xatom[2] < max_z:
                                if not flo:
                                    print ('low O detected!')
                                    flo = True
                                xatom = None
                        if xatom is not None and (self.surf is None or xatom[2] > min_z):
                            break
                        else:
                            xatom = None
                    if xatom is None:
                        print ("restarted second!")
                        cura = sx_nbak  # current atom
                        cure = 0  # current element or molecule
                        self.atoms = np.array(satoms_bak)
                        if self.surf is not None:
                            st = np.array([np.random.uniform(st_min[0], st_max[0]),
                                           np.random.uniform(st_min[1], st_max[1]), st_z])
                        continue
                    self.atoms[cura] = xatom
                    cura += 1
                    if self.surf is None:
                        self.atoms[0:cura] -= self.atoms[cura - 1] / cura
                    assert not np.isnan(self.atoms[0, 0])
            cure += 1
        self.atoms = np.array(self.atoms[indxr])
        if selected is not None:
            self.elems = np.array(elemr)[indxr]
        if hydro:
            hatoms = self.hydro_surf(mean, sigma)
            rem = def_hydro_num - len(hatoms)
            if xorder == -3:
                def_hydro_num -= rmhydro
            else:
                hatoms += self.hydro_cluster(rem, ['O'], hatoms, mean, sigma)
            if len(hatoms) != def_hydro_num:
                print ('hydro failed! added %d' % len(hatoms))
                self.atoms = sx_atomsbak
                self.elems = sx_elemsbak
                self.n = sx_nbak
                return None
            self.atoms = np.array(list(self.atoms) + hatoms)
            self.elems = np.array(list(self.elems) + ['H'] * len(hatoms))
            self.n += def_hydro_num
        return self


if __name__ == "__main__":
    from cluster.base import elem_char, update_stat
    from surface.base import read_surface
    surf = read_surface('../PGOPT/surfaces/mgo.xyz')
    surf.unit_cell = np.array([4.257, 4.257, 4.257])
    surf.space_group = "-F 4 2 3"
    xmean = {}
    xsigma = {}
    np.random.seed(100)
    k, _ = elem_char('Pt5(CH3OH)')
    kxx = [elem_char('Pt5')[0], elem_char('(CH3OH)')[0]]
    update_stat(np.array(list(k) + ["Mg", "O"]), xmean, xsigma, 0.0004)
    print (xmean, xsigma)
    app = False
    for i in range(100):
        print (i)
        c = CreateCluster()
        c.surf = surf
        c.create(kxx, xmean, xsigma, 2).to_cluster(
        ).write_xyz('test.xyz', append=app)
        app = True
