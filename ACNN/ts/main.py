
from __future__ import print_function
import numpy as np
import copy
import theano
import theano.tensor as T
from formod.at_comp import at_comp, at_comp_list
from formod.km_comp import km_comp
from formod.lbfgs import LBFGS
from cluster.base import Cluster, read_zip_clusters, elem_num
from cluster.coval import Covalent
from surface.base import read_surface, ClusterAtSurface
from surface.surf_comp import surface_align, surface_compare, apply_trans
from surface.surf_comp import surface_match, surface_adjust, cas_adjust
from surface.surf_symm import to_direct, to_cartesian, to_cellmat
from parallel.monte_carlo import PreOptimization
from parallel.main import RecordWritter
from utils.io import new_file_name


def rotation_g(atoms, rot_ang, dir_the, dir_phi):
    theta = rot_ang
    the = dir_the
    phi = dir_phi
    u = [T.sin(the) * T.cos(phi), T.sin(the) * T.sin(phi), T.cos(the)]
    va = T.stack([T.cos(theta) + u[0]**2 * (1 - T.cos(theta)),
                  u[0] * u[1] * (1 - T.cos(theta)) - u[2] * T.sin(theta),
                  u[0] * u[2] * (1 - T.cos(theta)) + u[1] * T.sin(theta)], axis=0)
    vb = T.stack([u[1] * u[0] * (1 - T.cos(theta)) + u[2] * T.sin(theta),
                  T.cos(theta) + u[1]**2 * (1 - T.cos(theta)),
                  u[1] * u[2] * (1 - T.cos(theta)) - u[0] * T.sin(theta)], axis=0)
    vc = T.stack([u[2] * u[0] * (1 - T.cos(theta)) - u[1] * T.sin(theta),
                  u[2] * u[1] * (1 - T.cos(theta)) + u[0] * T.sin(theta),
                  T.cos(theta) + u[2]**2 * (1 - T.cos(theta))], axis=0)
    r = T.stack([va, vb, vc], axis=0)
    return T.dot(atoms, r.T)

# only works for no periodic condition


def check(e, a, b, nmax, factor):
    for k in range(0, nmax + 1):
        x = a + (b - a) * k / nmax
        for i in range(0, len(x)):
            for j in range(0, i):
                rr = Covalent.x[e[i]] + Covalent.x[e[j]]
                r = np.linalg.norm(x[i] - x[j])
                if r < rr * factor:
                    return False
    return True


class Connection(object):
    RotFunc = []

    def __init__(self, clus, iprp, endp, idx, ref_structs=None):
        self.clus = clus
        self.endp = endp
        self.idx = idx
        self.nimages = iprp["nimages"]
        self.max_diff = iprp["max_diff"]
        self.acceptable_short_length_factor = 0.4
        self.adjusted_short_length_factor = 0.7
        self.max_trials = 20
        self.max_trials_surf = iprp.get("max_trials_surf", 8)
        # self.max_trials_surf = 2
        self.trail_tolerance = 0.5
        self.path_min_diff = 0.1
        self.ref_structs = ref_structs

    def init_func(self):
        print ("preparing rotation function...")
        x = T.dvector()
        atomsa = T.dmatrix()
        atomsb = T.dmatrix()
        xx = rotation_g(atomsa, x[0], x[1], x[2]) + x[3:6].reshape((1, 3))
        xxc = (xx - atomsb).norm(2, axis=1).mean()
        ff = theano.function([atomsa, atomsb, x], xxc)
        fd = theano.function([atomsa, atomsb, x], T.grad(xxc, x))
        fr = theano.function([atomsa, x], xx)
        self.RotFunc.extend([ff, fd, fr])

    def build(self):
        assert len(self.clus) != 0
        if isinstance(self.clus[0], ClusterAtSurface):
            if self.ref_structs is not None:
                return self.build_surf_ref()
            else:
                return self.interp_surf()
        elif len(self.clus) == 2:
            return self.build_gas_phase()
        else:
            return self.interp_gas_phase()

    def build_gas_phase(self):
        if len(self.RotFunc) == 0:
            self.init_func()
        ca = self.clus[0]
        cb = self.clus[1]
        eles, ne = elem_num(ca.elems)
        eles2, ne2 = elem_num(list(ca.elems) * 2)
        nmax = self.nimages + 1

        # initial measurement
        atomsa, atomsb = np.array(ca.atoms), np.array(cb.atoms)
        dd, _ = at_comp(atomsa, atomsb, ne, eles, self.max_diff + 0.1)
        dcl = np.array(self.max_trials)
        uph = self.trail_tolerance

        # build trial set
        dcf = None
        while dcl == self.max_trials and uph > 0.05:
            dcx = at_comp_list(atomsa, atomsb, ne, eles, dd + uph, dcl)
            dcx = dcx[:dcl - 1]
            if dcf is None or dcl == self.max_trials:
                dcf = np.array(dcx)
            uph /= np.sqrt(2)

        # solve rotation non-linear problem
        ntrial = len(dcf)
        nshort = 0
        nhigh = 0
        nsymm = 0
        mdiff = self.max_diff
        finals = []
        task = LBFGS(6)
        task.log_file = 0
        task.max_iter = 500
        for idx in dcf:
            caa, cab = np.array(atomsa), np.array(atomsb)
            caa = caa[idx - 1]
            xi = [0.0] * 6
            task.p.eval = lambda x: self.RotFunc[0](caa, cab, x)
            task.p.evald = lambda x: self.RotFunc[1](caa, cab, x)
            task.start(xi)
            task.opt()
            caf = self.RotFunc[2](caa, task.x)
            hdiff = self.RotFunc[0](caa, cab, task.x)
            if hdiff > self.max_diff:
                nhigh += 1
                continue
            if not check(ca.elems, caf, cab, nmax, self.acceptable_short_length_factor):
                nshort += 1
                continue
            finals.append([caf, cab, hdiff])
            if hdiff < mdiff:
                mdiff = hdiff
        finals.sort(key=lambda x: x[2])
        if len(finals) > 1:
            ufinals = [finals[0]]
            ufinalsrf = [np.concatenate(finals[0][:2])]
            for ff in finals[1:]:
                crg = np.concatenate(ff[:2])
                okay = True
                for crf in ufinalsrf:
                    crd, _ = at_comp(crf, crg, ne2, eles2,
                                     self.path_min_diff + 0.02)
                    if crd < self.path_min_diff:
                        okay = False
                        break
                if not okay:
                    nsymm += 1
                else:
                    ufinals.append(ff)
                    ufinalsrf.append(np.concatenate(ff[:2]))
            finals = ufinals

        # construct final images
        finalimgs = []
        for caa, cab, _ in finals:
            imgs = []
            for i in range(0, nmax + 1):
                cc = Cluster(ca.n, ca.elems)
                cc.atoms = caa + (cab - caa) * i / nmax
                cc.mag = None
                if i == 0:
                    cc.mag, cc.energy = ca.mag, ca.energy
                elif i == nmax:
                    cc.mag, cc.energy = cb.mag, cb.energy
                if cc.mag is None:
                    cc.mag = 0.0
                cc.label = "CONN-%d-%d:%%d.%d:%.2f" % (self.endp[0], self.endp[1],
                                                       i, cc.mag)
                if i != 0 and i != nmax:
                    pcc = PreOptimization(
                        cc, self.adjusted_short_length_factor)
                    pcc.opt()
                    cc = pcc.traj[-1]
                imgs.append(cc)
            finalimgs.append(imgs)
        return finalimgs, [self.endp[0], self.endp[1], len(finalimgs), ntrial, nshort, nhigh,
                           nsymm, mdiff]

    def interp_gas_phase(self):
        pass

    def adjust_surf(self, x, y, b):
        x.atoms = b[3]
        x.surf.atoms = b[4]
        _, sb = surface_match(x.surf, y.surf, 1.0)
        sbinv = np.zeros(len(sb), dtype=int)
        sbinv[sb - 1] = np.array(range(0, len(sb)))
        x.atoms = x.atoms[b[1] - 1]
        x.surf.atoms = x.surf.atoms[sbinv]
        cas_adjust(x, y)
        surface_adjust(x.surf, y.surf)

    def build_surf_ref(self):
        cax = self.clus
        surf = cax[0].surf
        # align atoms
        eles, ne = elem_num(cax[0].elems)
        cellmat = to_cellmat(surf.cell)
        for i in [0, -1]:
            _, b = surface_compare(cax[i], self.ref_structs[i], ne, eles, 0.5)
            self.adjust_surf(cax[i], self.ref_structs[i], b)
        mdiff, _ = km_comp(cax[0].atoms, cax[-1].atoms, ne, eles, cellmat)
        lscax = [cax]
        nmax = self.nimages + 1
        assert(len(self.ref_structs) == nmax + 1)
        uimgs = []
        mishort = 0
        for cax in lscax:
            imgs = []
            for i in range(0, nmax + 1):
                if i == 0:
                    cc = cax[0]
                elif i == nmax:
                    cc = cax[1]
                else:
                    cc = self.ref_structs[i]
                cc = copy.deepcopy(cc)
                cc.surf.forces = None
                if cc.mag is None:
                    cc.mag = 0.0
                cc.label = "CONN-%d-%d:%%d.%d:%.2f" % (self.endp[0], self.endp[1],
                                                       i, cc.mag)
                imgs.append(cc)
            uimgs.append(imgs)
        return uimgs, [self.endp[0], self.endp[1], len(uimgs), 1, mishort, 0, 0, mdiff]

    def interp_surf(self, ibest=None):
        if ibest is None:
            ibest = self.max_trials_surf
        cax = self.clus
        surf = cax[0].surf
        # align atoms
        eles, ne = elem_num(cax[0].elems)
        cellmat = to_cellmat(surf.cell)
        lscax = []
        mdiff = None
        mdmax = None
        misinglelong = 0
        if len(cax) == 2:
            bs = surface_compare(cax[0], cax[1], ne,
                                 eles, self.max_diff + 0.1, best=ibest)
            for md, b in bs:
                if len(lscax) != 0 and (md > self.max_diff or md > mdiff * 1.5):
                    break
                cxx = copy.deepcopy(cax)
                self.adjust_surf(cxx[0], cxx[1], b)
                mdx = np.linalg.norm(cxx[0].atoms - cxx[1].atoms, axis=1).max()
                if len(lscax) != 0 and (mdx > self.max_diff * 4):
                    misinglelong += 1
                    break
                lscax.append(cxx)
                if mdiff is None or md < mdiff:
                    mdiff = md
                if mdmax is None or mdx < mdmax:
                    mdmax = mdx
        else:
            for cai in range(len(cax) - 1, 0, -1):
                _, b = surface_compare(
                    cax[cai - 1], cax[cai], ne, eles, self.max_diff + 0.1)
                self.adjust_surf(cax[cai - 1], cax[cai], b)
            mdiff, _ = km_comp(cax[0].atoms, cax[-1].atoms, ne, eles, cellmat)
            lscax.append(cax)

        # construct final images
        nmax = self.nimages + 1
        dc = 1.0 * nmax / (len(cax) - 1)
        for cax in lscax:
            for cai in range(0, len(cax)):
                cax[cai].idx = cai * dc

        uimgs = []
        mishort = 0
        for cax in lscax:
            imgs = []
            ishort = False
            for i in range(0, nmax + 1):
                cc = ClusterAtSurface(cax[0].n, copy.deepcopy(cax[0].surf))
                cc.elems = cax[0].elems
                for cai in range(1, len(cax)):
                    if cax[cai].idx >= i - 1E-8:
                        break
                cx, cy = cax[cai - 1], cax[cai]
                rx = (cax[cai].idx - i) / dc
                cc.atoms = cx.atoms * rx + cy.atoms * (1 - rx)
                cc.surf.atoms = cx.surf.atoms * rx + cy.surf.atoms * (1 - rx)
                cc.mag = None
                if np.abs(rx - 1) < 1E-8:
                    cc.mag, cc.energy = cx.mag, cx.energy
                    cc.surf.forces = cx.surf.forces
                elif np.abs(rx) < 1E-8:
                    cc.mag, cc.energy = cy.mag, cy.energy
                    cc.surf.forces = cy.surf.forces
                else:
                    cc.surf.forces = None
                if cc.mag is None:
                    cc.mag = 0.0
                cc.label = "CONN-%d-%d:%%d.%d:%.2f" % (self.endp[0], self.endp[1],
                                                       i, cc.mag)
                if i != 0 and i != nmax:
                    pcc = PreOptimization(
                        cc, self.adjusted_short_length_factor)
                    pcc.maxiter = 12
                    lpp = pcc.opt()
                    if not lpp:
                        ishort = True
                        break
                    cc = pcc.traj[-1]
                imgs.append(cc)
            if not ishort:
                uimgs.append(imgs)
            else:
                mishort += 1
        return uimgs, [self.endp[0], self.endp[1], len(uimgs), ibest, mishort,
                       misinglelong, 0, mdiff]


class ConnectRun(object):

    def __init__(self, parent, ip):
        self.ip = ip
        self.p = parent
        self.iprp = self.ip["reaction-path"]
        self.surf = None
        self.init_surface()
        print (self.iprp["input_file"])
        self.clus = read_zip_clusters(self.iprp["input_file"])
        if self.surf is not None:
            self.clus = [ClusterAtSurface.from_cluster(
                sc, self.surf) for sc in self.clus]
            for sc in self.clus:
                surface_align(sc)
        self.nimages = self.iprp["nimages"]
        self.nisomers = self.iprp.get("nisomers", len(self.clus))
        self.clus = self.clus[:self.nisomers]
        self.max_diff = self.iprp["max_diff"]
        self.ref_neblocal = self.init_ref(self.iprp.get("path-ref-local", []))
        self.ref_nebmax = self.init_ref(self.iprp.get("path-ref-max", []))
        if len(self.ref_neblocal) != 0 or len(self.ref_nebmax) != 0:
            from meta.analysis import read_nebpaths
            strs = self.ref_neblocal + self.ref_nebmax
            self.ref_nebpaths, self.ref_nebmin = read_nebpaths(
                strs, len(self.ref_neblocal))
            # for neb path search
            self.ref_path_dict = {}
            for np in self.ref_nebpaths.values():
                # if np.focus and np.isdirect:
                if np.focus:
                    ft = "%d.%d" % tuple(np.from_to)
                    self.ref_path_dict[ft] = np
        else:
            self.ref_nebpaths, self.ref_nebmin = {}, {}

    def init_ref(self, x):
        aclus = []
        for xx in x:
            print ("transforming ref nebpaths ...")
            clus = read_zip_clusters(xx)
            if self.surf is not None:
                # here we did the surface replacement
                clus = [ClusterAtSurface.from_cluster(
                    sc, self.surf) for sc in clus]
            for c in clus:
                labels = c.label.split(":")
                if "." not in labels[1]:
                    c.tid = int(labels[1])
                    c.tidx = [int(labels[1]), 0]
                else:
                    c.tid = int(labels[1].split(".")[0])
                    c.tidx = [int(p) for p in labels[1].split(".")]
                c.stidx = ".".join([str(p) for p in c.tidx])
            # backup results
            if not "trans.xyz." in xx:
                xy = xx.replace(".xyz.", ".trans.xyz.")
                writter = RecordWritter()
                writter.write_clusters(xy, clus, False)
            aclus.extend(clus)
        return aclus

    def init_surface(self):
        if "creation-surface" in self.ip:
            ipcs = self.ip["creation-surface"]
            self.surf = read_surface(ipcs["surface"])
            self.surf.space_group = ipcs["space_group"]
            self.surf.space_group_ref = np.array(ipcs["space_group_ref"])
            self.surf.unit_cell = ipcs["unit_cell"]
            self.surf.fix = ipcs.get("fix", "")

    def run(self):
        # for surface case, need to solve the km_comp with cell periodic problem
        # interp needs to deal with crossing bounday problem
        conn = []
        clus_match = []
        # use ref paths to build
        if len(self.ref_nebmin) != 0:
            print ("## number of ref-isomers = %d" % len(self.ref_nebmin))
            eles, ne = elem_num(self.ref_nebmin.values()[0].elems)
            for i in range(len(self.clus)):
                kmm = -1
                for k, v in self.ref_nebmin.items():
                    aa = surface_compare(self.clus[i], v, ne, eles, 5)[0]
                    if kmm == -1 or aa < kmm:
                        kmm = aa
                        ki = (k, v)
                if kmm <= 0.15:
                    clus_match.append([ki, kmm])
                    print ("matched [ %4d -> %4s ] %.4f *" % (i, ki[0], kmm))
                else:
                    clus_match.append([("", None), -1])
                    print ("failed  [ %4d -> %4s ] %.4f" % (i, ki[0], kmm))
        else:
            for i in range(len(self.clus)):
                clus_match.append([("", None), -1])
        for i in range(len(self.clus)):
            eles, ne = elem_num(self.clus[i].elems)
            for j in range(i):
                if clus_match[i][1] != -1 and clus_match[j][1] != -1 and clus_match[i][0][0] != clus_match[j][0][0]:
                    xi, xj = clus_match[i][0][0], clus_match[j][0][0]
                    oft = "%s.%s" % (xi, xj)
                    rft = "%s.%s" % (xj, xi)
                    if oft in self.ref_path_dict or rft in self.ref_path_dict:
                        if oft in self.ref_path_dict:
                            rstr = self.ref_path_dict[oft].structs
                        elif rft in self.ref_path_dict:
                            rstr = self.ref_path_dict[rft].structs[::-1]
                        print (
                            "[%4d -> %4d ] projected to [%4s -> %4s ] *" % (i, j, xi, xj))
                        conn.append(Connection([self.clus[i], self.clus[j]], self.iprp,
                                               (i, j), len(conn), ref_structs=rstr))
                else:
                    # no ref paths
                    if self.surf is None:
                        d, _ = at_comp(self.clus[i].atoms, self.clus[j].atoms, ne, eles,
                                    self.max_diff + 0.1)
                    else:
                        d, _ = surface_compare(self.clus[i], self.clus[j], ne, eles,
                                            self.max_diff + 0.1)
                    if d < self.max_diff:
                        print ("[%4d -> %4d ] diff = %5.2f *" %
                            (i, j, d))
                        conn.append(Connection([self.clus[i], self.clus[j]], self.iprp,
                                            (i, j), len(conn)))
                    else:
                        print ("[%4d -> %4d ] diff = %5.2f" %
                            (i, j, d))
        if len(self.ref_nebmin) != 0:
            for cc in conn:
                cc.max_trials_surf = 4
        print ("## number of isomers = %d" % len(self.clus))
        print ("## number of pairs = %d" % len(conn))
        series = []
        title = "[from ->  to  ] accepted ntrial nshort nhighd  nsymm min-diff"
        holder = len(title)
        print (('%%%ds' % (holder)) % ('-' * holder, ))
        print (title)
        print (('%%%ds' % (holder)) % ('-' * holder, ))
        for c in conn:
            ser, prt = c.build()
            series.extend(ser)
            print ("[%4d -> %4d ] %8d %6d %6d %6d %6d %8.2f" % tuple(prt))
        print (('%%%ds' % (holder)) % ('-' * holder, ))
        print ("## number of final connections = %d" % len(series))
        conn_fn = new_file_name(self.p.conn_images_name)
        print ("write images: " + conn_fn)
        writter = RecordWritter()
        apd = False
        for it, s in enumerate(series):
            for c in s:
                c.label = c.label % it
            writter.write_clusters(conn_fn, s, apd)
            apd = True
