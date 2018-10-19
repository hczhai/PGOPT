
import copy, random
import numpy as np
from surface.base import ClusterAtSurface
from surface.create_cluster import random_direction, CreateCluster
from surface.surf_comp import DisjointSet
from surface.surf_symm import to_cellmat
from cluster.coval import Covalent, AtomicWeight, LightElements
from cluster.base import moles, Cluster, elem_char

def sort_elems(clu, elems):
    clu = copy.deepcopy(clu)
    dd = {}
    dii = {}
    assert len(elems) == clu.n
    for i in range(0, clu.n):
        if clu.elems[i] not in dd:
            dd[clu.elems[i]] = [np.array(clu.atoms[i])]
            dii[clu.elems[i]] = 0
        else:
            dd[clu.elems[i]].append(np.array(clu.atoms[i]))
    for i in range(0, clu.n):
        clu.elems[i] = elems[i]
        clu.atoms[i] = dd[elems[i]][dii[elems[i]]]
        dii[elems[i]] += 1
    return clu

def extract_molecule(clx, xelems, xatoms, main_elems, aux_elems, final_factor):
    ic = -1
    rgx = np.array(range(0, clx.n))
    np.random.shuffle(rgx)
    for i in rgx:
        if clx.elems[i] in main_elems:
            ic = i
            break
    if ic == -1:
        return
    ccelems = []
    ccatoms = []
    for i in range(0, clx.n):
        if i == ic:
            xelems.append(clx.elems[i])
            xatoms.append(clx.atoms[i])
        elif clx.elems[i] in aux_elems:
            rr = Covalent.x[clx.elems[i]] + Covalent.x[clx.elems[ic]]
            r = np.linalg.norm(clx.atoms[i] - clx.atoms[ic])
            if r < rr * final_factor:
                xelems.append(clx.elems[i])
                xatoms.append(clx.atoms[i])
            else:
                ccelems.append(clx.elems[i])
                ccatoms.append(clx.atoms[i])
        else:
            ccelems.append(clx.elems[i])
            ccatoms.append(clx.atoms[i])
    clx.n = len(ccelems)
    clx.elems = np.array(ccelems, dtype='|S20')
    clx.atoms = np.array(ccatoms, dtype=float)

def swap_site(clu, xelems, xatoms, mc):
    final_factor = mc.FinalFactor
    clx = copy.deepcopy(clu)
    xelems = list(xelems)
    xatoms = list(xatoms)
    main_elems = []
    aux_elems = []
    if mc.KeepCH3:
        if 'C' not in xelems:
            main_elems = ['C']
            aux_elems = ['H']
    elif mc.KeepCO:
        if 'C' not in xelems:
            main_elems = ['C']
            aux_elems = ['O']
    elif xelems == []:
        main_elems = LightElements
    extract_molecule(clx, xelems, xatoms, main_elems, aux_elems, final_factor)
    for ex in mc.SwapSiteMakeSpace:
        if ex in LightElements:
            if clx.elems.count(ex) == 0:
                continue
            while True:
                ri = random.randint(0, clx.n - 1)
                if clx.elems[ri] != ex:
                    continue
                xelems.append(ex)
                xatoms.append(clx.atoms[ri])
                clx.n -= 1
                clx.elems = np.array(list(clx.elems[0:ri]) + list(clx.elems[ri + 1:]))
                clx.atoms = np.array(list(clx.atoms[0:ri]) + list(clx.atoms[ri + 1:]))
                break
        else:
            exx = elem_char(ex)[0]
            main_elems = [exx[0]]
            aux_elems = exx[1:]
            extract_molecule(clx, xelems, xatoms, main_elems, aux_elems, final_factor)
    return clx, xelems, xatoms

def pull_back(clu, xelems, xatoms, mc, creopts):
    final_factor = mc.FinalFactor
    clx = copy.deepcopy(clu)
    xelems = list(xelems)
    xatoms = list(xatoms)
    use_temp = []
    assert not (mc.KeepCO and mc.KeepCH3)
    if mc.KeepCO:
        while xelems.count('C') != xelems.count('O'):
            while True:
                target = 'C' if xelems.count('C') < xelems.count('O') else 'O'
                ri = random.randint(0, clx.n - 1)
                if clx.elems[ri] != target:
                    continue
                xelems += [target]
                clx.n -= 1
                clx.elems = np.array(list(clx.elems[0:ri]) + list(clx.elems[ri + 1:]))
                clx.atoms = np.array(list(clx.atoms[0:ri]) + list(clx.atoms[ri + 1:]))
                break
        gelems = []
        for i in range(0, len(xelems)):
            if xelems[i] == 'C':
                gelems.append('CO')
            elif xelems[i] == 'O':
                continue
            else:
                gelems.append(xelems[i])
        xelems = gelems
    if mc.KeepCH3:
        ic = -1
        if 'C' not in xelems:
            for i in range(0, clx.n):
                if clx.elems[i] == 'C':
                    ic = i
                    break
            kkidx = []
            for i in range(0, clx.n):
                if clx.elems[i] != 'H':
                    continue
                rr = Covalent.x[clx.elems[i]] + Covalent.x[clx.elems[ic]]
                r = np.linalg.norm(clx.atoms[i] - clx.atoms[ic])
                if r < rr * final_factor:
                    kkidx.append(i)
            if len(kkidx) >= 3:
                ic = -1
            else:
                kkidx.append(ic)
                gelems = []
                gatoms = []
                for i in range(0, clx.n):
                    if i in kkidx:
                        if clx.elems[i] == 'C':
                            ic = len(xelems)
                        xelems.append(clx.elems[i])
                        xatoms.append(clx.atoms[i])
                    else:
                        gelems.append(clx.elems[i])
                        gatoms.append(clx.atoms[i])
                clx.n = len(gelems)
                clx.elems = np.array(gelems, dtype=clx.elems.dtype)
                clx.atoms = np.array(gatoms, dtype=clx.atoms.dtype)
        else:
            for i in range(0, len(xelems)):
                if xelems[i] == 'C':
                    ic = i
        if ic != -1:
            hofc_idx = []
            for i in range(0, len(xelems)):
                if xelems[i] != 'H':
                    continue
                rr = Covalent.x[xelems[i]] + Covalent.x[xelems[ic]]
                r = np.linalg.norm(xatoms[i] - xatoms[ic])
                if r < rr * final_factor:
                    hofc_idx.append(i)
            if len(hofc_idx) > 3:
                hofc_idx = hofc_idx[:3]
            if len(hofc_idx) == 3:
                moles["CH3TEMP"] = moles["CH3"].copy()
                moles["CH3TEMP"].mm.elems = np.array(['C', 'H', 'H', 'H'])
                moles["CH3TEMP"].mm.atoms = np.array([xatoms[ic]] + [xatoms[i] for i in hofc_idx])
                moles["CH3TEMP"].mm.atoms -= moles["CH3TEMP"].mm.atoms.mean(axis=0)
                use_temp.append('CH3')
                gelems = ["CH3"]
                for i in range(0, len(xelems)):
                    if i not in hofc_idx and i != ic:
                        gelems.append(xelems[i])
                xelems = gelems
            else:
                while len([i for i in xelems if i == 'H']) < 3:
                    while True:
                        ri = random.randint(0, clx.n - 1)
                        if clx.elems[ri] != 'H':
                            continue
                        xelems += ['H']
                        clx.n -= 1
                        clx.elems = np.array(list(clx.elems[0:ri]) + list(clx.elems[ri + 1:]))
                        clx.atoms = np.array(list(clx.atoms[0:ri]) + list(clx.atoms[ri + 1:]))
                        break
                gelems = ["CH3"]
                ih = len(hofc_idx)
                for i in range(0, len(xelems)):
                    if ih != 3 and i not in hofc_idx and xelems[i] == 'H':
                        ih += 1
                    elif i not in hofc_idx and i != ic:
                        gelems.append(xelems[i])
                xelems = gelems
        cck = ConnectivityCheck(clx)
        clx, cxelems, cxatoms = cck.separate()
        xelems = list(xelems) + list(cxelems)
        cxatoms = list(xatoms) + list(cxatoms)
    if mc.WriteSteps:
        clx.label = "after pull-back pre-treatment"
        clx.write_xyz("steps.xyz")
    cc = CreateCluster(clx.n, clx.atoms, clx.elems, clx.surf)
    creopts["elems"] = list(xelems)
    bak = []
    for ut in use_temp:
        bak.append(moles[ut].copy())
        moles[ut] = moles[ut + "TEMP"]
    cc.create(**creopts)
    for ut, bk in zip(use_temp, bak):
        moles[ut] = bk
    cg = cc.to_cluster()
    if mc.WriteSteps:
        cg.label = "after pull back"
        cg.write_xyz("steps.xyz")
    return cg

def make_movements(clu, mc, creopts):
    step = mc.CurrentStepLength
    lstep = mc.CurrentStepLength / mc.StepLength * mc.LightStepLength
    ffactor = mc.ShortDistanceFactor
    ref_sep = mc.RefuseSeparate
    if mc.WriteSteps:
        clu.label = "initial"
        clu.write_xyz("steps.xyz", append=False)
    if mc.LightShell and any(clu.elems[i] in LightElements for i in range(clu.n)):
        # detect separate atoms
        cck = ConnectivityCheck(clu)
        clx, cxelems, cxatoms = cck.separate()
        if mc.WriteSteps:
            clx.label = "after remove separate"
            clx.write_xyz("steps.xyz")
        # swap site
        if mc.SwapSite and random.random() < mc.SwapSiteRate:
            clx, cxelems, cxatoms = swap_site(clx, cxelems, cxatoms, mc)
            if mc.WriteSteps:
                clx.label = "after swap site"
                clx.write_xyz("steps.xyz")
        # pull back
        clu = pull_back(clx, cxelems, cxatoms, mc, creopts)
        if mc.WriteSteps:
            clu.label = "after initial pull back"
            clu.write_xyz("steps.xyz")
        # make movements of metal atoms
        main_atoms = []
        main_elems = []
        main_indices = []
        metal_atoms = []
        metal_elems = []
        for i in range(0, clu.n):
            if clu.elems[i] in LightElements:
                main_atoms.append(clu.atoms[i])
                main_elems.append(clu.elems[i])
            else:
                metal_atoms.append(clu.atoms[i])
                metal_elems.append(clu.elems[i])
        r = copy.deepcopy(clu)
        r.n = len(metal_elems)
        r.atoms = np.array(metal_atoms, r.atoms.dtype)
        r.elems = np.array(metal_elems, r.elems.dtype)
        if mc.WriteSteps:
            r.label = "metal only"
            r.write_xyz("steps.xyz")
        rr = make_movements_simple(r, step)
        if ref_sep:
            while True:
                cck = ConnectivityCheck(rr)
                if cck.check():
                    break
                rr = make_movements_simple(r, step)
        if mc.WriteSteps:
            rr.label = "after metal move"
            rr.write_xyz("steps.xyz")
        # pre-optimization of metal atoms
        preopt = PreOptimization(rr, ffactor)
        preopt.steplength = 0.2
        preopt.opt()
        rr = preopt.traj[-1]
        if mc.WriteSteps:
            rr.label = "after metal opt"
            rr.write_xyz("steps.xyz")
        # assign indices
        nmain = len(main_elems)
        for i in range(nmain):
            ix = -1
            dx = None
            for j in range(0, r.n):
                d = np.linalg.norm(metal_atoms[j] - main_atoms[i])
                if dx is None or d < dx:
                    dx = d
                    ix = j
            main_indices.append(ix)
        # co-move
        new_main_atoms = copy.deepcopy(main_atoms)
        new_main_elems = copy.deepcopy(main_elems)
        for i in range(nmain):
            # kdir = random_direction()
            # new_main_atoms[i] += lstep * Covalent.x[new_main_elems[i]] * 2 * kdir
            new_main_atoms[i] += rr.atoms[main_indices[i]] - metal_atoms[main_indices[i]]
        # main move
        rm = copy.deepcopy(clu)
        aux_elems = []
        aux_atoms = []
        aux_indices = []
        if mc.SolidMove != []:
            nx_main_atoms = []
            nx_main_elems = []
            main_list = []
            aux_list = []
            aux_dict = {}
            for sm in mc.SolidMove:
                smx = elem_char(sm)[0]
                main_list.append(smx[0])
                aux_list.extend(smx[1:])
                for sx in smx[1:]:
                    if sx not in aux_dict:
                        aux_dict[sx] = [smx[0]]
                    else:
                        aux_dict[sx].append(smx[0])
            assert len(set.intersection(set(aux_list), set(main_list))) == 0
            for i in range(0, len(new_main_elems)):
                if new_main_elems[i] not in aux_list:
                    nx_main_elems.append(new_main_elems[i])
                    nx_main_atoms.append(new_main_atoms[i])
                else:
                    aux_elems.append(new_main_elems[i])
                    aux_atoms.append(new_main_atoms[i])
            ad_nx_main_elems = []
            ad_nx_main_atoms = []
            for i in range(len(aux_elems)):
                ix = -1
                dx = None
                for j in range(0, len(nx_main_elems)):
                    d = np.linalg.norm(nx_main_atoms[j] - aux_atoms[i])
                    if dx is None or d < dx:
                        dx = d
                        ix = j
                if ix != -1 and nx_main_elems[ix] not in aux_dict[aux_elems[i]]:
                    ix = -1
                aux_indices.append(ix)
                if ix == -1:
                    ad_nx_main_elems.append(aux_elems[i])
                    ad_nx_main_atoms.append(aux_atoms[i])
            # if there are H and CH3, do not put H in aux
            aux_elems = [el for idx, el in zip(aux_indices, aux_elems) if idx != -1]
            aux_atoms = [am for idx, am in zip(aux_indices, aux_atoms) if idx != -1]
            new_main_elems = nx_main_elems + ad_nx_main_elems
            new_main_atoms = nx_main_atoms + ad_nx_main_atoms
        rm.n = len(new_main_elems)
        rm.atoms = np.array(new_main_atoms, rm.atoms.dtype)
        rm.elems = np.array(new_main_elems, rm.elems.dtype)
        if mc.WriteSteps:
            xrm = copy.deepcopy(rm)
            xrm.n = rm.n + rr.n
            xrm.atoms = np.array(list(xrm.atoms) + list(rr.atoms), xrm.atoms.dtype)
            xrm.elems = np.array(list(xrm.elems) + list(rr.elems), xrm.elems.dtype)
            xrm.label = "after put back main elems"
            xrm.write_xyz("steps.xyz")
        rmm = make_movements_simple(rm, lstep)
        # if ref_sep:
        #     while True:
        #         cck = ConnectivityCheck(rmm)
        #         if cck.check():
        #             break
        #         rmm = make_movements_simple(rm, lstep)
        # pre-optimization of main atoms, fix metal atoms
        if mc.WriteSteps:
            xrm = copy.deepcopy(rmm)
            xrm.n = rmm.n + rr.n
            xrm.atoms = np.array(list(xrm.atoms) + list(rr.atoms), xrm.atoms.dtype)
            xrm.elems = np.array(list(xrm.elems) + list(rr.elems), xrm.elems.dtype)
            xrm.label = "after main elems move"
            xrm.write_xyz("steps.xyz")
        if mc.SolidMove != []:
            # aux co-move
            new_aux_atoms = copy.deepcopy(aux_atoms)
            new_aux_elems = copy.deepcopy(aux_elems)
            for i in range(len(aux_elems)):
                new_aux_atoms[i] += rmm.atoms[aux_indices[i]] - rm.atoms[aux_indices[i]]
            rmm.n = rmm.n + len(aux_elems)
            rmm.atoms = np.array(list(rmm.atoms) + list(new_aux_atoms), rmm.atoms.dtype)
            rmm.elems = np.array(list(rmm.elems) + list(new_aux_elems), rmm.elems.dtype)
            if mc.WriteSteps:
                xrm = copy.deepcopy(rmm)
                xrm.n = rmm.n + rr.n
                xrm.atoms = np.array(list(xrm.atoms) + list(rr.atoms), xrm.atoms.dtype)
                xrm.elems = np.array(list(xrm.elems) + list(rr.elems), xrm.elems.dtype)
                xrm.label = "after aux elems co-move"
                xrm.write_xyz("steps.xyz")
        rf = copy.deepcopy(clu)
        rf.n = rmm.n
        rf.atoms = np.array(rmm.atoms, rf.atoms.dtype)
        rf.elems = np.array(rmm.elems, rf.elems.dtype)
        rf.surf.n = rf.surf.n + rr.n
        rf.surf.atoms = np.array(list(rf.surf.atoms) + list(rr.atoms), rf.atoms.dtype)
        rf.surf.elems = np.array(list(rf.surf.elems) + list(rr.elems), rf.elems.dtype)
        preopt = PreOptimization(rf, ffactor)
        preopt.steplength = 0.2
        preopt.opt()
        rf = preopt.traj[-1]
        # compile final structure
        ff = copy.deepcopy(clu)
        ima, ime = 0, 0
        for i in range(0, clu.n):
            if clu.elems[i] in LightElements:
                ff.atoms[i] = rf.atoms[ima]
                ff.elems[i] = rf.elems[ima]
                ima += 1
            else:
                ff.atoms[i] = rf.surf.atoms[ime + clu.surf.n]
                ff.elems[i] = rf.surf.elems[ime + clu.surf.n]
                ime += 1
        # final pull back
        cck = ConnectivityCheck(ff)
        clx, cxelems, cxatoms = cck.separate()
        ff = pull_back(clx, cxelems, cxatoms, mc, creopts)
        # sort elements
        ff = sort_elems(ff, clu.elems)
        if mc.WriteSteps:
            ff.label = "final"
            ff.write_xyz("steps.xyz")
        return ff
    else:
        r = make_movements_simple(clu, step)
        if ref_sep:
            while True:
                cck = ConnectivityCheck(r)
                if cck.check():
                    break
                r = make_movements_simple(clu, step)
        preopt = PreOptimization(r, ffactor)
        preopt.steplength = 0.2
        preopt.opt()
        if mc.WriteSteps:
            preopt.traj[-1].label = "final"
            preopt.traj[-1].write_xyz("steps.xyz")
        return preopt.traj[-1]

def make_movements_simple(clu, step):
    r = copy.deepcopy(clu)
    for i in range(0, r.n):
        kdir = random_direction()
        r.atoms[i] += step * Covalent.x[r.elems[i]] * 2 * kdir
    return r

class ConnectivityCheck(DisjointSet):
    def __init__(self, clu):
        self.clu = clu
        self.finalfactor = 1.2
        super(ConnectivityCheck, self).__init__(clu.n)
        self.sur = []
        if isinstance(clu, ClusterAtSurface):
            cellmat = to_cellmat(clu.surf.cell)
            for ix in [-1, 0, 1]:
                for iy in [-1, 0, 1]:
                    self.sur.append(ix * cellmat[0] + iy * cellmat[1])
        else:
            self.sur.append(np.array([0.0, 0.0, 0.0]))

    def separate(self):
        clx = copy.deepcopy(self.clu)
        idx_conn = []
        idx_sep = []
        self.clear()
        strx = self.clu
        for s in self.sur:
            for i in range(0, strx.n):
                for j in range(0, i):
                    rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                    r = np.linalg.norm(strx.atoms[i] - strx.atoms[j] - s)
                    if r < rr * self.finalfactor:
                        self.union(i, j)
        if strx.n > 1:
            idx_conn.append(0)
            for i in range(1, strx.n):
                if self.find(i) != self.find(0):
                    idx_sep.append(i)
                else:
                    idx_conn.append(i)
        if not isinstance(strx, ClusterAtSurface) or strx.n == 0:
            pass
        else:
            sok = False
            for s in self.sur + [np.zeros(3)]:
                for i in range(0, strx.n):
                    if self.find(i) == self.find(0):
                        for j in range(0, strx.surf.n):
                            rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.surf.elems[j]]
                            r = np.linalg.norm(strx.atoms[i] - strx.surf.atoms[j] - s)
                            if r < rr * self.finalfactor:
                                sok = True
                                break
            if not sok:
                idx_sep += idx_conn
                idx_conn = []
        xatoms = clx.atoms[np.array(idx_sep, dtype=int)]
        xelems = clx.elems[np.array(idx_sep, dtype=int)]
        clx.n = len(idx_conn)
        clx.atoms = clx.atoms[np.array(idx_conn, dtype=int)]
        clx.elems = clx.elems[np.array(idx_conn, dtype=int)]
        return (clx, xelems, xatoms)

    def check(self):
        return self.separate()[0].n == self.clu.n

class PreOptimization(object):
    def __init__(self, clu, ffactor):
        self.clu = clu
        self.finalfactor = ffactor
        self.steplength = 0.05
        self.traj = []
        self.maxiter = 1000
        self.sur = []
        if isinstance(clu, ClusterAtSurface):
            cellmat = to_cellmat(clu.surf.cell)
            for ix in [-1, 0, 1]:
                for iy in [-1, 0, 1]:
                    if ix == 0 and iy == 0:
                        continue
                    self.sur.append(ix * cellmat[0] + iy * cellmat[1])

    def check(self, strx):
        for i in range(0, strx.n):
            for j in range(0, i):
                rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                r = np.linalg.norm(strx.atoms[i] - strx.atoms[j])
                if r < rr * self.finalfactor:
                    return False
        if isinstance(strx, ClusterAtSurface):
            for s in self.sur + [np.zeros(3)]:
                for i in range(0, strx.n):
                    for j in range(0, strx.surf.n):
                        rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.surf.elems[j]]
                        r = np.linalg.norm(strx.atoms[i] - strx.surf.atoms[j] - s)
                        if r < rr * self.finalfactor:
                            return False
        for s in self.sur:
            for i in range(0, strx.n):
                for j in range(0, i):
                    rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                    r = np.linalg.norm(strx.atoms[i] - strx.atoms[j] - s)
                    if r < rr * self.finalfactor:
                        return False
        return True

    def opt(self):
        strx = self.clu
        self.traj.append(strx)
        it = 0
        while not self.check(strx) and it < self.maxiter:
            fce = self.force(strx)
            fcel = np.linalg.norm(fce, axis=1)
            fcelm = fcel.max()
            fce = fce / fcelm * self.steplength
            strt = copy.deepcopy(strx)
            strt.atoms += fce
            strx = strt
            self.traj.append(strx)
            it += 1
        return it < self.maxiter

    def energy(self, strx):
        ener = 0.0
        for i in range(0, strx.n):
            for j in range(0, i):
                rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                rr = rr / 2**(1.0/6.0)
                r = np.linalg.norm(strx.atoms[i] - strx.atoms[j])
                ener += (rr / r)**12 - (rr / r)**6
        if isinstance(strx, ClusterAtSurface):
            for s in self.sur + [np.zeros(3)]:
                for i in range(0, strx.n):
                    for j in range(0, strx.surf.n):
                        rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.surf.elems[j]]
                        rr = rr / 2**(1.0/6.0)
                        r = np.linalg.norm(strx.atoms[i] - strx.surf.atoms[j] - s)
                        ener += (rr / r)**12 - (rr / r)**6
        for s in self.sur:
            for i in range(0, strx.n):
                for j in range(0, strx.n):
                    rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                    rr = rr / 2**(1.0/6.0)
                    r = np.linalg.norm(strx.atoms[i] - strx.atoms[j] - s)
                    ener += (rr / r)**12 - (rr / r)**6
        return ener

    def force(self, strx):
        fce = np.zeros((strx.n, 3), dtype=np.float64)
        for i in range(0, strx.n):
            for j in range(0, i):
                rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                rr = rr / 2**(1.0/6.0)
                r = np.linalg.norm(strx.atoms[i] - strx.atoms[j])
                ed = -(12 * (rr / r)**13 - 6 * (rr / r)**7) / rr
                fce[i] += (strx.atoms[i] - strx.atoms[j]) / r * ed
                fce[j] += -(strx.atoms[i] - strx.atoms[j]) / r * ed
        if isinstance(strx, ClusterAtSurface):
            for s in self.sur + [np.zeros(3)]:
                for i in range(0, strx.n):
                    for j in range(0, strx.surf.n):
                        rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.surf.elems[j]]
                        rr = rr / 2**(1.0/6.0)
                        r = np.linalg.norm(strx.atoms[i] - strx.surf.atoms[j] - s)
                        ed = -(12 * (rr / r)**13 - 6 * (rr / r)**7) / rr
                        fce[i] += (strx.atoms[i] - strx.surf.atoms[j]) / r * ed
        for s in self.sur:
            for i in range(0, strx.n):
                for j in range(0, strx.n):
                    rr = Covalent.x[strx.elems[i]] + Covalent.x[strx.elems[j]]
                    rr = rr / 2**(1.0/6.0)
                    r = np.linalg.norm(strx.atoms[i] - strx.atoms[j] - s)
                    ed = -(12 * (rr / r)**13 - 6 * (rr / r)**7) / rr
                    fce[i] += (strx.atoms[i] - strx.atoms[j]) / r * ed
        return -fce
