
from __future__ import print_function
import os, time, random, numpy as np
from utils.base import reproducible, print_dict
from utils.io import read_json, write_json
from utils.io import load_data, dump_data, new_file_name, load_name
from parallel.main import ParallelRun
from ts.main import ConnectRun
from cluster.base import read_zip_clusters

class Driver(object):

    def __init__(self, ip, iprint=True):
        self.ip = ip
        self.iprint = iprint
        self.init_filenames()
        self.init_randomseed()

    def init_filenames(self):
        if "output_dir" not in self.ip:
            raise RuntimeError('No output dir given!')
        if not os.path.exists(self.ip["output_dir"]):
            os.mkdir(self.ip["output_dir"])
        self.fit_network_name = self.ip["output_dir"] + "/fit_network.dill"
        self.fit_test_name = self.ip["output_dir"] + "/fit_test.txt"
        self.fit_pertest_name = self.ip["output_dir"] + "/fit_pertest.txt"
        self.fit_error_name = self.ip["output_dir"] + "/fit_error.txt"
        self.opt_list_name = self.ip["output_dir"] + "/opt_list.txt"
        self.opt_structs_name = self.ip["output_dir"] + "/opt_structs.xyz"
        self.fil_list_name = self.ip["output_dir"] + "/fil_list.txt"
        self.fil_corr_name = self.ip["output_dir"] + "/fil_corr.txt"
        self.fil_structs_name = self.ip["output_dir"] + "/fil_structs.xyz"
        self.paraf_list_name = self.ip["output_dir"] + "/par_flist.txt"
        self.paraf_corr_name = self.ip["output_dir"] + "/par_fcorr.txt"
        self.paraf_structs_name = self.ip["output_dir"] + "/par_filter.xyz"
        self.summary_name = self.ip["output_dir"] + "/summary.json"
        self.para_run_name = self.ip["output_dir"] + "/par_run.json"
        self.para_data_name = self.ip["output_dir"] + "/par_data.dill"
        self.para_log_name = self.ip["output_dir"] + "/par_log.txt"
        self.para_shlog_name = self.ip["output_dir"] + "/par_shlog.txt"
        self.para_structs_name = self.ip["output_dir"] + "/par_structs.xyz"
        self.para_local_name = self.ip["output_dir"] + "/par_local.xyz"
        self.para_x_name = self.ip["output_dir"] + "/par_%s.%s"
        self.para_freqs_name = self.ip["output_dir"] + "/par_freqs.txt"
        self.para_disp_name = self.ip["output_dir"] + "/par_disp.txt"
        self.stat_dist_name = self.ip["output_dir"] + "/stat_dist.json"
        self.ali_structs_name = self.ip["output_dir"] + "/ali_structs.xyz"
        self.enl_structs_name = self.ip["output_dir"] + "/enl_structs.xyz"
        self.report_name = self.ip["output_dir"] + "/report.pdf"
        self.report_structs_name = self.ip["output_dir"] + "/report_structs_%s.xyz"
        self.report_props_name = self.ip["output_dir"] + "/report_props_%s.json"
        self.surfgen_json_name = self.ip["output_dir"] + "/surf.json"
        self.surfgen_xyz_name = self.ip["output_dir"] + "/surf.xyz"
        self.surfgen_orig_name = self.ip["output_dir"] + "/orig.xyz"
        self.conn_images_name = self.ip["output_dir"] + "/conn_images.xyz"

    def init_randomseed(self):
        if "random_seed" in self.ip:
            ips = self.ip["random_seed"]
            if ips == 0:
                ips = int(time.time() * 99991) % (2**32 - 1)
                self.ip["random_seed"] = ips
                if self.iprint:
                    print ('random seed = %d' % (ips, ))
            reproducible(seed=ips)

    def run_all(self):
        if "molecules" in self.ip:
            self.load_molecules()
        for task in self.ip['tasks']:
            if task == 'fit':
                self.fit()
            elif task == 'opt':
                self.opt()
            elif task == 'filter':
                self.filter()
            elif task == 'stat':
                self.stat()
            elif task == 'align':
                self.align()
            elif task == 'report':
                self.report()
            elif task == 'draw':
                self.draw()
            elif task == 'enlarge':
                self.enlarge()
            elif task == 'create':
                self.filter(create=True)
            elif task == 'global':
                ParallelRun(self, self.ip).run()
            elif task == 'connect':
                ConnectRun(self, self.ip).run()
            elif task == 'surfgen':
                self.surfgen()
            else:
                raise RuntimeError('Unknown task: ' + task)
        if 'write_summary' in self.ip and self.ip['write_summary']:
            print ('write summary ...')
            write_json(self.ip, self.summary_name)

    def load_molecules(self):
        from cluster.base import moles, MoleculePool, read_clusters
        for key, v in self.ip["molecules"].items():
            mm = read_clusters(v, iprint=False)
            keyx = key[:1].upper() + key[1:]
            if len(mm) == 1:
                moles[keyx] = MoleculePool(mm[0])
            else:
                moles[keyx] = MoleculePool(mm)

    def surfgen(self):
        from surface.base import read_surface, Surface, ClusterAtSurface
        from surface.surf_comp import surface_align
        from surface.surf_symm import read_contcar, num_list, rev_num_list, to_dcell
        from surface.surf_symm import determine_unit_cell, determine_space_group
        from surface.surf_symm import to_cartesian, to_direct, to_cellmat
        from surface.spacegroupdata import SymOpsHall
        from cluster.base import understand_name, elem_char
        from parallel.main import RecordWritter
        import copy
        ipcr = self.ip.get("creation", {"name": None})
        ipf = self.ip["surfgen"]["input_file"]
        cel = self.ip["surfgen"].get("cell", None)
        nospg = self.ip["surfgen"].get("no_space_group", False)
        fix = self.ip["surfgen"].get("fix", None)
        fixz = self.ip["surfgen"].get("fixz", None)
        gapz = self.ip["surfgen"].get("gapz", None)
        removez = self.ip["surfgen"].get("removez", None)
        supercell = self.ip["surfgen"].get("supercell", None)
        alldet = self.ip["surfgen"].get("alldet", False)
        gdx = self.ip["surfgen"].get("dx", 0.0)
        gdy = self.ip["surfgen"].get("dy", 0.0)
        surf = read_surface(ipf)
        clus = None
        ffix = []
        if surf is None:
            clus = read_zip_clusters(ipf)[0]
            if clus is not None:
                assert clus.surfnum != 0 and cel is not None
                surf = Surface(clus.surfnum)
                surf.cell = np.array(cel)
                clus = ClusterAtSurface.from_cluster(clus, surf)
                surf = clus.surf
        if surf is None:
            trclu, trcel, trfix = read_contcar(ipf)
            if ipcr["name"] is not None:
                _, elems, _ = understand_name(ipcr["name"])
            if ipcr["name"] is not None and len(elems) > 1 and isinstance(elems[0], \
                np.ndarray) and len(elems[0]) == 1:
                dc, _ = elem_char(elems[0][0])
                g = []
                for i in range(trclu.n):
                    g.append([i, trclu.atoms[i], trclu.elems[i], i in trfix,
                              trclu.elems[i] in dc and not i in trfix])
                gc = sorted([i for i in g if i[4]], key=lambda x: -x[1][2])
                gcf = []
                sdc = list(set(dc))
                for s in sdc:
                    gcl = [i[0] for i in gc if i[2] == s]
                    ld = len([d for d in dc if d == s])
                    assert len(gcl) >= ld
                    gcf += gcl[:ld]
                surf = Surface(trclu.n - len(gcf))
                surf.cell = np.array(list(trcel))
                isurf = 0
                ffix = []
                clus = ClusterAtSurface(len(gcf), surf)
                clus.label = trclu.label
                iclus = 0
                for i in range(trclu.n):
                    if i not in gcf:
                        surf.elems[isurf] = trclu.elems[i]
                        surf.atoms[isurf] = trclu.atoms[i]
                        if i in trfix:
                            ffix.append(isurf)
                        isurf += 1
                    else:
                        clus.elems[iclus] = trclu.elems[i]
                        clus.atoms[iclus] = trclu.atoms[i]
                        iclus += 1
                print (' SFTOT  = %d' % surf.n)
                print (' SFFIX  = %s' % num_list(ffix))
            else:
                surf = Surface(trclu.n)
                surf.cell = np.array(list(trcel))
                surf.atoms = trclu.atoms
                surf.elems = trclu.elems
                ffix = trfix
        if fix is not None:
            gfix = fix
            ffix = rev_num_list(fix)
        else:
            gfix = num_list(ffix)
        if cel is not None:
            surf.cell = np.array(cel)
        if surf.cellz == 0.0:
            surf.cellz = surf.cell[-1]
        # make sure surface atoms are not above cluster atoms
        if clus is not None:
            czmax = clus.atoms[:, 2].max(axis=0)
            for i in surf.atoms:
                if i[2] > czmax:
                    i[2] = i[2] - surf.cellz
        # gapz = 0 does not give bulk results!
        if gapz is not None:
            surf.cell[-1] = (surf.atoms.max(axis=0) - surf.atoms.min(axis=0))[2] + gapz
        # surface translation has influence on space group det
        csdcell = np.array(surf.cell)
        csdir = to_direct(surf.atoms, csdcell)
        csdmin = csdir.min(axis=0).reshape((1, 3))
        csdir -= csdmin
        surf.atoms = to_cartesian(csdir, csdcell)
        surf.atoms[:, 0:2] += np.array([gdx, gdy])
        # smin = surf.atoms.min(axis=0).reshape((1, 3))
        # surf.atoms -= smin
        # remove bottom layers/sort atoms
        sgaet = []
        for ii, (ia, ie) in enumerate(zip(surf.atoms, surf.elems)):
            if removez is not None and ia[2] <= removez:
                continue
            sgaet.append([ia, ie, ii])
        sgaet.sort(key=lambda x: [x[0][2], x[0][0], x[0][1]])
        nrem = surf.n - len(sgaet)
        surf.n = len(sgaet)
        surf.atoms = np.array([i[0] for i in sgaet])
        surf.elems = np.array([i[1] for i in sgaet])
        if removez is not None:
            print (' %d ATOMS REMOVED!' % nrem)
            surf.cell[-1] -= removez
            surf.cellz -= removez
        # revise fix list after z-sorting
        if len(ffix) != 0:
            fkix = ffix
            ffix = []
            for i in fkix:
                found = None
                for ij, j in enumerate(sgaet):
                    if j[2] == i:
                        found = ij
                        break
                if found is not None:
                    ffix.append(found)
            gfix = num_list(ffix)
            print (' SFFIX* = %s' % num_list(ffix))
            print (' %d ATOMS FIXED!' % len(ffix))
        # fix bottom layers
        if fixz is not None:
            ffix = [ii for ii, i in enumerate(surf.atoms) if i[2] <= fixz]
            gfix = num_list(ffix)
            print (' SFTOT  = %d' % surf.n)
            print (' SFFIX  = %s' % num_list(ffix))
            print (' %d ATOMS FIXED!' % len(ffix))
        if clus is not None:
            cdir = to_direct(clus.atoms, csdcell)
            cdir -= csdmin
            clus.atoms = to_cartesian(cdir, csdcell)
            # clus.atoms -= smin
        surf.fix = gfix
        if ffix != [] and not alldet:
            surfex = copy.deepcopy(surf)
            surfex.atoms = np.array([a for ia, a in enumerate(surf.atoms) if ia in ffix])
            surfex.elems = np.array([a for ia, a in enumerate(surf.elems) if ia in ffix])
            surfex.n = len(surfex.atoms)
        else:
            surfex = surf
        ucellmat = determine_unit_cell(surfex)
        ucellmat[2, 2] = surf.cellz
        surf.unit_cell = to_dcell(ucellmat)
        print (' UCELL = ' + ('%15.7f' * len(surf.unit_cell)) % tuple(surf.unit_cell))
        if not nospg:
            sgd = determine_space_group(surfex, ucellmat)[0]
            surf.space_group = sgd[0]
            surf.space_group_ref = sgd[3]
            print (' UCLAT = %d' % sgd[4])
        else:
            surf.space_group = "P 1"
            surf.space_group_ref = np.zeros((2, ))
        print (' SPGRP = %s' % surf.space_group)
        print (' SGREF = ' + '%15.7f%15.7f' % tuple(surf.space_group_ref))
        sgopt = " ".join(["[%s,%s]" % (i[0], i[1]) for i in SymOpsHall[surf.space_group]])
        print (' SYMMS = %s' % sgopt)
        if supercell is not None:
            icells = [int(x) for x in supercell.split(":")]
            assert len(icells) == 3
            cellmat = to_cellmat(surf.cell)
            cellmat[2, 2] = surf.cellz
            allz = surf.cell[-1]
            if icells[2] != 1:
                clus.atoms[:, 2] += (icells[2] - 1) * surf.cellz
                allz += (icells[2] - 1) * surf.cellz
            natoms = []
            nelems = []
            for i in range(icells[0]):
                for j in range(icells[1]):
                    for k in range(icells[2]):
                        natoms.extend(list(surf.atoms + cellmat[0] * i + \
                            + cellmat[1] * j + cellmat[2] * k))
                        nelems.extend(list(surf.elems))
            for i in range(3):
                cellmat[i] *= icells[i]
            ffix = rev_num_list(surf.fix)
            if ffix != [-1] and len(ffix) != 0:
                nf = 0
                nffix = []
                for i in range(icells[0]):
                    for j in range(icells[1]):
                        for k in range(icells[2]):
                            nffix.extend(list(np.array(ffix) + nf))
                            nf += surf.n
                surf.fix = num_list(nffix)
            surf.n = len(natoms)
            surf.atoms = np.array(natoms)
            surf.elems = np.array(nelems)
            ndcell = to_dcell(cellmat)
            surf.cellz = ndcell[-1]
            ndcell[-1] = allz
            surf.cell = ndcell
            print (' SPCEL = ' + ('%15.7f' * len(surf.cell)) % tuple(surf.cell))
        if clus is not None:
            surface_align(clus)
            ftn = self.surfgen_orig_name
            print ('write original struct file: ' + ftn)
            clus.write_xyz(fn=ftn, append=False)
        ftn = self.surfgen_xyz_name
        print ('write surface xyz file: ' + ftn)
        ftnj = self.surfgen_json_name
        print ('write surface json file: ' + ftnj)
        rw = RecordWritter(0)
        rw.write_surfs(ftn, ftnj, [ClusterAtSurface(0, surf)], append=False)

    def enlarge(self):
        from surface.create_periodic import enlarge
        ipf = self.ip["enlarge"]["input_file"]
        siz = self.ip["enlarge"]["size"]
        cel = self.ip["enlarge"]["cell"]
        finals = read_zip_clusters(ipf)
        for c in finals:
            enlarge(c, siz, cel)
        ftn = new_file_name(self.enl_structs_name)
        print ('write enlarged struct file: ' + ftn)
        for c in finals:
            c.write_xyz(fn=ftn)

    def draw(self):
        from meta.report import SimpleDraw
        ipf = self.ip["report"]["input_file"]
        finals = read_zip_clusters(ipf)
        if self.ip["report"].get("number", -1) != -1:
            finals = finals[:self.ip["report"]["number"]]
        if "input_props" in self.ip["report"]:
            ipp = self.ip["report"]["input_props"]
            pps = read_json(ipp)
            for c in finals:
                labels = c.label.split(":")
                stid = labels[1]
                if stid in pps:
                    c.props = pps[stid]
                    if "bader_charges" in c.props:
                        c.props["charges"] = c.props["bader_charges"]
        opf = self.ip["report"].get("output_file", self.report_name)
        report = SimpleDraw(filename=opf, finals=finals)
        copts = {}
        for key in ["surface_depth", "ratio", "plot_charge", "rotmat", "zdepth",
            "title", "plot_force", "force_factor", "force_image", "perspective"]:
            if key in self.ip["report"]:
                copts[key] = self.ip["report"][key]
        report.build(**copts)
        print ('write report file: ' + opf)

    def report(self):
        from meta.report import Report
        ipre = self.ip["report"]
        ipff = self.ip["filtering-report"]
        copts = {}
        if "surface_depth" in self.ip["report"]:
            copts["surface_depth"] = self.ip["report"]["surface_depth"]
            del ipre["surface_depth"]
        if "max_prop_count" in self.ip["report"]:
            copts["max_prop_count"] = self.ip["report"]["max_prop_count"]
            del ipre["max_prop_count"]
        report = Report(filename=self.report_name, ipff=ipff, **ipre)
        report.build(**copts)
        report.build_graph(**copts)

        from meta.energetics import multi_name, mini_name
        from meta.report import basis_short
        ftn = self.report_structs_name
        ftnj = self.report_props_name
        for dx in report.d.d:
            for g in dx.groups:
                gname = "%s.%s_%s_%s" % (dx.name, g.name, g.param['name'],
                    basis_short(g.param["basis"]))
                ftnx = ftn % gname
                print ('write reported struct file: ' + ftnx)
                for i, c in enumerate(g.finals):
                    c.label = "#%s:%s:%s" % (c.tname, multi_name(c.multiplicity), 
                        mini_name(c.minimum_type))
                    if i == 0: c.write_xyz(fn=ftnx, append=False)
                    else: c.write_xyz(fn=ftnx)
                ftnxj = ftnj % gname
                print ('write reported properties file: ' + ftnxj)
                jprops = { ("#%s" % c.tname) : c.props for c in g.finals }
                write_json(jprops, fn=ftnxj)

    def align(self):
        ipal = self.ip['align']
        if ipal["input_file"] != -1:
            if isinstance(ipal["input_file"], int):
                input_file = load_name(name=self.para_local_name, i=ipal["input_file"])
            else:
                input_file = load_name(name=ipal["input_file"])
        else:
            raise RuntimeError("The input_file must be set for stat!")

        if "creation-surface" in self.ip:
            from surface.base import ClusterAtSurface, read_surface
            from surface.surf_comp import surface_align
            ipcs = self.ip["creation-surface"]
            surf = read_surface(ipcs["surface"])
            surf.space_group = ipcs["space_group"]
            surf.space_group_ref = np.array(ipcs["space_group_ref"])
            surf.unit_cell = ipcs["unit_cell"]
            clus = read_zip_clusters(input_file)
            clus = [ClusterAtSurface.from_cluster(sc, surf) for sc in clus]
            for cc in clus:
                surface_align(cc)
        else:
            from cluster.align import Align
            align = Align(input_file=input_file)
            align.align()
            clus = align.clus

        ftn = new_file_name(self.ali_structs_name)
        print ('write aligned struct file: ' + ftn)
        for c in clus: c.write_xyz(fn=ftn)

    def stat(self):
        from cluster.stat import Stat
        ipst = self.ip['stat']
        if ipst["input_file"] != -1:
            if isinstance(ipst["input_file"], int):
                input_file = load_name(name=self.para_local_name, i=ipst["input_file"])
            else:
                input_file = load_name(name=ipst["input_file"])
        else:
            raise RuntimeError("The input_file must be set for stat!")
        ipx = {}
        for i, j in ipst.items():
            if i != 'input_file':
                ipx[i] = j
        stat = Stat(input_file, **ipx)
        distx, ema, emi = stat.stat()
        dist = { "dist": distx, "energy_max": ema, "energy_min": emi, 
            "sample_number": len(stat.clus) }
        
        print ('write statistical distribution ...')
        write_json(dist, self.stat_dist_name)

    def fit(self):
        from cluster.fitting import Fitting
        ip = self.ip
        ipx = {}
        pre_nets = []
        for i, j in ip['fitting'].items():
            if i != 'pre_nets':
                ipx[i] = j
            else:
                for net_id in j:
                    if isinstance(net_id, int):
                        k = load_data(name=self.fit_network_name, i=net_id)
                    else:
                        k = load_data(name=net_id)
                    pre_nets.append(k)
        fitting = Fitting(pre_nets=pre_nets, **ipx)

        if "outer_loop" in ip and ip["outer_loop"] == 0:
            return
        
        ipx = {}
        for i, j in ip['fitting_net'].items():
            if i != 'load_network' and i != 'dump_network':
                ipx[i] = j
        load_net = ip['fitting_net']['load_network']
        keeper = None
        if load_net != -1:
            if isinstance(load_net, int):
                keeper = load_data(name=self.fit_network_name, i=load_net)
            else:
                keeper = load_data(name=load_net)
        fitting.create_network(net_keep=keeper, **ipx)

        loop = 1
        if "outer_loop" in ip: loop = ip["outer_loop"]
        for il in range(loop):
            print ('== outer loop: %d ==' % (il + 1))
            fitting.trans_data_all()
            fitting.train()
            if fitting.net_eval.net.saved_parameters is not None:
                print ('saved parameters found.')
                fitting.net_eval.net.set_parameters(fitting.net_eval.net.saved_parameters)

        res = fitting.test()

        ip['data'] = { 'emax': fitting.emax, 'emin': fitting.emin, 
            'lmax': fitting.lmax, 'lmin': fitting.lmin, 'error': float(res) }
        
        dump_net = ip['fitting_net']['dump_network']

        ft = open(new_file_name(self.fit_error_name), 'w')
        print ('write fitting error file: ' + ft.name)
        ft.write('# %15.8f eV\n' % res)
        for i, x, y, z in zip(range(fitting.epochs), fitting.net_eval.net.train_errors, 
            fitting.net_eval.net.test_errors, fitting.net_eval.net.train_times):
            ft.write('%10d %15.8f %15.8f %7.2f\n' % (i, x, y, z))
        ft.close()

        ft = open(new_file_name(self.fit_test_name), 'w')
        print ('write test file: ' + ft.name)
        ft.write('#%7s %15s %15s %15s\n' % ('id', 'standard', 'predict', 'error'))
        for idx, std, pre, err in fitting.testings:
            ft.write('%8d %15.8f %15.8f %15.8f\n' % (idx, std, pre, err))
        ft.close()

        if "test_permutation" in ip and ip["test_permutation"]:
            fitting.test_permutation()
            ft = open(new_file_name(self.fit_pertest_name), 'w')
            print ('write permutation test file: ' + ft.name)
            ft.write('#%7s %15s %15s %15s %15s\n' % ('id', 'standard', 'predict-mean', 'std-pre error', 'predict-rmse'))
            for idx, std, pre, err, rmse in fitting.per_testings:
                ft.write('%8d %15.8f %15.8f %15.8f %15.8f\n' % (idx, std, pre, err, rmse))
            ft.close()

        if dump_net:
            print ('dump network ...')
            dump_data(name=self.fit_network_name, obj=fitting.net_eval.save_network())
    
    def opt(self):
        from cluster.opt import Optimization
        ip = self.ip

        ipop = ip["optimization"]
        print ('load network ...')
        if ipop["load_network"] != -1:
            if isinstance(ipop["load_network"], int):
                net_keep = load_data(name=self.fit_network_name, i=ipop["load_network"])
            else:
                net_keep = load_data(name=ipop["load_network"])
        else:
            raise RuntimeError("The fitted network is needed for optimization!")
        
        opt = Optimization(input_file=ipop["input_file"], net_keep=net_keep)
        ipx = {}
        # define a subset of structures to optimize
        if "opt_number" in ipop: ipx["nopt"] = ipop["opt_number"]
        opt.opt(**ipx)

        ft = open(new_file_name(self.opt_list_name), 'w')
        print ('write optimization list file: ' + ft.name)
        for l in opt.opts: ft.write(l + '\n')
        ft.close()

        if "sort" in ipop and ipop["sort"]:
            print ('sort optimized structs by energy.')
            opt.finals = sorted(opt.finals, key=lambda x: x.energy)

        ftn = new_file_name(self.opt_structs_name)
        print ('write optimized struct file: ' + ftn)
        for c in opt.finals: c.write_xyz(fn=ftn)
    
    @staticmethod
    def get_create_opts(ip):
        from cluster.base import understand_name, update_stat, Cluster
        if "creation-surface" in ip:
            from surface.base import ClusterAtSurface, read_surface
            ipcs = ip["creation-surface"]
            surf = read_surface(ipcs["surface"])
            surf.space_group = ipcs["space_group"]
            surf.space_group_ref = np.array(ipcs["space_group_ref"])
            surf.unit_cell = ipcs["unit_cell"]
        else:
            surf = None
        ipcr = ip["creation"]
        if "method" in ipcr and ipcr["method"] != 'blda':
            raise NotImplementedError('only blda generator available.')
        fl_elems, elems, et = understand_name(ipcr["name"])
        etx = et.replace(' ', '')
        defs = ipcr["default_sigma"] if "default_sigma" in ipcr else 0.04
        xorder = ipcr.get("order", 2)
        xmorder = ipcr.get("morder", 1)
        xtwod = ipcr.get("2d", 0.0)
        ext_range = ipcr.get("ext_range", None)
        xsampling = ipcr.get("sampling", "random")
        xmean, xsigma = {}, {}
        if "dist" in ip.keys(): 
            xmean, xsigma = ip["dist"]["mu"], ip["dist"]["sigma"]
        if "~hydro" in ipcr:
            fl_elems = list(fl_elems) + ['H']
        update_stat(fl_elems, xmean, xsigma, defs)
        copts = { "elems": elems, "mean": xmean, "sigma": xsigma, "order": xorder,
            "morder": xmorder, "twod": xtwod }
        for mk, mv in ipcr.items():
            if mk[0] == '~':
                copts[mk[1:]] = mv
        if ext_range is not None:
            copts["ext_range"] = ext_range
        # PES sampling, set number to any (such as 0). need "xyfix" in arg-opts
        if xsampling.startswith("grid"):
            copts["grid"] = map(int, xsampling.split(".")[1:])
        if "creation-periodic" in ip:
            del copts["order"], copts["morder"], copts["twod"]
            for mk, mv in ipcr.items():
                if mk[0] == '~':
                    del copts[mk[1:]]
            if "ext_range" in copts:
                del copts["ext_range"]
            if "max_height" in ipcs:
                copts["max_height"] = ipcs["max_height"]
            if "layer_height" in ipcs:
                copts["layer_height"] = ipcs["layer_height"]
            if "start_height" in ipcs:
                copts["start_height"] = ipcs["start_height"]
            update_stat(list(fl_elems) + list(surf.elems), xmean, xsigma, defs)
        else:
            if "creation-surface" in ip:
                del copts["twod"]
                if "ncore" in ipcs: copts["ncore"] = ipcs["ncore"]
                if "mix_rate" in ipcs: copts["mix_rate"] = ipcs["mix_rate"]
                update_stat(list(fl_elems) + list(surf.elems), xmean, xsigma, defs)
        return copts

    def filter(self, create=False):
        from cluster.opt import Filtering
        if create:
            from cluster.base import understand_name, update_stat, Cluster
        ip = self.ip
        if create:
            if "filtering-create" in self.ip:
                ipfl = self.ip["filtering-create"]
            else: ipfl = self.ip["filtering"]
        else:
            ipfl = ip["filtering"]

        if "creation-surface" in self.ip:
            from surface.base import ClusterAtSurface, read_surface
            ipcs = ip["creation-surface"]
            surf = read_surface(ipcs["surface"])
            surf.space_group = ipcs["space_group"]
            surf.space_group_ref = np.array(ipcs["space_group_ref"])
            surf.unit_cell = ipcs["unit_cell"]
        else:
            surf = None
        prej = None
        if create:
            ipcr = ip["creation"]
            copts = Driver.get_create_opts(self.ip)
            if "method" in ipcr and ipcr["method"] != 'blda':
                raise NotImplementedError('only blda generator available.')
            fl_elems, elems, et = understand_name(ipcr["name"])
            etx = et.replace(' ', '')
            clnum = ipcr["number"]
            xsampling = ipcr.get("sampling", "random")
            # PES sampling, set number to any (such as 0). need "xyfix" in arg-opts
            if xsampling.startswith("grid"):
                clnum = copts["grid"][0] * copts["grid"][1]
                ipfl["max_diff"] = 0.0
                ipfl["max_diff_report"] = 0.0
            if "creation-periodic" in self.ip:
                print ('Periodic Mode!!')
                from surface.create_periodic import CreatePeriodic
                xcreate = CreatePeriodic.generator(surf=surf, **copts)
            else:
                from surface.create_cluster import CreateCluster
                from cluster.base import moles
                ok = None
                if len(elems) > 1 and isinstance(elems[0], np.ndarray) and len(elems[0]) == 1 \
                    and elems[0][0] in moles:
                    if moles[elems[0][0]].surfnum == 0 and surf is not None:
                        ok = None
                    elif moles[elems[0][0]].surfnum == 0 and surf is None:
                        ok = moles[elems[0][0]]
                        copts["elems"] = elems[1:]
                    elif moles[elems[0][0]].surfnum != 0 and surf is None:
                        raise RuntimeError("Adding atoms to existing surface/cluster structure " + 
                            "requires surface information in input file!")
                    elif moles[elems[0][0]].surfnum != 0 and surf is not None:
                        from cluster.base import MoleculePool
                        if moles[elems[0][0]].ismulti:
                            cs = [ClusterAtSurface.from_cluster(csm, surf) for csm in moles[elems[0][0]].mm]
                            ok = MoleculePool(cs)
                        else:
                            ok = MoleculePool(ClusterAtSurface.from_cluster(moles[elems[0][0]], surf))
                        copts["elems"] = elems[1:]
                if ok is None:
                    xcreate = CreateCluster.generator(surf=surf, **copts)
                else:
                    if ok.ismulti and xsampling == "all":
                        ok.smode = "all"
                        ok.scur = 0
                        ok.snum = clnum
                        prej = ok.reject
                        clnum *= len(ok.mm)
                    elif ok.ismulti and xsampling == "random":
                        ok.smode = "random"
                    xcreate = CreateCluster.generator(surf=surf, xn=ok.n, xpool=ok, **copts)
            print ('blda parameters: ')
            print_dict(copts, group=1)
        else:
            if ipfl["input_file"] != -1:
                if isinstance(ipfl["input_file"], int):
                    input_file = load_name(name=self.opt_structs_name, i=ipfl["input_file"])
                else:
                    input_file = load_name(name=ipfl["input_file"])
            else:
                raise RuntimeError("The input_file must be set for filtering!")
        ipx = {}
        ipy = {}
        for i, j in ipfl.items():
            if i != 'input_file' and i != 'sort' and i != 'pre_sort' and i != 'align' and i != 'para':
                ipx[i] = j[-1] if isinstance(j, list) else j
            elif i == 'pre_sort' and not create:
                ipy[i] = j
        
        filt = Filtering(**ipx)
        if create:
            filt.init_from_create(citer=xcreate, cn=clnum, an=len(elems), prej=prej)
            filt.filter(prefix='INI')
        else:
            filt.init_from_file(filename=input_file, **ipy)
            if "creation-surface" in self.ip:
                filt.clus = [ClusterAtSurface.from_cluster(sc, surf) for sc in filt.clus]
            filt.filter()

        if "sort" in ipfl and ipfl["sort"] and not create:
            print ('sort filtered structs by energy.')
            filt.finals = sorted(filt.finals, key=lambda x: x.energy)

        if "align" in ipfl and ipfl["align"]:
            if "creation-surface" in self.ip:
                from surface.surf_comp import surface_align
                for cc in filt.finals:
                    surface_align(cc)
            else:
                from cluster.align import Align
                align = Align()
                align.clus = filt.finals
                align.align()

        if "para" in ipfl and ipfl["para"]:
            fidx = ipfl["input_file"].split(".")[-1]
            flist = self.paraf_list_name + "." + fidx
            fcorr = self.paraf_corr_name + "." + fidx
            fstru = self.paraf_structs_name + "." + fidx
        else:
            flist = new_file_name(self.fil_list_name)
            fcorr = new_file_name(self.fil_corr_name)
            fstru = new_file_name(self.fil_structs_name)

        ft = open(flist, 'w')
        print ('write filtering list file: ' + ft.name)
        ft.write('#%7s %10s %10s %8s %15s\n' % ('id', 'fil-id', 'orig-id', 'multi', 'fit-energy'))
        for idx, f in enumerate(filt.finals):
            ft.write('%8d %10s %10s %8d %15.8f\n' % (idx, f.new_label, f.label, f.multi, f.energy))
        ft.close()

        ft = open(fcorr, 'w')
        print ('write filtering correspondence file: ' + ft.name)
        for l in filt.corrs:
            ft.write(l + '\n')
        ft.close()

        ftn = fstru
        print ('write filtered struct file: ' + ftn)
        apd = False
        for c in filt.finals:
            c.write_xyz(fn=ftn, append=apd)
            apd = True

        return filt.finals
