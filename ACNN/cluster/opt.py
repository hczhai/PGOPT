
import itertools, theano, theano.tensor as T, numpy as np, time
from cluster.fitting import NetEvaluator
from cluster.fitting import get_length, trans_forward, trans_backward
from cluster.base import read_zip_clusters, Cluster, elem_num
from utils.base import cast_float
from formod.lbfgs import LBFGS
from formod.at_comp import at_comp
from surface.base import ClusterAtSurface
from surface.surf_comp import surface_compare

class Optimization(object):

  def __init__(self, input_file, net_keep):

    self.net_keep = net_keep

    # could be zip or normal serial xyz
    # auto detected
    self.clus = read_zip_clusters(input_file)
    self.cn = len(self.clus)
    self.an = self.clus[0].n
    self.net_eval = NetEvaluator(self.clus)
    self.net_eval.load_network(net_keep)
    self.init_derivatives()
    self.post_reject = None

    self.opt_data = [ None ] * 2
    self.opt_data[0] = np.array([c.atoms for c in self.clus])
    self.opt_data[1] = np.array([c.energy for c in self.clus])
  
  def init_derivatives(self):
    print ('initialize derivative functions ...')
    n = self.an
    num = self.net_keep.degree
    expl = self.net_keep.expl
    scale_l = self.net_keep.scale_l
    d_order = self.net_keep.d_order
    le = np.array([list(g) for g in itertools.combinations(range(0, n), num)])
    x = cast_float(T.dvector())
    xr = x.reshape((n, 3))
    lx = get_length(num)
    lena = le.shape[0] # number of all selections
    lenb = lx.shape[1] # for each selection, number of lengths
    yl = le[:, lx.T] # list of integers need to be selected for lengths
    lt = get_length(n)
    lenc = lt.shape[1] # for all atoms, number of lengths
    lti = np.zeros((n, n), dtype=np.int)
    for i, (j, k) in enumerate(lt.T): lti[j, k] = lti[k, j] = i
    if d_order: 
      raise NotImplementedError('second order has not been supported yet!')
    mt = (xr[lt[0]] - xr[lt[1]]).norm(2, axis=1) # all lengths
    if expl != 0: mt = np.exp(-mt / expl)
    if scale_l: mt = trans_forward(mt, *self.net_keep.max_min[0:2])
    ylt = lti[yl[:, :, 0], yl[:, :, 1]] # transfer the len indices to index in mt
    y = mt[ylt] # all length values after num'd, it is the input of net
    ylf = yl.reshape((lena * lenb, 2))
    yltf = lti[ylf[:, 0], ylf[:, 1]] # flatted ylt
    mtdl = [] # calculate the grad of all lengths wrt x
    for i in range(lenc): mtdl.append(T.grad(mt[i], x).reshape((1, n * 3)))
    mtd = T.concatenate(mtdl, axis=0) # the grad of mt wrt flat x
    yd = mtd[yltf] # the grad with sel-permu (flatted) indices
    xip = self.net_eval.net.variables['net_input']
    xop = self.net_eval.net.variables['prediction'][0][0]
    xop = trans_backward(xop, *self.net_keep.max_min[2:4])
    print ('## F: lengths -> energy ...')
    xener = theano.function([xip], xop) # from input lengths to energy
    print ('## dF: lengths -> energy ...')
    xenerd = theano.function([xip], T.grad(xop, xip).reshape((lena * lenb, ))) # energy gradient
    print ('## F: coords -> lengths ...')
    xprim = theano.function([x], y.reshape((1, lena, lenb))) # from atomic to input
    print ('## dF: coords -> lengths ...')
    xprimd = theano.function([x], yd) # from atomic to input, gradient
    self.opt_eval = lambda x, xener=xener, xprim=xprim: xener(xprim(x))
    self.opt_evald = (lambda x, xprim=xprim, xprimd=xprimd, xenerd=xenerd: 
      np.tensordot(xenerd(xprim(x)), xprimd(x), axes=(0, 0)))
  
  def test_opt(self, c):
    lt = get_length(c.n)
    mt = np.linalg.norm(c.atoms[lt[0]] - c.atoms[lt[1]], axis=1)
    rto = 0.05
    dr = rto / (1 + 2 * rto) * (self.net_keep.max_min[0] - self.net_keep.max_min[1])
    rmin = self.net_keep.max_min[1] + dr
    rmax = self.net_keep.max_min[0] - dr
    if mt.min() < rmin or mt.max() > rmax: return False
    else: return True

  def opt(self, nopt=-1):
    print ('optimization ...')
    task = LBFGS(self.an * 3)
    cast_forward = np.cast[np.float64]
    cast_backword = np.cast[theano.config.floatX]
    cast = lambda f: lambda x: cast_forward(f(cast_backword(x)))
    task.p.eval = cast(self.opt_eval)
    task.p.evald = cast(self.opt_evald)
    task.log_file = 0
    task.max_iter = 500
    self.opts = []
    self.finals = []
    elems = self.clus[0].elems
    opt_data = self.opt_data
    nopt = opt_data[0].shape[0] if nopt == -1 else min(opt_data[0].shape[0], nopt)
    holder = 65 + 16
    print ('Number of structs: %d' % nopt)
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    self.opts.append('#%9s %3s %15s %15s %15s %7s %7s  ' % ('struct #', '%', 'start E', 
      'final E', 'delta E', 'steps', 'time'))
    print (self.opts[-1])
    self.min_ener = None
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    init_time = start_time = time.time()
    for idx, x, y in zip(range(0, len(opt_data[0])), opt_data[0], opt_data[1]):
      if idx == nopt: break
      task.start(cast_forward(x.flatten()))
      task.opt()
      fo, ff = task.traj[0][1], task.traj[-1][1]
      if self.clus[idx].label == 'Energy': label = 'OPT:{}'.format(idx)
      else: label = 'OPT:{}'.format(self.clus[idx].label)
      cc = Cluster(n=self.an, atoms=cast_backword(task.x.reshape(self.an, 3)), 
        elems=elems, energy=ff, label=label)
      is_min = ' '
      if self.test_opt(cc):
        self.finals.append(cc)
        if self.min_ener is None or ff < self.min_ener:
          is_min = '*'
          self.min_ener = ff
      else:
        is_min = 'F'
      finish_time = time.time()
      self.opts.append('%10d %2d%% %15.8f %15.8f %15.8f %7d %7.2f%s' % (idx, idx * 100 / nopt, 
        fo, ff, ff - fo, len(task.traj), finish_time - start_time, ' %s' % is_min))
      print (self.opts[-1])
      start_time = finish_time
      if idx == 9:
        total_time = (finish_time - init_time) / 10 * nopt
        print (('%%%ds' % (holder)) % ('estimated total time: %d:%02d:%02d' % (total_time / 3600, 
          (total_time / 60) % 60, total_time % 60), ))
        print (('%%%ds' % (holder)) % ('-'*holder, ))
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    total_time = finish_time - init_time
    print (('%%%ds' % (holder)) % ('total time used: %d:%02d:%02d' % (total_time / 3600, 
      (total_time / 60) % 60, total_time % 60), ))
    print (('%%%ds' % (holder)) % ('-'*holder, ))

class Filtering(object):

  def __init__(self, max_diff=0.25, max_diff_report=1.00):

    self.max_diff = max_diff
    self.max_diff_report = max_diff_report
    self.post_reject = None
  
  def init_from_file(self, filename, pre_sort=False):
    self.clus = read_zip_clusters(filename)
    self.cn = len(self.clus)
    # self.an = self.clus[0].n
    if pre_sort:
      print ('pre-sort structs by energy for filtering.')
      self.clus = sorted(self.clus, key=lambda x: x.energy)
    self.flimit = False
  
  def init_from_create(self, citer, cn, an, prej=None):
    self.clus = citer
    self.cn = cn
    # self.an = an
    self.flimit = True
    self.post_reject = prej
  
  def filter(self, prefix='FIL', iprint=True):
    if iprint: print ('filtering ...')
    self.corrs = []
    n = self.cn
    dmax_rep = self.max_diff_report
    dmax = self.max_diff
    eles = None
    self.finals = []
    holder = 76
    if iprint: print ('Number of structs: %d' % n)
    if iprint: print (('%%%ds' % (holder)) % ('-'*holder, ))
    self.corrs.append('#%6s %3s %-8s %15s %-8s %15s %6s %5s  ' % 
      ('new #', '%', 'orig #', 'this E', 'sim. #', 'sim. E', 'mindm', 'time'))
    if iprint: print (self.corrs[-1])
    self.min_ener = None
    if iprint: print (('%%%ds' % (holder)) % ('-'*holder, ))
    init_time = start_time = time.time()
    for ci, c in enumerate(self.clus):
      if c.label == '': c.label = 'PRE:%d' % (ci, )
      if eles is None:
        eles, ne = elem_num(c.elems)
      # if n >= 10 and ci % (n / 10) == 0: print '{:.0f} %'.format(ci / (n / 100.0))
      min_dmax, mini = dmax_rep, None
      ok = True
      for di, d in enumerate(self.finals):
        if min_dmax == 0.0:
          v = 1.0
        elif isinstance(c, ClusterAtSurface):
          v, _ = surface_compare(d, c, ne, eles, min_dmax)
        else:
          v, _ = at_comp(d.atoms, c.atoms, ne, eles, min_dmax)
        if v < min_dmax:
          min_dmax, mini = v, d
        if v < dmax:
          ok = False
          if self.post_reject is not None:
            self.post_reject()
          break
      finish_time = time.time()
      self.corrs.append(('{:7} {:2d}% {:8} {:15.6f} {:24} {mindm:6.4f} ' + 
        '{time:5.2f} {star}').format(len(self.finals) if ok else '', 
        int((len(self.finals) if self.flimit else ci) / ((n + 1 if self.flimit else n) / 100.0)), 
        c.label, c.energy, '{:8} {:15.6f}'.format(mini.label, mini.energy) 
        if mini is not None else '{:8} {:>15}'.format('-', '-'), 
        mindm=min_dmax, time=finish_time - start_time, star='*' if ok else ' '))
      start_time = finish_time
      if iprint: print (self.corrs[-1])
      if ok: 
        self.finals.append(c)
        c.multi = 1
        c.new_label = prefix + ':' + str(len(self.finals) - 1)
      else:
        mini.multi += 1
        if hasattr(mini, "props") and "nrep" in mini.props:
          mini.props["nrep"] += c.props["nrep"]
      if len(self.finals) == n: break
      if ci >= 100000:
        raise RuntimeError("Cannot generate enough " + 
          "distinct structures after %d trials!" % ci)
    if iprint: print (('%%%ds' % (holder)) % ('-'*holder, ))
    total_time = finish_time - init_time
    if iprint: print (('%%%ds' % (holder)) % ('total time used: %d:%02d:%02d' % (total_time / 3600, 
      (total_time / 60) % 60, total_time % 60), ))
    if iprint: print (('%%%ds' % (holder)) % ('-'*holder, ))