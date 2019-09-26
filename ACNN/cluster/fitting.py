
import numpy as np
import random, itertools, copy
from cluster.base import read_zip_clusters
from acnn.trans_layers import trans_layers, LayerList
from acnn.learn import MomentumNet
from acnn.learn_adv import LevenbergMarquardtNet, LBFGSNet
from acnn.layers import CoLayerSize as C, ParameterLayer, Reshape, AxisPoolLayer, Diff
from utils.base import cast_float
from utils.base import reproducible, print_dict

# length indices
def get_length(n, selfinc=False):
  m = n * (n - 1) / 2
  zz = np.zeros((2, m), dtype=int)
  mi = 0
  for i in range(0, n):
    for j in range(0, i + 1 if selfinc else i):
      zz[:, mi] = i, j
      mi += 1
  return zz

def find_max_min(data, min_max_ext_ratio):
  ratio = min_max_ext_ratio
  dmax = np.ma.max(data)
  dmin = np.ma.min(data)
  dmm = dmax - dmin
  dmin -= dmm * ratio
  dmax += dmm * ratio
  return dmax, dmin

def trans_backward(ener, dmax, dmin):
  return ener * (dmax - dmin) + dmin

def trans_forward(ener, dmax, dmin):
  return (ener - dmin) / (dmax - dmin)

# size of input
def ipsize_inner(n, d_order=False):
  size = n * (n - 1) / 2
  if d_order: size = size + size * (size + 1) / 2
  return size

def get_ip_size(an, degree, d_order=False):
  if degree == 0: return (an, an - 1)
  lx = [list(g) for g in itertools.combinations(range(0, an), degree)]
  return (len(lx), ipsize_inner(degree, d_order))

# generate training/test length data
def trans_data(clul, n, degree, d_order, expl, diff, new_ener, 
  randorder=True, shuffle=True):
  num = degree
  print ('data preparing ...')
  xn = len(clul)
  gn = clul[0].n
  if num > 0:
    lx = [list(g) for g in itertools.combinations(range(0, gn), num)]
    gl = get_length(num)
    lend = len(gl[0])
    if d_order:
      gls = get_length(lend, selfinc=True)
    if not diff: x_d = np.zeros((n, len(lx), ipsize_inner(num, d_order)))
    else: x_d = np.zeros((n, 2, len(lx), ipsize_inner(num, d_order)))
  else:
    gl = np.array([[[i, j] for j in range(0, gn) if j != i]
      for i in range(0, gn)]).swapaxes(2, 0).swapaxes(1, 2)
    assert diff == False and d_order == False
    x_d = np.zeros((n, gn, gn - 1))
  y_d = np.zeros((n, ))
  for i in range(n):
    if n >= 10 and i % (n / 10) == 0: print '{:.0f} %'.format(i / (n / 100.0))
    if randorder:
      inds = [random.randrange(xn)]
    else:
      inds = [i / (n / xn)]
    if diff:
      if np.random.random() > 0.5: inds.append(inds[0])
      else: inds.append(random.randrange(xn))
    for j, ind in enumerate(inds):
      cl = clul[ind]
      if shuffle: cl.shuffle()
      if num == 0:
        if expl == 0: x_d[i] = np.linalg.norm(cl.atoms[gl[0]] - 
          cl.atoms[gl[1]], axis=2)
        else: x_d[i] = np.exp(-np.linalg.norm(cl.atoms[gl[0]] - 
          cl.atoms[gl[1]], axis=2) / expl)
        y_d[i] = cl.new_energy if new_ener else cl.energy
        continue
      xx = cl.atoms[np.array(lx)]
      if expl == 0:
        xx = np.linalg.norm(xx[:, gl[0]] - xx[:, gl[1]], axis=2)
      else:
        xx = np.exp(-np.linalg.norm(xx[:, gl[0]] - xx[:, gl[1]], axis=2) / expl)
      if not diff:
        x_d[i, :, 0:lend] = xx
        if d_order: 
          x_d[i, :, lend:] = xx[:, gls[0]] * xx[:, gls[1]]
        y_d[i] = cl.new_energy if new_ener else cl.energy
      else:
        x_d[i, j, :, 0:lend] = xx
        if d_order: 
          x_d[i, j, :, lend:] = xx[:, gls[0]] * xx[:, gls[1]]
        if j == 0: y_d[i] = cl.new_energy if new_ener else cl.energy
        else: 
          od = y_d[i]
          y_d[i] = cl.new_energy if new_ener else cl.energy
          y_d[i] -= od
  return x_d, y_d

# parameters needed for mixed layer
def mixed_param(elems, num):
  lelems = [list(g) for g in itertools.combinations(elems, num)]
  lelems = list(set(lelems))
  if len(lelems) != 1:
    leldic = {lelems[i]: i for i in range(len(lelems))}
    lxel = [leldic[[elems[g] for g in l]] for l in lxel]
  else:
    lxel = None
  return len(lelems), lxel

class NetBookkeeping(object):
  def __init__(self, layers_raw, params, max_min, elems, degree, d_order, expl, scale_l):
    self.params = params
    self.max_min = max_min
    self.degree = degree
    self.d_order = d_order
    self.expl = expl
    self.elems = elems
    self.layers_raw = layers_raw
    self.scale_l = scale_l

class NetEvaluator(object):
  def __init__(self, clus, diff=False, borrow=False):
    self.clus = clus if borrow else copy.deepcopy(clus)
    self.cn = len(clus)
    self.an = clus[0].n
    self.diff = diff
    self.net = None
    self.net_keep = None
    self.ener_pre = None
  
  def load_network(self, net_keep, **train_opts):
    ip_size = get_ip_size(len(net_keep.elems), net_keep.degree, net_keep.d_order)
    layers = trans_layers(net_keep.layers_raw, ip_size)
    if len(train_opts) == 0:
      self.net = MomentumNet(layers=layers)
    else:
      all_key = ["mu", "mu_update_factor", "nesterov", 
        "momentum", "method", "step", "batch_size", "function", 
        "h0_scale", "gradient_tol"]
      method_key = {"momentum": ["nesterov", "momentum", "step", "batch_size"], 
        "levenberg": ["mu", "mu_update_factor"], 
        "quasi_newton": ["function", "h0_scale", "gradient_tol"], 
        "lbfgs": [] }
      cc_opts = {}
      m = train_opts["method"] if "method" in train_opts else "momentum"
      for k, v in train_opts.items():
        if k not in all_key or k in method_key[m]: cc_opts[k] = v
      if m == "momentum":
        self.net = MomentumNet(layers=layers, **cc_opts)
      elif m == "levenberg":
        self.net = LevenbergMarquardtNet(layers=layers, **cc_opts)
      elif m == "lbfgs":
        self.net = LBFGSNet(layers=layers, **cc_opts)
    if net_keep.params is not None:
      self.net.set_parameters(net_keep.params)
    self.net_keep = net_keep
  
  def save_network(self):
    if self.net.param_save != 0 and self.net.saved_parameters is not None:
      params = self.net.saved_parameters
    else:
      params = self.net.get_parameters()
    self.net_keep.params = params
    return self.net_keep
  
  # this is only useful for determine max, min
  # not for determine the energy difference to train
  def predict_all(self):
    xd, yd = trans_data(self.clus, self.cn, degree=self.net_keep.degree, 
      d_order=self.net_keep.d_order, expl=self.net_keep.expl, diff=False, 
      randorder=False, shuffle=False, new_ener=False)
    if self.net_keep.scale_l:
      xd = trans_forward(xd, self.net_keep.max_min[0], self.net_keep.max_min[1])
    yp = self.net.predict(xd).reshape(self.cn)
    self.ener_pre = trans_backward(yp, self.net_keep.max_min[2], self.net_keep.max_min[3])

class Fitting(object):

  httoev = 27.21138505
  def __init__(self, input_file,  sample_number, sample_ratio, 
    energy_cut=0.0, exp_length=0, min_max_ext_ratio=0.05, 
    scale_lengths=True, second_order=False, energy_unit="hartree", 
    shuffle_input=False, degree=0, param_save=0, diff=False, orig_and_diff=0, pre_nets=[]):

    self.energy_cut = energy_cut
    self.sample_number = sample_number
    self.sample_ratio = sample_ratio
    self.exp_length = exp_length
    self.min_max_ext_ratio = min_max_ext_ratio
    self.scale_lengths = scale_lengths
    self.second_order = second_order
    self.energy_factor = Fitting.httoev if energy_unit == "hartree" else 1.0
    self.shuffle_input = shuffle_input
    self.param_save = param_save
    self.diff = diff
    self.degree = degree
    self.orig_and_diff = orig_and_diff
    self.pre_net_evals = []

    self.clus = read_zip_clusters(input_file)
    if isinstance(energy_cut, float):
      self.clus = [c for c in self.clus if c.energy < energy_cut]
    elif isinstance(energy_cut, basestring) and energy_cut.startswith("+"):
      ec = float(energy_cut[1:])
      emin = min([c.energy for c in self.clus])
      emax = max([c.energy for c in self.clus])
      print ('Original Energy: max = %15.6f, min = %15.6f' % (emax, emin))
      print ('Energy difference: %15.6f, scaled = %15.6f\n' % (emax - emin, 
        self.energy_factor * (emax - emin)))
      self.clus = [c for c in self.clus if c.energy < emin + ec]
      emin = min([c.energy for c in self.clus])
      emax = max([c.energy for c in self.clus])
      emean = np.array([c.energy for c in self.clus]).mean()
      estd = np.array([c.energy for c in self.clus]).std()
      print ("Structures selected: %d (%s)" % (len(self.clus), energy_cut))
      print ('Selected Energy: max = %15.6f, min = %15.6f' % (emax, emin))
      print ('Scaled: mean = %15.6f, std = %15.6f' % (self.energy_factor * emean, 
        self.energy_factor * estd))
      print ('Energy difference: %15.6f, scaled = %15.6f\n' % (emax - emin, 
        self.energy_factor * (emax - emin)))
    for c in self.clus:
      c.new_energy = c.energy
    if self.shuffle_input: random.shuffle(self.clus)
    self.cn = self.clus[0].n
    aclus = np.array([c.atoms for c in self.clus])
    plen = get_length(self.cn)
    al = np.linalg.norm(aclus[:, plen[0]] - aclus[:, plen[1]], axis=2)
    ae = np.array([c.energy for c in self.clus])
    self.emax, self.emin = find_max_min(ae, min_max_ext_ratio)
    self.lmax, self.lmin = find_max_min(al, min_max_ext_ratio)
    print ('Length: max = %15.6f, min = %15.6f' % (self.lmax, self.lmin))
    print ('Energy: max = %15.6f, min = %15.6f' % (self.emax, self.emin))
    print ('Energy difference: %15.6f, scaled = %15.6f' % (self.emax - self.emin, 
      self.energy_factor * (self.emax - self.emin)))
    if len(pre_nets) != 0:
      print ('Found %d pre-nets ...' % (len(pre_nets)))
      for net_keep in pre_nets:
        ne = NetEvaluator(self.clus)
        ne.load_network(net_keep)
        ne.predict_all()
        for ic, c in enumerate(self.clus):
          c.new_energy -= ne.ener_pre[ic]
        self.pre_net_evals.append(ne)
      nae = np.array([c.new_energy for c in self.clus])
      self.emax, self.emin = find_max_min(nae, min_max_ext_ratio)
      print ('New energy: max = %15.6f, min = %15.6f' % (self.emax, self.emin))
      print ('New energy difference: %15.6f, scaled = %15.6f' % (self.emax - self.emin, 
        self.energy_factor * (self.emax - self.emin)))
    self.net_eval = None

  def trans_data_all(self):
    ssr = sum(self.sample_ratio)
    nclus = len(self.clus)
    self.fit_data = []
    if self.orig_and_diff != 0 and self.diff:
      self.ori_fit_data = []
    rstart = 0.0
    for i in range(0, 3):
      ratio = self.sample_ratio[i] / ssr
      rend = rstart + ratio
      nstart, nend = int(nclus * rstart), int(nclus * rend)
      rstart = rend
      nseed = random.randrange(2**32)
      reproducible(seed=nseed)
      x, y = trans_data(self.clus[nstart:nend], self.sample_number[i], degree=self.degree, 
        d_order=self.second_order, expl=self.exp_length, diff=self.diff, 
        randorder=True, shuffle=True, new_ener=False)
      if self.scale_lengths:
        x = trans_forward(x, self.lmax, self.lmin)
      for ev in self.pre_net_evals:
        k = ev.net_keep
        reproducible(seed=nseed)
        xsr, ysr = trans_data(self.clus[nstart:nend], self.sample_number[i], degree=k.degree, 
          d_order=k.d_order, expl=k.expl, diff=self.diff, 
          randorder=True, shuffle=True, new_ener=False)
        print xsr.shape
        assert (ysr == y).all()
        if k.scale_l:
          xsr = trans_forward(xsr, k.max_min[0], k.max_min[1])
        ysr = ev.net.predict(xsr).reshape(xsr.shape[0])
        ysr = trans_backward(ysr, k.max_min[2], k.max_min[3])
        y -= ysr
      if not self.diff: y = trans_forward(y, self.emax, self.emin)
      else: y = trans_forward(y, self.emax - self.emin, 0)
      self.fit_data += [x, y]
      if self.orig_and_diff != 0 and self.diff:
        x_ori, y_ori = trans_data(self.clus[nstart:nend], self.sample_number[i], 
          diff=False, degree=self.degree, d_order=self.second_order, expl=self.exp_length, 
          randorder=True, shuffle=True, new_ener=True)
        y_ori = trans_forward(y_ori, self.emax, self.emin)
        if self.scale_lengths:
          x_ori = trans_forward(x_ori, self.lmax, self.lmin)
        self.ori_fit_data += [x_ori, y_ori]
    print ('##   training data shape: {} > {}'.format(self.fit_data[0].shape, self.fit_data[1].shape))
    print ('##    testing data shape: {} > {}'.format(self.fit_data[2].shape, self.fit_data[3].shape))
    print ('## validation data shape: {} > {}'.format(self.fit_data[4].shape, self.fit_data[5].shape))

  def create_network(self, layers=None, net_keep=None, epochs=100, **train_opts):

    if net_keep is None:
      net_keep = NetBookkeeping(layers_raw=layers, params=None, 
        max_min=[self.lmax, self.lmin, self.emax, self.emin], elems=self.clus[0].elems, 
        degree=self.degree, d_order=self.second_order, expl=self.exp_length, 
        scale_l=self.scale_lengths)
    else:
      # if these parameters changed, we should use new parameters
      net_keep.max_min = [self.lmax, self.lmin, self.emax, self.emin]
      net_keep.elems = self.clus[0].elems
      assert net_keep.degree == self.degree
      assert net_keep.d_order == self.second_order
      assert net_keep.expl == self.exp_length
      assert net_keep.scale_l == self.scale_lengths
    train_opts.update({'param_save': self.param_save})
    self.epochs = epochs
    
    self.net_eval = NetEvaluator(self.clus, diff=self.diff, borrow=True)
    self.net_eval.load_network(net_keep, **train_opts)
    if self.diff:
      self.net_ori_eval = NetEvaluator(self.clus, diff=False, borrow=True)
      self.net_ori_eval.load_network(net_keep, **train_opts)
    
    print ('network data parameters:')
    print_dict(net_keep.__dict__)
    print ('network training parameters:')
    print_dict(train_opts)
    print ('   EPOCHS = {}'.format(self.epochs))
    if self.diff: print ('!!! using difference network !!!')

  def train(self):
    scale = (self.emax - self.emin) * self.energy_factor
    if self.orig_and_diff == 0 or not self.diff:
      self.net_eval.net.train(*self.fit_data[0:4], epochs=self.epochs, scale=scale)
    else:
      dk = self.orig_and_diff
      self.net_eval.net.current_epoch = 0
      self.net_ori_eval.net.current_epoch = 0
      fil = np.array([x for x, y in zip(self.fit_data[2], self.fit_data[3]) if y == 0.0])
      self.filtered_test_data = [cast_float(fil), cast_float(np.zeros((fil.shape[0], 1)))]
      while self.net_eval.net.current_epoch < self.epochs:
        self.fit_data[0:4] = self.net_eval.net.train(*self.fit_data[0:4], 
          epochs=self.epochs, scale=scale, start_epoch=self.net_eval.net.current_epoch, 
          do_epochs=self.net_eval.net.current_epoch + dk)
        erx = np.sqrt(self.net_eval.net.methods['predict_error'](*self.filtered_test_data))
        print ('== filtered permutation error: %13.8f, scaled = %13.8f' % 
          (erx, erx * scale))
        param = self.net_eval.net.get_parameters()
        self.net_ori_eval.net.set_parameters(param)
        self.net_eval.net.saved_parameters = param

        self.ori_fit_data[0:4] = self.net_ori_eval.net.train(*self.ori_fit_data[0:4], 
          epochs=self.epochs, scale=scale, start_epoch=self.net_ori_eval.net.current_epoch, 
          do_epochs=self.net_ori_eval.net.current_epoch + dk)
        param = self.net_ori_eval.net.get_parameters()
        self.net_eval.net.set_parameters(param)
        self.net_ori_eval.net.saved_parameters = param
  
  def test(self):
    yp = self.net_eval.net.predict(self.fit_data[4])
    yp = yp.reshape(yp.shape[0])
    self.fit_data.append(yp)
    if self.diff:
      self.fit_data[5:] = [trans_backward(y, self.emax - self.emin, 0) 
        for y in self.fit_data[5:]]
    else:
      self.fit_data[5:] = [trans_backward(y, self.emax, self.emin) 
        for y in self.fit_data[5:]]
    self.testings = []
    res = 0.0
    rmax = 0.0
    n = len(self.fit_data[5])
    for i, std, pre in zip(range(n), self.fit_data[5], self.fit_data[6]):
      self.testings.append([i, std, pre, np.abs(std - pre)])
      terror = (std - pre) * self.energy_factor
      res += terror ** 2
      if np.abs(terror) > rmax: rmax = np.abs(terror)
    res = np.sqrt(res / n)
    print ('energy rmse = %15.8f, max = %15.8f' % (res, rmax))
    return res

  def test_permutation(self, m=20):
    print ('test permutation ...')
    if self.diff:
      param = self.net_eval.net.get_parameters()
      self.net_ori_eval.net.set_parameters(param)
      net = self.net_ori_eval.net
    else:
      net = self.net_eval.net
    permu, px = trans_data(self.clus, m * len(self.clus), degree=self.degree, 
      d_order=self.second_order, expl=self.exp_length, diff=False, 
      randorder=False, shuffle=True, new_ener=True)
    
    if self.scale_lengths:
      permu = trans_forward(permu, self.lmax, self.lmin)
    permu = permu.reshape((len(self.clus), m) + permu.shape[1:])
    px = px.reshape((len(self.clus), m))
    self.per_testings = []
    res = []
    for i, per in enumerate(permu):
      pred = net.predict(per).reshape(m)
      pxd = px[i]
      pred = trans_backward(pred, self.emax, self.emin)
      self.per_testings.append([i, pxd[0], pred.mean(), 
        np.abs(pred.mean() - pxd[0]) * self.energy_factor, pred.std() * self.energy_factor])
      res.append(pred.std() * self.energy_factor)
    print ('permutation energy rmse = %15.8f, max = %15.8f' % (
      np.array(res).mean(),np.array(res).max()))
    print (' (all data) energy rmse = %15.8f, max = %15.8f' % 
      (np.sqrt((np.array(self.per_testings)[:, 3] ** 2).mean()), 
      np.abs(np.array(self.per_testings)[:, 3]).max()))
    return np.array(res).mean()