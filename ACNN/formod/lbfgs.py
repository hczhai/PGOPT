import cp_lbfgs
import numpy as np

class LBFGS(object):
  class F(object):
    pass
  
  def longstr(self, x):
    return np.array(x + ' ' * (60 - len(x)))
  
  # cannot change float_dtype because it is fixed in fortran code
  def __init__(self, n, float_dtype=np.float64):
    self.float_type = float_dtype
    self.m = 4
    self.n = n
    self.stat = 0
    zb = np.zeros((n, ), dtype=self.float_type)
    zi = np.zeros((n, ), dtype=np.int32)
    self.lb = zb
    self.ub = zb
    self.kb = zi
    self.f = np.array(0.0, dtype=self.float_type)
    self.g = np.zeros((n, ), dtype=self.float_type)
    self.trust_radius = 0.0
    self.max_f_per_iter = 20
    self.iter = 0
    lenwa = 2 * self.m * n + 5 * n + 11 * self.m * self.m + 8 * self.m
    self.wa = np.zeros((lenwa, ), dtype=self.float_type)
    self.iwa = np.zeros((3 * n, ), dtype=np.int32)
    self.iprint = -1
    self.factr = 1E4
    self.lsave = np.array([ False ] * 4)
    self.isave = np.zeros(44, dtype=np.int32)
    self.dsave = np.zeros(29, dtype=self.float_type)
    self.csave = np.array(' ' * 60)
    self.grad_conv = 1E-8
    self.max_iter = 200
    self.log_file = 0
    self.p = self.F()
  
  def start(self, x):
    self.x = np.array(x, dtype=self.float_type)
    self.f = np.array(self.p.eval(self.x), dtype=self.float_type)
    self.traj = []
    self.traj.append([np.array(self.x, dtype=self.float_type), self.f])
    self.task = self.longstr('START')
    self._invoke()
  
  def _invoke(self):
    cp_lbfgs.setulb(self.m, self.x, self.lb, self.ub, self.kb, 
      self.f, self.g, self.factr, self.grad_conv, self.wa, self.iwa, 
      self.task, self.iprint, self.csave, self.lsave, self.isave, 
      self.dsave, self.trust_radius, n=self.n)
  
  def _step(self):
    just_entered = True
    while True:
      if self.task.tostring()[0:7] == 'RESTART':
        self.stat = 0
        self.task = self.longstr('START')
        self._invoke()
      if self.task.tostring()[0:2] == 'FG':
        if self.isave[35] > self.max_f_per_iter:
          self.task = longstr('STOP: CPU, hit max f eval in iter')
          self.stat = 5 # abnormal exit
          self._invoke()
        else:
          self.stat = 1 # evaluate f and g
      elif self.task.tostring()[0:5] == 'NEW_X':
        if just_entered:
          self.stat = 2
          self._invoke()
        else:
          self.stat = 3 # need to update the iter number
      elif self.task.tostring()[0:4] == 'CONV':
        self.stat = 4 # convergence
      elif self.task.tostring()[0:4] == 'STOP':
        self.stat = 5 # abnormal exit
      elif self.task.tostring()[0:5] == 'ERROR':
        self.stat = 5 # abnormal exit
      elif self.task.tostring()[0:8] == 'ABNORMAL':
        self.stat = 5 # abnormal exit
      if self.stat == 3:
        self.iter = self.isave[29]
      if self.stat == 1:
        self.f = np.array(self.p.eval(self.x), dtype=self.float_type)
        self.g = np.array(self.p.evald(self.x), dtype=self.float_type)
        self._invoke()
      elif self.stat >= 3:
        break
      just_entered = False
  
  def opt(self):
    ierr = 0
    if self.log_file != 0:
      print ('   ITER.   ENERGY(OLD)    ENERGY(NEW)      DE' + 
        '          GRADMAX     GRADNORM    GRADRMS     STEPMAX     STEPLEN' + 
        '     STEPRMS')
    xold = np.array(self.x, dtype=self.float_type)
    self.fold = np.array(self.p.eval(self.x), dtype=self.float_type)
    lmr = lambda x: [ np.linalg.norm(x), x.max(), np.std(x) ]
    while True:
      self._step()
      gl, gm, gr = lmr(self.g)
      xil, xim, xir = lmr(self.x - xold)
      if self.log_file != 0:
        print (' %5d%15.8f%15.8f%15.8f%12.6f%12.6f%12.6f%12.8f%12.8f%12.8f' % (
          self.iter, self.fold, self.f, self.f - self.fold, 
          gm, gl, gr, xim, xil, xir))
      if self.stat == 4 or self.stat == 5:
        break
      xold = np.array(self.x, dtype=self.float_type)
      self.fold = np.array(self.f, dtype=self.float_type)
      if self.iter >= self.max_iter:
        self.stat = 5
        break
      self.traj.append([np.array(self.x, dtype=self.float_type), self.f])
    if self.iter >= self.max_iter:
      if self.log_file != 0: print ('Maximum number of iteration reached.')
    if self.stat == 5: ierr = 1
    if self.log_file != 0:
      print ('Optimized variables')
      print (self.x)
    return ierr
