
from acnn.network import LearningNetwork
from utils.base import cast_float
from theano.ifelse import ifelse
import theano
import theano.tensor as T
from itertools import chain
import numpy as np
from formod.lbfgs import LBFGS

def iter_parameters(network):
  parameters = [layer.parameters for layer in network.layers]
  return chain(*parameters)

def parameters2vector(network):
  params = iter_parameters(network)
  return T.concatenate([param.flatten() for param in params])

class LBFGSNet(LearningNetwork):
  def __init__(self, **opts):
    super(LBFGSNet, self).__init__(**opts)
  
  def init_methods(self):
    super(LBFGSNet, self).init_methods()
    params = list(iter_parameters(self))
    xerr = self.variables['error']
    net_input = self.variables['net_input']
    net_output = self.variables['net_output']
    grad = T.grad(xerr, wrt=params)
    fgrad = T.concatenate([g.flatten() for g in grad])
    self.variables['grad'] = fgrad
  
  def pre_train(self, input_train, output_train):
    self.auto = True
    super(LBFGSNet, self).pre_train(input_train, output_train)
    xip = self.variables['net_input']
    xop = self.variables['net_output']
    xipt = self.variables['input_train']
    xopt = self.variables['output_train']
    fgrad = self.variables['grad']
    self.methods['auto_grad'] = theano.function(inputs=[], outputs=fgrad, 
      givens={xip: xipt, xop: xopt})
    x = parameters2vector(self).eval()
    self.task = LBFGS(x.shape[0])
    self.task.log_file = 0
    def opt_eval(x):
      set_parameter_value(self, x)
      return self.methods['auto_train']()
    self.task.p.eval = self.opt_eval
    self.task.p.evald = self.opt_evald
    self.task.log_file = 0
    self.task.start(x)

  def opt_eval(self, x):
    self.set_parameter_value(x)
    return self.methods['auto_train']()

  def opt_evald(self, x):
    self.set_parameter_value(x)
    return self.methods['auto_grad']()
  
  def train_epoch(self, input_data, output_data):
    self.task._step()
    return self.task.f
  
  def set_parameter_value(self, vector):
    vector = cast_float(vector)
    parameters = list(iter_parameters(self))
    start_position = 0
    for p in parameters:
      end_position = start_position + int(p.size.eval())
      pp = np.reshape(vector[start_position:end_position], p.shape.eval())
      p.set_value(pp)
      start_position = end_position
  
  def pre_epoch(self, input_train, output_train):
    pass

def setup_parameter_updates(parameters, vector):
  updates = []
  start_position = 0
  for p in parameters:
    end_position = start_position + p.size
    pp = T.reshape(vector[start_position:end_position], p.shape)
    updates.append((p, pp))
    start_position = end_position
  return updates

def compute_jaccobian(errors, parameters):
  n_samples = errors.shape[0]
  def find_jacobbian(i, errors, *params):
    return T.grad(T.sum(errors[i]), wrt=params)
  J, _ = theano.scan(find_jacobbian, sequences=T.arange(n_samples),
    non_sequences=[errors] + parameters)
  jaccobians = []
  for jaccobian, parameter in zip(J, parameters):
    jaccobian = jaccobian.reshape((n_samples, parameter.size))
    jaccobians.append(jaccobian)
  return T.concatenate(jaccobians, axis=1)

class LevenbergMarquardtNet(LearningNetwork):
  def __init__(self, mu=0.01, mu_update_factor=5, **opts):
    self.mu = mu
    self.mu_fac = mu_update_factor
    super(LevenbergMarquardtNet, self).__init__(**opts)
  
  def init_variables(self):
    super(LevenbergMarquardtNet, self).init_variables()
    self.variables['mu'] = theano.shared(name="mu", value=cast_float(self.mu))
    self.variables['last_error'] = theano.shared(name="last_error", 
      value=cast_float(np.nan))

  def init_updates(self):
    xop = self.variables['net_output']
    xpre = self.variables['prediction']
    xerr = self.variables['error']
    xlerr = self.variables['last_error']
    xmu = self.variables['mu']

    new_mu = ifelse(T.lt(xlerr, xerr), xmu * self.mu_fac, xmu / self.mu_fac)
    mse_for_each_sample = T.mean((xop - xpre) ** 2, axis=1)
    params = list(iter_parameters(self))
    param_vector = parameters2vector(self)
    J = compute_jaccobian(mse_for_each_sample, params)
    n_params = J.shape[1]
    updated_params = param_vector - T.nlinalg.matrix_inverse(
      J.T.dot(J) + new_mu * T.eye(n_params)).dot(J.T).dot(mse_for_each_sample)
    updates = [(xmu, new_mu)]
    parameter_updates = setup_parameter_updates(params, updated_params)
    updates += parameter_updates
    return updates

  def pre_epoch(self, input_train, output_train):
    super(LevenbergMarquardtNet, self).pre_epoch(input_train, output_train)
    last_error = self.train_errors[-1] if len(self.train_errors) != 0 else None
    if last_error is not None:
      self.variables['last_error'].set_value(last_error)


def bfgs(inverse_hessian, weight_delta, gradient_delta, maxrho=1e4):
  ident_matrix = cast_float(T.eye(inverse_hessian.shape[0]))
  maxrho = cast_float(maxrho)
  rho = cast_float(1.) / gradient_delta.dot(weight_delta)
  rho = ifelse(T.isinf(rho), maxrho * T.sgn(rho), rho)
  param1 = ident_matrix - T.outer(weight_delta, gradient_delta) * rho
  param2 = ident_matrix - T.outer(gradient_delta, weight_delta) * rho
  param3 = rho * T.outer(weight_delta, weight_delta)
  return param1.dot(inverse_hessian).dot(param2) + param3

def dfp(inverse_hessian, weight_delta, gradient_delta, maxnum=1e5):
  maxnum = cast_float(maxnum)
  quasi_dot_gradient = inverse_hessian.dot(gradient_delta)
  param1 = T.outer(weight_delta, weight_delta) / T.dot(gradient_delta, weight_delta)
  param2_numerator = T.clip(T.outer(quasi_dot_gradient, gradient_delta) * inverse_hessian,
    -maxnum, maxnum)
  param2_denominator = gradient_delta.dot(quasi_dot_gradient)
  param2 = param2_numerator / param2_denominator
  return inverse_hessian + param1 - param2

def psb(inverse_hessian, weight_delta, gradient_delta, **options):
    gradient_delta_t = gradient_delta.T
    param = weight_delta - inverse_hessian.dot(gradient_delta)
    devider = (1. / T.dot(gradient_delta, gradient_delta))
    param1 = T.outer(param, gradient_delta) + T.outer(gradient_delta, param)
    param2 = T.dot(gradient_delta, param) * T.outer(gradient_delta, gradient_delta_t)
    return inverse_hessian + param1 * devider - param2 * devider ** 2

def sr1(inverse_hessian, weight_delta, gradient_delta, epsilon=1e-8):
    epsilon = asfloat(epsilon)
    param = weight_delta - inverse_hessian.dot(gradient_delta)
    denominator = T.dot(param, gradient_delta)
    return ifelse(T.lt(T.abs_(denominator),
      epsilon * param.norm(L=2) * gradient_delta.norm(L=2)),
      inverse_hessian, inverse_hessian + T.outer(param, param) / denominator)
