
import theano
import theano.tensor as T
import math
import numpy as np
from acnn.network import LearningNetwork
from utils.base import cast_float

class GradientDescentNet(LearningNetwork):
  def __init__(self, step=0.1, step_decay=False, step_decay_factor=0, **opts):
    self.step = step
    self.step_decay = step_decay
    self.step_decay_factor = step_decay_factor
    super(GradientDescentNet, self).__init__(**opts)
  
  def init_param_updates(self, layer, param):
    grad = T.grad(self.variables['error'], wrt=param)
    if self.step_decay:
      this_step = (self.step * self.step_decay_factor / 
        (self.step_decay_factor + self.current_epoch))
    else:
      this_step = self.step
    return [(param, param - this_step * grad)]
  
  def post_train(self):
    if self.step_decay:
      this_step = (self.step * self.step_decay_factor / 
        (self.step_decay_factor + self.current_epoch))
      print ('final step = {}'.format(this_step))

class MinibatchGradientDescentNet(GradientDescentNet):
  def __init__(self, batch_size=1000, **opts):
    self.batch_size = batch_size
    super(MinibatchGradientDescentNet, self).__init__(**opts)
  
  def init_variables(self):
    super(MinibatchGradientDescentNet, self).init_variables()
    index = T.lscalar('index')
    self.variables['index'] = index
  
  def pre_train(self, input_train, output_train):
    if self.auto:
      batch_size = self.batch_size
      xerr = self.variables['error']
      xip = self.variables['net_input']
      xop = self.variables['net_output']
      xidx = self.variables['index']
      xipt = self.variables['input_train'] = LearningNetwork.create_shared('input_train', 
        input_train.shape)
      xopt = self.variables['output_train'] = LearningNetwork.create_shared('output_train', 
        output_train.shape)
      self.variables['input_train'].set_value(input_train)
      self.variables['output_train'].set_value(output_train)
      xupd = self.methods['updates']
      self.methods['auto_train'] = theano.function(inputs=[xidx], outputs=xerr, updates=xupd, 
        givens={xip: xipt[xidx * batch_size: (xidx + 1) * batch_size], 
        xop: xopt[xidx * batch_size: (xidx + 1) * batch_size]})
  
  def train_epoch(self, input_data, output_data):
    if self.batch_size is None or input_data.shape[0] <= self.batch_size:
      return super(MinibatchGradientDescentNet, self).train_epoch(input_data, output_data)
    batch_size = self.batch_size
    n_samples = input_data.shape[0]
    n_batches = n_samples / batch_size
    if n_batches * batch_size < n_samples: n_batches += 1
    last_size = n_samples % batch_size
    errors = []
    for batch_index in xrange(n_batches):
      if not self.auto:
        idx = slice(batch_index * batch_size, (batch_index + 1) * batch_size)
        error = super(MinibatchGradientDescentNet, self).train_epoch(input_data[idx], output_data[idx])
      else:
        error = self.methods['auto_train'](batch_index)
      errors.append(error)
    if last_size == 0:
      error = np.mean(errors)
    else:
      error = (np.sum(errors[:-1]) * batch_size + errors[-1] * last_size) / n_samples
    return error
  
  def predict(self, input_data):
    batch_size = self.batch_size
    n_samples = input_data.shape[0]
    n_batches = n_samples / batch_size
    if n_batches * batch_size < n_samples: n_batches += 1
    pred = cast_float(np.zeros((0, 1)))
    input_data = cast_float(input_data)
    for batch_index in xrange(n_batches):
      idx = slice(batch_index * batch_size, (batch_index + 1) * batch_size)
      predx = self.methods['predict'](input_data[idx])
      pred = np.concatenate([pred, predx])
    return pred

class MomentumNet(MinibatchGradientDescentNet):
  def __init__(self, momentum=0.0, nesterov=False, **opts):
    self.momentum = momentum
    self.nesterov = nesterov
    super(MomentumNet, self).__init__(**opts)

  def init_layers(self):
    super(MomentumNet, self).init_layers()
    for layer in self.layers:
      for parameter in layer.parameters:
        value = cast_float(np.zeros(T.shape(parameter).eval()))
        parameter.delta = theano.shared(name="delta_" + parameter.name, value=value)
  
  def init_param_updates(self, layer, param):
    step = self.step
    grad = T.grad(self.variables['error'], wrt=param)
    prev_delta = param.delta
    delta = self.momentum * prev_delta - step * grad
    if self.nesterov: delta = self.momentum * delta - step * grad
    return [(param, param + delta), (prev_delta, delta)]
