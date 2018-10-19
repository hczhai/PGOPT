
import theano
import theano.tensor as T
import numpy as np
from utils.base import cast_float

# abstract layer size
class LayerSize(object):
  def __init__(self, n):
    self.n = n
  def __repr__(self):
    return '{}({})'.format(self.__class__.__name__[0], self.n)

# basic input/output size
class BasicLayerSize(LayerSize):
  def __init__(self, n):
    super(BasicLayerSize, self).__init__(n)

# input size sharing common parameters
class CoLayerSize(LayerSize):
  def __init__(self, n):
    super(CoLayerSize, self).__init__(n)

# output size creating feature maps
class FeatureLayerSize(LayerSize):
  def __init__(self, n):
    super(FeatureLayerSize, self).__init__(n)

# mixed input/output size
# different feature sharing different common parameters
# output must have the same size as input
class MixedLayerSize(LayerSize):
  # n is the input/output size
  # m is the weight size
  def __init__(self, n, m):
    super(MixedLayerSize, self).__init__(n)

# basic neural network layer
class BasicLayer(object):
  base_id = 0
  def __init__(self):
    self.parameters = []
    self.calc_size = True
    BasicLayer.base_id += 1
    self.layerid = BasicLayer.base_id
  
  def initialize(self):
    pass

  def __repr__(self):
    return '{name}()'.format(name=self.__class__.__name__)

# layer with weight and bias
class ParameterLayer(BasicLayer):

  def __init__(self, input_size, output_size, mix_indices=None):
    super(ParameterLayer, self).__init__()
    self.mix_indices = mix_indices
    self.init_size(input_size, output_size)
  
  def init_size(self, input_size, output_size):
    if isinstance(input_size, int): input_size = [input_size]
    if isinstance(output_size, int): output_size = [output_size]
    self.input_size = [BasicLayerSize(i) if isinstance(i, int) else i for i in input_size]
    self.output_size = [BasicLayerSize(i) if isinstance(i, int) else i for i in output_size]
    self.weight = None
    self.bias = None
    self._weight_bias_shape()
    self.num_input_size = [n.n for n in self.input_size]
    self.num_output_size = [n.n for n in self.output_size]
    self.calc_size = False
  
  @staticmethod
  def xavier_normal(shape, isbias):
    if isbias: fanin, fanout = 1, shape[-1]
    else: fanin, fanout = shape[-2], shape[-1]
    std = np.sqrt(2.0 / (fanin + fanout))
    return np.random.normal(loc=0.0, scale=std, size=shape)

  @staticmethod
  def create_shared(value, name, shape, isbias):
    if isinstance(value, (T.sharedvar.SharedVariable, T.Variable)):
      return value
    if value is None:
      value = ParameterLayer.xavier_normal(shape, isbias)
    value = cast_float(value)
    return theano.shared(value=value, name=name, borrow=True)

  def _weight_bias_shape(self):
    # input order: (co, mix, basici) disordered
    # weight order: (mix,) feature, basici, basico
    # direct output order: (mix,) feature, basico, batch, co
    # bias order: (mix,) feature, basico
    # output order: batch, (feature, co, mix, basico) disordered
    ip_c, ip_m, ip_bi = [], [], [] # indices
    ip_cs, ip_ms = [], [] # sizes
    w_m, w_f, w_bi, w_bo = [], [], [], [] # sizes
    b_m = [] # sizes
    op_f, op_c, op_m, op_bo = [], [], [], [] # indices
    op_cs, op_ms = [], [] # sizes
    for i in range(len(self.input_size)):
      li = self.input_size[i]
      if isinstance(li, BasicLayerSize):
        ip_bi.append(i)
        w_bi.append(li.n)
      elif isinstance(li, CoLayerSize):
        ip_c.append(i)
        ip_cs.append(li.n)
      elif isinstance(li, MixedLayerSize):
        ip_m.append(i)
        ip_ms.append(li.n)
        w_m.append(li.m)
        b_m.append(li.n)
      else:
        raise ValueError('invalid input size! (layer #{}: {})'.format(i, li))
    for i in range(len(self.output_size)):
      li = self.output_size[i]
      if isinstance(li, BasicLayerSize):
        op_bo.append(i)
        w_bo.append(li.n)
      elif isinstance(li, CoLayerSize):
        op_c.append(i)
        op_cs.append(li.n)
      elif isinstance(li, MixedLayerSize):
        op_m.append(i)
        op_ms.append(li.n)
      elif isinstance(li, FeatureLayerSize):
        op_f.append(i)
        w_f.append(li.n)
      else:
        raise ValueError('invalid output size! (layer #{}: {})'.format(i, li))
    if ip_ms != op_ms or ip_cs != op_cs:
      raise ValueError('co/mixed layer size mismatch!')
    if len(w_m) > 1:
      raise ValueError('only support one mixed layer size!')
    if len(ip_bi) != 1 or len(op_bo) != 1:
      raise ValueError('must have one input/output basic layer size!')
    self.weight_shape = w_m + w_f + w_bi + w_bo
    self.bias_shape = b_m + w_f + w_bo
    self.bias_reshape = self.bias_shape + [1] * (len(ip_c) + 1)
    self.op_shuffle = [-1] * (len(self.output_size) + 1)
    ix = 0
    for k in (op_m, op_f, op_bo, [-1], op_c):
      for i in k:
        self.op_shuffle[i + 1] = ix
        ix += 1
    self.mixed_pair = [0, ip_m[0] + 1] if len(w_m) == 1 else None
    self.basic_pair = [len(w_m) + len(w_f), ip_bi[0] + 1]

  def initialize(self):
    super(ParameterLayer, self).initialize()
    self.weight = ParameterLayer.create_shared(value=self.weight, 
      name='weight_{}'.format(self.layerid), 
      shape=self.weight_shape, isbias=False)
    self.bias = ParameterLayer.create_shared(value=self.bias, 
      name='bias_{}'.format(self.layerid), 
      shape=self.bias_shape, isbias=True)
    self.parameters = [self.weight, self.bias]
  
  def __repr__(self):
    return '{name}(input = {ip}, output = {op})'.format(name=self.__class__.__name__, 
      ip=self.input_size, op=self.output_size)

class ActivationLayer(ParameterLayer):
  activations = []
  def __init__(self, *opts):
    self.activations = self.__class__.activations
    super(ActivationLayer, self).__init__(*opts)
  
  def initialize(self):
    super(ActivationLayer, self).initialize()
  
  def output(self, input_value):
    if self.mixed_pair != None:
      sn = slice(None)
      idw = [sn] * len(self.weight_shape)
      idip = [sn] * (len(self.input_size) + 1)
      opv = []
      for i, j in self.mix_indices:
        idw[self.mixed_pair[0]], idip[self.mixed_pair[1]] = i, j
        opvd = T.tensordot(self.weight[idw], input_value[idip], axes = self.basic_pair)
        opv.append(opvd.dimshuffle(['x'] + range(opvd.ndim)))
        idw[self.mixed_pair[0]], idip[self.mixed_pair[1]] = sn, sn
      input_value = T.concatenate(opv, axis=0)
    else:
      input_value = T.tensordot(self.weight, input_value, axes = self.basic_pair)
    input_value = input_value + self.bias.reshape(self.bias_reshape)
    input_value = input_value.dimshuffle(self.op_shuffle)
    for act in self.activations:
      input_value = act(input_value)
    return input_value

class Linear(ActivationLayer):
  activations = [lambda x: x]

class Sigmoid(ActivationLayer):
  activations = [T.nnet.sigmoid]

class Tanh(ActivationLayer):
  activations = [T.tanh]

class Relu(ActivationLayer):
  activations = [T.nnet.relu]

class Softplus(ActivationLayer):
  activations = [T.nnet.softplus]

class TanhSig(ActivationLayer):
  activations = [lambda x: T.tanh(x) + T.nnet.sigmoid(x)]

class Softmax(ActivationLayer):
  activations = [T.nnet.softmax]

class Reshape(BasicLayer):
  def __init__(self, shape=None):
    super(Reshape, self).__init__()
    self.shape = shape
  
  def output(self, input_value):
    shape = self.shape if self.shape != None else (T.prod(input_value.shape[1:]), )
    return T.reshape(input_value, (input_value.shape[0], ) + shape)
  
  def size_output(self, input_size):
    shape = self.shape if self.shape != None else (np.prod(input_size), )
    if shape[-1] == -1:
      shape = shape[:-1] + (np.prod(input_size) / np.prod(shape[:-1]), )
    return list(shape)
  
  def __repr__(self):
    return '{name}(shape = {sh})'.format(name=self.__class__.__name__, sh=self.shape)

class AxisPoolLayer(BasicLayer):
  methods = []
  def __init__(self, axis=None):
    self.methods = self.__class__.methods
    super(AxisPoolLayer, self).__init__()
    self.axis = axis
  
  def output(self, input_value):
    for method in self.methods:
      input_value = method(input_value)(axis=self.axis + 1)
    if input_value.ndim == 1:
      input_value = T.reshape(input_value, (input_value.shape[0], 1))
    return input_value
    
  def size_output(self, input_size):
    x = input_size[self.axis + 1:]
    if len(x) == 0: x = [1]
    return input_size[:self.axis] + x
  
  def __repr__(self):
    return '{name}(axis = {ax})'.format(name=self.__class__.__name__, ax=self.axis)

class Average(AxisPoolLayer):
  methods = [lambda x: x.mean]

class Maximum(AxisPoolLayer):
  methods = [lambda x: x.max]

class Minimum(AxisPoolLayer):
  methods = [lambda x: x.min]

class Diff(AxisPoolLayer):
  methods = [lambda x: (lambda axis, x=x: x[:, 1] - x[:, 0])]