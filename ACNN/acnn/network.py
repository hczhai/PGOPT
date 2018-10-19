
import time, numpy as np
import theano, theano.tensor as T
from utils.base import cast_float

class BasicNetwork(object):
  def __init__(self, layers=[], shuffle=False, shuffle_per=0, show_epoch=1, 
    param_save=0, **opts):
    for l in range(len(layers))[1:]:
      if layers[l].calc_size:
        layers[l].num_input_size = layers[l - 1].num_output_size
        layers[l].num_output_size = layers[l].size_output(layers[l].num_input_size)
      else:
        if layers[l].num_input_size != layers[l - 1].num_output_size:
          raise ValueError('layer size mismatch! {} > {}'.format(layers[l - 1], layers[l]) 
            + '\n' + '{} > {}'.format(layers[l - 1].num_output_size, layers[l].num_input_size))
    self.layers = layers
    self.shuffle = shuffle
    self.shuffle_per = shuffle_per
    self.show_epoch = show_epoch
    self.train_errors = []
    self.test_errors = []
    self.train_times = []
    self.input_size = tuple(self.layers[0].num_input_size)
    self.output_size = tuple(self.layers[-1].num_output_size)
    self.saved_parameters = None
    self.param_save = param_save
    self.current_epoch = 0
  
  def __repr__(self):
    return ' > '.join([l.__repr__() for l in self.layers])
  
  def train(self, input_train, output_train, input_test, output_test, 
    epochs=100, scale=None, start_epoch=0, do_epochs=None):
    if do_epochs is None or do_epochs > epochs: do_epochs = epochs
    holder = 52 if scale is None else 68
    train_ind = np.arange(input_train.shape[0])
    if start_epoch == 0:
      self.train_errors = []
      self.test_errors = []
      self.train_times = []
      input_train = cast_float(input_train)
      output_train = cast_float(output_train)
      input_test = cast_float(input_test)
      output_test = cast_float(output_test)
      if len(output_test.shape) == 1:
        output_test = output_test.reshape((output_test.shape[0], 1))
      if len(output_train.shape) == 1:
        output_train = output_train.reshape((output_train.shape[0], 1))
      if (input_train.shape[0] != output_train.shape[0] or 
        input_test.shape[0] != output_test.shape[0]):
        raise RuntimeError('train/test data first dimension mismatch!')
      if (input_train.shape[1:] != self.input_size or 
        input_test.shape[1:] != self.input_size or 
        output_train.shape[1:] != self.output_size or 
        output_test.shape[1:] != self.output_size):
        raise RuntimeError('train/test data size is not the size of network!\n' + 
          'input shape: {}, net input shape: {}\n'.format(input_train.shape[1:],self.input_size) +
          'ouput shape: {}, net ouput shape: {}'.format(output_train.shape[1:],self.output_size))
      print ('pre-training working ...')
      self.pre_train(input_train, output_train)
      self.pre_test(input_test, output_test)
      print (('%%%ds' % (holder)) % ('-'*holder, ))
      if scale is None:
        print ('%10s %15s %15s %7s  ' % ('epoch #', 'train error', 'test error', 'time'))
      else:
        print ('%10s %15s %15s %15s %7s  ' % ('epoch #', 'train error', 
          'test error', 'scaled error', 'time'))
      self.min_error = None
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    start_time_each = start_time = init_time = time.time()
    for epoch in xrange(start_epoch, do_epochs):
      self.current_epoch = epoch
      if self.shuffle_per != 0:
        if epoch % self.shuffle_per == 0:
          self.shuffle = True
        else:
          self.shuffle = False
      if self.shuffle:
        np.random.shuffle(train_ind)
        input_train = input_train[train_ind]
        output_train = output_train[train_ind]
      self.pre_epoch(input_train, output_train)
      mid_time = time.time()
      train_error = self.train_epoch(input_train, output_train)
      test_error = self.test_epoch(input_test, output_test)
      finish_time = time.time()
      train_time_each = finish_time - start_time_each
      start_time_each = finish_time
      self.train_errors.append(train_error)
      self.test_errors.append(test_error)
      self.train_times.append(train_time_each)
      is_min = False
      if self.param_save != 0 and epoch % self.param_save == 0:
        if self.min_error is None or test_error < self.min_error:
          self.min_error = test_error
          self.saved_parameters = self.get_parameters()
          is_min = True
      if epoch == 0 or (epoch + 1) % self.show_epoch == 0:
        train_time = finish_time - start_time
        start_time = finish_time
        if scale is None:
          print ('%10d %15.8f %15.8f %7.2f%s' % (epoch + 1, train_error, 
            test_error, train_time, ' *' if is_min else '  '))
        else:
          print ('%10d %15.8f %15.8f %15.8f %7.2f%s' % (epoch + 1, train_error, test_error, 
            np.sqrt(test_error) * scale, train_time, ' *' if is_min else '  '))
        if epoch == 0:
          if self.shuffle_per == 0:
            total_time = train_time * epochs
          else:
            shu_epochs = epochs / self.shuffle_per
            total_time = (finish_time - mid_time) * (epochs - shu_epochs) + train_time * shu_epochs
          print (('%%%ds' % (holder)) % ('estimated total time %d:%02d:%02d' % (total_time / 3600, 
            (total_time / 60) % 60, total_time % 60), ))
          print (('%%%ds' % (holder)) % ('-'*holder, ))
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    total_time = finish_time - init_time
    print (('%%%ds' % (holder)) % ('total time used: %d:%02d:%02d' % (total_time / 3600, 
      (total_time / 60) % 60, total_time % 60), ))
    print (('%%%ds' % (holder)) % ('-'*holder, ))
    self.current_epoch = do_epochs
    self.post_train()
    return input_train, output_train, input_test, output_test
  
  def pre_train(self, input_train, output_train):
    pass
  
  def pre_test(self, input_test, output_test):
    pass
  
  def pre_epoch(self, input_train, output_train):
    pass
  
  def post_train(self):
    pass
  
  def get_parameters(self):
    return [[p.get_value() for p in layer.parameters] for layer in self.layers]
  
  def set_parameters(self, params):
    for i, layer in enumerate(params):
      for j, p in enumerate(layer):
        if j == 0: self.layers[i].weight.set_value(p)
        elif j == 1: self.layers[i].bias.set_value(p)

class LearningNetwork(BasicNetwork):
  def __init__(self, auto=True, max_lambda=0, **opts):
    self.auto = auto
    self.max_lambda = max_lambda
    super(LearningNetwork, self).__init__(**opts)
    self.init_layers()
    self.init_variables()
    self.init_methods()
  
  def init_layers(self):
    print ('initialize layers ...')
    for layer in self.layers:
      layer.initialize()
  
  def init_variables(self):
    print ('initialize variables ...')
    dtype = theano.config.floatX
    net_input = T.TensorType(dtype, (False, ) * (len(self.input_size) + 1))('x')
    net_output = T.TensorType(dtype, (False, ) * (len(self.output_size) + 1))('y')
    prediction = net_input
    for layer in self.layers:
      prediction = layer.output(prediction)
    if self.max_lambda == 0:
      net_error = T.square(net_output - prediction).mean()
    else:
      net_error = (T.square(net_output - prediction).mean() + 
        self.max_lambda * T.square(net_output - prediction).max())
    self.variables = { 'net_input': net_input, 'net_output': net_output, 
      'prediction': prediction, 'error': net_error }
  
  def init_methods(self):
    print ('initialize methods ...')
    net_input = self.variables['net_input']
    net_output = self.variables['net_output']
    prediction = self.variables['prediction']
    error = self.variables['error']
    updates = self.init_updates()
    fpredict = theano.function([net_input], prediction)
    ftrain_epoch = theano.function([net_input, net_output], error, updates=updates)
    fpredict_error = theano.function([net_input, net_output], error)
    self.methods = { 'predict': fpredict, 'predict_error': fpredict_error, 
      'train_epoch': ftrain_epoch, 'updates': updates }
  
  def init_updates(self):
    updates = []
    for layer in self.layers: updates += self.init_layer_updates(layer)
    return updates
  
  def init_layer_updates(self, layer):
    updates = []
    for param in layer.parameters: updates += self.init_param_updates(layer, param)
    return updates
  
  def init_param_updates(self, layer, param):
    return []
  
  def predict(self, input_data):
    return self.methods['predict'](cast_float(input_data))
  
  def test_epoch(self, input_data, output_data):
    if self.auto:
      return self.methods['auto_test']()
    else:
      return self.methods['predict_error'](input_data, output_data)
  
  def train_epoch(self, input_data, output_data):
    if self.auto:
      return self.methods['auto_train']()
    else:
      return self.methods['train_epoch'](input_data, output_data)
  
  @staticmethod
  def create_shared(name, shape):
    value = cast_float(np.zeros(shape))
    return theano.shared(value=value, name=name, borrow=True)

  def pre_train(self, input_train, output_train):
    if self.auto:
      xerr = self.variables['error']
      xip = self.variables['net_input']
      xop = self.variables['net_output']
      xipt = self.variables['input_train'] = LearningNetwork.create_shared('input_train', 
        input_train.shape)
      xopt = self.variables['output_train'] = LearningNetwork.create_shared('output_train', 
        output_train.shape)
      self.variables['input_train'].set_value(input_train)
      self.variables['output_train'].set_value(output_train)
      xupd = self.methods['updates']
      self.methods['auto_train'] = theano.function(inputs=[], outputs=xerr, updates=xupd, 
        givens={xip: xipt, xop: xopt})
  
  def pre_test(self, input_test, output_test):
    if self.auto:
      xerr = self.variables['error']
      xip = self.variables['net_input']
      xop = self.variables['net_output']
      self.methods['auto_test'] = theano.function(inputs=[], outputs=xerr, 
        givens={xip: input_test, xop: output_test})
  
  def pre_epoch(self, input_train, output_train):
    if self.auto and self.shuffle:
      self.variables['input_train'].set_value(input_train)
      self.variables['output_train'].set_value(output_train)
  
  def release(self):
    if self.auto and self.shuffle:
      self.variables['input_train'].set_value([[]])
      self.variables['output_train'].set_value([[]])
    