#!/usr/bin/env python

# Parallel Global Optimization platform-independent part
# use PGOPT to generate submission scripts
# ACNN: Atomic Convolutional Neural Network
# theano, reportlab are optional packages
# go to formod subdirectory to compile fortran part
# L-BFGS and DFS-AM part is written in fortran

from __future__ import print_function
from sys import platform as _platform
import numpy as np, sys, os, time
sys.path.insert(0, sys.path[0] + '/..')

class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)

from utils.io import read_json, write_json
from utils.io import load_data, dump_data, new_file_name, load_name

if __name__ != "__main__":
    raise RuntimeError('Must run this file as the main program!')

if len(sys.argv) < 1 + 1:
    raise RuntimeError('Need input file!')

ip = read_json(sys.argv[1])

theano_flags = []
if _platform == 'darwin':
    theano_flags.append("cxx=/usr/local/bin/g++-5")
if 'OMP_NUM_THREADS' in os.environ:
    print ('number of openmp threads: {}'.format(os.environ['OMP_NUM_THREADS']))
    theano_flags.append("openmp=True")
    theano_flags.append("openmp_elemwise_minsize=100000")
if 'gpu' in ip and ip['gpu'] is not False:
    theano_flags.append("device=" + ip['gpu'])
    theano_flags.append("lib.cnmem=1")
    theano_flags.append("floatX=float32")
if len(theano_flags) != 0:
    if 'THEANO_FLAGS' in os.environ:
        os.environ['THEANO_FLAGS'] += ',' + ','.join(theano_flags)
    else:
        os.environ['THEANO_FLAGS'] = ','.join(theano_flags)

from utils.base import reproducible
reproducible(seed=0)
sys.setrecursionlimit(100000)

from utils.driver import Driver
Driver(ip).run_all()
