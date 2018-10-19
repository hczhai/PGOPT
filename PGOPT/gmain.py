#!/usr/bin/env python

# Parallel Global Optimization scripts part (GPU)
# this script is used to generate submit scripts
# in differnet platforms
# it must be used with ACNN package
# ACNN is platform independent while this script
# handles different DFT package and supercomputer systems

import os
from sys import argv
from gpu.driver import GPUDriver

if __name__ != "__main__":
  raise RuntimeError('Must run this file as the main program!')
if "PGOPTHOME" not in os.environ:
  raise RuntimeError("must set 'PGOPTHOME' environment variable!")

gpd = GPUDriver(argv)
gpd.run()
