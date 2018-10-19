#!/bin/bash

cat << EOF | python > bc.out
import ase
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.startgenerator import StartGenerator
from ase.data import covalent_radii
import numpy as np, random
from ase import Atoms

seed = @GRSEED
np.random.seed(seed)
random.seed(seed)

class ASECluster(object):
    @staticmethod
    def generator(name):
        c = Atoms(name).get_atomic_numbers()
        slab = fcc111('Au', size=(16, 16, 1), vacuum=10.0, orthogonal=True)
        slab.set_constraint(FixAtoms(mask=len(slab) * [True]))
        pos = slab.get_positions()
        cell = slab.get_cell()
        p0 = np.array([0., 0., max(pos[:, 2]) + 5.])
        v1 = np.array([0.0] * 3)
        v2 = np.array([0.0] * 3)
        v3 = np.array([0.0] * 3)
        v1[0] = 10.
        v2[1] = 10.
        v3[2] = 10.
        natom = len(c)
        atom_numbers = c
        unique_atom_types = get_all_atom_types(slab, atom_numbers)
        print 'cov: ', covalent_radii[atom_numbers[0]]
        cd = closest_distances_generator(atom_numbers=unique_atom_types,
                                         ratio_of_covalent_radii=0.8)
        sg = StartGenerator(slab=slab,
                            atom_numbers=atom_numbers,
                            closest_allowed_distances=cd,
                            box_to_place_in=[p0, [v1, v2, v3]])
        while True:
            a = sg.get_new_candidate().get_positions()[-natom:]
            a = a - a.mean(axis=0).reshape(1, 3)
            yield a, Atoms(name).get_chemical_symbols()

def write_clusters(fn, name, number):
    gen = ASECluster.generator(name)
    for i, (p, e) in enumerate(gen):
        if i % 10 == 0: print (i)
        if i == number: break
        if isinstance(fn, str): f = open(fn, 'a')
        else: f = fn
        f.write('{}\nSTR:{} = 0.0\n'.format(len(p), i))
        for j in range(len(p)):
            f.write('%5s%15.6f%15.6f%15.6f\n' % (e[j], p[j, 0], p[j, 1], p[j, 2]))
        if isinstance(fn, str): f.close()

write_clusters('./bc_structs.xyz', '@GNAME', @GNUMBER)
EOF
