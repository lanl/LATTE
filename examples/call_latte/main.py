import numpy as np
from sedacs.system import *
from sedacs.system import get_hindex, System
from sedacs.periodic_table import PeriodicTable
from sedacs.file_io import *
from latte_lib import latte_compute, latte_compute_hs #TODO

coords_file = "coords.xyz"
latticeVectors, symbols, atomTypes, coords0 = read_xyz_file(coords_file, lib="None")

field = np.zeros(3)
verb = True

#err, charges_out, forces_out, dipole_out, energy_out = latte_compute(latticeVectors, symbols, atomTypes, coords0, field, verb)

# Get hindex (the orbital index for each atom in the system)
orbs = [1, 4]
sy = System()
sy.nats = len(atomTypes)

sy.norbs, sy.orbs, hindex, sy.numel, sy.znuc = get_hindex(orbs, symbols, atomTypes)

err, ham, over = latte_compute_hs(sy.norbs, latticeVectors, symbols, atomTypes, coords0, field, verb)

print('ham', ham)


#print('force', forces_out)

#print(symbols)

#print('err', err)
