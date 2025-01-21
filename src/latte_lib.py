#!/usr/bin/env python3

import numpy as np
import ctypes as ct
import os
from periodic_table import PeriodicTable


def latte_compute_hs(norbs, latticeVectors, symbols, atomTypes, coords, field, verb):

    print(latticeVectors)
    # Import the shared library
    latteLibFileName = os.environ["LATTE_PATH"] + "/liblatte.so"

    latteLib = ct.CDLL(latteLibFileName)
    latte_compute_hs_f = latteLib.latte_compute_hs
    # Periodic table: We use this to pass the chemical atom types as integer instead of characters.
    pt = PeriodicTable()

    nats = len(coords[:, 1])

    # Getting atomic numbers
    nTypes = len(symbols)
    atomicNumbers = np.zeros((nTypes), dtype=np.int32)
    atomTypes32 = np.zeros((nats), dtype=np.int32)
    atomTypes32[:] = atomTypes
    for i in range(len(symbols)):
        atomicNumbers[i] = pt.get_atomic_number(symbols[i])

    print(atomicNumbers)

    # Vectorizing 2D arrays for C-Fortran interoperability
    coordsFlat_in = np.zeros(3 * nats)  # Vectorized coordinates
    hamFlat_out = np.zeros(norbs * norbs)  # Vectorized hamiltonian
    overFlat_out = np.zeros(norbs * norbs)  # Vectorized overlap
    fieldFlat_in = np.zeros(3)  # Same here

    for i in range(nats):
        coordsFlat_in[3 * i] = coords[i, 0]
        coordsFlat_in[3 * i + 1] = coords[i, 1]
        coordsFlat_in[3 * i + 2] = coords[i, 2]

    latticeVectorsFlat = np.zeros((9))
    latticeVectorsFlat[0] = latticeVectors[0, 0]
    latticeVectorsFlat[1] = latticeVectors[0, 1]
    latticeVectorsFlat[2] = latticeVectors[0, 2]

    latticeVectorsFlat[3] = latticeVectors[1, 0]
    latticeVectorsFlat[4] = latticeVectors[1, 1]
    latticeVectorsFlat[5] = latticeVectors[1, 2]

    latticeVectorsFlat[6] = latticeVectors[2, 0]
    latticeVectorsFlat[7] = latticeVectors[2, 1]
    latticeVectorsFlat[8] = latticeVectors[2, 2]

    fieldFlat_in[:] = field[:]

    # Specify arguments as a pointers to pass to Fortran
    latte_compute_hs_f.argtypes = [
        ct.c_int,
        ct.c_int,
        ct.c_int,
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.c_int,
    ]

    # Inputs
    coords_ptr = coordsFlat_in.ctypes.data_as(ct.POINTER(ct.c_double))
    atomTypes_ptr = atomTypes32.ctypes.data_as(ct.POINTER(ct.c_int))
    atomicNumbers_ptr = atomicNumbers.ctypes.data_as(ct.POINTER(ct.c_int))
    latticeVectors_ptr = latticeVectorsFlat.ctypes.data_as(ct.POINTER(ct.c_double))
    field_ptr = fieldFlat_in.ctypes.data_as(ct.POINTER(ct.c_double))

    # Outputs
    ham_ptr = hamFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    over_ptr = overFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))

    # Call to the fortran funtion
    err = latte_compute_hs_f(
        ct.c_int(norbs),
        ct.c_int(nats),
        ct.c_int(nTypes),
        coords_ptr,
        latticeVectors_ptr,
        atomTypes_ptr,
        atomicNumbers_ptr,
        field_ptr,
        ham_ptr,
        over_ptr,
        ct.c_int(verb),
    )

    # Back to a 2D array for the forces
    ham_out = np.zeros((norbs, norbs))
    over_out = np.zeros((norbs, norbs))

    for i in range(norbs):
        ham_out[:, i] = hamFlat_out[norbs * i : norbs + norbs * i]
        over_out[:, i] = overFlat_out[norbs * i : norbs + norbs * i]

    return err, ham_out, over_out


def latte_compute(latticeVectors, symbols, atomTypes, coords, field, verb):

    # Import the shared library
    latteLibFileName = os.environ["LATTE_PATH"] + "/liblatte.so"

    latteLib = ct.CDLL(latteLibFileName)
    latte_compute_f = latteLib.latte_compute
    # Periodic table: We use this to pass the chemical atom types as integer instead of characters.
    pt = PeriodicTable()

    nats = len(coords[:, 1])

    # Getting atomic numbers
    nTypes = len(symbols)
    atomicNumbers = np.zeros((nTypes), dtype=np.int32)
    atomTypes32 = np.zeros((nats), dtype=np.int32)
    atomTypes32[:] = atomTypes
    for i in range(len(symbols)):
        atomicNumbers[i] = pt.get_atomic_number(symbols[i])

    # Vectorizing 2D arrays for C-Fortran interoperability
    coordsFlat_in = np.zeros(3 * nats)  # Vectorized coordinates
    forcesFlat_out = np.zeros(3 * nats)  # Vectorized forces
    chargesFlat_out = np.zeros(nats)  # We call this one Flat just for consistency
    dipoleFlat_out = np.zeros(3)  # Same here
    fieldFlat_in = np.zeros(3)  # Same here
    energyFlat_out = np.zeros(1)

    for i in range(nats):
        coordsFlat_in[3 * i] = coords[i, 0]
        coordsFlat_in[3 * i + 1] = coords[i, 1]
        coordsFlat_in[3 * i + 2] = coords[i, 2]

    latticeVectorsFlat = np.zeros((9))
    latticeVectorsFlat[0] = latticeVectors[0, 0]
    latticeVectorsFlat[1] = latticeVectors[0, 1]
    latticeVectorsFlat[2] = latticeVectors[0, 2]

    latticeVectorsFlat[3] = latticeVectors[1, 0]
    latticeVectorsFlat[4] = latticeVectors[1, 1]
    latticeVectorsFlat[5] = latticeVectors[1, 2]

    latticeVectorsFlat[6] = latticeVectors[2, 0]
    latticeVectorsFlat[7] = latticeVectors[2, 1]
    latticeVectorsFlat[8] = latticeVectors[2, 2]

    fieldFlat_in[:] = field[:]

    # Specify arguments as a pointers to pass to Fortran
    latte_compute_f.argtypes = [
        ct.c_int,
        ct.c_int,
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_double),
        ct.c_int,
    ]

    # Inputs
    coords_ptr = coordsFlat_in.ctypes.data_as(ct.POINTER(ct.c_double))
    atomTypes_ptr = atomTypes32.ctypes.data_as(ct.POINTER(ct.c_int))
    atomicNumbers_ptr = atomicNumbers.ctypes.data_as(ct.POINTER(ct.c_int))
    latticeVectors_ptr = latticeVectorsFlat.ctypes.data_as(ct.POINTER(ct.c_double))
    field_ptr = fieldFlat_in.ctypes.data_as(ct.POINTER(ct.c_double))

    # Outputs
    charges_ptr = chargesFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    forces_ptr = forcesFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    dipole_ptr = dipoleFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    energy_ptr = energyFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    energy_out = 0.0

    # Call to the fortran funtion
    err = latte_compute_f(
        ct.c_int(nats),
        ct.c_int(nTypes),
        coords_ptr,
        latticeVectors_ptr,
        atomTypes_ptr,
        atomicNumbers_ptr,
        field_ptr,
        charges_ptr,
        forces_ptr,
        dipole_ptr,
        energy_ptr,
        ct.c_int(verb),
    )
    charges_out = np.zeros(nats)
    charges_out[:] = chargesFlat_out[:]
    dipole_out = np.zeros(3)
    dipole_out[:] = dipoleFlat_out[:]
    energy_out = energyFlat_out[0]

    # Back to a 2D array for the forces
    forces_out = np.zeros((nats, 3))
    for i in range(nats):
        forces_out[i, 0] = forcesFlat_out[i * 3 + 0]
        forces_out[i, 1] = forcesFlat_out[i * 3 + 1]
        forces_out[i, 2] = forcesFlat_out[i * 3 + 2]

    return err, charges_out, forces_out, dipole_out, energy_out
