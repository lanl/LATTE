#!/bin/bash -e
rm -rf bml; git clone git@github.com:lanl/bml.git
rm -rf qmd-progress; git clone git@github.com:lanl/qmd-progress.git
rm -rf LATTE; git clone -b ECP git@github.com:lanl/LATTE.git
rm -rf lammps; git clone https://github.com/lammps/lammps.git

