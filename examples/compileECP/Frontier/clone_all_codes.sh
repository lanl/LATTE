#!/bin/bash -e
rm -rf bml; git clone git@github.com:lanl/bml.git
rm -rf qmd-progress; git clone git@github.com:lanl/qmd-progress.git
rm -rf LATTE; git clone -b ECP git@github.com:lanl/LATTE.git
wget https://github.com/lammps/lammps/archive/refs/tags/patch_23Jun2022_update2.tar.gz
rm -rf lammps; tar -xvzf patch_23Jun2022_update2.tar.gz; mv patch_23Jun2022_update2 lammps


