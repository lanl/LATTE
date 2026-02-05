INSTALATION INSTRUCTIONS
========================

Load an appropriate environment with the mpif90 compiler, oenblas, cublas, magma, and metis library.

Clone the repo: 

	git clone -b ECP  git@github.com:lanl/LATTE.git LATTE_ECP

Get into the example ecp compilation

	cd $HOME/LATTE_ECP/examples/compileECP/Workstation_AMD_NVIDIA

Clone all module (this will clone LATTE again)

	./clone_all_codes.sh

Build and install bml

	./build_bml.sh

Build and install progress

	./build_progress.sh 

Build and install LATTE

	./build_latte.sh

Build and install LAMMPS

	./build_lammps
