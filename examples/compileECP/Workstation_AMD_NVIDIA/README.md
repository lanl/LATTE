INSTALLATION INSTRUCTIONS
========================
Here we assume that we are in `$HOME` as the base directory (`cd `) 

Git clone spack to machine:

	git clone -b v1.1.1 https://github.com/spack/spack.git spack_latte_lammps 

Run:

	. ~/spack_latte_lammps/share/spack/setup-env.sh  

Check current compilers that spack knows about:

	spack compiler find 

View all these compilers: 

	vi ~/.spack/packages.yaml

Pick a system compiler you typically use and remove all the rest of the 
blocks of codes for the other compilers.

Now install the compiler you want (for example, pick your favorite version):

	spack install -j 64 gcc@12.3.0 

Now create environment:

	spack env create latte_lammps

Enter environment:

	spack env activate latte_lammps -p

Build bml:

        spack add bml@master+cusolver+magma%gcc@12.3.0 ^magma+cuda cuda_arch=89

        spack concretize

        spack install -v -j 64    (-v = verbose, -j = build in parallel w 64 cpus)

Clone the repo: 

	git clone -b ECP  git@github.com:lanl/LATTE.git LATTE_ECP

Get into the example ecp compilation:

	cd $HOME/LATTE_ECP/examples/compileECP/Workstation_AMD_NVIDIA

Clone all module (this will clone LATTE again):

	./clone_all_codes.sh

Create a source file with all the relevant path:

	./get_spack_paths ; source paths.sh 

Build and install progress:

	./build_progress.sh 

Build and install LATTE:

	./build_latte.sh

Build and install LAMMPS:

	./build_lammps


## Running an example via LAMMPS 

Get into the example folders:

	cd $HOME/LATTE_ECP/examples/compileECP/Workstation_AMD_NVIDIA/lammps/examples/latte

	../../src/lmp_serial < in.latte.sucrose

From another terminal, one can type `nvidia-smi` to see the GPU load and double check that the code is using the available GPUS.

