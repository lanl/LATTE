Example using latte with lammps 
===============================

Copy any of the latte.*.in files into latte.in 

The wat.lmp file is generated using a PROGRESS tool as: 
  
    ~/qmd-progress/build/changecoords wat.pdb wat.lmp    

Run lammps doing: 

    ~/lammps/src/lmp_serial < in.file 

latte.fullscf.in: To perform QMD with full SCF at each time step.

latte.0SCF.in: To perform XLBOMD with only one diagonalization per time step.

latte.Pulay.in: Includes the use of the Pulay mixing scheme (Only when compiled with PROGRESS/BML libraries)
