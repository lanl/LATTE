Runnng latte as an MDI engine 
=============================

All the MDI imlementations are done in the `/src/latte_mdi.F90` file.

# Compile LATTE

Source the variables in `vars`

  $ source vars 

Use the `makefile.CHOICHES.mdi` to compile. 

  $ cp makefiles/makefile.CHOICHES.mdi makefile.CHOICES ; make 

# Compile dummy driver

  $ cd ./examples/MDI ; sh compile.sh 

# Run 

Run the driver as follws: 

  $ ./driver



