# To be used with the latte-lib input file.  

units		metal
atom_style	full
atom_modify    sort 0 0.0  # This is to avoid sorting the coordinates

read_data data.wat_24.lmp

velocity	all create 0.0 87287 loop geom
#velocity all zero linear

pair_style	lj/cut/coul/cut 8.0
pair_style zero 8.0
pair_coeff	* *  

neighbor	1.0 bin
neigh_modify every 1 delay 0 check yes 

timestep 0.00025

compute coulomb all pe/atom 

# Compute temperature, pressure, potential energy, kinetic energy and total energy
#compute   cT  all temp
#compute   cP  all pressure thermo_temp

fix		1 all nve

fix   2 all latte NULL

# To print some useful data
#variable latteE equal "(ke + f_2)/23.06035"
#variable kinE equal "ke/23.06035"
#variable potE equal "f_2/23.06035"
#variable myT equal "temp"

variable latteE equal "(ke + f_2)"
variable kinE equal "ke"
variable potE equal "f_2"
variable myT equal "temp"


fix 3 all print 1 "Total Energy = ${latteE}"
fix 4 all print 1 "Kin = ${kinE}"
fix 5 all print 1 "Pot = ${potE}"
fix 6 all print 1 "Temp = ${myT}"

#dump dmpvtk all custom 1 dump*.myforce.vtk id type vx fx

run_style verlet

run		100

