# To be used with the latte-lib input file.  

units		metal
atom_style	full
atom_modify    sort 0 0.0  # This is to avoid sorting the coordinates

read_data coords.lmp

velocity	all create 300.0 87287 loop geom

pair_style zero 1.0
pair_coeff	* *  

neighbor	1.0 bin
neigh_modify every 1 delay 0 check yes 

timestep 0.00025

compute coulomb all pe/atom 

fix		1 all nvt temp 300 300  10

fix   2 all latte NULL
fix_modify      2 energy yes

variable latteE equal "(ke + f_2)"
variable kinE equal "ke"
variable potE equal "f_2"
variable myT equal "temp"
variable myP equal "press"

fix 3 all print 1 "Total Energy = ${latteE}"
fix 4 all print 1 "Kin = ${kinE}"
fix 5 all print 1 "Pot = ${potE}"
fix 6 all print 1 "Temp = ${myT}"
fix 7 all print 1 "Press = ${myP}"

run_style verlet

thermo 1
thermo_modify flush yes
run		10000

