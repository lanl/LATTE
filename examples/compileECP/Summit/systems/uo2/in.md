#LAMMPS (23 Jun 2022)
# AIMD test of two UO2 molecules with LATTE in MDI stand-alone mode

units           metal
atom_style      full
atom_modify     sort 0 0.0

read_data       uo2.lmp

variable  x index 2
variable  y index 1
variable  z index 1

variable        nrep equal v_x*v_y*v_z
if              "${nrep} > 1" then "replicate $x $y $z"

velocity        all create 300.0 87287 loop geom

neighbor        0.0 bin
neigh_modify    every 1 delay 0 check yes

pair_style      zero 1.0
pair_coeff  * *

timestep        0.00025

fix             1 all nve
fix   2 all latte NULL
fix_modify 2 energy yes

thermo 1
thermo_style    custom step temp pe etotal press
thermo_modify flush yes

dump            1 all custom 1 dump.mdi.aimd                 id type x y z vx vy vz fx fy fz

run             14

