#!/bin/sh

export OMP_NUM_THREADS=24

../../LATTE-skim/LATTE_DOUBLE < latte.in | tee out
