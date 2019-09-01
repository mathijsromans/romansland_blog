#!/bin/bash

set -E

mkdir -p build/part5
octave part5.m
cd build/part5
convert -delay 20 -loop 0 chart_*.png all_potentials.gif
cd ../..
