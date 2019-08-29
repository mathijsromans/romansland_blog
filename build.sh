#!/bin/bash

set -E

mkdir -p part5
cd part5
octave ../part5.m
convert -delay 20 -loop 0 chart_*.png myimage.gif
cd ..
