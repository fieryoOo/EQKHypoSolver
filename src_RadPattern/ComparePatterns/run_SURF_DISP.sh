#!/bin/bash

#fmod=/projects/yeti4009/Model/US/41.25_245.00
fmod=./41.25_245.00
/projects/yeti4009/code/Programs/SURF_DISP/SURF_DISP $fmod 245_41.25_broad R 0 2 2 60 1 -a -f
/projects/yeti4009/code/Programs/SURF_DISP/SURF_DISP $fmod 245_41.25_broad L 0 2 2 60 1 -a -f
