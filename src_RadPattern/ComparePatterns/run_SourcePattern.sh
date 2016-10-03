#!/bin/bash


### main ###
feigR=245_41.25_dense.R; feigL=245_41.25_dense.L
#feigR=245_41.25_broad.R; feigL=245_41.25_broad.L

#stk1=247.8; dip1=33.3; rak1=-59.6; dep1=5.2
#stk2=247.8; dip2=33.3; rak2=120.4; dep1=5.2

stk1=180; dip1=45; rak1=90; dep1=10; M01=1
stk2=60.431; dip2=72.865; rak2=-118.431; dep2=0.003; M02=0.55



../PredSourcePattern R $feigR perlst ${dep1} ${M01} ${stk1} ${dip1} ${rak1} solution1_R
../PredSourcePattern L $feigL perlst ${dep1} ${M01} ${stk1} ${dip1} ${rak1} solution1_L
../PredSourcePattern R $feigR perlst ${dep2} ${M02} ${stk2} ${dip2} ${rak2} solution2_R
../PredSourcePattern L $feigL perlst ${dep2} ${M02} ${stk2} ${dip2} ${rak2} solution2_L

#../PredSourcePattern R $feigR perlst 10 1 0 -1 1 0 0 0 solution1_R
#../PredSourcePattern L $feigL perlst 10 1 0 -1 1 0 0 0 solution1_L
#../PredSourcePattern R $feigR perlst 0.003 0.55 0.76514 -0.26996 -0.49518 0.02085 -0.56286 0.48062 solution2_R
#../PredSourcePattern L $feigL perlst 0.003 0.55 0.76514 -0.26996 -0.49518 0.02085 -0.56286 0.48062 solution2_L
