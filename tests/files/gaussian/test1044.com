%mem=1600mb
%chk=test1044
#p b3lyp/sto-3g extrabasis 5d 7f td=(singlets,maxdav=1000,nstates=6,conver=6,eqsolv) scrf int=(acc2e=12,grid=ultrafine) test

Gaussian Test Job 1044 (Part 1):
Singlet part of 50-50 test

1 1
 C                 -1.71345031    0.47368420    0.00000000
 C                 -0.31829031    0.47368420    0.00000000
 C                  0.37924769    1.68143520    0.00000000
 C                 -0.31840631    2.88994420   -0.00119900
 C                 -1.71323131    2.88986620   -0.00167800
 C                 -2.41083231    1.68166020   -0.00068200
 H                  0.23121769   -0.47882880    0.00131500
 H                  1.47892769    1.68151520    0.00063400
 H                 -2.26335331    3.84214720   -0.00263100
 H                 -3.51043631    1.68184320   -0.00086200
 C                 -2.11359106   -1.68562799    1.19745661
 H                 -2.24052229   -1.40337812   -0.94478398
 H                 -3.57689783   -0.63148543    0.00052223
 H                 -2.67190382   -2.65313022    1.19832897
 H                 -1.01998863   -1.91390717    1.19725366
 H                 -2.35634965   -1.14191707    2.14271910
 N                 -2.48338567   -0.86003209    0.00063022
 O                  0.39706176    4.12809017   -0.00127572
 H                  1.34083210    3.95231320   -0.00126483

c n o 0
d 1 1.0
0.8 1.0
****

--Link1--
%mem=1600mb
%chk=test1044
#p b3lyp/chkbas 5d 7f td=(triplets,maxdav=1000,nstates=6,conver=6,eqsolv) scrf int=(acc2e=12,grid=ultrafine) geom=check guess=read test

Gaussian Test Job 1044 (Part 3):
Triplet part of 50-50 test

1 1

--Link1--
%mem=1600mb
%chk=test1044
%nosave
#p b3lyp/chkbas 5d 7f td=(50-50,maxdav=1000,nstates=6,conver=6,eqsolv) scrf int=(acc2e=12,grid=ultrafine) geom=check guess=read test

Gaussian Test Job 1044 (Part 3):
50-50 part of 50-50 test

1 1

