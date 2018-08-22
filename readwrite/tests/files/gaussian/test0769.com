%chk=test0769
#p freq=noraman oniom(external="mygau o g09 b3lyp 6-311G\(d\,p\)":uff=qeq) nosymm geom=connectivity test

Gaussian Test job 769 (Part 1):
Test of External procedure using Gaussian as the external program.

0 1 0 1 0 1
 C-C_3            0    0.000000    0.000000    0.000000 H
 C-C_3            0    0.000000    0.000000    1.566200 H
 C-C_R            0    1.473182    0.000000   -0.492419 L H-H_      1
 C-C_R            0   -0.736591   -1.275813   -0.492419 L H-H_      1
 C-C_R            0   -0.736591    1.275813   -0.492419 L H-H_      1
 C-C_R            0    2.412644   -0.869263    0.097427 L
 C-C_R            0    3.741978   -0.900299   -0.333237 L
 C-C_R            0    4.167062   -0.063464   -1.372126 L
 C-C_R            0    3.245068    0.795956   -1.975494 L
 C-C_R            0    1.913283    0.825913   -1.540967 L
 H-H_             0    2.104823   -1.538161    0.893372 L
 H-H_             0    4.444197   -1.578833    0.141800 L
 H-H_             0    5.199639   -0.084597   -1.706264 L
 H-H_             0    3.556903    1.447801   -2.785937 L
 H-H_             0    1.215276    1.500721   -2.021098 L
 C-C_R            0   -0.453518    2.524042    0.097427 L
 C-C_R            0   -1.091307    3.690798   -0.333237 L
 C-C_R            0   -2.028569    3.640513   -1.372126 L
 C-C_R            0   -2.311852    2.412333   -1.975494 L
 C-C_R            0   -1.671903    1.243996   -1.540967 L
 H-H_             0    0.279675    2.591911    0.893372 L
 H-H_             0   -0.854789    4.638204    0.141800 L
 H-H_             0   -2.526556    4.545318   -1.706264 L
 H-H_             0   -3.032284    2.356467   -2.785937 L
 H-H_             0   -1.907300    0.302099   -2.021098 L
 C-C_R            0   -1.959126   -1.654779    0.097427 L
 C-C_R            0   -2.650671   -2.790499   -0.333237 L
 C-C_R            0   -2.138492   -3.577049   -1.372126 L
 C-C_R            0   -0.933216   -3.208289   -1.975494 L
 C-C_R            0   -0.241380   -2.069908   -1.540967 L
 H-H_             0   -2.384498   -1.053750    0.893372 L
 H-H_             0   -3.589408   -3.059371    0.141800 L
 H-H_             0   -2.673082   -4.460721   -1.706264 L
 H-H_             0   -0.524619   -3.804269   -2.785937 L
 H-H_             0    0.692024   -1.802820   -2.021098 L
 H-H_             0   -0.996673    0.221850    1.956448 H
 H-H_             0    0.690464    0.752219    1.956448 H
 H-H_             0    0.306209   -0.974069    1.956448 H

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2 36 1.0 37 1.0 38 1.0
 3 6 1.5 10 1.5
 4 26 1.5 30 1.5
 5 16 1.5 20 1.5
 6 7 1.5 11 1.0
 7 8 1.5 12 1.0
 8 9 1.5 13 1.0
 9 10 1.5 14 1.0
 10 15 1.0
 11
 12
 13
 14
 15
 16 17 1.5 21 1.0
 17 18 1.5 22 1.0
 18 19 1.5 23 1.0
 19 20 1.5 24 1.0
 20 25 1.0
 21
 22
 23
 24
 25
 26 27 1.5 31 1.0
 27 28 1.5 32 1.0
 28 29 1.5 33 1.0
 29 30 1.5 34 1.0
 30 35 1.0
 31
 32
 33
 34
 35
 36
 37
 38

--Link1--
%chk=test0769
#p freq=noraman oniom(external=("mygau c g09 b3lyp 6-311G\(d\,p\)",iofchk)/6-311g**:uff=qeq) nosymm geom=check test

Gaussian Test job 769 (Part 2):
Test of External procedure using Gaussian as the external program,
results from an fchk file.

0 1 0 1 0 1

--Link1--
%chk=test0769
#p freq=noraman oniom(external=("mygau c g09 b3lyp 6-311G\(d\,p\)",iofchk,1elect)/6-311g**:uff=qeq) nosymm geom=check test

Gaussian Test job 769 (Part 3):
Test of External procedure using Gaussian as the external program,
results from an fchk file, writing 1e integral file

0 1 0 1 0 1

--Link1--
%chk=test0769
%nosave
#p freq=noraman oniom(external=("mygau c g09 b3lyp 6-311G\(d\,p\)",iofchk,2elect)/6-311g**:uff=qeq) nosymm geom=check test

Gaussian Test job 769 (Part 4):
Test of External procedure using Gaussian as the external program,
results from an fchk file, writing 1e/2e integral file

0 1 0 1 0 1

