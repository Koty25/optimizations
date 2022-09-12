#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par
i=36

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 36 $CAMINHO/Res_36_1

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 36 1 $CAMINHO/Res_36_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 2 18 $CAMINHO/Res_36_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 18 2 $CAMINHO/Res_36_1


mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 3 12 $CAMINHO/Res_36_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 12 3 $CAMINHO/Res_36_1


mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 4 9 $CAMINHO/Res_36_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 9 4 $CAMINHO/Res_36_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 6 6 $CAMINHO/Res_36_1


