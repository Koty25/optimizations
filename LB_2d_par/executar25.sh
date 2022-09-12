#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par
i=25

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 25 $CAMINHO/Res_25_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 25 1 $CAMINHO/Res_25_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 5 5 $CAMINHO/Res_25_1

