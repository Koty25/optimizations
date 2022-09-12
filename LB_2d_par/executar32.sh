#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par
i=32

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 32 $CAMINHO/Res_32_1

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 32 1 $CAMINHO/Res_32_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 2 16 $CAMINHO/Res_32_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 16 2 $CAMINHO/Res_32_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 4 8 $CAMINHO/Res_32_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 8 4 $CAMINHO/Res_32_1

