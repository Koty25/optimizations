#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par
i=40

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 40 $CAMINHO/Res_40_1

#mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 40 1 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 2 20 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 20 2 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 4 10 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 10 4 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 5 8 $CAMINHO/Res_40_1

mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 8 5 $CAMINHO/Res_40_1

