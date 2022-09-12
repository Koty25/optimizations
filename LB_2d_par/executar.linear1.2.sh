#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par

##for i in `seq -f %0g 5 8`
for i in 6
do
	mpirun -machinefile $CAMINHO/../arquivo_nos_labtec9 -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res $i 1 $CAMINHO/Res_linha_seq_1
done

#for i in `seq -f %0g 5 8`
for i in 6
do
	mpirun -machinefile $CAMINHO/../arquivo_nos_labtec9 -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 $i $CAMINHO/Res_coluna_seq_1
done
