#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/gppd/cschepke/LB_2d_par

##for i in `seq -f %0g 17 40`
for i in 8 #40 38 36 34 32 30 28 26 24 22 20 18
do
	mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res $i 1 $CAMINHO/Res_linha_seq_1
done

##for i in `seq -f %0g 17 40`
#for i in 40 #38 36 34 32 30 28 26 24 22 20 18 #20 22 24 26 28 30 32 34 36 38 40
#do
#	mpirun -machinefile $CAMINHO/../arquivo_nos_labtec -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 $i $CAMINHO/Res_coluna_seq_1
#done
