#/bin/bash

#mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence

#/usr/local/mpich-gm/bin/

CAMINHO=/home/users/fpmjunior/ccpe/LB_3d_par


for i in `seq -f %0g 37 40`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res $i 1 1 $CAMINHO/Res40_labtec_linha_seq_10
done

for i in `seq -f %0g 37 40`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 $i 1 $CAMINHO/Res40_labtec_coluna_seq_10
done

for i in `seq -f %0g 37 40`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 1 $i $CAMINHO/Res40_labtec_profundidade_seq_10
done


for i in `seq -f %0g 25 36`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res $i 1 1 $CAMINHO/Res40_labtec_linha_seq_10
done

for i in `seq -f %0g 25 36`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 $i 1 $CAMINHO/Res40_labtec_coluna_seq_10
done

for i in `seq -f %0g 25 36`
do
	mpirun -machinefile $CAMINHO/../arquivo_nos -np $i $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs41 $CAMINHO/anb.res 1 1 $i $CAMINHO/Res40_labtec_profundidade_seq_10
done



