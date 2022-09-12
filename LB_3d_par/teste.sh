CAMINHO=/home/gppd/cschepke/LB_3d_par

mpirun -machinefile /home/gppd/cschepke/arquivo_nos_labtec -np 2 $CAMINHO/lb $CAMINHO/anb.par $CAMINHO/anb.obs11 $CAMINHO/anb.res 1 1 2 $CAMINHO/Res
