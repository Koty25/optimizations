#/bin/bash

CAMINHO=/home/gppd/cschepke/LB_3d_par

for i in `seq -f %0g 1 5`
do
	sh $CAMINHO/executar.sh
done
