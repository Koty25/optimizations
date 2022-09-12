#/bin/bash

CAMINHO=/home/gppd/cschepke/LB_2d_par

for i in `seq -f %0g 1 3`
do
	sh $CAMINHO/executar36.sh
done
