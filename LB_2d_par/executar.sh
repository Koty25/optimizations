#/bin/bash
#for i in 20000 30000 40000 50000 60000 70000 80000 90000
for i in `seq -f %0g 1 10`
do
	for j in `seq -f %0g 1 16`
	do
		#/usr/local/mpich-gm/bin/
		mpirun -machinefile ~/arquivo_nos -np $j ~/LB/lb ~/LB/anb.par ~/LB/anb.obs3 >> ~/LB/C1eth200x50_$j.txt
	done
done

