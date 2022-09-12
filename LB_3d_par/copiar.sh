for i in `seq -f %0g 1 1`
do
	a=1
	b=`echo "128*($i - 1) + 3 - 1"|bc -l`
	c=2
	d=`echo "128*128*128 - $b" | bc -l`
	cat anb.res.bak | sed '$a,$b' | sed '$c,$d' >> 3d44.txt
done


#cat anb.res | sed -r '/^44/!d'
