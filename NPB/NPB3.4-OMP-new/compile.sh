#!/bin/csh


logfile=npb-make.log
touch $logfile
outf=npb-make.out
touch $outf

echo "Date: `date`" >> $logfile
echo "Host: `hostname`" >> $logfile
echo "" >> $logfile

cnt=0
cntf=0

aps=()
spv=()
c="C"

for cf in ncc intel; do

bindir=bin/bin_$cf
if [ ! -d $bindir ]; then
 mkdir $bindir
fi
cp -f config/NAS.samples/make.def_$cf config/make.def
make clean >> $outf

for ap in bt cg ep ft is lu mg sp ua; do
   make $ap CLASS=$c &>> $outf
   pgm=${ap}.${c}.x
   pgmx=bin/$pgm
   ((cnt++))
   if [ -e $pgmx ]; then
      mv $pgmx $bindir
      echo ">>> make $cf/$pgm - successful" | tee -a $logfile
   else
      echo "*** make $cf/$pgm - FAILED" | tee -a $logfile
      ((cntf++))
   fi
done
done

echo "" >> $logfile
echo "Date: `date`" >> $logfile
echo "Total number of cases: $cnt" | tee -a $logfile
echo "Total number of FAILED cases: $cntf" | tee -a $logfile
echo "" >> $logfile
