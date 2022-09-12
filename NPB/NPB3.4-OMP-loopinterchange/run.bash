#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=tsubasa
#SBATCH --time=2:00:00
#SBATCH --job-name=cmp270-job
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

export PATH=/opt/nec/ve/nfort/3.0.4/bin:$PATH
export PATH=/opt/nec/ve/ncc/3.0.4/bin:$PATH
source /home/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

bash ./compile.sh

logfile=npb-run.log
touch $logfile
tmpf=npb.tmp.$$

timestamp=`date +"%Y%m%d_%I%M"`
expdir=experiments/exp_$timestamp
if [ ! -d $expdir ]; then
 mkdir $expdir
fi

echo "Date: `date`" >> $logfile
echo "Host: `hostname`" >> $logfile
echo "" >> $logfile

cnt=0
cntf=0
cntp=0

aps=()
spv=()
c="B"
export NPB_TIMER_FLAG=1

for nt in 4; do
for cf in ncc intel; do

bindir=bin/bin_$cf
outdir=bin/out_$cf
if [ ! -d $outdir ]; then
 mkdir -p $outdir
fi

for ap in bt cg ep ft mg sp ua; do
for i in {0..4}; do

   pgm=${ap}.${c}.x
   pgmx=$bindir/$pgm
   ((cnt++))
   if [[ -e $pgmx ]]; then
      outf=$outdir/${ap}.${c}.out.$nt
      touch $outf
      export OMP_NUM_THREADS=$nt
      $pgmx !&>> $tmpf
      
      grep -i successful $tmpf &>> /dev/null
      if [[ $? = 0 ]]; then
         echo ">>> run $cf/$pgm nt=$nt - successful" | tee -a $logfile
      else
         echo "*** run $cf/$pgm nt=$nt - FAILED" | tee -a $logfile
         ((cntf++))
      fi
      cat $tmpf &> $outf
      cat $outf | grep "NAS Parallel Benchmarks (NPB3.4-OMP)" | awk '{printf "%s,", $6}' >> times.${cf}.csv
      cat $outf | grep "Time in seconds =" | awk '{print $5}' >> times.${cf}.csv
      rm $tmpf
   else
      echo "... run $cf/$pgm nt=$nt - not present" | tee -a $logfile
      ((cntp++))
   fi
done
done
mv $outdir $expdir
Rscript trace.r "${cf}" &>> $logfile
done
done



echo "" >> $logfile
echo "Date: `date`" >> $logfile
echo "Total number of cases: $cnt" | tee -a $logfile
echo "Total number of FAILED cases: $cntf" | tee -a $logfile
echo "Total number of not present cases: $cntp" | tee -a $logfile
echo "" >> $logfile
mv $logfile $expdir
mv npb-make.log $expdir
mv npb-make.out $expdir
mv times*.csv $expdir
cd $expdir
Rscript ~/NPB/NPB3.4-OMP/combinedtrace.r &>> $logfile

