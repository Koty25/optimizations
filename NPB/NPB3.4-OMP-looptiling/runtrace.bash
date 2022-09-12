#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=tsubasa
#SBATCH --time=2:00:00
#SBATCH --job-name=cmp270-jobtrace
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# QUANDO RODAR ESSE SCRIPT, POR FAVOR, ADICIONAR A FLAG -ftrace NO ARQUIVO ./config/NAS.samples/make.def_ncc


export PATH=/opt/nec/ve/nfort/3.0.4/bin:$PATH
export PATH=/opt/nec/ve/ncc/3.0.4/bin:$PATH
export VE_PROGINF=DETAIL
export VE_PERF_MODE=VECTOR-OP  #VECTOR-OP or VECTOR-MEM
source /home/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

bash ./compile.sh


step=1
timestamp=`date +"%Y%m%d_%I%M"`
expdir=experiments/exp_2_og_classB-trace1-noftrace-fopenmp-loopunroll-inline
if [ ! -d $expdir ]; then
 mkdir $expdir
fi

logfile=$expdir/npb-run.$step.log
touch $logfile
tmpf=$expdir/npb.tmp.$$

echo "Date: `date`" >> $logfile
echo "Host: `hostname`" >> $logfile
echo "" >> $logfile

cnt=0
cntf=0
cntp=0

aps=(bt sp lu lu ua ua)
spv=(blk blk hp doac au rd)
c="B"
export NPB_TIMER_FLAG=1

for nt in 8; do
for cf in ncc intel; do

bindir=bin/bin_$cf
outdir=$expdir/out_$cf
if [ ! -d $outdir ]; then
 mkdir -p $outdir
fi

for i in {1..1}; do
for ap in bt cg ep ft is mg sp ua; do

   pgm=${ap}.${c}.x
   pgmx=$bindir/$pgm
   ((cnt++))
   if [[ -e $pgmx ]]; then
      outf=$outdir/${ap}.${c}.out.$nt
      touch $outf
      export OMP_NUM_THREADS=$nt
	res1=$(date +"%s.%3N")
      $pgmx !&>> $tmpf
	res2=$(date +"%s.%3N")
	aptime=$(echo "$res2-$res1" | bc)
      
      grep -i successful $tmpf &>> /dev/null
      if [[ $? = 0 ]]; then
         echo ">>> run $cf/$pgm nt=$nt - successful" | tee -a $logfile
      else
         echo "*** run $cf/$pgm nt=$nt - FAILED" | tee -a $logfile
         ((cntf++))
      fi
      cat $tmpf &> $outf
      cat $outf | grep "NAS Parallel Benchmarks (NPB3.4-OMP)" | awk '{printf "%s,", $6}' >> $expdir/times.${cf}.csv
      cat $outf | grep "Time in seconds =" | awk '{print $5}' >> $expdir/times.${cf}.csv
	echo $ap,$aptime,$cf,$i >> $expdir/timest.${cf}.csv #application, runtime of application, architecture, step.

      rm $tmpf
   else
      echo "... run $cf/$pgm nt=$nt - not present" | tee -a $logfile
      ((cntp++))
   fi
   mv ./ftrace.out $outdir/ftrace.${ap}.out
done
done

#mv $outdir $expdir
#Rscript trace.r "${cf}" &>> $logfile
#Rscript trace_ana.r.R "${cf}" all &>> $logfile
#pode-se usar o trace_ana.r.R para criar plots dos benchamarks especificos. Somente mude o argumento "all" para um dos benchmarks!
done
done


echo "" >> $logfile
echo "Date: `date`" >> $logfile
echo "Total number of cases: $cnt" | tee -a $logfile
echo "Total number of FAILED cases: $cntf" | tee -a $logfile
echo "Total number of not present cases: $cntp" | tee -a $logfile
echo "" >> $logfile
#mv $logfile $expdir
mv npb-make.log $expdir
mv npb-make.out $expdir
#mv times*.csv $expdir
cd $expdir
#Rscript ~/NPB/NPB3.4-OMP/combinedtrace.r &>> $logfile


