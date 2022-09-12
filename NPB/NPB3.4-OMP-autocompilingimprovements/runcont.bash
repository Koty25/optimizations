#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=tsubasa
#SBATCH --time=2:00:00
#SBATCH --job-name=cmp270-job2a
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

export PATH=/opt/nec/ve/nfort/3.0.4/bin:$PATH
export PATH=/opt/nec/ve/ncc/3.0.4/bin:$PATH
export VE_PROGINF=DETAIL
export VE_PERF_MODE=VECTOR-OP  #VECTOR-OP or VECTOR-MEM

source /home/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

bash ./compile.sh


step=1
timestamp=`date +"%Y%m%d_%I%M"`
expdir=experiments/exp_1
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

for i in {1..5}; do
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
      cat $outf | grep "MFLOPS                                  :" | awk '{printf "%s,", $3}' >> $expdir/times.${cf}.csv
      cat $outf | grep "V. Op. Ratio (%)                        :" | awk '{printf "%s,", $3}' >> $expdir/times.${cf}.csv
      cat $outf | grep "VLD LLC Hit Element Ratio (%)           :" | awk '{printf "%s,", $3}' >> $expdir/times.${cf}.csv
      cat $outf | grep "Time in seconds =" | awk '{print $5}' >> $expdir/times.${cf}.csv
   echo $ap,$aptime,$cf,$i >> $expdir/timest.${cf}.csv #application, runtime of application, architecture, step.

   nlines=`sed -n '/SECTION/,$p' $outf | wc -l`
   nlines="$(($nlines-1))"

   for line in  $(seq 1 $nlines); do #[[line=1; line<nlines; line++]]; do # lines for each ap (bt cg ep ft is mg sp ua)

      echo -n $ap, >> $expdir/test.${cf}.csv
      cat $outf | grep -A20 SECTION | tail -n+2 | awk -v line="$line" 'NR==line {printf "%s,", $1}' >> $expdir/test.${cf}.csv
      cat $outf | grep -A20 SECTION | tail -n+2 | awk -v line="$line" 'NR==line {printf "%s,", $3}' >> $expdir/test.${cf}.csv
      cat $outf | grep -A20 SECTION | tail -n+2 | awk -v line="$line" 'NR==line {printf $5}' >> $expdir/test.${cf}.csv
      echo >> $expdir/test.${cf}.csv

   done

      rm $tmpf
   else
      echo "... run $cf/$pgm nt=$nt - not present" | tee -a $logfile
      ((cntp++))
   fi
done
done

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

