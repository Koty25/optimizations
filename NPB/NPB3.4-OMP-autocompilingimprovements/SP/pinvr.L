NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME : pinvr.f

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fassume-contiguous                      : disable
        -fcopyin-intent-out                      : enable
        -fcse-after-vectorization                : disable
        -ffast-formatted-io                      : enable
        -ffast-math                              : enable
        -fignore-asynchronous                    : disable
        -fignore-induction-variable-overflow     : enable
        -fignore-volatile                        : enable
        -fivdep                                  : disable
        -floop-collapse                          : enable
        -floop-count                             : 5000
        -floop-fusion                            : enable
        -floop-interchange                       : enable
        -floop-normalize                         : enable
        -floop-split                             : disable
        -floop-strip-mine                        : enable
        -floop-unroll                            : enable
        -floop-unroll-complete                   : 4
        -floop-unroll-complete-nest              : 3
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : enable
        -fmove-nested-loop-invariants-outer      : enable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fnamed-noalias-aggressive               : disable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freorder-logical-expression             : enable
        -freplace-loop-equation                  : enable
        -marray-io                               : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mstack-arrays                           : enable
        -mno-stack-arrays                        : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : enable
        -mvector-loop-count-test                 : disable
        -mvector-low-precise-divide-function     : disable
        -mvector-merge-conditional               : enable
        -mvector-packed                          : disable
        -mvector-power-to-explog                 : disable
        -mvector-power-to-sqrt                   : enable
        -mvector-reduction                       : enable
        -mvector-shortloop-reduction             : disable
        -mvector-sqrt-instruction                : disable
        -mvector-threshold                       : 5
        -mwork-vector-kind=none                  : disable

      Parallelization Options : 
        -fopenmp                                 : enable
        -mparallel                               : disable
        -mparallel-innerloop                     : disable
        -mparallel-omp-routine                   : enable
        -mparallel-outerloop-strip-mine          : disable
        -mparallel-sections                      : disable
        -mparallel-threshold                     : 2000
        -mschedule-dynamic                       : disable
        -mschedule-runtime                       : disable
        -mschedule-static                        : disable
        -mschedule-chunk-size                    : 0
        -pthread                                 : disable

      Inlining Options : 
        -finline-copy-arguments                  : disable
        -finline-functions                       : enable
        -finline-max-depth                       : 2
        -finline-max-function-size               : 50
        -finline-max-times                       : 6

      Code Generation Options : 
        -finstrument-functions                   : disable
        -fpic | -fPIC                            : disable
        -ftrace                                  : disable
        -p                                       : disable

      Debugging Options : 
        -fbounds-check                           : disable
        -fcheck=all                              : disable
        -fcheck=alias                            : disable
        -fcheck=bits                             : disable
        -fcheck=bounds                           : disable
        -fcheck=dangling                         : disable
        -fcheck=do                               : disable
        -fcheck=iovf                             : disable
        -fcheck=pointer                          : disable
        -fcheck=present                          : disable
        -fcheck=recursion                        : disable
        -g                                       : disable
        -minit-stack                             : not specified
        -mmemory-trace                           : disable
        -mmemory-trace-full                      : disable
        -traceback                               : disable

      Language Options : 
        -bss                                     : disable
        -fdefault-integer                        : 4
        -fdefault-double                         : 8
        -fdefault-real                           : 4
        -fextend-source                          : disable
        -ffree-form                              : disable
        -ffixed-form                             : enable
        -ff90-sign                               : disable
        -fmax-continuation-lines                 : 1023
        -frealloc-lhs                            : enable
        -frealloc-lhs-array                      : enable
        -frealloc-lhs-scalar                     : enable
        -save                                    : disable
        -std                                     : Fortran 2008

      Message Options : 
        -Wall                                    : disable
        -Werror                                  : disable
        -Wextension                              : enable
        -Wobsolescent                            : enable
        -Woverflow                               : disable
        -Woverflow-errors                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -pedantic-errors                         : disable
        -w                                       : disable

      List Output Options : 
        -report-file                             : not specified
        -report-append-mode                      : disable
        -report-all                              : enable
        -report-cg                               : enable
        -report-diagnostics                      : enable
        -report-format                           : enable
        -report-inline                           : enable
        -report-option                           : enable
        -report-vector                           : enable

      Directory Options : 
        -fintrinsic-modules-path                 : /opt/nec/ve/nfort/3.2.0/include
        -module                                  : not specified

      Miscellaneous Options : 
        -v                                       : disable


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: pinvr.f

PROCEDURE NAME: PINVR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    20: inl(1212): Source for routine not found.: TIMER_START
    21: par(1801): Parallel routine generated.: PINVR$1
    23: par(1803): Parallelized by "do".
    23: opt(2495): Loop nest collapsed by "collapse" clause.
    23: vec( 101): Vectorized loop.
    25: vec( 101): Vectorized loop.
    39: vec( 128): Fused multiply-add operation applied.
    44: inl(1212): Source for routine not found.: TIMER_STOP


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: pinvr.f

PROCEDURE NAME: PINVR
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

     5:                  subroutine pinvr
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   block-diagonal matrix-vector multiplication
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  integer i, j, k
    18:                  double precision r1, r2, r3, r4, r5, t1, t2
    19:           
    20:                  if (timeron) call timer_start(t_pinvr)
    21:           !$omp parallel do default(shared) private(i,j,k,r1,r2,r3,r4,r5,t1,t2)
    22:           !$omp&  collapse(2)
    23: Y------>         do   k = 1, nz2
    24: |+----->            do   j = 1, ny2
    25: ||V---->               do   i = 1, nx2
    26: |||       
    27: |||                       r1 = rhs(1,i,j,k)
    28: |||                       r2 = rhs(2,i,j,k)
    29: |||                       r3 = rhs(3,i,j,k)
    30: |||                       r4 = rhs(4,i,j,k)
    31: |||                       r5 = rhs(5,i,j,k)
    32: |||       
    33: |||                       t1 = bt * r1
    34: |||                       t2 = 0.5d0 * ( r4 + r5 )
    35: |||       
    36: |||                       rhs(1,i,j,k) =  bt * ( r4 - r5 )
    37: |||                       rhs(2,i,j,k) = -r3
    38: |||                       rhs(3,i,j,k) =  r2
    39: |||     F                 rhs(4,i,j,k) = -t1 + t2
    40: |||                       rhs(5,i,j,k) =  t1 + t2
    41: ||V----                end do
    42: |+-----             end do
    43: Y------          end do
    44:                  if (timeron) call timer_stop(t_pinvr)
    45:           
    46:                  return
    47:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: pinvr.f

PROCEDURE NAME: PINVR
INLINE LIST

  ROOT: PINVR (pinvr.f:5)
  -> NOINLINE: TIMER_START (pinvr.f:20)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (pinvr.f:44)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: pinvr.f

PROCEDURE NAME: PINVR
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (pinvr.f:23)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (pinvr.f:23)
    *** The number of VLOAD, VSTORE. :  0,  2. (pinvr.f:23)

    LOOP BEGIN: (pinvr.f:25)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (pinvr.f:25)
      *** The number of VLOAD, VSTORE. :  5,  5. (pinvr.f:25)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: pinvr.f

PROCEDURE NAME: PINVR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 43 [s0-s1 s8-s12 s15-s16 s18 s23-s28 s37-s63]
      Vector registers         : 26 [v38-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   8 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :   8 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (pinvr.f:23)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (pinvr.f:23)
    *** Estimated execution cycle                       : 23
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (pinvr.f:25)
      *** Estimated execution cycle                     : 460
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


