NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME : txinvr.f

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
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
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
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
        -fmove-loop-invariants-unsafe            : disable
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
        -freplace-loop-equation                  : disable
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
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    23: par(1801): Parallel routine generated.: TXINVR$1
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    26: vec( 101): Vectorized loop.
    28: vec( 101): Vectorized loop.
    43: vec( 128): Fused multiply-add operation applied.
    45: vec( 128): Fused multiply-add operation applied.
    49: vec( 128): Fused multiply-add operation applied.
    50: vec( 128): Fused multiply-add operation applied.
    51: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
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

     5:                  subroutine  txinvr
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c block-diagonal matrix-vector multiplication
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  integer i, j, k
    18:                  double precision t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3,
    19:                >                  r4, r5, ac2inv
    20:           
    21:           
    22:                  if (timeron) call timer_start(t_txinvr)
    23:           !$omp parallel do default(shared)
    24:           !$omp& private(i,j,k,t1,t2,t3,ac,ru1,uu,vv,ww,r1,r2,r3,r4,r5,ac2inv)
    25:           !$omp&  collapse(2)
    26: Y------>         do    k = 1, nz2
    27: |+----->            do    j = 1, ny2
    28: ||V---->               do    i = 1, nx2
    29: |||       
    30: |||                       ru1 = rho_i(i,j,k)
    31: |||                       uu = us(i,j,k)
    32: |||                       vv = vs(i,j,k)
    33: |||                       ww = ws(i,j,k)
    34: |||                       ac = speed(i,j,k)
    35: |||                       ac2inv = ac*ac
    36: |||       
    37: |||                       r1 = rhs(1,i,j,k)
    38: |||                       r2 = rhs(2,i,j,k)
    39: |||                       r3 = rhs(3,i,j,k)
    40: |||                       r4 = rhs(4,i,j,k)
    41: |||                       r5 = rhs(5,i,j,k)
    42: |||       
    43: |||     F                 t1 = c2 / ac2inv * ( qs(i,j,k)*r1 - uu*r2  -
    44: |||            >                  vv*r3 - ww*r4 + r5 )
    45: |||     F                 t2 = bt * ru1 * ( uu * r1 - r2 )
    46: |||                       t3 = ( bt * ru1 * ac ) * t1
    47: |||       
    48: |||                       rhs(1,i,j,k) = r1 - t1
    49: |||     F                 rhs(2,i,j,k) = - ru1 * ( ww*r1 - r4 )
    50: |||     F                 rhs(3,i,j,k) =   ru1 * ( vv*r1 - r3 )
    51: |||     F                 rhs(4,i,j,k) = - t2 + t3
    52: |||                       rhs(5,i,j,k) =   t2 + t3
    53: |||       
    54: ||V----                end do
    55: |+-----             end do
    56: Y------          end do
    57:                  if (timeron) call timer_stop(t_txinvr)
    58:           
    59:                  return
    60:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (txinvr.f:26)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (txinvr.f:26)
    *** The number of VLOAD, VSTORE. :  0,  2. (txinvr.f:26)

    LOOP BEGIN: (txinvr.f:28)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (txinvr.f:28)
      *** The number of VLOAD, VSTORE. : 11,  5. (txinvr.f:28)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 50 [v14-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 464 bytes
      Register spill area      : 136 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           : 136 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (txinvr.f:26)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (txinvr.f:26)
    *** Estimated execution cycle                       : 74
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 16
            Over basic blocks                           : 16
    *** The number of SCALAR REGISTER TRANSFER          : 7

    LOOP BEGIN: (txinvr.f:28)
      *** Estimated execution cycle                     : 1392
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


