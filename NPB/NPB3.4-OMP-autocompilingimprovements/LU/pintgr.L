NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME : pintgr.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME: pintgr.f

PROCEDURE NAME: PINTGR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    33: par(1801): Parallel routine generated.: PINTGR$1
    37: par(1803): Parallelized by "do".
    37: opt(2495): Loop nest collapsed by "collapse" clause.
    37: vec( 101): Vectorized loop.
    42: vec( 128): Fused multiply-add operation applied.
    50: vec( 128): Fused multiply-add operation applied.
    60: par(1806): Single-construct generated.
    65: par(1807): Critical-construct generated.
    65: par(1809): Barrier synchronization.
    65: par(1803): Parallelized by "do".
    65: opt(2495): Loop nest collapsed by "collapse" clause.
    65: vec( 101): Vectorized loop.
    67: vec( 126): Idiom detected.: SUM
    80: par(1805): Master-construct generated.
    86: par(1803): Parallelized by "do".
    86: opt(2495): Loop nest collapsed by "collapse" clause.
    86: vec( 101): Vectorized loop.
    88: vec( 128): Fused multiply-add operation applied.
    98: par(1803): Parallelized by "do".
    98: opt(2495): Loop nest collapsed by "collapse" clause.
    98: vec( 101): Vectorized loop.
   100: vec( 128): Fused multiply-add operation applied.
   110: par(1806): Single-construct generated.
   115: par(1807): Critical-construct generated.
   115: par(1809): Barrier synchronization.
   115: par(1803): Parallelized by "do".
   115: opt(2495): Loop nest collapsed by "collapse" clause.
   115: vec( 101): Vectorized loop.
   117: vec( 126): Idiom detected.: SUM
   130: par(1805): Master-construct generated.
   136: par(1803): Parallelized by "do".
   136: opt(2495): Loop nest collapsed by "collapse" clause.
   136: vec( 101): Vectorized loop.
   138: vec( 128): Fused multiply-add operation applied.
   148: par(1803): Parallelized by "do".
   148: opt(2495): Loop nest collapsed by "collapse" clause.
   148: vec( 101): Vectorized loop.
   150: vec( 128): Fused multiply-add operation applied.
   160: par(1806): Single-construct generated.
   165: par(1807): Critical-construct generated.
   165: par(1809): Barrier synchronization.
   165: par(1803): Parallelized by "do".
   165: opt(2495): Loop nest collapsed by "collapse" clause.
   165: vec( 101): Vectorized loop.
   167: vec( 126): Idiom detected.: SUM
   180: par(1805): Master-construct generated.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME: pintgr.f

PROCEDURE NAME: PINTGR
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

     5:                 subroutine pintgr
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:                 use lu_data
    11:                 implicit none
    12:           
    13:           c---------------------------------------------------------------------
    14:           c  local variables
    15:           c---------------------------------------------------------------------
    16:                 integer i, j, k
    17:                 integer ibeg, ifin, ifin1
    18:                 integer jbeg, jfin, jfin1
    19:                 double precision frc1, frc2, frc3
    20:           
    21:           
    22:           
    23:           c---------------------------------------------------------------------
    24:           c   set up the sub-domains for integeration in each processor
    25:           c---------------------------------------------------------------------
    26:                 ibeg = ii1
    27:                 ifin = ii2
    28:                 jbeg = ji1
    29:                 jfin = ji2
    30:                 ifin1 = ifin - 1
    31:                 jfin1 = jfin - 1
    32:           
    33:           !$omp parallel default(shared) private(i,j,k)
    34:           !$omp&  shared(ki1,ki2,ifin,ibeg,jfin,jbeg,ifin1,jfin1)
    35:           
    36:           !$omp do schedule(static) collapse(2)
    37: Y------>        do j = jbeg,jfin
    38: |+----->           do i = ibeg,ifin
    39: ||        
    40: ||                    k = ki1
    41: ||        
    42: ||      C             phi1(i,j) = c2*(  u(5,i,j,k)
    43: ||             >           - 0.50d+00 * (  u(2,i,j,k) ** 2
    44: ||             >                         + u(3,i,j,k) ** 2
    45: ||             >                         + u(4,i,j,k) ** 2 )
    46: ||             >                        / u(1,i,j,k) )
    47: ||        
    48: ||                    k = ki2
    49: ||        
    50: ||      C             phi2(i,j) = c2*(  u(5,i,j,k)
    51: ||             >           - 0.50d+00 * (  u(2,i,j,k) ** 2
    52: ||             >                         + u(3,i,j,k) ** 2
    53: ||             >                         + u(4,i,j,k) ** 2 )
    54: ||             >                        / u(1,i,j,k) )
    55: |+-----            end do
    56: Y------         end do
    57:           !$omp end do nowait
    58:           
    59:           
    60:           !$omp single
    61:                 frc1 = 0.0d+00
    62:           !$omp end single
    63:           
    64:           !$omp do schedule(static) reduction(+:frc1) collapse(2)
    65: Y------>        do j = jbeg,jfin1
    66: |+----->           do i = ibeg, ifin1
    67: ||      G             frc1 = frc1 + (  phi1(i,j)
    68: ||             >                     + phi1(i+1,j)
    69: ||             >                     + phi1(i,j+1)
    70: ||             >                     + phi1(i+1,j+1)
    71: ||             >                     + phi2(i,j)
    72: ||             >                     + phi2(i+1,j)
    73: ||             >                     + phi2(i,j+1)
    74: ||             >                     + phi2(i+1,j+1) )
    75: |+-----            end do
    76: Y------         end do
    77:           !$omp end do
    78:           
    79:           
    80:           !$omp master
    81:                 frc1 = dxi * deta * frc1
    82:           !$omp end master
    83:           
    84:           
    85:           !$omp do schedule(static) collapse(2)
    86: Y------>        do k = ki1, ki2
    87: |+----->           do i = ibeg, ifin
    88: ||      C             phi1(i,k) = c2*(  u(5,i,jbeg,k)
    89: ||             >           - 0.50d+00 * (  u(2,i,jbeg,k) ** 2
    90: ||             >                         + u(3,i,jbeg,k) ** 2
    91: ||             >                         + u(4,i,jbeg,k) ** 2 )
    92: ||             >                        / u(1,i,jbeg,k) )
    93: |+-----            end do
    94: Y------         end do
    95:           !$omp end do nowait
    96:           
    97:           !$omp do schedule(static) collapse(2)
    98: Y------>        do k = ki1, ki2
    99: |+----->           do i = ibeg, ifin
   100: ||      C             phi2(i,k) = c2*(  u(5,i,jfin,k)
   101: ||             >           - 0.50d+00 * (  u(2,i,jfin,k) ** 2
   102: ||             >                         + u(3,i,jfin,k) ** 2
   103: ||             >                         + u(4,i,jfin,k) ** 2 )
   104: ||             >                        / u(1,i,jfin,k) )
   105: |+-----            end do
   106: Y------         end do
   107:           !$omp end do nowait
   108:           
   109:           
   110:           !$omp single
   111:                 frc2 = 0.0d+00
   112:           !$omp end single
   113:           
   114:           !$omp do schedule(static) reduction(+:frc2) collapse(2)
   115: Y------>        do k = ki1, ki2-1
   116: |+----->           do i = ibeg, ifin1
   117: ||      G             frc2 = frc2 + (  phi1(i,k)
   118: ||             >                     + phi1(i+1,k)
   119: ||             >                     + phi1(i,k+1)
   120: ||             >                     + phi1(i+1,k+1)
   121: ||             >                     + phi2(i,k)
   122: ||             >                     + phi2(i+1,k)
   123: ||             >                     + phi2(i,k+1)
   124: ||             >                     + phi2(i+1,k+1) )
   125: |+-----            end do
   126: Y------         end do
   127:           !$omp end do
   128:           
   129:           
   130:           !$omp master
   131:                 frc2 = dxi * dzeta * frc2
   132:           !$omp end master
   133:           
   134:           
   135:           !$omp do schedule(static) collapse(2)
   136: Y------>        do k = ki1, ki2
   137: |+----->           do j = jbeg, jfin
   138: ||      C             phi1(j,k) = c2*(  u(5,ibeg,j,k)
   139: ||             >           - 0.50d+00 * (  u(2,ibeg,j,k) ** 2
   140: ||             >                         + u(3,ibeg,j,k) ** 2
   141: ||             >                         + u(4,ibeg,j,k) ** 2 )
   142: ||             >                        / u(1,ibeg,j,k) )
   143: |+-----            end do
   144: Y------         end do
   145:           !$omp end do nowait
   146:           
   147:           !$omp do schedule(static) collapse(2)
   148: Y------>        do k = ki1, ki2
   149: |+----->           do j = jbeg, jfin
   150: ||      C             phi2(j,k) = c2*(  u(5,ifin,j,k)
   151: ||             >           - 0.50d+00 * (  u(2,ifin,j,k) ** 2
   152: ||             >                         + u(3,ifin,j,k) ** 2
   153: ||             >                         + u(4,ifin,j,k) ** 2 )
   154: ||             >                        / u(1,ifin,j,k) )
   155: |+-----            end do
   156: Y------         end do
   157:           !$omp end do nowait
   158:           
   159:           
   160:           !$omp single
   161:                 frc3 = 0.0d+00
   162:           !$omp end single
   163:           
   164:           !$omp do schedule(static) reduction(+:frc3) collapse(2)
   165: Y------>        do k = ki1, ki2-1
   166: |+----->           do j = jbeg, jfin1
   167: ||      G             frc3 = frc3 + (  phi1(j,k)
   168: ||             >                     + phi1(j+1,k)
   169: ||             >                     + phi1(j,k+1)
   170: ||             >                     + phi1(j+1,k+1)
   171: ||             >                     + phi2(j,k)
   172: ||             >                     + phi2(j+1,k)
   173: ||             >                     + phi2(j,k+1)
   174: ||             >                     + phi2(j+1,k+1) )
   175: |+-----            end do
   176: Y------         end do
   177:           !$omp end do
   178:           
   179:           
   180:           !$omp master
   181:                 frc3 = deta * dzeta * frc3
   182:           !$omp end master
   183:           !$omp end parallel
   184:           
   185:                 frc = 0.25d+00 * ( frc1 + frc2 + frc3 )
   186:           c      write (*,1001) frc
   187:           
   188:                 return
   189:           
   190:           c 1001 format (//5x,'surface integral = ',1pe12.5//)
   191:           
   192:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME: pintgr.f

PROCEDURE NAME: PINTGR
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME: pintgr.f

PROCEDURE NAME: PINTGR
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (pintgr.f:37)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    : 10,  2. (pintgr.f:37)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:37)
    *** VGT generated (pintgr.f:42)
    *** VGT generated (pintgr.f:50)
    *** VSC generated (pintgr.f:42)
    *** VSC generated (pintgr.f:50)
  LOOP END

  LOOP BEGIN: (pintgr.f:65)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  8,  0. (pintgr.f:65)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:65)
    *** VGT generated (pintgr.f:67)
    *** Idiom detected. : SUM (pintgr.f:67)
  LOOP END

  LOOP BEGIN: (pintgr.f:86)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  5,  1. (pintgr.f:86)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:86)
    *** VGT generated (pintgr.f:88)
    *** VSC generated (pintgr.f:88)
  LOOP END

  LOOP BEGIN: (pintgr.f:98)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  5,  1. (pintgr.f:98)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:98)
    *** VGT generated (pintgr.f:100)
    *** VSC generated (pintgr.f:100)
  LOOP END

  LOOP BEGIN: (pintgr.f:115)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  8,  0. (pintgr.f:115)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:115)
    *** VGT generated (pintgr.f:117)
    *** Idiom detected. : SUM (pintgr.f:117)
  LOOP END

  LOOP BEGIN: (pintgr.f:136)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  5,  1. (pintgr.f:136)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:136)
    *** VGT generated (pintgr.f:138)
    *** VSC generated (pintgr.f:138)
  LOOP END

  LOOP BEGIN: (pintgr.f:148)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  5,  1. (pintgr.f:148)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:148)
    *** VGT generated (pintgr.f:150)
    *** VSC generated (pintgr.f:150)
  LOOP END

  LOOP BEGIN: (pintgr.f:165)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  8,  0. (pintgr.f:165)
    *** The number of VLOAD, VSTORE. :  0,  0. (pintgr.f:165)
    *** VGT generated (pintgr.f:167)
    *** Idiom detected. : SUM (pintgr.f:167)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:54 2021
FILE NAME: pintgr.f

PROCEDURE NAME: PINTGR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 47 [s0-s2 s8-s12 s15-s16 s23-s34 s39-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 45424 bytes
      Register spill area      : 22576 bytes
      Parameter area           :    16 bytes
      Register save area       :   176 bytes
      User data area           : 22656 bytes
      Others                   :     0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (pintgr.f:37)
    *** Estimated execution cycle                       : 3292
  LOOP END

  LOOP BEGIN: (pintgr.f:65)
    *** Estimated execution cycle                       : 2212
  LOOP END

  LOOP BEGIN: (pintgr.f:86)
    *** Estimated execution cycle                       : 2055
  LOOP END

  LOOP BEGIN: (pintgr.f:98)
    *** Estimated execution cycle                       : 2055
  LOOP END

  LOOP BEGIN: (pintgr.f:115)
    *** Estimated execution cycle                       : 2212
  LOOP END

  LOOP BEGIN: (pintgr.f:136)
    *** Estimated execution cycle                       : 2055
  LOOP END

  LOOP BEGIN: (pintgr.f:148)
    *** Estimated execution cycle                       : 2055
  LOOP END

  LOOP BEGIN: (pintgr.f:165)
    *** Estimated execution cycle                       : 2212
  LOOP END


