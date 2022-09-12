NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME : utils.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    12: par(1801): Parallel routine generated.: RECIPROCAL$1
    13: par(1803): Parallelized by "do".
    13: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
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

     2:                 subroutine reciprocal (a, n)
     3:           c------------------------------------------------------------------
     4:           c     initialize double precision array a with length of n
     5:           c------------------------------------------------------------------
     6:           
     7:                 implicit none
     8:           
     9:                 integer n, i
    10:                 double precision a(n)
    11:           
    12:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    13: Y------>        do i = 1, n
    14: |                 a(i) = 1.d0/a(i)
    15: Y------         end do
    16:           c$OMP END PARALLEL DO
    17:           
    18:                 return
    19:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:13)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:13)
    *** The number of VLOAD, VSTORE. :  1,  1. (utils.f:13)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 20 [s0-s2 s8-s12 s15-s16 s23-s25 s57-s63]
      Vector registers         :  7 [v57-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:13)
    *** Estimated execution cycle                       : 393
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    31: par(1801): Parallel routine generated.: R_INIT_OMP$1
    32: par(1803): Parallelized by "do".
    32: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
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

    21:                 subroutine r_init_omp (a, n, const)
    22:           c------------------------------------------------------------------
    23:           c     initialize double precision array a with length of n
    24:           c------------------------------------------------------------------
    25:           
    26:                 implicit none
    27:           
    28:                 integer n, i
    29:                 double precision a(n), const
    30:           
    31:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    32: Y------>        do i = 1, n
    33: |                 a(i) = const
    34: Y------         end do
    35:           c$OMP END PARALLEL DO
    36:           
    37:                 return
    38:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:32)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:32)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:32)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s25 s56-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:32)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    50: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
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

    40:                 subroutine r_init (a, n, const)
    41:           c------------------------------------------------------------------
    42:           c     initialize double precision array a with length of n
    43:           c------------------------------------------------------------------
    44:           
    45:                 implicit none
    46:           
    47:                 integer n, i
    48:                 double precision a(n), const
    49:           
    50: V------>        do i = 1, n
    51: |                 a(i) = const
    52: V------         end do
    53:           
    54:                 return
    55:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:50)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:50)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:50)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0-s2 s8-s11 s58-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:50)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: par(1801): Parallel routine generated.: NR_INIT_OMP$1
    67: par(1803): Parallelized by "do".
    67: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
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

    57:                 subroutine nr_init_omp (a, n, const)
    58:           c------------------------------------------------------------------
    59:           c     initialize integer array a with length of n
    60:           c------------------------------------------------------------------
    61:           
    62:                 implicit none
    63:           
    64:                 integer n, i, a(n), const
    65:           
    66:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    67: Y------>        do i = 1, n
    68: |                 a(i) = const
    69: Y------         end do
    70:           c$OMP END PARALLEL DO
    71:           
    72:                 return
    73:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:67)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:67)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:67)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s25 s56-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:67)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    85: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
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

    76:                 subroutine nr_init (a, n, const)
    77:           c------------------------------------------------------------------
    78:           c     initialize integer array a with length of n
    79:           c------------------------------------------------------------------
    80:           
    81:                 implicit none
    82:           
    83:                 integer n, i, a(n), const
    84:           
    85: V------>        do i = 1, n
    86: |                 a(i) = const
    87: V------         end do
    88:           
    89:                 return
    90:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:85)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:85)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:85)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0-s2 s8-s11 s58-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:85)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   101: par(1801): Parallel routine generated.: L_INIT_OMP$1
   102: par(1803): Parallelized by "do".
   102: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
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

    92:                 subroutine l_init_omp (a, n, const)
    93:           c------------------------------------------------------------------
    94:           c     initialize integer array a with length of n
    95:           c------------------------------------------------------------------
    96:           
    97:                 implicit none
    98:                 integer n, i
    99:                 logical a(n), const
   100:           
   101:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
   102: Y------>        do i = 1, n
   103: |                 a(i) = const
   104: Y------         end do
   105:           c$OMP END PARALLEL DO
   106:           
   107:                 return
   108:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:102)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:102)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:102)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s25 s56-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:102)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   121: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
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

   111:                 subroutine ncopy (a,b,n)
   112:           c------------------------------------------------------------------
   113:           c     copy array of integers b to a, the length of array is n
   114:           c------------------------------------------------------------------
   115:           
   116:                 implicit none
   117:           
   118:                 integer n,i
   119:                 integer a(n), b(n)
   120:           
   121: V------>        do i = 1, n
   122: |                 a(i) = b(i)
   123: V------         end do
   124:           
   125:                 return
   126:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:121)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:121)
    *** The number of VLOAD, VSTORE. :  1,  1. (utils.f:121)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 15 [s0-s2 s8-s11 s56-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:121)
    *** Estimated execution cycle                       : 41
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   139: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
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

   129:                 subroutine copy (a,b,n)
   130:           c------------------------------------------------------------------
   131:           c     copy double precision array b to a, the length of array is n
   132:           c------------------------------------------------------------------
   133:           
   134:                 implicit none
   135:           
   136:                 integer n,i
   137:                 double precision a(n), b(n)
   138:           
   139: V------>        do i = 1, n
   140: |                  a(i) = b(i)
   141: V------         end do
   142:           
   143:                 return
   144:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:139)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:139)
    *** The number of VLOAD, VSTORE. :  1,  1. (utils.f:139)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 15 [s0-s2 s8-s11 s56-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:139)
    *** Estimated execution cycle                       : 41
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   155: par(1801): Parallel routine generated.: ADDS2M1$1
   156: par(1803): Parallelized by "do".
   156: vec( 101): Vectorized loop.
   157: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
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

   147:                 subroutine adds2m1(a,b,c1,n)
   148:           c-----------------------------------------------------------------
   149:           c     a=b*c1
   150:           c-----------------------------------------------------------------
   151:                 implicit none
   152:           
   153:                 integer n,i
   154:                 double precision a(n),b(n),c1
   155:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   156: Y------>        do i=1,n
   157: |       F         a(i)=a(i)+c1*b(i)
   158: Y------         end do
   159:           c$OMP END PARALLEL DO
   160:           
   161:                 return
   162:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:156)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:156)
    *** The number of VLOAD, VSTORE. :  2,  1. (utils.f:156)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 23 [s0-s2 s8-s12 s15-s16 s23-s25 s54-s63]
      Vector registers         :  3 [v61-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:156)
    *** Estimated execution cycle                       : 106
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   174: par(1801): Parallel routine generated.: ADDS1M1$1
   175: par(1803): Parallelized by "do".
   175: vec( 101): Vectorized loop.
   176: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
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

   165:                 subroutine adds1m1(a,b,c1,n )
   166:           c-----------------------------------------------------------------
   167:           c     a=c1*a+b
   168:           c-----------------------------------------------------------------
   169:           
   170:                 implicit none
   171:           
   172:                 integer n,i
   173:                 double precision a(n),b(n),c1
   174:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   175: Y------>        do i=1,n
   176: |       F         a(i)=c1*a(i)+b(i)
   177: Y------         end do
   178:           c$OMP END PARALLEL DO
   179:           
   180:                 return
   181:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:175)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:175)
    *** The number of VLOAD, VSTORE. :  2,  1. (utils.f:175)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 23 [s0-s2 s8-s12 s15-s16 s23-s25 s54-s63]
      Vector registers         :  3 [v61-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:175)
    *** Estimated execution cycle                       : 106
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   194: par(1801): Parallel routine generated.: COL2$1
   195: par(1803): Parallelized by "do".
   195: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
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

   184:                 subroutine col2(a,b,n)
   185:           c------------------------------------------------------------------
   186:           c     a=a*b
   187:           c------------------------------------------------------------------
   188:           
   189:                 implicit none
   190:           
   191:                 integer n,i
   192:                 double precision a(n),b(n)
   193:           
   194:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   195: Y------>        do i=1,n
   196: |                 a(i)=a(i)*b(i)
   197: Y------         end do
   198:           c$OMP END PARALLEL DO
   199:           
   200:                 return
   201:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:195)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:195)
    *** The number of VLOAD, VSTORE. :  2,  1. (utils.f:195)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s25 s56-s63]
      Vector registers         :  3 [v61-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:195)
    *** Estimated execution cycle                       : 107
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   213: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
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

   204:                 subroutine nrzero (na,n)
   205:           c------------------------------------------------------------------
   206:           c     zero out array of integers
   207:           c------------------------------------------------------------------
   208:           
   209:                 implicit none
   210:           
   211:                 integer n,i,na(n)
   212:           
   213: V------>        do i = 1, n
   214: |                 na(i ) = 0
   215: V------         end do
   216:           
   217:                 return
   218:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:213)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:213)
    *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:213)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 11 [s0-s1 s8-s11 s59-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:213)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   230: par(1801): Parallel routine generated.: ADD2$1
   231: par(1803): Parallelized by "do".
   231: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
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

   221:                 subroutine add2(a,b,n)
   222:           c------------------------------------------------------------------
   223:           c     a=a+b
   224:           c------------------------------------------------------------------
   225:           
   226:                 implicit none
   227:           
   228:                 integer n,i
   229:                 double precision  a(n),b(n)
   230:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   231: Y------>        do i=1,n
   232: |                 a(i)=a(i)+b(i)
   233: Y------         end do
   234:           c$OMP END PARALLEL DO
   235:           
   236:                 return
   237:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:231)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:231)
    *** The number of VLOAD, VSTORE. :  2,  1. (utils.f:231)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s25 s56-s63]
      Vector registers         :  3 [v61-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:231)
    *** Estimated execution cycle                       : 107
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   252: par(1801): Parallel routine generated.: CALC_NORM$1
   255: par(1803): Parallelized by "do".
   255: vec( 103): Unvectorized loop.
   255: vec( 113): Overhead of loop division is too large.
   258: vec( 101): Vectorized loop.
   259: opt(1418): Constant-length loop is expanded.
   260: opt(1418): Constant-length loop is expanded.
   261: vec( 126): Idiom detected.: SUM
   261: vec( 128): Fused multiply-add operation applied.
   268: par(1807): Critical-construct generated.
   268: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
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

   240:                 double precision function calc_norm()
   241:           c------------------------------------------------------------------
   242:           c     calculate the integral of ta1 over the whole domain
   243:           c------------------------------------------------------------------
   244:           
   245:                 use ua_data
   246:                 implicit none
   247:           
   248:                 double precision total,ieltotal
   249:                 integer iel,k,j,i,isize
   250:           
   251:                 total=0.d0
   252:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,isize,ieltotal,iel)
   253:           c$OMP& REDUCTION(+:total)
   254:           
   255: P------>        do iel=1,nelt
   256: |                 ieltotal=0.d0
   257: |                 isize=size_e(iel)
   258: |V----->          do k=1,lx1
   259: ||*---->            do j=1,lx1
   260: |||*--->              do i=1,lx1
   261: ||||    F               ieltotal=ieltotal+ta1(i,j,k,iel)*w3m1(i,j,k)
   262: ||||           &                               *jacm1_s(i,j,k,isize)
   263: |||*---               end do
   264: ||*----             end do
   265: |V-----           end do
   266: |               total=total+ieltotal
   267: P------         end do
   268:           c$OMP END PARALLEL DO
   269:           
   270:                 calc_norm = total
   271:           
   272:                 return
   273:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:255)
    <Unvectorized loop.>
    *** Overhead of loop division is too large. (utils.f:255)

    LOOP BEGIN: (utils.f:258)
      <Vectorized loop.>
      **  Runtime-tested short loop. (utils.f:258)
      *** The number of VGT,   VSC.    :  0,  0. (utils.f:258)
      *** The number of VLOAD, VSTORE. : 75,  0. (utils.f:258)
      *** Idiom detected. : SUM (utils.f:261)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 45296 bytes
      Register spill area      : 22536 bytes
      Parameter area           :     8 bytes
      Register save area       :   176 bytes
      User data area           : 22568 bytes
      Others                   :     8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:255)
    *** Estimated execution cycle                       : 1451
    *** The number of VECTOR REGISTER SPILL
          Total                                         : 1
            Not enough registers                        : 1
    *** The number of VECTOR REGISTER RESTORE
          Total                                         : 7
            Over basic blocks                           : 7
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   275: vec( 103): Unvectorized loop.
   275: vec( 108): Unvectorizable loop structure.
   293: vec( 109): Vectorization obstructive statement.
   293: vec( 103): Unvectorized loop.
   295: par(1801): Parallel routine generated.: PARALLEL_ADD$1
   296: par(1803): Parallelized by "do".
   298: vec( 101): Vectorized loop.
   308: par(1801): Parallel routine generated.: PARALLEL_ADD$2
   309: par(1803): Parallelized by "do".
   309: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
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

   275:                 subroutine parallel_add(frontier)
   276:           c-----------------------------------------------------------------
   277:           c     input array frontier, perform (potentially) parallel add so that
   278:           c     the output frontier(i) has sum of frontier(1)+frontier(2)+...+frontier(i)
   279:           c-----------------------------------------------------------------
   280:           
   281:                 use ua_data
   282:                 implicit none
   283:           
   284:                 integer nellog,i,ahead,ii,ntemp,n1,ntemp1,frontier(lelt),iel
   285:           
   286:                 nellog=0
   287:                 iel=1
   288:              10 iel=iel*2
   289:                 nellog=nellog+1
   290:                 if (iel.lt.nelt) goto 10
   291:           
   292:                 ntemp=1
   293: +------>        do i=1,nellog
   294: |                 n1=ntemp*2
   295: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ahead,ii,iel)
   296: |P----->          do iel=n1, nelt,n1
   297: ||                  ahead=frontier(iel-ntemp)
   298: ||V---->            do ii=ntemp-1,0,-1
   299: |||                   frontier(iel-ii)=frontier(iel-ii)+ahead
   300: ||V----             end do
   301: |P-----           end do
   302: |         c$OMP END PARALLEL DO
   303: |         
   304: |                 iel=(nelt/n1+1)*n1
   305: |                 ntemp1=iel-nelt
   306: |                 if(ntemp1.lt.ntemp)then
   307: |                   ahead=frontier(iel-ntemp)
   308: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)
   309: |Y----->            do ii=ntemp-1,ntemp1,-1
   310: ||                    frontier(iel-ii)=frontier(iel-ii)+ahead
   311: |Y-----             end do
   312: |         c$OMP END PARALLEL DO
   313: |                 end if
   314: |         
   315: |                 ntemp=n1
   316: +------         end do
   317:           
   318:                 return
   319:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:309)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (utils.f:309)
    *** The number of VLOAD, VSTORE. :  1,  1. (utils.f:309)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 23 [s0-s2 s8-s12 s15-s16 s23-s25 s54-s63]
      Vector registers         :  2 [v62-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 240 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  48 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:309)
    *** Estimated execution cycle                       : 72
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
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

   322:                 subroutine dssum
   323:           
   324:           c------------------------------------------------------------------
   325:           c     Perform stiffness summation: element-mortar-element mapping
   326:           c------------------------------------------------------------------
   327:           
   328:                 use ua_data
   329:                 implicit none
   330:           
   331:                 call transfb(dpcmor,dpcelm)
   332:                 call transf (dpcmor,dpcelm)
   333:           
   334:                 return
   335:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 14 [s0-s1 s8-s12 s15-s16 s23-s24 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   362: opt(1589): Outer loop moved inside inner loop(s).: IX
   362: vec( 101): Vectorized loop.
   364: opt(1592): Outer loop unrolled inside inner loop.: IZ


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
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

   338:                 subroutine facev(a,iface,val)
   339:           c------------------------------------------------------------------
   340:           c     assign the value val to face(iface,iel) of array a.
   341:           c------------------------------------------------------------------
   342:           
   343:                 use ua_data
   344:                 implicit none
   345:           
   346:                 double precision a(lx1,lx1,lx1), val
   347:                 integer iface, kx1, kx2, ky1, ky2, kz1, kz2, ix, iy, iz
   348:           
   349:                 kx1=1
   350:                 ky1=1
   351:                 kz1=1
   352:                 kx2=lx1
   353:                 ky2=lx1
   354:                 kz2=lx1
   355:                 if (iface.eq.1) kx1=lx1
   356:                 if (iface.eq.2) kx2=1
   357:                 if (iface.eq.3) ky1=lx1
   358:                 if (iface.eq.4) ky2=1
   359:                 if (iface.eq.5) kz1=lx1
   360:                 if (iface.eq.6) kz2=1
   361:           
   362: X------>        do ix = kx1, kx2
   363: |*----->          do iy = ky1, ky2
   364: ||U---->            do iz = kz1, kz2
   365: |||                   a(ix,iy,iz)=val
   366: ||U----             end do
   367: |*-----           end do
   368: X------         end do
   369:           
   370:                 return
   371:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (utils.f:363)
    <Unvectorized loop.>

    LOOP BEGIN: (utils.f:364)
      <Unvectorized loop.>
      **  Copied loop. (utils.f:364)

      LOOP BEGIN: (utils.f:362)
        <Vectorized loop.>
        **  Copied loop. (utils.f:362)
        **  Interchanged loop. (utils.f:362)
        **  Runtime-tested short loop. (utils.f:362)
        *** The number of VGT,   VSC.    :  0,  0. (utils.f:362)
        *** The number of VLOAD, VSTORE. :  0,  1. (utils.f:362)
      LOOP END
    LOOP END

    LOOP BEGIN: (utils.f:364)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (utils.f:364)

      LOOP BEGIN: (utils.f:362)
        <Vectorized loop.>
        **  Interchanged loop. (utils.f:362)
        **  Runtime-tested short loop. (utils.f:362)
        *** The number of VGT,   VSC.    :  0,  0. (utils.f:362)
        *** The number of VLOAD, VSTORE. :  0,  4. (utils.f:362)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 36 [s0-s11 s18-s21 s44-s63]
      Vector registers         :  5 [v59-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (utils.f:363)
    *** Estimated execution cycle                       : 29
    *** The number of SCALAR REGISTER TRANSFER          : 20

    LOOP BEGIN: (utils.f:364)
      *** Estimated execution cycle                     : 41
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END

    LOOP BEGIN: (utils.f:364)
      *** Estimated execution cycle                     : 137
      *** The number of SCALAR REGISTER TRANSFER        : 11
    LOOP END
  LOOP END


