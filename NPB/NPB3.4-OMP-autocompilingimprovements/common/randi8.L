NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME : randi8.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
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

     1:                 double precision function randlc(x, a)
     2:           
     3:           c---------------------------------------------------------------------
     4:           c
     5:           c   This routine returns a uniform pseudorandom double precision number in the
     6:           c   range (0, 1) by using the linear congruential generator
     7:           c
     8:           c   x_{k+1} = a x_k  (mod 2^46)
     9:           c
    10:           c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
    11:           c   before repeating.  The argument A is the same as 'a' in the above formula,
    12:           c   and X is the same as x_0.  A and X must be odd double precision integers
    13:           c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
    14:           c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
    15:           c   the new seed x_1, so that subsequent calls to RANDLC using the same
    16:           c   arguments will generate a continuous sequence.
    17:           
    18:                 implicit none
    19:                 double precision x, a
    20:                 integer*8 i246m1, Lx, La
    21:                 double precision d2m46
    22:           
    23:                 parameter(d2m46=0.5d0**46)
    24:           
    25:                 save i246m1
    26:                 data i246m1/X'00003FFFFFFFFFFF'/
    27:           
    28:                 Lx = X
    29:                 La = A
    30:           
    31:                 Lx   = iand(Lx*La,i246m1)
    32:                 randlc = d2m46*dble(Lx)
    33:                 x    = dble(Lx)
    34:                 return
    35:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         :  8 [s0-s1 s8-s11 s62-s63]
      Vector registers         :  0
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    71: vec( 102): Partially vectorized loop.
    72: opt(1019): Feedback of scalar value from one loop pass to another.: LX
    72: vec( 121): Unvectorizable dependency.
    74: opt(1512): Loop unrolled.: I


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
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

    42:                 SUBROUTINE VRANLC (N, X, A, Y)
    43:                 implicit none
    44:                 integer n, i
    45:                 double precision x, a, y(*)
    46:                 integer*8 i246m1, Lx, La
    47:                 double precision d2m46
    48:           
    49:           c This doesn't work, because the compiler does the calculation in 32
    50:           c bits and overflows. No standard way (without f90 stuff) to specify
    51:           c that the rhs should be done in 64 bit arithmetic.
    52:           c      parameter(i246m1=2**46-1)
    53:           
    54:                 parameter(d2m46=0.5d0**46)
    55:           
    56:                 save i246m1
    57:                 data i246m1/X'00003FFFFFFFFFFF'/
    58:           
    59:           c Note that the v6 compiler on an R8000 does something stupid with
    60:           c the above. Using the following instead (or various other things)
    61:           c makes the calculation run almost 10 times as fast.
    62:           c
    63:           c      save d2m46
    64:           c      data d2m46/0.0d0/
    65:           c      if (d2m46 .eq. 0.0d0) then
    66:           c         d2m46 = 0.5d0**46
    67:           c      endif
    68:           
    69:                 Lx = X
    70:                 La = A
    71: S------>        do i = 1, N
    72: |                  Lx   = iand(Lx*La,i246m1)
    73: |                  y(i) = d2m46*dble(Lx)
    74: S------         end do
    75:                 x    = dble(Lx)
    76:           
    77:                 return
    78:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (randi8.f:71)
    <Partially vectorized loop.>
    **  Unrolled loop. : THE UNROLL TIME IS 2. (randi8.f:71)
    *** The number of VGT,   VSC.    :  0,  0. (randi8.f:71)
    *** The number of VLOAD, VSTORE. :  1,  1. (randi8.f:71)
    *** Unvectorizable dependency. (randi8.f:72)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:48 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 29 [s0-s3 s8-s12 s15-s16 s18 s23-s25 s33 s51-s63]
      Vector registers         :  3 [v61-v63]
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

  LOOP BEGIN: (randi8.f:71)
    *** Estimated execution cycle                       : 18
  LOOP END

  LOOP BEGIN: (randi8.f:71)
    *** Estimated execution cycle                       : 106
  LOOP END


