NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME : error.f

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
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    21: opt(1418): Constant-length loop is expanded.
    25: par(1801): Parallel routine generated.: ERROR_NORM$1
    25: opt(1418): Constant-length loop is expanded.
    29: par(1803): Parallelized by "do".
    29: opt(2495): Loop nest collapsed by "collapse" clause.
    33: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    33: vec( 103): Unvectorized loop.
    35: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    37: opt(1418): Constant-length loop is expanded.
    44: par(1807): Critical-construct generated.
    44: par(1809): Barrier synchronization.
    44: opt(1418): Constant-length loop is expanded.
    46: opt(1418): Constant-length loop is expanded.
    47: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
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

     5:                  subroutine error_norm(rms)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function computes the norm of the difference between the
    12:           c computed solution and the exact solution
    13:           c---------------------------------------------------------------------
    14:           
    15:                  use sp_data
    16:                  implicit none
    17:           
    18:                  integer i, j, k, m, d
    19:                  double precision xi, eta, zeta, u_exact(5), rms(5), add
    20:           
    21: *------>         do m = 1, 5
    22: |                   rms(m) = 0.0d0
    23: *------          enddo
    24:           
    25:           !$omp parallel do default(shared)
    26:           !$omp&   private(i,j,k,m,zeta,eta,xi,add,u_exact)
    27:           !$omp&   reduction(+:rms)
    28:           !$omp&   schedule(static) collapse(2)
    29: P------>         do   k = 0, grid_points(3)-1
    30: |+----->            do   j = 0, grid_points(2)-1
    31: ||                     zeta = dble(k) * dnzm1
    32: ||                     eta = dble(j) * dnym1
    33: ||+---->               do   i = 0, grid_points(1)-1
    34: |||                       xi = dble(i) * dnxm1
    35: |||                       call exact_solution(xi, eta, zeta, u_exact)
    36: |||       
    37: |||*--->                  do   m = 1, 5
    38: ||||                         add = u(m,i,j,k)-u_exact(m)
    39: ||||                         rms(m) = rms(m) + add*add
    40: |||*---                   end do
    41: ||+----                end do
    42: |+-----             end do
    43: P------          end do
    44:           !$omp end parallel do
    45:           
    46: *------>         do    m = 1, 5
    47: |*----->            do    d = 1, 3
    48: ||                     rms(m) = rms(m) / dble(grid_points(d)-2)
    49: |*-----             end do
    50: |                   rms(m) = dsqrt(rms(m))
    51: *------          end do
    52:           
    53:                  return
    54:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
VECTORIZATION LIST

  LOOP BEGIN: (error.f:29)
    <Unvectorized loop.>

    LOOP BEGIN: (error.f:33)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (error.f:33)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 40 [s0-s3 s8-s12 s15-s16 s18-s33 s51-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 336 bytes
      Register spill area      :  24 bytes
      Parameter area           :  32 bytes
      Register save area       : 176 bytes
      User data area           : 104 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (error.f:29)
    *** Estimated execution cycle                       : 64
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 8

    LOOP BEGIN: (error.f:33)
      *** Estimated execution cycle                     : 99
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 17
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: opt(1418): Constant-length loop is expanded.
    70: par(1801): Parallel routine generated.: RHS_NORM$1
    70: opt(1418): Constant-length loop is expanded.
    73: opt(3008): Reference within a conditional branch moved outside loop.
    73: par(1803): Parallelized by "do".
    73: opt(2495): Loop nest collapsed by "collapse" clause.
    73: vec( 108): Unvectorizable loop structure.
    75: vec( 101): Vectorized loop.
    76: opt(1418): Constant-length loop is expanded.
    78: vec( 126): Idiom detected.: SUM
    78: vec( 128): Fused multiply-add operation applied.
    83: par(1807): Critical-construct generated.
    83: par(1809): Barrier synchronization.
    83: opt(1418): Constant-length loop is expanded.
    85: opt(1418): Constant-length loop is expanded.
    86: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
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

    58:                  subroutine rhs_norm(rms)
    59:           
    60:                  use sp_data
    61:                  implicit none
    62:           
    63:                  integer i, j, k, d, m
    64:                  double precision rms(5), add
    65:           
    66: *------>         do m = 1, 5
    67: |                   rms(m) = 0.0d0
    68: *------          enddo
    69:           
    70:           !$omp parallel do default(shared) private(i,j,k,m,add)
    71:           !$omp&   reduction(+:rms)
    72:           !$omp&   schedule(static) collapse(2)
    73: P------>         do k = 1, nz2
    74: |+----->            do j = 1, ny2
    75: ||V---->               do i = 1, nx2
    76: |||*--->                  do m = 1, 5
    77: ||||                         add = rhs(m,i,j,k)
    78: ||||    F                    rms(m) = rms(m) + add*add
    79: |||*---                   end do
    80: ||V----                end do
    81: |+-----             end do
    82: P------          end do
    83:           !$omp end parallel do
    84:           
    85: *------>         do   m = 1, 5
    86: |*----->            do   d = 1, 3
    87: ||                     rms(m) = rms(m) / dble(grid_points(d)-2)
    88: |*-----             end do
    89: |                   rms(m) = dsqrt(rms(m))
    90: *------          end do
    91:           
    92:                  return
    93:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (error.f:73)
    <Unvectorized loop.>
    *** Unvectorizable loop structure. (error.f:73)

    LOOP BEGIN: (error.f:75)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (error.f:75)
      *** The number of VLOAD, VSTORE. :  5,  0. (error.f:75)
      *** Idiom detected. : SUM (error.f:78)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:06 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 44 [s0-s2 s8-s12 s15-s16 s18-s19 s23-s29 s39-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (error.f:73)
    *** Estimated execution cycle                       : 352
    *** The number of SCALAR REGISTER TRANSFER          : 13

    LOOP BEGIN: (error.f:75)
      *** Estimated execution cycle                     : 328
    LOOP END
  LOOP END


