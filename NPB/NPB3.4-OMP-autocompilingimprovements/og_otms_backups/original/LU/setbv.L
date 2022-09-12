NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME : setbv.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    26: par(1801): Parallel routine generated.: SETBV$1
    29: vec( 110): Vectorization obstructive procedure reference.: EXACT
    29: par(1803): Parallelized by "do".
    29: opt(2495): Loop nest collapsed by "collapse" clause.
    29: vec( 103): Unvectorized loop.
    33: opt(1418): Constant-length loop is expanded.
    45: vec( 110): Vectorization obstructive procedure reference.: EXACT
    45: par(1803): Parallelized by "do".
    45: opt(2495): Loop nest collapsed by "collapse" clause.
    45: vec( 103): Unvectorized loop.
    49: opt(1418): Constant-length loop is expanded.
    61: vec( 110): Vectorization obstructive procedure reference.: EXACT
    61: par(1803): Parallelized by "do".
    61: opt(2495): Loop nest collapsed by "collapse" clause.
    61: vec( 103): Unvectorized loop.
    65: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
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

     5:                 subroutine setbv
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   set the boundary values of dependent variables
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use lu_data
    15:                 implicit none
    16:           
    17:           c---------------------------------------------------------------------
    18:           c   local variables
    19:           c---------------------------------------------------------------------
    20:                 integer i, j, k, m
    21:                 double precision temp1(5), temp2(5)
    22:           
    23:           c---------------------------------------------------------------------
    24:           c   set the dependent variable values along the top and bottom faces
    25:           c---------------------------------------------------------------------
    26:           !$omp parallel default(shared) private(i,j,k,m,temp1,temp2)
    27:           !$omp& shared(nx,ny,nz)
    28:           !$omp do schedule(static) collapse(2)
    29: P------>        do j = 1, ny
    30: |+----->           do i = 1, nx
    31: ||                    call exact( i, j, 1, temp1 )
    32: ||                    call exact( i, j, nz, temp2 )
    33: ||*---->              do m = 1, 5
    34: |||                      u( m, i, j, 1 ) = temp1(m)
    35: |||                      u( m, i, j, nz ) = temp2(m)
    36: ||*----               end do
    37: |+-----            end do
    38: P------         end do
    39:           !$omp end do
    40:           
    41:           c---------------------------------------------------------------------
    42:           c   set the dependent variable values along north and south faces
    43:           c---------------------------------------------------------------------
    44:           !$omp do schedule(static) collapse(2)
    45: P------>        do k = 1, nz
    46: |+----->           do i = 1, nx
    47: ||                    call exact( i, 1, k, temp1 )
    48: ||                    call exact( i, ny, k, temp2 )
    49: ||*---->              do m = 1, 5
    50: |||                      u( m, i, 1, k ) = temp1(m)
    51: |||                      u( m, i, ny, k ) = temp2(m)
    52: ||*----               end do
    53: |+-----            end do
    54: P------         end do
    55:           !$omp end do
    56:           
    57:           c---------------------------------------------------------------------
    58:           c   set the dependent variable values along east and west faces
    59:           c---------------------------------------------------------------------
    60:           !$omp do schedule(static) collapse(2)
    61: P------>        do k = 1, nz
    62: |+----->           do j = 1, ny
    63: ||                    call exact( 1, j, k, temp1 )
    64: ||                    call exact( nx, j, k, temp2 )
    65: ||*---->              do m = 1, 5
    66: |||                      u( m, 1, j, k ) = temp1(m)
    67: |||                      u( m, nx, j, k ) = temp2(m)
    68: ||*----               end do
    69: |+-----            end do
    70: P------         end do
    71:           !$omp end do nowait
    72:           !$omp end parallel
    73:           
    74:                 return
    75:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
VECTORIZATION LIST

  LOOP BEGIN: (setbv.f:29)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : EXACT (setbv.f:29)
  LOOP END

  LOOP BEGIN: (setbv.f:45)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : EXACT (setbv.f:45)
  LOOP END

  LOOP BEGIN: (setbv.f:61)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : EXACT (setbv.f:61)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 45 [s0-s3 s8-s12 s15-s16 s18-s37 s50-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 448 bytes
      Register spill area      :  64 bytes
      Parameter area           :  32 bytes
      Register save area       : 176 bytes
      User data area           : 176 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (setbv.f:29)
    *** Estimated execution cycle                       : 135
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER TRANSFER          : 14
  LOOP END

  LOOP BEGIN: (setbv.f:45)
    *** Estimated execution cycle                       : 137
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 4
            Across calls                                : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER TRANSFER          : 14
  LOOP END

  LOOP BEGIN: (setbv.f:61)
    *** Estimated execution cycle                       : 135
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Across calls                                : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER TRANSFER          : 11
  LOOP END


