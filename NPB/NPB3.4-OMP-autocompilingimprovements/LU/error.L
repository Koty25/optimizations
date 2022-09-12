NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME : error.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    26: opt(1418): Constant-length loop is expanded.
    30: par(1801): Parallel routine generated.: ERROR$1
    30: opt(1418): Constant-length loop is expanded.
    32: par(1803): Parallelized by "do".
    32: opt(2495): Loop nest collapsed by "collapse" clause.
    34: vec( 110): Vectorization obstructive procedure reference.: EXACT
    34: vec( 103): Unvectorized loop.
    35: inl(1212): Source for routine not found.: EXACT
    36: opt(1418): Constant-length loop is expanded.
    43: par(1807): Critical-construct generated.
    43: par(1809): Barrier synchronization.
    43: opt(1418): Constant-length loop is expanded.
    45: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR
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

     4:                 subroutine error
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c
    11:           c   compute the solution error
    12:           c
    13:           c---------------------------------------------------------------------
    14:           
    15:                 use lu_data
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  local variables
    20:           c---------------------------------------------------------------------
    21:                 integer i, j, k, m
    22:                 double precision  tmp
    23:                 double precision  u000ijk(5)
    24:           
    25:           
    26: *------>        do m = 1, 5
    27: |                  errnm(m) = 0.0d+00
    28: *------         end do
    29:           
    30:           !$omp parallel do schedule(static) collapse(2) default(shared)
    31:           !$omp&  private(i,j,k,m,tmp,u000ijk) reduction(+: errnm)
    32: P------>        do k = 2, nz-1
    33: |+----->           do j = jst, jend
    34: ||+---->              do i = ist, iend
    35: |||                      call exact( i, j, k, u000ijk )
    36: |||*--->                 do m = 1, 5
    37: ||||                        tmp = ( u000ijk(m) - u(m,i,j,k) )
    38: ||||                        errnm(m) = errnm(m) + tmp * tmp
    39: |||*---                  end do
    40: ||+----               end do
    41: |+-----            end do
    42: P------         end do
    43:           !$omp end parallel do
    44:           
    45: *------>        do m = 1, 5
    46: |                  errnm(m) = sqrt ( errnm(m) / ( dble(nx0-2)*(ny0-2)*(nz0-2) ) )
    47: *------         end do
    48:           
    49:           c        write (*,1002) ( errnm(m), m = 1, 5 )
    50:           
    51:            1002 format (1x/1x,'RMS-norm of error in soln. to ',
    52:                > 'first pde  = ',1pe12.5/,
    53:                > 1x,'RMS-norm of error in soln. to ',
    54:                > 'second pde = ',1pe12.5/,
    55:                > 1x,'RMS-norm of error in soln. to ',
    56:                > 'third pde  = ',1pe12.5/,
    57:                > 1x,'RMS-norm of error in soln. to ',
    58:                > 'fourth pde = ',1pe12.5/,
    59:                > 1x,'RMS-norm of error in soln. to ',
    60:                > 'fifth pde  = ',1pe12.5)
    61:           
    62:                 return
    63:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR
INLINE LIST

  ROOT: ERROR (error.f:4)
  -> NOINLINE: EXACT (error.f:35)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR
VECTORIZATION LIST

  LOOP BEGIN: (error.f:32)
    <Unvectorized loop.>

    LOOP BEGIN: (error.f:34)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT (error.f:34)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 37 [s0-s3 s8-s12 s15-s16 s18-s31 s52-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 304 bytes
      Register spill area      :   0 bytes
      Parameter area           :  32 bytes
      Register save area       : 176 bytes
      User data area           :  96 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (error.f:32)
    *** Estimated execution cycle                       : 55
    *** The number of SCALAR REGISTER TRANSFER          : 12

    LOOP BEGIN: (error.f:34)
      *** Estimated execution cycle                     : 106
      *** The number of SCALAR REGISTER TRANSFER        : 14
    LOOP END
  LOOP END


