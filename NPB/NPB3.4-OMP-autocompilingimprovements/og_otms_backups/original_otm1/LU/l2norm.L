NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:01 2021
FILE NAME : l2norm.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:01 2021
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: opt(1418): Constant-length loop is expanded.
    41: par(1801): Parallel routine generated.: L2NORM$1
    41: opt(1418): Constant-length loop is expanded.
    43: opt(3008): Reference within a conditional branch moved outside loop.
    43: par(1803): Parallelized by "do".
    43: opt(2495): Loop nest collapsed by "collapse" clause.
    43: vec( 108): Unvectorizable loop structure.
    45: vec( 101): Vectorized loop.
    46: opt(1418): Constant-length loop is expanded.
    47: vec( 126): Idiom detected.: SUM
    47: vec( 128): Fused multiply-add operation applied.
    52: par(1807): Critical-construct generated.
    52: par(1809): Barrier synchronization.
    52: opt(1418): Constant-length loop is expanded.
    54: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:01 2021
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
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

     4:                 subroutine l2norm ( ldx, ldy, ldz,
     5:                >                    nx0, ny0, nz0,
     6:                >                    ist, iend,
     7:                >                    jst, jend,
     8:                >                    v, sum )
     9:           c---------------------------------------------------------------------
    10:           c---------------------------------------------------------------------
    11:           
    12:           c---------------------------------------------------------------------
    13:           c   to compute the l2-norm of vector v.
    14:           c---------------------------------------------------------------------
    15:           
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  input parameters
    20:           c---------------------------------------------------------------------
    21:                 integer ldx, ldy, ldz
    22:                 integer nx0, ny0, nz0
    23:                 integer ist, iend
    24:                 integer jst, jend
    25:           c---------------------------------------------------------------------
    26:           c   To improve cache performance, second two dimensions padded by 1
    27:           c   for even number sizes only.  Only needed in v.
    28:           c---------------------------------------------------------------------
    29:                 double precision  v(5,ldx/2*2+1,ldy/2*2+1,*), sum(5)
    30:           
    31:           c---------------------------------------------------------------------
    32:           c  local variables
    33:           c---------------------------------------------------------------------
    34:                 integer i, j, k, m
    35:           
    36:           
    37: *------>        do m = 1, 5
    38: |                  sum(m) = 0.0d+00
    39: *------         end do
    40:           
    41:           !$omp parallel do schedule(static) collapse(2) default(shared)
    42:           !$omp&  private(i,j,k,m) reduction(+: sum)
    43: P------>        do k = 2, nz0-1
    44: |+----->           do j = jst, jend
    45: ||V---->              do i = ist, iend
    46: |||*--->                 do m = 1, 5
    47: ||||    F                   sum(m) = sum(m) + v(m,i,j,k)*v(m,i,j,k)
    48: |||*---                  end do
    49: ||V----               end do
    50: |+-----            end do
    51: P------         end do
    52:           !$omp end parallel do
    53:           
    54: *------>        do m = 1, 5
    55: |                  sum(m) = sqrt ( sum(m) / ( dble(nx0-2)*(ny0-2)*(nz0-2) ) )
    56: *------         end do
    57:           
    58:                 return
    59:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:01 2021
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (l2norm.f:43)
    <Unvectorized loop.>
    *** Unvectorizable loop structure. (l2norm.f:43)

    LOOP BEGIN: (l2norm.f:45)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (l2norm.f:45)
      *** The number of VLOAD, VSTORE. :  5,  0. (l2norm.f:45)
      *** Idiom detected. : SUM (l2norm.f:47)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:01 2021
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 48 [s0-s2 s8-s12 s15-s16 s18-s19 s23-s32 s38-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 272 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  80 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (l2norm.f:43)
    *** Estimated execution cycle                       : 352
    *** The number of SCALAR REGISTER TRANSFER          : 13

    LOOP BEGIN: (l2norm.f:45)
      *** Estimated execution cycle                     : 328
    LOOP END
  LOOP END


