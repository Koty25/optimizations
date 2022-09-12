NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME : work_lhs.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: vec( 101): Vectorized loop.
    38: opt(1418): Constant-length loop is expanded.
    46: opt(1112): Loop fused with previous loop.: M
    47: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
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

    24:                 subroutine lhsinit(lhs, ni)
    25:           
    26:           c---------------------------------------------------------------------
    27:           c---------------------------------------------------------------------
    28:            
    29:                 integer i, m, n, ni
    30:                 double precision lhs(5,5,3,0:ni)
    31:           
    32:           c---------------------------------------------------------------------
    33:           c     zero the whole left hand side for starters
    34:           c     set all diagonal values to 1. This is overkill, but convenient
    35:           c---------------------------------------------------------------------
    36:                 i = 0
    37: V------>        do m = 1, 5
    38: |*----->           do n = 1, 5
    39: ||                    lhs(m,n,1,i) = 0.0d0
    40: ||                    lhs(m,n,2,i) = 0.0d0
    41: ||                    lhs(m,n,3,i) = 0.0d0
    42: |*-----            end do
    43: |                  lhs(m,m,2,i) = 1.0d0
    44: |               end do
    45: |               i = ni
    46: |               do m = 1, 5
    47: |                  do n = 1, 5
    48: ||                    lhs(m,n,1,i) = 0.0d0
    49: ||                    lhs(m,n,2,i) = 0.0d0
    50: ||                    lhs(m,n,3,i) = 0.0d0
    51: |*-----            end do
    52: |                  lhs(m,m,2,i) = 1.0d0
    53: V------         end do
    54:           
    55:                 return
    56:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (work_lhs.f:37)
    <Vectorized loop.>
    **  Fused loop. (work_lhs.f:37)
    **  Runtime-tested short loop. (work_lhs.f:37)
    *** The number of VGT,   VSC.    :  0,  0. (work_lhs.f:37)
    *** The number of VLOAD, VSTORE. :  0, 32. (work_lhs.f:37)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:46 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 50 [s0-s11 s26-s63]
      Vector registers         : 32 [v32-v63]
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


