NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME : y_solve.f

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
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    29: inl(1212): Source for routine not found.: TIMER_START
    30: par(1801): Parallel routine generated.: Y_SOLVE$1
    32: opt(3008): Reference within a conditional branch moved outside loop.
    32: vec( 110): Vectorization obstructive procedure reference.: LHSINIT
    32: par(1803): Parallelized by "do".
    32: opt(2495): Loop nest collapsed by "collapse" clause.
    35: inl(1212): Source for routine not found.: LHSINIT
    45: vec( 101): Vectorized loop.
    45: vec( 124): Iteration count is assumed. Iteration count.: 162
    48: vec( 128): Fused multiply-add operation applied.
    54: vec( 101): Vectorized loop.
    54: vec( 124): Iteration count is assumed. Iteration count.: 162
    56: vec( 128): Fused multiply-add operation applied.
    57: vec( 128): Fused multiply-add operation applied.
    58: vec( 128): Fused multiply-add operation applied.
    76: vec( 101): Vectorized loop.
    76: vec( 124): Iteration count is assumed. Iteration count.: 163
    97: vec( 101): Vectorized loop.
    97: vec( 124): Iteration count is assumed. Iteration count.: 163
    99: vec( 128): Fused multiply-add operation applied.
   102: vec( 128): Fused multiply-add operation applied.
   106: vec( 128): Fused multiply-add operation applied.
   109: vec( 128): Fused multiply-add operation applied.
   119: vec( 103): Unvectorized loop.
   119: vec( 113): Overhead of loop division is too large.
   122: opt(1037): Feedback of array elements.: LHS
   123: opt(1037): Feedback of array elements.: LHS
   125: opt(1418): Constant-length loop is expanded.
   126: opt(1037): Feedback of array elements.: RHS
   126: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   128: opt(1037): Feedback of array elements.: LHS
   128: vec( 120): Unvectorizable dependency.: LHS
   132: opt(1418): Constant-length loop is expanded.
   133: opt(1037): Feedback of array elements.: RHS
   133: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   136: vec( 120): Unvectorizable dependency.: LHS
   138: vec( 120): Unvectorizable dependency.: LHS
   140: opt(1418): Constant-length loop is expanded.
   141: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   157: opt(1418): Constant-length loop is expanded.
   164: opt(1418): Constant-length loop is expanded.
   172: opt(1418): Constant-length loop is expanded.
   179: vec( 103): Unvectorized loop.
   179: vec( 113): Overhead of loop division is too large.
   183: opt(1037): Feedback of array elements.: LHSP
   184: opt(1037): Feedback of array elements.: LHSP
   186: opt(1037): Feedback of array elements.: RHS
   186: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   187: opt(1037): Feedback of array elements.: LHSP
   187: vec( 120): Unvectorizable dependency.: LHSP
   191: opt(1037): Feedback of array elements.: RHS
   191: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   193: vec( 120): Unvectorizable dependency.: LHSP
   195: vec( 120): Unvectorizable dependency.: LHSP
   197: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   200: opt(1037): Feedback of array elements.: LHSM
   201: opt(1037): Feedback of array elements.: LHSM
   203: opt(1037): Feedback of array elements.: RHS
   203: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   204: opt(1037): Feedback of array elements.: LHSM
   204: vec( 120): Unvectorizable dependency.: LHSM
   208: opt(1037): Feedback of array elements.: RHS
   208: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   210: vec( 120): Unvectorizable dependency.: LHSM
   212: vec( 120): Unvectorizable dependency.: LHSM
   214: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   258: opt(1418): Constant-length loop is expanded.
   271: vec( 103): Unvectorized loop.
   271: vec( 113): Overhead of loop division is too large.
   274: opt(1418): Constant-length loop is expanded.
   275: opt(1037): Feedback of array elements.: RHS
   275: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   283: opt(1037): Feedback of array elements.: RHS
   283: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   286: opt(1037): Feedback of array elements.: RHS
   286: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   293: inl(1212): Source for routine not found.: TIMER_STOP
   296: inl(1212): Source for routine not found.: PINVR


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
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

     5:                  subroutine y_solve
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function performs the solution of the approximate factorization
    12:           c step in the y-direction for all five matrix components
    13:           c simultaneously. The Thomas algorithm is employed to solve the
    14:           c systems for the y-lines. Boundary conditions are non-periodic
    15:           c---------------------------------------------------------------------
    16:           
    17:                  use sp_data
    18:                  use work_lhs
    19:           
    20:                  implicit none
    21:           
    22:                  integer i, j, k, j1, j2, m
    23:                  double precision ru1, fac1, fac2
    24:           
    25:           
    26:           c---------------------------------------------------------------------
    27:           c---------------------------------------------------------------------
    28:           
    29:                  if (timeron) call timer_start(t_ysolve)
    30:           !$omp parallel do default(shared) private(i,j,k,j1,j2,m,
    31:           !$omp&    ru1,fac1,fac2) collapse(2)
    32: P------>         do  k = 1, nz2
    33: |+----->            do  i = 1, grid_points(1)-2
    34: ||        
    35: ||                    call lhsinit(ny2+1, lhs, lhsp, lhsm)
    36: ||        
    37: ||        c---------------------------------------------------------------------
    38: ||        c Computes the left hand side for the three y-factors
    39: ||        c---------------------------------------------------------------------
    40: ||        
    41: ||        c---------------------------------------------------------------------
    42: ||        c      first fill the lhs for the u-eigenvalue
    43: ||        c---------------------------------------------------------------------
    44: ||        
    45: ||V---->               do  j = 0, grid_points(2)-1
    46: |||                       ru1 = c3c4*rho_i(i,j,k)
    47: |||                       cv(j) = vs(i,j,k)
    48: |||     F                 rhov(j) = dmax1( dy3 + con43 * ru1,
    49: |||            >                           dy5 + c1c5*ru1,
    50: |||            >                           dymax + ru1,
    51: |||            >                           dy1)
    52: ||V----                end do
    53: ||         
    54: ||V---->               do  j = 1, grid_points(2)-2
    55: |||                       lhs(1,j) =  0.0d0
    56: |||     F                 lhs(2,j) = -dtty2 * cv(j-1) - dtty1 * rhov(j-1)
    57: |||     F                 lhs(3,j) =  1.0 + c2dtty1 * rhov(j)
    58: |||     F                 lhs(4,j) =  dtty2 * cv(j+1) - dtty1 * rhov(j+1)
    59: |||                       lhs(5,j) =  0.0d0
    60: ||V----                end do
    61: ||        
    62: ||        c---------------------------------------------------------------------
    63: ||        c      add fourth order dissipation
    64: ||        c---------------------------------------------------------------------
    65: ||        
    66: ||                     j = 1
    67: ||                     lhs(3,j) = lhs(3,j) + comz5
    68: ||                     lhs(4,j) = lhs(4,j) - comz4
    69: ||                     lhs(5,j) = lhs(5,j) + comz1
    70: ||         
    71: ||                     lhs(2,j+1) = lhs(2,j+1) - comz4
    72: ||                     lhs(3,j+1) = lhs(3,j+1) + comz6
    73: ||                     lhs(4,j+1) = lhs(4,j+1) - comz4
    74: ||                     lhs(5,j+1) = lhs(5,j+1) + comz1
    75: ||        
    76: ||V---->               do   j=3, grid_points(2)-4
    77: |||                       lhs(1,j) = lhs(1,j) + comz1
    78: |||                       lhs(2,j) = lhs(2,j) - comz4
    79: |||                       lhs(3,j) = lhs(3,j) + comz6
    80: |||                       lhs(4,j) = lhs(4,j) - comz4
    81: |||                       lhs(5,j) = lhs(5,j) + comz1
    82: ||V----                end do
    83: ||        
    84: ||                     j = grid_points(2)-3
    85: ||                     lhs(1,j) = lhs(1,j) + comz1
    86: ||                     lhs(2,j) = lhs(2,j) - comz4
    87: ||                     lhs(3,j) = lhs(3,j) + comz6
    88: ||                     lhs(4,j) = lhs(4,j) - comz4
    89: ||        
    90: ||                     lhs(1,j+1) = lhs(1,j+1) + comz1
    91: ||                     lhs(2,j+1) = lhs(2,j+1) - comz4
    92: ||                     lhs(3,j+1) = lhs(3,j+1) + comz5
    93: ||        
    94: ||        c---------------------------------------------------------------------
    95: ||        c      subsequently, do the other two factors
    96: ||        c---------------------------------------------------------------------
    97: ||V---->               do    j = 1, grid_points(2)-2
    98: |||                       lhsp(1,j) = lhs(1,j)
    99: |||     F                 lhsp(2,j) = lhs(2,j) -
   100: |||            >                            dtty2 * speed(i,j-1,k)
   101: |||                       lhsp(3,j) = lhs(3,j)
   102: |||     F                 lhsp(4,j) = lhs(4,j) +
   103: |||            >                            dtty2 * speed(i,j+1,k)
   104: |||                       lhsp(5,j) = lhs(5,j)
   105: |||                       lhsm(1,j) = lhs(1,j)
   106: |||     F                 lhsm(2,j) = lhs(2,j) +
   107: |||            >                            dtty2 * speed(i,j-1,k)
   108: |||                       lhsm(3,j) = lhs(3,j)
   109: |||     F                 lhsm(4,j) = lhs(4,j) -
   110: |||            >                            dtty2 * speed(i,j+1,k)
   111: |||                       lhsm(5,j) = lhs(5,j)
   112: ||V----                end do
   113: ||        
   114: ||        
   115: ||        c---------------------------------------------------------------------
   116: ||        c                          FORWARD ELIMINATION
   117: ||        c---------------------------------------------------------------------
   118: ||        
   119: ||+---->               do    j = 0, grid_points(2)-3
   120: |||                       j1 = j  + 1
   121: |||                       j2 = j  + 2
   122: |||                       fac1      = 1.d0/lhs(3,j)
   123: |||                       lhs(4,j)  = fac1*lhs(4,j)
   124: |||                       lhs(5,j)  = fac1*lhs(5,j)
   125: |||*--->                  do    m = 1, 3
   126: ||||                         rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   127: |||*---                   end do
   128: |||                       lhs(3,j1) = lhs(3,j1) -
   129: |||            >                         lhs(2,j1)*lhs(4,j)
   130: |||                       lhs(4,j1) = lhs(4,j1) -
   131: |||            >                         lhs(2,j1)*lhs(5,j)
   132: |||*--->                  do    m = 1, 3
   133: ||||                         rhs(m,i,j1,k) = rhs(m,i,j1,k) -
   134: ||||           >                         lhs(2,j1)*rhs(m,i,j,k)
   135: |||*---                   end do
   136: |||                       lhs(2,j2) = lhs(2,j2) -
   137: |||            >                         lhs(1,j2)*lhs(4,j)
   138: |||                       lhs(3,j2) = lhs(3,j2) -
   139: |||            >                         lhs(1,j2)*lhs(5,j)
   140: |||*--->                  do    m = 1, 3
   141: ||||                         rhs(m,i,j2,k) = rhs(m,i,j2,k) -
   142: ||||           >                         lhs(1,j2)*rhs(m,i,j,k)
   143: |||*---                   end do
   144: ||+----                end do
   145: ||        
   146: ||        c---------------------------------------------------------------------
   147: ||        c      The last two rows in this grid block are a bit different,
   148: ||        c      since they do not have two more rows available for the
   149: ||        c      elimination of off-diagonal entries
   150: ||        c---------------------------------------------------------------------
   151: ||        
   152: ||                     j  = grid_points(2)-2
   153: ||                     j1 = grid_points(2)-1
   154: ||                     fac1      = 1.d0/lhs(3,j)
   155: ||                     lhs(4,j)  = fac1*lhs(4,j)
   156: ||                     lhs(5,j)  = fac1*lhs(5,j)
   157: ||*---->               do    m = 1, 3
   158: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   159: ||*----                end do
   160: ||                     lhs(3,j1) = lhs(3,j1) -
   161: ||             >                      lhs(2,j1)*lhs(4,j)
   162: ||                     lhs(4,j1) = lhs(4,j1) -
   163: ||             >                      lhs(2,j1)*lhs(5,j)
   164: ||*---->               do    m = 1, 3
   165: |||                       rhs(m,i,j1,k) = rhs(m,i,j1,k) -
   166: |||            >                      lhs(2,j1)*rhs(m,i,j,k)
   167: ||*----                end do
   168: ||        c---------------------------------------------------------------------
   169: ||        c            scale the last row immediately
   170: ||        c---------------------------------------------------------------------
   171: ||                     fac2      = 1.d0/lhs(3,j1)
   172: ||*---->               do    m = 1, 3
   173: |||                       rhs(m,i,j1,k) = fac2*rhs(m,i,j1,k)
   174: ||*----                end do
   175: ||        
   176: ||        c---------------------------------------------------------------------
   177: ||        c      do the u+c and the u-c factors
   178: ||        c---------------------------------------------------------------------
   179: ||+---->               do    j = 0, grid_points(2)-3
   180: |||                       j1 = j  + 1
   181: |||                       j2 = j  + 2
   182: |||                       m = 4
   183: |||                       fac1       = 1.d0/lhsp(3,j)
   184: |||                       lhsp(4,j)  = fac1*lhsp(4,j)
   185: |||                       lhsp(5,j)  = fac1*lhsp(5,j)
   186: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   187: |||                       lhsp(3,j1) = lhsp(3,j1) -
   188: |||            >                       lhsp(2,j1)*lhsp(4,j)
   189: |||                       lhsp(4,j1) = lhsp(4,j1) -
   190: |||            >                       lhsp(2,j1)*lhsp(5,j)
   191: |||                       rhs(m,i,j1,k) = rhs(m,i,j1,k) -
   192: |||            >                       lhsp(2,j1)*rhs(m,i,j,k)
   193: |||                       lhsp(2,j2) = lhsp(2,j2) -
   194: |||            >                       lhsp(1,j2)*lhsp(4,j)
   195: |||                       lhsp(3,j2) = lhsp(3,j2) -
   196: |||            >                       lhsp(1,j2)*lhsp(5,j)
   197: |||                       rhs(m,i,j2,k) = rhs(m,i,j2,k) -
   198: |||            >                       lhsp(1,j2)*rhs(m,i,j,k)
   199: |||                       m = 5
   200: |||                       fac1       = 1.d0/lhsm(3,j)
   201: |||                       lhsm(4,j)  = fac1*lhsm(4,j)
   202: |||                       lhsm(5,j)  = fac1*lhsm(5,j)
   203: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   204: |||                       lhsm(3,j1) = lhsm(3,j1) -
   205: |||            >                       lhsm(2,j1)*lhsm(4,j)
   206: |||                       lhsm(4,j1) = lhsm(4,j1) -
   207: |||            >                       lhsm(2,j1)*lhsm(5,j)
   208: |||                       rhs(m,i,j1,k) = rhs(m,i,j1,k) -
   209: |||            >                       lhsm(2,j1)*rhs(m,i,j,k)
   210: |||                       lhsm(2,j2) = lhsm(2,j2) -
   211: |||            >                       lhsm(1,j2)*lhsm(4,j)
   212: |||                       lhsm(3,j2) = lhsm(3,j2) -
   213: |||            >                       lhsm(1,j2)*lhsm(5,j)
   214: |||                       rhs(m,i,j2,k) = rhs(m,i,j2,k) -
   215: |||            >                       lhsm(1,j2)*rhs(m,i,j,k)
   216: ||+----                end do
   217: ||        
   218: ||        c---------------------------------------------------------------------
   219: ||        c         And again the last two rows separately
   220: ||        c---------------------------------------------------------------------
   221: ||                     j  = grid_points(2)-2
   222: ||                     j1 = grid_points(2)-1
   223: ||                     m = 4
   224: ||                     fac1       = 1.d0/lhsp(3,j)
   225: ||                     lhsp(4,j)  = fac1*lhsp(4,j)
   226: ||                     lhsp(5,j)  = fac1*lhsp(5,j)
   227: ||                     rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   228: ||                     lhsp(3,j1) = lhsp(3,j1) -
   229: ||             >                    lhsp(2,j1)*lhsp(4,j)
   230: ||                     lhsp(4,j1) = lhsp(4,j1) -
   231: ||             >                    lhsp(2,j1)*lhsp(5,j)
   232: ||                     rhs(m,i,j1,k)   = rhs(m,i,j1,k) -
   233: ||             >                    lhsp(2,j1)*rhs(m,i,j,k)
   234: ||                     m = 5
   235: ||                     fac1       = 1.d0/lhsm(3,j)
   236: ||                     lhsm(4,j)  = fac1*lhsm(4,j)
   237: ||                     lhsm(5,j)  = fac1*lhsm(5,j)
   238: ||                     rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   239: ||                     lhsm(3,j1) = lhsm(3,j1) -
   240: ||             >                    lhsm(2,j1)*lhsm(4,j)
   241: ||                     lhsm(4,j1) = lhsm(4,j1) -
   242: ||             >                    lhsm(2,j1)*lhsm(5,j)
   243: ||                     rhs(m,i,j1,k)   = rhs(m,i,j1,k) -
   244: ||             >                    lhsm(2,j1)*rhs(m,i,j,k)
   245: ||        c---------------------------------------------------------------------
   246: ||        c               Scale the last row immediately
   247: ||        c---------------------------------------------------------------------
   248: ||                     rhs(4,i,j1,k)   = rhs(4,i,j1,k)/lhsp(3,j1)
   249: ||                     rhs(5,i,j1,k)   = rhs(5,i,j1,k)/lhsm(3,j1)
   250: ||        
   251: ||        
   252: ||        c---------------------------------------------------------------------
   253: ||        c                         BACKSUBSTITUTION
   254: ||        c---------------------------------------------------------------------
   255: ||        
   256: ||                     j  = grid_points(2)-2
   257: ||                     j1 = grid_points(2)-1
   258: ||*---->               do   m = 1, 3
   259: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) -
   260: |||            >                           lhs(4,j)*rhs(m,i,j1,k)
   261: ||*----                end do
   262: ||        
   263: ||                     rhs(4,i,j,k) = rhs(4,i,j,k) -
   264: ||             >                           lhsp(4,j)*rhs(4,i,j1,k)
   265: ||                     rhs(5,i,j,k) = rhs(5,i,j,k) -
   266: ||             >                           lhsm(4,j)*rhs(5,i,j1,k)
   267: ||        
   268: ||        c---------------------------------------------------------------------
   269: ||        c      The first three factors
   270: ||        c---------------------------------------------------------------------
   271: ||+---->               do   j = grid_points(2)-3, 0, -1
   272: |||                       j1 = j  + 1
   273: |||                       j2 = j  + 2
   274: |||*--->                  do   m = 1, 3
   275: ||||                         rhs(m,i,j,k) = rhs(m,i,j,k) -
   276: ||||           >                          lhs(4,j)*rhs(m,i,j1,k) -
   277: ||||           >                          lhs(5,j)*rhs(m,i,j2,k)
   278: |||*---                   end do
   279: |||       
   280: |||       c---------------------------------------------------------------------
   281: |||       c      And the remaining two
   282: |||       c---------------------------------------------------------------------
   283: |||                       rhs(4,i,j,k) = rhs(4,i,j,k) -
   284: |||            >                          lhsp(4,j)*rhs(4,i,j1,k) -
   285: |||            >                          lhsp(5,j)*rhs(4,i,j2,k)
   286: |||                       rhs(5,i,j,k) = rhs(5,i,j,k) -
   287: |||            >                          lhsm(4,j)*rhs(5,i,j1,k) -
   288: |||            >                          lhsm(5,j)*rhs(5,i,j2,k)
   289: ||+----                end do
   290: ||        
   291: |+-----             end do
   292: P------          end do
   293:                  if (timeron) call timer_stop(t_ysolve)
   294:           
   295:           
   296:                  call pinvr
   297:           
   298:                  return
   299:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
INLINE LIST

  ROOT: Y_SOLVE (y_solve.f:5)
  -> NOINLINE: TIMER_START (y_solve.f:29)
     *** Source for routine not found.
  -> NOINLINE: LHSINIT (y_solve.f:35)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (y_solve.f:293)
     *** Source for routine not found.
  -> NOINLINE: PINVR (y_solve.f:296)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (y_solve.f:32)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : LHSINIT (y_solve.f:32)

    LOOP BEGIN: (y_solve.f:45)
      <Vectorized loop.>
      **  Runtime-tested short loop. (y_solve.f:45)
      *** The number of VGT,   VSC.    :  0,  0. (y_solve.f:45)
      *** The number of VLOAD, VSTORE. :  2,  2. (y_solve.f:45)
      *** Iteration count is assumed. Iteration count. : 162 (y_solve.f:45)
    LOOP END

    LOOP BEGIN: (y_solve.f:54)
      <Vectorized loop.>
      **  Runtime-tested short loop. (y_solve.f:54)
      *** The number of VGT,   VSC.    :  0,  0. (y_solve.f:54)
      *** The number of VLOAD, VSTORE. :  5,  5. (y_solve.f:54)
      *** Iteration count is assumed. Iteration count. : 162 (y_solve.f:54)
    LOOP END

    LOOP BEGIN: (y_solve.f:76)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (y_solve.f:76)
      *** The number of VLOAD, VSTORE. :  5,  5. (y_solve.f:76)
      *** Iteration count is assumed. Iteration count. : 163 (y_solve.f:76)
    LOOP END

    LOOP BEGIN: (y_solve.f:97)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (y_solve.f:97)
      *** The number of VLOAD, VSTORE. :  7, 10. (y_solve.f:97)
      *** Iteration count is assumed. Iteration count. : 163 (y_solve.f:97)
    LOOP END

    LOOP BEGIN: (y_solve.f:119)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (y_solve.f:119)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:141)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:133)
      *** Unvectorizable dependency. : LHS (y_solve.f:138)
      *** Unvectorizable dependency. : LHS (y_solve.f:136)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:126)
      *** Unvectorizable dependency. : LHS (y_solve.f:128)
    LOOP END

    LOOP BEGIN: (y_solve.f:179)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (y_solve.f:179)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:214)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:208)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:203)
      *** Unvectorizable dependency. : LHSP (y_solve.f:195)
      *** Unvectorizable dependency. : LHSM (y_solve.f:212)
      *** Unvectorizable dependency. : LHSP (y_solve.f:193)
      *** Unvectorizable dependency. : LHSM (y_solve.f:210)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:191)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:186)
      *** Unvectorizable dependency. : LHSP (y_solve.f:187)
      *** Unvectorizable dependency. : LHSM (y_solve.f:204)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:197)
    LOOP END

    LOOP BEGIN: (y_solve.f:271)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (y_solve.f:271)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:286)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:283)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (y_solve.f:275)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 60 [s0-s5 s7-s12 s15-s16 s18-s63]
      Vector registers         : 46 [v18-v63]
      Vector mask registers    :  1 [vm15]
      VREG directive           :  0

  Routine stack
    Total size                 : 37712 bytes
      Register spill area      : 18744 bytes
      Parameter area           :    32 bytes
      Register save area       :   176 bytes
      User data area           : 18760 bytes
      Others                   :     0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (y_solve.f:32)
    *** Estimated execution cycle                       : 492
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 44
            Over basic blocks                           : 44
    *** The number of SCALAR REGISTER TRANSFER          : 17

    LOOP BEGIN: (y_solve.f:119)
      *** Estimated execution cycle                     : 146
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 2
              Over basic blocks                         : 2
    LOOP END

    LOOP BEGIN: (y_solve.f:179)
      *** Estimated execution cycle                     : 189
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
    LOOP END

    LOOP BEGIN: (y_solve.f:271)
      *** Estimated execution cycle                     : 90
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END


