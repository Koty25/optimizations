NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME : z_solve.f

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
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    33: inl(1212): Source for routine not found.: TIMER_START
    34: par(1801): Parallel routine generated.: Z_SOLVE$1
    36: opt(3008): Reference within a conditional branch moved outside loop.
    36: vec( 110): Vectorization obstructive procedure reference.: LHSINIT
    36: par(1803): Parallelized by "do".
    36: opt(2495): Loop nest collapsed by "collapse" clause.
    39: inl(1212): Source for routine not found.: LHSINIT
    49: vec( 101): Vectorized loop.
    49: vec( 124): Iteration count is assumed. Iteration count.: 162
    52: vec( 128): Fused multiply-add operation applied.
    58: vec( 101): Vectorized loop.
    58: vec( 124): Iteration count is assumed. Iteration count.: 162
    60: vec( 128): Fused multiply-add operation applied.
    61: vec( 128): Fused multiply-add operation applied.
    62: vec( 128): Fused multiply-add operation applied.
    81: vec( 101): Vectorized loop.
    81: vec( 124): Iteration count is assumed. Iteration count.: 163
   104: vec( 101): Vectorized loop.
   104: vec( 124): Iteration count is assumed. Iteration count.: 163
   106: vec( 128): Fused multiply-add operation applied.
   109: vec( 128): Fused multiply-add operation applied.
   113: vec( 128): Fused multiply-add operation applied.
   116: vec( 128): Fused multiply-add operation applied.
   126: vec( 103): Unvectorized loop.
   126: vec( 113): Overhead of loop division is too large.
   129: opt(1037): Feedback of array elements.: LHS
   130: opt(1037): Feedback of array elements.: LHS
   132: opt(1418): Constant-length loop is expanded.
   133: opt(1037): Feedback of array elements.: RHS
   133: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   135: opt(1037): Feedback of array elements.: LHS
   135: vec( 120): Unvectorizable dependency.: LHS
   139: opt(1418): Constant-length loop is expanded.
   140: opt(1037): Feedback of array elements.: RHS
   140: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   143: vec( 120): Unvectorizable dependency.: LHS
   145: vec( 120): Unvectorizable dependency.: LHS
   147: opt(1418): Constant-length loop is expanded.
   148: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   163: opt(1418): Constant-length loop is expanded.
   170: opt(1418): Constant-length loop is expanded.
   178: opt(1418): Constant-length loop is expanded.
   185: vec( 103): Unvectorized loop.
   185: vec( 113): Overhead of loop division is too large.
   189: opt(1037): Feedback of array elements.: LHSP
   190: opt(1037): Feedback of array elements.: LHSP
   192: opt(1037): Feedback of array elements.: RHS
   192: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   193: opt(1037): Feedback of array elements.: LHSP
   193: vec( 120): Unvectorizable dependency.: LHSP
   197: opt(1037): Feedback of array elements.: RHS
   197: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   199: vec( 120): Unvectorizable dependency.: LHSP
   201: vec( 120): Unvectorizable dependency.: LHSP
   203: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   206: opt(1037): Feedback of array elements.: LHSM
   207: opt(1037): Feedback of array elements.: LHSM
   209: opt(1037): Feedback of array elements.: RHS
   209: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   210: opt(1037): Feedback of array elements.: LHSM
   210: vec( 120): Unvectorizable dependency.: LHSM
   214: opt(1037): Feedback of array elements.: RHS
   214: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   216: vec( 120): Unvectorizable dependency.: LHSM
   218: vec( 120): Unvectorizable dependency.: LHSM
   220: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   265: opt(1418): Constant-length loop is expanded.
   283: vec( 103): Unvectorized loop.
   283: vec( 113): Overhead of loop division is too large.
   286: opt(1418): Constant-length loop is expanded.
   287: opt(1037): Feedback of array elements.: RHS
   287: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   295: opt(1037): Feedback of array elements.: RHS
   295: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   298: opt(1037): Feedback of array elements.: RHS
   298: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   305: inl(1212): Source for routine not found.: TIMER_STOP
   307: inl(1212): Source for routine not found.: TZETAR


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
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

     5:                  subroutine z_solve
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function performs the solution of the approximate factorization
    12:           c step in the z-direction for all five matrix components
    13:           c simultaneously. The Thomas algorithm is employed to solve the
    14:           c systems for the z-lines. Boundary conditions are non-periodic
    15:           c---------------------------------------------------------------------
    16:           
    17:                  use sp_data
    18:                  use work_lhs
    19:           
    20:                  implicit none
    21:           
    22:                  integer i, j, k, k1, k2, m
    23:                  double precision ru1, fac1, fac2
    24:           
    25:           
    26:           c---------------------------------------------------------------------
    27:           c---------------------------------------------------------------------
    28:           
    29:           c---------------------------------------------------------------------
    30:           c Prepare for z-solve, array redistribution
    31:           c---------------------------------------------------------------------
    32:           
    33:                  if (timeron) call timer_start(t_zsolve)
    34:           !$omp parallel do default(shared) private(i,j,k,k1,k2,m,
    35:           !$omp&    ru1,fac1,fac2) collapse(2)
    36: P------>         do   j = 1, ny2
    37: |+----->            do   i = 1, nx2
    38: ||        
    39: ||                    call lhsinit(nz2+1, lhs, lhsp, lhsm)
    40: ||        
    41: ||        c---------------------------------------------------------------------
    42: ||        c Computes the left hand side for the three z-factors
    43: ||        c---------------------------------------------------------------------
    44: ||        
    45: ||        c---------------------------------------------------------------------
    46: ||        c first fill the lhs for the u-eigenvalue
    47: ||        c---------------------------------------------------------------------
    48: ||        
    49: ||V---->               do   k = 0, nz2 + 1
    50: |||                       ru1 = c3c4*rho_i(i,j,k)
    51: |||                       cv(k) = ws(i,j,k)
    52: |||     F                 rhov(k) = dmax1(dz4 + con43 * ru1,
    53: |||            >                          dz5 + c1c5 * ru1,
    54: |||            >                          dzmax + ru1,
    55: |||            >                          dz1)
    56: ||V----                end do
    57: ||        
    58: ||V---->               do   k =  1, nz2
    59: |||                       lhs(1,k) =  0.0d0
    60: |||     F                 lhs(2,k) = -dttz2 * cv(k-1) - dttz1 * rhov(k-1)
    61: |||     F                 lhs(3,k) =  1.0 + c2dttz1 * rhov(k)
    62: |||     F                 lhs(4,k) =  dttz2 * cv(k+1) - dttz1 * rhov(k+1)
    63: |||                       lhs(5,k) =  0.0d0
    64: ||V----                end do
    65: ||        
    66: ||        c---------------------------------------------------------------------
    67: ||        c      add fourth order dissipation
    68: ||        c---------------------------------------------------------------------
    69: ||        
    70: ||                     k = 1
    71: ||                     lhs(3,k) = lhs(3,k) + comz5
    72: ||                     lhs(4,k) = lhs(4,k) - comz4
    73: ||                     lhs(5,k) = lhs(5,k) + comz1
    74: ||        
    75: ||                     k = 2
    76: ||                     lhs(2,k) = lhs(2,k) - comz4
    77: ||                     lhs(3,k) = lhs(3,k) + comz6
    78: ||                     lhs(4,k) = lhs(4,k) - comz4
    79: ||                     lhs(5,k) = lhs(5,k) + comz1
    80: ||        
    81: ||V---->               do    k = 3, nz2-2
    82: |||                       lhs(1,k) = lhs(1,k) + comz1
    83: |||                       lhs(2,k) = lhs(2,k) - comz4
    84: |||                       lhs(3,k) = lhs(3,k) + comz6
    85: |||                       lhs(4,k) = lhs(4,k) - comz4
    86: |||                       lhs(5,k) = lhs(5,k) + comz1
    87: ||V----                end do
    88: ||        
    89: ||                     k = nz2-1
    90: ||                     lhs(1,k) = lhs(1,k) + comz1
    91: ||                     lhs(2,k) = lhs(2,k) - comz4
    92: ||                     lhs(3,k) = lhs(3,k) + comz6
    93: ||                     lhs(4,k) = lhs(4,k) - comz4
    94: ||        
    95: ||                     k = nz2
    96: ||                     lhs(1,k) = lhs(1,k) + comz1
    97: ||                     lhs(2,k) = lhs(2,k) - comz4
    98: ||                     lhs(3,k) = lhs(3,k) + comz5
    99: ||        
   100: ||        
   101: ||        c---------------------------------------------------------------------
   102: ||        c      subsequently, fill the other factors (u+c), (u-c)
   103: ||        c---------------------------------------------------------------------
   104: ||V---->               do    k = 1, nz2
   105: |||                       lhsp(1,k) = lhs(1,k)
   106: |||     F                 lhsp(2,k) = lhs(2,k) -
   107: |||            >                            dttz2 * speed(i,j,k-1)
   108: |||                       lhsp(3,k) = lhs(3,k)
   109: |||     F                 lhsp(4,k) = lhs(4,k) +
   110: |||            >                            dttz2 * speed(i,j,k+1)
   111: |||                       lhsp(5,k) = lhs(5,k)
   112: |||                       lhsm(1,k) = lhs(1,k)
   113: |||     F                 lhsm(2,k) = lhs(2,k) +
   114: |||            >                            dttz2 * speed(i,j,k-1)
   115: |||                       lhsm(3,k) = lhs(3,k)
   116: |||     F                 lhsm(4,k) = lhs(4,k) -
   117: |||            >                            dttz2 * speed(i,j,k+1)
   118: |||                       lhsm(5,k) = lhs(5,k)
   119: ||V----                end do
   120: ||        
   121: ||        
   122: ||        c---------------------------------------------------------------------
   123: ||        c                          FORWARD ELIMINATION
   124: ||        c---------------------------------------------------------------------
   125: ||        
   126: ||+---->               do    k = 0, grid_points(3)-3
   127: |||                       k1 = k  + 1
   128: |||                       k2 = k  + 2
   129: |||                       fac1      = 1.d0/lhs(3,k)
   130: |||                       lhs(4,k)  = fac1*lhs(4,k)
   131: |||                       lhs(5,k)  = fac1*lhs(5,k)
   132: |||*--->                  do    m = 1, 3
   133: ||||                         rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   134: |||*---                   end do
   135: |||                       lhs(3,k1) = lhs(3,k1) -
   136: |||            >                         lhs(2,k1)*lhs(4,k)
   137: |||                       lhs(4,k1) = lhs(4,k1) -
   138: |||            >                         lhs(2,k1)*lhs(5,k)
   139: |||*--->                  do    m = 1, 3
   140: ||||                         rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   141: ||||           >                         lhs(2,k1)*rhs(m,i,j,k)
   142: |||*---                   end do
   143: |||                       lhs(2,k2) = lhs(2,k2) -
   144: |||            >                         lhs(1,k2)*lhs(4,k)
   145: |||                       lhs(3,k2) = lhs(3,k2) -
   146: |||            >                         lhs(1,k2)*lhs(5,k)
   147: |||*--->                  do    m = 1, 3
   148: ||||                         rhs(m,i,j,k2) = rhs(m,i,j,k2) -
   149: ||||           >                         lhs(1,k2)*rhs(m,i,j,k)
   150: |||*---                   end do
   151: ||+----                end do
   152: ||        
   153: ||        c---------------------------------------------------------------------
   154: ||        c      The last two rows in this grid block are a bit different,
   155: ||        c      since they do not have two more rows available for the
   156: ||        c      elimination of off-diagonal entries
   157: ||        c---------------------------------------------------------------------
   158: ||                     k  = grid_points(3)-2
   159: ||                     k1 = grid_points(3)-1
   160: ||                     fac1      = 1.d0/lhs(3,k)
   161: ||                     lhs(4,k)  = fac1*lhs(4,k)
   162: ||                     lhs(5,k)  = fac1*lhs(5,k)
   163: ||*---->               do    m = 1, 3
   164: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   165: ||*----                end do
   166: ||                     lhs(3,k1) = lhs(3,k1) -
   167: ||             >                      lhs(2,k1)*lhs(4,k)
   168: ||                     lhs(4,k1) = lhs(4,k1) -
   169: ||             >                      lhs(2,k1)*lhs(5,k)
   170: ||*---->               do    m = 1, 3
   171: |||                       rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   172: |||            >                      lhs(2,k1)*rhs(m,i,j,k)
   173: ||*----                end do
   174: ||        c---------------------------------------------------------------------
   175: ||        c               scale the last row immediately
   176: ||        c---------------------------------------------------------------------
   177: ||                     fac2      = 1.d0/lhs(3,k1)
   178: ||*---->               do    m = 1, 3
   179: |||                       rhs(m,i,j,k1) = fac2*rhs(m,i,j,k1)
   180: ||*----                end do
   181: ||        
   182: ||        c---------------------------------------------------------------------
   183: ||        c      do the u+c and the u-c factors
   184: ||        c---------------------------------------------------------------------
   185: ||+---->               do    k = 0, grid_points(3)-3
   186: |||                       k1 = k  + 1
   187: |||                       k2 = k  + 2
   188: |||                       m = 4
   189: |||                       fac1       = 1.d0/lhsp(3,k)
   190: |||                       lhsp(4,k)  = fac1*lhsp(4,k)
   191: |||                       lhsp(5,k)  = fac1*lhsp(5,k)
   192: |||                       rhs(m,i,j,k)  = fac1*rhs(m,i,j,k)
   193: |||                       lhsp(3,k1) = lhsp(3,k1) -
   194: |||            >                       lhsp(2,k1)*lhsp(4,k)
   195: |||                       lhsp(4,k1) = lhsp(4,k1) -
   196: |||            >                       lhsp(2,k1)*lhsp(5,k)
   197: |||                       rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   198: |||            >                       lhsp(2,k1)*rhs(m,i,j,k)
   199: |||                       lhsp(2,k2) = lhsp(2,k2) -
   200: |||            >                       lhsp(1,k2)*lhsp(4,k)
   201: |||                       lhsp(3,k2) = lhsp(3,k2) -
   202: |||            >                       lhsp(1,k2)*lhsp(5,k)
   203: |||                       rhs(m,i,j,k2) = rhs(m,i,j,k2) -
   204: |||            >                       lhsp(1,k2)*rhs(m,i,j,k)
   205: |||                       m = 5
   206: |||                       fac1       = 1.d0/lhsm(3,k)
   207: |||                       lhsm(4,k)  = fac1*lhsm(4,k)
   208: |||                       lhsm(5,k)  = fac1*lhsm(5,k)
   209: |||                       rhs(m,i,j,k)  = fac1*rhs(m,i,j,k)
   210: |||                       lhsm(3,k1) = lhsm(3,k1) -
   211: |||            >                       lhsm(2,k1)*lhsm(4,k)
   212: |||                       lhsm(4,k1) = lhsm(4,k1) -
   213: |||            >                       lhsm(2,k1)*lhsm(5,k)
   214: |||                       rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   215: |||            >                       lhsm(2,k1)*rhs(m,i,j,k)
   216: |||                       lhsm(2,k2) = lhsm(2,k2) -
   217: |||            >                       lhsm(1,k2)*lhsm(4,k)
   218: |||                       lhsm(3,k2) = lhsm(3,k2) -
   219: |||            >                       lhsm(1,k2)*lhsm(5,k)
   220: |||                       rhs(m,i,j,k2) = rhs(m,i,j,k2) -
   221: |||            >                       lhsm(1,k2)*rhs(m,i,j,k)
   222: ||+----                end do
   223: ||        
   224: ||        c---------------------------------------------------------------------
   225: ||        c         And again the last two rows separately
   226: ||        c---------------------------------------------------------------------
   227: ||                     k  = grid_points(3)-2
   228: ||                     k1 = grid_points(3)-1
   229: ||                     m = 4
   230: ||                     fac1       = 1.d0/lhsp(3,k)
   231: ||                     lhsp(4,k)  = fac1*lhsp(4,k)
   232: ||                     lhsp(5,k)  = fac1*lhsp(5,k)
   233: ||                     rhs(m,i,j,k)  = fac1*rhs(m,i,j,k)
   234: ||                     lhsp(3,k1) = lhsp(3,k1) -
   235: ||             >                    lhsp(2,k1)*lhsp(4,k)
   236: ||                     lhsp(4,k1) = lhsp(4,k1) -
   237: ||             >                    lhsp(2,k1)*lhsp(5,k)
   238: ||                     rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   239: ||             >                    lhsp(2,k1)*rhs(m,i,j,k)
   240: ||                     m = 5
   241: ||                     fac1       = 1.d0/lhsm(3,k)
   242: ||                     lhsm(4,k)  = fac1*lhsm(4,k)
   243: ||                     lhsm(5,k)  = fac1*lhsm(5,k)
   244: ||                     rhs(m,i,j,k)  = fac1*rhs(m,i,j,k)
   245: ||                     lhsm(3,k1) = lhsm(3,k1) -
   246: ||             >                    lhsm(2,k1)*lhsm(4,k)
   247: ||                     lhsm(4,k1) = lhsm(4,k1) -
   248: ||             >                    lhsm(2,k1)*lhsm(5,k)
   249: ||                     rhs(m,i,j,k1) = rhs(m,i,j,k1) -
   250: ||             >                    lhsm(2,k1)*rhs(m,i,j,k)
   251: ||        c---------------------------------------------------------------------
   252: ||        c               Scale the last row immediately (some of this is overkill
   253: ||        c               if this is the last cell)
   254: ||        c---------------------------------------------------------------------
   255: ||                     rhs(4,i,j,k1) = rhs(4,i,j,k1)/lhsp(3,k1)
   256: ||                     rhs(5,i,j,k1) = rhs(5,i,j,k1)/lhsm(3,k1)
   257: ||        
   258: ||        
   259: ||        c---------------------------------------------------------------------
   260: ||        c                         BACKSUBSTITUTION
   261: ||        c---------------------------------------------------------------------
   262: ||        
   263: ||                     k  = grid_points(3)-2
   264: ||                     k1 = grid_points(3)-1
   265: ||*---->               do   m = 1, 3
   266: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) -
   267: |||            >                             lhs(4,k)*rhs(m,i,j,k1)
   268: ||*----                end do
   269: ||        
   270: ||                     rhs(4,i,j,k) = rhs(4,i,j,k) -
   271: ||             >                             lhsp(4,k)*rhs(4,i,j,k1)
   272: ||                     rhs(5,i,j,k) = rhs(5,i,j,k) -
   273: ||             >                             lhsm(4,k)*rhs(5,i,j,k1)
   274: ||        
   275: ||        c---------------------------------------------------------------------
   276: ||        c      Whether or not this is the last processor, we always have
   277: ||        c      to complete the back-substitution
   278: ||        c---------------------------------------------------------------------
   279: ||        
   280: ||        c---------------------------------------------------------------------
   281: ||        c      The first three factors
   282: ||        c---------------------------------------------------------------------
   283: ||+---->               do   k = grid_points(3)-3, 0, -1
   284: |||                       k1 = k  + 1
   285: |||                       k2 = k  + 2
   286: |||*--->                  do   m = 1, 3
   287: ||||                         rhs(m,i,j,k) = rhs(m,i,j,k) -
   288: ||||           >                          lhs(4,k)*rhs(m,i,j,k1) -
   289: ||||           >                          lhs(5,k)*rhs(m,i,j,k2)
   290: |||*---                   end do
   291: |||       
   292: |||       c---------------------------------------------------------------------
   293: |||       c      And the remaining two
   294: |||       c---------------------------------------------------------------------
   295: |||                       rhs(4,i,j,k) = rhs(4,i,j,k) -
   296: |||            >                          lhsp(4,k)*rhs(4,i,j,k1) -
   297: |||            >                          lhsp(5,k)*rhs(4,i,j,k2)
   298: |||                       rhs(5,i,j,k) = rhs(5,i,j,k) -
   299: |||            >                          lhsm(4,k)*rhs(5,i,j,k1) -
   300: |||            >                          lhsm(5,k)*rhs(5,i,j,k2)
   301: ||+----                end do
   302: ||        
   303: |+-----             end do
   304: P------          end do
   305:                  if (timeron) call timer_stop(t_zsolve)
   306:           
   307:                  call tzetar
   308:           
   309:                  return
   310:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
INLINE LIST

  ROOT: Z_SOLVE (z_solve.f:5)
  -> NOINLINE: TIMER_START (z_solve.f:33)
     *** Source for routine not found.
  -> NOINLINE: LHSINIT (z_solve.f:39)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (z_solve.f:305)
     *** Source for routine not found.
  -> NOINLINE: TZETAR (z_solve.f:307)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (z_solve.f:36)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : LHSINIT (z_solve.f:36)

    LOOP BEGIN: (z_solve.f:49)
      <Vectorized loop.>
      **  Runtime-tested short loop. (z_solve.f:49)
      *** The number of VGT,   VSC.    :  0,  0. (z_solve.f:49)
      *** The number of VLOAD, VSTORE. :  2,  2. (z_solve.f:49)
      *** Iteration count is assumed. Iteration count. : 162 (z_solve.f:49)
    LOOP END

    LOOP BEGIN: (z_solve.f:58)
      <Vectorized loop.>
      **  Runtime-tested short loop. (z_solve.f:58)
      *** The number of VGT,   VSC.    :  0,  0. (z_solve.f:58)
      *** The number of VLOAD, VSTORE. :  5,  5. (z_solve.f:58)
      *** Iteration count is assumed. Iteration count. : 162 (z_solve.f:58)
    LOOP END

    LOOP BEGIN: (z_solve.f:81)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (z_solve.f:81)
      *** The number of VLOAD, VSTORE. :  5,  5. (z_solve.f:81)
      *** Iteration count is assumed. Iteration count. : 163 (z_solve.f:81)
    LOOP END

    LOOP BEGIN: (z_solve.f:104)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (z_solve.f:104)
      *** The number of VLOAD, VSTORE. :  7, 10. (z_solve.f:104)
      *** Iteration count is assumed. Iteration count. : 163 (z_solve.f:104)
    LOOP END

    LOOP BEGIN: (z_solve.f:126)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (z_solve.f:126)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:148)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:140)
      *** Unvectorizable dependency. : LHS (z_solve.f:145)
      *** Unvectorizable dependency. : LHS (z_solve.f:143)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:133)
      *** Unvectorizable dependency. : LHS (z_solve.f:135)
    LOOP END

    LOOP BEGIN: (z_solve.f:185)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (z_solve.f:185)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:220)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:214)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:209)
      *** Unvectorizable dependency. : LHSP (z_solve.f:201)
      *** Unvectorizable dependency. : LHSM (z_solve.f:218)
      *** Unvectorizable dependency. : LHSP (z_solve.f:199)
      *** Unvectorizable dependency. : LHSM (z_solve.f:216)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:197)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:192)
      *** Unvectorizable dependency. : LHSP (z_solve.f:193)
      *** Unvectorizable dependency. : LHSM (z_solve.f:210)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:203)
    LOOP END

    LOOP BEGIN: (z_solve.f:283)
      <Unvectorized loop.>
      *** Overhead of loop division is too large. (z_solve.f:283)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:298)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:295)
      *** Dependency unknown. Unvectorizable dependency is assumed. : RHS (z_solve.f:287)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:57 2021
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 58 [s0-s4 s8-s12 s15-s16 s18-s63]
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

  LOOP BEGIN: (z_solve.f:36)
    *** Estimated execution cycle                       : 494
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 42
            Over basic blocks                           : 42
    *** The number of SCALAR REGISTER TRANSFER          : 19

    LOOP BEGIN: (z_solve.f:126)
      *** Estimated execution cycle                     : 146
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 2
              Over basic blocks                         : 2
    LOOP END

    LOOP BEGIN: (z_solve.f:185)
      *** Estimated execution cycle                     : 189
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
    LOOP END

    LOOP BEGIN: (z_solve.f:283)
      *** Estimated execution cycle                     : 90
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END


