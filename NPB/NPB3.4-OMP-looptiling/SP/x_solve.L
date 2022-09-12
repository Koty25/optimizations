NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: x_solve.f

PROCEDURE NAME: X_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    29: inl(1212): Source for routine not found.: TIMER_START
    30: par(1801): Parallel routine generated.: X_SOLVE$1
    32: vec( 110): Vectorization obstructive procedure reference.: LHSINIT
    32: par(1803): Parallelized by "do".
    32: opt(2495): Loop nest collapsed by "collapse" clause.
    35: inl(1212): Source for routine not found.: LHSINIT
    44: vec( 101): Vectorized loop.
    47: vec( 128): Fused multiply-add operation applied.
    53: vec( 124): Iteration count is assumed. Iteration count.: 102
    53: vec( 101): Vectorized loop.
    55: vec( 128): Fused multiply-add operation applied.
    56: vec( 128): Fused multiply-add operation applied.
    57: vec( 128): Fused multiply-add operation applied.
    75: vec( 124): Iteration count is assumed. Iteration count.: 103
    75: vec( 101): Vectorized loop.
    97: vec( 124): Iteration count is assumed. Iteration count.: 103
    97: vec( 101): Vectorized loop.
    99: vec( 128): Fused multiply-add operation applied.
   102: vec( 128): Fused multiply-add operation applied.
   106: vec( 128): Fused multiply-add operation applied.
   109: vec( 128): Fused multiply-add operation applied.
   122: vec( 103): Unvectorized loop.
   122: vec( 113): Overhead of loop division is too large.
   125: opt(1037): Feedback of array elements.: LHS
   126: opt(1037): Feedback of array elements.: LHS
   128: opt(1418): Constant-length loop is expanded.
   129: opt(1037): Feedback of array elements.: RHS
   129: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   131: opt(1037): Feedback of array elements.: LHS
   131: vec( 120): Unvectorizable dependency.: LHS
   135: opt(1418): Constant-length loop is expanded.
   136: opt(1037): Feedback of array elements.: RHS
   136: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   139: vec( 120): Unvectorizable dependency.: LHS
   141: vec( 120): Unvectorizable dependency.: LHS
   143: opt(1418): Constant-length loop is expanded.
   144: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   160: opt(1418): Constant-length loop is expanded.
   167: opt(1418): Constant-length loop is expanded.
   175: opt(1418): Constant-length loop is expanded.
   183: vec( 103): Unvectorized loop.
   183: vec( 113): Overhead of loop division is too large.
   187: opt(1037): Feedback of array elements.: LHSP
   188: opt(1037): Feedback of array elements.: LHSP
   190: opt(1037): Feedback of array elements.: RHS
   190: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   191: opt(1037): Feedback of array elements.: LHSP
   191: vec( 120): Unvectorizable dependency.: LHSP
   195: opt(1037): Feedback of array elements.: RHS
   195: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   197: vec( 120): Unvectorizable dependency.: LHSP
   199: vec( 120): Unvectorizable dependency.: LHSP
   201: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   204: opt(1037): Feedback of array elements.: LHSM
   205: opt(1037): Feedback of array elements.: LHSM
   207: opt(1037): Feedback of array elements.: RHS
   207: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   208: opt(1037): Feedback of array elements.: LHSM
   208: vec( 120): Unvectorizable dependency.: LHSM
   212: opt(1037): Feedback of array elements.: RHS
   212: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   214: vec( 120): Unvectorizable dependency.: LHSM
   216: vec( 120): Unvectorizable dependency.: LHSM
   218: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   263: opt(1418): Constant-length loop is expanded.
   276: vec( 103): Unvectorized loop.
   276: vec( 113): Overhead of loop division is too large.
   279: opt(1418): Constant-length loop is expanded.
   280: opt(1037): Feedback of array elements.: RHS
   280: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   288: opt(1037): Feedback of array elements.: RHS
   288: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   291: opt(1037): Feedback of array elements.: RHS
   291: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS
   298: inl(1212): Source for routine not found.: TIMER_STOP
   303: inl(1212): Source for routine not found.: NINVR


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: x_solve.f

PROCEDURE NAME: X_SOLVE
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine x_solve
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function performs the solution of the approximate factorization
    12:           c step in the x-direction for all five matrix components
    13:           c simultaneously. The Thomas algorithm is employed to solve the
    14:           c systems for the x-lines. Boundary conditions are non-periodic
    15:           c---------------------------------------------------------------------
    16:           
    17:                  use sp_data
    18:                  use work_lhs
    19:           
    20:                  implicit none
    21:           
    22:                  integer i, j, k, i1, i2, m
    23:                  double precision  ru1, fac1, fac2
    24:           
    25:           
    26:           c---------------------------------------------------------------------
    27:           c---------------------------------------------------------------------
    28:           
    29:                  if (timeron) call timer_start(t_xsolve)
    30:           !$omp parallel do default(shared) private(i,j,k,i1,i2,m,
    31:           !$omp&    ru1,fac1,fac2) collapse(2)
    32: P------>         do  k = 1, nz2
    33: |+----->            do  j = 1, ny2
    34: ||        
    35: ||                    call lhsinit(nx2+1, lhs, lhsp, lhsm)
    36: ||        
    37: ||        c---------------------------------------------------------------------
    38: ||        c Computes the left hand side for the three x-factors
    39: ||        c---------------------------------------------------------------------
    40: ||        
    41: ||        c---------------------------------------------------------------------
    42: ||        c      first fill the lhs for the u-eigenvalue
    43: ||        c---------------------------------------------------------------------
    44: ||V---->               do  i = 0, grid_points(1)-1
    45: |||                       ru1 = c3c4*rho_i(i,j,k)
    46: |||                       cv(i) = us(i,j,k)
    47: |||     F                 rhov(i) = dmax1(dx2+con43*ru1,
    48: |||            >                          dx5+c1c5*ru1,
    49: |||            >                          dxmax+ru1,
    50: |||            >                          dx1)
    51: ||V----                end do
    52: ||        
    53: ||V---->               do  i = 1, nx2
    54: |||                       lhs(1,i) =  0.0d0
    55: |||     F                 lhs(2,i) = -dttx2 * cv(i-1) - dttx1 * rhov(i-1)
    56: |||     F                 lhs(3,i) =  1.0d0 + c2dttx1 * rhov(i)
    57: |||     F                 lhs(4,i) =  dttx2 * cv(i+1) - dttx1 * rhov(i+1)
    58: |||                       lhs(5,i) =  0.0d0
    59: ||V----                end do
    60: ||        
    61: ||        c---------------------------------------------------------------------
    62: ||        c      add fourth order dissipation
    63: ||        c---------------------------------------------------------------------
    64: ||        
    65: ||                     i = 1
    66: ||                     lhs(3,i) = lhs(3,i) + comz5
    67: ||                     lhs(4,i) = lhs(4,i) - comz4
    68: ||                     lhs(5,i) = lhs(5,i) + comz1
    69: ||         
    70: ||                     lhs(2,i+1) = lhs(2,i+1) - comz4
    71: ||                     lhs(3,i+1) = lhs(3,i+1) + comz6
    72: ||                     lhs(4,i+1) = lhs(4,i+1) - comz4
    73: ||                     lhs(5,i+1) = lhs(5,i+1) + comz1
    74: ||        
    75: ||V---->               do   i=3, grid_points(1)-4
    76: |||                       lhs(1,i) = lhs(1,i) + comz1
    77: |||                       lhs(2,i) = lhs(2,i) - comz4
    78: |||                       lhs(3,i) = lhs(3,i) + comz6
    79: |||                       lhs(4,i) = lhs(4,i) - comz4
    80: |||                       lhs(5,i) = lhs(5,i) + comz1
    81: ||V----                end do
    82: ||        
    83: ||                     i = grid_points(1)-3
    84: ||                     lhs(1,i) = lhs(1,i) + comz1
    85: ||                     lhs(2,i) = lhs(2,i) - comz4
    86: ||                     lhs(3,i) = lhs(3,i) + comz6
    87: ||                     lhs(4,i) = lhs(4,i) - comz4
    88: ||        
    89: ||                     lhs(1,i+1) = lhs(1,i+1) + comz1
    90: ||                     lhs(2,i+1) = lhs(2,i+1) - comz4
    91: ||                     lhs(3,i+1) = lhs(3,i+1) + comz5
    92: ||        
    93: ||        c---------------------------------------------------------------------
    94: ||        c      subsequently, fill the other factors (u+c), (u-c) by adding to
    95: ||        c      the first
    96: ||        c---------------------------------------------------------------------
    97: ||V---->               do   i = 1, nx2
    98: |||                       lhsp(1,i) = lhs(1,i)
    99: |||     F                 lhsp(2,i) = lhs(2,i) -
   100: |||            >                            dttx2 * speed(i-1,j,k)
   101: |||                       lhsp(3,i) = lhs(3,i)
   102: |||     F                 lhsp(4,i) = lhs(4,i) +
   103: |||            >                            dttx2 * speed(i+1,j,k)
   104: |||                       lhsp(5,i) = lhs(5,i)
   105: |||                       lhsm(1,i) = lhs(1,i)
   106: |||     F                 lhsm(2,i) = lhs(2,i) +
   107: |||            >                            dttx2 * speed(i-1,j,k)
   108: |||                       lhsm(3,i) = lhs(3,i)
   109: |||     F                 lhsm(4,i) = lhs(4,i) -
   110: |||            >                            dttx2 * speed(i+1,j,k)
   111: |||                       lhsm(5,i) = lhs(5,i)
   112: ||V----                end do
   113: ||        
   114: ||        c---------------------------------------------------------------------
   115: ||        c                          FORWARD ELIMINATION
   116: ||        c---------------------------------------------------------------------
   117: ||        
   118: ||        c---------------------------------------------------------------------
   119: ||        c      perform the Thomas algorithm; first, FORWARD ELIMINATION
   120: ||        c---------------------------------------------------------------------
   121: ||        
   122: ||+---->               do    i = 0, grid_points(1)-3
   123: |||                       i1 = i  + 1
   124: |||                       i2 = i  + 2
   125: |||                       fac1      = 1.d0/lhs(3,i)
   126: |||                       lhs(4,i)  = fac1*lhs(4,i)
   127: |||                       lhs(5,i)  = fac1*lhs(5,i)
   128: |||*--->                  do    m = 1, 3
   129: ||||                         rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   130: |||*---                   end do
   131: |||                       lhs(3,i1) = lhs(3,i1) -
   132: |||            >                         lhs(2,i1)*lhs(4,i)
   133: |||                       lhs(4,i1) = lhs(4,i1) -
   134: |||            >                         lhs(2,i1)*lhs(5,i)
   135: |||*--->                  do    m = 1, 3
   136: ||||                         rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   137: ||||           >                         lhs(2,i1)*rhs(m,i,j,k)
   138: |||*---                   end do
   139: |||                       lhs(2,i2) = lhs(2,i2) -
   140: |||            >                         lhs(1,i2)*lhs(4,i)
   141: |||                       lhs(3,i2) = lhs(3,i2) -
   142: |||            >                         lhs(1,i2)*lhs(5,i)
   143: |||*--->                  do    m = 1, 3
   144: ||||                         rhs(m,i2,j,k) = rhs(m,i2,j,k) -
   145: ||||           >                         lhs(1,i2)*rhs(m,i,j,k)
   146: |||*---                   end do
   147: ||+----                end do
   148: ||        
   149: ||        c---------------------------------------------------------------------
   150: ||        c      The last two rows in this grid block are a bit different,
   151: ||        c      since they do not have two more rows available for the
   152: ||        c      elimination of off-diagonal entries
   153: ||        c---------------------------------------------------------------------
   154: ||        
   155: ||                     i  = grid_points(1)-2
   156: ||                     i1 = grid_points(1)-1
   157: ||                     fac1      = 1.d0/lhs(3,i)
   158: ||                     lhs(4,i)  = fac1*lhs(4,i)
   159: ||                     lhs(5,i)  = fac1*lhs(5,i)
   160: ||*---->               do    m = 1, 3
   161: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   162: ||*----                end do
   163: ||                     lhs(3,i1) = lhs(3,i1) -
   164: ||             >                      lhs(2,i1)*lhs(4,i)
   165: ||                     lhs(4,i1) = lhs(4,i1) -
   166: ||             >                      lhs(2,i1)*lhs(5,i)
   167: ||*---->               do    m = 1, 3
   168: |||                       rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   169: |||            >                      lhs(2,i1)*rhs(m,i,j,k)
   170: ||*----                end do
   171: ||        c---------------------------------------------------------------------
   172: ||        c            scale the last row immediately
   173: ||        c---------------------------------------------------------------------
   174: ||                     fac2             = 1.d0/lhs(3,i1)
   175: ||*---->               do    m = 1, 3
   176: |||                       rhs(m,i1,j,k) = fac2*rhs(m,i1,j,k)
   177: ||*----                end do
   178: ||        
   179: ||        c---------------------------------------------------------------------
   180: ||        c      do the u+c and the u-c factors
   181: ||        c---------------------------------------------------------------------
   182: ||        
   183: ||+---->               do    i = 0, grid_points(1)-3
   184: |||                       i1 = i  + 1
   185: |||                       i2 = i  + 2
   186: |||                       m = 4
   187: |||                       fac1       = 1.d0/lhsp(3,i)
   188: |||                       lhsp(4,i)  = fac1*lhsp(4,i)
   189: |||                       lhsp(5,i)  = fac1*lhsp(5,i)
   190: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   191: |||                       lhsp(3,i1) = lhsp(3,i1) -
   192: |||            >                        lhsp(2,i1)*lhsp(4,i)
   193: |||                       lhsp(4,i1) = lhsp(4,i1) -
   194: |||            >                        lhsp(2,i1)*lhsp(5,i)
   195: |||                       rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   196: |||            >                        lhsp(2,i1)*rhs(m,i,j,k)
   197: |||                       lhsp(2,i2) = lhsp(2,i2) -
   198: |||            >                        lhsp(1,i2)*lhsp(4,i)
   199: |||                       lhsp(3,i2) = lhsp(3,i2) -
   200: |||            >                        lhsp(1,i2)*lhsp(5,i)
   201: |||                       rhs(m,i2,j,k) = rhs(m,i2,j,k) -
   202: |||            >                        lhsp(1,i2)*rhs(m,i,j,k)
   203: |||                       m = 5
   204: |||                       fac1       = 1.d0/lhsm(3,i)
   205: |||                       lhsm(4,i)  = fac1*lhsm(4,i)
   206: |||                       lhsm(5,i)  = fac1*lhsm(5,i)
   207: |||                       rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   208: |||                       lhsm(3,i1) = lhsm(3,i1) -
   209: |||            >                        lhsm(2,i1)*lhsm(4,i)
   210: |||                       lhsm(4,i1) = lhsm(4,i1) -
   211: |||            >                        lhsm(2,i1)*lhsm(5,i)
   212: |||                       rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   213: |||            >                        lhsm(2,i1)*rhs(m,i,j,k)
   214: |||                       lhsm(2,i2) = lhsm(2,i2) -
   215: |||            >                        lhsm(1,i2)*lhsm(4,i)
   216: |||                       lhsm(3,i2) = lhsm(3,i2) -
   217: |||            >                        lhsm(1,i2)*lhsm(5,i)
   218: |||                       rhs(m,i2,j,k) = rhs(m,i2,j,k) -
   219: |||            >                        lhsm(1,i2)*rhs(m,i,j,k)
   220: ||+----                end do
   221: ||        
   222: ||        c---------------------------------------------------------------------
   223: ||        c         And again the last two rows separately
   224: ||        c---------------------------------------------------------------------
   225: ||                     i  = grid_points(1)-2
   226: ||                     i1 = grid_points(1)-1
   227: ||                     m = 4
   228: ||                     fac1       = 1.d0/lhsp(3,i)
   229: ||                     lhsp(4,i)  = fac1*lhsp(4,i)
   230: ||                     lhsp(5,i)  = fac1*lhsp(5,i)
   231: ||                     rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   232: ||                     lhsp(3,i1) = lhsp(3,i1) -
   233: ||             >                      lhsp(2,i1)*lhsp(4,i)
   234: ||                     lhsp(4,i1) = lhsp(4,i1) -
   235: ||             >                      lhsp(2,i1)*lhsp(5,i)
   236: ||                     rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   237: ||             >                      lhsp(2,i1)*rhs(m,i,j,k)
   238: ||                     m = 5
   239: ||                     fac1       = 1.d0/lhsm(3,i)
   240: ||                     lhsm(4,i)  = fac1*lhsm(4,i)
   241: ||                     lhsm(5,i)  = fac1*lhsm(5,i)
   242: ||                     rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
   243: ||                     lhsm(3,i1) = lhsm(3,i1) -
   244: ||             >                      lhsm(2,i1)*lhsm(4,i)
   245: ||                     lhsm(4,i1) = lhsm(4,i1) -
   246: ||             >                      lhsm(2,i1)*lhsm(5,i)
   247: ||                     rhs(m,i1,j,k) = rhs(m,i1,j,k) -
   248: ||             >                      lhsm(2,i1)*rhs(m,i,j,k)
   249: ||        c---------------------------------------------------------------------
   250: ||        c               Scale the last row immediately
   251: ||        c---------------------------------------------------------------------
   252: ||                     rhs(4,i1,j,k) = rhs(4,i1,j,k)/lhsp(3,i1)
   253: ||                     rhs(5,i1,j,k) = rhs(5,i1,j,k)/lhsm(3,i1)
   254: ||        
   255: ||        
   256: ||        c---------------------------------------------------------------------
   257: ||        c                         BACKSUBSTITUTION
   258: ||        c---------------------------------------------------------------------
   259: ||        
   260: ||        
   261: ||                     i  = grid_points(1)-2
   262: ||                     i1 = grid_points(1)-1
   263: ||*---->               do   m = 1, 3
   264: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) -
   265: |||            >                             lhs(4,i)*rhs(m,i1,j,k)
   266: ||*----                end do
   267: ||        
   268: ||                     rhs(4,i,j,k) = rhs(4,i,j,k) -
   269: ||             >                          lhsp(4,i)*rhs(4,i1,j,k)
   270: ||                     rhs(5,i,j,k) = rhs(5,i,j,k) -
   271: ||             >                          lhsm(4,i)*rhs(5,i1,j,k)
   272: ||        
   273: ||        c---------------------------------------------------------------------
   274: ||        c      The first three factors
   275: ||        c---------------------------------------------------------------------
   276: ||+---->               do    i = grid_points(1)-3, 0, -1
   277: |||                       i1 = i  + 1
   278: |||                       i2 = i  + 2
   279: |||*--->                  do   m = 1, 3
   280: ||||                         rhs(m,i,j,k) = rhs(m,i,j,k) -
   281: ||||           >                          lhs(4,i)*rhs(m,i1,j,k) -
   282: ||||           >                          lhs(5,i)*rhs(m,i2,j,k)
   283: |||*---                   end do
   284: |||       
   285: |||       c---------------------------------------------------------------------
   286: |||       c      And the remaining two
   287: |||       c---------------------------------------------------------------------
   288: |||                       rhs(4,i,j,k) = rhs(4,i,j,k) -
   289: |||            >                          lhsp(4,i)*rhs(4,i1,j,k) -
   290: |||            >                          lhsp(5,i)*rhs(4,i2,j,k)
   291: |||                       rhs(5,i,j,k) = rhs(5,i,j,k) -
   292: |||            >                          lhsm(4,i)*rhs(5,i1,j,k) -
   293: |||            >                          lhsm(5,i)*rhs(5,i2,j,k)
   294: ||+----                end do
   295: |+-----             end do
   296: |         
   297: P------          end do
   298:                  if (timeron) call timer_stop(t_xsolve)
   299:           
   300:           c---------------------------------------------------------------------
   301:           c      Do the block-diagonal inversion
   302:           c---------------------------------------------------------------------
   303:                  call ninvr
   304:           
   305:                  return
   306:                  end


