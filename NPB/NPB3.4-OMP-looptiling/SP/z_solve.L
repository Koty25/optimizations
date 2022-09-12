NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    33: inl(1212): Source for routine not found.: TIMER_START
    34: par(1801): Parallel routine generated.: Z_SOLVE$1
    36: vec( 110): Vectorization obstructive procedure reference.: LHSINIT
    36: par(1803): Parallelized by "do".
    36: opt(2495): Loop nest collapsed by "collapse" clause.
    39: inl(1212): Source for routine not found.: LHSINIT
    49: vec( 101): Vectorized loop.
    52: vec( 128): Fused multiply-add operation applied.
    58: vec( 124): Iteration count is assumed. Iteration count.: 102
    58: vec( 101): Vectorized loop.
    60: vec( 128): Fused multiply-add operation applied.
    61: vec( 128): Fused multiply-add operation applied.
    62: vec( 128): Fused multiply-add operation applied.
    81: vec( 124): Iteration count is assumed. Iteration count.: 103
    81: vec( 101): Vectorized loop.
   104: vec( 124): Iteration count is assumed. Iteration count.: 103
   104: vec( 101): Vectorized loop.
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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
FORMAT LIST

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


