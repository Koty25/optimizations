NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: y_solve_blk.f

PROCEDURE NAME: Y_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    46: par(1801): Parallel routine generated.: Y_SOLVE$1
    52: par(1803): Parallelized by "do".
    52: opt(2495): Loop nest collapsed by "collapse" clause.
    57: vec( 101): Vectorized loop.
    59: opt(1034): Multiple store conflict.: RHSX
    60: opt(1034): Multiple store conflict.: RHSX
    61: opt(1034): Multiple store conflict.: RHSX
    62: opt(1034): Multiple store conflict.: RHSX
    71: opt(1135): Outer loop conditionally executes inner loop.
    71: vec( 110): Vectorization obstructive procedure reference.: BINVRHS
    71: vec( 110): Vectorization obstructive procedure reference.: MATMUL_SUB
    71: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
    71: vec( 110): Vectorization obstructive procedure reference.: BINVCRHS
    71: vec( 110): Vectorization obstructive procedure reference.: MATVEC_SUB
    71: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
    78: vec( 101): Vectorized loop.
    85: opt(1034): Multiple store conflict.: FJAC
    98: vec( 128): Fused multiply-add operation applied.
   113: vec( 128): Fused multiply-add operation applied.
   118: vec( 128): Fused multiply-add operation applied.
   150: vec( 128): Fused multiply-add operation applied.
   164: opt(1062): Too many data dependency problems.
   172: vec( 101): Vectorized loop.
   177: vec( 128): Fused multiply-add operation applied.
   180: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   184: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   189: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   194: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   198: vec( 128): Fused multiply-add operation applied.
   201: vec( 128): Fused multiply-add operation applied.
   203: vec( 128): Fused multiply-add operation applied.
   205: vec( 128): Fused multiply-add operation applied.
   208: vec( 128): Fused multiply-add operation applied.
   210: vec( 128): Fused multiply-add operation applied.
   213: vec( 128): Fused multiply-add operation applied.
   215: vec( 128): Fused multiply-add operation applied.
   217: vec( 128): Fused multiply-add operation applied.
   219: vec( 128): Fused multiply-add operation applied.
   222: vec( 128): Fused multiply-add operation applied.
   225: vec( 128): Fused multiply-add operation applied.
   227: vec( 128): Fused multiply-add operation applied.
   229: vec( 128): Fused multiply-add operation applied.
   231: vec( 128): Fused multiply-add operation applied.
   233: vec( 128): Fused multiply-add operation applied.
   237: vec( 128): Fused multiply-add operation applied.
   246: vec( 128): Fused multiply-add operation applied.
   255: vec( 128): Fused multiply-add operation applied.
   264: vec( 128): Fused multiply-add operation applied.
   273: vec( 128): Fused multiply-add operation applied.
   277: vec( 128): Fused multiply-add operation applied.
   280: vec( 128): Fused multiply-add operation applied.
   282: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   286: vec( 128): Fused multiply-add operation applied.
   289: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   294: vec( 128): Fused multiply-add operation applied.
   296: vec( 128): Fused multiply-add operation applied.
   298: vec( 128): Fused multiply-add operation applied.
   301: vec( 128): Fused multiply-add operation applied.
   303: vec( 128): Fused multiply-add operation applied.
   305: vec( 128): Fused multiply-add operation applied.
   308: vec( 128): Fused multiply-add operation applied.
   310: vec( 128): Fused multiply-add operation applied.
   313: vec( 128): Fused multiply-add operation applied.
   315: vec( 128): Fused multiply-add operation applied.
   317: vec( 128): Fused multiply-add operation applied.
   319: vec( 128): Fused multiply-add operation applied.
   322: vec( 128): Fused multiply-add operation applied.
   325: vec( 128): Fused multiply-add operation applied.
   327: vec( 128): Fused multiply-add operation applied.
   329: vec( 128): Fused multiply-add operation applied.
   331: vec( 128): Fused multiply-add operation applied.
   333: vec( 128): Fused multiply-add operation applied.
   430: vec( 101): Vectorized loop.
   432: opt(1418): Constant-length loop is expanded.
   433: opt(1036): Potential feedback - use directive or compiler option if OK.: RHSX
   433: opt(1037): Feedback of array elements.: RHSX
   433: vec( 128): Fused multiply-add operation applied.
   441: opt(1135): Outer loop conditionally executes inner loop.
   444: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: y_solve_blk.f

PROCEDURE NAME: Y_SOLVE
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine y_solve
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     Performs line solves in Y direction by first factoring
    11:           c     the block-tridiagonal matrix into an upper triangular matrix,
    12:           c     and then performing back substitution to solve for the unknow
    13:           c     vectors of each line.
    14:           c
    15:           c     Make sure we treat elements zero to cell_size in the direction
    16:           c     of the sweep.
    17:           c---------------------------------------------------------------------
    18:           
    19:                 use bt_data
    20:                 use work_lhs
    21:           
    22:                 implicit none
    23:           
    24:                 integer i, j, k, m, n, jsize
    25:                 integer ii,ib,jj,jm,jp,jb
    26:                 double precision tmp1, tmp2, tmp3
    27:           
    28:           c---------------------------------------------------------------------
    29:           c---------------------------------------------------------------------
    30:           
    31:                 if (timeron) call timer_start(t_ysolve)
    32:           
    33:           c---------------------------------------------------------------------
    34:           c---------------------------------------------------------------------
    35:           
    36:           c---------------------------------------------------------------------
    37:           c     This function computes the left hand side for the three y-factors
    38:           c---------------------------------------------------------------------
    39:           
    40:                 jsize = grid_points(2)-1
    41:           
    42:           c---------------------------------------------------------------------
    43:           c     Compute the indices for storing the tri-diagonal matrix;
    44:           c     determine a (labeled f) and n jacobians for cell c
    45:           c---------------------------------------------------------------------
    46:           !$omp parallel default(shared) shared(jsize)
    47:           !$omp& private(i,j,k,m,n,ii,ib,jj,jm,jp,jb,tmp1,tmp2,tmp3)
    48:           
    49:                 call lhsinit(jsize)
    50:           
    51:           !$omp do collapse(2)
    52: P------>        do k = 1, grid_points(3)-2
    53: |+----->           do ii = 1, grid_points(1)-2, bsize
    54: ||        
    55: ||                    if (timeron) call timer_start(t_rdis1)
    56: ||+---->              do j=0,jsize
    57: |||V--->              do ib = 1, bsize
    58: ||||                     i = min(ii+ib-1, grid_points(1)-2)
    59: ||||    G                rhsx(ib,1,j) = rhs(1,i,j,k)
    60: ||||    G                rhsx(ib,2,j) = rhs(2,i,j,k)
    61: ||||    G                rhsx(ib,3,j) = rhs(3,i,j,k)
    62: ||||    G                rhsx(ib,4,j) = rhs(4,i,j,k)
    63: ||||    G                rhsx(ib,5,j) = rhs(5,i,j,k)
    64: |||V---               end do
    65: ||+----               end do
    66: ||                    if (timeron) call timer_stop(t_rdis1)
    67: ||        
    68: ||                    call lhsinit(0)
    69: ||        
    70: ||                    jb = 0
    71: ||+---->              do jj = 1, jsize-1
    72: |||                   jb = mod(jb + 1, 3)
    73: |||                   jm = min(2*jj - 3, 1)     ! -1 or 1
    74: |||                   jp = mod(jb + jm, 3) - 1
    75: |||       
    76: |||+--->              do j = jj+jm, jj+1
    77: ||||                  jp = jp + 1
    78: ||||V-->              do ib = 1, bsize
    79: |||||                    i = min(ii+ib-1, grid_points(1)-2)
    80: |||||     
    81: |||||   G                tmp1 = rho_i(i,j,k)
    82: |||||                    tmp2 = tmp1 * tmp1
    83: |||||                    tmp3 = tmp1 * tmp2
    84: |||||     
    85: |||||                    fjac(ib,1,1,jp) = 0.0d+00
    86: |||||                    fjac(ib,1,2,jp) = 0.0d+00
    87: |||||                    fjac(ib,1,3,jp) = 1.0d+00
    88: |||||                    fjac(ib,1,4,jp) = 0.0d+00
    89: |||||                    fjac(ib,1,5,jp) = 0.0d+00
    90: |||||     
    91: |||||   G                fjac(ib,2,1,jp) = - ( u(2,i,j,k)*u(3,i,j,k) )
    92: |||||          >              * tmp2
    93: |||||                    fjac(ib,2,2,jp) = u(3,i,j,k) * tmp1
    94: |||||                    fjac(ib,2,3,jp) = u(2,i,j,k) * tmp1
    95: |||||                    fjac(ib,2,4,jp) = 0.0d+00
    96: |||||                    fjac(ib,2,5,jp) = 0.0d+00
    97: |||||     
    98: |||||   G                fjac(ib,3,1,jp) = - ( u(3,i,j,k)*u(3,i,j,k)*tmp2)
    99: |||||          >              + c2 * qs(i,j,k)
   100: |||||                    fjac(ib,3,2,jp) = - c2 *  u(2,i,j,k) * tmp1
   101: |||||                    fjac(ib,3,3,jp) = ( 2.0d+00 - c2 )
   102: |||||          >              *  u(3,i,j,k) * tmp1
   103: |||||   G                fjac(ib,3,4,jp) = - c2 * u(4,i,j,k) * tmp1
   104: |||||                    fjac(ib,3,5,jp) = c2
   105: |||||     
   106: |||||                    fjac(ib,4,1,jp) = - ( u(3,i,j,k)*u(4,i,j,k) )
   107: |||||          >              * tmp2
   108: |||||                    fjac(ib,4,2,jp) = 0.0d+00
   109: |||||                    fjac(ib,4,3,jp) = u(4,i,j,k) * tmp1
   110: |||||                    fjac(ib,4,4,jp) = u(3,i,j,k) * tmp1
   111: |||||                    fjac(ib,4,5,jp) = 0.0d+00
   112: |||||     
   113: |||||   G                fjac(ib,5,1,jp) = ( c2 * 2.0d0 * square(i,j,k)
   114: |||||          >              - c1 * u(5,i,j,k) )
   115: |||||          >              * u(3,i,j,k) * tmp2
   116: |||||                    fjac(ib,5,2,jp) = - c2 * u(2,i,j,k)*u(3,i,j,k)
   117: |||||          >              * tmp2
   118: |||||   F                fjac(ib,5,3,jp) = c1 * u(5,i,j,k) * tmp1
   119: |||||          >              - c2
   120: |||||          >              * ( qs(i,j,k)
   121: |||||          >              + u(3,i,j,k)*u(3,i,j,k) * tmp2 )
   122: |||||                    fjac(ib,5,4,jp) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
   123: |||||          >              * tmp2
   124: |||||                    fjac(ib,5,5,jp) = c1 * u(3,i,j,k) * tmp1
   125: |||||     
   126: |||||                    njac(ib,1,1,jp) = 0.0d+00
   127: |||||                    njac(ib,1,2,jp) = 0.0d+00
   128: |||||                    njac(ib,1,3,jp) = 0.0d+00
   129: |||||                    njac(ib,1,4,jp) = 0.0d+00
   130: |||||                    njac(ib,1,5,jp) = 0.0d+00
   131: |||||     
   132: |||||                    njac(ib,2,1,jp) = - c3c4 * tmp2 * u(2,i,j,k)
   133: |||||                    njac(ib,2,2,jp) =   c3c4 * tmp1
   134: |||||                    njac(ib,2,3,jp) =   0.0d+00
   135: |||||                    njac(ib,2,4,jp) =   0.0d+00
   136: |||||                    njac(ib,2,5,jp) =   0.0d+00
   137: |||||     
   138: |||||                    njac(ib,3,1,jp) = - con43 * c3c4 * tmp2 * u(3,i,j,k)
   139: |||||                    njac(ib,3,2,jp) =   0.0d+00
   140: |||||                    njac(ib,3,3,jp) =   con43 * c3c4 * tmp1
   141: |||||                    njac(ib,3,4,jp) =   0.0d+00
   142: |||||                    njac(ib,3,5,jp) =   0.0d+00
   143: |||||     
   144: |||||                    njac(ib,4,1,jp) = - c3c4 * tmp2 * u(4,i,j,k)
   145: |||||                    njac(ib,4,2,jp) =   0.0d+00
   146: |||||                    njac(ib,4,3,jp) =   0.0d+00
   147: |||||                    njac(ib,4,4,jp) =   c3c4 * tmp1
   148: |||||                    njac(ib,4,5,jp) =   0.0d+00
   149: |||||     
   150: |||||   F                njac(ib,5,1,jp) = - (  c3c4
   151: |||||          >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   152: |||||          >              - ( con43 * c3c4
   153: |||||          >              - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   154: |||||          >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   155: |||||          >              - c1345 * tmp2 * u(5,i,j,k)
   156: |||||     
   157: |||||                    njac(ib,5,2,jp) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
   158: |||||                    njac(ib,5,3,jp) = ( con43 * c3c4
   159: |||||          >              - c1345 ) * tmp2 * u(3,i,j,k)
   160: |||||                    njac(ib,5,4,jp) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
   161: |||||                    njac(ib,5,5,jp) = ( c1345 ) * tmp1
   162: |||||     
   163: ||||V--               enddo
   164: |||+---               enddo
   165: |||       
   166: |||       c---------------------------------------------------------------------
   167: |||       c     now joacobians set, so form left hand side in y direction
   168: |||       c---------------------------------------------------------------------
   169: |||                   jm = mod(jb + 2, 3)
   170: |||                   j = jj
   171: |||       !dir$ vector always
   172: |||V--->              do ib = 1, bsize
   173: ||||      
   174: ||||                     tmp1 = dt * ty1
   175: ||||                     tmp2 = dt * ty2
   176: ||||      
   177: ||||    F                lhsa(ib,1,1,1) = - tmp2 * fjac(ib,1,1,jm)
   178: ||||           >              - tmp1 * njac(ib,1,1,jm)
   179: ||||           >              - tmp1 * dy1
   180: ||||    F                lhsa(ib,1,2,1) = - tmp2 * fjac(ib,1,2,jm)
   181: ||||           >              - tmp1 * njac(ib,1,2,jm)
   182: ||||    F                lhsa(ib,1,3,1) = - tmp2 * fjac(ib,1,3,jm)
   183: ||||           >              - tmp1 * njac(ib,1,3,jm)
   184: ||||    F                lhsa(ib,1,4,1) = - tmp2 * fjac(ib,1,4,jm)
   185: ||||           >              - tmp1 * njac(ib,1,4,jm)
   186: ||||    F                lhsa(ib,1,5,1) = - tmp2 * fjac(ib,1,5,jm)
   187: ||||           >              - tmp1 * njac(ib,1,5,jm)
   188: ||||      
   189: ||||    F                lhsa(ib,2,1,1) = - tmp2 * fjac(ib,2,1,jm)
   190: ||||           >              - tmp1 * njac(ib,2,1,jm)
   191: ||||    F                lhsa(ib,2,2,1) = - tmp2 * fjac(ib,2,2,jm)
   192: ||||           >              - tmp1 * njac(ib,2,2,jm)
   193: ||||           >              - tmp1 * dy2
   194: ||||    F                lhsa(ib,2,3,1) = - tmp2 * fjac(ib,2,3,jm)
   195: ||||           >              - tmp1 * njac(ib,2,3,jm)
   196: ||||    F                lhsa(ib,2,4,1) = - tmp2 * fjac(ib,2,4,jm)
   197: ||||           >              - tmp1 * njac(ib,2,4,jm)
   198: ||||    F                lhsa(ib,2,5,1) = - tmp2 * fjac(ib,2,5,jm)
   199: ||||           >              - tmp1 * njac(ib,2,5,jm)
   200: ||||      
   201: ||||    F                lhsa(ib,3,1,1) = - tmp2 * fjac(ib,3,1,jm)
   202: ||||           >              - tmp1 * njac(ib,3,1,jm)
   203: ||||    F                lhsa(ib,3,2,1) = - tmp2 * fjac(ib,3,2,jm)
   204: ||||           >              - tmp1 * njac(ib,3,2,jm)
   205: ||||    F                lhsa(ib,3,3,1) = - tmp2 * fjac(ib,3,3,jm)
   206: ||||           >              - tmp1 * njac(ib,3,3,jm)
   207: ||||           >              - tmp1 * dy3
   208: ||||    F                lhsa(ib,3,4,1) = - tmp2 * fjac(ib,3,4,jm)
   209: ||||           >              - tmp1 * njac(ib,3,4,jm)
   210: ||||    F                lhsa(ib,3,5,1) = - tmp2 * fjac(ib,3,5,jm)
   211: ||||           >              - tmp1 * njac(ib,3,5,jm)
   212: ||||      
   213: ||||    F                lhsa(ib,4,1,1) = - tmp2 * fjac(ib,4,1,jm)
   214: ||||           >              - tmp1 * njac(ib,4,1,jm)
   215: ||||    F                lhsa(ib,4,2,1) = - tmp2 * fjac(ib,4,2,jm)
   216: ||||           >              - tmp1 * njac(ib,4,2,jm)
   217: ||||    F                lhsa(ib,4,3,1) = - tmp2 * fjac(ib,4,3,jm)
   218: ||||           >              - tmp1 * njac(ib,4,3,jm)
   219: ||||    F                lhsa(ib,4,4,1) = - tmp2 * fjac(ib,4,4,jm)
   220: ||||           >              - tmp1 * njac(ib,4,4,jm)
   221: ||||           >              - tmp1 * dy4
   222: ||||    F                lhsa(ib,4,5,1) = - tmp2 * fjac(ib,4,5,jm)
   223: ||||           >              - tmp1 * njac(ib,4,5,jm)
   224: ||||      
   225: ||||    F                lhsa(ib,5,1,1) = - tmp2 * fjac(ib,5,1,jm)
   226: ||||           >              - tmp1 * njac(ib,5,1,jm)
   227: ||||    F                lhsa(ib,5,2,1) = - tmp2 * fjac(ib,5,2,jm)
   228: ||||           >              - tmp1 * njac(ib,5,2,jm)
   229: ||||    F                lhsa(ib,5,3,1) = - tmp2 * fjac(ib,5,3,jm)
   230: ||||           >              - tmp1 * njac(ib,5,3,jm)
   231: ||||    F                lhsa(ib,5,4,1) = - tmp2 * fjac(ib,5,4,jm)
   232: ||||           >              - tmp1 * njac(ib,5,4,jm)
   233: ||||    F                lhsa(ib,5,5,1) = - tmp2 * fjac(ib,5,5,jm)
   234: ||||           >              - tmp1 * njac(ib,5,5,jm)
   235: ||||           >              - tmp1 * dy5
   236: ||||      
   237: ||||    F                lhsb(ib,1,1,1) = 1.0d+00
   238: ||||           >              + tmp1 * 2.0d+00 * njac(ib,1,1,jb)
   239: ||||           >              + tmp1 * 2.0d+00 * dy1
   240: ||||                     lhsb(ib,1,2,1) = tmp1 * 2.0d+00 * njac(ib,1,2,jb)
   241: ||||                     lhsb(ib,1,3,1) = tmp1 * 2.0d+00 * njac(ib,1,3,jb)
   242: ||||                     lhsb(ib,1,4,1) = tmp1 * 2.0d+00 * njac(ib,1,4,jb)
   243: ||||                     lhsb(ib,1,5,1) = tmp1 * 2.0d+00 * njac(ib,1,5,jb)
   244: ||||      
   245: ||||                     lhsb(ib,2,1,1) = tmp1 * 2.0d+00 * njac(ib,2,1,jb)
   246: ||||    F                lhsb(ib,2,2,1) = 1.0d+00
   247: ||||           >              + tmp1 * 2.0d+00 * njac(ib,2,2,jb)
   248: ||||           >              + tmp1 * 2.0d+00 * dy2
   249: ||||                     lhsb(ib,2,3,1) = tmp1 * 2.0d+00 * njac(ib,2,3,jb)
   250: ||||                     lhsb(ib,2,4,1) = tmp1 * 2.0d+00 * njac(ib,2,4,jb)
   251: ||||                     lhsb(ib,2,5,1) = tmp1 * 2.0d+00 * njac(ib,2,5,jb)
   252: ||||      
   253: ||||                     lhsb(ib,3,1,1) = tmp1 * 2.0d+00 * njac(ib,3,1,jb)
   254: ||||                     lhsb(ib,3,2,1) = tmp1 * 2.0d+00 * njac(ib,3,2,jb)
   255: ||||    F                lhsb(ib,3,3,1) = 1.0d+00
   256: ||||           >              + tmp1 * 2.0d+00 * njac(ib,3,3,jb)
   257: ||||           >              + tmp1 * 2.0d+00 * dy3
   258: ||||                     lhsb(ib,3,4,1) = tmp1 * 2.0d+00 * njac(ib,3,4,jb)
   259: ||||                     lhsb(ib,3,5,1) = tmp1 * 2.0d+00 * njac(ib,3,5,jb)
   260: ||||      
   261: ||||                     lhsb(ib,4,1,1) = tmp1 * 2.0d+00 * njac(ib,4,1,jb)
   262: ||||                     lhsb(ib,4,2,1) = tmp1 * 2.0d+00 * njac(ib,4,2,jb)
   263: ||||                     lhsb(ib,4,3,1) = tmp1 * 2.0d+00 * njac(ib,4,3,jb)
   264: ||||    F                lhsb(ib,4,4,1) = 1.0d+00
   265: ||||           >              + tmp1 * 2.0d+00 * njac(ib,4,4,jb)
   266: ||||           >              + tmp1 * 2.0d+00 * dy4
   267: ||||                     lhsb(ib,4,5,1) = tmp1 * 2.0d+00 * njac(ib,4,5,jb)
   268: ||||      
   269: ||||                     lhsb(ib,5,1,1) = tmp1 * 2.0d+00 * njac(ib,5,1,jb)
   270: ||||                     lhsb(ib,5,2,1) = tmp1 * 2.0d+00 * njac(ib,5,2,jb)
   271: ||||                     lhsb(ib,5,3,1) = tmp1 * 2.0d+00 * njac(ib,5,3,jb)
   272: ||||                     lhsb(ib,5,4,1) = tmp1 * 2.0d+00 * njac(ib,5,4,jb)
   273: ||||    F                lhsb(ib,5,5,1) = 1.0d+00
   274: ||||           >              + tmp1 * 2.0d+00 * njac(ib,5,5,jb)
   275: ||||           >              + tmp1 * 2.0d+00 * dy5
   276: ||||      
   277: ||||    F                lhsc(ib,1,1,j) =  tmp2 * fjac(ib,1,1,jp)
   278: ||||           >              - tmp1 * njac(ib,1,1,jp)
   279: ||||           >              - tmp1 * dy1
   280: ||||    F                lhsc(ib,1,2,j) =  tmp2 * fjac(ib,1,2,jp)
   281: ||||           >              - tmp1 * njac(ib,1,2,jp)
   282: ||||    F                lhsc(ib,1,3,j) =  tmp2 * fjac(ib,1,3,jp)
   283: ||||           >              - tmp1 * njac(ib,1,3,jp)
   284: ||||    F                lhsc(ib,1,4,j) =  tmp2 * fjac(ib,1,4,jp)
   285: ||||           >              - tmp1 * njac(ib,1,4,jp)
   286: ||||    F                lhsc(ib,1,5,j) =  tmp2 * fjac(ib,1,5,jp)
   287: ||||           >              - tmp1 * njac(ib,1,5,jp)
   288: ||||      
   289: ||||    F                lhsc(ib,2,1,j) =  tmp2 * fjac(ib,2,1,jp)
   290: ||||           >              - tmp1 * njac(ib,2,1,jp)
   291: ||||    F                lhsc(ib,2,2,j) =  tmp2 * fjac(ib,2,2,jp)
   292: ||||           >              - tmp1 * njac(ib,2,2,jp)
   293: ||||           >              - tmp1 * dy2
   294: ||||    F                lhsc(ib,2,3,j) =  tmp2 * fjac(ib,2,3,jp)
   295: ||||           >              - tmp1 * njac(ib,2,3,jp)
   296: ||||    F                lhsc(ib,2,4,j) =  tmp2 * fjac(ib,2,4,jp)
   297: ||||           >              - tmp1 * njac(ib,2,4,jp)
   298: ||||    F                lhsc(ib,2,5,j) =  tmp2 * fjac(ib,2,5,jp)
   299: ||||           >              - tmp1 * njac(ib,2,5,jp)
   300: ||||      
   301: ||||    F                lhsc(ib,3,1,j) =  tmp2 * fjac(ib,3,1,jp)
   302: ||||           >              - tmp1 * njac(ib,3,1,jp)
   303: ||||    F                lhsc(ib,3,2,j) =  tmp2 * fjac(ib,3,2,jp)
   304: ||||           >              - tmp1 * njac(ib,3,2,jp)
   305: ||||    F                lhsc(ib,3,3,j) =  tmp2 * fjac(ib,3,3,jp)
   306: ||||           >              - tmp1 * njac(ib,3,3,jp)
   307: ||||           >              - tmp1 * dy3
   308: ||||    F                lhsc(ib,3,4,j) =  tmp2 * fjac(ib,3,4,jp)
   309: ||||           >              - tmp1 * njac(ib,3,4,jp)
   310: ||||    F                lhsc(ib,3,5,j) =  tmp2 * fjac(ib,3,5,jp)
   311: ||||           >              - tmp1 * njac(ib,3,5,jp)
   312: ||||      
   313: ||||    F                lhsc(ib,4,1,j) =  tmp2 * fjac(ib,4,1,jp)
   314: ||||           >              - tmp1 * njac(ib,4,1,jp)
   315: ||||    F                lhsc(ib,4,2,j) =  tmp2 * fjac(ib,4,2,jp)
   316: ||||           >              - tmp1 * njac(ib,4,2,jp)
   317: ||||    F                lhsc(ib,4,3,j) =  tmp2 * fjac(ib,4,3,jp)
   318: ||||           >              - tmp1 * njac(ib,4,3,jp)
   319: ||||    F                lhsc(ib,4,4,j) =  tmp2 * fjac(ib,4,4,jp)
   320: ||||           >              - tmp1 * njac(ib,4,4,jp)
   321: ||||           >              - tmp1 * dy4
   322: ||||    F                lhsc(ib,4,5,j) =  tmp2 * fjac(ib,4,5,jp)
   323: ||||           >              - tmp1 * njac(ib,4,5,jp)
   324: ||||      
   325: ||||    F                lhsc(ib,5,1,j) =  tmp2 * fjac(ib,5,1,jp)
   326: ||||           >              - tmp1 * njac(ib,5,1,jp)
   327: ||||    F                lhsc(ib,5,2,j) =  tmp2 * fjac(ib,5,2,jp)
   328: ||||           >              - tmp1 * njac(ib,5,2,jp)
   329: ||||    F                lhsc(ib,5,3,j) =  tmp2 * fjac(ib,5,3,jp)
   330: ||||           >              - tmp1 * njac(ib,5,3,jp)
   331: ||||    F                lhsc(ib,5,4,j) =  tmp2 * fjac(ib,5,4,jp)
   332: ||||           >              - tmp1 * njac(ib,5,4,jp)
   333: ||||    F                lhsc(ib,5,5,j) =  tmp2 * fjac(ib,5,5,jp)
   334: ||||           >              - tmp1 * njac(ib,5,5,jp)
   335: ||||           >              - tmp1 * dy5
   336: ||||      
   337: |||V---               enddo
   338: |||       
   339: |||       c---------------------------------------------------------------------
   340: |||       c---------------------------------------------------------------------
   341: |||       
   342: |||       c---------------------------------------------------------------------
   343: |||       c     performs guaussian elimination on this cell.
   344: |||       c
   345: |||       c     assumes that unpacking routines for non-first cells
   346: |||       c     preload C' and rhs' from previous cell.
   347: |||       c
   348: |||       c     assumed send happens outside this routine, but that
   349: |||       c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
   350: |||       c---------------------------------------------------------------------
   351: |||       
   352: |||                   if (timeron) call timer_start(t_solsub)
   353: |||       c---------------------------------------------------------------------
   354: |||       c     multiply c(i,0,k) by b_inverse and copy back to c
   355: |||       c     multiply rhs(0) by b_inverse(0) and copy to rhs
   356: |||       c---------------------------------------------------------------------
   357: |||                   if (jj .eq. 1) then
   358: |||                   call binvcrhs( lhsb(1,1,1,0),
   359: |||            >                        lhsc(1,1,1,0),
   360: |||            >                        rhsx(1,1,0) )
   361: |||                   endif
   362: |||       
   363: |||       c---------------------------------------------------------------------
   364: |||       c     begin inner most do loop
   365: |||       c     do all the elements of the cell unless last
   366: |||       c---------------------------------------------------------------------
   367: |||       
   368: |||       c---------------------------------------------------------------------
   369: |||       c     subtract A*lhs_vector(j-1) from lhs_vector(j)
   370: |||       c
   371: |||       c     rhs(j) = rhs(j) - A*rhs(j-1)
   372: |||       c---------------------------------------------------------------------
   373: |||                      call matvec_sub(lhsa(1,1,1,1),
   374: |||            >                         rhsx(1,1,j-1),rhsx(1,1,j))
   375: |||       
   376: |||       c---------------------------------------------------------------------
   377: |||       c     B(j) = B(j) - C(j-1)*A(j)
   378: |||       c---------------------------------------------------------------------
   379: |||                      call matmul_sub(lhsa(1,1,1,1),
   380: |||            >                         lhsc(1,1,1,j-1),
   381: |||            >                         lhsb(1,1,1,1))
   382: |||       
   383: |||       c---------------------------------------------------------------------
   384: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   385: |||       c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
   386: |||       c---------------------------------------------------------------------
   387: |||                      call binvcrhs( lhsb(1,1,1,1),
   388: |||            >                        lhsc(1,1,1,j),
   389: |||            >                        rhsx(1,1,j) )
   390: |||       
   391: |||       
   392: |||                   if (jj .eq. jsize-1) then
   393: |||       c---------------------------------------------------------------------
   394: |||       c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
   395: |||       c---------------------------------------------------------------------
   396: |||                   call matvec_sub(lhsa(1,1,1,2),
   397: |||            >                         rhsx(1,1,jsize-1),rhsx(1,1,jsize))
   398: |||       
   399: |||       c---------------------------------------------------------------------
   400: |||       c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
   401: |||       c     call matmul_sub(aa,i,jsize,k,c,
   402: |||       c     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
   403: |||       c---------------------------------------------------------------------
   404: |||                   call matmul_sub(lhsa(1,1,1,2),
   405: |||            >                         lhsc(1,1,1,jsize-1),
   406: |||            >                         lhsb(1,1,1,2))
   407: |||       
   408: |||       c---------------------------------------------------------------------
   409: |||       c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
   410: |||       c---------------------------------------------------------------------
   411: |||                   call binvrhs( lhsb(1,1,1,2),
   412: |||            >                       rhsx(1,1,jsize) )
   413: |||                   endif
   414: |||       
   415: |||                   if (timeron) call timer_stop(t_solsub)
   416: |||       
   417: ||+----               enddo
   418: ||        
   419: ||        c---------------------------------------------------------------------
   420: ||        c     back solve: if last cell, then generate U(jsize)=rhs(jsize)
   421: ||        c     else assume U(jsize) is loaded in un pack backsub_info
   422: ||        c     so just use it
   423: ||        c     after call u(jstart) will be sent to next cell
   424: ||        c---------------------------------------------------------------------
   425: ||         
   426: ||+---->              do j=jsize-1,0,-1
   427: |||       !dir$ vector always
   428: |||+--->              do ib=1,bsize
   429: ||||      !dir$ unroll
   430: ||||V-->                 do m=1,BLOCK_SIZE
   431: |||||     !dir$ unroll
   432: |||||*->                    do n=1,BLOCK_SIZE
   433: ||||||  F                      rhsx(ib,m,j) = rhsx(ib,m,j)
   434: ||||||         >                    - lhsc(ib,m,n,j)*rhsx(ib,n,j+1)
   435: |||||*-                     enddo
   436: ||||V--                  enddo
   437: |||+---               enddo
   438: ||+----               enddo
   439: ||        
   440: ||                    if (timeron) call timer_start(t_rdis1)
   441: ||+---->              do ib = 1, bsize
   442: |||                      i = ii+ib-1
   443: |||                      if (i .lt. grid_points(1)-1) then
   444: |||V--->                 do j=0,jsize
   445: ||||                        rhs(1,i,j,k) = rhsx(ib,1,j)
   446: ||||                        rhs(2,i,j,k) = rhsx(ib,2,j)
   447: ||||                        rhs(3,i,j,k) = rhsx(ib,3,j)
   448: ||||                        rhs(4,i,j,k) = rhsx(ib,4,j)
   449: ||||                        rhs(5,i,j,k) = rhsx(ib,5,j)
   450: |||V---                  end do
   451: |||                      endif
   452: ||+----               end do
   453: ||                    if (timeron) call timer_stop(t_rdis1)
   454: ||        
   455: |+-----            enddo
   456: P------         enddo
   457:           !$omp end parallel
   458:                 if (timeron) call timer_stop(t_ysolve)
   459:           
   460:                 return
   461:                 end


