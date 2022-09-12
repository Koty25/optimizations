NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:58 2020
FILE NAME: z_solve_blk.f

PROCEDURE NAME: Z_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    46: par(1801): Parallel routine generated.: Z_SOLVE$1
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
   105: vec( 128): Fused multiply-add operation applied.
   113: vec( 128): Fused multiply-add operation applied.
   120: vec( 128): Fused multiply-add operation applied.
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
   443: vec( 101): Vectorized loop.
   445: opt(1418): Constant-length loop is expanded.
   446: opt(1036): Potential feedback - use directive or compiler option if OK.: RHSX
   446: opt(1037): Feedback of array elements.: RHSX
   446: vec( 128): Fused multiply-add operation applied.
   454: opt(1135): Outer loop conditionally executes inner loop.
   457: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:58 2020
FILE NAME: z_solve_blk.f

PROCEDURE NAME: Z_SOLVE
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine z_solve
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     Performs line solves in Z direction by first factoring
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
    24:                 integer i, j, k, m, n, ksize
    25:                 integer ii, ib,kk,km,kp,kb
    26:                 double precision tmp1, tmp2, tmp3
    27:           
    28:           c---------------------------------------------------------------------
    29:           c---------------------------------------------------------------------
    30:           
    31:                 if (timeron) call timer_start(t_zsolve)
    32:           
    33:           c---------------------------------------------------------------------
    34:           c---------------------------------------------------------------------
    35:           
    36:           c---------------------------------------------------------------------
    37:           c     This function computes the left hand side for the three z-factors
    38:           c---------------------------------------------------------------------
    39:           
    40:                 ksize = grid_points(3)-1
    41:           
    42:           c---------------------------------------------------------------------
    43:           c     Compute the indices for storing the block-diagonal matrix;
    44:           c     determine c (labeled f) and s jacobians
    45:           c---------------------------------------------------------------------
    46:           !$omp parallel default(shared) shared(ksize)
    47:           !$omp& private(i,j,k,m,n,ii,ib,kk,km,kp,kb,tmp1,tmp2,tmp3)
    48:           
    49:                 call lhsinit(ksize)
    50:           
    51:           !$omp do collapse(2)
    52: P------>        do j = 1, grid_points(2)-2
    53: |+----->           do ii = 1, grid_points(1)-2, bsize
    54: ||        
    55: ||                    if (timeron) call timer_start(t_rdis1)
    56: ||+---->              do k=0,ksize
    57: |||V--->              do ib = 1, bsize
    58: ||||                     i = min(ii+ib-1, grid_points(1)-2)
    59: ||||    G                rhsx(ib,1,k) = rhs(1,i,j,k)
    60: ||||    G                rhsx(ib,2,k) = rhs(2,i,j,k)
    61: ||||    G                rhsx(ib,3,k) = rhs(3,i,j,k)
    62: ||||    G                rhsx(ib,4,k) = rhs(4,i,j,k)
    63: ||||    G                rhsx(ib,5,k) = rhs(5,i,j,k)
    64: |||V---               end do
    65: ||+----               end do
    66: ||                    if (timeron) call timer_stop(t_rdis1)
    67: ||        
    68: ||                    call lhsinit(0)
    69: ||        
    70: ||                    kb = 0
    71: ||+---->              do kk = 1, ksize-1
    72: |||                   kb = mod(kb + 1, 3)
    73: |||                   km = min(2*kk - 3, 1)     ! -1 or 1
    74: |||                   kp = mod(kb + km, 3) - 1
    75: |||       
    76: |||+--->              do k = kk+km, kk+1
    77: ||||                  kp = kp + 1
    78: ||||V-->              do ib = 1, bsize
    79: |||||                    i = min(ii+ib-1, grid_points(1)-2)
    80: |||||     
    81: |||||   G                tmp1 = 1.0d+00 / u(1,i,j,k)
    82: |||||                    tmp2 = tmp1 * tmp1
    83: |||||                    tmp3 = tmp1 * tmp2
    84: |||||     
    85: |||||                    fjac(ib,1,1,kp) = 0.0d+00
    86: |||||                    fjac(ib,1,2,kp) = 0.0d+00
    87: |||||                    fjac(ib,1,3,kp) = 0.0d+00
    88: |||||                    fjac(ib,1,4,kp) = 1.0d+00
    89: |||||                    fjac(ib,1,5,kp) = 0.0d+00
    90: |||||     
    91: |||||   G                fjac(ib,2,1,kp) = - ( u(2,i,j,k)*u(4,i,j,k) )
    92: |||||          >              * tmp2
    93: |||||                    fjac(ib,2,2,kp) = u(4,i,j,k) * tmp1
    94: |||||                    fjac(ib,2,3,kp) = 0.0d+00
    95: |||||                    fjac(ib,2,4,kp) = u(2,i,j,k) * tmp1
    96: |||||                    fjac(ib,2,5,kp) = 0.0d+00
    97: |||||     
    98: |||||   G                fjac(ib,3,1,kp) = - ( u(3,i,j,k)*u(4,i,j,k) )
    99: |||||          >              * tmp2
   100: |||||                    fjac(ib,3,2,kp) = 0.0d+00
   101: |||||                    fjac(ib,3,3,kp) = u(4,i,j,k) * tmp1
   102: |||||                    fjac(ib,3,4,kp) = u(3,i,j,k) * tmp1
   103: |||||                    fjac(ib,3,5,kp) = 0.0d+00
   104: |||||     
   105: |||||   G                fjac(ib,4,1,kp) = - (u(4,i,j,k)*u(4,i,j,k) * tmp2 )
   106: |||||          >              + c2 * qs(i,j,k)
   107: |||||                    fjac(ib,4,2,kp) = - c2 *  u(2,i,j,k) * tmp1
   108: |||||                    fjac(ib,4,3,kp) = - c2 *  u(3,i,j,k) * tmp1
   109: |||||                    fjac(ib,4,4,kp) = ( 2.0d+00 - c2 )
   110: |||||          >              *  u(4,i,j,k) * tmp1
   111: |||||                    fjac(ib,4,5,kp) = c2
   112: |||||     
   113: |||||   G                fjac(ib,5,1,kp) = ( c2 * 2.0d0 * square(i,j,k)
   114: |||||          >              - c1 * u(5,i,j,k) )
   115: |||||          >              * u(4,i,j,k) * tmp2
   116: |||||                    fjac(ib,5,2,kp) = - c2 * ( u(2,i,j,k)*u(4,i,j,k) )
   117: |||||          >              * tmp2
   118: |||||                    fjac(ib,5,3,kp) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
   119: |||||          >              * tmp2
   120: |||||   F                fjac(ib,5,4,kp) = c1 * ( u(5,i,j,k) * tmp1 )
   121: |||||          >              - c2
   122: |||||          >              * ( qs(i,j,k)
   123: |||||          >              + u(4,i,j,k)*u(4,i,j,k) * tmp2 )
   124: |||||                    fjac(ib,5,5,kp) = c1 * u(4,i,j,k) * tmp1
   125: |||||     
   126: |||||                    njac(ib,1,1,kp) = 0.0d+00
   127: |||||                    njac(ib,1,2,kp) = 0.0d+00
   128: |||||                    njac(ib,1,3,kp) = 0.0d+00
   129: |||||                    njac(ib,1,4,kp) = 0.0d+00
   130: |||||                    njac(ib,1,5,kp) = 0.0d+00
   131: |||||     
   132: |||||                    njac(ib,2,1,kp) = - c3c4 * tmp2 * u(2,i,j,k)
   133: |||||                    njac(ib,2,2,kp) =   c3c4 * tmp1
   134: |||||                    njac(ib,2,3,kp) =   0.0d+00
   135: |||||                    njac(ib,2,4,kp) =   0.0d+00
   136: |||||                    njac(ib,2,5,kp) =   0.0d+00
   137: |||||     
   138: |||||                    njac(ib,3,1,kp) = - c3c4 * tmp2 * u(3,i,j,k)
   139: |||||                    njac(ib,3,2,kp) =   0.0d+00
   140: |||||                    njac(ib,3,3,kp) =   c3c4 * tmp1
   141: |||||                    njac(ib,3,4,kp) =   0.0d+00
   142: |||||                    njac(ib,3,5,kp) =   0.0d+00
   143: |||||     
   144: |||||                    njac(ib,4,1,kp) = - con43 * c3c4 * tmp2 * u(4,i,j,k)
   145: |||||                    njac(ib,4,2,kp) =   0.0d+00
   146: |||||                    njac(ib,4,3,kp) =   0.0d+00
   147: |||||                    njac(ib,4,4,kp) =   con43 * c3 * c4 * tmp1
   148: |||||                    njac(ib,4,5,kp) =   0.0d+00
   149: |||||     
   150: |||||   F                njac(ib,5,1,kp) = - (  c3c4
   151: |||||          >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   152: |||||          >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   153: |||||          >              - ( con43 * c3c4
   154: |||||          >              - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   155: |||||          >              - c1345 * tmp2 * u(5,i,j,k)
   156: |||||     
   157: |||||                    njac(ib,5,2,kp) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
   158: |||||                    njac(ib,5,3,kp) = (  c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
   159: |||||                    njac(ib,5,4,kp) = ( con43 * c3c4
   160: |||||          >              - c1345 ) * tmp2 * u(4,i,j,k)
   161: |||||                    njac(ib,5,5,kp) = ( c1345 )* tmp1
   162: |||||     
   163: ||||V--               enddo
   164: |||+---               enddo
   165: |||       
   166: |||       c---------------------------------------------------------------------
   167: |||       c     now jacobians set, so form left hand side in z direction
   168: |||       c---------------------------------------------------------------------
   169: |||                   km = mod(kb + 2, 3)
   170: |||                   k = kk
   171: |||       !dir$ vector always
   172: |||V--->              do ib = 1, bsize
   173: ||||      
   174: ||||                     tmp1 = dt * tz1
   175: ||||                     tmp2 = dt * tz2
   176: ||||      
   177: ||||    F                lhsa(ib,1,1,1) = - tmp2 * fjac(ib,1,1,km)
   178: ||||           >              - tmp1 * njac(ib,1,1,km)
   179: ||||           >              - tmp1 * dz1
   180: ||||    F                lhsa(ib,1,2,1) = - tmp2 * fjac(ib,1,2,km)
   181: ||||           >              - tmp1 * njac(ib,1,2,km)
   182: ||||    F                lhsa(ib,1,3,1) = - tmp2 * fjac(ib,1,3,km)
   183: ||||           >              - tmp1 * njac(ib,1,3,km)
   184: ||||    F                lhsa(ib,1,4,1) = - tmp2 * fjac(ib,1,4,km)
   185: ||||           >              - tmp1 * njac(ib,1,4,km)
   186: ||||    F                lhsa(ib,1,5,1) = - tmp2 * fjac(ib,1,5,km)
   187: ||||           >              - tmp1 * njac(ib,1,5,km)
   188: ||||      
   189: ||||    F                lhsa(ib,2,1,1) = - tmp2 * fjac(ib,2,1,km)
   190: ||||           >              - tmp1 * njac(ib,2,1,km)
   191: ||||    F                lhsa(ib,2,2,1) = - tmp2 * fjac(ib,2,2,km)
   192: ||||           >              - tmp1 * njac(ib,2,2,km)
   193: ||||           >              - tmp1 * dz2
   194: ||||    F                lhsa(ib,2,3,1) = - tmp2 * fjac(ib,2,3,km)
   195: ||||           >              - tmp1 * njac(ib,2,3,km)
   196: ||||    F                lhsa(ib,2,4,1) = - tmp2 * fjac(ib,2,4,km)
   197: ||||           >              - tmp1 * njac(ib,2,4,km)
   198: ||||    F                lhsa(ib,2,5,1) = - tmp2 * fjac(ib,2,5,km)
   199: ||||           >              - tmp1 * njac(ib,2,5,km)
   200: ||||      
   201: ||||    F                lhsa(ib,3,1,1) = - tmp2 * fjac(ib,3,1,km)
   202: ||||           >              - tmp1 * njac(ib,3,1,km)
   203: ||||    F                lhsa(ib,3,2,1) = - tmp2 * fjac(ib,3,2,km)
   204: ||||           >              - tmp1 * njac(ib,3,2,km)
   205: ||||    F                lhsa(ib,3,3,1) = - tmp2 * fjac(ib,3,3,km)
   206: ||||           >              - tmp1 * njac(ib,3,3,km)
   207: ||||           >              - tmp1 * dz3
   208: ||||    F                lhsa(ib,3,4,1) = - tmp2 * fjac(ib,3,4,km)
   209: ||||           >              - tmp1 * njac(ib,3,4,km)
   210: ||||    F                lhsa(ib,3,5,1) = - tmp2 * fjac(ib,3,5,km)
   211: ||||           >              - tmp1 * njac(ib,3,5,km)
   212: ||||      
   213: ||||    F                lhsa(ib,4,1,1) = - tmp2 * fjac(ib,4,1,km)
   214: ||||           >              - tmp1 * njac(ib,4,1,km)
   215: ||||    F                lhsa(ib,4,2,1) = - tmp2 * fjac(ib,4,2,km)
   216: ||||           >              - tmp1 * njac(ib,4,2,km)
   217: ||||    F                lhsa(ib,4,3,1) = - tmp2 * fjac(ib,4,3,km)
   218: ||||           >              - tmp1 * njac(ib,4,3,km)
   219: ||||    F                lhsa(ib,4,4,1) = - tmp2 * fjac(ib,4,4,km)
   220: ||||           >              - tmp1 * njac(ib,4,4,km)
   221: ||||           >              - tmp1 * dz4
   222: ||||    F                lhsa(ib,4,5,1) = - tmp2 * fjac(ib,4,5,km)
   223: ||||           >              - tmp1 * njac(ib,4,5,km)
   224: ||||      
   225: ||||    F                lhsa(ib,5,1,1) = - tmp2 * fjac(ib,5,1,km)
   226: ||||           >              - tmp1 * njac(ib,5,1,km)
   227: ||||    F                lhsa(ib,5,2,1) = - tmp2 * fjac(ib,5,2,km)
   228: ||||           >              - tmp1 * njac(ib,5,2,km)
   229: ||||    F                lhsa(ib,5,3,1) = - tmp2 * fjac(ib,5,3,km)
   230: ||||           >              - tmp1 * njac(ib,5,3,km)
   231: ||||    F                lhsa(ib,5,4,1) = - tmp2 * fjac(ib,5,4,km)
   232: ||||           >              - tmp1 * njac(ib,5,4,km)
   233: ||||    F                lhsa(ib,5,5,1) = - tmp2 * fjac(ib,5,5,km)
   234: ||||           >              - tmp1 * njac(ib,5,5,km)
   235: ||||           >              - tmp1 * dz5
   236: ||||      
   237: ||||    F                lhsb(ib,1,1,1) = 1.0d+00
   238: ||||           >              + tmp1 * 2.0d+00 * njac(ib,1,1,kb)
   239: ||||           >              + tmp1 * 2.0d+00 * dz1
   240: ||||                     lhsb(ib,1,2,1) = tmp1 * 2.0d+00 * njac(ib,1,2,kb)
   241: ||||                     lhsb(ib,1,3,1) = tmp1 * 2.0d+00 * njac(ib,1,3,kb)
   242: ||||                     lhsb(ib,1,4,1) = tmp1 * 2.0d+00 * njac(ib,1,4,kb)
   243: ||||                     lhsb(ib,1,5,1) = tmp1 * 2.0d+00 * njac(ib,1,5,kb)
   244: ||||      
   245: ||||                     lhsb(ib,2,1,1) = tmp1 * 2.0d+00 * njac(ib,2,1,kb)
   246: ||||    F                lhsb(ib,2,2,1) = 1.0d+00
   247: ||||           >              + tmp1 * 2.0d+00 * njac(ib,2,2,kb)
   248: ||||           >              + tmp1 * 2.0d+00 * dz2
   249: ||||                     lhsb(ib,2,3,1) = tmp1 * 2.0d+00 * njac(ib,2,3,kb)
   250: ||||                     lhsb(ib,2,4,1) = tmp1 * 2.0d+00 * njac(ib,2,4,kb)
   251: ||||                     lhsb(ib,2,5,1) = tmp1 * 2.0d+00 * njac(ib,2,5,kb)
   252: ||||      
   253: ||||                     lhsb(ib,3,1,1) = tmp1 * 2.0d+00 * njac(ib,3,1,kb)
   254: ||||                     lhsb(ib,3,2,1) = tmp1 * 2.0d+00 * njac(ib,3,2,kb)
   255: ||||    F                lhsb(ib,3,3,1) = 1.0d+00
   256: ||||           >              + tmp1 * 2.0d+00 * njac(ib,3,3,kb)
   257: ||||           >              + tmp1 * 2.0d+00 * dz3
   258: ||||                     lhsb(ib,3,4,1) = tmp1 * 2.0d+00 * njac(ib,3,4,kb)
   259: ||||                     lhsb(ib,3,5,1) = tmp1 * 2.0d+00 * njac(ib,3,5,kb)
   260: ||||      
   261: ||||                     lhsb(ib,4,1,1) = tmp1 * 2.0d+00 * njac(ib,4,1,kb)
   262: ||||                     lhsb(ib,4,2,1) = tmp1 * 2.0d+00 * njac(ib,4,2,kb)
   263: ||||                     lhsb(ib,4,3,1) = tmp1 * 2.0d+00 * njac(ib,4,3,kb)
   264: ||||    F                lhsb(ib,4,4,1) = 1.0d+00
   265: ||||           >              + tmp1 * 2.0d+00 * njac(ib,4,4,kb)
   266: ||||           >              + tmp1 * 2.0d+00 * dz4
   267: ||||                     lhsb(ib,4,5,1) = tmp1 * 2.0d+00 * njac(ib,4,5,kb)
   268: ||||      
   269: ||||                     lhsb(ib,5,1,1) = tmp1 * 2.0d+00 * njac(ib,5,1,kb)
   270: ||||                     lhsb(ib,5,2,1) = tmp1 * 2.0d+00 * njac(ib,5,2,kb)
   271: ||||                     lhsb(ib,5,3,1) = tmp1 * 2.0d+00 * njac(ib,5,3,kb)
   272: ||||                     lhsb(ib,5,4,1) = tmp1 * 2.0d+00 * njac(ib,5,4,kb)
   273: ||||    F                lhsb(ib,5,5,1) = 1.0d+00
   274: ||||           >              + tmp1 * 2.0d+00 * njac(ib,5,5,kb)
   275: ||||           >              + tmp1 * 2.0d+00 * dz5
   276: ||||      
   277: ||||    F                lhsc(ib,1,1,k) =  tmp2 * fjac(ib,1,1,kp)
   278: ||||           >              - tmp1 * njac(ib,1,1,kp)
   279: ||||           >              - tmp1 * dz1
   280: ||||    F                lhsc(ib,1,2,k) =  tmp2 * fjac(ib,1,2,kp)
   281: ||||           >              - tmp1 * njac(ib,1,2,kp)
   282: ||||    F                lhsc(ib,1,3,k) =  tmp2 * fjac(ib,1,3,kp)
   283: ||||           >              - tmp1 * njac(ib,1,3,kp)
   284: ||||    F                lhsc(ib,1,4,k) =  tmp2 * fjac(ib,1,4,kp)
   285: ||||           >              - tmp1 * njac(ib,1,4,kp)
   286: ||||    F                lhsc(ib,1,5,k) =  tmp2 * fjac(ib,1,5,kp)
   287: ||||           >              - tmp1 * njac(ib,1,5,kp)
   288: ||||      
   289: ||||    F                lhsc(ib,2,1,k) =  tmp2 * fjac(ib,2,1,kp)
   290: ||||           >              - tmp1 * njac(ib,2,1,kp)
   291: ||||    F                lhsc(ib,2,2,k) =  tmp2 * fjac(ib,2,2,kp)
   292: ||||           >              - tmp1 * njac(ib,2,2,kp)
   293: ||||           >              - tmp1 * dz2
   294: ||||    F                lhsc(ib,2,3,k) =  tmp2 * fjac(ib,2,3,kp)
   295: ||||           >              - tmp1 * njac(ib,2,3,kp)
   296: ||||    F                lhsc(ib,2,4,k) =  tmp2 * fjac(ib,2,4,kp)
   297: ||||           >              - tmp1 * njac(ib,2,4,kp)
   298: ||||    F                lhsc(ib,2,5,k) =  tmp2 * fjac(ib,2,5,kp)
   299: ||||           >              - tmp1 * njac(ib,2,5,kp)
   300: ||||      
   301: ||||    F                lhsc(ib,3,1,k) =  tmp2 * fjac(ib,3,1,kp)
   302: ||||           >              - tmp1 * njac(ib,3,1,kp)
   303: ||||    F                lhsc(ib,3,2,k) =  tmp2 * fjac(ib,3,2,kp)
   304: ||||           >              - tmp1 * njac(ib,3,2,kp)
   305: ||||    F                lhsc(ib,3,3,k) =  tmp2 * fjac(ib,3,3,kp)
   306: ||||           >              - tmp1 * njac(ib,3,3,kp)
   307: ||||           >              - tmp1 * dz3
   308: ||||    F                lhsc(ib,3,4,k) =  tmp2 * fjac(ib,3,4,kp)
   309: ||||           >              - tmp1 * njac(ib,3,4,kp)
   310: ||||    F                lhsc(ib,3,5,k) =  tmp2 * fjac(ib,3,5,kp)
   311: ||||           >              - tmp1 * njac(ib,3,5,kp)
   312: ||||      
   313: ||||    F                lhsc(ib,4,1,k) =  tmp2 * fjac(ib,4,1,kp)
   314: ||||           >              - tmp1 * njac(ib,4,1,kp)
   315: ||||    F                lhsc(ib,4,2,k) =  tmp2 * fjac(ib,4,2,kp)
   316: ||||           >              - tmp1 * njac(ib,4,2,kp)
   317: ||||    F                lhsc(ib,4,3,k) =  tmp2 * fjac(ib,4,3,kp)
   318: ||||           >              - tmp1 * njac(ib,4,3,kp)
   319: ||||    F                lhsc(ib,4,4,k) =  tmp2 * fjac(ib,4,4,kp)
   320: ||||           >              - tmp1 * njac(ib,4,4,kp)
   321: ||||           >              - tmp1 * dz4
   322: ||||    F                lhsc(ib,4,5,k) =  tmp2 * fjac(ib,4,5,kp)
   323: ||||           >              - tmp1 * njac(ib,4,5,kp)
   324: ||||      
   325: ||||    F                lhsc(ib,5,1,k) =  tmp2 * fjac(ib,5,1,kp)
   326: ||||           >              - tmp1 * njac(ib,5,1,kp)
   327: ||||    F                lhsc(ib,5,2,k) =  tmp2 * fjac(ib,5,2,kp)
   328: ||||           >              - tmp1 * njac(ib,5,2,kp)
   329: ||||    F                lhsc(ib,5,3,k) =  tmp2 * fjac(ib,5,3,kp)
   330: ||||           >              - tmp1 * njac(ib,5,3,kp)
   331: ||||    F                lhsc(ib,5,4,k) =  tmp2 * fjac(ib,5,4,kp)
   332: ||||           >              - tmp1 * njac(ib,5,4,kp)
   333: ||||    F                lhsc(ib,5,5,k) =  tmp2 * fjac(ib,5,5,kp)
   334: ||||           >              - tmp1 * njac(ib,5,5,kp)
   335: ||||           >              - tmp1 * dz5
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
   349: |||       c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
   350: |||       c---------------------------------------------------------------------
   351: |||       
   352: |||       c---------------------------------------------------------------------
   353: |||       c     outer most do loops - sweeping in i direction
   354: |||       c---------------------------------------------------------------------
   355: |||       
   356: |||                   if (timeron) call timer_start(t_solsub)
   357: |||       c---------------------------------------------------------------------
   358: |||       c     multiply c(i,j,0) by b_inverse and copy back to c
   359: |||       c     multiply rhs(0) by b_inverse(0) and copy to rhs
   360: |||       c---------------------------------------------------------------------
   361: |||                   if (kk .eq. 1) then
   362: |||                   call binvcrhs( lhsb(1,1,1,0),
   363: |||            >                        lhsc(1,1,1,0),
   364: |||            >                        rhsx(1,1,0) )
   365: |||                   endif
   366: |||       
   367: |||       
   368: |||       c---------------------------------------------------------------------
   369: |||       c     begin inner most do loop
   370: |||       c     do all the elements of the cell unless last
   371: |||       c---------------------------------------------------------------------
   372: |||       
   373: |||       c---------------------------------------------------------------------
   374: |||       c     subtract A*lhs_vector(k-1) from lhs_vector(k)
   375: |||       c
   376: |||       c     rhs(k) = rhs(k) - A*rhs(k-1)
   377: |||       c---------------------------------------------------------------------
   378: |||                      call matvec_sub(lhsa(1,1,1,1),
   379: |||            >                         rhsx(1,1,k-1),rhsx(1,1,k))
   380: |||       
   381: |||       c---------------------------------------------------------------------
   382: |||       c     B(k) = B(k) - C(k-1)*A(k)
   383: |||       c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
   384: |||       c---------------------------------------------------------------------
   385: |||                      call matmul_sub(lhsa(1,1,1,1),
   386: |||            >                         lhsc(1,1,1,k-1),
   387: |||            >                         lhsb(1,1,1,1))
   388: |||       
   389: |||       c---------------------------------------------------------------------
   390: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   391: |||       c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
   392: |||       c---------------------------------------------------------------------
   393: |||                      call binvcrhs( lhsb(1,1,1,1),
   394: |||            >                        lhsc(1,1,1,k),
   395: |||            >                        rhsx(1,1,k) )
   396: |||       
   397: |||       
   398: |||       c---------------------------------------------------------------------
   399: |||       c     Now finish up special cases for last cell
   400: |||       c---------------------------------------------------------------------
   401: |||       
   402: |||                   if (kk .eq. ksize-1) then
   403: |||       c---------------------------------------------------------------------
   404: |||       c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
   405: |||       c---------------------------------------------------------------------
   406: |||                   call matvec_sub(lhsa(1,1,1,2),
   407: |||            >                         rhsx(1,1,ksize-1),rhsx(1,1,ksize))
   408: |||       
   409: |||       c---------------------------------------------------------------------
   410: |||       c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
   411: |||       c     call matmul_sub(aa,i,j,ksize,c,
   412: |||       c     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
   413: |||       c---------------------------------------------------------------------
   414: |||                   call matmul_sub(lhsa(1,1,1,2),
   415: |||            >                         lhsc(1,1,1,ksize-1),
   416: |||            >                         lhsb(1,1,1,2))
   417: |||       
   418: |||       c---------------------------------------------------------------------
   419: |||       c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
   420: |||       c---------------------------------------------------------------------
   421: |||                   call binvrhs( lhsb(1,1,1,2),
   422: |||            >                       rhsx(1,1,ksize) )
   423: |||                   endif
   424: |||       
   425: |||                   if (timeron) call timer_stop(t_solsub)
   426: |||       
   427: ||+----               enddo
   428: ||        
   429: ||        c---------------------------------------------------------------------
   430: ||        c---------------------------------------------------------------------
   431: ||        
   432: ||        c---------------------------------------------------------------------
   433: ||        c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
   434: ||        c     else assume U(ksize) is loaded in un pack backsub_info
   435: ||        c     so just use it
   436: ||        c     after call u(kstart) will be sent to next cell
   437: ||        c---------------------------------------------------------------------
   438: ||        
   439: ||+---->              do k=ksize-1,0,-1
   440: |||       !dir$ vector always
   441: |||+--->              do ib = 1, bsize
   442: ||||      !dir$ unroll
   443: ||||V-->                 do m=1,BLOCK_SIZE
   444: |||||     !dir$ unroll
   445: |||||*->                    do n=1,BLOCK_SIZE
   446: ||||||  F                      rhsx(ib,m,k) = rhsx(ib,m,k)
   447: ||||||         >                    - lhsc(ib,m,n,k)*rhsx(ib,n,k+1)
   448: |||||*-                     enddo
   449: ||||V--                  enddo
   450: |||+---               enddo
   451: ||+----               enddo
   452: ||        
   453: ||                    if (timeron) call timer_start(t_rdis1)
   454: ||+---->              do ib = 1, bsize
   455: |||                      i = ii+ib-1
   456: |||                      if (i < grid_points(1)-1) then
   457: |||V--->                 do k=0,ksize
   458: ||||                        rhs(1,i,j,k) = rhsx(ib,1,k)
   459: ||||                        rhs(2,i,j,k) = rhsx(ib,2,k)
   460: ||||                        rhs(3,i,j,k) = rhsx(ib,3,k)
   461: ||||                        rhs(4,i,j,k) = rhsx(ib,4,k)
   462: ||||                        rhs(5,i,j,k) = rhsx(ib,5,k)
   463: |||V---                  end do
   464: |||                      endif
   465: ||+----               end do
   466: ||                    if (timeron) call timer_stop(t_rdis1)
   467: ||        
   468: |+-----            enddo
   469: P------         enddo
   470:           !$omp end parallel
   471:                 if (timeron) call timer_stop(t_zsolve)
   472:           
   473:                 return
   474:                 end


