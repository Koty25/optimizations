NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: x_solve_blk.f

PROCEDURE NAME: X_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    48: par(1801): Parallel routine generated.: X_SOLVE$1
    54: par(1803): Parallelized by "do".
    54: opt(2495): Loop nest collapsed by "collapse" clause.
    59: vec( 101): Vectorized loop.
    61: opt(1034): Multiple store conflict.: RHSX
    62: opt(1034): Multiple store conflict.: RHSX
    63: opt(1034): Multiple store conflict.: RHSX
    64: opt(1034): Multiple store conflict.: RHSX
    73: opt(1135): Outer loop conditionally executes inner loop.
    73: vec( 110): Vectorization obstructive procedure reference.: BINVRHS
    73: vec( 110): Vectorization obstructive procedure reference.: MATMUL_SUB
    73: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
    73: vec( 110): Vectorization obstructive procedure reference.: BINVCRHS
    73: vec( 110): Vectorization obstructive procedure reference.: MATVEC_SUB
    73: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
    80: vec( 101): Vectorized loop.
    89: opt(1034): Multiple store conflict.: FJAC
    95: vec( 128): Fused multiply-add operation applied.
   116: vec( 128): Fused multiply-add operation applied.
   119: vec( 128): Fused multiply-add operation applied.
   153: vec( 128): Fused multiply-add operation applied.
   166: opt(1062): Too many data dependency problems.
   174: vec( 101): Vectorized loop.
   179: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   184: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   188: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   193: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   198: vec( 128): Fused multiply-add operation applied.
   200: vec( 128): Fused multiply-add operation applied.
   203: vec( 128): Fused multiply-add operation applied.
   205: vec( 128): Fused multiply-add operation applied.
   207: vec( 128): Fused multiply-add operation applied.
   210: vec( 128): Fused multiply-add operation applied.
   212: vec( 128): Fused multiply-add operation applied.
   215: vec( 128): Fused multiply-add operation applied.
   217: vec( 128): Fused multiply-add operation applied.
   219: vec( 128): Fused multiply-add operation applied.
   221: vec( 128): Fused multiply-add operation applied.
   224: vec( 128): Fused multiply-add operation applied.
   227: vec( 128): Fused multiply-add operation applied.
   229: vec( 128): Fused multiply-add operation applied.
   231: vec( 128): Fused multiply-add operation applied.
   233: vec( 128): Fused multiply-add operation applied.
   235: vec( 128): Fused multiply-add operation applied.
   239: vec( 128): Fused multiply-add operation applied.
   248: vec( 128): Fused multiply-add operation applied.
   257: vec( 128): Fused multiply-add operation applied.
   266: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   279: vec( 128): Fused multiply-add operation applied.
   282: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   286: vec( 128): Fused multiply-add operation applied.
   288: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   293: vec( 128): Fused multiply-add operation applied.
   296: vec( 128): Fused multiply-add operation applied.
   298: vec( 128): Fused multiply-add operation applied.
   300: vec( 128): Fused multiply-add operation applied.
   303: vec( 128): Fused multiply-add operation applied.
   305: vec( 128): Fused multiply-add operation applied.
   307: vec( 128): Fused multiply-add operation applied.
   310: vec( 128): Fused multiply-add operation applied.
   312: vec( 128): Fused multiply-add operation applied.
   315: vec( 128): Fused multiply-add operation applied.
   317: vec( 128): Fused multiply-add operation applied.
   319: vec( 128): Fused multiply-add operation applied.
   321: vec( 128): Fused multiply-add operation applied.
   324: vec( 128): Fused multiply-add operation applied.
   327: vec( 128): Fused multiply-add operation applied.
   329: vec( 128): Fused multiply-add operation applied.
   331: vec( 128): Fused multiply-add operation applied.
   333: vec( 128): Fused multiply-add operation applied.
   335: vec( 128): Fused multiply-add operation applied.
   432: vec( 101): Vectorized loop.
   434: opt(1418): Constant-length loop is expanded.
   435: opt(1036): Potential feedback - use directive or compiler option if OK.: RHSX
   435: opt(1037): Feedback of array elements.: RHSX
   435: vec( 128): Fused multiply-add operation applied.
   443: opt(1135): Outer loop conditionally executes inner loop.
   446: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: x_solve_blk.f

PROCEDURE NAME: X_SOLVE
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine x_solve
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c
    12:           c     Performs line solves in X direction by first factoring
    13:           c     the block-tridiagonal matrix into an upper triangular matrix,
    14:           c     and then performing back substitution to solve for the unknow
    15:           c     vectors of each line.
    16:           c
    17:           c     Make sure we treat elements zero to cell_size in the direction
    18:           c     of the sweep.
    19:           c
    20:           c---------------------------------------------------------------------
    21:           
    22:                 use bt_data
    23:                 use work_lhs
    24:           
    25:                 implicit none
    26:           
    27:                 integer i,j,k,m,n,isize
    28:                 integer ii,im,ip,ib,jj,jb
    29:                 double precision tmp1, tmp2, tmp3
    30:           
    31:           c---------------------------------------------------------------------
    32:           c---------------------------------------------------------------------
    33:           
    34:                 if (timeron) call timer_start(t_xsolve)
    35:           
    36:           c---------------------------------------------------------------------
    37:           c---------------------------------------------------------------------
    38:           
    39:           c---------------------------------------------------------------------
    40:           c     This function computes the left hand side in the xi-direction
    41:           c---------------------------------------------------------------------
    42:           
    43:                 isize = grid_points(1)-1
    44:           
    45:           c---------------------------------------------------------------------
    46:           c     determine a (labeled f) and n jacobians
    47:           c---------------------------------------------------------------------
    48:           !$omp parallel default(shared) shared(isize)
    49:           !$omp& private(i,j,k,m,n,ii,im,ip,ib,jj,jb,tmp1,tmp2,tmp3)
    50:           
    51:                 call lhsinit(isize)
    52:           
    53:           !$omp do collapse(2)
    54: P------>        do k = 1, grid_points(3)-2
    55: |+----->           do jj = 1, grid_points(2)-2, bsize
    56: ||        
    57: ||                    if (timeron) call timer_start(t_rdis1)
    58: ||+---->              do i=0,isize
    59: |||V--->              do jb = 1, bsize
    60: ||||                     j = min(jj+jb-1, grid_points(2)-2)
    61: ||||    G                rhsx(jb,1,i) = rhs(1,i,j,k)
    62: ||||    G                rhsx(jb,2,i) = rhs(2,i,j,k)
    63: ||||    G                rhsx(jb,3,i) = rhs(3,i,j,k)
    64: ||||    G                rhsx(jb,4,i) = rhs(4,i,j,k)
    65: ||||    G                rhsx(jb,5,i) = rhs(5,i,j,k)
    66: |||V---               end do
    67: ||+----               end do
    68: ||                    if (timeron) call timer_stop(t_rdis1)
    69: ||        
    70: ||                    call lhsinit(0)
    71: ||        
    72: ||                    ib = 0
    73: ||+---->              do ii = 1, isize-1
    74: |||                   ib = mod(ib + 1, 3)
    75: |||                   im = min(2*ii - 3, 1)     ! -1 or 1
    76: |||                   ip = mod(ib + im, 3) - 1
    77: |||       
    78: |||+--->              do i = ii+im, ii+1
    79: ||||                  ip = ip + 1
    80: ||||V-->              do jb = 1, bsize
    81: |||||                    j = min(jj+jb-1, grid_points(2)-2)
    82: |||||     
    83: |||||   G                tmp1 = rho_i(i,j,k)
    84: |||||                    tmp2 = tmp1 * tmp1
    85: |||||                    tmp3 = tmp1 * tmp2
    86: |||||     c---------------------------------------------------------------------
    87: |||||     c
    88: |||||     c---------------------------------------------------------------------
    89: |||||                    fjac(jb,1,1,ip) = 0.0d+00
    90: |||||                    fjac(jb,1,2,ip) = 1.0d+00
    91: |||||                    fjac(jb,1,3,ip) = 0.0d+00
    92: |||||                    fjac(jb,1,4,ip) = 0.0d+00
    93: |||||                    fjac(jb,1,5,ip) = 0.0d+00
    94: |||||     
    95: |||||   G                fjac(jb,2,1,ip) = -(u(2,i,j,k) * tmp2 *
    96: |||||          >              u(2,i,j,k))
    97: |||||          >              + c2 * qs(i,j,k)
    98: |||||   G                fjac(jb,2,2,ip) = ( 2.0d+00 - c2 )
    99: |||||          >              * ( u(2,i,j,k) / u(1,i,j,k) )
   100: |||||   G                fjac(jb,2,3,ip) = - c2 * ( u(3,i,j,k) * tmp1 )
   101: |||||   G                fjac(jb,2,4,ip) = - c2 * ( u(4,i,j,k) * tmp1 )
   102: |||||                    fjac(jb,2,5,ip) = c2
   103: |||||     
   104: |||||                    fjac(jb,3,1,ip) = - ( u(2,i,j,k)*u(3,i,j,k) ) * tmp2
   105: |||||                    fjac(jb,3,2,ip) = u(3,i,j,k) * tmp1
   106: |||||                    fjac(jb,3,3,ip) = u(2,i,j,k) * tmp1
   107: |||||                    fjac(jb,3,4,ip) = 0.0d+00
   108: |||||                    fjac(jb,3,5,ip) = 0.0d+00
   109: |||||     
   110: |||||                    fjac(jb,4,1,ip) = - ( u(2,i,j,k)*u(4,i,j,k) ) * tmp2
   111: |||||                    fjac(jb,4,2,ip) = u(4,i,j,k) * tmp1
   112: |||||                    fjac(jb,4,3,ip) = 0.0d+00
   113: |||||                    fjac(jb,4,4,ip) = u(2,i,j,k) * tmp1
   114: |||||                    fjac(jb,4,5,ip) = 0.0d+00
   115: |||||     
   116: |||||   G                fjac(jb,5,1,ip) = ( c2 * 2.0d0 * square(i,j,k)
   117: |||||          >              - c1 * u(5,i,j,k) )
   118: |||||          >              * ( u(2,i,j,k) * tmp2 )
   119: |||||   F                fjac(jb,5,2,ip) = c1 *  u(5,i,j,k) * tmp1
   120: |||||          >              - c2
   121: |||||          >              * ( u(2,i,j,k)*u(2,i,j,k) * tmp2
   122: |||||          >              + qs(i,j,k) )
   123: |||||                    fjac(jb,5,3,ip) = - c2 * ( u(3,i,j,k)*u(2,i,j,k) )
   124: |||||          >              * tmp2
   125: |||||                    fjac(jb,5,4,ip) = - c2 * ( u(4,i,j,k)*u(2,i,j,k) )
   126: |||||          >              * tmp2
   127: |||||                    fjac(jb,5,5,ip) = c1 * ( u(2,i,j,k) * tmp1 )
   128: |||||     
   129: |||||                    njac(jb,1,1,ip) = 0.0d+00
   130: |||||                    njac(jb,1,2,ip) = 0.0d+00
   131: |||||                    njac(jb,1,3,ip) = 0.0d+00
   132: |||||                    njac(jb,1,4,ip) = 0.0d+00
   133: |||||                    njac(jb,1,5,ip) = 0.0d+00
   134: |||||     
   135: |||||                    njac(jb,2,1,ip) = - con43 * c3c4 * tmp2 * u(2,i,j,k)
   136: |||||                    njac(jb,2,2,ip) =   con43 * c3c4 * tmp1
   137: |||||                    njac(jb,2,3,ip) =   0.0d+00
   138: |||||                    njac(jb,2,4,ip) =   0.0d+00
   139: |||||                    njac(jb,2,5,ip) =   0.0d+00
   140: |||||     
   141: |||||                    njac(jb,3,1,ip) = - c3c4 * tmp2 * u(3,i,j,k)
   142: |||||                    njac(jb,3,2,ip) =   0.0d+00
   143: |||||                    njac(jb,3,3,ip) =   c3c4 * tmp1
   144: |||||                    njac(jb,3,4,ip) =   0.0d+00
   145: |||||                    njac(jb,3,5,ip) =   0.0d+00
   146: |||||     
   147: |||||                    njac(jb,4,1,ip) = - c3c4 * tmp2 * u(4,i,j,k)
   148: |||||                    njac(jb,4,2,ip) =   0.0d+00
   149: |||||                    njac(jb,4,3,ip) =   0.0d+00
   150: |||||                    njac(jb,4,4,ip) =   c3c4 * tmp1
   151: |||||                    njac(jb,4,5,ip) =   0.0d+00
   152: |||||     
   153: |||||   F                njac(jb,5,1,ip) = - ( con43 * c3c4
   154: |||||          >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   155: |||||          >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   156: |||||          >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   157: |||||          >              - c1345 * tmp2 * u(5,i,j,k)
   158: |||||     
   159: |||||                    njac(jb,5,2,ip) = ( con43 * c3c4
   160: |||||          >              - c1345 ) * tmp2 * u(2,i,j,k)
   161: |||||                    njac(jb,5,3,ip) = ( c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
   162: |||||                    njac(jb,5,4,ip) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
   163: |||||                    njac(jb,5,5,ip) = ( c1345 ) * tmp1
   164: |||||     
   165: ||||V--               enddo
   166: |||+---               enddo
   167: |||       
   168: |||       c---------------------------------------------------------------------
   169: |||       c     now jacobians set, so form left hand side in x direction
   170: |||       c---------------------------------------------------------------------
   171: |||                   im = mod(ib + 2, 3)
   172: |||                   i = ii
   173: |||       !dir$ vector always
   174: |||V--->              do jb = 1, bsize
   175: ||||      
   176: ||||                     tmp1 = dt * tx1
   177: ||||                     tmp2 = dt * tx2
   178: ||||      
   179: ||||    F                lhsa(jb,1,1,1) = - tmp2 * fjac(jb,1,1,im)
   180: ||||           >              - tmp1 * njac(jb,1,1,im)
   181: ||||           >              - tmp1 * dx1
   182: ||||    F                lhsa(jb,1,2,1) = - tmp2 * fjac(jb,1,2,im)
   183: ||||           >              - tmp1 * njac(jb,1,2,im)
   184: ||||    F                lhsa(jb,1,3,1) = - tmp2 * fjac(jb,1,3,im)
   185: ||||           >              - tmp1 * njac(jb,1,3,im)
   186: ||||    F                lhsa(jb,1,4,1) = - tmp2 * fjac(jb,1,4,im)
   187: ||||           >              - tmp1 * njac(jb,1,4,im)
   188: ||||    F                lhsa(jb,1,5,1) = - tmp2 * fjac(jb,1,5,im)
   189: ||||           >              - tmp1 * njac(jb,1,5,im)
   190: ||||      
   191: ||||    F                lhsa(jb,2,1,1) = - tmp2 * fjac(jb,2,1,im)
   192: ||||           >              - tmp1 * njac(jb,2,1,im)
   193: ||||    F                lhsa(jb,2,2,1) = - tmp2 * fjac(jb,2,2,im)
   194: ||||           >              - tmp1 * njac(jb,2,2,im)
   195: ||||           >              - tmp1 * dx2
   196: ||||    F                lhsa(jb,2,3,1) = - tmp2 * fjac(jb,2,3,im)
   197: ||||           >              - tmp1 * njac(jb,2,3,im)
   198: ||||    F                lhsa(jb,2,4,1) = - tmp2 * fjac(jb,2,4,im)
   199: ||||           >              - tmp1 * njac(jb,2,4,im)
   200: ||||    F                lhsa(jb,2,5,1) = - tmp2 * fjac(jb,2,5,im)
   201: ||||           >              - tmp1 * njac(jb,2,5,im)
   202: ||||      
   203: ||||    F                lhsa(jb,3,1,1) = - tmp2 * fjac(jb,3,1,im)
   204: ||||           >              - tmp1 * njac(jb,3,1,im)
   205: ||||    F                lhsa(jb,3,2,1) = - tmp2 * fjac(jb,3,2,im)
   206: ||||           >              - tmp1 * njac(jb,3,2,im)
   207: ||||    F                lhsa(jb,3,3,1) = - tmp2 * fjac(jb,3,3,im)
   208: ||||           >              - tmp1 * njac(jb,3,3,im)
   209: ||||           >              - tmp1 * dx3
   210: ||||    F                lhsa(jb,3,4,1) = - tmp2 * fjac(jb,3,4,im)
   211: ||||           >              - tmp1 * njac(jb,3,4,im)
   212: ||||    F                lhsa(jb,3,5,1) = - tmp2 * fjac(jb,3,5,im)
   213: ||||           >              - tmp1 * njac(jb,3,5,im)
   214: ||||      
   215: ||||    F                lhsa(jb,4,1,1) = - tmp2 * fjac(jb,4,1,im)
   216: ||||           >              - tmp1 * njac(jb,4,1,im)
   217: ||||    F                lhsa(jb,4,2,1) = - tmp2 * fjac(jb,4,2,im)
   218: ||||           >              - tmp1 * njac(jb,4,2,im)
   219: ||||    F                lhsa(jb,4,3,1) = - tmp2 * fjac(jb,4,3,im)
   220: ||||           >              - tmp1 * njac(jb,4,3,im)
   221: ||||    F                lhsa(jb,4,4,1) = - tmp2 * fjac(jb,4,4,im)
   222: ||||           >              - tmp1 * njac(jb,4,4,im)
   223: ||||           >              - tmp1 * dx4
   224: ||||    F                lhsa(jb,4,5,1) = - tmp2 * fjac(jb,4,5,im)
   225: ||||           >              - tmp1 * njac(jb,4,5,im)
   226: ||||      
   227: ||||    F                lhsa(jb,5,1,1) = - tmp2 * fjac(jb,5,1,im)
   228: ||||           >              - tmp1 * njac(jb,5,1,im)
   229: ||||    F                lhsa(jb,5,2,1) = - tmp2 * fjac(jb,5,2,im)
   230: ||||           >              - tmp1 * njac(jb,5,2,im)
   231: ||||    F                lhsa(jb,5,3,1) = - tmp2 * fjac(jb,5,3,im)
   232: ||||           >              - tmp1 * njac(jb,5,3,im)
   233: ||||    F                lhsa(jb,5,4,1) = - tmp2 * fjac(jb,5,4,im)
   234: ||||           >              - tmp1 * njac(jb,5,4,im)
   235: ||||    F                lhsa(jb,5,5,1) = - tmp2 * fjac(jb,5,5,im)
   236: ||||           >              - tmp1 * njac(jb,5,5,im)
   237: ||||           >              - tmp1 * dx5
   238: ||||      
   239: ||||    F                lhsb(jb,1,1,1) = 1.0d+00
   240: ||||           >              + tmp1 * 2.0d+00 * njac(jb,1,1,ib)
   241: ||||           >              + tmp1 * 2.0d+00 * dx1
   242: ||||                     lhsb(jb,1,2,1) = tmp1 * 2.0d+00 * njac(jb,1,2,ib)
   243: ||||                     lhsb(jb,1,3,1) = tmp1 * 2.0d+00 * njac(jb,1,3,ib)
   244: ||||                     lhsb(jb,1,4,1) = tmp1 * 2.0d+00 * njac(jb,1,4,ib)
   245: ||||                     lhsb(jb,1,5,1) = tmp1 * 2.0d+00 * njac(jb,1,5,ib)
   246: ||||      
   247: ||||                     lhsb(jb,2,1,1) = tmp1 * 2.0d+00 * njac(jb,2,1,ib)
   248: ||||    F                lhsb(jb,2,2,1) = 1.0d+00
   249: ||||           >              + tmp1 * 2.0d+00 * njac(jb,2,2,ib)
   250: ||||           >              + tmp1 * 2.0d+00 * dx2
   251: ||||                     lhsb(jb,2,3,1) = tmp1 * 2.0d+00 * njac(jb,2,3,ib)
   252: ||||                     lhsb(jb,2,4,1) = tmp1 * 2.0d+00 * njac(jb,2,4,ib)
   253: ||||                     lhsb(jb,2,5,1) = tmp1 * 2.0d+00 * njac(jb,2,5,ib)
   254: ||||      
   255: ||||                     lhsb(jb,3,1,1) = tmp1 * 2.0d+00 * njac(jb,3,1,ib)
   256: ||||                     lhsb(jb,3,2,1) = tmp1 * 2.0d+00 * njac(jb,3,2,ib)
   257: ||||    F                lhsb(jb,3,3,1) = 1.0d+00
   258: ||||           >              + tmp1 * 2.0d+00 * njac(jb,3,3,ib)
   259: ||||           >              + tmp1 * 2.0d+00 * dx3
   260: ||||                     lhsb(jb,3,4,1) = tmp1 * 2.0d+00 * njac(jb,3,4,ib)
   261: ||||                     lhsb(jb,3,5,1) = tmp1 * 2.0d+00 * njac(jb,3,5,ib)
   262: ||||      
   263: ||||                     lhsb(jb,4,1,1) = tmp1 * 2.0d+00 * njac(jb,4,1,ib)
   264: ||||                     lhsb(jb,4,2,1) = tmp1 * 2.0d+00 * njac(jb,4,2,ib)
   265: ||||                     lhsb(jb,4,3,1) = tmp1 * 2.0d+00 * njac(jb,4,3,ib)
   266: ||||    F                lhsb(jb,4,4,1) = 1.0d+00
   267: ||||           >              + tmp1 * 2.0d+00 * njac(jb,4,4,ib)
   268: ||||           >              + tmp1 * 2.0d+00 * dx4
   269: ||||                     lhsb(jb,4,5,1) = tmp1 * 2.0d+00 * njac(jb,4,5,ib)
   270: ||||      
   271: ||||                     lhsb(jb,5,1,1) = tmp1 * 2.0d+00 * njac(jb,5,1,ib)
   272: ||||                     lhsb(jb,5,2,1) = tmp1 * 2.0d+00 * njac(jb,5,2,ib)
   273: ||||                     lhsb(jb,5,3,1) = tmp1 * 2.0d+00 * njac(jb,5,3,ib)
   274: ||||                     lhsb(jb,5,4,1) = tmp1 * 2.0d+00 * njac(jb,5,4,ib)
   275: ||||    F                lhsb(jb,5,5,1) = 1.0d+00
   276: ||||           >              + tmp1 * 2.0d+00 * njac(jb,5,5,ib)
   277: ||||           >              + tmp1 * 2.0d+00 * dx5
   278: ||||      
   279: ||||    F                lhsc(jb,1,1,i) =  tmp2 * fjac(jb,1,1,ip)
   280: ||||           >              - tmp1 * njac(jb,1,1,ip)
   281: ||||           >              - tmp1 * dx1
   282: ||||    F                lhsc(jb,1,2,i) =  tmp2 * fjac(jb,1,2,ip)
   283: ||||           >              - tmp1 * njac(jb,1,2,ip)
   284: ||||    F                lhsc(jb,1,3,i) =  tmp2 * fjac(jb,1,3,ip)
   285: ||||           >              - tmp1 * njac(jb,1,3,ip)
   286: ||||    F                lhsc(jb,1,4,i) =  tmp2 * fjac(jb,1,4,ip)
   287: ||||           >              - tmp1 * njac(jb,1,4,ip)
   288: ||||    F                lhsc(jb,1,5,i) =  tmp2 * fjac(jb,1,5,ip)
   289: ||||           >              - tmp1 * njac(jb,1,5,ip)
   290: ||||      
   291: ||||    F                lhsc(jb,2,1,i) =  tmp2 * fjac(jb,2,1,ip)
   292: ||||           >              - tmp1 * njac(jb,2,1,ip)
   293: ||||    F                lhsc(jb,2,2,i) =  tmp2 * fjac(jb,2,2,ip)
   294: ||||           >              - tmp1 * njac(jb,2,2,ip)
   295: ||||           >              - tmp1 * dx2
   296: ||||    F                lhsc(jb,2,3,i) =  tmp2 * fjac(jb,2,3,ip)
   297: ||||           >              - tmp1 * njac(jb,2,3,ip)
   298: ||||    F                lhsc(jb,2,4,i) =  tmp2 * fjac(jb,2,4,ip)
   299: ||||           >              - tmp1 * njac(jb,2,4,ip)
   300: ||||    F                lhsc(jb,2,5,i) =  tmp2 * fjac(jb,2,5,ip)
   301: ||||           >              - tmp1 * njac(jb,2,5,ip)
   302: ||||      
   303: ||||    F                lhsc(jb,3,1,i) =  tmp2 * fjac(jb,3,1,ip)
   304: ||||           >              - tmp1 * njac(jb,3,1,ip)
   305: ||||    F                lhsc(jb,3,2,i) =  tmp2 * fjac(jb,3,2,ip)
   306: ||||           >              - tmp1 * njac(jb,3,2,ip)
   307: ||||    F                lhsc(jb,3,3,i) =  tmp2 * fjac(jb,3,3,ip)
   308: ||||           >              - tmp1 * njac(jb,3,3,ip)
   309: ||||           >              - tmp1 * dx3
   310: ||||    F                lhsc(jb,3,4,i) =  tmp2 * fjac(jb,3,4,ip)
   311: ||||           >              - tmp1 * njac(jb,3,4,ip)
   312: ||||    F                lhsc(jb,3,5,i) =  tmp2 * fjac(jb,3,5,ip)
   313: ||||           >              - tmp1 * njac(jb,3,5,ip)
   314: ||||      
   315: ||||    F                lhsc(jb,4,1,i) =  tmp2 * fjac(jb,4,1,ip)
   316: ||||           >              - tmp1 * njac(jb,4,1,ip)
   317: ||||    F                lhsc(jb,4,2,i) =  tmp2 * fjac(jb,4,2,ip)
   318: ||||           >              - tmp1 * njac(jb,4,2,ip)
   319: ||||    F                lhsc(jb,4,3,i) =  tmp2 * fjac(jb,4,3,ip)
   320: ||||           >              - tmp1 * njac(jb,4,3,ip)
   321: ||||    F                lhsc(jb,4,4,i) =  tmp2 * fjac(jb,4,4,ip)
   322: ||||           >              - tmp1 * njac(jb,4,4,ip)
   323: ||||           >              - tmp1 * dx4
   324: ||||    F                lhsc(jb,4,5,i) =  tmp2 * fjac(jb,4,5,ip)
   325: ||||           >              - tmp1 * njac(jb,4,5,ip)
   326: ||||      
   327: ||||    F                lhsc(jb,5,1,i) =  tmp2 * fjac(jb,5,1,ip)
   328: ||||           >              - tmp1 * njac(jb,5,1,ip)
   329: ||||    F                lhsc(jb,5,2,i) =  tmp2 * fjac(jb,5,2,ip)
   330: ||||           >              - tmp1 * njac(jb,5,2,ip)
   331: ||||    F                lhsc(jb,5,3,i) =  tmp2 * fjac(jb,5,3,ip)
   332: ||||           >              - tmp1 * njac(jb,5,3,ip)
   333: ||||    F                lhsc(jb,5,4,i) =  tmp2 * fjac(jb,5,4,ip)
   334: ||||           >              - tmp1 * njac(jb,5,4,ip)
   335: ||||    F                lhsc(jb,5,5,i) =  tmp2 * fjac(jb,5,5,ip)
   336: ||||           >              - tmp1 * njac(jb,5,5,ip)
   337: ||||           >              - tmp1 * dx5
   338: ||||      
   339: |||V---               enddo
   340: |||       
   341: |||       c---------------------------------------------------------------------
   342: |||       c---------------------------------------------------------------------
   343: |||       
   344: |||       c---------------------------------------------------------------------
   345: |||       c     performs guaussian elimination on this cell.
   346: |||       c
   347: |||       c     assumes that unpacking routines for non-first cells
   348: |||       c     preload C' and rhs' from previous cell.
   349: |||       c
   350: |||       c     assumed send happens outside this routine, but that
   351: |||       c     c'(IMAX) and rhs'(IMAX) will be sent to next cell
   352: |||       c---------------------------------------------------------------------
   353: |||       
   354: |||       c---------------------------------------------------------------------
   355: |||       c     outer most do loops - sweeping in i direction
   356: |||       c---------------------------------------------------------------------
   357: |||       
   358: |||                   if (timeron) call timer_start(t_solsub)
   359: |||       c---------------------------------------------------------------------
   360: |||       c     multiply c(0,j,k) by b_inverse and copy back to c
   361: |||       c     multiply rhs(0) by b_inverse(0) and copy to rhs
   362: |||       c---------------------------------------------------------------------
   363: |||                   if (ii .eq. 1) then
   364: |||                      call binvcrhs( lhsb(1,1,1,0),
   365: |||            >                        lhsc(1,1,1,0),
   366: |||            >                        rhsx(1,1,0) )
   367: |||                   endif
   368: |||       
   369: |||       c---------------------------------------------------------------------
   370: |||       c     begin inner most do loop
   371: |||       c     do all the elements of the cell unless last
   372: |||       c---------------------------------------------------------------------
   373: |||       
   374: |||       c---------------------------------------------------------------------
   375: |||       c     rhs(i) = rhs(i) - A*rhs(i-1)
   376: |||       c---------------------------------------------------------------------
   377: |||                      call matvec_sub(lhsa(1,1,1,1),
   378: |||            >                         rhsx(1,1,i-1),rhsx(1,1,i))
   379: |||       
   380: |||       c---------------------------------------------------------------------
   381: |||       c     B(i) = B(i) - C(i-1)*A(i)
   382: |||       c---------------------------------------------------------------------
   383: |||                      call matmul_sub(lhsa(1,1,1,1),
   384: |||            >                         lhsc(1,1,1,i-1),
   385: |||            >                         lhsb(1,1,1,1))
   386: |||       
   387: |||       
   388: |||       c---------------------------------------------------------------------
   389: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   390: |||       c     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
   391: |||       c---------------------------------------------------------------------
   392: |||                      call binvcrhs( lhsb(1,1,1,1),
   393: |||            >                        lhsc(1,1,1,i),
   394: |||            >                        rhsx(1,1,i) )
   395: |||       
   396: |||       
   397: |||                   if (ii .eq. isize-1) then
   398: |||       c---------------------------------------------------------------------
   399: |||       c     rhs(isize) = rhs(isize) - A*rhs(isize-1)
   400: |||       c---------------------------------------------------------------------
   401: |||                      call matvec_sub(lhsa(1,1,1,2),
   402: |||            >                         rhsx(1,1,isize-1),rhsx(1,1,isize))
   403: |||       
   404: |||       c---------------------------------------------------------------------
   405: |||       c     B(isize) = B(isize) - C(isize-1)*A(isize)
   406: |||       c---------------------------------------------------------------------
   407: |||                      call matmul_sub(lhsa(1,1,1,2),
   408: |||            >                         lhsc(1,1,1,isize-1),
   409: |||            >                         lhsb(1,1,1,2))
   410: |||       
   411: |||       c---------------------------------------------------------------------
   412: |||       c     multiply rhs() by b_inverse() and copy to rhs
   413: |||       c---------------------------------------------------------------------
   414: |||                      call binvrhs( lhsb(1,1,1,2),
   415: |||            >                       rhsx(1,1,isize) )
   416: |||                   endif
   417: |||                   if (timeron) call timer_stop(t_solsub)
   418: |||       
   419: ||+----               enddo
   420: ||        
   421: ||        c---------------------------------------------------------------------
   422: ||        c     back solve: if last cell, then generate U(isize)=rhs(isize)
   423: ||        c     else assume U(isize) is loaded in un pack backsub_info
   424: ||        c     so just use it
   425: ||        c     after call u(istart) will be sent to next cell
   426: ||        c---------------------------------------------------------------------
   427: ||        
   428: ||+---->              do i=isize-1,0,-1
   429: |||       !dir$ vector always
   430: |||+--->              do jb=1,bsize
   431: ||||      !dir$ unroll
   432: ||||V-->                 do m=1,BLOCK_SIZE
   433: |||||     !dir$ unroll
   434: |||||*->                    do n=1,BLOCK_SIZE
   435: ||||||  F                      rhsx(jb,m,i) = rhsx(jb,m,i)
   436: ||||||         >                    - lhsc(jb,m,n,i)*rhsx(jb,n,i+1)
   437: |||||*-                     enddo
   438: ||||V--                  enddo
   439: |||+---               enddo
   440: ||+----               enddo
   441: ||        
   442: ||                    if (timeron) call timer_start(t_rdis1)
   443: ||+---->              do jb = 1, bsize
   444: |||                      j = jj+jb-1
   445: |||                      if (j .lt. grid_points(2)-1) then
   446: |||V--->                 do i=0,isize
   447: ||||                        rhs(1,i,j,k) = rhsx(jb,1,i)
   448: ||||                        rhs(2,i,j,k) = rhsx(jb,2,i)
   449: ||||                        rhs(3,i,j,k) = rhsx(jb,3,i)
   450: ||||                        rhs(4,i,j,k) = rhsx(jb,4,i)
   451: ||||                        rhs(5,i,j,k) = rhsx(jb,5,i)
   452: |||V---                  end do
   453: |||                      endif
   454: ||+----               end do
   455: ||                    if (timeron) call timer_stop(t_rdis1)
   456: ||        
   457: |+-----            enddo
   458: P------         enddo
   459:           !$omp end parallel
   460:                 if (timeron) call timer_stop(t_xsolve)
   461:           
   462:                 return
   463:                 end


