NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:17 2020
FILE NAME: x_solve.f

PROCEDURE NAME: X_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    33: inl(1212): Source for routine not found.: TIMER_START
    47: par(1801): Parallel routine generated.: X_SOLVE$1
    49: par(1803): Parallelized by "do".
    49: opt(2495): Loop nest collapsed by "collapse" clause.
    51: vec( 124): Iteration count is assumed. Iteration count.: 103
    51: vec( 101): Vectorized loop.
    65: vec( 128): Fused multiply-add operation applied.
    86: vec( 128): Fused multiply-add operation applied.
    89: vec( 128): Fused multiply-add operation applied.
   123: vec( 128): Fused multiply-add operation applied.
   140: inl(1212): Source for routine not found.: LHSINIT
   141: vec( 124): Iteration count is assumed. Iteration count.: 103
   141: vec( 101): Vectorized loop.
   146: vec( 128): Fused multiply-add operation applied.
   149: vec( 128): Fused multiply-add operation applied.
   151: vec( 128): Fused multiply-add operation applied.
   153: vec( 128): Fused multiply-add operation applied.
   155: vec( 128): Fused multiply-add operation applied.
   158: vec( 128): Fused multiply-add operation applied.
   160: vec( 128): Fused multiply-add operation applied.
   163: vec( 128): Fused multiply-add operation applied.
   165: vec( 128): Fused multiply-add operation applied.
   167: vec( 128): Fused multiply-add operation applied.
   170: vec( 128): Fused multiply-add operation applied.
   172: vec( 128): Fused multiply-add operation applied.
   174: vec( 128): Fused multiply-add operation applied.
   177: vec( 128): Fused multiply-add operation applied.
   179: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   184: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   188: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   194: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   198: vec( 128): Fused multiply-add operation applied.
   200: vec( 128): Fused multiply-add operation applied.
   202: vec( 128): Fused multiply-add operation applied.
   206: vec( 128): Fused multiply-add operation applied.
   215: vec( 128): Fused multiply-add operation applied.
   224: vec( 128): Fused multiply-add operation applied.
   233: vec( 128): Fused multiply-add operation applied.
   242: vec( 128): Fused multiply-add operation applied.
   246: vec( 128): Fused multiply-add operation applied.
   249: vec( 128): Fused multiply-add operation applied.
   251: vec( 128): Fused multiply-add operation applied.
   253: vec( 128): Fused multiply-add operation applied.
   255: vec( 128): Fused multiply-add operation applied.
   258: vec( 128): Fused multiply-add operation applied.
   260: vec( 128): Fused multiply-add operation applied.
   263: vec( 128): Fused multiply-add operation applied.
   265: vec( 128): Fused multiply-add operation applied.
   267: vec( 128): Fused multiply-add operation applied.
   270: vec( 128): Fused multiply-add operation applied.
   272: vec( 128): Fused multiply-add operation applied.
   274: vec( 128): Fused multiply-add operation applied.
   277: vec( 128): Fused multiply-add operation applied.
   279: vec( 128): Fused multiply-add operation applied.
   282: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   286: vec( 128): Fused multiply-add operation applied.
   288: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   294: vec( 128): Fused multiply-add operation applied.
   296: vec( 128): Fused multiply-add operation applied.
   298: vec( 128): Fused multiply-add operation applied.
   300: vec( 128): Fused multiply-add operation applied.
   302: vec( 128): Fused multiply-add operation applied.
   330: inl(1212): Source for routine not found.: BINVCRHS
   338: vec( 110): Vectorization obstructive procedure reference.: MATMUL_SUB
   338: vec( 110): Vectorization obstructive procedure reference.: BINVCRHS
   338: vec( 110): Vectorization obstructive procedure reference.: MATVEC_SUB
   338: vec( 103): Unvectorized loop.
   343: inl(1212): Source for routine not found.: MATVEC_SUB
   343: opt(1025): Reference to this procedure inhibits optimization.: MATVEC_SUB
   349: inl(1212): Source for routine not found.: MATMUL_SUB
   380: inl(1212): Source for routine not found.: BINVRHS
   382: inl(1212): Source for routine not found.: TIMER_STOP
   392: vec( 103): Unvectorized loop.
   392: vec( 113): Overhead of loop division is too large.
   393: opt(1418): Constant-length loop is expanded.
   394: opt(1418): Constant-length loop is expanded.
   395: opt(1037): Feedback of array elements.: RHS
   395: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:17 2020
FILE NAME: x_solve.f

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
    28:                 double precision tmp1, tmp2, tmp3
    29:           
    30:           c---------------------------------------------------------------------
    31:           c---------------------------------------------------------------------
    32:           
    33:                 if (timeron) call timer_start(t_xsolve)
    34:           
    35:           c---------------------------------------------------------------------
    36:           c---------------------------------------------------------------------
    37:           
    38:           c---------------------------------------------------------------------
    39:           c     This function computes the left hand side in the xi-direction
    40:           c---------------------------------------------------------------------
    41:           
    42:                 isize = grid_points(1)-1
    43:           
    44:           c---------------------------------------------------------------------
    45:           c     determine a (labeled f) and n jacobians
    46:           c---------------------------------------------------------------------
    47:           !$omp parallel do default(shared) shared(isize) collapse(2)
    48:           !$omp& private(i,j,k,m,n,tmp1,tmp2,tmp3)
    49: P------>        do k = 1, grid_points(3)-2
    50: |+----->           do j = 1, grid_points(2)-2
    51: ||V---->              do i = 0, isize
    52: |||       
    53: |||                      tmp1 = rho_i(i,j,k)
    54: |||                      tmp2 = tmp1 * tmp1
    55: |||                      tmp3 = tmp1 * tmp2
    56: |||       c---------------------------------------------------------------------
    57: |||       c
    58: |||       c---------------------------------------------------------------------
    59: |||                      fjac(1,1,i) = 0.0d+00
    60: |||                      fjac(1,2,i) = 1.0d+00
    61: |||                      fjac(1,3,i) = 0.0d+00
    62: |||                      fjac(1,4,i) = 0.0d+00
    63: |||                      fjac(1,5,i) = 0.0d+00
    64: |||       
    65: |||     F                fjac(2,1,i) = -(u(2,i,j,k) * tmp2 *
    66: |||            >              u(2,i,j,k))
    67: |||            >              + c2 * qs(i,j,k)
    68: |||                      fjac(2,2,i) = ( 2.0d+00 - c2 )
    69: |||            >              * ( u(2,i,j,k) / u(1,i,j,k) )
    70: |||                      fjac(2,3,i) = - c2 * ( u(3,i,j,k) * tmp1 )
    71: |||                      fjac(2,4,i) = - c2 * ( u(4,i,j,k) * tmp1 )
    72: |||                      fjac(2,5,i) = c2
    73: |||       
    74: |||                      fjac(3,1,i) = - ( u(2,i,j,k)*u(3,i,j,k) ) * tmp2
    75: |||                      fjac(3,2,i) = u(3,i,j,k) * tmp1
    76: |||                      fjac(3,3,i) = u(2,i,j,k) * tmp1
    77: |||                      fjac(3,4,i) = 0.0d+00
    78: |||                      fjac(3,5,i) = 0.0d+00
    79: |||       
    80: |||                      fjac(4,1,i) = - ( u(2,i,j,k)*u(4,i,j,k) ) * tmp2
    81: |||                      fjac(4,2,i) = u(4,i,j,k) * tmp1
    82: |||                      fjac(4,3,i) = 0.0d+00
    83: |||                      fjac(4,4,i) = u(2,i,j,k) * tmp1
    84: |||                      fjac(4,5,i) = 0.0d+00
    85: |||       
    86: |||     F                fjac(5,1,i) = ( c2 * 2.0d0 * square(i,j,k)
    87: |||            >              - c1 * u(5,i,j,k) )
    88: |||            >              * ( u(2,i,j,k) * tmp2 )
    89: |||     F                fjac(5,2,i) = c1 *  u(5,i,j,k) * tmp1
    90: |||            >              - c2
    91: |||            >              * ( u(2,i,j,k)*u(2,i,j,k) * tmp2
    92: |||            >              + qs(i,j,k) )
    93: |||                      fjac(5,3,i) = - c2 * ( u(3,i,j,k)*u(2,i,j,k) )
    94: |||            >              * tmp2
    95: |||                      fjac(5,4,i) = - c2 * ( u(4,i,j,k)*u(2,i,j,k) )
    96: |||            >              * tmp2
    97: |||                      fjac(5,5,i) = c1 * ( u(2,i,j,k) * tmp1 )
    98: |||       
    99: |||                      njac(1,1,i) = 0.0d+00
   100: |||                      njac(1,2,i) = 0.0d+00
   101: |||                      njac(1,3,i) = 0.0d+00
   102: |||                      njac(1,4,i) = 0.0d+00
   103: |||                      njac(1,5,i) = 0.0d+00
   104: |||       
   105: |||                      njac(2,1,i) = - con43 * c3c4 * tmp2 * u(2,i,j,k)
   106: |||                      njac(2,2,i) =   con43 * c3c4 * tmp1
   107: |||                      njac(2,3,i) =   0.0d+00
   108: |||                      njac(2,4,i) =   0.0d+00
   109: |||                      njac(2,5,i) =   0.0d+00
   110: |||       
   111: |||                      njac(3,1,i) = - c3c4 * tmp2 * u(3,i,j,k)
   112: |||                      njac(3,2,i) =   0.0d+00
   113: |||                      njac(3,3,i) =   c3c4 * tmp1
   114: |||                      njac(3,4,i) =   0.0d+00
   115: |||                      njac(3,5,i) =   0.0d+00
   116: |||       
   117: |||                      njac(4,1,i) = - c3c4 * tmp2 * u(4,i,j,k)
   118: |||                      njac(4,2,i) =   0.0d+00
   119: |||                      njac(4,3,i) =   0.0d+00
   120: |||                      njac(4,4,i) =   c3c4 * tmp1
   121: |||                      njac(4,5,i) =   0.0d+00
   122: |||       
   123: |||     F                njac(5,1,i) = - ( con43 * c3c4
   124: |||            >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   125: |||            >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   126: |||            >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   127: |||            >              - c1345 * tmp2 * u(5,i,j,k)
   128: |||       
   129: |||                      njac(5,2,i) = ( con43 * c3c4
   130: |||            >              - c1345 ) * tmp2 * u(2,i,j,k)
   131: |||                      njac(5,3,i) = ( c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
   132: |||                      njac(5,4,i) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
   133: |||                      njac(5,5,i) = ( c1345 ) * tmp1
   134: |||       
   135: ||V----               enddo
   136: ||        
   137: ||        c---------------------------------------------------------------------
   138: ||        c     now jacobians set, so form left hand side in x direction
   139: ||        c---------------------------------------------------------------------
   140: ||                    call lhsinit(lhs, isize)
   141: ||V---->              do i = 1, isize-1
   142: |||       
   143: |||                      tmp1 = dt * tx1
   144: |||                      tmp2 = dt * tx2
   145: |||       
   146: |||     F                lhs(1,1,aa,i) = - tmp2 * fjac(1,1,i-1)
   147: |||            >              - tmp1 * njac(1,1,i-1)
   148: |||            >              - tmp1 * dx1
   149: |||     F                lhs(1,2,aa,i) = - tmp2 * fjac(1,2,i-1)
   150: |||            >              - tmp1 * njac(1,2,i-1)
   151: |||     F                lhs(1,3,aa,i) = - tmp2 * fjac(1,3,i-1)
   152: |||            >              - tmp1 * njac(1,3,i-1)
   153: |||     F                lhs(1,4,aa,i) = - tmp2 * fjac(1,4,i-1)
   154: |||            >              - tmp1 * njac(1,4,i-1)
   155: |||     F                lhs(1,5,aa,i) = - tmp2 * fjac(1,5,i-1)
   156: |||            >              - tmp1 * njac(1,5,i-1)
   157: |||       
   158: |||     F                lhs(2,1,aa,i) = - tmp2 * fjac(2,1,i-1)
   159: |||            >              - tmp1 * njac(2,1,i-1)
   160: |||     F                lhs(2,2,aa,i) = - tmp2 * fjac(2,2,i-1)
   161: |||            >              - tmp1 * njac(2,2,i-1)
   162: |||            >              - tmp1 * dx2
   163: |||     F                lhs(2,3,aa,i) = - tmp2 * fjac(2,3,i-1)
   164: |||            >              - tmp1 * njac(2,3,i-1)
   165: |||     F                lhs(2,4,aa,i) = - tmp2 * fjac(2,4,i-1)
   166: |||            >              - tmp1 * njac(2,4,i-1)
   167: |||     F                lhs(2,5,aa,i) = - tmp2 * fjac(2,5,i-1)
   168: |||            >              - tmp1 * njac(2,5,i-1)
   169: |||       
   170: |||     F                lhs(3,1,aa,i) = - tmp2 * fjac(3,1,i-1)
   171: |||            >              - tmp1 * njac(3,1,i-1)
   172: |||     F                lhs(3,2,aa,i) = - tmp2 * fjac(3,2,i-1)
   173: |||            >              - tmp1 * njac(3,2,i-1)
   174: |||     F                lhs(3,3,aa,i) = - tmp2 * fjac(3,3,i-1)
   175: |||            >              - tmp1 * njac(3,3,i-1)
   176: |||            >              - tmp1 * dx3
   177: |||     F                lhs(3,4,aa,i) = - tmp2 * fjac(3,4,i-1)
   178: |||            >              - tmp1 * njac(3,4,i-1)
   179: |||     F                lhs(3,5,aa,i) = - tmp2 * fjac(3,5,i-1)
   180: |||            >              - tmp1 * njac(3,5,i-1)
   181: |||       
   182: |||     F                lhs(4,1,aa,i) = - tmp2 * fjac(4,1,i-1)
   183: |||            >              - tmp1 * njac(4,1,i-1)
   184: |||     F                lhs(4,2,aa,i) = - tmp2 * fjac(4,2,i-1)
   185: |||            >              - tmp1 * njac(4,2,i-1)
   186: |||     F                lhs(4,3,aa,i) = - tmp2 * fjac(4,3,i-1)
   187: |||            >              - tmp1 * njac(4,3,i-1)
   188: |||     F                lhs(4,4,aa,i) = - tmp2 * fjac(4,4,i-1)
   189: |||            >              - tmp1 * njac(4,4,i-1)
   190: |||            >              - tmp1 * dx4
   191: |||     F                lhs(4,5,aa,i) = - tmp2 * fjac(4,5,i-1)
   192: |||            >              - tmp1 * njac(4,5,i-1)
   193: |||       
   194: |||     F                lhs(5,1,aa,i) = - tmp2 * fjac(5,1,i-1)
   195: |||            >              - tmp1 * njac(5,1,i-1)
   196: |||     F                lhs(5,2,aa,i) = - tmp2 * fjac(5,2,i-1)
   197: |||            >              - tmp1 * njac(5,2,i-1)
   198: |||     F                lhs(5,3,aa,i) = - tmp2 * fjac(5,3,i-1)
   199: |||            >              - tmp1 * njac(5,3,i-1)
   200: |||     F                lhs(5,4,aa,i) = - tmp2 * fjac(5,4,i-1)
   201: |||            >              - tmp1 * njac(5,4,i-1)
   202: |||     F                lhs(5,5,aa,i) = - tmp2 * fjac(5,5,i-1)
   203: |||            >              - tmp1 * njac(5,5,i-1)
   204: |||            >              - tmp1 * dx5
   205: |||       
   206: |||     F                lhs(1,1,bb,i) = 1.0d+00
   207: |||            >              + tmp1 * 2.0d+00 * njac(1,1,i)
   208: |||            >              + tmp1 * 2.0d+00 * dx1
   209: |||                      lhs(1,2,bb,i) = tmp1 * 2.0d+00 * njac(1,2,i)
   210: |||                      lhs(1,3,bb,i) = tmp1 * 2.0d+00 * njac(1,3,i)
   211: |||                      lhs(1,4,bb,i) = tmp1 * 2.0d+00 * njac(1,4,i)
   212: |||                      lhs(1,5,bb,i) = tmp1 * 2.0d+00 * njac(1,5,i)
   213: |||       
   214: |||                      lhs(2,1,bb,i) = tmp1 * 2.0d+00 * njac(2,1,i)
   215: |||     F                lhs(2,2,bb,i) = 1.0d+00
   216: |||            >              + tmp1 * 2.0d+00 * njac(2,2,i)
   217: |||            >              + tmp1 * 2.0d+00 * dx2
   218: |||                      lhs(2,3,bb,i) = tmp1 * 2.0d+00 * njac(2,3,i)
   219: |||                      lhs(2,4,bb,i) = tmp1 * 2.0d+00 * njac(2,4,i)
   220: |||                      lhs(2,5,bb,i) = tmp1 * 2.0d+00 * njac(2,5,i)
   221: |||       
   222: |||                      lhs(3,1,bb,i) = tmp1 * 2.0d+00 * njac(3,1,i)
   223: |||                      lhs(3,2,bb,i) = tmp1 * 2.0d+00 * njac(3,2,i)
   224: |||     F                lhs(3,3,bb,i) = 1.0d+00
   225: |||            >              + tmp1 * 2.0d+00 * njac(3,3,i)
   226: |||            >              + tmp1 * 2.0d+00 * dx3
   227: |||                      lhs(3,4,bb,i) = tmp1 * 2.0d+00 * njac(3,4,i)
   228: |||                      lhs(3,5,bb,i) = tmp1 * 2.0d+00 * njac(3,5,i)
   229: |||       
   230: |||                      lhs(4,1,bb,i) = tmp1 * 2.0d+00 * njac(4,1,i)
   231: |||                      lhs(4,2,bb,i) = tmp1 * 2.0d+00 * njac(4,2,i)
   232: |||                      lhs(4,3,bb,i) = tmp1 * 2.0d+00 * njac(4,3,i)
   233: |||     F                lhs(4,4,bb,i) = 1.0d+00
   234: |||            >              + tmp1 * 2.0d+00 * njac(4,4,i)
   235: |||            >              + tmp1 * 2.0d+00 * dx4
   236: |||                      lhs(4,5,bb,i) = tmp1 * 2.0d+00 * njac(4,5,i)
   237: |||       
   238: |||                      lhs(5,1,bb,i) = tmp1 * 2.0d+00 * njac(5,1,i)
   239: |||                      lhs(5,2,bb,i) = tmp1 * 2.0d+00 * njac(5,2,i)
   240: |||                      lhs(5,3,bb,i) = tmp1 * 2.0d+00 * njac(5,3,i)
   241: |||                      lhs(5,4,bb,i) = tmp1 * 2.0d+00 * njac(5,4,i)
   242: |||     F                lhs(5,5,bb,i) = 1.0d+00
   243: |||            >              + tmp1 * 2.0d+00 * njac(5,5,i)
   244: |||            >              + tmp1 * 2.0d+00 * dx5
   245: |||       
   246: |||     F                lhs(1,1,cc,i) =  tmp2 * fjac(1,1,i+1)
   247: |||            >              - tmp1 * njac(1,1,i+1)
   248: |||            >              - tmp1 * dx1
   249: |||     F                lhs(1,2,cc,i) =  tmp2 * fjac(1,2,i+1)
   250: |||            >              - tmp1 * njac(1,2,i+1)
   251: |||     F                lhs(1,3,cc,i) =  tmp2 * fjac(1,3,i+1)
   252: |||            >              - tmp1 * njac(1,3,i+1)
   253: |||     F                lhs(1,4,cc,i) =  tmp2 * fjac(1,4,i+1)
   254: |||            >              - tmp1 * njac(1,4,i+1)
   255: |||     F                lhs(1,5,cc,i) =  tmp2 * fjac(1,5,i+1)
   256: |||            >              - tmp1 * njac(1,5,i+1)
   257: |||       
   258: |||     F                lhs(2,1,cc,i) =  tmp2 * fjac(2,1,i+1)
   259: |||            >              - tmp1 * njac(2,1,i+1)
   260: |||     F                lhs(2,2,cc,i) =  tmp2 * fjac(2,2,i+1)
   261: |||            >              - tmp1 * njac(2,2,i+1)
   262: |||            >              - tmp1 * dx2
   263: |||     F                lhs(2,3,cc,i) =  tmp2 * fjac(2,3,i+1)
   264: |||            >              - tmp1 * njac(2,3,i+1)
   265: |||     F                lhs(2,4,cc,i) =  tmp2 * fjac(2,4,i+1)
   266: |||            >              - tmp1 * njac(2,4,i+1)
   267: |||     F                lhs(2,5,cc,i) =  tmp2 * fjac(2,5,i+1)
   268: |||            >              - tmp1 * njac(2,5,i+1)
   269: |||       
   270: |||     F                lhs(3,1,cc,i) =  tmp2 * fjac(3,1,i+1)
   271: |||            >              - tmp1 * njac(3,1,i+1)
   272: |||     F                lhs(3,2,cc,i) =  tmp2 * fjac(3,2,i+1)
   273: |||            >              - tmp1 * njac(3,2,i+1)
   274: |||     F                lhs(3,3,cc,i) =  tmp2 * fjac(3,3,i+1)
   275: |||            >              - tmp1 * njac(3,3,i+1)
   276: |||            >              - tmp1 * dx3
   277: |||     F                lhs(3,4,cc,i) =  tmp2 * fjac(3,4,i+1)
   278: |||            >              - tmp1 * njac(3,4,i+1)
   279: |||     F                lhs(3,5,cc,i) =  tmp2 * fjac(3,5,i+1)
   280: |||            >              - tmp1 * njac(3,5,i+1)
   281: |||       
   282: |||     F                lhs(4,1,cc,i) =  tmp2 * fjac(4,1,i+1)
   283: |||            >              - tmp1 * njac(4,1,i+1)
   284: |||     F                lhs(4,2,cc,i) =  tmp2 * fjac(4,2,i+1)
   285: |||            >              - tmp1 * njac(4,2,i+1)
   286: |||     F                lhs(4,3,cc,i) =  tmp2 * fjac(4,3,i+1)
   287: |||            >              - tmp1 * njac(4,3,i+1)
   288: |||     F                lhs(4,4,cc,i) =  tmp2 * fjac(4,4,i+1)
   289: |||            >              - tmp1 * njac(4,4,i+1)
   290: |||            >              - tmp1 * dx4
   291: |||     F                lhs(4,5,cc,i) =  tmp2 * fjac(4,5,i+1)
   292: |||            >              - tmp1 * njac(4,5,i+1)
   293: |||       
   294: |||     F                lhs(5,1,cc,i) =  tmp2 * fjac(5,1,i+1)
   295: |||            >              - tmp1 * njac(5,1,i+1)
   296: |||     F                lhs(5,2,cc,i) =  tmp2 * fjac(5,2,i+1)
   297: |||            >              - tmp1 * njac(5,2,i+1)
   298: |||     F                lhs(5,3,cc,i) =  tmp2 * fjac(5,3,i+1)
   299: |||            >              - tmp1 * njac(5,3,i+1)
   300: |||     F                lhs(5,4,cc,i) =  tmp2 * fjac(5,4,i+1)
   301: |||            >              - tmp1 * njac(5,4,i+1)
   302: |||     F                lhs(5,5,cc,i) =  tmp2 * fjac(5,5,i+1)
   303: |||            >              - tmp1 * njac(5,5,i+1)
   304: |||            >              - tmp1 * dx5
   305: |||       
   306: ||V----               enddo
   307: ||        
   308: ||        c---------------------------------------------------------------------
   309: ||        c---------------------------------------------------------------------
   310: ||        
   311: ||        c---------------------------------------------------------------------
   312: ||        c     performs guaussian elimination on this cell.
   313: ||        c
   314: ||        c     assumes that unpacking routines for non-first cells
   315: ||        c     preload C' and rhs' from previous cell.
   316: ||        c
   317: ||        c     assumed send happens outside this routine, but that
   318: ||        c     c'(IMAX) and rhs'(IMAX) will be sent to next cell
   319: ||        c---------------------------------------------------------------------
   320: ||        
   321: ||        c---------------------------------------------------------------------
   322: ||        c     outer most do loops - sweeping in i direction
   323: ||        c---------------------------------------------------------------------
   324: ||        
   325: ||                    if (timeron) call timer_start(t_solsub)
   326: ||        c---------------------------------------------------------------------
   327: ||        c     multiply c(0,j,k) by b_inverse and copy back to c
   328: ||        c     multiply rhs(0) by b_inverse(0) and copy to rhs
   329: ||        c---------------------------------------------------------------------
   330: ||                    call binvcrhs( lhs(1,1,bb,0),
   331: ||             >                        lhs(1,1,cc,0),
   332: ||             >                        rhs(1,0,j,k) )
   333: ||        
   334: ||        c---------------------------------------------------------------------
   335: ||        c     begin inner most do loop
   336: ||        c     do all the elements of the cell unless last
   337: ||        c---------------------------------------------------------------------
   338: ||+---->              do i=1,isize-1
   339: |||       
   340: |||       c---------------------------------------------------------------------
   341: |||       c     rhs(i) = rhs(i) - A*rhs(i-1)
   342: |||       c---------------------------------------------------------------------
   343: |||                      call matvec_sub(lhs(1,1,aa,i),
   344: |||            >                         rhs(1,i-1,j,k),rhs(1,i,j,k))
   345: |||       
   346: |||       c---------------------------------------------------------------------
   347: |||       c     B(i) = B(i) - C(i-1)*A(i)
   348: |||       c---------------------------------------------------------------------
   349: |||                      call matmul_sub(lhs(1,1,aa,i),
   350: |||            >                         lhs(1,1,cc,i-1),
   351: |||            >                         lhs(1,1,bb,i))
   352: |||       
   353: |||       
   354: |||       c---------------------------------------------------------------------
   355: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   356: |||       c     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
   357: |||       c---------------------------------------------------------------------
   358: |||                      call binvcrhs( lhs(1,1,bb,i),
   359: |||            >                        lhs(1,1,cc,i),
   360: |||            >                        rhs(1,i,j,k) )
   361: |||       
   362: ||+----               enddo
   363: ||        
   364: ||        c---------------------------------------------------------------------
   365: ||        c     rhs(isize) = rhs(isize) - A*rhs(isize-1)
   366: ||        c---------------------------------------------------------------------
   367: ||                    call matvec_sub(lhs(1,1,aa,isize),
   368: ||             >                         rhs(1,isize-1,j,k),rhs(1,isize,j,k))
   369: ||        
   370: ||        c---------------------------------------------------------------------
   371: ||        c     B(isize) = B(isize) - C(isize-1)*A(isize)
   372: ||        c---------------------------------------------------------------------
   373: ||                    call matmul_sub(lhs(1,1,aa,isize),
   374: ||             >                         lhs(1,1,cc,isize-1),
   375: ||             >                         lhs(1,1,bb,isize))
   376: ||        
   377: ||        c---------------------------------------------------------------------
   378: ||        c     multiply rhs() by b_inverse() and copy to rhs
   379: ||        c---------------------------------------------------------------------
   380: ||                    call binvrhs( lhs(1,1,bb,isize),
   381: ||             >                       rhs(1,isize,j,k) )
   382: ||                    if (timeron) call timer_stop(t_solsub)
   383: ||        
   384: ||        
   385: ||        c---------------------------------------------------------------------
   386: ||        c     back solve: if last cell, then generate U(isize)=rhs(isize)
   387: ||        c     else assume U(isize) is loaded in un pack backsub_info
   388: ||        c     so just use it
   389: ||        c     after call u(istart) will be sent to next cell
   390: ||        c---------------------------------------------------------------------
   391: ||        
   392: ||+---->              do i=isize-1,0,-1
   393: |||*--->                 do m=1,BLOCK_SIZE
   394: ||||*-->                    do n=1,BLOCK_SIZE
   395: |||||                          rhs(m,i,j,k) = rhs(m,i,j,k)
   396: |||||          >                    - lhs(m,n,cc,i)*rhs(n,i+1,j,k)
   397: ||||*--                     enddo
   398: |||*---                  enddo
   399: ||+----               enddo
   400: ||        
   401: |+-----            enddo
   402: P------         enddo
   403:                 if (timeron) call timer_stop(t_xsolve)
   404:           
   405:                 return
   406:                 end


