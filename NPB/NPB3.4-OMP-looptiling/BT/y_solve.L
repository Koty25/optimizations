NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:17 2020
FILE NAME: y_solve.f

PROCEDURE NAME: Y_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    30: inl(1212): Source for routine not found.: TIMER_START
    45: par(1801): Parallel routine generated.: Y_SOLVE$1
    47: par(1803): Parallelized by "do".
    47: opt(2495): Loop nest collapsed by "collapse" clause.
    49: vec( 124): Iteration count is assumed. Iteration count.: 103
    49: vec( 101): Vectorized loop.
    68: vec( 128): Fused multiply-add operation applied.
    83: vec( 128): Fused multiply-add operation applied.
    88: vec( 128): Fused multiply-add operation applied.
   120: vec( 128): Fused multiply-add operation applied.
   138: inl(1212): Source for routine not found.: LHSINIT
   139: vec( 124): Iteration count is assumed. Iteration count.: 103
   139: vec( 101): Vectorized loop.
   144: vec( 128): Fused multiply-add operation applied.
   147: vec( 128): Fused multiply-add operation applied.
   149: vec( 128): Fused multiply-add operation applied.
   151: vec( 128): Fused multiply-add operation applied.
   153: vec( 128): Fused multiply-add operation applied.
   156: vec( 128): Fused multiply-add operation applied.
   158: vec( 128): Fused multiply-add operation applied.
   161: vec( 128): Fused multiply-add operation applied.
   163: vec( 128): Fused multiply-add operation applied.
   165: vec( 128): Fused multiply-add operation applied.
   168: vec( 128): Fused multiply-add operation applied.
   170: vec( 128): Fused multiply-add operation applied.
   172: vec( 128): Fused multiply-add operation applied.
   175: vec( 128): Fused multiply-add operation applied.
   177: vec( 128): Fused multiply-add operation applied.
   180: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   184: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   189: vec( 128): Fused multiply-add operation applied.
   192: vec( 128): Fused multiply-add operation applied.
   194: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   198: vec( 128): Fused multiply-add operation applied.
   200: vec( 128): Fused multiply-add operation applied.
   204: vec( 128): Fused multiply-add operation applied.
   213: vec( 128): Fused multiply-add operation applied.
   222: vec( 128): Fused multiply-add operation applied.
   231: vec( 128): Fused multiply-add operation applied.
   240: vec( 128): Fused multiply-add operation applied.
   244: vec( 128): Fused multiply-add operation applied.
   247: vec( 128): Fused multiply-add operation applied.
   249: vec( 128): Fused multiply-add operation applied.
   251: vec( 128): Fused multiply-add operation applied.
   253: vec( 128): Fused multiply-add operation applied.
   256: vec( 128): Fused multiply-add operation applied.
   258: vec( 128): Fused multiply-add operation applied.
   261: vec( 128): Fused multiply-add operation applied.
   263: vec( 128): Fused multiply-add operation applied.
   265: vec( 128): Fused multiply-add operation applied.
   268: vec( 128): Fused multiply-add operation applied.
   270: vec( 128): Fused multiply-add operation applied.
   272: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   277: vec( 128): Fused multiply-add operation applied.
   280: vec( 128): Fused multiply-add operation applied.
   282: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   286: vec( 128): Fused multiply-add operation applied.
   289: vec( 128): Fused multiply-add operation applied.
   292: vec( 128): Fused multiply-add operation applied.
   294: vec( 128): Fused multiply-add operation applied.
   296: vec( 128): Fused multiply-add operation applied.
   298: vec( 128): Fused multiply-add operation applied.
   300: vec( 128): Fused multiply-add operation applied.
   324: inl(1212): Source for routine not found.: BINVCRHS
   332: vec( 110): Vectorization obstructive procedure reference.: MATMUL_SUB
   332: vec( 110): Vectorization obstructive procedure reference.: BINVCRHS
   332: vec( 110): Vectorization obstructive procedure reference.: MATVEC_SUB
   332: vec( 103): Unvectorized loop.
   339: inl(1212): Source for routine not found.: MATVEC_SUB
   339: opt(1025): Reference to this procedure inhibits optimization.: MATVEC_SUB
   345: inl(1212): Source for routine not found.: MATMUL_SUB
   378: inl(1212): Source for routine not found.: BINVRHS
   380: inl(1212): Source for routine not found.: TIMER_STOP
   390: vec( 103): Unvectorized loop.
   390: vec( 113): Overhead of loop division is too large.
   391: opt(1418): Constant-length loop is expanded.
   392: opt(1418): Constant-length loop is expanded.
   393: opt(1037): Feedback of array elements.: RHS
   393: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:17 2020
FILE NAME: y_solve.f

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
    25:                 double precision tmp1, tmp2, tmp3
    26:           
    27:           c---------------------------------------------------------------------
    28:           c---------------------------------------------------------------------
    29:           
    30:                 if (timeron) call timer_start(t_ysolve)
    31:           
    32:           c---------------------------------------------------------------------
    33:           c---------------------------------------------------------------------
    34:           
    35:           c---------------------------------------------------------------------
    36:           c     This function computes the left hand side for the three y-factors
    37:           c---------------------------------------------------------------------
    38:           
    39:                 jsize = grid_points(2)-1
    40:           
    41:           c---------------------------------------------------------------------
    42:           c     Compute the indices for storing the tri-diagonal matrix;
    43:           c     determine a (labeled f) and n jacobians for cell c
    44:           c---------------------------------------------------------------------
    45:           !$omp parallel do default(shared) shared(jsize) collapse(2)
    46:           !$omp& private(i,j,k,m,n,tmp1,tmp2,tmp3)
    47: P------>        do k = 1, grid_points(3)-2
    48: |+----->           do i = 1, grid_points(1)-2
    49: ||V---->              do j = 0, jsize
    50: |||       
    51: |||                      tmp1 = rho_i(i,j,k)
    52: |||                      tmp2 = tmp1 * tmp1
    53: |||                      tmp3 = tmp1 * tmp2
    54: |||       
    55: |||                      fjac(1,1,j) = 0.0d+00
    56: |||                      fjac(1,2,j) = 0.0d+00
    57: |||                      fjac(1,3,j) = 1.0d+00
    58: |||                      fjac(1,4,j) = 0.0d+00
    59: |||                      fjac(1,5,j) = 0.0d+00
    60: |||       
    61: |||                      fjac(2,1,j) = - ( u(2,i,j,k)*u(3,i,j,k) )
    62: |||            >              * tmp2
    63: |||                      fjac(2,2,j) = u(3,i,j,k) * tmp1
    64: |||                      fjac(2,3,j) = u(2,i,j,k) * tmp1
    65: |||                      fjac(2,4,j) = 0.0d+00
    66: |||                      fjac(2,5,j) = 0.0d+00
    67: |||       
    68: |||     F                fjac(3,1,j) = - ( u(3,i,j,k)*u(3,i,j,k)*tmp2)
    69: |||            >              + c2 * qs(i,j,k)
    70: |||                      fjac(3,2,j) = - c2 *  u(2,i,j,k) * tmp1
    71: |||                      fjac(3,3,j) = ( 2.0d+00 - c2 )
    72: |||            >              *  u(3,i,j,k) * tmp1
    73: |||                      fjac(3,4,j) = - c2 * u(4,i,j,k) * tmp1
    74: |||                      fjac(3,5,j) = c2
    75: |||       
    76: |||                      fjac(4,1,j) = - ( u(3,i,j,k)*u(4,i,j,k) )
    77: |||            >              * tmp2
    78: |||                      fjac(4,2,j) = 0.0d+00
    79: |||                      fjac(4,3,j) = u(4,i,j,k) * tmp1
    80: |||                      fjac(4,4,j) = u(3,i,j,k) * tmp1
    81: |||                      fjac(4,5,j) = 0.0d+00
    82: |||       
    83: |||     F                fjac(5,1,j) = ( c2 * 2.0d0 * square(i,j,k)
    84: |||            >              - c1 * u(5,i,j,k) )
    85: |||            >              * u(3,i,j,k) * tmp2
    86: |||                      fjac(5,2,j) = - c2 * u(2,i,j,k)*u(3,i,j,k)
    87: |||            >              * tmp2
    88: |||     F                fjac(5,3,j) = c1 * u(5,i,j,k) * tmp1
    89: |||            >              - c2
    90: |||            >              * ( qs(i,j,k)
    91: |||            >              + u(3,i,j,k)*u(3,i,j,k) * tmp2 )
    92: |||                      fjac(5,4,j) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
    93: |||            >              * tmp2
    94: |||                      fjac(5,5,j) = c1 * u(3,i,j,k) * tmp1
    95: |||       
    96: |||                      njac(1,1,j) = 0.0d+00
    97: |||                      njac(1,2,j) = 0.0d+00
    98: |||                      njac(1,3,j) = 0.0d+00
    99: |||                      njac(1,4,j) = 0.0d+00
   100: |||                      njac(1,5,j) = 0.0d+00
   101: |||       
   102: |||                      njac(2,1,j) = - c3c4 * tmp2 * u(2,i,j,k)
   103: |||                      njac(2,2,j) =   c3c4 * tmp1
   104: |||                      njac(2,3,j) =   0.0d+00
   105: |||                      njac(2,4,j) =   0.0d+00
   106: |||                      njac(2,5,j) =   0.0d+00
   107: |||       
   108: |||                      njac(3,1,j) = - con43 * c3c4 * tmp2 * u(3,i,j,k)
   109: |||                      njac(3,2,j) =   0.0d+00
   110: |||                      njac(3,3,j) =   con43 * c3c4 * tmp1
   111: |||                      njac(3,4,j) =   0.0d+00
   112: |||                      njac(3,5,j) =   0.0d+00
   113: |||       
   114: |||                      njac(4,1,j) = - c3c4 * tmp2 * u(4,i,j,k)
   115: |||                      njac(4,2,j) =   0.0d+00
   116: |||                      njac(4,3,j) =   0.0d+00
   117: |||                      njac(4,4,j) =   c3c4 * tmp1
   118: |||                      njac(4,5,j) =   0.0d+00
   119: |||       
   120: |||     F                njac(5,1,j) = - (  c3c4
   121: |||            >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   122: |||            >              - ( con43 * c3c4
   123: |||            >              - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   124: |||            >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   125: |||            >              - c1345 * tmp2 * u(5,i,j,k)
   126: |||       
   127: |||                      njac(5,2,j) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
   128: |||                      njac(5,3,j) = ( con43 * c3c4
   129: |||            >              - c1345 ) * tmp2 * u(3,i,j,k)
   130: |||                      njac(5,4,j) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
   131: |||                      njac(5,5,j) = ( c1345 ) * tmp1
   132: |||       
   133: ||V----               enddo
   134: ||        
   135: ||        c---------------------------------------------------------------------
   136: ||        c     now joacobians set, so form left hand side in y direction
   137: ||        c---------------------------------------------------------------------
   138: ||                    call lhsinit(lhs, jsize)
   139: ||V---->              do j = 1, jsize-1
   140: |||       
   141: |||                      tmp1 = dt * ty1
   142: |||                      tmp2 = dt * ty2
   143: |||       
   144: |||     F                lhs(1,1,aa,j) = - tmp2 * fjac(1,1,j-1)
   145: |||            >              - tmp1 * njac(1,1,j-1)
   146: |||            >              - tmp1 * dy1
   147: |||     F                lhs(1,2,aa,j) = - tmp2 * fjac(1,2,j-1)
   148: |||            >              - tmp1 * njac(1,2,j-1)
   149: |||     F                lhs(1,3,aa,j) = - tmp2 * fjac(1,3,j-1)
   150: |||            >              - tmp1 * njac(1,3,j-1)
   151: |||     F                lhs(1,4,aa,j) = - tmp2 * fjac(1,4,j-1)
   152: |||            >              - tmp1 * njac(1,4,j-1)
   153: |||     F                lhs(1,5,aa,j) = - tmp2 * fjac(1,5,j-1)
   154: |||            >              - tmp1 * njac(1,5,j-1)
   155: |||       
   156: |||     F                lhs(2,1,aa,j) = - tmp2 * fjac(2,1,j-1)
   157: |||            >              - tmp1 * njac(2,1,j-1)
   158: |||     F                lhs(2,2,aa,j) = - tmp2 * fjac(2,2,j-1)
   159: |||            >              - tmp1 * njac(2,2,j-1)
   160: |||            >              - tmp1 * dy2
   161: |||     F                lhs(2,3,aa,j) = - tmp2 * fjac(2,3,j-1)
   162: |||            >              - tmp1 * njac(2,3,j-1)
   163: |||     F                lhs(2,4,aa,j) = - tmp2 * fjac(2,4,j-1)
   164: |||            >              - tmp1 * njac(2,4,j-1)
   165: |||     F                lhs(2,5,aa,j) = - tmp2 * fjac(2,5,j-1)
   166: |||            >              - tmp1 * njac(2,5,j-1)
   167: |||       
   168: |||     F                lhs(3,1,aa,j) = - tmp2 * fjac(3,1,j-1)
   169: |||            >              - tmp1 * njac(3,1,j-1)
   170: |||     F                lhs(3,2,aa,j) = - tmp2 * fjac(3,2,j-1)
   171: |||            >              - tmp1 * njac(3,2,j-1)
   172: |||     F                lhs(3,3,aa,j) = - tmp2 * fjac(3,3,j-1)
   173: |||            >              - tmp1 * njac(3,3,j-1)
   174: |||            >              - tmp1 * dy3
   175: |||     F                lhs(3,4,aa,j) = - tmp2 * fjac(3,4,j-1)
   176: |||            >              - tmp1 * njac(3,4,j-1)
   177: |||     F                lhs(3,5,aa,j) = - tmp2 * fjac(3,5,j-1)
   178: |||            >              - tmp1 * njac(3,5,j-1)
   179: |||       
   180: |||     F                lhs(4,1,aa,j) = - tmp2 * fjac(4,1,j-1)
   181: |||            >              - tmp1 * njac(4,1,j-1)
   182: |||     F                lhs(4,2,aa,j) = - tmp2 * fjac(4,2,j-1)
   183: |||            >              - tmp1 * njac(4,2,j-1)
   184: |||     F                lhs(4,3,aa,j) = - tmp2 * fjac(4,3,j-1)
   185: |||            >              - tmp1 * njac(4,3,j-1)
   186: |||     F                lhs(4,4,aa,j) = - tmp2 * fjac(4,4,j-1)
   187: |||            >              - tmp1 * njac(4,4,j-1)
   188: |||            >              - tmp1 * dy4
   189: |||     F                lhs(4,5,aa,j) = - tmp2 * fjac(4,5,j-1)
   190: |||            >              - tmp1 * njac(4,5,j-1)
   191: |||       
   192: |||     F                lhs(5,1,aa,j) = - tmp2 * fjac(5,1,j-1)
   193: |||            >              - tmp1 * njac(5,1,j-1)
   194: |||     F                lhs(5,2,aa,j) = - tmp2 * fjac(5,2,j-1)
   195: |||            >              - tmp1 * njac(5,2,j-1)
   196: |||     F                lhs(5,3,aa,j) = - tmp2 * fjac(5,3,j-1)
   197: |||            >              - tmp1 * njac(5,3,j-1)
   198: |||     F                lhs(5,4,aa,j) = - tmp2 * fjac(5,4,j-1)
   199: |||            >              - tmp1 * njac(5,4,j-1)
   200: |||     F                lhs(5,5,aa,j) = - tmp2 * fjac(5,5,j-1)
   201: |||            >              - tmp1 * njac(5,5,j-1)
   202: |||            >              - tmp1 * dy5
   203: |||       
   204: |||     F                lhs(1,1,bb,j) = 1.0d+00
   205: |||            >              + tmp1 * 2.0d+00 * njac(1,1,j)
   206: |||            >              + tmp1 * 2.0d+00 * dy1
   207: |||                      lhs(1,2,bb,j) = tmp1 * 2.0d+00 * njac(1,2,j)
   208: |||                      lhs(1,3,bb,j) = tmp1 * 2.0d+00 * njac(1,3,j)
   209: |||                      lhs(1,4,bb,j) = tmp1 * 2.0d+00 * njac(1,4,j)
   210: |||                      lhs(1,5,bb,j) = tmp1 * 2.0d+00 * njac(1,5,j)
   211: |||       
   212: |||                      lhs(2,1,bb,j) = tmp1 * 2.0d+00 * njac(2,1,j)
   213: |||     F                lhs(2,2,bb,j) = 1.0d+00
   214: |||            >              + tmp1 * 2.0d+00 * njac(2,2,j)
   215: |||            >              + tmp1 * 2.0d+00 * dy2
   216: |||                      lhs(2,3,bb,j) = tmp1 * 2.0d+00 * njac(2,3,j)
   217: |||                      lhs(2,4,bb,j) = tmp1 * 2.0d+00 * njac(2,4,j)
   218: |||                      lhs(2,5,bb,j) = tmp1 * 2.0d+00 * njac(2,5,j)
   219: |||       
   220: |||                      lhs(3,1,bb,j) = tmp1 * 2.0d+00 * njac(3,1,j)
   221: |||                      lhs(3,2,bb,j) = tmp1 * 2.0d+00 * njac(3,2,j)
   222: |||     F                lhs(3,3,bb,j) = 1.0d+00
   223: |||            >              + tmp1 * 2.0d+00 * njac(3,3,j)
   224: |||            >              + tmp1 * 2.0d+00 * dy3
   225: |||                      lhs(3,4,bb,j) = tmp1 * 2.0d+00 * njac(3,4,j)
   226: |||                      lhs(3,5,bb,j) = tmp1 * 2.0d+00 * njac(3,5,j)
   227: |||       
   228: |||                      lhs(4,1,bb,j) = tmp1 * 2.0d+00 * njac(4,1,j)
   229: |||                      lhs(4,2,bb,j) = tmp1 * 2.0d+00 * njac(4,2,j)
   230: |||                      lhs(4,3,bb,j) = tmp1 * 2.0d+00 * njac(4,3,j)
   231: |||     F                lhs(4,4,bb,j) = 1.0d+00
   232: |||            >              + tmp1 * 2.0d+00 * njac(4,4,j)
   233: |||            >              + tmp1 * 2.0d+00 * dy4
   234: |||                      lhs(4,5,bb,j) = tmp1 * 2.0d+00 * njac(4,5,j)
   235: |||       
   236: |||                      lhs(5,1,bb,j) = tmp1 * 2.0d+00 * njac(5,1,j)
   237: |||                      lhs(5,2,bb,j) = tmp1 * 2.0d+00 * njac(5,2,j)
   238: |||                      lhs(5,3,bb,j) = tmp1 * 2.0d+00 * njac(5,3,j)
   239: |||                      lhs(5,4,bb,j) = tmp1 * 2.0d+00 * njac(5,4,j)
   240: |||     F                lhs(5,5,bb,j) = 1.0d+00
   241: |||            >              + tmp1 * 2.0d+00 * njac(5,5,j)
   242: |||            >              + tmp1 * 2.0d+00 * dy5
   243: |||       
   244: |||     F                lhs(1,1,cc,j) =  tmp2 * fjac(1,1,j+1)
   245: |||            >              - tmp1 * njac(1,1,j+1)
   246: |||            >              - tmp1 * dy1
   247: |||     F                lhs(1,2,cc,j) =  tmp2 * fjac(1,2,j+1)
   248: |||            >              - tmp1 * njac(1,2,j+1)
   249: |||     F                lhs(1,3,cc,j) =  tmp2 * fjac(1,3,j+1)
   250: |||            >              - tmp1 * njac(1,3,j+1)
   251: |||     F                lhs(1,4,cc,j) =  tmp2 * fjac(1,4,j+1)
   252: |||            >              - tmp1 * njac(1,4,j+1)
   253: |||     F                lhs(1,5,cc,j) =  tmp2 * fjac(1,5,j+1)
   254: |||            >              - tmp1 * njac(1,5,j+1)
   255: |||       
   256: |||     F                lhs(2,1,cc,j) =  tmp2 * fjac(2,1,j+1)
   257: |||            >              - tmp1 * njac(2,1,j+1)
   258: |||     F                lhs(2,2,cc,j) =  tmp2 * fjac(2,2,j+1)
   259: |||            >              - tmp1 * njac(2,2,j+1)
   260: |||            >              - tmp1 * dy2
   261: |||     F                lhs(2,3,cc,j) =  tmp2 * fjac(2,3,j+1)
   262: |||            >              - tmp1 * njac(2,3,j+1)
   263: |||     F                lhs(2,4,cc,j) =  tmp2 * fjac(2,4,j+1)
   264: |||            >              - tmp1 * njac(2,4,j+1)
   265: |||     F                lhs(2,5,cc,j) =  tmp2 * fjac(2,5,j+1)
   266: |||            >              - tmp1 * njac(2,5,j+1)
   267: |||       
   268: |||     F                lhs(3,1,cc,j) =  tmp2 * fjac(3,1,j+1)
   269: |||            >              - tmp1 * njac(3,1,j+1)
   270: |||     F                lhs(3,2,cc,j) =  tmp2 * fjac(3,2,j+1)
   271: |||            >              - tmp1 * njac(3,2,j+1)
   272: |||     F                lhs(3,3,cc,j) =  tmp2 * fjac(3,3,j+1)
   273: |||            >              - tmp1 * njac(3,3,j+1)
   274: |||            >              - tmp1 * dy3
   275: |||     F                lhs(3,4,cc,j) =  tmp2 * fjac(3,4,j+1)
   276: |||            >              - tmp1 * njac(3,4,j+1)
   277: |||     F                lhs(3,5,cc,j) =  tmp2 * fjac(3,5,j+1)
   278: |||            >              - tmp1 * njac(3,5,j+1)
   279: |||       
   280: |||     F                lhs(4,1,cc,j) =  tmp2 * fjac(4,1,j+1)
   281: |||            >              - tmp1 * njac(4,1,j+1)
   282: |||     F                lhs(4,2,cc,j) =  tmp2 * fjac(4,2,j+1)
   283: |||            >              - tmp1 * njac(4,2,j+1)
   284: |||     F                lhs(4,3,cc,j) =  tmp2 * fjac(4,3,j+1)
   285: |||            >              - tmp1 * njac(4,3,j+1)
   286: |||     F                lhs(4,4,cc,j) =  tmp2 * fjac(4,4,j+1)
   287: |||            >              - tmp1 * njac(4,4,j+1)
   288: |||            >              - tmp1 * dy4
   289: |||     F                lhs(4,5,cc,j) =  tmp2 * fjac(4,5,j+1)
   290: |||            >              - tmp1 * njac(4,5,j+1)
   291: |||       
   292: |||     F                lhs(5,1,cc,j) =  tmp2 * fjac(5,1,j+1)
   293: |||            >              - tmp1 * njac(5,1,j+1)
   294: |||     F                lhs(5,2,cc,j) =  tmp2 * fjac(5,2,j+1)
   295: |||            >              - tmp1 * njac(5,2,j+1)
   296: |||     F                lhs(5,3,cc,j) =  tmp2 * fjac(5,3,j+1)
   297: |||            >              - tmp1 * njac(5,3,j+1)
   298: |||     F                lhs(5,4,cc,j) =  tmp2 * fjac(5,4,j+1)
   299: |||            >              - tmp1 * njac(5,4,j+1)
   300: |||     F                lhs(5,5,cc,j) =  tmp2 * fjac(5,5,j+1)
   301: |||            >              - tmp1 * njac(5,5,j+1)
   302: |||            >              - tmp1 * dy5
   303: |||       
   304: ||V----               enddo
   305: ||        
   306: ||        c---------------------------------------------------------------------
   307: ||        c---------------------------------------------------------------------
   308: ||        
   309: ||        c---------------------------------------------------------------------
   310: ||        c     performs guaussian elimination on this cell.
   311: ||        c
   312: ||        c     assumes that unpacking routines for non-first cells
   313: ||        c     preload C' and rhs' from previous cell.
   314: ||        c
   315: ||        c     assumed send happens outside this routine, but that
   316: ||        c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
   317: ||        c---------------------------------------------------------------------
   318: ||        
   319: ||                    if (timeron) call timer_start(t_solsub)
   320: ||        c---------------------------------------------------------------------
   321: ||        c     multiply c(i,0,k) by b_inverse and copy back to c
   322: ||        c     multiply rhs(0) by b_inverse(0) and copy to rhs
   323: ||        c---------------------------------------------------------------------
   324: ||                    call binvcrhs( lhs(1,1,bb,0),
   325: ||             >                        lhs(1,1,cc,0),
   326: ||             >                        rhs(1,i,0,k) )
   327: ||        
   328: ||        c---------------------------------------------------------------------
   329: ||        c     begin inner most do loop
   330: ||        c     do all the elements of the cell unless last
   331: ||        c---------------------------------------------------------------------
   332: ||+---->              do j=1,jsize-1
   333: |||       
   334: |||       c---------------------------------------------------------------------
   335: |||       c     subtract A*lhs_vector(j-1) from lhs_vector(j)
   336: |||       c
   337: |||       c     rhs(j) = rhs(j) - A*rhs(j-1)
   338: |||       c---------------------------------------------------------------------
   339: |||                      call matvec_sub(lhs(1,1,aa,j),
   340: |||            >                         rhs(1,i,j-1,k),rhs(1,i,j,k))
   341: |||       
   342: |||       c---------------------------------------------------------------------
   343: |||       c     B(j) = B(j) - C(j-1)*A(j)
   344: |||       c---------------------------------------------------------------------
   345: |||                      call matmul_sub(lhs(1,1,aa,j),
   346: |||            >                         lhs(1,1,cc,j-1),
   347: |||            >                         lhs(1,1,bb,j))
   348: |||       
   349: |||       c---------------------------------------------------------------------
   350: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   351: |||       c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
   352: |||       c---------------------------------------------------------------------
   353: |||                      call binvcrhs( lhs(1,1,bb,j),
   354: |||            >                        lhs(1,1,cc,j),
   355: |||            >                        rhs(1,i,j,k) )
   356: |||       
   357: ||+----               enddo
   358: ||        
   359: ||        
   360: ||        c---------------------------------------------------------------------
   361: ||        c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
   362: ||        c---------------------------------------------------------------------
   363: ||                    call matvec_sub(lhs(1,1,aa,jsize),
   364: ||             >                         rhs(1,i,jsize-1,k),rhs(1,i,jsize,k))
   365: ||        
   366: ||        c---------------------------------------------------------------------
   367: ||        c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
   368: ||        c     call matmul_sub(aa,i,jsize,k,c,
   369: ||        c     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
   370: ||        c---------------------------------------------------------------------
   371: ||                    call matmul_sub(lhs(1,1,aa,jsize),
   372: ||             >                         lhs(1,1,cc,jsize-1),
   373: ||             >                         lhs(1,1,bb,jsize))
   374: ||        
   375: ||        c---------------------------------------------------------------------
   376: ||        c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
   377: ||        c---------------------------------------------------------------------
   378: ||                    call binvrhs( lhs(1,1,bb,jsize),
   379: ||             >                       rhs(1,i,jsize,k) )
   380: ||                    if (timeron) call timer_stop(t_solsub)
   381: ||        
   382: ||        
   383: ||        c---------------------------------------------------------------------
   384: ||        c     back solve: if last cell, then generate U(jsize)=rhs(jsize)
   385: ||        c     else assume U(jsize) is loaded in un pack backsub_info
   386: ||        c     so just use it
   387: ||        c     after call u(jstart) will be sent to next cell
   388: ||        c---------------------------------------------------------------------
   389: ||         
   390: ||+---->              do j=jsize-1,0,-1
   391: |||*--->                 do m=1,BLOCK_SIZE
   392: ||||*-->                    do n=1,BLOCK_SIZE
   393: |||||                          rhs(m,i,j,k) = rhs(m,i,j,k)
   394: |||||          >                    - lhs(m,n,cc,j)*rhs(n,i,j+1,k)
   395: ||||*--                     enddo
   396: |||*---                  enddo
   397: ||+----               enddo
   398: ||        
   399: |+-----            enddo
   400: P------         enddo
   401:                 if (timeron) call timer_stop(t_ysolve)
   402:           
   403:                 return
   404:                 end


