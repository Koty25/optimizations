NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: z_solve.f

PROCEDURE NAME: Z_SOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    45: par(1801): Parallel routine generated.: Z_SOLVE$1
    47: par(1803): Parallelized by "do".
    47: opt(2495): Loop nest collapsed by "collapse" clause.
    49: vec( 124): Iteration count is assumed. Iteration count.: 163
    49: vec( 101): Vectorized loop.
    75: vec( 128): Fused multiply-add operation applied.
    83: vec( 128): Fused multiply-add operation applied.
    90: vec( 128): Fused multiply-add operation applied.
   120: vec( 128): Fused multiply-add operation applied.
   139: vec( 124): Iteration count is assumed. Iteration count.: 163
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
   337: vec( 110): Vectorization obstructive procedure reference.: MATMUL_SUB
   337: vec( 110): Vectorization obstructive procedure reference.: BINVCRHS
   337: vec( 110): Vectorization obstructive procedure reference.: MATVEC_SUB
   337: vec( 103): Unvectorized loop.
   344: opt(1025): Reference to this procedure inhibits optimization.: MATVEC_SUB
   402: vec( 103): Unvectorized loop.
   402: vec( 113): Overhead of loop division is too large.
   403: opt(1418): Constant-length loop is expanded.
   404: opt(1418): Constant-length loop is expanded.
   405: opt(1037): Feedback of array elements.: RHS
   405: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: RHS


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: z_solve.f

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
    25:                 double precision tmp1, tmp2, tmp3
    26:            
    27:           c---------------------------------------------------------------------
    28:           c---------------------------------------------------------------------
    29:           
    30:                 if (timeron) call timer_start(t_zsolve)
    31:           
    32:           c---------------------------------------------------------------------
    33:           c---------------------------------------------------------------------
    34:           
    35:           c---------------------------------------------------------------------
    36:           c     This function computes the left hand side for the three z-factors
    37:           c---------------------------------------------------------------------
    38:           
    39:                 ksize = grid_points(3)-1
    40:           
    41:           c---------------------------------------------------------------------
    42:           c     Compute the indices for storing the block-diagonal matrix;
    43:           c     determine c (labeled f) and s jacobians
    44:           c---------------------------------------------------------------------
    45:           !$omp parallel do default(shared) shared(ksize) collapse(2)
    46:           !$omp& private(i,j,k,m,n,tmp1,tmp2,tmp3)
    47: P------>        do j = 1, grid_points(2)-2
    48: |+----->           do i = 1, grid_points(1)-2
    49: ||V---->              do k = 0, ksize
    50: |||       
    51: |||                      tmp1 = 1.0d+00 / u(1,i,j,k)
    52: |||                      tmp2 = tmp1 * tmp1
    53: |||                      tmp3 = tmp1 * tmp2
    54: |||       
    55: |||                      fjac(1,1,k) = 0.0d+00
    56: |||                      fjac(1,2,k) = 0.0d+00
    57: |||                      fjac(1,3,k) = 0.0d+00
    58: |||                      fjac(1,4,k) = 1.0d+00
    59: |||                      fjac(1,5,k) = 0.0d+00
    60: |||       
    61: |||                      fjac(2,1,k) = - ( u(2,i,j,k)*u(4,i,j,k) )
    62: |||            >              * tmp2
    63: |||                      fjac(2,2,k) = u(4,i,j,k) * tmp1
    64: |||                      fjac(2,3,k) = 0.0d+00
    65: |||                      fjac(2,4,k) = u(2,i,j,k) * tmp1
    66: |||                      fjac(2,5,k) = 0.0d+00
    67: |||       
    68: |||                      fjac(3,1,k) = - ( u(3,i,j,k)*u(4,i,j,k) )
    69: |||            >              * tmp2
    70: |||                      fjac(3,2,k) = 0.0d+00
    71: |||                      fjac(3,3,k) = u(4,i,j,k) * tmp1
    72: |||                      fjac(3,4,k) = u(3,i,j,k) * tmp1
    73: |||                      fjac(3,5,k) = 0.0d+00
    74: |||       
    75: |||     F                fjac(4,1,k) = - (u(4,i,j,k)*u(4,i,j,k) * tmp2 )
    76: |||            >              + c2 * qs(i,j,k)
    77: |||                      fjac(4,2,k) = - c2 *  u(2,i,j,k) * tmp1
    78: |||                      fjac(4,3,k) = - c2 *  u(3,i,j,k) * tmp1
    79: |||                      fjac(4,4,k) = ( 2.0d+00 - c2 )
    80: |||            >              *  u(4,i,j,k) * tmp1
    81: |||                      fjac(4,5,k) = c2
    82: |||       
    83: |||     F                fjac(5,1,k) = ( c2 * 2.0d0 * square(i,j,k)
    84: |||            >              - c1 * u(5,i,j,k) )
    85: |||            >              * u(4,i,j,k) * tmp2
    86: |||                      fjac(5,2,k) = - c2 * ( u(2,i,j,k)*u(4,i,j,k) )
    87: |||            >              * tmp2
    88: |||                      fjac(5,3,k) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
    89: |||            >              * tmp2
    90: |||     F                fjac(5,4,k) = c1 * ( u(5,i,j,k) * tmp1 )
    91: |||            >              - c2
    92: |||            >              * ( qs(i,j,k)
    93: |||            >              + u(4,i,j,k)*u(4,i,j,k) * tmp2 )
    94: |||                      fjac(5,5,k) = c1 * u(4,i,j,k) * tmp1
    95: |||       
    96: |||                      njac(1,1,k) = 0.0d+00
    97: |||                      njac(1,2,k) = 0.0d+00
    98: |||                      njac(1,3,k) = 0.0d+00
    99: |||                      njac(1,4,k) = 0.0d+00
   100: |||                      njac(1,5,k) = 0.0d+00
   101: |||       
   102: |||                      njac(2,1,k) = - c3c4 * tmp2 * u(2,i,j,k)
   103: |||                      njac(2,2,k) =   c3c4 * tmp1
   104: |||                      njac(2,3,k) =   0.0d+00
   105: |||                      njac(2,4,k) =   0.0d+00
   106: |||                      njac(2,5,k) =   0.0d+00
   107: |||       
   108: |||                      njac(3,1,k) = - c3c4 * tmp2 * u(3,i,j,k)
   109: |||                      njac(3,2,k) =   0.0d+00
   110: |||                      njac(3,3,k) =   c3c4 * tmp1
   111: |||                      njac(3,4,k) =   0.0d+00
   112: |||                      njac(3,5,k) =   0.0d+00
   113: |||       
   114: |||                      njac(4,1,k) = - con43 * c3c4 * tmp2 * u(4,i,j,k)
   115: |||                      njac(4,2,k) =   0.0d+00
   116: |||                      njac(4,3,k) =   0.0d+00
   117: |||                      njac(4,4,k) =   con43 * c3 * c4 * tmp1
   118: |||                      njac(4,5,k) =   0.0d+00
   119: |||       
   120: |||     F                njac(5,1,k) = - (  c3c4
   121: |||            >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
   122: |||            >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
   123: |||            >              - ( con43 * c3c4
   124: |||            >              - c1345 ) * tmp3 * (u(4,i,j,k)**2)
   125: |||            >              - c1345 * tmp2 * u(5,i,j,k)
   126: |||       
   127: |||                      njac(5,2,k) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
   128: |||                      njac(5,3,k) = (  c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
   129: |||                      njac(5,4,k) = ( con43 * c3c4
   130: |||            >              - c1345 ) * tmp2 * u(4,i,j,k)
   131: |||                      njac(5,5,k) = ( c1345 )* tmp1
   132: |||       
   133: ||V----               enddo
   134: ||        
   135: ||        c---------------------------------------------------------------------
   136: ||        c     now jacobians set, so form left hand side in z direction
   137: ||        c---------------------------------------------------------------------
   138: ||                    call lhsinit(lhs, ksize)
   139: ||V---->              do k = 1, ksize-1
   140: |||       
   141: |||                      tmp1 = dt * tz1
   142: |||                      tmp2 = dt * tz2
   143: |||       
   144: |||     F                lhs(1,1,aa,k) = - tmp2 * fjac(1,1,k-1)
   145: |||            >              - tmp1 * njac(1,1,k-1)
   146: |||            >              - tmp1 * dz1
   147: |||     F                lhs(1,2,aa,k) = - tmp2 * fjac(1,2,k-1)
   148: |||            >              - tmp1 * njac(1,2,k-1)
   149: |||     F                lhs(1,3,aa,k) = - tmp2 * fjac(1,3,k-1)
   150: |||            >              - tmp1 * njac(1,3,k-1)
   151: |||     F                lhs(1,4,aa,k) = - tmp2 * fjac(1,4,k-1)
   152: |||            >              - tmp1 * njac(1,4,k-1)
   153: |||     F                lhs(1,5,aa,k) = - tmp2 * fjac(1,5,k-1)
   154: |||            >              - tmp1 * njac(1,5,k-1)
   155: |||       
   156: |||     F                lhs(2,1,aa,k) = - tmp2 * fjac(2,1,k-1)
   157: |||            >              - tmp1 * njac(2,1,k-1)
   158: |||     F                lhs(2,2,aa,k) = - tmp2 * fjac(2,2,k-1)
   159: |||            >              - tmp1 * njac(2,2,k-1)
   160: |||            >              - tmp1 * dz2
   161: |||     F                lhs(2,3,aa,k) = - tmp2 * fjac(2,3,k-1)
   162: |||            >              - tmp1 * njac(2,3,k-1)
   163: |||     F                lhs(2,4,aa,k) = - tmp2 * fjac(2,4,k-1)
   164: |||            >              - tmp1 * njac(2,4,k-1)
   165: |||     F                lhs(2,5,aa,k) = - tmp2 * fjac(2,5,k-1)
   166: |||            >              - tmp1 * njac(2,5,k-1)
   167: |||       
   168: |||     F                lhs(3,1,aa,k) = - tmp2 * fjac(3,1,k-1)
   169: |||            >              - tmp1 * njac(3,1,k-1)
   170: |||     F                lhs(3,2,aa,k) = - tmp2 * fjac(3,2,k-1)
   171: |||            >              - tmp1 * njac(3,2,k-1)
   172: |||     F                lhs(3,3,aa,k) = - tmp2 * fjac(3,3,k-1)
   173: |||            >              - tmp1 * njac(3,3,k-1)
   174: |||            >              - tmp1 * dz3
   175: |||     F                lhs(3,4,aa,k) = - tmp2 * fjac(3,4,k-1)
   176: |||            >              - tmp1 * njac(3,4,k-1)
   177: |||     F                lhs(3,5,aa,k) = - tmp2 * fjac(3,5,k-1)
   178: |||            >              - tmp1 * njac(3,5,k-1)
   179: |||       
   180: |||     F                lhs(4,1,aa,k) = - tmp2 * fjac(4,1,k-1)
   181: |||            >              - tmp1 * njac(4,1,k-1)
   182: |||     F                lhs(4,2,aa,k) = - tmp2 * fjac(4,2,k-1)
   183: |||            >              - tmp1 * njac(4,2,k-1)
   184: |||     F                lhs(4,3,aa,k) = - tmp2 * fjac(4,3,k-1)
   185: |||            >              - tmp1 * njac(4,3,k-1)
   186: |||     F                lhs(4,4,aa,k) = - tmp2 * fjac(4,4,k-1)
   187: |||            >              - tmp1 * njac(4,4,k-1)
   188: |||            >              - tmp1 * dz4
   189: |||     F                lhs(4,5,aa,k) = - tmp2 * fjac(4,5,k-1)
   190: |||            >              - tmp1 * njac(4,5,k-1)
   191: |||       
   192: |||     F                lhs(5,1,aa,k) = - tmp2 * fjac(5,1,k-1)
   193: |||            >              - tmp1 * njac(5,1,k-1)
   194: |||     F                lhs(5,2,aa,k) = - tmp2 * fjac(5,2,k-1)
   195: |||            >              - tmp1 * njac(5,2,k-1)
   196: |||     F                lhs(5,3,aa,k) = - tmp2 * fjac(5,3,k-1)
   197: |||            >              - tmp1 * njac(5,3,k-1)
   198: |||     F                lhs(5,4,aa,k) = - tmp2 * fjac(5,4,k-1)
   199: |||            >              - tmp1 * njac(5,4,k-1)
   200: |||     F                lhs(5,5,aa,k) = - tmp2 * fjac(5,5,k-1)
   201: |||            >              - tmp1 * njac(5,5,k-1)
   202: |||            >              - tmp1 * dz5
   203: |||       
   204: |||     F                lhs(1,1,bb,k) = 1.0d+00
   205: |||            >              + tmp1 * 2.0d+00 * njac(1,1,k)
   206: |||            >              + tmp1 * 2.0d+00 * dz1
   207: |||                      lhs(1,2,bb,k) = tmp1 * 2.0d+00 * njac(1,2,k)
   208: |||                      lhs(1,3,bb,k) = tmp1 * 2.0d+00 * njac(1,3,k)
   209: |||                      lhs(1,4,bb,k) = tmp1 * 2.0d+00 * njac(1,4,k)
   210: |||                      lhs(1,5,bb,k) = tmp1 * 2.0d+00 * njac(1,5,k)
   211: |||       
   212: |||                      lhs(2,1,bb,k) = tmp1 * 2.0d+00 * njac(2,1,k)
   213: |||     F                lhs(2,2,bb,k) = 1.0d+00
   214: |||            >              + tmp1 * 2.0d+00 * njac(2,2,k)
   215: |||            >              + tmp1 * 2.0d+00 * dz2
   216: |||                      lhs(2,3,bb,k) = tmp1 * 2.0d+00 * njac(2,3,k)
   217: |||                      lhs(2,4,bb,k) = tmp1 * 2.0d+00 * njac(2,4,k)
   218: |||                      lhs(2,5,bb,k) = tmp1 * 2.0d+00 * njac(2,5,k)
   219: |||       
   220: |||                      lhs(3,1,bb,k) = tmp1 * 2.0d+00 * njac(3,1,k)
   221: |||                      lhs(3,2,bb,k) = tmp1 * 2.0d+00 * njac(3,2,k)
   222: |||     F                lhs(3,3,bb,k) = 1.0d+00
   223: |||            >              + tmp1 * 2.0d+00 * njac(3,3,k)
   224: |||            >              + tmp1 * 2.0d+00 * dz3
   225: |||                      lhs(3,4,bb,k) = tmp1 * 2.0d+00 * njac(3,4,k)
   226: |||                      lhs(3,5,bb,k) = tmp1 * 2.0d+00 * njac(3,5,k)
   227: |||       
   228: |||                      lhs(4,1,bb,k) = tmp1 * 2.0d+00 * njac(4,1,k)
   229: |||                      lhs(4,2,bb,k) = tmp1 * 2.0d+00 * njac(4,2,k)
   230: |||                      lhs(4,3,bb,k) = tmp1 * 2.0d+00 * njac(4,3,k)
   231: |||     F                lhs(4,4,bb,k) = 1.0d+00
   232: |||            >              + tmp1 * 2.0d+00 * njac(4,4,k)
   233: |||            >              + tmp1 * 2.0d+00 * dz4
   234: |||                      lhs(4,5,bb,k) = tmp1 * 2.0d+00 * njac(4,5,k)
   235: |||       
   236: |||                      lhs(5,1,bb,k) = tmp1 * 2.0d+00 * njac(5,1,k)
   237: |||                      lhs(5,2,bb,k) = tmp1 * 2.0d+00 * njac(5,2,k)
   238: |||                      lhs(5,3,bb,k) = tmp1 * 2.0d+00 * njac(5,3,k)
   239: |||                      lhs(5,4,bb,k) = tmp1 * 2.0d+00 * njac(5,4,k)
   240: |||     F                lhs(5,5,bb,k) = 1.0d+00
   241: |||            >              + tmp1 * 2.0d+00 * njac(5,5,k)
   242: |||            >              + tmp1 * 2.0d+00 * dz5
   243: |||       
   244: |||     F                lhs(1,1,cc,k) =  tmp2 * fjac(1,1,k+1)
   245: |||            >              - tmp1 * njac(1,1,k+1)
   246: |||            >              - tmp1 * dz1
   247: |||     F                lhs(1,2,cc,k) =  tmp2 * fjac(1,2,k+1)
   248: |||            >              - tmp1 * njac(1,2,k+1)
   249: |||     F                lhs(1,3,cc,k) =  tmp2 * fjac(1,3,k+1)
   250: |||            >              - tmp1 * njac(1,3,k+1)
   251: |||     F                lhs(1,4,cc,k) =  tmp2 * fjac(1,4,k+1)
   252: |||            >              - tmp1 * njac(1,4,k+1)
   253: |||     F                lhs(1,5,cc,k) =  tmp2 * fjac(1,5,k+1)
   254: |||            >              - tmp1 * njac(1,5,k+1)
   255: |||       
   256: |||     F                lhs(2,1,cc,k) =  tmp2 * fjac(2,1,k+1)
   257: |||            >              - tmp1 * njac(2,1,k+1)
   258: |||     F                lhs(2,2,cc,k) =  tmp2 * fjac(2,2,k+1)
   259: |||            >              - tmp1 * njac(2,2,k+1)
   260: |||            >              - tmp1 * dz2
   261: |||     F                lhs(2,3,cc,k) =  tmp2 * fjac(2,3,k+1)
   262: |||            >              - tmp1 * njac(2,3,k+1)
   263: |||     F                lhs(2,4,cc,k) =  tmp2 * fjac(2,4,k+1)
   264: |||            >              - tmp1 * njac(2,4,k+1)
   265: |||     F                lhs(2,5,cc,k) =  tmp2 * fjac(2,5,k+1)
   266: |||            >              - tmp1 * njac(2,5,k+1)
   267: |||       
   268: |||     F                lhs(3,1,cc,k) =  tmp2 * fjac(3,1,k+1)
   269: |||            >              - tmp1 * njac(3,1,k+1)
   270: |||     F                lhs(3,2,cc,k) =  tmp2 * fjac(3,2,k+1)
   271: |||            >              - tmp1 * njac(3,2,k+1)
   272: |||     F                lhs(3,3,cc,k) =  tmp2 * fjac(3,3,k+1)
   273: |||            >              - tmp1 * njac(3,3,k+1)
   274: |||            >              - tmp1 * dz3
   275: |||     F                lhs(3,4,cc,k) =  tmp2 * fjac(3,4,k+1)
   276: |||            >              - tmp1 * njac(3,4,k+1)
   277: |||     F                lhs(3,5,cc,k) =  tmp2 * fjac(3,5,k+1)
   278: |||            >              - tmp1 * njac(3,5,k+1)
   279: |||       
   280: |||     F                lhs(4,1,cc,k) =  tmp2 * fjac(4,1,k+1)
   281: |||            >              - tmp1 * njac(4,1,k+1)
   282: |||     F                lhs(4,2,cc,k) =  tmp2 * fjac(4,2,k+1)
   283: |||            >              - tmp1 * njac(4,2,k+1)
   284: |||     F                lhs(4,3,cc,k) =  tmp2 * fjac(4,3,k+1)
   285: |||            >              - tmp1 * njac(4,3,k+1)
   286: |||     F                lhs(4,4,cc,k) =  tmp2 * fjac(4,4,k+1)
   287: |||            >              - tmp1 * njac(4,4,k+1)
   288: |||            >              - tmp1 * dz4
   289: |||     F                lhs(4,5,cc,k) =  tmp2 * fjac(4,5,k+1)
   290: |||            >              - tmp1 * njac(4,5,k+1)
   291: |||       
   292: |||     F                lhs(5,1,cc,k) =  tmp2 * fjac(5,1,k+1)
   293: |||            >              - tmp1 * njac(5,1,k+1)
   294: |||     F                lhs(5,2,cc,k) =  tmp2 * fjac(5,2,k+1)
   295: |||            >              - tmp1 * njac(5,2,k+1)
   296: |||     F                lhs(5,3,cc,k) =  tmp2 * fjac(5,3,k+1)
   297: |||            >              - tmp1 * njac(5,3,k+1)
   298: |||     F                lhs(5,4,cc,k) =  tmp2 * fjac(5,4,k+1)
   299: |||            >              - tmp1 * njac(5,4,k+1)
   300: |||     F                lhs(5,5,cc,k) =  tmp2 * fjac(5,5,k+1)
   301: |||            >              - tmp1 * njac(5,5,k+1)
   302: |||            >              - tmp1 * dz5
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
   316: ||        c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
   317: ||        c---------------------------------------------------------------------
   318: ||        
   319: ||        c---------------------------------------------------------------------
   320: ||        c     outer most do loops - sweeping in i direction
   321: ||        c---------------------------------------------------------------------
   322: ||        
   323: ||                    if (timeron) call timer_start(t_solsub)
   324: ||        c---------------------------------------------------------------------
   325: ||        c     multiply c(i,j,0) by b_inverse and copy back to c
   326: ||        c     multiply rhs(0) by b_inverse(0) and copy to rhs
   327: ||        c---------------------------------------------------------------------
   328: ||                    call binvcrhs( lhs(1,1,bb,0),
   329: ||             >                        lhs(1,1,cc,0),
   330: ||             >                        rhs(1,i,j,0) )
   331: ||        
   332: ||        
   333: ||        c---------------------------------------------------------------------
   334: ||        c     begin inner most do loop
   335: ||        c     do all the elements of the cell unless last
   336: ||        c---------------------------------------------------------------------
   337: ||+---->              do k=1,ksize-1
   338: |||       
   339: |||       c---------------------------------------------------------------------
   340: |||       c     subtract A*lhs_vector(k-1) from lhs_vector(k)
   341: |||       c
   342: |||       c     rhs(k) = rhs(k) - A*rhs(k-1)
   343: |||       c---------------------------------------------------------------------
   344: |||                      call matvec_sub(lhs(1,1,aa,k),
   345: |||            >                         rhs(1,i,j,k-1),rhs(1,i,j,k))
   346: |||       
   347: |||       c---------------------------------------------------------------------
   348: |||       c     B(k) = B(k) - C(k-1)*A(k)
   349: |||       c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
   350: |||       c---------------------------------------------------------------------
   351: |||                      call matmul_sub(lhs(1,1,aa,k),
   352: |||            >                         lhs(1,1,cc,k-1),
   353: |||            >                         lhs(1,1,bb,k))
   354: |||       
   355: |||       c---------------------------------------------------------------------
   356: |||       c     multiply c(i,j,k) by b_inverse and copy back to c
   357: |||       c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
   358: |||       c---------------------------------------------------------------------
   359: |||                      call binvcrhs( lhs(1,1,bb,k),
   360: |||            >                        lhs(1,1,cc,k),
   361: |||            >                        rhs(1,i,j,k) )
   362: |||       
   363: ||+----               enddo
   364: ||        
   365: ||        c---------------------------------------------------------------------
   366: ||        c     Now finish up special cases for last cell
   367: ||        c---------------------------------------------------------------------
   368: ||        
   369: ||        c---------------------------------------------------------------------
   370: ||        c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
   371: ||        c---------------------------------------------------------------------
   372: ||                    call matvec_sub(lhs(1,1,aa,ksize),
   373: ||             >                         rhs(1,i,j,ksize-1),rhs(1,i,j,ksize))
   374: ||        
   375: ||        c---------------------------------------------------------------------
   376: ||        c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
   377: ||        c     call matmul_sub(aa,i,j,ksize,c,
   378: ||        c     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
   379: ||        c---------------------------------------------------------------------
   380: ||                    call matmul_sub(lhs(1,1,aa,ksize),
   381: ||             >                         lhs(1,1,cc,ksize-1),
   382: ||             >                         lhs(1,1,bb,ksize))
   383: ||        
   384: ||        c---------------------------------------------------------------------
   385: ||        c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
   386: ||        c---------------------------------------------------------------------
   387: ||                    call binvrhs( lhs(1,1,bb,ksize),
   388: ||             >                       rhs(1,i,j,ksize) )
   389: ||                    if (timeron) call timer_stop(t_solsub)
   390: ||        
   391: ||        
   392: ||        c---------------------------------------------------------------------
   393: ||        c---------------------------------------------------------------------
   394: ||        
   395: ||        c---------------------------------------------------------------------
   396: ||        c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
   397: ||        c     else assume U(ksize) is loaded in un pack backsub_info
   398: ||        c     so just use it
   399: ||        c     after call u(kstart) will be sent to next cell
   400: ||        c---------------------------------------------------------------------
   401: ||        
   402: ||+---->              do k=ksize-1,0,-1
   403: |||*--->                 do m=1,BLOCK_SIZE
   404: ||||*-->                    do n=1,BLOCK_SIZE
   405: |||||                          rhs(m,i,j,k) = rhs(m,i,j,k)
   406: |||||          >                    - lhs(m,n,cc,k)*rhs(n,i,j,k+1)
   407: ||||*--                     enddo
   408: |||*---                  enddo
   409: ||+----               enddo
   410: ||        
   411: |+-----            enddo
   412: P------         enddo
   413:                 if (timeron) call timer_stop(t_zsolve)
   414:           
   415:                 return
   416:                 end


