NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME : exact_rhs.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    20: par(1801): Parallel routine generated.: EXACT_RHS$1
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    26: vec( 101): Vectorized loop.
    28: vec( 101): Vectorized loop.
    29: opt(1418): Constant-length loop is expanded.
    40: opt(3008): Reference within a conditional branch moved outside loop.
    40: par(1803): Parallelized by "do".
    40: opt(2495): Loop nest collapsed by "collapse" clause.
    45: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    45: vec( 103): Unvectorized loop.
    48: inl(1212): Source for routine not found.: EXACT_SOLUTION
    48: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    49: opt(1418): Constant-length loop is expanded.
    55: opt(1418): Constant-length loop is expanded.
    67: vec( 101): Vectorized loop.
    67: vec( 124): Iteration count is assumed. Iteration count.: 163
    71: vec( 128): Fused multiply-add operation applied.
    75: vec( 128): Fused multiply-add operation applied.
    81: vec( 128): Fused multiply-add operation applied.
    86: vec( 128): Fused multiply-add operation applied.
    91: vec( 128): Fused multiply-add operation applied.
   105: opt(1418): Constant-length loop is expanded.
   115: vec( 101): Vectorized loop.
   116: opt(1418): Constant-length loop is expanded.
   117: vec( 128): Fused multiply-add operation applied.
   123: opt(1418): Constant-length loop is expanded.
   141: opt(3008): Reference within a conditional branch moved outside loop.
   141: par(1803): Parallelized by "do".
   141: opt(2495): Loop nest collapsed by "collapse" clause.
   146: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   146: vec( 103): Unvectorized loop.
   149: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   150: opt(1418): Constant-length loop is expanded.
   156: opt(1418): Constant-length loop is expanded.
   167: vec( 101): Vectorized loop.
   167: vec( 124): Iteration count is assumed. Iteration count.: 163
   171: vec( 128): Fused multiply-add operation applied.
   175: vec( 128): Fused multiply-add operation applied.
   180: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   204: opt(1418): Constant-length loop is expanded.
   214: vec( 101): Vectorized loop.
   215: opt(1418): Constant-length loop is expanded.
   216: vec( 128): Fused multiply-add operation applied.
   222: opt(1418): Constant-length loop is expanded.
   241: opt(3008): Reference within a conditional branch moved outside loop.
   241: par(1803): Parallelized by "do".
   241: opt(2495): Loop nest collapsed by "collapse" clause.
   246: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   246: vec( 103): Unvectorized loop.
   249: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   250: opt(1418): Constant-length loop is expanded.
   256: opt(1418): Constant-length loop is expanded.
   267: vec( 101): Vectorized loop.
   267: vec( 124): Iteration count is assumed. Iteration count.: 163
   271: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   280: vec( 128): Fused multiply-add operation applied.
   285: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   304: opt(1418): Constant-length loop is expanded.
   314: vec( 101): Vectorized loop.
   315: opt(1418): Constant-length loop is expanded.
   316: vec( 128): Fused multiply-add operation applied.
   322: opt(1418): Constant-length loop is expanded.
   340: par(1803): Parallelized by "do".
   340: opt(2495): Loop nest collapsed by "collapse" clause.
   340: vec( 101): Vectorized loop.
   342: vec( 101): Vectorized loop.
   343: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
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

     5:                 subroutine exact_rhs
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c     compute the right hand side based on exact solution
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use bt_data
    15:                 implicit none
    16:           
    17:                 double precision dtemp(5), xi, eta, zeta, dtpp
    18:                 integer m, i, j, k, ip1, im1, jp1, jm1, km1, kp1
    19:           
    20:           !$omp parallel default(shared) private(i,j,k,m,zeta,eta,xi,
    21:           !$omp&  dtpp,im1,ip1,jm1,jp1,km1,kp1,dtemp)
    22:           c---------------------------------------------------------------------
    23:           c     initialize
    24:           c---------------------------------------------------------------------
    25:           !$omp do schedule(static) collapse(2)
    26: Y------>        do k= 0, grid_points(3)-1
    27: |+----->           do j = 0, grid_points(2)-1
    28: ||V---->              do i = 0, grid_points(1)-1
    29: |||*--->                 do m = 1, 5
    30: ||||                        forcing(m,i,j,k) = 0.0d0
    31: |||*---                  enddo
    32: ||V----               enddo
    33: |+-----            enddo
    34: Y------         enddo
    35:           
    36:           c---------------------------------------------------------------------
    37:           c     xi-direction flux differences
    38:           c---------------------------------------------------------------------
    39:           !$omp do schedule(static) collapse(2)
    40: P------>        do k = 1, grid_points(3)-2
    41: |+----->           do j = 1, grid_points(2)-2
    42: ||                    zeta = dble(k) * dnzm1
    43: ||                    eta = dble(j) * dnym1
    44: ||        
    45: ||+---->              do i=0, grid_points(1)-1
    46: |||                      xi = dble(i) * dnxm1
    47: |||       
    48: |||                      call exact_solution(xi, eta, zeta, dtemp)
    49: |||*--->                 do m = 1, 5
    50: ||||                        ue(i,m) = dtemp(m)
    51: |||*---                  enddo
    52: |||       
    53: |||                      dtpp = 1.0d0 / dtemp(1)
    54: |||       
    55: |||*--->                 do m = 2, 5
    56: ||||                        buf(i,m) = dtpp * dtemp(m)
    57: |||*---                  enddo
    58: |||       
    59: |||                      cuf(i)   = buf(i,2) * buf(i,2)
    60: |||                      buf(i,1) = cuf(i) + buf(i,3) * buf(i,3) +
    61: |||            >                 buf(i,4) * buf(i,4)
    62: |||                      q(i) = 0.5d0*(buf(i,2)*ue(i,2) + buf(i,3)*ue(i,3) +
    63: |||            >                 buf(i,4)*ue(i,4))
    64: |||       
    65: ||+----               enddo
    66: ||         
    67: ||V---->              do i = 1, grid_points(1)-2
    68: |||                      im1 = i-1
    69: |||                      ip1 = i+1
    70: |||       
    71: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
    72: |||            >                 tx2*( ue(ip1,2)-ue(im1,2) )+
    73: |||            >                 dx1tx1*(ue(ip1,1)-2.0d0*ue(i,1)+ue(im1,1))
    74: |||       
    75: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - tx2 * (
    76: |||            >                 (ue(ip1,2)*buf(ip1,2)+c2*(ue(ip1,5)-q(ip1)))-
    77: |||            >                 (ue(im1,2)*buf(im1,2)+c2*(ue(im1,5)-q(im1))))+
    78: |||            >                 xxcon1*(buf(ip1,2)-2.0d0*buf(i,2)+buf(im1,2))+
    79: |||            >                 dx2tx1*( ue(ip1,2)-2.0d0* ue(i,2)+ue(im1,2))
    80: |||       
    81: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - tx2 * (
    82: |||            >                 ue(ip1,3)*buf(ip1,2)-ue(im1,3)*buf(im1,2))+
    83: |||            >                 xxcon2*(buf(ip1,3)-2.0d0*buf(i,3)+buf(im1,3))+
    84: |||            >                 dx3tx1*( ue(ip1,3)-2.0d0*ue(i,3) +ue(im1,3))
    85: |||        
    86: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - tx2*(
    87: |||            >                 ue(ip1,4)*buf(ip1,2)-ue(im1,4)*buf(im1,2))+
    88: |||            >                 xxcon2*(buf(ip1,4)-2.0d0*buf(i,4)+buf(im1,4))+
    89: |||            >                 dx4tx1*( ue(ip1,4)-2.0d0* ue(i,4)+ ue(im1,4))
    90: |||       
    91: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - tx2*(
    92: |||            >                 buf(ip1,2)*(c1*ue(ip1,5)-c2*q(ip1))-
    93: |||            >                 buf(im1,2)*(c1*ue(im1,5)-c2*q(im1)))+
    94: |||            >                 0.5d0*xxcon3*(buf(ip1,1)-2.0d0*buf(i,1)+
    95: |||            >                 buf(im1,1))+
    96: |||            >                 xxcon4*(cuf(ip1)-2.0d0*cuf(i)+cuf(im1))+
    97: |||            >                 xxcon5*(buf(ip1,5)-2.0d0*buf(i,5)+buf(im1,5))+
    98: |||            >                 dx5tx1*( ue(ip1,5)-2.0d0* ue(i,5)+ ue(im1,5))
    99: ||V----               enddo
   100: ||        
   101: ||        c---------------------------------------------------------------------
   102: ||        c     Fourth-order dissipation
   103: ||        c---------------------------------------------------------------------
   104: ||        
   105: ||*---->              do m = 1, 5
   106: |||                      i = 1
   107: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   108: |||            >                    (5.0d0*ue(i,m) - 4.0d0*ue(i+1,m) +ue(i+2,m))
   109: |||                      i = 2
   110: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   111: |||            >                    (-4.0d0*ue(i-1,m) + 6.0d0*ue(i,m) -
   112: |||            >                    4.0d0*ue(i+1,m) +       ue(i+2,m))
   113: ||*----               enddo
   114: ||        
   115: ||V---->              do i = 3, grid_points(1)-4
   116: |||*--->                 do m = 1, 5
   117: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   118: ||||           >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   119: ||||           >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m) + ue(i+2,m))
   120: |||*---                  enddo
   121: ||V----               enddo
   122: ||        
   123: ||*---->              do m = 1, 5
   124: |||                      i = grid_points(1)-3
   125: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   126: |||            >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   127: |||            >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m))
   128: |||                      i = grid_points(1)-2
   129: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   130: |||            >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) + 5.0d0*ue(i,m))
   131: ||*----               enddo
   132: ||        
   133: |+-----            enddo
   134: P------         enddo
   135:           !$omp end do
   136:           
   137:           c---------------------------------------------------------------------
   138:           c     eta-direction flux differences
   139:           c---------------------------------------------------------------------
   140:           !$omp do schedule(static) collapse(2)
   141: P------>        do k = 1, grid_points(3)-2
   142: |+----->           do i=1, grid_points(1)-2
   143: ||                    zeta = dble(k) * dnzm1
   144: ||                    xi = dble(i) * dnxm1
   145: ||        
   146: ||+---->              do j=0, grid_points(2)-1
   147: |||                      eta = dble(j) * dnym1
   148: |||       
   149: |||                      call exact_solution(xi, eta, zeta, dtemp)
   150: |||*--->                 do m = 1, 5
   151: ||||                        ue(j,m) = dtemp(m)
   152: |||*---                  enddo
   153: |||        
   154: |||                      dtpp = 1.0d0/dtemp(1)
   155: |||       
   156: |||*--->                 do m = 2, 5
   157: ||||                        buf(j,m) = dtpp * dtemp(m)
   158: |||*---                  enddo
   159: |||       
   160: |||                      cuf(j)   = buf(j,3) * buf(j,3)
   161: |||                      buf(j,1) = cuf(j) + buf(j,2) * buf(j,2) +
   162: |||            >                 buf(j,4) * buf(j,4)
   163: |||                      q(j) = 0.5d0*(buf(j,2)*ue(j,2) + buf(j,3)*ue(j,3) +
   164: |||            >                 buf(j,4)*ue(j,4))
   165: ||+----               enddo
   166: ||        
   167: ||V---->              do j = 1, grid_points(2)-2
   168: |||                      jm1 = j-1
   169: |||                      jp1 = j+1
   170: |||        
   171: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
   172: |||            >                 ty2*( ue(jp1,3)-ue(jm1,3) )+
   173: |||            >                 dy1ty1*(ue(jp1,1)-2.0d0*ue(j,1)+ue(jm1,1))
   174: |||       
   175: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - ty2*(
   176: |||            >                 ue(jp1,2)*buf(jp1,3)-ue(jm1,2)*buf(jm1,3))+
   177: |||            >                 yycon2*(buf(jp1,2)-2.0d0*buf(j,2)+buf(jm1,2))+
   178: |||            >                 dy2ty1*( ue(jp1,2)-2.0* ue(j,2)+ ue(jm1,2))
   179: |||       
   180: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - ty2*(
   181: |||            >                 (ue(jp1,3)*buf(jp1,3)+c2*(ue(jp1,5)-q(jp1)))-
   182: |||            >                 (ue(jm1,3)*buf(jm1,3)+c2*(ue(jm1,5)-q(jm1))))+
   183: |||            >                 yycon1*(buf(jp1,3)-2.0d0*buf(j,3)+buf(jm1,3))+
   184: |||            >                 dy3ty1*( ue(jp1,3)-2.0d0*ue(j,3) +ue(jm1,3))
   185: |||       
   186: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - ty2*(
   187: |||            >                 ue(jp1,4)*buf(jp1,3)-ue(jm1,4)*buf(jm1,3))+
   188: |||            >                 yycon2*(buf(jp1,4)-2.0d0*buf(j,4)+buf(jm1,4))+
   189: |||            >                 dy4ty1*( ue(jp1,4)-2.0d0*ue(j,4)+ ue(jm1,4))
   190: |||       
   191: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - ty2*(
   192: |||            >                 buf(jp1,3)*(c1*ue(jp1,5)-c2*q(jp1))-
   193: |||            >                 buf(jm1,3)*(c1*ue(jm1,5)-c2*q(jm1)))+
   194: |||            >                 0.5d0*yycon3*(buf(jp1,1)-2.0d0*buf(j,1)+
   195: |||            >                 buf(jm1,1))+
   196: |||            >                 yycon4*(cuf(jp1)-2.0d0*cuf(j)+cuf(jm1))+
   197: |||            >                 yycon5*(buf(jp1,5)-2.0d0*buf(j,5)+buf(jm1,5))+
   198: |||            >                 dy5ty1*(ue(jp1,5)-2.0d0*ue(j,5)+ue(jm1,5))
   199: ||V----               enddo
   200: ||        
   201: ||        c---------------------------------------------------------------------
   202: ||        c     Fourth-order dissipation
   203: ||        c---------------------------------------------------------------------
   204: ||*---->              do m = 1, 5
   205: |||                      j = 1
   206: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   207: |||            >                    (5.0d0*ue(j,m) - 4.0d0*ue(j+1,m) +ue(j+2,m))
   208: |||                      j = 2
   209: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   210: |||            >                    (-4.0d0*ue(j-1,m) + 6.0d0*ue(j,m) -
   211: |||            >                    4.0d0*ue(j+1,m) +       ue(j+2,m))
   212: ||*----               enddo
   213: ||        
   214: ||V---->              do j = 3, grid_points(2)-4
   215: |||*--->                 do m = 1, 5
   216: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   217: ||||           >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   218: ||||           >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m) + ue(j+2,m))
   219: |||*---                  enddo
   220: ||V----               enddo
   221: ||        
   222: ||*---->              do m = 1, 5
   223: |||                      j = grid_points(2)-3
   224: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   225: |||            >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   226: |||            >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m))
   227: |||                      j = grid_points(2)-2
   228: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   229: |||            >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) + 5.0d0*ue(j,m))
   230: |||       
   231: ||*----               enddo
   232: ||        
   233: |+-----            enddo
   234: P------         enddo
   235:           !$omp end do
   236:           
   237:           c---------------------------------------------------------------------
   238:           c     zeta-direction flux differences
   239:           c---------------------------------------------------------------------
   240:           !$omp do schedule(static) collapse(2)
   241: P------>        do j=1, grid_points(2)-2
   242: |+----->           do i = 1, grid_points(1)-2
   243: ||                    eta = dble(j) * dnym1
   244: ||                    xi = dble(i) * dnxm1
   245: ||        
   246: ||+---->              do k=0, grid_points(3)-1
   247: |||                      zeta = dble(k) * dnzm1
   248: |||       
   249: |||                      call exact_solution(xi, eta, zeta, dtemp)
   250: |||*--->                 do m = 1, 5
   251: ||||                        ue(k,m) = dtemp(m)
   252: |||*---                  enddo
   253: |||       
   254: |||                      dtpp = 1.0d0/dtemp(1)
   255: |||       
   256: |||*--->                 do m = 2, 5
   257: ||||                        buf(k,m) = dtpp * dtemp(m)
   258: |||*---                  enddo
   259: |||       
   260: |||                      cuf(k)   = buf(k,4) * buf(k,4)
   261: |||                      buf(k,1) = cuf(k) + buf(k,2) * buf(k,2) +
   262: |||            >                 buf(k,3) * buf(k,3)
   263: |||                      q(k) = 0.5d0*(buf(k,2)*ue(k,2) + buf(k,3)*ue(k,3) +
   264: |||            >                 buf(k,4)*ue(k,4))
   265: ||+----               enddo
   266: ||        
   267: ||V---->              do k=1, grid_points(3)-2
   268: |||                      km1 = k-1
   269: |||                      kp1 = k+1
   270: |||        
   271: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
   272: |||            >                 tz2*( ue(kp1,4)-ue(km1,4) )+
   273: |||            >                 dz1tz1*(ue(kp1,1)-2.0d0*ue(k,1)+ue(km1,1))
   274: |||       
   275: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - tz2 * (
   276: |||            >                 ue(kp1,2)*buf(kp1,4)-ue(km1,2)*buf(km1,4))+
   277: |||            >                 zzcon2*(buf(kp1,2)-2.0d0*buf(k,2)+buf(km1,2))+
   278: |||            >                 dz2tz1*( ue(kp1,2)-2.0d0* ue(k,2)+ ue(km1,2))
   279: |||       
   280: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - tz2 * (
   281: |||            >                 ue(kp1,3)*buf(kp1,4)-ue(km1,3)*buf(km1,4))+
   282: |||            >                 zzcon2*(buf(kp1,3)-2.0d0*buf(k,3)+buf(km1,3))+
   283: |||            >                 dz3tz1*(ue(kp1,3)-2.0d0*ue(k,3)+ue(km1,3))
   284: |||       
   285: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - tz2 * (
   286: |||            >                 (ue(kp1,4)*buf(kp1,4)+c2*(ue(kp1,5)-q(kp1)))-
   287: |||            >                 (ue(km1,4)*buf(km1,4)+c2*(ue(km1,5)-q(km1))))+
   288: |||            >                 zzcon1*(buf(kp1,4)-2.0d0*buf(k,4)+buf(km1,4))+
   289: |||            >                 dz4tz1*( ue(kp1,4)-2.0d0*ue(k,4) +ue(km1,4))
   290: |||       
   291: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - tz2 * (
   292: |||            >                 buf(kp1,4)*(c1*ue(kp1,5)-c2*q(kp1))-
   293: |||            >                 buf(km1,4)*(c1*ue(km1,5)-c2*q(km1)))+
   294: |||            >                 0.5d0*zzcon3*(buf(kp1,1)-2.0d0*buf(k,1)
   295: |||            >                 +buf(km1,1))+
   296: |||            >                 zzcon4*(cuf(kp1)-2.0d0*cuf(k)+cuf(km1))+
   297: |||            >                 zzcon5*(buf(kp1,5)-2.0d0*buf(k,5)+buf(km1,5))+
   298: |||            >                 dz5tz1*( ue(kp1,5)-2.0d0*ue(k,5)+ ue(km1,5))
   299: ||V----               enddo
   300: ||        
   301: ||        c---------------------------------------------------------------------
   302: ||        c     Fourth-order dissipation
   303: ||        c---------------------------------------------------------------------
   304: ||*---->              do m = 1, 5
   305: |||                      k = 1
   306: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   307: |||            >                    (5.0d0*ue(k,m) - 4.0d0*ue(k+1,m) +ue(k+2,m))
   308: |||                      k = 2
   309: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   310: |||            >                    (-4.0d0*ue(k-1,m) + 6.0d0*ue(k,m) -
   311: |||            >                    4.0d0*ue(k+1,m) +       ue(k+2,m))
   312: ||*----               enddo
   313: ||        
   314: ||V---->              do k = 3, grid_points(3)-4
   315: |||*--->                 do m = 1, 5
   316: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   317: ||||           >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   318: ||||           >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m) + ue(k+2,m))
   319: |||*---                  enddo
   320: ||V----               enddo
   321: ||        
   322: ||*---->              do m = 1, 5
   323: |||                      k = grid_points(3)-3
   324: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   325: |||            >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   326: |||            >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m))
   327: |||                      k = grid_points(3)-2
   328: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   329: |||            >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) + 5.0d0*ue(k,m))
   330: ||*----               enddo
   331: ||        
   332: |+-----            enddo
   333: P------         enddo
   334:           !$omp end do
   335:           
   336:           c---------------------------------------------------------------------
   337:           c     now change the sign of the forcing function,
   338:           c---------------------------------------------------------------------
   339:           !$omp do schedule(static) collapse(2)
   340: Y------>        do k = 1, grid_points(3)-2
   341: |+----->           do j = 1, grid_points(2)-2
   342: ||V---->              do i = 1, grid_points(1)-2
   343: |||*--->                 do m = 1, 5
   344: ||||                        forcing(m,i,j,k) = -1.d0 * forcing(m,i,j,k)
   345: |||*---                  enddo
   346: ||V----               enddo
   347: |+-----            enddo
   348: Y------         enddo
   349:           !$omp end do nowait
   350:           !$omp end parallel
   351:           
   352:                 return
   353:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
INLINE LIST

  ROOT: EXACT_RHS (exact_rhs.f:5)
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:48)
     *** Source for routine not found.
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:149)
     *** Source for routine not found.
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:249)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (exact_rhs.f:26)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:26)
    *** The number of VLOAD, VSTORE. :  0,  2. (exact_rhs.f:26)

    LOOP BEGIN: (exact_rhs.f:28)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:28)
      *** The number of VLOAD, VSTORE. :  0,  5. (exact_rhs.f:28)
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:40)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:45)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:45)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:67)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:67)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:67)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:67)
      *** Iteration count is assumed. Iteration count. : 163 (exact_rhs.f:67)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:115)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:115)
      *** The number of VLOAD, VSTORE. : 30,  5. (exact_rhs.f:115)
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:141)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:146)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:146)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:167)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:167)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:167)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:167)
      *** Iteration count is assumed. Iteration count. : 163 (exact_rhs.f:167)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:214)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:214)
      *** The number of VLOAD, VSTORE. : 30,  5. (exact_rhs.f:214)
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:241)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:246)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:246)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:267)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:267)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:267)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:267)
      *** Iteration count is assumed. Iteration count. : 163 (exact_rhs.f:267)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:314)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:314)
      *** The number of VLOAD, VSTORE. : 30,  5. (exact_rhs.f:314)
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:340)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:340)
    *** The number of VLOAD, VSTORE. :  0,  2. (exact_rhs.f:340)

    LOOP BEGIN: (exact_rhs.f:342)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:342)
      *** The number of VLOAD, VSTORE. :  5,  5. (exact_rhs.f:342)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:45 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 2288 bytes
      Register spill area      : 1000 bytes
      Parameter area           :   32 bytes
      Register save area       :  176 bytes
      User data area           : 1080 bytes
      Others                   :    0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (exact_rhs.f:26)
    *** Estimated execution cycle                       : 511
  LOOP END

  LOOP BEGIN: (exact_rhs.f:26)
    *** Estimated execution cycle                       : 21
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (exact_rhs.f:28)
      *** Estimated execution cycle                     : 167
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:40)
    *** Estimated execution cycle                       : 1911
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 30
            Not enough registers                        :  2
            Over basic blocks                           : 28
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 36
            Over basic blocks                           : 36
    *** The number of SCALAR REGISTER TRANSFER          : 56

    LOOP BEGIN: (exact_rhs.f:45)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 4
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 12
              Over basic blocks                         : 12
      *** The number of SCALAR REGISTER TRANSFER        : 12
    LOOP END

    LOOP BEGIN: (exact_rhs.f:115)
      *** Estimated execution cycle                     : 1224
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 9
              Over basic blocks                         : 9
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:141)
    *** Estimated execution cycle                       : 1909
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 29
            Not enough registers                        :  1
            Over basic blocks                           : 28
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 35
            Over basic blocks                           : 35
    *** The number of SCALAR REGISTER TRANSFER          : 62

    LOOP BEGIN: (exact_rhs.f:146)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 12
              Over basic blocks                         : 12
      *** The number of SCALAR REGISTER TRANSFER        : 12
    LOOP END

    LOOP BEGIN: (exact_rhs.f:214)
      *** Estimated execution cycle                     : 1224
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 9
              Over basic blocks                         : 9
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:241)
    *** Estimated execution cycle                       : 1913
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 28
            Not enough registers                        :  1
            Over basic blocks                           : 27
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 34
            Over basic blocks                           : 34
    *** The number of SCALAR REGISTER TRANSFER          : 61

    LOOP BEGIN: (exact_rhs.f:246)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 12
              Over basic blocks                         : 12
      *** The number of SCALAR REGISTER TRANSFER        : 12
    LOOP END

    LOOP BEGIN: (exact_rhs.f:314)
      *** Estimated execution cycle                     : 1224
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 9
              Over basic blocks                         : 9
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:340)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (exact_rhs.f:340)
    *** Estimated execution cycle                       : 21
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (exact_rhs.f:342)
      *** Estimated execution cycle                     : 328
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 4
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


