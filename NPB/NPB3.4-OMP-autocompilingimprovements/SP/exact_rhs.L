NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    21: par(1801): Parallel routine generated.: EXACT_RHS$1
    28: par(1803): Parallelized by "do".
    28: opt(2495): Loop nest collapsed by "collapse" clause.
    28: vec( 101): Vectorized loop.
    30: vec( 101): Vectorized loop.
    31: opt(1418): Constant-length loop is expanded.
    42: opt(3008): Reference within a conditional branch moved outside loop.
    42: par(1803): Parallelized by "do".
    42: opt(2495): Loop nest collapsed by "collapse" clause.
    47: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    47: vec( 103): Unvectorized loop.
    50: inl(1212): Source for routine not found.: EXACT_SOLUTION
    50: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    51: opt(1418): Constant-length loop is expanded.
    57: opt(1418): Constant-length loop is expanded.
    69: vec( 101): Vectorized loop.
    69: vec( 124): Iteration count is assumed. Iteration count.: 162
    73: vec( 128): Fused multiply-add operation applied.
    77: vec( 128): Fused multiply-add operation applied.
    83: vec( 128): Fused multiply-add operation applied.
    88: vec( 128): Fused multiply-add operation applied.
    93: vec( 128): Fused multiply-add operation applied.
   106: opt(1418): Constant-length loop is expanded.
   116: opt(1592): Outer loop unrolled inside inner loop.: M
   117: vec( 101): Vectorized loop.
   118: vec( 128): Fused multiply-add operation applied.
   124: opt(1418): Constant-length loop is expanded.
   142: opt(3008): Reference within a conditional branch moved outside loop.
   142: par(1803): Parallelized by "do".
   142: opt(2495): Loop nest collapsed by "collapse" clause.
   147: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   147: vec( 103): Unvectorized loop.
   150: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   151: opt(1418): Constant-length loop is expanded.
   156: opt(1418): Constant-length loop is expanded.
   167: vec( 101): Vectorized loop.
   167: vec( 124): Iteration count is assumed. Iteration count.: 162
   171: vec( 128): Fused multiply-add operation applied.
   175: vec( 128): Fused multiply-add operation applied.
   180: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   204: opt(1418): Constant-length loop is expanded.
   214: opt(1592): Outer loop unrolled inside inner loop.: M
   215: vec( 101): Vectorized loop.
   216: vec( 128): Fused multiply-add operation applied.
   222: opt(1418): Constant-length loop is expanded.
   240: opt(3008): Reference within a conditional branch moved outside loop.
   240: par(1803): Parallelized by "do".
   240: opt(2495): Loop nest collapsed by "collapse" clause.
   245: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   245: vec( 103): Unvectorized loop.
   248: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   249: opt(1418): Constant-length loop is expanded.
   255: opt(1418): Constant-length loop is expanded.
   266: vec( 101): Vectorized loop.
   266: vec( 124): Iteration count is assumed. Iteration count.: 162
   270: vec( 128): Fused multiply-add operation applied.
   274: vec( 128): Fused multiply-add operation applied.
   279: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   290: vec( 128): Fused multiply-add operation applied.
   303: opt(1418): Constant-length loop is expanded.
   313: opt(1592): Outer loop unrolled inside inner loop.: M
   314: vec( 101): Vectorized loop.
   315: vec( 128): Fused multiply-add operation applied.
   321: opt(1418): Constant-length loop is expanded.
   338: par(1803): Parallelized by "do".
   338: opt(2495): Loop nest collapsed by "collapse" clause.
   338: vec( 101): Vectorized loop.
   340: vec( 101): Vectorized loop.
   341: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
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

     5:                  subroutine exact_rhs
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c compute the right hand side based on exact solution
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  double precision dtemp(5), xi, eta, zeta, dtpp
    18:                  integer          m, i, j, k, ip1, im1, jp1,
    19:                >                  jm1, km1, kp1
    20:           
    21:           !$omp parallel default(shared)
    22:           !$omp& private(i,j,k,m,zeta,eta,xi,dtpp,im1,ip1,
    23:           !$omp&         jm1,jp1,km1,kp1,dtemp)
    24:           c---------------------------------------------------------------------
    25:           c      initialize
    26:           c---------------------------------------------------------------------
    27:           !$omp do schedule(static) collapse(2)
    28: Y------>         do   k= 0, grid_points(3)-1
    29: |+----->            do   j = 0, grid_points(2)-1
    30: ||V---->               do   i = 0, grid_points(1)-1
    31: |||*--->                  do   m = 1, 5
    32: ||||                         forcing(m,i,j,k) = 0.0d0
    33: |||*---                   end do
    34: ||V----                end do
    35: |+-----             end do
    36: Y------          end do
    37:           
    38:           c---------------------------------------------------------------------
    39:           c      xi-direction flux differences
    40:           c---------------------------------------------------------------------
    41:           !$omp do schedule(static) collapse(2)
    42: P------>         do   k = 1, grid_points(3)-2
    43: |+----->            do   j = 1, grid_points(2)-2
    44: ||                  zeta = dble(k) * dnzm1
    45: ||                     eta = dble(j) * dnym1
    46: ||        
    47: ||+---->               do  i=0, grid_points(1)-1
    48: |||                       xi = dble(i) * dnxm1
    49: |||       
    50: |||                       call exact_solution(xi, eta, zeta, dtemp)
    51: |||*--->                  do  m = 1, 5
    52: ||||                         ue(i,m) = dtemp(m)
    53: |||*---                   end do
    54: |||       
    55: |||                       dtpp = 1.0d0 / dtemp(1)
    56: |||       
    57: |||*--->                  do  m = 2, 5
    58: ||||                         buf(i,m) = dtpp * dtemp(m)
    59: |||*---                   end do
    60: |||       
    61: |||                       cuf(i)   = buf(i,2) * buf(i,2)
    62: |||                       buf(i,1) = cuf(i) + buf(i,3) * buf(i,3) +
    63: |||            >                     buf(i,4) * buf(i,4)
    64: |||                       q(i) = 0.5d0*(buf(i,2)*ue(i,2) + buf(i,3)*ue(i,3) +
    65: |||            >                        buf(i,4)*ue(i,4))
    66: |||       
    67: ||+----                end do
    68: ||         
    69: ||V---->               do  i = 1, grid_points(1)-2
    70: |||                       im1 = i-1
    71: |||                       ip1 = i+1
    72: |||       
    73: |||     F                 forcing(1,i,j,k) = forcing(1,i,j,k) -
    74: |||            >                 tx2*( ue(ip1,2)-ue(im1,2) )+
    75: |||            >                 dx1tx1*(ue(ip1,1)-2.0d0*ue(i,1)+ue(im1,1))
    76: |||       
    77: |||     F                 forcing(2,i,j,k) = forcing(2,i,j,k) - tx2 * (
    78: |||            >                (ue(ip1,2)*buf(ip1,2)+c2*(ue(ip1,5)-q(ip1)))-
    79: |||            >                (ue(im1,2)*buf(im1,2)+c2*(ue(im1,5)-q(im1))))+
    80: |||            >                 xxcon1*(buf(ip1,2)-2.0d0*buf(i,2)+buf(im1,2))+
    81: |||            >                 dx2tx1*( ue(ip1,2)-2.0d0* ue(i,2)+ue(im1,2))
    82: |||       
    83: |||     F                 forcing(3,i,j,k) = forcing(3,i,j,k) - tx2 * (
    84: |||            >                 ue(ip1,3)*buf(ip1,2)-ue(im1,3)*buf(im1,2))+
    85: |||            >                 xxcon2*(buf(ip1,3)-2.0d0*buf(i,3)+buf(im1,3))+
    86: |||            >                 dx3tx1*( ue(ip1,3)-2.0d0*ue(i,3) +ue(im1,3))
    87: |||        
    88: |||     F                 forcing(4,i,j,k) = forcing(4,i,j,k) - tx2*(
    89: |||            >                 ue(ip1,4)*buf(ip1,2)-ue(im1,4)*buf(im1,2))+
    90: |||            >                 xxcon2*(buf(ip1,4)-2.0d0*buf(i,4)+buf(im1,4))+
    91: |||            >                 dx4tx1*( ue(ip1,4)-2.0d0* ue(i,4)+ ue(im1,4))
    92: |||       
    93: |||     F                 forcing(5,i,j,k) = forcing(5,i,j,k) - tx2*(
    94: |||            >                 buf(ip1,2)*(c1*ue(ip1,5)-c2*q(ip1))-
    95: |||            >                 buf(im1,2)*(c1*ue(im1,5)-c2*q(im1)))+
    96: |||            >                 0.5d0*xxcon3*(buf(ip1,1)-2.0d0*buf(i,1)+
    97: |||            >                               buf(im1,1))+
    98: |||            >                 xxcon4*(cuf(ip1)-2.0d0*cuf(i)+cuf(im1))+
    99: |||            >                 xxcon5*(buf(ip1,5)-2.0d0*buf(i,5)+buf(im1,5))+
   100: |||            >                 dx5tx1*( ue(ip1,5)-2.0d0* ue(i,5)+ ue(im1,5))
   101: ||V----                end do
   102: ||        
   103: ||        c---------------------------------------------------------------------
   104: ||        c            Fourth-order dissipation
   105: ||        c---------------------------------------------------------------------
   106: ||*---->               do   m = 1, 5
   107: |||                       i = 1
   108: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   109: |||            >                    (5.0d0*ue(i,m) - 4.0d0*ue(i+1,m) +ue(i+2,m))
   110: |||                       i = 2
   111: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   112: |||            >                   (-4.0d0*ue(i-1,m) + 6.0d0*ue(i,m) -
   113: |||            >                     4.0d0*ue(i+1,m) +       ue(i+2,m))
   114: ||*----                end do
   115: ||        
   116: ||U---->               do   m = 1, 5
   117: |||V--->                  do  i = 3, grid_points(1)-4
   118: ||||    F                    forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   119: ||||           >                   (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   120: ||||           >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m) + ue(i+2,m))
   121: |||V---                   end do
   122: ||U----                end do
   123: ||        
   124: ||*---->               do   m = 1, 5
   125: |||                       i = grid_points(1)-3
   126: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   127: |||            >                   (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   128: |||            >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m))
   129: |||                       i = grid_points(1)-2
   130: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   131: |||            >                   (ue(i-2,m) - 4.0d0*ue(i-1,m) + 5.0d0*ue(i,m))
   132: ||*----                end do
   133: ||        
   134: |+-----             end do
   135: P------          end do
   136:           !$omp end do
   137:           
   138:           c---------------------------------------------------------------------
   139:           c  eta-direction flux differences
   140:           c---------------------------------------------------------------------
   141:           !$omp do schedule(static) collapse(2)
   142: P------>         do   k = 1, grid_points(3)-2
   143: |+----->            do   i=1, grid_points(1)-2
   144: ||                  zeta = dble(k) * dnzm1
   145: ||                     xi = dble(i) * dnxm1
   146: ||        
   147: ||+---->               do  j=0, grid_points(2)-1
   148: |||                       eta = dble(j) * dnym1
   149: |||       
   150: |||                       call exact_solution(xi, eta, zeta, dtemp)
   151: |||*--->                  do   m = 1, 5
   152: ||||                         ue(j,m) = dtemp(m)
   153: |||*---                   end do
   154: |||                       dtpp = 1.0d0/dtemp(1)
   155: |||       
   156: |||*--->                  do  m = 2, 5
   157: ||||                         buf(j,m) = dtpp * dtemp(m)
   158: |||*---                   end do
   159: |||       
   160: |||                       cuf(j)   = buf(j,3) * buf(j,3)
   161: |||                       buf(j,1) = cuf(j) + buf(j,2) * buf(j,2) +
   162: |||            >                     buf(j,4) * buf(j,4)
   163: |||                       q(j) = 0.5d0*(buf(j,2)*ue(j,2) + buf(j,3)*ue(j,3) +
   164: |||            >                        buf(j,4)*ue(j,4))
   165: ||+----                end do
   166: ||        
   167: ||V---->               do  j = 1, grid_points(2)-2
   168: |||                       jm1 = j-1
   169: |||                       jp1 = j+1
   170: |||        
   171: |||     F                 forcing(1,i,j,k) = forcing(1,i,j,k) -
   172: |||            >                ty2*( ue(jp1,3)-ue(jm1,3) )+
   173: |||            >                dy1ty1*(ue(jp1,1)-2.0d0*ue(j,1)+ue(jm1,1))
   174: |||       
   175: |||     F                 forcing(2,i,j,k) = forcing(2,i,j,k) - ty2*(
   176: |||            >                ue(jp1,2)*buf(jp1,3)-ue(jm1,2)*buf(jm1,3))+
   177: |||            >                yycon2*(buf(jp1,2)-2.0d0*buf(j,2)+buf(jm1,2))+
   178: |||            >                dy2ty1*( ue(jp1,2)-2.0* ue(j,2)+ ue(jm1,2))
   179: |||       
   180: |||     F                 forcing(3,i,j,k) = forcing(3,i,j,k) - ty2*(
   181: |||            >                (ue(jp1,3)*buf(jp1,3)+c2*(ue(jp1,5)-q(jp1)))-
   182: |||            >                (ue(jm1,3)*buf(jm1,3)+c2*(ue(jm1,5)-q(jm1))))+
   183: |||            >                yycon1*(buf(jp1,3)-2.0d0*buf(j,3)+buf(jm1,3))+
   184: |||            >                dy3ty1*( ue(jp1,3)-2.0d0*ue(j,3) +ue(jm1,3))
   185: |||       
   186: |||     F                 forcing(4,i,j,k) = forcing(4,i,j,k) - ty2*(
   187: |||            >                ue(jp1,4)*buf(jp1,3)-ue(jm1,4)*buf(jm1,3))+
   188: |||            >                yycon2*(buf(jp1,4)-2.0d0*buf(j,4)+buf(jm1,4))+
   189: |||            >                dy4ty1*( ue(jp1,4)-2.0d0*ue(j,4)+ ue(jm1,4))
   190: |||       
   191: |||     F                 forcing(5,i,j,k) = forcing(5,i,j,k) - ty2*(
   192: |||            >                buf(jp1,3)*(c1*ue(jp1,5)-c2*q(jp1))-
   193: |||            >                buf(jm1,3)*(c1*ue(jm1,5)-c2*q(jm1)))+
   194: |||            >                0.5d0*yycon3*(buf(jp1,1)-2.0d0*buf(j,1)+
   195: |||            >                              buf(jm1,1))+
   196: |||            >                yycon4*(cuf(jp1)-2.0d0*cuf(j)+cuf(jm1))+
   197: |||            >                yycon5*(buf(jp1,5)-2.0d0*buf(j,5)+buf(jm1,5))+
   198: |||            >                dy5ty1*(ue(jp1,5)-2.0d0*ue(j,5)+ue(jm1,5))
   199: ||V----                end do
   200: ||        
   201: ||        c---------------------------------------------------------------------
   202: ||        c            Fourth-order dissipation
   203: ||        c---------------------------------------------------------------------
   204: ||*---->               do   m = 1, 5
   205: |||                       j = 1
   206: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   207: |||            >                    (5.0d0*ue(j,m) - 4.0d0*ue(j+1,m) +ue(j+2,m))
   208: |||                       j = 2
   209: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   210: |||            >                   (-4.0d0*ue(j-1,m) + 6.0d0*ue(j,m) -
   211: |||            >                     4.0d0*ue(j+1,m) +       ue(j+2,m))
   212: ||*----                end do
   213: ||        
   214: ||U---->               do   m = 1, 5
   215: |||V--->                  do  j = 3, grid_points(2)-4
   216: ||||    F                    forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   217: ||||           >                   (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   218: ||||           >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m) + ue(j+2,m))
   219: |||V---                   end do
   220: ||U----                end do
   221: ||        
   222: ||*---->               do   m = 1, 5
   223: |||                       j = grid_points(2)-3
   224: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   225: |||            >                   (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   226: |||            >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m))
   227: |||                       j = grid_points(2)-2
   228: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   229: |||            >                   (ue(j-2,m) - 4.0d0*ue(j-1,m) + 5.0d0*ue(j,m))
   230: |||       
   231: ||*----                end do
   232: ||        
   233: |+-----             end do
   234: P------          end do
   235:           
   236:           c---------------------------------------------------------------------
   237:           c      zeta-direction flux differences
   238:           c---------------------------------------------------------------------
   239:           !$omp do schedule(static) collapse(2)
   240: P------>         do  j=1, grid_points(2)-2
   241: |+----->            do   i = 1, grid_points(1)-2
   242: ||                  eta = dble(j) * dnym1
   243: ||                     xi = dble(i) * dnxm1
   244: ||        
   245: ||+---->               do k=0, grid_points(3)-1
   246: |||                       zeta = dble(k) * dnzm1
   247: |||       
   248: |||                       call exact_solution(xi, eta, zeta, dtemp)
   249: |||*--->                  do   m = 1, 5
   250: ||||                         ue(k,m) = dtemp(m)
   251: |||*---                   end do
   252: |||       
   253: |||                       dtpp = 1.0d0/dtemp(1)
   254: |||       
   255: |||*--->                  do   m = 2, 5
   256: ||||                         buf(k,m) = dtpp * dtemp(m)
   257: |||*---                   end do
   258: |||       
   259: |||                       cuf(k)   = buf(k,4) * buf(k,4)
   260: |||                       buf(k,1) = cuf(k) + buf(k,2) * buf(k,2) +
   261: |||            >                     buf(k,3) * buf(k,3)
   262: |||                       q(k) = 0.5d0*(buf(k,2)*ue(k,2) + buf(k,3)*ue(k,3) +
   263: |||            >                        buf(k,4)*ue(k,4))
   264: ||+----                end do
   265: ||        
   266: ||V---->               do    k=1, grid_points(3)-2
   267: |||                       km1 = k-1
   268: |||                       kp1 = k+1
   269: |||       
   270: |||     F                 forcing(1,i,j,k) = forcing(1,i,j,k) -
   271: |||            >                 tz2*( ue(kp1,4)-ue(km1,4) )+
   272: |||            >                 dz1tz1*(ue(kp1,1)-2.0d0*ue(k,1)+ue(km1,1))
   273: |||       
   274: |||     F                 forcing(2,i,j,k) = forcing(2,i,j,k) - tz2 * (
   275: |||            >                 ue(kp1,2)*buf(kp1,4)-ue(km1,2)*buf(km1,4))+
   276: |||            >                 zzcon2*(buf(kp1,2)-2.0d0*buf(k,2)+buf(km1,2))+
   277: |||            >                 dz2tz1*( ue(kp1,2)-2.0d0* ue(k,2)+ ue(km1,2))
   278: |||       
   279: |||     F                 forcing(3,i,j,k) = forcing(3,i,j,k) - tz2 * (
   280: |||            >                 ue(kp1,3)*buf(kp1,4)-ue(km1,3)*buf(km1,4))+
   281: |||            >                 zzcon2*(buf(kp1,3)-2.0d0*buf(k,3)+buf(km1,3))+
   282: |||            >                 dz3tz1*(ue(kp1,3)-2.0d0*ue(k,3)+ue(km1,3))
   283: |||       
   284: |||     F                 forcing(4,i,j,k) = forcing(4,i,j,k) - tz2 * (
   285: |||            >                (ue(kp1,4)*buf(kp1,4)+c2*(ue(kp1,5)-q(kp1)))-
   286: |||            >                (ue(km1,4)*buf(km1,4)+c2*(ue(km1,5)-q(km1))))+
   287: |||            >                zzcon1*(buf(kp1,4)-2.0d0*buf(k,4)+buf(km1,4))+
   288: |||            >                dz4tz1*( ue(kp1,4)-2.0d0*ue(k,4) +ue(km1,4))
   289: |||       
   290: |||     F                 forcing(5,i,j,k) = forcing(5,i,j,k) - tz2 * (
   291: |||            >                 buf(kp1,4)*(c1*ue(kp1,5)-c2*q(kp1))-
   292: |||            >                 buf(km1,4)*(c1*ue(km1,5)-c2*q(km1)))+
   293: |||            >                 0.5d0*zzcon3*(buf(kp1,1)-2.0d0*buf(k,1)
   294: |||            >                              +buf(km1,1))+
   295: |||            >                 zzcon4*(cuf(kp1)-2.0d0*cuf(k)+cuf(km1))+
   296: |||            >                 zzcon5*(buf(kp1,5)-2.0d0*buf(k,5)+buf(km1,5))+
   297: |||            >                 dz5tz1*( ue(kp1,5)-2.0d0*ue(k,5)+ ue(km1,5))
   298: ||V----                end do
   299: ||        
   300: ||        c---------------------------------------------------------------------
   301: ||        c            Fourth-order dissipation
   302: ||        c---------------------------------------------------------------------
   303: ||*---->               do   m = 1, 5
   304: |||                       k = 1
   305: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   306: |||            >                    (5.0d0*ue(k,m) - 4.0d0*ue(k+1,m) +ue(k+2,m))
   307: |||                       k = 2
   308: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   309: |||            >                   (-4.0d0*ue(k-1,m) + 6.0d0*ue(k,m) -
   310: |||            >                     4.0d0*ue(k+1,m) +       ue(k+2,m))
   311: ||*----                end do
   312: ||        
   313: ||U---->               do   m = 1, 5
   314: |||V--->                  do  k = 3, grid_points(3)-4
   315: ||||    F                    forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   316: ||||           >                   (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   317: ||||           >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m) + ue(k+2,m))
   318: |||V---                   end do
   319: ||U----                end do
   320: ||        
   321: ||*---->               do    m = 1, 5
   322: |||                       k = grid_points(3)-3
   323: |||                       forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   324: |||            >                   (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   325: |||            >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m))
   326: |||                          k = grid_points(3)-2
   327: |||                          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   328: |||            >                   (ue(k-2,m) - 4.0d0*ue(k-1,m) + 5.0d0*ue(k,m))
   329: ||*----                   end do
   330: ||        
   331: |+-----             end do
   332: P------          end do
   333:           
   334:           c---------------------------------------------------------------------
   335:           c now change the sign of the forcing function,
   336:           c---------------------------------------------------------------------
   337:           !$omp do schedule(static) collapse(2)
   338: Y------>         do   k = 1, grid_points(3)-2
   339: |+----->            do   j = 1, grid_points(2)-2
   340: ||V---->               do   i = 1, grid_points(1)-2
   341: |||*--->                  do   m = 1, 5
   342: ||||                         forcing(m,i,j,k) = -1.d0 * forcing(m,i,j,k)
   343: |||*---                   end do
   344: ||V----                end do
   345: |+-----             end do
   346: Y------          end do
   347:           !$omp end do nowait
   348:           !$omp end parallel
   349:           
   350:                  return
   351:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
INLINE LIST

  ROOT: EXACT_RHS (exact_rhs.f:5)
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:50)
     *** Source for routine not found.
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:150)
     *** Source for routine not found.
  -> NOINLINE: EXACT_SOLUTION (exact_rhs.f:248)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (exact_rhs.f:28)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:28)
    *** The number of VLOAD, VSTORE. :  0,  2. (exact_rhs.f:28)

    LOOP BEGIN: (exact_rhs.f:30)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:30)
      *** The number of VLOAD, VSTORE. :  0,  5. (exact_rhs.f:30)
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:42)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:47)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:47)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:69)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:69)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:69)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:69)
      *** Iteration count is assumed. Iteration count. : 162 (exact_rhs.f:69)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:117)
      <Vectorized loop.>
      **  Copied loop. (exact_rhs.f:117)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:117)
      *** The number of VLOAD, VSTORE. :  6,  1. (exact_rhs.f:117)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:116)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (exact_rhs.f:116)

      LOOP BEGIN: (exact_rhs.f:117)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:117)
        *** The number of VLOAD, VSTORE. : 24,  4. (exact_rhs.f:117)
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:142)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:147)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:147)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:167)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:167)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:167)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:167)
      *** Iteration count is assumed. Iteration count. : 162 (exact_rhs.f:167)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:215)
      <Vectorized loop.>
      **  Copied loop. (exact_rhs.f:215)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:215)
      *** The number of VLOAD, VSTORE. :  6,  1. (exact_rhs.f:215)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:214)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (exact_rhs.f:214)

      LOOP BEGIN: (exact_rhs.f:215)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:215)
        *** The number of VLOAD, VSTORE. : 24,  4. (exact_rhs.f:215)
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:240)
    <Unvectorized loop.>

    LOOP BEGIN: (exact_rhs.f:245)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT_SOLUTION (exact_rhs.f:245)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:266)
      <Vectorized loop.>
      **  Runtime-tested short loop. (exact_rhs.f:266)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:266)
      *** The number of VLOAD, VSTORE. : 40,  5. (exact_rhs.f:266)
      *** Iteration count is assumed. Iteration count. : 162 (exact_rhs.f:266)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:314)
      <Vectorized loop.>
      **  Copied loop. (exact_rhs.f:314)
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:314)
      *** The number of VLOAD, VSTORE. :  6,  1. (exact_rhs.f:314)
    LOOP END

    LOOP BEGIN: (exact_rhs.f:313)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (exact_rhs.f:313)

      LOOP BEGIN: (exact_rhs.f:314)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:314)
        *** The number of VLOAD, VSTORE. : 24,  4. (exact_rhs.f:314)
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:338)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:338)
    *** The number of VLOAD, VSTORE. :  0,  2. (exact_rhs.f:338)

    LOOP BEGIN: (exact_rhs.f:340)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (exact_rhs.f:340)
      *** The number of VLOAD, VSTORE. :  5,  5. (exact_rhs.f:340)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:56 2021
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
    Total size                 : 3520 bytes
      Register spill area      : 1616 bytes
      Parameter area           :   32 bytes
      Register save area       :  176 bytes
      User data area           : 1696 bytes
      Others                   :    0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (exact_rhs.f:28)
    *** Estimated execution cycle                       : 511
  LOOP END

  LOOP BEGIN: (exact_rhs.f:28)
    *** Estimated execution cycle                       : 21
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (exact_rhs.f:30)
      *** Estimated execution cycle                     : 167
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:42)
    *** Estimated execution cycle                       : 1911
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 29
            Over basic blocks                           : 29
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 42
            Over basic blocks                           : 42
    *** The number of SCALAR REGISTER TRANSFER          : 46

    LOOP BEGIN: (exact_rhs.f:47)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 7
              Over basic blocks                         : 7
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 10
              Over basic blocks                         : 10
      *** The number of SCALAR REGISTER TRANSFER        : 11
    LOOP END

    LOOP BEGIN: (exact_rhs.f:117)
      *** Estimated execution cycle                     : 360
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END

    LOOP BEGIN: (exact_rhs.f:116)
      *** Estimated execution cycle                     : 46
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 7
              Over basic blocks                         : 7
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER TRANSFER        : 11

      LOOP BEGIN: (exact_rhs.f:117)
        *** Estimated execution cycle                   : 1417
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 15
                Not enough registers                    : 15
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 37
                Over basic blocks                       : 37
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:142)
    *** Estimated execution cycle                       : 1909
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 29
            Over basic blocks                           : 29
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 42
            Over basic blocks                           : 42
    *** The number of SCALAR REGISTER TRANSFER          : 46

    LOOP BEGIN: (exact_rhs.f:147)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 11
              Over basic blocks                         : 11
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 10
              Over basic blocks                         : 10
      *** The number of SCALAR REGISTER TRANSFER        : 11
    LOOP END

    LOOP BEGIN: (exact_rhs.f:215)
      *** Estimated execution cycle                     : 360
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END

    LOOP BEGIN: (exact_rhs.f:214)
      *** Estimated execution cycle                     : 46
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 7
              Over basic blocks                         : 7
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER TRANSFER        : 11

      LOOP BEGIN: (exact_rhs.f:215)
        *** Estimated execution cycle                   : 1417
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 15
                Not enough registers                    : 15
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 37
                Over basic blocks                       : 37
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:240)
    *** Estimated execution cycle                       : 1913
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 28
            Over basic blocks                           : 28
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 41
            Over basic blocks                           : 41
    *** The number of SCALAR REGISTER TRANSFER          : 46

    LOOP BEGIN: (exact_rhs.f:245)
      *** Estimated execution cycle                     : 78
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 11
              Over basic blocks                         : 11
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 10
              Over basic blocks                         : 10
      *** The number of SCALAR REGISTER TRANSFER        : 11
    LOOP END

    LOOP BEGIN: (exact_rhs.f:314)
      *** Estimated execution cycle                     : 360
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END

    LOOP BEGIN: (exact_rhs.f:313)
      *** Estimated execution cycle                     : 46
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 7
              Over basic blocks                         : 7
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER TRANSFER        : 11

      LOOP BEGIN: (exact_rhs.f:314)
        *** Estimated execution cycle                   : 1417
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 15
                Not enough registers                    : 15
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 37
                Over basic blocks                       : 37
      LOOP END
    LOOP END
  LOOP END

  LOOP BEGIN: (exact_rhs.f:338)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (exact_rhs.f:338)
    *** Estimated execution cycle                       : 20
    *** The number of SCALAR REGISTER TRANSFER          : 4

    LOOP BEGIN: (exact_rhs.f:340)
      *** Estimated execution cycle                     : 328
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 4
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END


