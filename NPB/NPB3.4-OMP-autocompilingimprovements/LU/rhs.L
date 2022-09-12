NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME : rhs.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: rhs.f

PROCEDURE NAME: RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    32: inl(1212): Source for routine not found.: TIMER_START
    34: par(1801): Parallel routine generated.: RHS$1
    39: par(1803): Parallelized by "do".
    39: opt(2495): Loop nest collapsed by "collapse" clause.
    39: vec( 101): Vectorized loop.
    41: vec( 101): Vectorized loop.
    42: opt(1418): Constant-length loop is expanded.
    47: vec( 128): Fused multiply-add operation applied.
    56: par(1805): Master-construct generated.
    64: opt(3008): Reference within a conditional branch moved outside loop.
    64: par(1803): Parallelized by "do".
    64: opt(2495): Loop nest collapsed by "collapse" clause.
    64: vec( 102): Partially vectorized loop.
    66: vec( 101): Vectorized loop.
    66: vec( 124): Iteration count is assumed. Iteration count.: 162
    72: vec( 128): Fused multiply-add operation applied.
    76: vec( 128): Fused multiply-add operation applied.
    79: vec( 101): Vectorized loop.
    80: opt(1418): Constant-length loop is expanded.
    81: vec( 128): Fused multiply-add operation applied.
    86: vec( 101): Vectorized loop.
   104: vec( 128): Fused multiply-add operation applied.
   112: vec( 101): Vectorized loop.
   113: vec( 128): Fused multiply-add operation applied.
   117: vec( 128): Fused multiply-add operation applied.
   122: vec( 128): Fused multiply-add operation applied.
   127: vec( 128): Fused multiply-add operation applied.
   132: vec( 128): Fused multiply-add operation applied.
   142: opt(1418): Constant-length loop is expanded.
   154: vec( 101): Vectorized loop.
   155: opt(1418): Constant-length loop is expanded.
   156: vec( 128): Fused multiply-add operation applied.
   166: opt(1418): Constant-length loop is expanded.
   181: par(1805): Master-construct generated.
   182: inl(1212): Source for routine not found.: TIMER_STOP
   190: opt(3008): Reference within a conditional branch moved outside loop.
   190: par(1803): Parallelized by "do".
   190: opt(2495): Loop nest collapsed by "collapse" clause.
   190: vec( 101): Vectorized loop.
   192: vec( 101): Vectorized loop.
   192: vec( 124): Iteration count is assumed. Iteration count.: 162
   199: vec( 128): Fused multiply-add operation applied.
   201: vec( 128): Fused multiply-add operation applied.
   204: vec( 101): Vectorized loop.
   205: opt(1418): Constant-length loop is expanded.
   206: vec( 128): Fused multiply-add operation applied.
   211: vec( 101): Vectorized loop.
   228: vec( 128): Fused multiply-add operation applied.
   236: vec( 101): Vectorized loop.
   238: vec( 128): Fused multiply-add operation applied.
   243: vec( 128): Fused multiply-add operation applied.
   249: vec( 128): Fused multiply-add operation applied.
   255: vec( 128): Fused multiply-add operation applied.
   261: vec( 128): Fused multiply-add operation applied.
   277: par(1803): Parallelized by "do".
   277: opt(2495): Loop nest collapsed by "collapse" clause.
   277: vec( 102): Partially vectorized loop.
   280: vec( 101): Vectorized loop.
   281: opt(1418): Constant-length loop is expanded.
   282: vec( 128): Fused multiply-add operation applied.
   290: vec( 101): Vectorized loop.
   291: opt(1418): Constant-length loop is expanded.
   292: vec( 128): Fused multiply-add operation applied.
   301: vec( 101): Vectorized loop.
   302: opt(1418): Constant-length loop is expanded.
   303: vec( 128): Fused multiply-add operation applied.
   312: vec( 101): Vectorized loop.
   313: opt(1418): Constant-length loop is expanded.
   314: vec( 128): Fused multiply-add operation applied.
   322: vec( 101): Vectorized loop.
   323: opt(1418): Constant-length loop is expanded.
   324: vec( 128): Fused multiply-add operation applied.
   337: par(1805): Master-construct generated.
   346: opt(3008): Reference within a conditional branch moved outside loop.
   346: par(1803): Parallelized by "do".
   346: opt(2495): Loop nest collapsed by "collapse" clause.
   346: vec( 102): Partially vectorized loop.
   348: opt(1112): Loop fused with previous loop.
   348: vec( 101): Vectorized loop.
   348: vec( 124): Iteration count is assumed. Iteration count.: 162
   356: opt(1112): Loop fused with previous loop.: K
   364: vec( 128): Fused multiply-add operation applied.
   365: vec( 128): Fused multiply-add operation applied.
   368: vec( 101): Vectorized loop.
   368: vec( 124): Iteration count is assumed. Iteration count.: 162
   369: opt(1418): Constant-length loop is expanded.
   370: vec( 128): Fused multiply-add operation applied.
   375: vec( 101): Vectorized loop.
   375: vec( 124): Iteration count is assumed. Iteration count.: 162
   393: vec( 128): Fused multiply-add operation applied.
   401: vec( 101): Vectorized loop.
   401: vec( 124): Iteration count is assumed. Iteration count.: 162
   402: vec( 128): Fused multiply-add operation applied.
   406: vec( 128): Fused multiply-add operation applied.
   411: vec( 128): Fused multiply-add operation applied.
   416: vec( 128): Fused multiply-add operation applied.
   421: vec( 128): Fused multiply-add operation applied.
   431: opt(1418): Constant-length loop is expanded.
   443: vec( 101): Vectorized loop.
   443: vec( 124): Iteration count is assumed. Iteration count.: 162
   444: opt(1418): Constant-length loop is expanded.
   445: vec( 128): Fused multiply-add operation applied.
   454: opt(1418): Constant-length loop is expanded.
   468: par(1805): Master-construct generated.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: rhs.f

PROCEDURE NAME: RHS
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

     5:                 subroutine rhs
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   compute the right hand sides
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use lu_data
    15:                 implicit none
    16:           
    17:           c---------------------------------------------------------------------
    18:           c  local variables
    19:           c---------------------------------------------------------------------
    20:                 integer i, j, k, m
    21:                 double precision  q
    22:                 double precision  tmp, utmp(6,isiz3), rtmp(5,isiz3)
    23:                 double precision  u21, u31, u41
    24:                 double precision  u21i, u31i, u41i, u51i
    25:                 double precision  u21j, u31j, u41j, u51j
    26:                 double precision  u21k, u31k, u41k, u51k
    27:                 double precision  u21im1, u31im1, u41im1, u51im1
    28:                 double precision  u21jm1, u31jm1, u41jm1, u51jm1
    29:                 double precision  u21km1, u31km1, u41km1, u51km1
    30:           
    31:           
    32:                 if (timeron) call timer_start(t_rhs)
    33:           
    34:           !$omp parallel default(shared) private(i,j,k,m,q,tmp,utmp,rtmp,
    35:           !$omp& u51im1,u41im1,u31im1,u21im1,u51i,u41i,u31i,u21i,u21,
    36:           !$omp& u51jm1,u41jm1,u31jm1,u21jm1,u51j,u41j,u31j,u21j,u31,
    37:           !$omp& u51km1,u41km1,u31km1,u21km1,u51k,u41k,u31k,u21k,u41)
    38:           !$omp do schedule(static) collapse(2)
    39: Y------>        do k = 1, nz
    40: |+----->           do j = 1, ny
    41: ||V---->              do i = 1, nx
    42: |||*--->                 do m = 1, 5
    43: ||||                        rsd(m,i,j,k) = - frct(m,i,j,k)
    44: |||*---                  end do
    45: |||                      tmp = 1.0d+00 / u(1,i,j,k)
    46: |||                      rho_i(i,j,k) = tmp
    47: |||     F                qs(i,j,k) = 0.50d+00 * (  u(2,i,j,k) * u(2,i,j,k)
    48: |||            >                         + u(3,i,j,k) * u(3,i,j,k)
    49: |||            >                         + u(4,i,j,k) * u(4,i,j,k) )
    50: |||            >                      * tmp
    51: ||V----               end do
    52: |+-----            end do
    53: Y------         end do
    54:           !$omp end do
    55:           
    56:           !$omp master
    57:                 if (timeron) call timer_start(t_rhsx)
    58:           !$omp end master
    59:           c---------------------------------------------------------------------
    60:           c   xi-direction flux differences
    61:           c---------------------------------------------------------------------
    62:           
    63:           !$omp do schedule(static) collapse(2)
    64: Y------>        do k = 2, nz - 1
    65: |+----->           do j = jst, jend
    66: ||V---->              do i = 1, nx
    67: |||                      flux(1,i) = u(2,i,j,k)
    68: |||                      u21 = u(2,i,j,k) * rho_i(i,j,k)
    69: |||       
    70: |||                      q = qs(i,j,k)
    71: |||       
    72: |||     F                flux(2,i) = u(2,i,j,k) * u21 + c2 *
    73: |||            >                        ( u(5,i,j,k) - q )
    74: |||                      flux(3,i) = u(3,i,j,k) * u21
    75: |||                      flux(4,i) = u(4,i,j,k) * u21
    76: |||     F                flux(5,i) = ( c1 * u(5,i,j,k) - c2 * q ) * u21
    77: ||V----               end do
    78: ||        
    79: ||V---->              do i = ist, iend
    80: |||*--->                 do m = 1, 5
    81: ||||    F                   rsd(m,i,j,k) =  rsd(m,i,j,k)
    82: ||||           >                 - tx2 * ( flux(m,i+1) - flux(m,i-1) )
    83: |||*---                  end do
    84: ||V----               end do
    85: ||        
    86: ||V---->              do i = ist, nx
    87: |||                      tmp = rho_i(i,j,k)
    88: |||       
    89: |||                      u21i = tmp * u(2,i,j,k)
    90: |||                      u31i = tmp * u(3,i,j,k)
    91: |||                      u41i = tmp * u(4,i,j,k)
    92: |||                      u51i = tmp * u(5,i,j,k)
    93: |||       
    94: |||                      tmp = rho_i(i-1,j,k)
    95: |||       
    96: |||                      u21im1 = tmp * u(2,i-1,j,k)
    97: |||                      u31im1 = tmp * u(3,i-1,j,k)
    98: |||                      u41im1 = tmp * u(4,i-1,j,k)
    99: |||                      u51im1 = tmp * u(5,i-1,j,k)
   100: |||       
   101: |||                      flux(2,i) = (4.0d+00/3.0d+00) * tx3 * (u21i-u21im1)
   102: |||                      flux(3,i) = tx3 * ( u31i - u31im1 )
   103: |||                      flux(4,i) = tx3 * ( u41i - u41im1 )
   104: |||     F                flux(5,i) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   105: |||            >              * tx3 * ( ( u21i  **2 + u31i  **2 + u41i  **2 )
   106: |||            >                      - ( u21im1**2 + u31im1**2 + u41im1**2 ) )
   107: |||            >              + (1.0d+00/6.0d+00)
   108: |||            >              * tx3 * ( u21i**2 - u21im1**2 )
   109: |||            >              + c1 * c5 * tx3 * ( u51i - u51im1 )
   110: ||V----               end do
   111: ||        
   112: ||V---->              do i = ist, iend
   113: |||     F                rsd(1,i,j,k) = rsd(1,i,j,k)
   114: |||            >              + dx1 * tx1 * (            u(1,i-1,j,k)
   115: |||            >                             - 2.0d+00 * u(1,i,j,k)
   116: |||            >                             +           u(1,i+1,j,k) )
   117: |||     F                rsd(2,i,j,k) = rsd(2,i,j,k)
   118: |||            >          + tx3 * c3 * c4 * ( flux(2,i+1) - flux(2,i) )
   119: |||            >              + dx2 * tx1 * (            u(2,i-1,j,k)
   120: |||            >                             - 2.0d+00 * u(2,i,j,k)
   121: |||            >                             +           u(2,i+1,j,k) )
   122: |||     F                rsd(3,i,j,k) = rsd(3,i,j,k)
   123: |||            >          + tx3 * c3 * c4 * ( flux(3,i+1) - flux(3,i) )
   124: |||            >              + dx3 * tx1 * (            u(3,i-1,j,k)
   125: |||            >                             - 2.0d+00 * u(3,i,j,k)
   126: |||            >                             +           u(3,i+1,j,k) )
   127: |||     F                rsd(4,i,j,k) = rsd(4,i,j,k)
   128: |||            >          + tx3 * c3 * c4 * ( flux(4,i+1) - flux(4,i) )
   129: |||            >              + dx4 * tx1 * (            u(4,i-1,j,k)
   130: |||            >                             - 2.0d+00 * u(4,i,j,k)
   131: |||            >                             +           u(4,i+1,j,k) )
   132: |||     F                rsd(5,i,j,k) = rsd(5,i,j,k)
   133: |||            >          + tx3 * c3 * c4 * ( flux(5,i+1) - flux(5,i) )
   134: |||            >              + dx5 * tx1 * (            u(5,i-1,j,k)
   135: |||            >                             - 2.0d+00 * u(5,i,j,k)
   136: |||            >                             +           u(5,i+1,j,k) )
   137: ||V----               end do
   138: ||        
   139: ||        c---------------------------------------------------------------------
   140: ||        c   Fourth-order dissipation
   141: ||        c---------------------------------------------------------------------
   142: ||*---->              do m = 1, 5
   143: |||                      rsd(m,2,j,k) = rsd(m,2,j,k)
   144: |||            >           - dssp * ( + 5.0d+00 * u(m,2,j,k)
   145: |||            >                      - 4.0d+00 * u(m,3,j,k)
   146: |||            >                      +           u(m,4,j,k) )
   147: |||                      rsd(m,3,j,k) = rsd(m,3,j,k)
   148: |||            >           - dssp * ( - 4.0d+00 * u(m,2,j,k)
   149: |||            >                      + 6.0d+00 * u(m,3,j,k)
   150: |||            >                      - 4.0d+00 * u(m,4,j,k)
   151: |||            >                      +           u(m,5,j,k) )
   152: ||*----               end do
   153: ||        
   154: ||V---->              do i = 4, nx - 3
   155: |||*--->                 do m = 1, 5
   156: ||||    F                   rsd(m,i,j,k) = rsd(m,i,j,k)
   157: ||||           >              - dssp * (            u(m,i-2,j,k)
   158: ||||           >                        - 4.0d+00 * u(m,i-1,j,k)
   159: ||||           >                        + 6.0d+00 * u(m,i,j,k)
   160: ||||           >                        - 4.0d+00 * u(m,i+1,j,k)
   161: ||||           >                        +           u(m,i+2,j,k) )
   162: |||*---                  end do
   163: ||V----               end do
   164: ||        
   165: ||        
   166: ||*---->              do m = 1, 5
   167: |||                      rsd(m,nx-2,j,k) = rsd(m,nx-2,j,k)
   168: |||            >           - dssp * (             u(m,nx-4,j,k)
   169: |||            >                      - 4.0d+00 * u(m,nx-3,j,k)
   170: |||            >                      + 6.0d+00 * u(m,nx-2,j,k)
   171: |||            >                      - 4.0d+00 * u(m,nx-1,j,k)  )
   172: |||                      rsd(m,nx-1,j,k) = rsd(m,nx-1,j,k)
   173: |||            >           - dssp * (             u(m,nx-3,j,k)
   174: |||            >                      - 4.0d+00 * u(m,nx-2,j,k)
   175: |||            >                      + 5.0d+00 * u(m,nx-1,j,k) )
   176: ||*----               end do
   177: ||        
   178: |+-----            end do
   179: Y------         end do
   180:           !$omp end do
   181:           !$omp master
   182:                 if (timeron) call timer_stop(t_rhsx)
   183:           
   184:                 if (timeron) call timer_start(t_rhsy)
   185:           !$omp end master
   186:           c---------------------------------------------------------------------
   187:           c   eta-direction flux differences
   188:           c---------------------------------------------------------------------
   189:           !$omp do schedule(static) collapse(2)
   190: Y------>        do k = 2, nz - 1
   191: |+----->           do i = ist, iend
   192: ||V---->              do j = 1, ny
   193: |||                      flux(1,j) = u(3,i,j,k)
   194: |||                      u31 = u(3,i,j,k) * rho_i(i,j,k)
   195: |||       
   196: |||                      q = qs(i,j,k)
   197: |||       
   198: |||                      flux(2,j) = u(2,i,j,k) * u31
   199: |||     F                flux(3,j) = u(3,i,j,k) * u31 + c2 * (u(5,i,j,k)-q)
   200: |||                      flux(4,j) = u(4,i,j,k) * u31
   201: |||     F                flux(5,j) = ( c1 * u(5,i,j,k) - c2 * q ) * u31
   202: ||V----               end do
   203: ||        
   204: ||V---->              do j = jst, jend
   205: |||*--->                 do m = 1, 5
   206: ||||    F                   rsd(m,i,j,k) =  rsd(m,i,j,k)
   207: ||||           >                   - ty2 * ( flux(m,j+1) - flux(m,j-1) )
   208: |||*---                  end do
   209: ||V----               end do
   210: ||        
   211: ||V---->              do j = jst, ny
   212: |||                      tmp = rho_i(i,j,k)
   213: |||       
   214: |||                      u21j = tmp * u(2,i,j,k)
   215: |||                      u31j = tmp * u(3,i,j,k)
   216: |||                      u41j = tmp * u(4,i,j,k)
   217: |||                      u51j = tmp * u(5,i,j,k)
   218: |||       
   219: |||                      tmp = rho_i(i,j-1,k)
   220: |||                      u21jm1 = tmp * u(2,i,j-1,k)
   221: |||                      u31jm1 = tmp * u(3,i,j-1,k)
   222: |||                      u41jm1 = tmp * u(4,i,j-1,k)
   223: |||                      u51jm1 = tmp * u(5,i,j-1,k)
   224: |||       
   225: |||                      flux(2,j) = ty3 * ( u21j - u21jm1 )
   226: |||                      flux(3,j) = (4.0d+00/3.0d+00) * ty3 * (u31j-u31jm1)
   227: |||                      flux(4,j) = ty3 * ( u41j - u41jm1 )
   228: |||     F                flux(5,j) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   229: |||            >              * ty3 * ( ( u21j  **2 + u31j  **2 + u41j  **2 )
   230: |||            >                      - ( u21jm1**2 + u31jm1**2 + u41jm1**2 ) )
   231: |||            >              + (1.0d+00/6.0d+00)
   232: |||            >              * ty3 * ( u31j**2 - u31jm1**2 )
   233: |||            >              + c1 * c5 * ty3 * ( u51j - u51jm1 )
   234: ||V----               end do
   235: ||        
   236: ||V---->              do j = jst, jend
   237: |||       
   238: |||     F                rsd(1,i,j,k) = rsd(1,i,j,k)
   239: |||            >              + dy1 * ty1 * (            u(1,i,j-1,k)
   240: |||            >                             - 2.0d+00 * u(1,i,j,k)
   241: |||            >                             +           u(1,i,j+1,k) )
   242: |||       
   243: |||     F                rsd(2,i,j,k) = rsd(2,i,j,k)
   244: |||            >          + ty3 * c3 * c4 * ( flux(2,j+1) - flux(2,j) )
   245: |||            >              + dy2 * ty1 * (            u(2,i,j-1,k)
   246: |||            >                             - 2.0d+00 * u(2,i,j,k)
   247: |||            >                             +           u(2,i,j+1,k) )
   248: |||       
   249: |||     F                rsd(3,i,j,k) = rsd(3,i,j,k)
   250: |||            >          + ty3 * c3 * c4 * ( flux(3,j+1) - flux(3,j) )
   251: |||            >              + dy3 * ty1 * (            u(3,i,j-1,k)
   252: |||            >                             - 2.0d+00 * u(3,i,j,k)
   253: |||            >                             +           u(3,i,j+1,k) )
   254: |||       
   255: |||     F                rsd(4,i,j,k) = rsd(4,i,j,k)
   256: |||            >          + ty3 * c3 * c4 * ( flux(4,j+1) - flux(4,j) )
   257: |||            >              + dy4 * ty1 * (            u(4,i,j-1,k)
   258: |||            >                             - 2.0d+00 * u(4,i,j,k)
   259: |||            >                             +           u(4,i,j+1,k) )
   260: |||       
   261: |||     F                rsd(5,i,j,k) = rsd(5,i,j,k)
   262: |||            >          + ty3 * c3 * c4 * ( flux(5,j+1) - flux(5,j) )
   263: |||            >              + dy5 * ty1 * (            u(5,i,j-1,k)
   264: |||            >                             - 2.0d+00 * u(5,i,j,k)
   265: |||            >                             +           u(5,i,j+1,k) )
   266: |||       
   267: ||V----               end do
   268: ||        
   269: |+-----            end do
   270: Y------         end do
   271:           !$omp end do
   272:           
   273:           c---------------------------------------------------------------------
   274:           c   fourth-order dissipation
   275:           c---------------------------------------------------------------------
   276:           !$omp do schedule(static) collapse(2)
   277: Y------>        do k = 2, nz - 1
   278: |+----->           do j = jst, jend
   279: ||                    if (j .eq. 2) then
   280: ||V---->                 do i = ist, iend
   281: |||*--->                 do m = 1, 5
   282: ||||    F                   rsd(m,i,2,k) = rsd(m,i,2,k)
   283: ||||           >              - dssp * ( + 5.0d+00 * u(m,i,2,k)
   284: ||||           >                      - 4.0d+00 * u(m,i,3,k)
   285: ||||           >                      +           u(m,i,4,k) )
   286: |||*---                  end do
   287: ||V----                  end do
   288: ||        
   289: ||                    else if (j .eq. 3) then
   290: ||V---->                 do i = ist, iend
   291: |||*--->                 do m = 1, 5
   292: ||||    F                   rsd(m,i,3,k) = rsd(m,i,3,k)
   293: ||||           >              - dssp * ( - 4.0d+00 * u(m,i,2,k)
   294: ||||           >                      + 6.0d+00 * u(m,i,3,k)
   295: ||||           >                      - 4.0d+00 * u(m,i,4,k)
   296: ||||           >                      +           u(m,i,5,k) )
   297: |||*---                  end do
   298: ||V----                  end do
   299: ||        
   300: ||                    else if (j .eq. ny-2) then
   301: ||V---->                 do i = ist, iend
   302: |||*--->                 do m = 1, 5
   303: ||||    F                   rsd(m,i,ny-2,k) = rsd(m,i,ny-2,k)
   304: ||||           >              - dssp * (          u(m,i,ny-4,k)
   305: ||||           >                      - 4.0d+00 * u(m,i,ny-3,k)
   306: ||||           >                      + 6.0d+00 * u(m,i,ny-2,k)
   307: ||||           >                      - 4.0d+00 * u(m,i,ny-1,k)  )
   308: |||*---                  end do
   309: ||V----                  end do
   310: ||        
   311: ||                    else if (j .eq. ny-1) then
   312: ||V---->                 do i = ist, iend
   313: |||*--->                 do m = 1, 5
   314: ||||    F                   rsd(m,i,ny-1,k) = rsd(m,i,ny-1,k)
   315: ||||           >              - dssp * (          u(m,i,ny-3,k)
   316: ||||           >                      - 4.0d+00 * u(m,i,ny-2,k)
   317: ||||           >                      + 5.0d+00 * u(m,i,ny-1,k) )
   318: |||*---                  end do
   319: ||V----                  end do
   320: ||        
   321: ||                    else
   322: ||V---->                 do i = ist, iend
   323: |||*--->                 do m = 1, 5
   324: ||||    F                   rsd(m,i,j,k) = rsd(m,i,j,k)
   325: ||||           >              - dssp * (            u(m,i,j-2,k)
   326: ||||           >                        - 4.0d+00 * u(m,i,j-1,k)
   327: ||||           >                        + 6.0d+00 * u(m,i,j,k)
   328: ||||           >                        - 4.0d+00 * u(m,i,j+1,k)
   329: ||||           >                        +           u(m,i,j+2,k) )
   330: |||*---                  end do
   331: ||V----                  end do
   332: ||                    endif
   333: ||        
   334: |+-----            end do
   335: Y------         end do
   336:           !$omp end do
   337:           !$omp master
   338:                 if (timeron) call timer_stop(t_rhsy)
   339:           
   340:                 if (timeron) call timer_start(t_rhsz)
   341:           !$omp end master
   342:           c---------------------------------------------------------------------
   343:           c   zeta-direction flux differences
   344:           c---------------------------------------------------------------------
   345:           !$omp do schedule(static) collapse(2)
   346: Y------>        do j = jst, jend
   347: |+----->           do i = ist, iend
   348: ||V---->              do k = 1, nz
   349: |||                      utmp(1,k) = u(1,i,j,k)
   350: |||                      utmp(2,k) = u(2,i,j,k)
   351: |||                      utmp(3,k) = u(3,i,j,k)
   352: |||                      utmp(4,k) = u(4,i,j,k)
   353: |||                      utmp(5,k) = u(5,i,j,k)
   354: |||                      utmp(6,k) = rho_i(i,j,k)
   355: |||                   end do
   356: |||                   do k = 1, nz
   357: |||                      flux(1,k) = utmp(4,k)
   358: |||                      u41 = utmp(4,k) * utmp(6,k)
   359: |||       
   360: |||                      q = qs(i,j,k)
   361: |||       
   362: |||                      flux(2,k) = utmp(2,k) * u41
   363: |||                      flux(3,k) = utmp(3,k) * u41
   364: |||     F                flux(4,k) = utmp(4,k) * u41 + c2 * (utmp(5,k)-q)
   365: |||     F                flux(5,k) = ( c1 * utmp(5,k) - c2 * q ) * u41
   366: ||V----               end do
   367: ||        
   368: ||V---->              do k = 2, nz - 1
   369: |||*--->                 do m = 1, 5
   370: ||||    F                   rtmp(m,k) =  rsd(m,i,j,k)
   371: ||||           >                - tz2 * ( flux(m,k+1) - flux(m,k-1) )
   372: |||*---                  end do
   373: ||V----               end do
   374: ||        
   375: ||V---->              do k = 2, nz
   376: |||                      tmp = utmp(6,k)
   377: |||       
   378: |||                      u21k = tmp * utmp(2,k)
   379: |||                      u31k = tmp * utmp(3,k)
   380: |||                      u41k = tmp * utmp(4,k)
   381: |||                      u51k = tmp * utmp(5,k)
   382: |||       
   383: |||                      tmp = utmp(6,k-1)
   384: |||       
   385: |||                      u21km1 = tmp * utmp(2,k-1)
   386: |||                      u31km1 = tmp * utmp(3,k-1)
   387: |||                      u41km1 = tmp * utmp(4,k-1)
   388: |||                      u51km1 = tmp * utmp(5,k-1)
   389: |||       
   390: |||                      flux(2,k) = tz3 * ( u21k - u21km1 )
   391: |||                      flux(3,k) = tz3 * ( u31k - u31km1 )
   392: |||                      flux(4,k) = (4.0d+00/3.0d+00) * tz3 * (u41k-u41km1)
   393: |||     F                flux(5,k) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   394: |||            >              * tz3 * ( ( u21k  **2 + u31k  **2 + u41k  **2 )
   395: |||            >                      - ( u21km1**2 + u31km1**2 + u41km1**2 ) )
   396: |||            >              + (1.0d+00/6.0d+00)
   397: |||            >              * tz3 * ( u41k**2 - u41km1**2 )
   398: |||            >              + c1 * c5 * tz3 * ( u51k - u51km1 )
   399: ||V----               end do
   400: ||        
   401: ||V---->              do k = 2, nz - 1
   402: |||     F                rtmp(1,k) = rtmp(1,k)
   403: |||            >              + dz1 * tz1 * (            utmp(1,k-1)
   404: |||            >                             - 2.0d+00 * utmp(1,k)
   405: |||            >                             +           utmp(1,k+1) )
   406: |||     F                rtmp(2,k) = rtmp(2,k)
   407: |||            >          + tz3 * c3 * c4 * ( flux(2,k+1) - flux(2,k) )
   408: |||            >              + dz2 * tz1 * (            utmp(2,k-1)
   409: |||            >                             - 2.0d+00 * utmp(2,k)
   410: |||            >                             +           utmp(2,k+1) )
   411: |||     F                rtmp(3,k) = rtmp(3,k)
   412: |||            >          + tz3 * c3 * c4 * ( flux(3,k+1) - flux(3,k) )
   413: |||            >              + dz3 * tz1 * (            utmp(3,k-1)
   414: |||            >                             - 2.0d+00 * utmp(3,k)
   415: |||            >                             +           utmp(3,k+1) )
   416: |||     F                rtmp(4,k) = rtmp(4,k)
   417: |||            >          + tz3 * c3 * c4 * ( flux(4,k+1) - flux(4,k) )
   418: |||            >              + dz4 * tz1 * (            utmp(4,k-1)
   419: |||            >                             - 2.0d+00 * utmp(4,k)
   420: |||            >                             +           utmp(4,k+1) )
   421: |||     F                rtmp(5,k) = rtmp(5,k)
   422: |||            >          + tz3 * c3 * c4 * ( flux(5,k+1) - flux(5,k) )
   423: |||            >              + dz5 * tz1 * (            utmp(5,k-1)
   424: |||            >                             - 2.0d+00 * utmp(5,k)
   425: |||            >                             +           utmp(5,k+1) )
   426: ||V----               end do
   427: ||        
   428: ||        c---------------------------------------------------------------------
   429: ||        c   fourth-order dissipation
   430: ||        c---------------------------------------------------------------------
   431: ||*---->              do m = 1, 5
   432: |||                      rsd(m,i,j,2) = rtmp(m,2)
   433: |||            >           - dssp * ( + 5.0d+00 * utmp(m,2)
   434: |||            >                      - 4.0d+00 * utmp(m,3)
   435: |||            >                      +           utmp(m,4) )
   436: |||                      rsd(m,i,j,3) = rtmp(m,3)
   437: |||            >           - dssp * ( - 4.0d+00 * utmp(m,2)
   438: |||            >                      + 6.0d+00 * utmp(m,3)
   439: |||            >                      - 4.0d+00 * utmp(m,4)
   440: |||            >                      +           utmp(m,5) )
   441: ||*----               end do
   442: ||        
   443: ||V---->              do k = 4, nz - 3
   444: |||*--->                 do m = 1, 5
   445: ||||    F                   rsd(m,i,j,k) = rtmp(m,k)
   446: ||||           >              - dssp * (            utmp(m,k-2)
   447: ||||           >                        - 4.0d+00 * utmp(m,k-1)
   448: ||||           >                        + 6.0d+00 * utmp(m,k)
   449: ||||           >                        - 4.0d+00 * utmp(m,k+1)
   450: ||||           >                        +           utmp(m,k+2) )
   451: |||*---                  end do
   452: ||V----               end do
   453: ||        
   454: ||*---->              do m = 1, 5
   455: |||                      rsd(m,i,j,nz-2) = rtmp(m,nz-2)
   456: |||            >           - dssp * (             utmp(m,nz-4)
   457: |||            >                      - 4.0d+00 * utmp(m,nz-3)
   458: |||            >                      + 6.0d+00 * utmp(m,nz-2)
   459: |||            >                      - 4.0d+00 * utmp(m,nz-1)  )
   460: |||                      rsd(m,i,j,nz-1) = rtmp(m,nz-1)
   461: |||            >           - dssp * (             utmp(m,nz-3)
   462: |||            >                      - 4.0d+00 * utmp(m,nz-2)
   463: |||            >                      + 5.0d+00 * utmp(m,nz-1) )
   464: ||*----               end do
   465: |+-----            end do
   466: Y------         end do
   467:           !$omp end do nowait
   468:           !$omp master
   469:                 if (timeron) call timer_stop(t_rhsz)
   470:           !$omp end master
   471:           !$omp end parallel
   472:           
   473:                 if (timeron) call timer_stop(t_rhs)
   474:           
   475:                 return
   476:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: rhs.f

PROCEDURE NAME: RHS
INLINE LIST

  ROOT: RHS (rhs.f:5)
  -> NOINLINE: TIMER_START (rhs.f:32)
     *** Source for routine not found.
  -> NOINLINE: TIMER_START (rhs.f:57)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (rhs.f:182)
     *** Source for routine not found.
  -> NOINLINE: TIMER_START (rhs.f:184)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (rhs.f:338)
     *** Source for routine not found.
  -> NOINLINE: TIMER_START (rhs.f:340)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (rhs.f:469)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (rhs.f:473)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: rhs.f

PROCEDURE NAME: RHS
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (rhs.f:39)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (rhs.f:39)
    *** The number of VLOAD, VSTORE. :  0,  2. (rhs.f:39)

    LOOP BEGIN: (rhs.f:41)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:41)
      *** The number of VLOAD, VSTORE. :  9,  7. (rhs.f:41)
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:64)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (rhs.f:64)
    *** The number of VLOAD, VSTORE. :  0,  2. (rhs.f:64)

    LOOP BEGIN: (rhs.f:66)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:66)
      *** The number of VLOAD, VSTORE. :  6,  5. (rhs.f:66)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:66)
    LOOP END

    LOOP BEGIN: (rhs.f:79)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:79)
      *** The number of VLOAD, VSTORE. : 15,  5. (rhs.f:79)
    LOOP END

    LOOP BEGIN: (rhs.f:86)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:86)
      *** The number of VLOAD, VSTORE. : 10,  4. (rhs.f:86)
    LOOP END

    LOOP BEGIN: (rhs.f:112)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:112)
      *** The number of VLOAD, VSTORE. : 28,  5. (rhs.f:112)
    LOOP END

    LOOP BEGIN: (rhs.f:154)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:154)
      *** The number of VLOAD, VSTORE. : 30,  5. (rhs.f:154)
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:190)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (rhs.f:190)
    *** The number of VLOAD, VSTORE. :  0,  2. (rhs.f:190)

    LOOP BEGIN: (rhs.f:192)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:192)
      *** The number of VLOAD, VSTORE. :  6,  5. (rhs.f:192)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:192)
    LOOP END

    LOOP BEGIN: (rhs.f:204)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:204)
      *** The number of VLOAD, VSTORE. : 15,  5. (rhs.f:204)
    LOOP END

    LOOP BEGIN: (rhs.f:211)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:211)
      *** The number of VLOAD, VSTORE. : 10,  4. (rhs.f:211)
    LOOP END

    LOOP BEGIN: (rhs.f:236)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:236)
      *** The number of VLOAD, VSTORE. : 28,  5. (rhs.f:236)
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:277)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (rhs.f:277)
    *** The number of VLOAD, VSTORE. :  0,  2. (rhs.f:277)

    LOOP BEGIN: (rhs.f:322)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:322)
      *** The number of VLOAD, VSTORE. : 30,  5. (rhs.f:322)
    LOOP END

    LOOP BEGIN: (rhs.f:312)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:312)
      *** The number of VLOAD, VSTORE. : 20,  5. (rhs.f:312)
    LOOP END

    LOOP BEGIN: (rhs.f:301)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:301)
      *** The number of VLOAD, VSTORE. : 25,  5. (rhs.f:301)
    LOOP END

    LOOP BEGIN: (rhs.f:290)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:290)
      *** The number of VLOAD, VSTORE. : 25,  5. (rhs.f:290)
    LOOP END

    LOOP BEGIN: (rhs.f:280)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:280)
      *** The number of VLOAD, VSTORE. : 20,  5. (rhs.f:280)
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:346)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (rhs.f:346)
    *** The number of VLOAD, VSTORE. :  0,  2. (rhs.f:346)

    LOOP BEGIN: (rhs.f:348)
      <Vectorized loop.>
      **  Fused loop. (rhs.f:348)
      **  Runtime-tested short loop. (rhs.f:348)
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:348)
      *** The number of VLOAD, VSTORE. :  7, 11. (rhs.f:348)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:348)
    LOOP END

    LOOP BEGIN: (rhs.f:368)
      <Vectorized loop.>
      **  Runtime-tested short loop. (rhs.f:368)
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:368)
      *** The number of VLOAD, VSTORE. : 15,  5. (rhs.f:368)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:368)
    LOOP END

    LOOP BEGIN: (rhs.f:375)
      <Vectorized loop.>
      **  Runtime-tested short loop. (rhs.f:375)
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:375)
      *** The number of VLOAD, VSTORE. : 10,  4. (rhs.f:375)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:375)
    LOOP END

    LOOP BEGIN: (rhs.f:401)
      <Vectorized loop.>
      **  Runtime-tested short loop. (rhs.f:401)
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:401)
      *** The number of VLOAD, VSTORE. : 28,  5. (rhs.f:401)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:401)
    LOOP END

    LOOP BEGIN: (rhs.f:443)
      <Vectorized loop.>
      **  Runtime-tested short loop. (rhs.f:443)
      *** The number of VGT,   VSC.    :  0,  0. (rhs.f:443)
      *** The number of VLOAD, VSTORE. : 30,  5. (rhs.f:443)
      *** Iteration count is assumed. Iteration count. : 162 (rhs.f:443)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: rhs.f

PROCEDURE NAME: RHS
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
    Total size                 : 18960 bytes
      Register spill area      :  2240 bytes
      Parameter area           :    16 bytes
      Register save area       :   176 bytes
      User data area           : 16528 bytes
      Others                   :     0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (rhs.f:39)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (rhs.f:39)
    *** Estimated execution cycle                       : 60
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 12
            Over basic blocks                           : 12
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (rhs.f:41)
      *** Estimated execution cycle                     : 979
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:64)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (rhs.f:64)
    *** Estimated execution cycle                       : 508
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 51
            Over basic blocks                           : 51
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 67
            Over basic blocks                           : 67
    *** The number of SCALAR REGISTER TRANSFER          : 100

    LOOP BEGIN: (rhs.f:79)
      *** Estimated execution cycle                     : 813
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 18
    LOOP END

    LOOP BEGIN: (rhs.f:86)
      *** Estimated execution cycle                     : 1199
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER TRANSFER        : 19
    LOOP END

    LOOP BEGIN: (rhs.f:112)
      *** Estimated execution cycle                     : 1510
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 16
              Over basic blocks                         : 16
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (rhs.f:154)
      *** Estimated execution cycle                     : 1287
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 3
              Over basic blocks                         : 3
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 16
              Over basic blocks                         : 16
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:190)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (rhs.f:190)
    *** Estimated execution cycle                       : 516
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 12
            Over basic blocks                           : 12
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 35
            Over basic blocks                           : 35
    *** The number of SCALAR REGISTER TRANSFER          : 38

    LOOP BEGIN: (rhs.f:204)
      *** Estimated execution cycle                     : 813
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (rhs.f:211)
      *** Estimated execution cycle                     : 1168
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER TRANSFER        : 8
    LOOP END

    LOOP BEGIN: (rhs.f:236)
      *** Estimated execution cycle                     : 1510
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 22
              Over basic blocks                         : 22
      *** The number of SCALAR REGISTER TRANSFER        : 5
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:277)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (rhs.f:277)
    *** Estimated execution cycle                       : 121
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 63
            Over basic blocks                           : 63
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 15
            Over basic blocks                           : 15
    *** The number of SCALAR REGISTER TRANSFER          : 90

    LOOP BEGIN: (rhs.f:322)
      *** Estimated execution cycle                     : 1287
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 18
              Over basic blocks                         : 18
      *** The number of SCALAR REGISTER TRANSFER        : 1
    LOOP END

    LOOP BEGIN: (rhs.f:312)
      *** Estimated execution cycle                     : 1127
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 10
              Over basic blocks                         : 10
      *** The number of SCALAR REGISTER TRANSFER        : 9
    LOOP END

    LOOP BEGIN: (rhs.f:301)
      *** Estimated execution cycle                     : 1190
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 15
              Over basic blocks                         : 15
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END

    LOOP BEGIN: (rhs.f:290)
      *** Estimated execution cycle                     : 1389
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 17
              Over basic blocks                         : 17
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END

    LOOP BEGIN: (rhs.f:280)
      *** Estimated execution cycle                     : 1291
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 9
              Over basic blocks                         : 9
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 10
              Over basic blocks                         : 10
      *** The number of SCALAR REGISTER TRANSFER        : 7
    LOOP END
  LOOP END

  LOOP BEGIN: (rhs.f:346)
    *** Estimated execution cycle                       : 575
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
  LOOP END

  LOOP BEGIN: (rhs.f:346)
    *** Estimated execution cycle                       : 1299
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 30
            Not enough registers                        :  9
            Over basic blocks                           : 21
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 73
            Over basic blocks                           : 73
    *** The number of SCALAR REGISTER TRANSFER          : 23
  LOOP END


