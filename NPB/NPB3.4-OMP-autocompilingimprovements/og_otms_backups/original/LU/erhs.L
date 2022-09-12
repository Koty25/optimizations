NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME : erhs.f

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
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
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
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
        -fmove-loop-invariants-unsafe            : disable
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
        -freplace-loop-equation                  : disable
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
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: erhs.f

PROCEDURE NAME: ERHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    34: par(1801): Parallel routine generated.: ERHS$1
    39: par(1803): Parallelized by "do".
    39: opt(2495): Loop nest collapsed by "collapse" clause.
    39: vec( 101): Vectorized loop.
    41: vec( 101): Vectorized loop.
    42: opt(1418): Constant-length loop is expanded.
    51: par(1803): Parallelized by "do".
    51: opt(2495): Loop nest collapsed by "collapse" clause.
    51: vec( 101): Vectorized loop.
    55: vec( 101): Vectorized loop.
    57: opt(1418): Constant-length loop is expanded.
    58: vec( 128): Fused multiply-add operation applied.
    82: opt(3008): Reference within a conditional branch moved outside loop.
    82: par(1803): Parallelized by "do".
    82: opt(2495): Loop nest collapsed by "collapse" clause.
    82: vec( 102): Partially vectorized loop.
    84: vec( 101): Vectorized loop.
    84: vec( 124): Iteration count is assumed. Iteration count.: 162
    87: vec( 128): Fused multiply-add operation applied.
    91: vec( 128): Fused multiply-add operation applied.
    95: vec( 128): Fused multiply-add operation applied.
    98: vec( 101): Vectorized loop.
    99: opt(1418): Constant-length loop is expanded.
   100: vec( 128): Fused multiply-add operation applied.
   104: vec( 101): Vectorized loop.
   123: vec( 128): Fused multiply-add operation applied.
   131: vec( 101): Vectorized loop.
   132: vec( 128): Fused multiply-add operation applied.
   136: vec( 128): Fused multiply-add operation applied.
   141: vec( 128): Fused multiply-add operation applied.
   146: vec( 128): Fused multiply-add operation applied.
   151: vec( 128): Fused multiply-add operation applied.
   161: opt(1418): Constant-length loop is expanded.
   173: vec( 101): Vectorized loop.
   174: opt(1418): Constant-length loop is expanded.
   175: vec( 128): Fused multiply-add operation applied.
   184: opt(1418): Constant-length loop is expanded.
   205: opt(3008): Reference within a conditional branch moved outside loop.
   205: par(1803): Parallelized by "do".
   205: opt(2495): Loop nest collapsed by "collapse" clause.
   205: vec( 102): Partially vectorized loop.
   207: vec( 101): Vectorized loop.
   207: vec( 124): Iteration count is assumed. Iteration count.: 162
   210: vec( 128): Fused multiply-add operation applied.
   215: vec( 128): Fused multiply-add operation applied.
   218: vec( 128): Fused multiply-add operation applied.
   221: vec( 101): Vectorized loop.
   222: opt(1418): Constant-length loop is expanded.
   223: vec( 128): Fused multiply-add operation applied.
   228: vec( 101): Vectorized loop.
   247: vec( 128): Fused multiply-add operation applied.
   255: vec( 101): Vectorized loop.
   256: vec( 128): Fused multiply-add operation applied.
   260: vec( 128): Fused multiply-add operation applied.
   265: vec( 128): Fused multiply-add operation applied.
   270: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   285: opt(1418): Constant-length loop is expanded.
   297: vec( 101): Vectorized loop.
   298: opt(1418): Constant-length loop is expanded.
   299: vec( 128): Fused multiply-add operation applied.
   308: opt(1418): Constant-length loop is expanded.
   328: opt(3008): Reference within a conditional branch moved outside loop.
   328: par(1803): Parallelized by "do".
   328: opt(2495): Loop nest collapsed by "collapse" clause.
   328: vec( 102): Partially vectorized loop.
   330: vec( 101): Vectorized loop.
   330: vec( 124): Iteration count is assumed. Iteration count.: 162
   333: vec( 128): Fused multiply-add operation applied.
   339: vec( 128): Fused multiply-add operation applied.
   341: vec( 128): Fused multiply-add operation applied.
   344: vec( 101): Vectorized loop.
   344: vec( 124): Iteration count is assumed. Iteration count.: 162
   345: opt(1418): Constant-length loop is expanded.
   346: vec( 128): Fused multiply-add operation applied.
   351: vec( 101): Vectorized loop.
   351: vec( 124): Iteration count is assumed. Iteration count.: 162
   370: vec( 128): Fused multiply-add operation applied.
   378: vec( 101): Vectorized loop.
   378: vec( 124): Iteration count is assumed. Iteration count.: 162
   379: vec( 128): Fused multiply-add operation applied.
   383: vec( 128): Fused multiply-add operation applied.
   388: vec( 128): Fused multiply-add operation applied.
   393: vec( 128): Fused multiply-add operation applied.
   398: vec( 128): Fused multiply-add operation applied.
   408: opt(1418): Constant-length loop is expanded.
   420: vec( 101): Vectorized loop.
   421: opt(1418): Constant-length loop is expanded.
   422: vec( 128): Fused multiply-add operation applied.
   431: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: erhs.f

PROCEDURE NAME: ERHS
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

     4:                 subroutine erhs
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c
    11:           c   compute the right hand side based on exact solution
    12:           c
    13:           c---------------------------------------------------------------------
    14:           
    15:                 use lu_data
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  local variables
    20:           c---------------------------------------------------------------------
    21:                 integer i, j, k, m
    22:                 double precision  xi, eta, zeta
    23:                 double precision  q
    24:                 double precision  u21, u31, u41
    25:                 double precision  tmp
    26:                 double precision  u21i, u31i, u41i, u51i
    27:                 double precision  u21j, u31j, u41j, u51j
    28:                 double precision  u21k, u31k, u41k, u51k
    29:                 double precision  u21im1, u31im1, u41im1, u51im1
    30:                 double precision  u21jm1, u31jm1, u41jm1, u51jm1
    31:                 double precision  u21km1, u31km1, u41km1, u51km1
    32:           
    33:           
    34:           !$omp parallel default(shared) private(i,j,k,m,xi,eta,zeta,tmp,q,
    35:           !$omp&   u51im1,u41im1,u31im1,u21im1,u51i,u41i,u31i,u21i,u21,
    36:           !$omp&   u51jm1,u41jm1,u31jm1,u21jm1,u51j,u41j,u31j,u21j,u31,
    37:           !$omp&   u51km1,u41km1,u31km1,u21km1,u51k,u41k,u31k,u21k,u41)
    38:           !$omp do schedule(static) collapse(2)
    39: Y------>        do k = 1, nz
    40: |+----->           do j = 1, ny
    41: ||V---->              do i = 1, nx
    42: |||*--->                 do m = 1, 5
    43: ||||                        frct( m, i, j, k ) = 0.0d+00
    44: |||*---                  end do
    45: ||V----               end do
    46: |+-----            end do
    47: Y------         end do
    48:           !$omp end do nowait
    49:           
    50:           !$omp do schedule(static) collapse(2)
    51: Y------>        do k = 1, nz
    52: |+----->           do j = 1, ny
    53: ||                    zeta = ( dble(k-1) ) / ( nz - 1 )
    54: ||                    eta = ( dble(j-1) ) / ( ny0 - 1 )
    55: ||V---->              do i = 1, nx
    56: |||                      xi = ( dble(i-1) ) / ( nx0 - 1 )
    57: |||*--->                 do m = 1, 5
    58: ||||    F                   rsd(m,i,j,k) =  ce(m,1)
    59: ||||           >                 + (ce(m,2)
    60: ||||           >                 + (ce(m,5)
    61: ||||           >                 + (ce(m,8)
    62: ||||           >                 +  ce(m,11) * xi) * xi) * xi) * xi
    63: ||||           >                 + (ce(m,3)
    64: ||||           >                 + (ce(m,6)
    65: ||||           >                 + (ce(m,9)
    66: ||||           >                 +  ce(m,12) * eta) * eta) * eta) * eta
    67: ||||           >                 + (ce(m,4)
    68: ||||           >                 + (ce(m,7)
    69: ||||           >                 + (ce(m,10)
    70: ||||           >                 +  ce(m,13) * zeta) * zeta) * zeta) * zeta
    71: |||*---                  end do
    72: ||V----               end do
    73: |+-----            end do
    74: Y------         end do
    75:           !$omp end do
    76:           
    77:           c---------------------------------------------------------------------
    78:           c   xi-direction flux differences
    79:           c---------------------------------------------------------------------
    80:           
    81:           !$omp do schedule(static) collapse(2)
    82: Y------>        do k = 2, nz - 1
    83: |+----->           do j = jst, jend
    84: ||V---->              do i = 1, nx
    85: |||                      flux(1,i) = rsd(2,i,j,k)
    86: |||                      u21 = rsd(2,i,j,k) / rsd(1,i,j,k)
    87: |||     F                q = 0.50d+00 * (  rsd(2,i,j,k) * rsd(2,i,j,k)
    88: |||            >                         + rsd(3,i,j,k) * rsd(3,i,j,k)
    89: |||            >                         + rsd(4,i,j,k) * rsd(4,i,j,k) )
    90: |||            >                      / rsd(1,i,j,k)
    91: |||     F                flux(2,i) = rsd(2,i,j,k) * u21 + c2 *
    92: |||            >                         ( rsd(5,i,j,k) - q )
    93: |||                      flux(3,i) = rsd(3,i,j,k) * u21
    94: |||                      flux(4,i) = rsd(4,i,j,k) * u21
    95: |||     F                flux(5,i) = ( c1 * rsd(5,i,j,k) - c2 * q ) * u21
    96: ||V----               end do
    97: ||        
    98: ||V---->              do i = ist, iend
    99: |||*--->                 do m = 1, 5
   100: ||||    F                   frct(m,i,j,k) =  frct(m,i,j,k)
   101: ||||           >                   - tx2 * ( flux(m,i+1) - flux(m,i-1) )
   102: |||*---                  end do
   103: ||V----               end do
   104: ||V---->              do i = ist, nx
   105: |||                      tmp = 1.0d+00 / rsd(1,i,j,k)
   106: |||       
   107: |||                      u21i = tmp * rsd(2,i,j,k)
   108: |||                      u31i = tmp * rsd(3,i,j,k)
   109: |||                      u41i = tmp * rsd(4,i,j,k)
   110: |||                      u51i = tmp * rsd(5,i,j,k)
   111: |||       
   112: |||                      tmp = 1.0d+00 / rsd(1,i-1,j,k)
   113: |||       
   114: |||                      u21im1 = tmp * rsd(2,i-1,j,k)
   115: |||                      u31im1 = tmp * rsd(3,i-1,j,k)
   116: |||                      u41im1 = tmp * rsd(4,i-1,j,k)
   117: |||                      u51im1 = tmp * rsd(5,i-1,j,k)
   118: |||       
   119: |||                      flux(2,i) = (4.0d+00/3.0d+00) * tx3 *
   120: |||            >                        ( u21i - u21im1 )
   121: |||                      flux(3,i) = tx3 * ( u31i - u31im1 )
   122: |||                      flux(4,i) = tx3 * ( u41i - u41im1 )
   123: |||     F                flux(5,i) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   124: |||            >              * tx3 * ( ( u21i  **2 + u31i  **2 + u41i  **2 )
   125: |||            >                      - ( u21im1**2 + u31im1**2 + u41im1**2 ) )
   126: |||            >              + (1.0d+00/6.0d+00)
   127: |||            >              * tx3 * ( u21i**2 - u21im1**2 )
   128: |||            >              + c1 * c5 * tx3 * ( u51i - u51im1 )
   129: ||V----               end do
   130: ||        
   131: ||V---->              do i = ist, iend
   132: |||     F                frct(1,i,j,k) = frct(1,i,j,k)
   133: |||            >              + dx1 * tx1 * (            rsd(1,i-1,j,k)
   134: |||            >                             - 2.0d+00 * rsd(1,i,j,k)
   135: |||            >                             +           rsd(1,i+1,j,k) )
   136: |||     F                frct(2,i,j,k) = frct(2,i,j,k)
   137: |||            >           + tx3 * c3 * c4 * ( flux(2,i+1) - flux(2,i) )
   138: |||            >              + dx2 * tx1 * (            rsd(2,i-1,j,k)
   139: |||            >                             - 2.0d+00 * rsd(2,i,j,k)
   140: |||            >                             +           rsd(2,i+1,j,k) )
   141: |||     F                frct(3,i,j,k) = frct(3,i,j,k)
   142: |||            >           + tx3 * c3 * c4 * ( flux(3,i+1) - flux(3,i) )
   143: |||            >              + dx3 * tx1 * (            rsd(3,i-1,j,k)
   144: |||            >                             - 2.0d+00 * rsd(3,i,j,k)
   145: |||            >                             +           rsd(3,i+1,j,k) )
   146: |||     F                frct(4,i,j,k) = frct(4,i,j,k)
   147: |||            >            + tx3 * c3 * c4 * ( flux(4,i+1) - flux(4,i) )
   148: |||            >              + dx4 * tx1 * (            rsd(4,i-1,j,k)
   149: |||            >                             - 2.0d+00 * rsd(4,i,j,k)
   150: |||            >                             +           rsd(4,i+1,j,k) )
   151: |||     F                frct(5,i,j,k) = frct(5,i,j,k)
   152: |||            >           + tx3 * c3 * c4 * ( flux(5,i+1) - flux(5,i) )
   153: |||            >              + dx5 * tx1 * (            rsd(5,i-1,j,k)
   154: |||            >                             - 2.0d+00 * rsd(5,i,j,k)
   155: |||            >                             +           rsd(5,i+1,j,k) )
   156: ||V----               end do
   157: ||        
   158: ||        c---------------------------------------------------------------------
   159: ||        c   Fourth-order dissipation
   160: ||        c---------------------------------------------------------------------
   161: ||*---->              do m = 1, 5
   162: |||                      frct(m,2,j,k) = frct(m,2,j,k)
   163: |||            >           - dssp * ( + 5.0d+00 * rsd(m,2,j,k)
   164: |||            >                       - 4.0d+00 * rsd(m,3,j,k)
   165: |||            >                       +           rsd(m,4,j,k) )
   166: |||                      frct(m,3,j,k) = frct(m,3,j,k)
   167: |||            >           - dssp * ( - 4.0d+00 * rsd(m,2,j,k)
   168: |||            >                       + 6.0d+00 * rsd(m,3,j,k)
   169: |||            >                       - 4.0d+00 * rsd(m,4,j,k)
   170: |||            >                       +           rsd(m,5,j,k) )
   171: ||*----               end do
   172: ||        
   173: ||V---->              do i = 4, nx - 3
   174: |||*--->                 do m = 1, 5
   175: ||||    F                   frct(m,i,j,k) = frct(m,i,j,k)
   176: ||||           >              - dssp * (            rsd(m,i-2,j,k)
   177: ||||           >                         - 4.0d+00 * rsd(m,i-1,j,k)
   178: ||||           >                         + 6.0d+00 * rsd(m,i,j,k)
   179: ||||           >                         - 4.0d+00 * rsd(m,i+1,j,k)
   180: ||||           >                         +           rsd(m,i+2,j,k) )
   181: |||*---                  end do
   182: ||V----               end do
   183: ||        
   184: ||*---->              do m = 1, 5
   185: |||                      frct(m,nx-2,j,k) = frct(m,nx-2,j,k)
   186: |||            >           - dssp * (             rsd(m,nx-4,j,k)
   187: |||            >                       - 4.0d+00 * rsd(m,nx-3,j,k)
   188: |||            >                       + 6.0d+00 * rsd(m,nx-2,j,k)
   189: |||            >                       - 4.0d+00 * rsd(m,nx-1,j,k)  )
   190: |||                      frct(m,nx-1,j,k) = frct(m,nx-1,j,k)
   191: |||            >           - dssp * (             rsd(m,nx-3,j,k)
   192: |||            >                       - 4.0d+00 * rsd(m,nx-2,j,k)
   193: |||            >                       + 5.0d+00 * rsd(m,nx-1,j,k) )
   194: ||*----               end do
   195: ||        
   196: |+-----            end do
   197: Y------         end do
   198:           !$omp end do
   199:           
   200:           c---------------------------------------------------------------------
   201:           c   eta-direction flux differences
   202:           c---------------------------------------------------------------------
   203:           
   204:           !$omp do schedule(static) collapse(2)
   205: Y------>        do k = 2, nz - 1
   206: |+----->           do i = ist, iend
   207: ||V---->              do j = 1, ny
   208: |||                      flux(1,j) = rsd(3,i,j,k)
   209: |||                      u31 = rsd(3,i,j,k) / rsd(1,i,j,k)
   210: |||     F                q = 0.50d+00 * (  rsd(2,i,j,k) * rsd(2,i,j,k)
   211: |||            >                         + rsd(3,i,j,k) * rsd(3,i,j,k)
   212: |||            >                         + rsd(4,i,j,k) * rsd(4,i,j,k) )
   213: |||            >                      / rsd(1,i,j,k)
   214: |||                      flux(2,j) = rsd(2,i,j,k) * u31
   215: |||     F                flux(3,j) = rsd(3,i,j,k) * u31 + c2 *
   216: |||            >                       ( rsd(5,i,j,k) - q )
   217: |||                      flux(4,j) = rsd(4,i,j,k) * u31
   218: |||     F                flux(5,j) = ( c1 * rsd(5,i,j,k) - c2 * q ) * u31
   219: ||V----               end do
   220: ||        
   221: ||V---->              do j = jst, jend
   222: |||*--->                 do m = 1, 5
   223: ||||    F                   frct(m,i,j,k) =  frct(m,i,j,k)
   224: ||||           >                 - ty2 * ( flux(m,j+1) - flux(m,j-1) )
   225: |||*---                  end do
   226: ||V----               end do
   227: ||        
   228: ||V---->              do j = jst, ny
   229: |||                      tmp = 1.0d+00 / rsd(1,i,j,k)
   230: |||       
   231: |||                      u21j = tmp * rsd(2,i,j,k)
   232: |||                      u31j = tmp * rsd(3,i,j,k)
   233: |||                      u41j = tmp * rsd(4,i,j,k)
   234: |||                      u51j = tmp * rsd(5,i,j,k)
   235: |||       
   236: |||                      tmp = 1.0d+00 / rsd(1,i,j-1,k)
   237: |||       
   238: |||                      u21jm1 = tmp * rsd(2,i,j-1,k)
   239: |||                      u31jm1 = tmp * rsd(3,i,j-1,k)
   240: |||                      u41jm1 = tmp * rsd(4,i,j-1,k)
   241: |||                      u51jm1 = tmp * rsd(5,i,j-1,k)
   242: |||       
   243: |||                      flux(2,j) = ty3 * ( u21j - u21jm1 )
   244: |||                      flux(3,j) = (4.0d+00/3.0d+00) * ty3 *
   245: |||            >                       ( u31j - u31jm1 )
   246: |||                      flux(4,j) = ty3 * ( u41j - u41jm1 )
   247: |||     F                flux(5,j) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   248: |||            >              * ty3 * ( ( u21j  **2 + u31j  **2 + u41j  **2 )
   249: |||            >                      - ( u21jm1**2 + u31jm1**2 + u41jm1**2 ) )
   250: |||            >              + (1.0d+00/6.0d+00)
   251: |||            >              * ty3 * ( u31j**2 - u31jm1**2 )
   252: |||            >              + c1 * c5 * ty3 * ( u51j - u51jm1 )
   253: ||V----               end do
   254: ||        
   255: ||V---->              do j = jst, jend
   256: |||     F                frct(1,i,j,k) = frct(1,i,j,k)
   257: |||            >              + dy1 * ty1 * (            rsd(1,i,j-1,k)
   258: |||            >                             - 2.0d+00 * rsd(1,i,j,k)
   259: |||            >                             +           rsd(1,i,j+1,k) )
   260: |||     F                frct(2,i,j,k) = frct(2,i,j,k)
   261: |||            >          + ty3 * c3 * c4 * ( flux(2,j+1) - flux(2,j) )
   262: |||            >              + dy2 * ty1 * (            rsd(2,i,j-1,k)
   263: |||            >                             - 2.0d+00 * rsd(2,i,j,k)
   264: |||            >                             +           rsd(2,i,j+1,k) )
   265: |||     F                frct(3,i,j,k) = frct(3,i,j,k)
   266: |||            >          + ty3 * c3 * c4 * ( flux(3,j+1) - flux(3,j) )
   267: |||            >              + dy3 * ty1 * (            rsd(3,i,j-1,k)
   268: |||            >                             - 2.0d+00 * rsd(3,i,j,k)
   269: |||            >                             +           rsd(3,i,j+1,k) )
   270: |||     F                frct(4,i,j,k) = frct(4,i,j,k)
   271: |||            >          + ty3 * c3 * c4 * ( flux(4,j+1) - flux(4,j) )
   272: |||            >              + dy4 * ty1 * (            rsd(4,i,j-1,k)
   273: |||            >                             - 2.0d+00 * rsd(4,i,j,k)
   274: |||            >                             +           rsd(4,i,j+1,k) )
   275: |||     F                frct(5,i,j,k) = frct(5,i,j,k)
   276: |||            >          + ty3 * c3 * c4 * ( flux(5,j+1) - flux(5,j) )
   277: |||            >              + dy5 * ty1 * (            rsd(5,i,j-1,k)
   278: |||            >                             - 2.0d+00 * rsd(5,i,j,k)
   279: |||            >                             +           rsd(5,i,j+1,k) )
   280: ||V----               end do
   281: ||        
   282: ||        c---------------------------------------------------------------------
   283: ||        c   fourth-order dissipation
   284: ||        c---------------------------------------------------------------------
   285: ||*---->              do m = 1, 5
   286: |||                      frct(m,i,2,k) = frct(m,i,2,k)
   287: |||            >           - dssp * ( + 5.0d+00 * rsd(m,i,2,k)
   288: |||            >                       - 4.0d+00 * rsd(m,i,3,k)
   289: |||            >                       +           rsd(m,i,4,k) )
   290: |||                      frct(m,i,3,k) = frct(m,i,3,k)
   291: |||            >           - dssp * ( - 4.0d+00 * rsd(m,i,2,k)
   292: |||            >                       + 6.0d+00 * rsd(m,i,3,k)
   293: |||            >                       - 4.0d+00 * rsd(m,i,4,k)
   294: |||            >                       +           rsd(m,i,5,k) )
   295: ||*----               end do
   296: ||        
   297: ||V---->              do j = 4, ny - 3
   298: |||*--->                 do m = 1, 5
   299: ||||    F                   frct(m,i,j,k) = frct(m,i,j,k)
   300: ||||           >              - dssp * (            rsd(m,i,j-2,k)
   301: ||||           >                        - 4.0d+00 * rsd(m,i,j-1,k)
   302: ||||           >                        + 6.0d+00 * rsd(m,i,j,k)
   303: ||||           >                        - 4.0d+00 * rsd(m,i,j+1,k)
   304: ||||           >                        +           rsd(m,i,j+2,k) )
   305: |||*---                  end do
   306: ||V----               end do
   307: ||        
   308: ||*---->              do m = 1, 5
   309: |||                      frct(m,i,ny-2,k) = frct(m,i,ny-2,k)
   310: |||            >           - dssp * (             rsd(m,i,ny-4,k)
   311: |||            >                       - 4.0d+00 * rsd(m,i,ny-3,k)
   312: |||            >                       + 6.0d+00 * rsd(m,i,ny-2,k)
   313: |||            >                       - 4.0d+00 * rsd(m,i,ny-1,k)  )
   314: |||                      frct(m,i,ny-1,k) = frct(m,i,ny-1,k)
   315: |||            >           - dssp * (             rsd(m,i,ny-3,k)
   316: |||            >                       - 4.0d+00 * rsd(m,i,ny-2,k)
   317: |||            >                       + 5.0d+00 * rsd(m,i,ny-1,k)  )
   318: ||*----               end do
   319: ||        
   320: |+-----            end do
   321: Y------         end do
   322:           !$omp end do
   323:           
   324:           c---------------------------------------------------------------------
   325:           c   zeta-direction flux differences
   326:           c---------------------------------------------------------------------
   327:           !$omp do schedule(static) collapse(2)
   328: Y------>        do j = jst, jend
   329: |+----->           do i = ist, iend
   330: ||V---->              do k = 1, nz
   331: |||                      flux(1,k) = rsd(4,i,j,k)
   332: |||                      u41 = rsd(4,i,j,k) / rsd(1,i,j,k)
   333: |||     F                q = 0.50d+00 * (  rsd(2,i,j,k) * rsd(2,i,j,k)
   334: |||            >                         + rsd(3,i,j,k) * rsd(3,i,j,k)
   335: |||            >                         + rsd(4,i,j,k) * rsd(4,i,j,k) )
   336: |||            >                      / rsd(1,i,j,k)
   337: |||                      flux(2,k) = rsd(2,i,j,k) * u41
   338: |||                      flux(3,k) = rsd(3,i,j,k) * u41
   339: |||     F                flux(4,k) = rsd(4,i,j,k) * u41 + c2 *
   340: |||            >                         ( rsd(5,i,j,k) - q )
   341: |||     F                flux(5,k) = ( c1 * rsd(5,i,j,k) - c2 * q ) * u41
   342: ||V----               end do
   343: ||        
   344: ||V---->              do k = 2, nz - 1
   345: |||*--->                 do m = 1, 5
   346: ||||    F                   frct(m,i,j,k) =  frct(m,i,j,k)
   347: ||||           >                  - tz2 * ( flux(m,k+1) - flux(m,k-1) )
   348: |||*---                  end do
   349: ||V----               end do
   350: ||        
   351: ||V---->              do k = 2, nz
   352: |||                      tmp = 1.0d+00 / rsd(1,i,j,k)
   353: |||       
   354: |||                      u21k = tmp * rsd(2,i,j,k)
   355: |||                      u31k = tmp * rsd(3,i,j,k)
   356: |||                      u41k = tmp * rsd(4,i,j,k)
   357: |||                      u51k = tmp * rsd(5,i,j,k)
   358: |||       
   359: |||                      tmp = 1.0d+00 / rsd(1,i,j,k-1)
   360: |||       
   361: |||                      u21km1 = tmp * rsd(2,i,j,k-1)
   362: |||                      u31km1 = tmp * rsd(3,i,j,k-1)
   363: |||                      u41km1 = tmp * rsd(4,i,j,k-1)
   364: |||                      u51km1 = tmp * rsd(5,i,j,k-1)
   365: |||       
   366: |||                      flux(2,k) = tz3 * ( u21k - u21km1 )
   367: |||                      flux(3,k) = tz3 * ( u31k - u31km1 )
   368: |||                      flux(4,k) = (4.0d+00/3.0d+00) * tz3 * ( u41k
   369: |||            >                       - u41km1 )
   370: |||     F                flux(5,k) = 0.50d+00 * ( 1.0d+00 - c1*c5 )
   371: |||            >              * tz3 * ( ( u21k  **2 + u31k  **2 + u41k  **2 )
   372: |||            >                      - ( u21km1**2 + u31km1**2 + u41km1**2 ) )
   373: |||            >              + (1.0d+00/6.0d+00)
   374: |||            >              * tz3 * ( u41k**2 - u41km1**2 )
   375: |||            >              + c1 * c5 * tz3 * ( u51k - u51km1 )
   376: ||V----               end do
   377: ||        
   378: ||V---->              do k = 2, nz - 1
   379: |||     F                frct(1,i,j,k) = frct(1,i,j,k)
   380: |||            >              + dz1 * tz1 * (            rsd(1,i,j,k+1)
   381: |||            >                             - 2.0d+00 * rsd(1,i,j,k)
   382: |||            >                             +           rsd(1,i,j,k-1) )
   383: |||     F                frct(2,i,j,k) = frct(2,i,j,k)
   384: |||            >          + tz3 * c3 * c4 * ( flux(2,k+1) - flux(2,k) )
   385: |||            >              + dz2 * tz1 * (            rsd(2,i,j,k+1)
   386: |||            >                             - 2.0d+00 * rsd(2,i,j,k)
   387: |||            >                             +           rsd(2,i,j,k-1) )
   388: |||     F                frct(3,i,j,k) = frct(3,i,j,k)
   389: |||            >          + tz3 * c3 * c4 * ( flux(3,k+1) - flux(3,k) )
   390: |||            >              + dz3 * tz1 * (            rsd(3,i,j,k+1)
   391: |||            >                             - 2.0d+00 * rsd(3,i,j,k)
   392: |||            >                             +           rsd(3,i,j,k-1) )
   393: |||     F                frct(4,i,j,k) = frct(4,i,j,k)
   394: |||            >          + tz3 * c3 * c4 * ( flux(4,k+1) - flux(4,k) )
   395: |||            >              + dz4 * tz1 * (            rsd(4,i,j,k+1)
   396: |||            >                             - 2.0d+00 * rsd(4,i,j,k)
   397: |||            >                             +           rsd(4,i,j,k-1) )
   398: |||     F                frct(5,i,j,k) = frct(5,i,j,k)
   399: |||            >          + tz3 * c3 * c4 * ( flux(5,k+1) - flux(5,k) )
   400: |||            >              + dz5 * tz1 * (            rsd(5,i,j,k+1)
   401: |||            >                             - 2.0d+00 * rsd(5,i,j,k)
   402: |||            >                             +           rsd(5,i,j,k-1) )
   403: ||V----               end do
   404: ||        
   405: ||        c---------------------------------------------------------------------
   406: ||        c   fourth-order dissipation
   407: ||        c---------------------------------------------------------------------
   408: ||*---->              do m = 1, 5
   409: |||                      frct(m,i,j,2) = frct(m,i,j,2)
   410: |||            >           - dssp * ( + 5.0d+00 * rsd(m,i,j,2)
   411: |||            >                       - 4.0d+00 * rsd(m,i,j,3)
   412: |||            >                       +           rsd(m,i,j,4) )
   413: |||                      frct(m,i,j,3) = frct(m,i,j,3)
   414: |||            >           - dssp * (- 4.0d+00 * rsd(m,i,j,2)
   415: |||            >                      + 6.0d+00 * rsd(m,i,j,3)
   416: |||            >                      - 4.0d+00 * rsd(m,i,j,4)
   417: |||            >                      +           rsd(m,i,j,5) )
   418: ||*----               end do
   419: ||        
   420: ||V---->              do k = 4, nz - 3
   421: |||*--->                 do m = 1, 5
   422: ||||    F                   frct(m,i,j,k) = frct(m,i,j,k)
   423: ||||           >              - dssp * (           rsd(m,i,j,k-2)
   424: ||||           >                        - 4.0d+00 * rsd(m,i,j,k-1)
   425: ||||           >                        + 6.0d+00 * rsd(m,i,j,k)
   426: ||||           >                        - 4.0d+00 * rsd(m,i,j,k+1)
   427: ||||           >                        +           rsd(m,i,j,k+2) )
   428: |||*---                  end do
   429: ||V----               end do
   430: ||        
   431: ||*---->              do m = 1, 5
   432: |||                      frct(m,i,j,nz-2) = frct(m,i,j,nz-2)
   433: |||            >           - dssp * (            rsd(m,i,j,nz-4)
   434: |||            >                      - 4.0d+00 * rsd(m,i,j,nz-3)
   435: |||            >                      + 6.0d+00 * rsd(m,i,j,nz-2)
   436: |||            >                      - 4.0d+00 * rsd(m,i,j,nz-1)  )
   437: |||                      frct(m,i,j,nz-1) = frct(m,i,j,nz-1)
   438: |||            >           - dssp * (             rsd(m,i,j,nz-3)
   439: |||            >                       - 4.0d+00 * rsd(m,i,j,nz-2)
   440: |||            >                       + 5.0d+00 * rsd(m,i,j,nz-1)  )
   441: ||*----               end do
   442: |+-----            end do
   443: Y------         end do
   444:           !$omp end do nowait
   445:           !$omp end parallel
   446:           
   447:                 return
   448:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: erhs.f

PROCEDURE NAME: ERHS
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (erhs.f:39)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (erhs.f:39)
    *** The number of VLOAD, VSTORE. :  0,  2. (erhs.f:39)

    LOOP BEGIN: (erhs.f:41)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:41)
      *** The number of VLOAD, VSTORE. :  0,  5. (erhs.f:41)
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:51)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (erhs.f:51)
    *** The number of VLOAD, VSTORE. :  0,  4. (erhs.f:51)

    LOOP BEGIN: (erhs.f:55)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:55)
      *** The number of VLOAD, VSTORE. :  0,  5. (erhs.f:55)
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:82)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (erhs.f:82)
    *** The number of VLOAD, VSTORE. :  0,  2. (erhs.f:82)

    LOOP BEGIN: (erhs.f:84)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:84)
      *** The number of VLOAD, VSTORE. :  5,  5. (erhs.f:84)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:84)
    LOOP END

    LOOP BEGIN: (erhs.f:98)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:98)
      *** The number of VLOAD, VSTORE. : 15,  5. (erhs.f:98)
    LOOP END

    LOOP BEGIN: (erhs.f:104)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:104)
      *** The number of VLOAD, VSTORE. : 10,  4. (erhs.f:104)
    LOOP END

    LOOP BEGIN: (erhs.f:131)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:131)
      *** The number of VLOAD, VSTORE. : 28,  5. (erhs.f:131)
    LOOP END

    LOOP BEGIN: (erhs.f:173)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:173)
      *** The number of VLOAD, VSTORE. : 30,  5. (erhs.f:173)
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:205)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (erhs.f:205)
    *** The number of VLOAD, VSTORE. :  0,  2. (erhs.f:205)

    LOOP BEGIN: (erhs.f:207)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:207)
      *** The number of VLOAD, VSTORE. :  5,  5. (erhs.f:207)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:207)
    LOOP END

    LOOP BEGIN: (erhs.f:221)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:221)
      *** The number of VLOAD, VSTORE. : 15,  5. (erhs.f:221)
    LOOP END

    LOOP BEGIN: (erhs.f:228)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:228)
      *** The number of VLOAD, VSTORE. : 10,  4. (erhs.f:228)
    LOOP END

    LOOP BEGIN: (erhs.f:255)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:255)
      *** The number of VLOAD, VSTORE. : 28,  5. (erhs.f:255)
    LOOP END

    LOOP BEGIN: (erhs.f:297)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:297)
      *** The number of VLOAD, VSTORE. : 30,  5. (erhs.f:297)
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:328)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (erhs.f:328)
    *** The number of VLOAD, VSTORE. :  0,  2. (erhs.f:328)

    LOOP BEGIN: (erhs.f:330)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:330)
      *** The number of VLOAD, VSTORE. :  5,  5. (erhs.f:330)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:330)
    LOOP END

    LOOP BEGIN: (erhs.f:344)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:344)
      *** The number of VLOAD, VSTORE. : 15,  5. (erhs.f:344)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:344)
    LOOP END

    LOOP BEGIN: (erhs.f:351)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:351)
      *** The number of VLOAD, VSTORE. : 10,  4. (erhs.f:351)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:351)
    LOOP END

    LOOP BEGIN: (erhs.f:378)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:378)
      *** The number of VLOAD, VSTORE. : 28,  5. (erhs.f:378)
      *** Iteration count is assumed. Iteration count. : 162 (erhs.f:378)
    LOOP END

    LOOP BEGIN: (erhs.f:420)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (erhs.f:420)
      *** The number of VLOAD, VSTORE. : 30,  5. (erhs.f:420)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: erhs.f

PROCEDURE NAME: ERHS
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
    Total size                 : 5104 bytes
      Register spill area      : 2448 bytes
      Parameter area           :   16 bytes
      Register save area       :  176 bytes
      User data area           : 2464 bytes
      Others                   :    0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (erhs.f:39)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (erhs.f:39)
    *** Estimated execution cycle                       : 24
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (erhs.f:41)
      *** Estimated execution cycle                     : 167
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:51)
    *** Estimated execution cycle                       : 737
  LOOP END

  LOOP BEGIN: (erhs.f:51)
    *** Estimated execution cycle                       : 252
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 43
            Not enough registers                        :  2
            Over basic blocks                           : 41
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 35
            Not enough registers                        :  2
            Over basic blocks                           : 33
    *** The number of SCALAR REGISTER TRANSFER          : 1

    LOOP BEGIN: (erhs.f:55)
      *** Estimated execution cycle                     : 903
      *** The number of SCALAR REGISTER TRANSFER        : 7
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:82)
    *** Estimated execution cycle                       : 575
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
  LOOP END

  LOOP BEGIN: (erhs.f:82)
    *** Estimated execution cycle                       : 990
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 57
            Over basic blocks                           : 57
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 52
            Over basic blocks                           : 52
    *** The number of SCALAR REGISTER TRANSFER          : 95

    LOOP BEGIN: (erhs.f:98)
      *** Estimated execution cycle                     : 813
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 18
    LOOP END

    LOOP BEGIN: (erhs.f:104)
      *** Estimated execution cycle                     : 1904
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 7
              Over basic blocks                         : 7
      *** The number of SCALAR REGISTER TRANSFER        : 20
    LOOP END

    LOOP BEGIN: (erhs.f:131)
      *** Estimated execution cycle                     : 1510
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 18
              Over basic blocks                         : 18
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (erhs.f:173)
      *** Estimated execution cycle                     : 1287
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 4
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 19
              Over basic blocks                         : 19
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:205)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (erhs.f:205)
    *** Estimated execution cycle                       : 993
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 53
            Over basic blocks                           : 53
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 54
            Over basic blocks                           : 54
    *** The number of SCALAR REGISTER TRANSFER          : 86

    LOOP BEGIN: (erhs.f:221)
      *** Estimated execution cycle                     : 813
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 19
    LOOP END

    LOOP BEGIN: (erhs.f:228)
      *** Estimated execution cycle                     : 1873
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 6
              Over basic blocks                         : 6
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 5
              Over basic blocks                         : 5
      *** The number of SCALAR REGISTER TRANSFER        : 20
    LOOP END

    LOOP BEGIN: (erhs.f:255)
      *** Estimated execution cycle                     : 1510
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 16
              Over basic blocks                         : 16
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (erhs.f:297)
      *** Estimated execution cycle                     : 1287
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 4
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 17
              Over basic blocks                         : 17
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (erhs.f:328)
    *** Estimated execution cycle                       : 575
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
  LOOP END

  LOOP BEGIN: (erhs.f:328)
    *** Estimated execution cycle                       : 1921
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 43
            Over basic blocks                           : 43
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 59
            Over basic blocks                           : 59
    *** The number of SCALAR REGISTER TRANSFER          : 37

    LOOP BEGIN: (erhs.f:420)
      *** Estimated execution cycle                     : 1287
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 16
              Over basic blocks                         : 16
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END


