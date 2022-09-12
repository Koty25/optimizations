NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME : jacu.f

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
FILE NAME: jacu.f

PROCEDURE NAME: JACU
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: vec( 101): Vectorized loop.
    58: vec( 128): Fused multiply-add operation applied.
    72: vec( 128): Fused multiply-add operation applied.
    86: vec( 128): Fused multiply-add operation applied.
    94: vec( 128): Fused multiply-add operation applied.
   119: vec( 128): Fused multiply-add operation applied.
   138: vec( 128): Fused multiply-add operation applied.
   142: vec( 128): Fused multiply-add operation applied.
   152: vec( 128): Fused multiply-add operation applied.
   156: vec( 128): Fused multiply-add operation applied.
   162: vec( 128): Fused multiply-add operation applied.
   172: vec( 128): Fused multiply-add operation applied.
   181: vec( 128): Fused multiply-add operation applied.
   188: vec( 128): Fused multiply-add operation applied.
   192: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   214: vec( 128): Fused multiply-add operation applied.
   217: vec( 128): Fused multiply-add operation applied.
   224: vec( 128): Fused multiply-add operation applied.
   230: vec( 128): Fused multiply-add operation applied.
   238: vec( 128): Fused multiply-add operation applied.
   248: vec( 128): Fused multiply-add operation applied.
   257: vec( 128): Fused multiply-add operation applied.
   261: vec( 128): Fused multiply-add operation applied.
   268: vec( 128): Fused multiply-add operation applied.
   271: vec( 128): Fused multiply-add operation applied.
   289: vec( 128): Fused multiply-add operation applied.
   292: vec( 128): Fused multiply-add operation applied.
   299: vec( 128): Fused multiply-add operation applied.
   309: vec( 128): Fused multiply-add operation applied.
   317: vec( 128): Fused multiply-add operation applied.
   323: vec( 128): Fused multiply-add operation applied.
   332: vec( 128): Fused multiply-add operation applied.
   335: vec( 128): Fused multiply-add operation applied.
   338: vec( 128): Fused multiply-add operation applied.
   344: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacu.f

PROCEDURE NAME: JACU
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

     5:                 subroutine jacu(j, k)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   compute the upper triangular part of the jacobian matrix
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use lu_data
    15:                 implicit none
    16:           
    17:           c---------------------------------------------------------------------
    18:           c  input parameters
    19:           c---------------------------------------------------------------------
    20:                 integer j, k
    21:           
    22:           c---------------------------------------------------------------------
    23:           c  local variables
    24:           c---------------------------------------------------------------------
    25:                 integer i
    26:                 double precision  r43
    27:                 double precision  c1345
    28:                 double precision  c34
    29:                 double precision  tmp1, tmp2, tmp3
    30:           
    31:           
    32:           
    33:                 r43 = ( 4.0d+00 / 3.0d+00 )
    34:                 c1345 = c1 * c3 * c4 * c5
    35:                 c34 = c3 * c4
    36:           
    37: V------>              do i = iend, ist, -1
    38: |         
    39: |         c---------------------------------------------------------------------
    40: |         c   form the block daigonal
    41: |         c---------------------------------------------------------------------
    42: |                        tmp1 = rho_i(i,j,k)
    43: |                        tmp2 = tmp1 * tmp1
    44: |                        tmp3 = tmp1 * tmp2
    45: |         
    46: |                        d(1,1,i) =  1.0d+00
    47: |              >                       + dt * 2.0d+00 * (   tx1 * dx1
    48: |              >                                          + ty1 * dy1
    49: |              >                                          + tz1 * dz1 )
    50: |                        d(1,2,i) =  0.0d+00
    51: |                        d(1,3,i) =  0.0d+00
    52: |                        d(1,4,i) =  0.0d+00
    53: |                        d(1,5,i) =  0.0d+00
    54: |         
    55: |                        d(2,1,i) =  dt * 2.0d+00
    56: |              >           * ( - tx1 * r43 - ty1 - tz1 )
    57: |              >           * ( c34 * tmp2 * u(2,i,j,k) )
    58: |       F                d(2,2,i) =  1.0d+00
    59: |              >          + dt * 2.0d+00 * c34 * tmp1
    60: |              >          * (  tx1 * r43 + ty1 + tz1 )
    61: |              >          + dt * 2.0d+00 * (   tx1 * dx2
    62: |              >                             + ty1 * dy2
    63: |              >                             + tz1 * dz2  )
    64: |                        d(2,3,i) = 0.0d+00
    65: |                        d(2,4,i) = 0.0d+00
    66: |                        d(2,5,i) = 0.0d+00
    67: |         
    68: |                        d(3,1,i) = dt * 2.0d+00
    69: |              >           * ( - tx1 - ty1 * r43 - tz1 )
    70: |              >           * ( c34 * tmp2 * u(3,i,j,k) )
    71: |                        d(3,2,i) = 0.0d+00
    72: |       F                d(3,3,i) = 1.0d+00
    73: |              >         + dt * 2.0d+00 * c34 * tmp1
    74: |              >              * (  tx1 + ty1 * r43 + tz1 )
    75: |              >         + dt * 2.0d+00 * (  tx1 * dx3
    76: |              >                           + ty1 * dy3
    77: |              >                           + tz1 * dz3 )
    78: |                        d(3,4,i) = 0.0d+00
    79: |                        d(3,5,i) = 0.0d+00
    80: |         
    81: |                        d(4,1,i) = dt * 2.0d+00
    82: |              >           * ( - tx1 - ty1 - tz1 * r43 )
    83: |              >           * ( c34 * tmp2 * u(4,i,j,k) )
    84: |                        d(4,2,i) = 0.0d+00
    85: |                        d(4,3,i) = 0.0d+00
    86: |       F                d(4,4,i) = 1.0d+00
    87: |              >         + dt * 2.0d+00 * c34 * tmp1
    88: |              >              * (  tx1 + ty1 + tz1 * r43 )
    89: |              >         + dt * 2.0d+00 * (  tx1 * dx4
    90: |              >                           + ty1 * dy4
    91: |              >                           + tz1 * dz4 )
    92: |                        d(4,5,i) = 0.0d+00
    93: |         
    94: |       F                d(5,1,i) = -dt * 2.0d+00
    95: |              >  * ( ( ( tx1 * ( r43*c34 - c1345 )
    96: |              >     + ty1 * ( c34 - c1345 )
    97: |              >     + tz1 * ( c34 - c1345 ) ) * ( u(2,i,j,k) ** 2 )
    98: |              >   + ( tx1 * ( c34 - c1345 )
    99: |              >     + ty1 * ( r43*c34 - c1345 )
   100: |              >     + tz1 * ( c34 - c1345 ) ) * ( u(3,i,j,k) ** 2 )
   101: |              >   + ( tx1 * ( c34 - c1345 )
   102: |              >     + ty1 * ( c34 - c1345 )
   103: |              >     + tz1 * ( r43*c34 - c1345 ) ) * ( u(4,i,j,k) ** 2 )
   104: |              >      ) * tmp3
   105: |              >   + ( tx1 + ty1 + tz1 ) * c1345 * tmp2 * u(5,i,j,k) )
   106: |         
   107: |                        d(5,2,i) = dt * 2.0d+00
   108: |              > * ( tx1 * ( r43*c34 - c1345 )
   109: |              >   + ty1 * (     c34 - c1345 )
   110: |              >   + tz1 * (     c34 - c1345 ) ) * tmp2 * u(2,i,j,k)
   111: |                        d(5,3,i) = dt * 2.0d+00
   112: |              > * ( tx1 * ( c34 - c1345 )
   113: |              >   + ty1 * ( r43*c34 -c1345 )
   114: |              >   + tz1 * ( c34 - c1345 ) ) * tmp2 * u(3,i,j,k)
   115: |                        d(5,4,i) = dt * 2.0d+00
   116: |              > * ( tx1 * ( c34 - c1345 )
   117: |              >   + ty1 * ( c34 - c1345 )
   118: |              >   + tz1 * ( r43*c34 - c1345 ) ) * tmp2 * u(4,i,j,k)
   119: |       F                d(5,5,i) = 1.0d+00
   120: |              >   + dt * 2.0d+00 * ( tx1 + ty1 + tz1 ) * c1345 * tmp1
   121: |              >   + dt * 2.0d+00 * (  tx1 * dx5
   122: |              >                    +  ty1 * dy5
   123: |              >                    +  tz1 * dz5 )
   124: |         
   125: |         c---------------------------------------------------------------------
   126: |         c   form the first block sub-diagonal
   127: |         c---------------------------------------------------------------------
   128: |                        tmp1 = rho_i(i+1,j,k)
   129: |                        tmp2 = tmp1 * tmp1
   130: |                        tmp3 = tmp1 * tmp2
   131: |         
   132: |                        a(1,1,i) = - dt * tx1 * dx1
   133: |                        a(1,2,i) =   dt * tx2
   134: |                        a(1,3,i) =   0.0d+00
   135: |                        a(1,4,i) =   0.0d+00
   136: |                        a(1,5,i) =   0.0d+00
   137: |         
   138: |       F                a(2,1,i) =  dt * tx2
   139: |              >          * ( - ( u(2,i+1,j,k) * tmp1 ) ** 2
   140: |              >     + c2 * qs(i+1,j,k) * tmp1 )
   141: |              >          - dt * tx1 * ( - r43 * c34 * tmp2 * u(2,i+1,j,k) )
   142: |       F                a(2,2,i) =  dt * tx2
   143: |              >          * ( ( 2.0d+00 - c2 ) * ( u(2,i+1,j,k) * tmp1 ) )
   144: |              >          - dt * tx1 * ( r43 * c34 * tmp1 )
   145: |              >          - dt * tx1 * dx2
   146: |                        a(2,3,i) =  dt * tx2
   147: |              >              * ( - c2 * ( u(3,i+1,j,k) * tmp1 ) )
   148: |                        a(2,4,i) =  dt * tx2
   149: |              >              * ( - c2 * ( u(4,i+1,j,k) * tmp1 ) )
   150: |                        a(2,5,i) =  dt * tx2 * c2
   151: |         
   152: |       F                a(3,1,i) =  dt * tx2
   153: |              >              * ( - ( u(2,i+1,j,k) * u(3,i+1,j,k) ) * tmp2 )
   154: |              >         - dt * tx1 * ( - c34 * tmp2 * u(3,i+1,j,k) )
   155: |                        a(3,2,i) =  dt * tx2 * ( u(3,i+1,j,k) * tmp1 )
   156: |       F                a(3,3,i) =  dt * tx2 * ( u(2,i+1,j,k) * tmp1 )
   157: |              >          - dt * tx1 * ( c34 * tmp1 )
   158: |              >          - dt * tx1 * dx3
   159: |                        a(3,4,i) = 0.0d+00
   160: |                        a(3,5,i) = 0.0d+00
   161: |         
   162: |       F                a(4,1,i) = dt * tx2
   163: |              >          * ( - ( u(2,i+1,j,k)*u(4,i+1,j,k) ) * tmp2 )
   164: |              >          - dt * tx1 * ( - c34 * tmp2 * u(4,i+1,j,k) )
   165: |                        a(4,2,i) = dt * tx2 * ( u(4,i+1,j,k) * tmp1 )
   166: |                        a(4,3,i) = 0.0d+00
   167: |                        a(4,4,i) = dt * tx2 * ( u(2,i+1,j,k) * tmp1 )
   168: |              >          - dt * tx1 * ( c34 * tmp1 )
   169: |              >          - dt * tx1 * dx4
   170: |                        a(4,5,i) = 0.0d+00
   171: |         
   172: |       F                a(5,1,i) = dt * tx2
   173: |              >          * ( ( c2 * 2.0d0 * qs(i+1,j,k)
   174: |              >              - c1 * u(5,i+1,j,k) )
   175: |              >          * ( u(2,i+1,j,k) * tmp2 ) )
   176: |              >          - dt * tx1
   177: |              >          * ( - ( r43*c34 - c1345 ) * tmp3 * ( u(2,i+1,j,k)**2 )
   178: |              >              - (     c34 - c1345 ) * tmp3 * ( u(3,i+1,j,k)**2 )
   179: |              >              - (     c34 - c1345 ) * tmp3 * ( u(4,i+1,j,k)**2 )
   180: |              >              - c1345 * tmp2 * u(5,i+1,j,k) )
   181: |       F                a(5,2,i) = dt * tx2
   182: |              >          * ( c1 * ( u(5,i+1,j,k) * tmp1 )
   183: |              >             - c2
   184: |              >             * (  u(2,i+1,j,k)*u(2,i+1,j,k) * tmp2
   185: |              >                  + qs(i+1,j,k) * tmp1 ) )
   186: |              >           - dt * tx1
   187: |              >           * ( r43*c34 - c1345 ) * tmp2 * u(2,i+1,j,k)
   188: |       F                a(5,3,i) = dt * tx2
   189: |              >           * ( - c2 * ( u(3,i+1,j,k)*u(2,i+1,j,k) ) * tmp2 )
   190: |              >           - dt * tx1
   191: |              >           * (  c34 - c1345 ) * tmp2 * u(3,i+1,j,k)
   192: |       F                a(5,4,i) = dt * tx2
   193: |              >           * ( - c2 * ( u(4,i+1,j,k)*u(2,i+1,j,k) ) * tmp2 )
   194: |              >           - dt * tx1
   195: |              >           * (  c34 - c1345 ) * tmp2 * u(4,i+1,j,k)
   196: |       F                a(5,5,i) = dt * tx2
   197: |              >           * ( c1 * ( u(2,i+1,j,k) * tmp1 ) )
   198: |              >           - dt * tx1 * c1345 * tmp1
   199: |              >           - dt * tx1 * dx5
   200: |         
   201: |         c---------------------------------------------------------------------
   202: |         c   form the second block sub-diagonal
   203: |         c---------------------------------------------------------------------
   204: |                        tmp1 = rho_i(i,j+1,k)
   205: |                        tmp2 = tmp1 * tmp1
   206: |                        tmp3 = tmp1 * tmp2
   207: |         
   208: |                        b(1,1,i) = - dt * ty1 * dy1
   209: |                        b(1,2,i) =   0.0d+00
   210: |                        b(1,3,i) =  dt * ty2
   211: |                        b(1,4,i) =   0.0d+00
   212: |                        b(1,5,i) =   0.0d+00
   213: |         
   214: |       F                b(2,1,i) =  dt * ty2
   215: |              >           * ( - ( u(2,i,j+1,k)*u(3,i,j+1,k) ) * tmp2 )
   216: |              >           - dt * ty1 * ( - c34 * tmp2 * u(2,i,j+1,k) )
   217: |       F                b(2,2,i) =  dt * ty2 * ( u(3,i,j+1,k) * tmp1 )
   218: |              >          - dt * ty1 * ( c34 * tmp1 )
   219: |              >          - dt * ty1 * dy2
   220: |                        b(2,3,i) =  dt * ty2 * ( u(2,i,j+1,k) * tmp1 )
   221: |                        b(2,4,i) = 0.0d+00
   222: |                        b(2,5,i) = 0.0d+00
   223: |         
   224: |       F                b(3,1,i) =  dt * ty2
   225: |              >           * ( - ( u(3,i,j+1,k) * tmp1 ) ** 2
   226: |              >      + c2 * ( qs(i,j+1,k) * tmp1 ) )
   227: |              >       - dt * ty1 * ( - r43 * c34 * tmp2 * u(3,i,j+1,k) )
   228: |                        b(3,2,i) =  dt * ty2
   229: |              >                   * ( - c2 * ( u(2,i,j+1,k) * tmp1 ) )
   230: |       F                b(3,3,i) =  dt * ty2 * ( ( 2.0d+00 - c2 )
   231: |              >                   * ( u(3,i,j+1,k) * tmp1 ) )
   232: |              >       - dt * ty1 * ( r43 * c34 * tmp1 )
   233: |              >       - dt * ty1 * dy3
   234: |                        b(3,4,i) =  dt * ty2
   235: |              >                   * ( - c2 * ( u(4,i,j+1,k) * tmp1 ) )
   236: |                        b(3,5,i) =  dt * ty2 * c2
   237: |         
   238: |       F                b(4,1,i) =  dt * ty2
   239: |              >              * ( - ( u(3,i,j+1,k)*u(4,i,j+1,k) ) * tmp2 )
   240: |              >       - dt * ty1 * ( - c34 * tmp2 * u(4,i,j+1,k) )
   241: |                        b(4,2,i) = 0.0d+00
   242: |                        b(4,3,i) =  dt * ty2 * ( u(4,i,j+1,k) * tmp1 )
   243: |                        b(4,4,i) =  dt * ty2 * ( u(3,i,j+1,k) * tmp1 )
   244: |              >                        - dt * ty1 * ( c34 * tmp1 )
   245: |              >                        - dt * ty1 * dy4
   246: |                        b(4,5,i) = 0.0d+00
   247: |         
   248: |       F                b(5,1,i) =  dt * ty2
   249: |              >          * ( ( c2 * 2.0d0 * qs(i,j+1,k)
   250: |              >               - c1 * u(5,i,j+1,k) )
   251: |              >          * ( u(3,i,j+1,k) * tmp2 ) )
   252: |              >          - dt * ty1
   253: |              >          * ( - (     c34 - c1345 )*tmp3*(u(2,i,j+1,k)**2)
   254: |              >              - ( r43*c34 - c1345 )*tmp3*(u(3,i,j+1,k)**2)
   255: |              >              - (     c34 - c1345 )*tmp3*(u(4,i,j+1,k)**2)
   256: |              >              - c1345*tmp2*u(5,i,j+1,k) )
   257: |       F                b(5,2,i) =  dt * ty2
   258: |              >          * ( - c2 * ( u(2,i,j+1,k)*u(3,i,j+1,k) ) * tmp2 )
   259: |              >          - dt * ty1
   260: |              >          * ( c34 - c1345 ) * tmp2 * u(2,i,j+1,k)
   261: |       F                b(5,3,i) =  dt * ty2
   262: |              >          * ( c1 * ( u(5,i,j+1,k) * tmp1 )
   263: |              >          - c2
   264: |              >          * ( qs(i,j+1,k) * tmp1
   265: |              >               + u(3,i,j+1,k)*u(3,i,j+1,k) * tmp2 ) )
   266: |              >          - dt * ty1
   267: |              >          * ( r43*c34 - c1345 ) * tmp2 * u(3,i,j+1,k)
   268: |       F                b(5,4,i) =  dt * ty2
   269: |              >          * ( - c2 * ( u(3,i,j+1,k)*u(4,i,j+1,k) ) * tmp2 )
   270: |              >          - dt * ty1 * ( c34 - c1345 ) * tmp2 * u(4,i,j+1,k)
   271: |       F                b(5,5,i) =  dt * ty2
   272: |              >          * ( c1 * ( u(3,i,j+1,k) * tmp1 ) )
   273: |              >          - dt * ty1 * c1345 * tmp1
   274: |              >          - dt * ty1 * dy5
   275: |         
   276: |         c---------------------------------------------------------------------
   277: |         c   form the third block sub-diagonal
   278: |         c---------------------------------------------------------------------
   279: |                        tmp1 = rho_i(i,j,k+1)
   280: |                        tmp2 = tmp1 * tmp1
   281: |                        tmp3 = tmp1 * tmp2
   282: |         
   283: |                        c(1,1,i) = - dt * tz1 * dz1
   284: |                        c(1,2,i) =   0.0d+00
   285: |                        c(1,3,i) =   0.0d+00
   286: |                        c(1,4,i) = dt * tz2
   287: |                        c(1,5,i) =   0.0d+00
   288: |         
   289: |       F                c(2,1,i) = dt * tz2
   290: |              >           * ( - ( u(2,i,j,k+1)*u(4,i,j,k+1) ) * tmp2 )
   291: |              >           - dt * tz1 * ( - c34 * tmp2 * u(2,i,j,k+1) )
   292: |       F                c(2,2,i) = dt * tz2 * ( u(4,i,j,k+1) * tmp1 )
   293: |              >           - dt * tz1 * c34 * tmp1
   294: |              >           - dt * tz1 * dz2
   295: |                        c(2,3,i) = 0.0d+00
   296: |                        c(2,4,i) = dt * tz2 * ( u(2,i,j,k+1) * tmp1 )
   297: |                        c(2,5,i) = 0.0d+00
   298: |         
   299: |       F                c(3,1,i) = dt * tz2
   300: |              >           * ( - ( u(3,i,j,k+1)*u(4,i,j,k+1) ) * tmp2 )
   301: |              >           - dt * tz1 * ( - c34 * tmp2 * u(3,i,j,k+1) )
   302: |                        c(3,2,i) = 0.0d+00
   303: |                        c(3,3,i) = dt * tz2 * ( u(4,i,j,k+1) * tmp1 )
   304: |              >           - dt * tz1 * ( c34 * tmp1 )
   305: |              >           - dt * tz1 * dz3
   306: |                        c(3,4,i) = dt * tz2 * ( u(3,i,j,k+1) * tmp1 )
   307: |                        c(3,5,i) = 0.0d+00
   308: |         
   309: |       F                c(4,1,i) = dt * tz2
   310: |              >        * ( - ( u(4,i,j,k+1) * tmp1 ) ** 2
   311: |              >            + c2 * ( qs(i,j,k+1) * tmp1 ) )
   312: |              >        - dt * tz1 * ( - r43 * c34 * tmp2 * u(4,i,j,k+1) )
   313: |                        c(4,2,i) = dt * tz2
   314: |              >             * ( - c2 * ( u(2,i,j,k+1) * tmp1 ) )
   315: |                        c(4,3,i) = dt * tz2
   316: |              >             * ( - c2 * ( u(3,i,j,k+1) * tmp1 ) )
   317: |       F                c(4,4,i) = dt * tz2 * ( 2.0d+00 - c2 )
   318: |              >             * ( u(4,i,j,k+1) * tmp1 )
   319: |              >             - dt * tz1 * ( r43 * c34 * tmp1 )
   320: |              >             - dt * tz1 * dz4
   321: |                        c(4,5,i) = dt * tz2 * c2
   322: |         
   323: |       F                c(5,1,i) = dt * tz2
   324: |              >     * ( ( c2 * 2.0d0 * qs(i,j,k+1)
   325: |              >       - c1 * u(5,i,j,k+1) )
   326: |              >            * ( u(4,i,j,k+1) * tmp2 ) )
   327: |              >       - dt * tz1
   328: |              >       * ( - ( c34 - c1345 ) * tmp3 * (u(2,i,j,k+1)**2)
   329: |              >           - ( c34 - c1345 ) * tmp3 * (u(3,i,j,k+1)**2)
   330: |              >           - ( r43*c34 - c1345 )* tmp3 * (u(4,i,j,k+1)**2)
   331: |              >          - c1345 * tmp2 * u(5,i,j,k+1) )
   332: |       F                c(5,2,i) = dt * tz2
   333: |              >       * ( - c2 * ( u(2,i,j,k+1)*u(4,i,j,k+1) ) * tmp2 )
   334: |              >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * u(2,i,j,k+1)
   335: |       F                c(5,3,i) = dt * tz2
   336: |              >       * ( - c2 * ( u(3,i,j,k+1)*u(4,i,j,k+1) ) * tmp2 )
   337: |              >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * u(3,i,j,k+1)
   338: |       F                c(5,4,i) = dt * tz2
   339: |              >       * ( c1 * ( u(5,i,j,k+1) * tmp1 )
   340: |              >       - c2
   341: |              >       * ( qs(i,j,k+1) * tmp1
   342: |              >            + u(4,i,j,k+1)*u(4,i,j,k+1) * tmp2 ) )
   343: |              >       - dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u(4,i,j,k+1)
   344: |       F                c(5,5,i) = dt * tz2
   345: |              >       * ( c1 * ( u(4,i,j,k+1) * tmp1 ) )
   346: |              >       - dt * tz1 * c1345 * tmp1
   347: |              >       - dt * tz1 * dz5
   348: |         
   349: V------               end do
   350:           
   351:           
   352:                 return
   353:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacu.f

PROCEDURE NAME: JACU
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacu.f

PROCEDURE NAME: JACU
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (jacu.f:37)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (jacu.f:37)
    *** The number of VLOAD, VSTORE. : 23, 100. (jacu.f:37)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacu.f

PROCEDURE NAME: JACU
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 58 [s0-s11 s18-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1167216 bytes
      Register spill area      :  583512 bytes
      Parameter area           :       0 bytes
      Register save area       :     176 bytes
      User data area           :  583512 bytes
      Others                   :       8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (jacu.f:37)
    *** Estimated execution cycle                       : 2100
    *** The number of VECTOR REGISTER SPILL
          Total                                         : 5
            Not enough registers                        : 5
    *** The number of VECTOR REGISTER RESTORE
          Total                                         : 7
            Over basic blocks                           : 7
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 18
            Not enough registers                        : 17
            Over basic blocks                           :  1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 186
            Over basic blocks                           : 186
    *** The number of SCALAR REGISTER TRANSFER          : 6
  LOOP END


