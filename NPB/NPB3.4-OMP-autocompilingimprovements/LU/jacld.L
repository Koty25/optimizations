NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME : jacld.f

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
FILE NAME: jacld.f

PROCEDURE NAME: JACLD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    38: vec( 101): Vectorized loop.
    59: vec( 128): Fused multiply-add operation applied.
    73: vec( 128): Fused multiply-add operation applied.
    87: vec( 128): Fused multiply-add operation applied.
    95: vec( 128): Fused multiply-add operation applied.
   120: vec( 128): Fused multiply-add operation applied.
   139: vec( 128): Fused multiply-add operation applied.
   142: vec( 128): Fused multiply-add operation applied.
   149: vec( 128): Fused multiply-add operation applied.
   159: vec( 128): Fused multiply-add operation applied.
   167: vec( 128): Fused multiply-add operation applied.
   173: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   185: vec( 128): Fused multiply-add operation applied.
   188: vec( 128): Fused multiply-add operation applied.
   194: vec( 128): Fused multiply-add operation applied.
   212: vec( 128): Fused multiply-add operation applied.
   215: vec( 128): Fused multiply-add operation applied.
   222: vec( 128): Fused multiply-add operation applied.
   228: vec( 128): Fused multiply-add operation applied.
   236: vec( 128): Fused multiply-add operation applied.
   246: vec( 128): Fused multiply-add operation applied.
   255: vec( 128): Fused multiply-add operation applied.
   259: vec( 128): Fused multiply-add operation applied.
   266: vec( 128): Fused multiply-add operation applied.
   269: vec( 128): Fused multiply-add operation applied.
   287: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   301: vec( 128): Fused multiply-add operation applied.
   305: vec( 128): Fused multiply-add operation applied.
   311: vec( 128): Fused multiply-add operation applied.
   321: vec( 128): Fused multiply-add operation applied.
   330: vec( 128): Fused multiply-add operation applied.
   337: vec( 128): Fused multiply-add operation applied.
   341: vec( 128): Fused multiply-add operation applied.
   345: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacld.f

PROCEDURE NAME: JACLD
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

     5:                 subroutine jacld(j, k)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           
    11:           c---------------------------------------------------------------------
    12:           c   compute the lower triangular part of the jacobian matrix
    13:           c---------------------------------------------------------------------
    14:           
    15:                 use lu_data
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  input parameters
    20:           c---------------------------------------------------------------------
    21:                 integer j, k
    22:           
    23:           c---------------------------------------------------------------------
    24:           c  local variables
    25:           c---------------------------------------------------------------------
    26:                 integer i
    27:                 double precision  r43
    28:                 double precision  c1345
    29:                 double precision  c34
    30:                 double precision  tmp1, tmp2, tmp3
    31:           
    32:           
    33:           
    34:                 r43 = ( 4.0d+00 / 3.0d+00 )
    35:                 c1345 = c1 * c3 * c4 * c5
    36:                 c34 = c3 * c4
    37:           
    38: V------>              do i = ist, iend
    39: |         
    40: |         c---------------------------------------------------------------------
    41: |         c   form the block daigonal
    42: |         c---------------------------------------------------------------------
    43: |                        tmp1 = rho_i(i,j,k)
    44: |                        tmp2 = tmp1 * tmp1
    45: |                        tmp3 = tmp1 * tmp2
    46: |         
    47: |                        d(1,1,i) =  1.0d+00
    48: |              >                       + dt * 2.0d+00 * (   tx1 * dx1
    49: |              >                                          + ty1 * dy1
    50: |              >                                          + tz1 * dz1 )
    51: |                        d(1,2,i) =  0.0d+00
    52: |                        d(1,3,i) =  0.0d+00
    53: |                        d(1,4,i) =  0.0d+00
    54: |                        d(1,5,i) =  0.0d+00
    55: |         
    56: |                        d(2,1,i) = -dt * 2.0d+00
    57: |              >          * (  tx1 * r43 + ty1 + tz1  )
    58: |              >          * c34 * tmp2 * u(2,i,j,k)
    59: |       F                d(2,2,i) =  1.0d+00
    60: |              >          + dt * 2.0d+00 * c34 * tmp1
    61: |              >          * (  tx1 * r43 + ty1 + tz1 )
    62: |              >          + dt * 2.0d+00 * (   tx1 * dx2
    63: |              >                             + ty1 * dy2
    64: |              >                             + tz1 * dz2  )
    65: |                        d(2,3,i) = 0.0d+00
    66: |                        d(2,4,i) = 0.0d+00
    67: |                        d(2,5,i) = 0.0d+00
    68: |         
    69: |                        d(3,1,i) = -dt * 2.0d+00
    70: |              >           * (  tx1 + ty1 * r43 + tz1  )
    71: |              >           * c34 * tmp2 * u(3,i,j,k)
    72: |                        d(3,2,i) = 0.0d+00
    73: |       F                d(3,3,i) = 1.0d+00
    74: |              >         + dt * 2.0d+00 * c34 * tmp1
    75: |              >              * (  tx1 + ty1 * r43 + tz1 )
    76: |              >         + dt * 2.0d+00 * (  tx1 * dx3
    77: |              >                           + ty1 * dy3
    78: |              >                           + tz1 * dz3 )
    79: |                        d(3,4,i) = 0.0d+00
    80: |                        d(3,5,i) = 0.0d+00
    81: |         
    82: |                        d(4,1,i) = -dt * 2.0d+00
    83: |              >           * (  tx1 + ty1 + tz1 * r43  )
    84: |              >           * c34 * tmp2 * u(4,i,j,k)
    85: |                        d(4,2,i) = 0.0d+00
    86: |                        d(4,3,i) = 0.0d+00
    87: |       F                d(4,4,i) = 1.0d+00
    88: |              >         + dt * 2.0d+00 * c34 * tmp1
    89: |              >              * (  tx1 + ty1 + tz1 * r43 )
    90: |              >         + dt * 2.0d+00 * (  tx1 * dx4
    91: |              >                           + ty1 * dy4
    92: |              >                           + tz1 * dz4 )
    93: |                        d(4,5,i) = 0.0d+00
    94: |         
    95: |       F                d(5,1,i) = -dt * 2.0d+00
    96: |              >  * ( ( ( tx1 * ( r43*c34 - c1345 )
    97: |              >     + ty1 * ( c34 - c1345 )
    98: |              >     + tz1 * ( c34 - c1345 ) ) * ( u(2,i,j,k) ** 2 )
    99: |              >   + ( tx1 * ( c34 - c1345 )
   100: |              >     + ty1 * ( r43*c34 - c1345 )
   101: |              >     + tz1 * ( c34 - c1345 ) ) * ( u(3,i,j,k) ** 2 )
   102: |              >   + ( tx1 * ( c34 - c1345 )
   103: |              >     + ty1 * ( c34 - c1345 )
   104: |              >     + tz1 * ( r43*c34 - c1345 ) ) * ( u(4,i,j,k) ** 2 )
   105: |              >      ) * tmp3
   106: |              >   + ( tx1 + ty1 + tz1 ) * c1345 * tmp2 * u(5,i,j,k) )
   107: |         
   108: |                        d(5,2,i) = dt * 2.0d+00 * tmp2 * u(2,i,j,k)
   109: |              > * ( tx1 * ( r43*c34 - c1345 )
   110: |              >   + ty1 * (     c34 - c1345 )
   111: |              >   + tz1 * (     c34 - c1345 ) )
   112: |                        d(5,3,i) = dt * 2.0d+00 * tmp2 * u(3,i,j,k)
   113: |              > * ( tx1 * ( c34 - c1345 )
   114: |              >   + ty1 * ( r43*c34 -c1345 )
   115: |              >   + tz1 * ( c34 - c1345 ) )
   116: |                        d(5,4,i) = dt * 2.0d+00 * tmp2 * u(4,i,j,k)
   117: |              > * ( tx1 * ( c34 - c1345 )
   118: |              >   + ty1 * ( c34 - c1345 )
   119: |              >   + tz1 * ( r43*c34 - c1345 ) )
   120: |       F                d(5,5,i) = 1.0d+00
   121: |              >   + dt * 2.0d+00 * ( tx1  + ty1 + tz1 ) * c1345 * tmp1
   122: |              >   + dt * 2.0d+00 * (  tx1 * dx5
   123: |              >                    +  ty1 * dy5
   124: |              >                    +  tz1 * dz5 )
   125: |         
   126: |         c---------------------------------------------------------------------
   127: |         c   form the first block sub-diagonal
   128: |         c---------------------------------------------------------------------
   129: |                        tmp1 = rho_i(i,j,k-1)
   130: |                        tmp2 = tmp1 * tmp1
   131: |                        tmp3 = tmp1 * tmp2
   132: |         
   133: |                        a(1,1,i) = - dt * tz1 * dz1
   134: |                        a(1,2,i) =   0.0d+00
   135: |                        a(1,3,i) =   0.0d+00
   136: |                        a(1,4,i) = - dt * tz2
   137: |                        a(1,5,i) =   0.0d+00
   138: |         
   139: |       F                a(2,1,i) = - dt * tz2
   140: |              >           * ( - ( u(2,i,j,k-1)*u(4,i,j,k-1) ) * tmp2 )
   141: |              >           - dt * tz1 * ( - c34 * tmp2 * u(2,i,j,k-1) )
   142: |       F                a(2,2,i) = - dt * tz2 * ( u(4,i,j,k-1) * tmp1 )
   143: |              >           - dt * tz1 * c34 * tmp1
   144: |              >           - dt * tz1 * dz2
   145: |                        a(2,3,i) = 0.0d+00
   146: |                        a(2,4,i) = - dt * tz2 * ( u(2,i,j,k-1) * tmp1 )
   147: |                        a(2,5,i) = 0.0d+00
   148: |         
   149: |       F                a(3,1,i) = - dt * tz2
   150: |              >           * ( - ( u(3,i,j,k-1)*u(4,i,j,k-1) ) * tmp2 )
   151: |              >           - dt * tz1 * ( - c34 * tmp2 * u(3,i,j,k-1) )
   152: |                        a(3,2,i) = 0.0d+00
   153: |                        a(3,3,i) = - dt * tz2 * ( u(4,i,j,k-1) * tmp1 )
   154: |              >           - dt * tz1 * ( c34 * tmp1 )
   155: |              >           - dt * tz1 * dz3
   156: |                        a(3,4,i) = - dt * tz2 * ( u(3,i,j,k-1) * tmp1 )
   157: |                        a(3,5,i) = 0.0d+00
   158: |         
   159: |       F                a(4,1,i) = - dt * tz2
   160: |              >        * ( - ( u(4,i,j,k-1) * tmp1 ) ** 2
   161: |              >            + c2 * qs(i,j,k-1) * tmp1 )
   162: |              >        - dt * tz1 * ( - r43 * c34 * tmp2 * u(4,i,j,k-1) )
   163: |                        a(4,2,i) = - dt * tz2
   164: |              >             * ( - c2 * ( u(2,i,j,k-1) * tmp1 ) )
   165: |                        a(4,3,i) = - dt * tz2
   166: |              >             * ( - c2 * ( u(3,i,j,k-1) * tmp1 ) )
   167: |       F                a(4,4,i) = - dt * tz2 * ( 2.0d+00 - c2 )
   168: |              >             * ( u(4,i,j,k-1) * tmp1 )
   169: |              >             - dt * tz1 * ( r43 * c34 * tmp1 )
   170: |              >             - dt * tz1 * dz4
   171: |                        a(4,5,i) = - dt * tz2 * c2
   172: |         
   173: |       F                a(5,1,i) = - dt * tz2
   174: |              >       * ( ( c2 * 2.0d0 * qs(i,j,k-1)
   175: |              >       - c1 * u(5,i,j,k-1) )
   176: |              >            * u(4,i,j,k-1) * tmp2 )
   177: |              >       - dt * tz1
   178: |              >       * ( - ( c34 - c1345 ) * tmp3 * (u(2,i,j,k-1)**2)
   179: |              >           - ( c34 - c1345 ) * tmp3 * (u(3,i,j,k-1)**2)
   180: |              >           - ( r43*c34 - c1345 )* tmp3 * (u(4,i,j,k-1)**2)
   181: |              >          - c1345 * tmp2 * u(5,i,j,k-1) )
   182: |       F                a(5,2,i) = - dt * tz2
   183: |              >       * ( - c2 * ( u(2,i,j,k-1)*u(4,i,j,k-1) ) * tmp2 )
   184: |              >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * u(2,i,j,k-1)
   185: |       F                a(5,3,i) = - dt * tz2
   186: |              >       * ( - c2 * ( u(3,i,j,k-1)*u(4,i,j,k-1) ) * tmp2 )
   187: |              >       - dt * tz1 * ( c34 - c1345 ) * tmp2 * u(3,i,j,k-1)
   188: |       F                a(5,4,i) = - dt * tz2
   189: |              >       * ( c1 * ( u(5,i,j,k-1) * tmp1 )
   190: |              >       - c2
   191: |              >       * ( qs(i,j,k-1) * tmp1
   192: |              >            + u(4,i,j,k-1)*u(4,i,j,k-1) * tmp2 ) )
   193: |              >       - dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u(4,i,j,k-1)
   194: |       F                a(5,5,i) = - dt * tz2
   195: |              >       * ( c1 * ( u(4,i,j,k-1) * tmp1 ) )
   196: |              >       - dt * tz1 * c1345 * tmp1
   197: |              >       - dt * tz1 * dz5
   198: |         
   199: |         c---------------------------------------------------------------------
   200: |         c   form the second block sub-diagonal
   201: |         c---------------------------------------------------------------------
   202: |                        tmp1 = rho_i(i,j-1,k)
   203: |                        tmp2 = tmp1 * tmp1
   204: |                        tmp3 = tmp1 * tmp2
   205: |         
   206: |                        b(1,1,i) = - dt * ty1 * dy1
   207: |                        b(1,2,i) =   0.0d+00
   208: |                        b(1,3,i) = - dt * ty2
   209: |                        b(1,4,i) =   0.0d+00
   210: |                        b(1,5,i) =   0.0d+00
   211: |         
   212: |       F                b(2,1,i) = - dt * ty2
   213: |              >           * ( - ( u(2,i,j-1,k)*u(3,i,j-1,k) ) * tmp2 )
   214: |              >           - dt * ty1 * ( - c34 * tmp2 * u(2,i,j-1,k) )
   215: |       F                b(2,2,i) = - dt * ty2 * ( u(3,i,j-1,k) * tmp1 )
   216: |              >          - dt * ty1 * ( c34 * tmp1 )
   217: |              >          - dt * ty1 * dy2
   218: |                        b(2,3,i) = - dt * ty2 * ( u(2,i,j-1,k) * tmp1 )
   219: |                        b(2,4,i) = 0.0d+00
   220: |                        b(2,5,i) = 0.0d+00
   221: |         
   222: |       F                b(3,1,i) = - dt * ty2
   223: |              >           * ( - ( u(3,i,j-1,k) * tmp1 ) ** 2
   224: |              >       + c2 * ( qs(i,j-1,k) * tmp1 ) )
   225: |              >       - dt * ty1 * ( - r43 * c34 * tmp2 * u(3,i,j-1,k) )
   226: |                        b(3,2,i) = - dt * ty2
   227: |              >                   * ( - c2 * ( u(2,i,j-1,k) * tmp1 ) )
   228: |       F                b(3,3,i) = - dt * ty2 * ( ( 2.0d+00 - c2 )
   229: |              >                   * ( u(3,i,j-1,k) * tmp1 ) )
   230: |              >       - dt * ty1 * ( r43 * c34 * tmp1 )
   231: |              >       - dt * ty1 * dy3
   232: |                        b(3,4,i) = - dt * ty2
   233: |              >                   * ( - c2 * ( u(4,i,j-1,k) * tmp1 ) )
   234: |                        b(3,5,i) = - dt * ty2 * c2
   235: |         
   236: |       F                b(4,1,i) = - dt * ty2
   237: |              >              * ( - ( u(3,i,j-1,k)*u(4,i,j-1,k) ) * tmp2 )
   238: |              >       - dt * ty1 * ( - c34 * tmp2 * u(4,i,j-1,k) )
   239: |                        b(4,2,i) = 0.0d+00
   240: |                        b(4,3,i) = - dt * ty2 * ( u(4,i,j-1,k) * tmp1 )
   241: |                        b(4,4,i) = - dt * ty2 * ( u(3,i,j-1,k) * tmp1 )
   242: |              >                        - dt * ty1 * ( c34 * tmp1 )
   243: |              >                        - dt * ty1 * dy4
   244: |                        b(4,5,i) = 0.0d+00
   245: |         
   246: |       F                b(5,1,i) = - dt * ty2
   247: |              >          * ( ( c2 * 2.0d0 * qs(i,j-1,k)
   248: |              >               - c1 * u(5,i,j-1,k) )
   249: |              >          * ( u(3,i,j-1,k) * tmp2 ) )
   250: |              >          - dt * ty1
   251: |              >          * ( - (     c34 - c1345 )*tmp3*(u(2,i,j-1,k)**2)
   252: |              >              - ( r43*c34 - c1345 )*tmp3*(u(3,i,j-1,k)**2)
   253: |              >              - (     c34 - c1345 )*tmp3*(u(4,i,j-1,k)**2)
   254: |              >              - c1345*tmp2*u(5,i,j-1,k) )
   255: |       F                b(5,2,i) = - dt * ty2
   256: |              >          * ( - c2 * ( u(2,i,j-1,k)*u(3,i,j-1,k) ) * tmp2 )
   257: |              >          - dt * ty1
   258: |              >          * ( c34 - c1345 ) * tmp2 * u(2,i,j-1,k)
   259: |       F                b(5,3,i) = - dt * ty2
   260: |              >          * ( c1 * ( u(5,i,j-1,k) * tmp1 )
   261: |              >          - c2
   262: |              >          * ( qs(i,j-1,k) * tmp1
   263: |              >               + u(3,i,j-1,k)*u(3,i,j-1,k) * tmp2 ) )
   264: |              >          - dt * ty1
   265: |              >          * ( r43*c34 - c1345 ) * tmp2 * u(3,i,j-1,k)
   266: |       F                b(5,4,i) = - dt * ty2
   267: |              >          * ( - c2 * ( u(3,i,j-1,k)*u(4,i,j-1,k) ) * tmp2 )
   268: |              >          - dt * ty1 * ( c34 - c1345 ) * tmp2 * u(4,i,j-1,k)
   269: |       F                b(5,5,i) = - dt * ty2
   270: |              >          * ( c1 * ( u(3,i,j-1,k) * tmp1 ) )
   271: |              >          - dt * ty1 * c1345 * tmp1
   272: |              >          - dt * ty1 * dy5
   273: |         
   274: |         c---------------------------------------------------------------------
   275: |         c   form the third block sub-diagonal
   276: |         c---------------------------------------------------------------------
   277: |                        tmp1 = rho_i(i-1,j,k)
   278: |                        tmp2 = tmp1 * tmp1
   279: |                        tmp3 = tmp1 * tmp2
   280: |         
   281: |                        c(1,1,i) = - dt * tx1 * dx1
   282: |                        c(1,2,i) = - dt * tx2
   283: |                        c(1,3,i) =   0.0d+00
   284: |                        c(1,4,i) =   0.0d+00
   285: |                        c(1,5,i) =   0.0d+00
   286: |         
   287: |       F                c(2,1,i) = - dt * tx2
   288: |              >          * ( - ( u(2,i-1,j,k) * tmp1 ) ** 2
   289: |              >       + c2 * qs(i-1,j,k) * tmp1 )
   290: |              >          - dt * tx1 * ( - r43 * c34 * tmp2 * u(2,i-1,j,k) )
   291: |       F                c(2,2,i) = - dt * tx2
   292: |              >          * ( ( 2.0d+00 - c2 ) * ( u(2,i-1,j,k) * tmp1 ) )
   293: |              >          - dt * tx1 * ( r43 * c34 * tmp1 )
   294: |              >          - dt * tx1 * dx2
   295: |                        c(2,3,i) = - dt * tx2
   296: |              >              * ( - c2 * ( u(3,i-1,j,k) * tmp1 ) )
   297: |                        c(2,4,i) = - dt * tx2
   298: |              >              * ( - c2 * ( u(4,i-1,j,k) * tmp1 ) )
   299: |                        c(2,5,i) = - dt * tx2 * c2
   300: |         
   301: |       F                c(3,1,i) = - dt * tx2
   302: |              >              * ( - ( u(2,i-1,j,k) * u(3,i-1,j,k) ) * tmp2 )
   303: |              >         - dt * tx1 * ( - c34 * tmp2 * u(3,i-1,j,k) )
   304: |                        c(3,2,i) = - dt * tx2 * ( u(3,i-1,j,k) * tmp1 )
   305: |       F                c(3,3,i) = - dt * tx2 * ( u(2,i-1,j,k) * tmp1 )
   306: |              >          - dt * tx1 * ( c34 * tmp1 )
   307: |              >          - dt * tx1 * dx3
   308: |                        c(3,4,i) = 0.0d+00
   309: |                        c(3,5,i) = 0.0d+00
   310: |         
   311: |       F                c(4,1,i) = - dt * tx2
   312: |              >          * ( - ( u(2,i-1,j,k)*u(4,i-1,j,k) ) * tmp2 )
   313: |              >          - dt * tx1 * ( - c34 * tmp2 * u(4,i-1,j,k) )
   314: |                        c(4,2,i) = - dt * tx2 * ( u(4,i-1,j,k) * tmp1 )
   315: |                        c(4,3,i) = 0.0d+00
   316: |                        c(4,4,i) = - dt * tx2 * ( u(2,i-1,j,k) * tmp1 )
   317: |              >          - dt * tx1 * ( c34 * tmp1 )
   318: |              >          - dt * tx1 * dx4
   319: |                        c(4,5,i) = 0.0d+00
   320: |         
   321: |       F                c(5,1,i) = - dt * tx2
   322: |              >          * ( ( c2 * 2.0d0 * qs(i-1,j,k)
   323: |              >              - c1 * u(5,i-1,j,k) )
   324: |              >          * u(2,i-1,j,k) * tmp2 )
   325: |              >          - dt * tx1
   326: |              >          * ( - ( r43*c34 - c1345 ) * tmp3 * ( u(2,i-1,j,k)**2 )
   327: |              >              - (     c34 - c1345 ) * tmp3 * ( u(3,i-1,j,k)**2 )
   328: |              >              - (     c34 - c1345 ) * tmp3 * ( u(4,i-1,j,k)**2 )
   329: |              >              - c1345 * tmp2 * u(5,i-1,j,k) )
   330: |       F                c(5,2,i) = - dt * tx2
   331: |              >          * ( c1 * ( u(5,i-1,j,k) * tmp1 )
   332: |              >             - c2
   333: |              >             * ( u(2,i-1,j,k)*u(2,i-1,j,k) * tmp2
   334: |              >                  + qs(i-1,j,k) * tmp1 ) )
   335: |              >           - dt * tx1
   336: |              >           * ( r43*c34 - c1345 ) * tmp2 * u(2,i-1,j,k)
   337: |       F                c(5,3,i) = - dt * tx2
   338: |              >           * ( - c2 * ( u(3,i-1,j,k)*u(2,i-1,j,k) ) * tmp2 )
   339: |              >           - dt * tx1
   340: |              >           * (  c34 - c1345 ) * tmp2 * u(3,i-1,j,k)
   341: |       F                c(5,4,i) = - dt * tx2
   342: |              >           * ( - c2 * ( u(4,i-1,j,k)*u(2,i-1,j,k) ) * tmp2 )
   343: |              >           - dt * tx1
   344: |              >           * (  c34 - c1345 ) * tmp2 * u(4,i-1,j,k)
   345: |       F                c(5,5,i) = - dt * tx2
   346: |              >           * ( c1 * ( u(2,i-1,j,k) * tmp1 ) )
   347: |              >           - dt * tx1 * c1345 * tmp1
   348: |              >           - dt * tx1 * dx5
   349: |         
   350: V------               end do
   351:           
   352:           
   353:                 return
   354:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacld.f

PROCEDURE NAME: JACLD
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacld.f

PROCEDURE NAME: JACLD
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (jacld.f:38)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (jacld.f:38)
    *** The number of VLOAD, VSTORE. : 23, 100. (jacld.f:38)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:53 2021
FILE NAME: jacld.f

PROCEDURE NAME: JACLD
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
    Total size                 : 1394480 bytes
      Register spill area      :  697144 bytes
      Parameter area           :       0 bytes
      Register save area       :     176 bytes
      User data area           :  697144 bytes
      Others                   :       8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (jacld.f:38)
    *** Estimated execution cycle                       : 1895
    *** The number of VECTOR REGISTER SPILL
          Total                                         : 2
            Not enough registers                        : 2
    *** The number of VECTOR REGISTER RESTORE
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 7
            Not enough registers                        : 3
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 88
            Over basic blocks                           : 88
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END


