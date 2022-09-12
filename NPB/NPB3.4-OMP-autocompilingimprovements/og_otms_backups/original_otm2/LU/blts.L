NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME : blts.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: blts.f

PROCEDURE NAME: BLTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    50: vec( 101): Vectorized loop.
    52: opt(1037): Feedback of array elements.
    52: vec( 128): Fused multiply-add operation applied.
    59: opt(1037): Feedback of array elements.
    59: vec( 128): Fused multiply-add operation applied.
    78: opt(1418): Constant-length loop is expanded.
    96: opt(1019): Feedback of scalar value from one loop pass to another.
   108: opt(1019): Feedback of scalar value from one loop pass to another.
   120: opt(1019): Feedback of scalar value from one loop pass to another.
   132: opt(1019): Feedback of scalar value from one loop pass to another.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: blts.f

PROCEDURE NAME: BLTS
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

     4:                 subroutine blts ( ldmx, ldmy, ldmz,
     5:                >                  nx, ny, nz,
     6:                >                  omega,
     7:                >                  v,
     8:                >                  ldz, ldy, ldx, d,
     9:                >                  ist, iend, j, k )
    10:           
    11:           c---------------------------------------------------------------------
    12:           c---------------------------------------------------------------------
    13:           
    14:           c---------------------------------------------------------------------
    15:           c
    16:           c   compute the regular-sparse, block lower triangular solution:
    17:           c
    18:           c                     v <-- ( L-inv ) * v
    19:           c
    20:           c---------------------------------------------------------------------
    21:           
    22:                 implicit none
    23:           
    24:           c---------------------------------------------------------------------
    25:           c  input parameters
    26:           c---------------------------------------------------------------------
    27:                 integer ldmx, ldmy, ldmz
    28:                 integer nx, ny, nz
    29:                 double precision  omega
    30:           c---------------------------------------------------------------------
    31:           c   To improve cache performance, second two dimensions padded by 1
    32:           c   for even number sizes only.  Only needed in v.
    33:           c---------------------------------------------------------------------
    34:                 double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz),
    35:                >        ldz( 5, 5, ldmx ),
    36:                >        ldy( 5, 5, ldmx ),
    37:                >        ldx( 5, 5, ldmx ),
    38:                >        d  ( 5, 5, ldmx )
    39:                 integer ist, iend, j, k
    40:           
    41:           c---------------------------------------------------------------------
    42:           c  local variables
    43:           c---------------------------------------------------------------------
    44:                 integer i, m
    45:                 double precision  tmp, tmp1
    46:                 double precision  tmat(5,5), tv(5)
    47:           
    48:           
    49: +------>           do i = ist, iend
    50: |V----->              do m = 1, 5
    51: ||        
    52: ||      F                   tv( m ) =  v( m, i, j, k )
    53: ||             >    - omega * (  ldz( m, 1, i ) * v( 1, i, j, k-1 )
    54: ||             >               + ldz( m, 2, i ) * v( 2, i, j, k-1 )
    55: ||             >               + ldz( m, 3, i ) * v( 3, i, j, k-1 )
    56: ||             >               + ldz( m, 4, i ) * v( 4, i, j, k-1 )
    57: ||             >               + ldz( m, 5, i ) * v( 5, i, j, k-1 )  )
    58: ||        
    59: ||      F                   tv( m ) =  tv( m )
    60: ||             > - omega * ( ldy( m, 1, i ) * v( 1, i, j-1, k )
    61: ||             >           + ldx( m, 1, i ) * v( 1, i-1, j, k )
    62: ||             >           + ldy( m, 2, i ) * v( 2, i, j-1, k )
    63: ||             >           + ldx( m, 2, i ) * v( 2, i-1, j, k )
    64: ||             >           + ldy( m, 3, i ) * v( 3, i, j-1, k )
    65: ||             >           + ldx( m, 3, i ) * v( 3, i-1, j, k )
    66: ||             >           + ldy( m, 4, i ) * v( 4, i, j-1, k )
    67: ||             >           + ldx( m, 4, i ) * v( 4, i-1, j, k )
    68: ||             >           + ldy( m, 5, i ) * v( 5, i, j-1, k )
    69: ||             >           + ldx( m, 5, i ) * v( 5, i-1, j, k ) )
    70: ||        
    71: |V-----               end do
    72: |          
    73: |         c---------------------------------------------------------------------
    74: |         c   diagonal block inversion
    75: |         c
    76: |         c   forward elimination
    77: |         c---------------------------------------------------------------------
    78: |*----->              do m = 1, 5
    79: ||                       tmat( m, 1 ) = d( m, 1, i )
    80: ||                       tmat( m, 2 ) = d( m, 2, i )
    81: ||                       tmat( m, 3 ) = d( m, 3, i )
    82: ||                       tmat( m, 4 ) = d( m, 4, i )
    83: ||                       tmat( m, 5 ) = d( m, 5, i )
    84: |*-----               end do
    85: |         
    86: |                     tmp1 = 1.0d+00 / tmat( 1, 1 )
    87: |                     tmp = tmp1 * tmat( 2, 1 )
    88: |                     tmat( 2, 2 ) =  tmat( 2, 2 )
    89: |              >           - tmp * tmat( 1, 2 )
    90: |                     tmat( 2, 3 ) =  tmat( 2, 3 )
    91: |              >           - tmp * tmat( 1, 3 )
    92: |                     tmat( 2, 4 ) =  tmat( 2, 4 )
    93: |              >           - tmp * tmat( 1, 4 )
    94: |                     tmat( 2, 5 ) =  tmat( 2, 5 )
    95: |              >           - tmp * tmat( 1, 5 )
    96: |                     tv( 2 ) = tv( 2 )
    97: |              >        - tv( 1 ) * tmp
    98: |         
    99: |                     tmp = tmp1 * tmat( 3, 1 )
   100: |                     tmat( 3, 2 ) =  tmat( 3, 2 )
   101: |              >           - tmp * tmat( 1, 2 )
   102: |                     tmat( 3, 3 ) =  tmat( 3, 3 )
   103: |              >           - tmp * tmat( 1, 3 )
   104: |                     tmat( 3, 4 ) =  tmat( 3, 4 )
   105: |              >           - tmp * tmat( 1, 4 )
   106: |                     tmat( 3, 5 ) =  tmat( 3, 5 )
   107: |              >           - tmp * tmat( 1, 5 )
   108: |                     tv( 3 ) = tv( 3 )
   109: |              >        - tv( 1 ) * tmp
   110: |         
   111: |                     tmp = tmp1 * tmat( 4, 1 )
   112: |                     tmat( 4, 2 ) =  tmat( 4, 2 )
   113: |              >           - tmp * tmat( 1, 2 )
   114: |                     tmat( 4, 3 ) =  tmat( 4, 3 )
   115: |              >           - tmp * tmat( 1, 3 )
   116: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   117: |              >           - tmp * tmat( 1, 4 )
   118: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   119: |              >           - tmp * tmat( 1, 5 )
   120: |                     tv( 4 ) = tv( 4 )
   121: |              >        - tv( 1 ) * tmp
   122: |         
   123: |                     tmp = tmp1 * tmat( 5, 1 )
   124: |                     tmat( 5, 2 ) =  tmat( 5, 2 )
   125: |              >           - tmp * tmat( 1, 2 )
   126: |                     tmat( 5, 3 ) =  tmat( 5, 3 )
   127: |              >           - tmp * tmat( 1, 3 )
   128: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   129: |              >           - tmp * tmat( 1, 4 )
   130: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   131: |              >           - tmp * tmat( 1, 5 )
   132: |                     tv( 5 ) = tv( 5 )
   133: |              >        - tv( 1 ) * tmp
   134: |         
   135: |         
   136: |         
   137: |                     tmp1 = 1.0d+00 / tmat( 2, 2 )
   138: |                     tmp = tmp1 * tmat( 3, 2 )
   139: |                     tmat( 3, 3 ) =  tmat( 3, 3 )
   140: |              >           - tmp * tmat( 2, 3 )
   141: |                     tmat( 3, 4 ) =  tmat( 3, 4 )
   142: |              >           - tmp * tmat( 2, 4 )
   143: |                     tmat( 3, 5 ) =  tmat( 3, 5 )
   144: |              >           - tmp * tmat( 2, 5 )
   145: |                     tv( 3 ) = tv( 3 )
   146: |              >        - tv( 2 ) * tmp
   147: |         
   148: |                     tmp = tmp1 * tmat( 4, 2 )
   149: |                     tmat( 4, 3 ) =  tmat( 4, 3 )
   150: |              >           - tmp * tmat( 2, 3 )
   151: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   152: |              >           - tmp * tmat( 2, 4 )
   153: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   154: |              >           - tmp * tmat( 2, 5 )
   155: |                     tv( 4 ) = tv( 4 )
   156: |              >        - tv( 2 ) * tmp
   157: |         
   158: |                     tmp = tmp1 * tmat( 5, 2 )
   159: |                     tmat( 5, 3 ) =  tmat( 5, 3 )
   160: |              >           - tmp * tmat( 2, 3 )
   161: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   162: |              >           - tmp * tmat( 2, 4 )
   163: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   164: |              >           - tmp * tmat( 2, 5 )
   165: |                     tv( 5 ) = tv( 5 )
   166: |              >        - tv( 2 ) * tmp
   167: |         
   168: |         
   169: |         
   170: |                     tmp1 = 1.0d+00 / tmat( 3, 3 )
   171: |                     tmp = tmp1 * tmat( 4, 3 )
   172: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   173: |              >           - tmp * tmat( 3, 4 )
   174: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   175: |              >           - tmp * tmat( 3, 5 )
   176: |                     tv( 4 ) = tv( 4 )
   177: |              >        - tv( 3 ) * tmp
   178: |         
   179: |                     tmp = tmp1 * tmat( 5, 3 )
   180: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   181: |              >           - tmp * tmat( 3, 4 )
   182: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   183: |              >           - tmp * tmat( 3, 5 )
   184: |                     tv( 5 ) = tv( 5 )
   185: |              >        - tv( 3 ) * tmp
   186: |         
   187: |         
   188: |         
   189: |                     tmp1 = 1.0d+00 / tmat( 4, 4 )
   190: |                     tmp = tmp1 * tmat( 5, 4 )
   191: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   192: |              >           - tmp * tmat( 4, 5 )
   193: |                     tv( 5 ) = tv( 5 )
   194: |              >        - tv( 4 ) * tmp
   195: |         
   196: |         c---------------------------------------------------------------------
   197: |         c   back substitution
   198: |         c---------------------------------------------------------------------
   199: |                     v( 5, i, j, k ) = tv( 5 )
   200: |              >                      / tmat( 5, 5 )
   201: |         
   202: |                     tv( 4 ) = tv( 4 )
   203: |              >           - tmat( 4, 5 ) * v( 5, i, j, k )
   204: |                     v( 4, i, j, k ) = tv( 4 )
   205: |              >                      / tmat( 4, 4 )
   206: |         
   207: |                     tv( 3 ) = tv( 3 )
   208: |              >           - tmat( 3, 4 ) * v( 4, i, j, k )
   209: |              >           - tmat( 3, 5 ) * v( 5, i, j, k )
   210: |                     v( 3, i, j, k ) = tv( 3 )
   211: |              >                      / tmat( 3, 3 )
   212: |         
   213: |                     tv( 2 ) = tv( 2 )
   214: |              >           - tmat( 2, 3 ) * v( 3, i, j, k )
   215: |              >           - tmat( 2, 4 ) * v( 4, i, j, k )
   216: |              >           - tmat( 2, 5 ) * v( 5, i, j, k )
   217: |                     v( 2, i, j, k ) = tv( 2 )
   218: |              >                      / tmat( 2, 2 )
   219: |         
   220: |                     tv( 1 ) = tv( 1 )
   221: |              >           - tmat( 1, 2 ) * v( 2, i, j, k )
   222: |              >           - tmat( 1, 3 ) * v( 3, i, j, k )
   223: |              >           - tmat( 1, 4 ) * v( 4, i, j, k )
   224: |              >           - tmat( 1, 5 ) * v( 5, i, j, k )
   225: |                     v( 1, i, j, k ) = tv( 1 )
   226: |              >                      / tmat( 1, 1 )
   227: |         
   228: |         
   229: +------           enddo
   230:           
   231:           
   232:                 return
   233:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: blts.f

PROCEDURE NAME: BLTS
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (blts.f:49)
    <Unvectorized loop.>

    LOOP BEGIN: (blts.f:50)
      <Vectorized loop.>
      **  Runtime-tested short loop. (blts.f:50)
      *** The number of VGT,   VSC.    :  0,  0. (blts.f:50)
      *** The number of VLOAD, VSTORE. : 16,  1. (blts.f:50)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: blts.f

PROCEDURE NAME: BLTS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 58 [s0-s11 s18-s63]
      Vector registers         : 33 [v31-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 432 bytes
      Register spill area      :  96 bytes
      Parameter area           :   0 bytes
      Register save area       :   0 bytes
      User data area           : 336 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (blts.f:49)
    *** Estimated execution cycle                       : 1082
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 12
            Not enough registers                        :  8
            Over basic blocks                           :  4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 15
            Over basic blocks                           : 15
    *** The number of SCALAR REGISTER TRANSFER          : 3
  LOOP END


