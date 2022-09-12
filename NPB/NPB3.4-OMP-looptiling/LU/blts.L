NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: blts.f

PROCEDURE NAME: BLTS
FORMAT LIST

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


