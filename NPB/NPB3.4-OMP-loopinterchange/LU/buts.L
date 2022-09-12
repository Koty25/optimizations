NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:40 2020
FILE NAME: buts.f

PROCEDURE NAME: BUTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    50: vec( 101): Vectorized loop.
    51: opt(1037): Feedback of array elements.
    51: vec( 128): Fused multiply-add operation applied.
    58: opt(1037): Feedback of array elements.
    58: vec( 128): Fused multiply-add operation applied.
    74: opt(1418): Constant-length loop is expanded.
    92: opt(1019): Feedback of scalar value from one loop pass to another.
   104: opt(1019): Feedback of scalar value from one loop pass to another.
   116: opt(1019): Feedback of scalar value from one loop pass to another.
   128: opt(1019): Feedback of scalar value from one loop pass to another.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:40 2020
FILE NAME: buts.f

PROCEDURE NAME: BUTS
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine buts( ldmx, ldmy, ldmz,
     5:                >                 nx, ny, nz,
     6:                >                 omega,
     7:                >                 v,
     8:                >                 d, udx, udy, udz,
     9:                >                 ist, iend, j, k )
    10:           
    11:           c---------------------------------------------------------------------
    12:           c---------------------------------------------------------------------
    13:           
    14:           c---------------------------------------------------------------------
    15:           c
    16:           c   compute the regular-sparse, block upper triangular solution:
    17:           c
    18:           c                     v <-- ( U-inv ) * v
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
    34:                 double precision  v( 5,ldmx/2*2+1, ldmy/2*2+1, ldmz),
    35:                >        d  ( 5, 5, ldmx ),
    36:                >        udx( 5, 5, ldmx ),
    37:                >        udy( 5, 5, ldmx ),
    38:                >        udz( 5, 5, ldmx )
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
    49: +------>           do i = iend, ist, -1
    50: |V----->              do m = 1, 5
    51: ||      F                   tv( m ) =
    52: ||             >      omega * (  udz( m, 1, i ) * v( 1, i, j, k+1 )
    53: ||             >               + udz( m, 2, i ) * v( 2, i, j, k+1 )
    54: ||             >               + udz( m, 3, i ) * v( 3, i, j, k+1 )
    55: ||             >               + udz( m, 4, i ) * v( 4, i, j, k+1 )
    56: ||             >               + udz( m, 5, i ) * v( 5, i, j, k+1 ) )
    57: ||        
    58: ||      F                   tv( m ) = tv( m )
    59: ||             > + omega * ( udy( m, 1, i ) * v( 1, i, j+1, k )
    60: ||             >           + udx( m, 1, i ) * v( 1, i+1, j, k )
    61: ||             >           + udy( m, 2, i ) * v( 2, i, j+1, k )
    62: ||             >           + udx( m, 2, i ) * v( 2, i+1, j, k )
    63: ||             >           + udy( m, 3, i ) * v( 3, i, j+1, k )
    64: ||             >           + udx( m, 3, i ) * v( 3, i+1, j, k )
    65: ||             >           + udy( m, 4, i ) * v( 4, i, j+1, k )
    66: ||             >           + udx( m, 4, i ) * v( 4, i+1, j, k )
    67: ||             >           + udy( m, 5, i ) * v( 5, i, j+1, k )
    68: ||             >           + udx( m, 5, i ) * v( 5, i+1, j, k ) )
    69: |V-----               end do
    70: |         
    71: |         c---------------------------------------------------------------------
    72: |         c   diagonal block inversion
    73: |         c---------------------------------------------------------------------
    74: |*----->              do m = 1, 5
    75: ||                       tmat( m, 1 ) = d( m, 1, i )
    76: ||                       tmat( m, 2 ) = d( m, 2, i )
    77: ||                       tmat( m, 3 ) = d( m, 3, i )
    78: ||                       tmat( m, 4 ) = d( m, 4, i )
    79: ||                       tmat( m, 5 ) = d( m, 5, i )
    80: |*-----               end do
    81: |         
    82: |                     tmp1 = 1.0d+00 / tmat( 1, 1 )
    83: |                     tmp = tmp1 * tmat( 2, 1 )
    84: |                     tmat( 2, 2 ) =  tmat( 2, 2 )
    85: |              >           - tmp * tmat( 1, 2 )
    86: |                     tmat( 2, 3 ) =  tmat( 2, 3 )
    87: |              >           - tmp * tmat( 1, 3 )
    88: |                     tmat( 2, 4 ) =  tmat( 2, 4 )
    89: |              >           - tmp * tmat( 1, 4 )
    90: |                     tmat( 2, 5 ) =  tmat( 2, 5 )
    91: |              >           - tmp * tmat( 1, 5 )
    92: |                     tv( 2 ) = tv( 2 )
    93: |              >        - tv( 1 ) * tmp
    94: |         
    95: |                     tmp = tmp1 * tmat( 3, 1 )
    96: |                     tmat( 3, 2 ) =  tmat( 3, 2 )
    97: |              >           - tmp * tmat( 1, 2 )
    98: |                     tmat( 3, 3 ) =  tmat( 3, 3 )
    99: |              >           - tmp * tmat( 1, 3 )
   100: |                     tmat( 3, 4 ) =  tmat( 3, 4 )
   101: |              >           - tmp * tmat( 1, 4 )
   102: |                     tmat( 3, 5 ) =  tmat( 3, 5 )
   103: |              >           - tmp * tmat( 1, 5 )
   104: |                     tv( 3 ) = tv( 3 )
   105: |              >        - tv( 1 ) * tmp
   106: |         
   107: |                     tmp = tmp1 * tmat( 4, 1 )
   108: |                     tmat( 4, 2 ) =  tmat( 4, 2 )
   109: |              >           - tmp * tmat( 1, 2 )
   110: |                     tmat( 4, 3 ) =  tmat( 4, 3 )
   111: |              >           - tmp * tmat( 1, 3 )
   112: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   113: |              >           - tmp * tmat( 1, 4 )
   114: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   115: |              >           - tmp * tmat( 1, 5 )
   116: |                     tv( 4 ) = tv( 4 )
   117: |              >        - tv( 1 ) * tmp
   118: |         
   119: |                     tmp = tmp1 * tmat( 5, 1 )
   120: |                     tmat( 5, 2 ) =  tmat( 5, 2 )
   121: |              >           - tmp * tmat( 1, 2 )
   122: |                     tmat( 5, 3 ) =  tmat( 5, 3 )
   123: |              >           - tmp * tmat( 1, 3 )
   124: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   125: |              >           - tmp * tmat( 1, 4 )
   126: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   127: |              >           - tmp * tmat( 1, 5 )
   128: |                     tv( 5 ) = tv( 5 )
   129: |              >        - tv( 1 ) * tmp
   130: |         
   131: |         
   132: |         
   133: |                     tmp1 = 1.0d+00 / tmat( 2, 2 )
   134: |                     tmp = tmp1 * tmat( 3, 2 )
   135: |                     tmat( 3, 3 ) =  tmat( 3, 3 )
   136: |              >           - tmp * tmat( 2, 3 )
   137: |                     tmat( 3, 4 ) =  tmat( 3, 4 )
   138: |              >           - tmp * tmat( 2, 4 )
   139: |                     tmat( 3, 5 ) =  tmat( 3, 5 )
   140: |              >           - tmp * tmat( 2, 5 )
   141: |                     tv( 3 ) = tv( 3 )
   142: |              >        - tv( 2 ) * tmp
   143: |         
   144: |                     tmp = tmp1 * tmat( 4, 2 )
   145: |                     tmat( 4, 3 ) =  tmat( 4, 3 )
   146: |              >           - tmp * tmat( 2, 3 )
   147: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   148: |              >           - tmp * tmat( 2, 4 )
   149: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   150: |              >           - tmp * tmat( 2, 5 )
   151: |                     tv( 4 ) = tv( 4 )
   152: |              >        - tv( 2 ) * tmp
   153: |         
   154: |                     tmp = tmp1 * tmat( 5, 2 )
   155: |                     tmat( 5, 3 ) =  tmat( 5, 3 )
   156: |              >           - tmp * tmat( 2, 3 )
   157: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   158: |              >           - tmp * tmat( 2, 4 )
   159: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   160: |              >           - tmp * tmat( 2, 5 )
   161: |                     tv( 5 ) = tv( 5 )
   162: |              >        - tv( 2 ) * tmp
   163: |         
   164: |         
   165: |         
   166: |                     tmp1 = 1.0d+00 / tmat( 3, 3 )
   167: |                     tmp = tmp1 * tmat( 4, 3 )
   168: |                     tmat( 4, 4 ) =  tmat( 4, 4 )
   169: |              >           - tmp * tmat( 3, 4 )
   170: |                     tmat( 4, 5 ) =  tmat( 4, 5 )
   171: |              >           - tmp * tmat( 3, 5 )
   172: |                     tv( 4 ) = tv( 4 )
   173: |              >        - tv( 3 ) * tmp
   174: |         
   175: |                     tmp = tmp1 * tmat( 5, 3 )
   176: |                     tmat( 5, 4 ) =  tmat( 5, 4 )
   177: |              >           - tmp * tmat( 3, 4 )
   178: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   179: |              >           - tmp * tmat( 3, 5 )
   180: |                     tv( 5 ) = tv( 5 )
   181: |              >        - tv( 3 ) * tmp
   182: |         
   183: |         
   184: |         
   185: |                     tmp1 = 1.0d+00 / tmat( 4, 4 )
   186: |                     tmp = tmp1 * tmat( 5, 4 )
   187: |                     tmat( 5, 5 ) =  tmat( 5, 5 )
   188: |              >           - tmp * tmat( 4, 5 )
   189: |                     tv( 5 ) = tv( 5 )
   190: |              >        - tv( 4 ) * tmp
   191: |         
   192: |         c---------------------------------------------------------------------
   193: |         c   back substitution
   194: |         c---------------------------------------------------------------------
   195: |                     tv( 5 ) = tv( 5 )
   196: |              >                      / tmat( 5, 5 )
   197: |         
   198: |                     tv( 4 ) = tv( 4 )
   199: |              >           - tmat( 4, 5 ) * tv( 5 )
   200: |                     tv( 4 ) = tv( 4 )
   201: |              >                      / tmat( 4, 4 )
   202: |         
   203: |                     tv( 3 ) = tv( 3 )
   204: |              >           - tmat( 3, 4 ) * tv( 4 )
   205: |              >           - tmat( 3, 5 ) * tv( 5 )
   206: |                     tv( 3 ) = tv( 3 )
   207: |              >                      / tmat( 3, 3 )
   208: |         
   209: |                     tv( 2 ) = tv( 2 )
   210: |              >           - tmat( 2, 3 ) * tv( 3 )
   211: |              >           - tmat( 2, 4 ) * tv( 4 )
   212: |              >           - tmat( 2, 5 ) * tv( 5 )
   213: |                     tv( 2 ) = tv( 2 )
   214: |              >                      / tmat( 2, 2 )
   215: |         
   216: |                     tv( 1 ) = tv( 1 )
   217: |              >           - tmat( 1, 2 ) * tv( 2 )
   218: |              >           - tmat( 1, 3 ) * tv( 3 )
   219: |              >           - tmat( 1, 4 ) * tv( 4 )
   220: |              >           - tmat( 1, 5 ) * tv( 5 )
   221: |                     tv( 1 ) = tv( 1 )
   222: |              >                      / tmat( 1, 1 )
   223: |         
   224: |                     v( 1, i, j, k ) = v( 1, i, j, k ) - tv( 1 )
   225: |                     v( 2, i, j, k ) = v( 2, i, j, k ) - tv( 2 )
   226: |                     v( 3, i, j, k ) = v( 3, i, j, k ) - tv( 3 )
   227: |                     v( 4, i, j, k ) = v( 4, i, j, k ) - tv( 4 )
   228: |                     v( 5, i, j, k ) = v( 5, i, j, k ) - tv( 5 )
   229: |         
   230: +------           enddo
   231:           
   232:            
   233:                 return
   234:                 end


