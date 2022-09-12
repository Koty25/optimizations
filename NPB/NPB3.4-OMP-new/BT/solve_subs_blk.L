NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: MATVEC_SUB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    27: vec( 101): Vectorized loop.
    28: vec( 128): Fused multiply-add operation applied.
    33: vec( 128): Fused multiply-add operation applied.
    38: vec( 128): Fused multiply-add operation applied.
    43: vec( 128): Fused multiply-add operation applied.
    48: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: MATVEC_SUB
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine matvec_sub(ablock,avec,bvec)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c     subtracts bvec=bvec - ablock*avec
    12:           c---------------------------------------------------------------------
    13:           
    14:                 implicit none
    15:                 include 'blk_par.h'
    16:           
    17:                 double precision ablock,avec,bvec
    18:                 dimension ablock(blkdim,5,5),avec(blkdim,5),bvec(blkdim,5)
    19:           
    20:                 integer i
    21:           
    22:           c---------------------------------------------------------------------
    23:           c            rhs(i,ic,jc,kc) = rhs(i,ic,jc,kc)
    24:           c     $           - lhs(i,1,ablock,ia)*
    25:           c---------------------------------------------------------------------
    26:           !dir$ vector always
    27: V------>        do i = 1, bsize
    28: |       F          bvec(i,1) = bvec(i,1) - ablock(i,1,1)*avec(i,1)
    29: |              >                     - ablock(i,1,2)*avec(i,2)
    30: |              >                     - ablock(i,1,3)*avec(i,3)
    31: |              >                     - ablock(i,1,4)*avec(i,4)
    32: |              >                     - ablock(i,1,5)*avec(i,5)
    33: |       F          bvec(i,2) = bvec(i,2) - ablock(i,2,1)*avec(i,1)
    34: |              >                     - ablock(i,2,2)*avec(i,2)
    35: |              >                     - ablock(i,2,3)*avec(i,3)
    36: |              >                     - ablock(i,2,4)*avec(i,4)
    37: |              >                     - ablock(i,2,5)*avec(i,5)
    38: |       F          bvec(i,3) = bvec(i,3) - ablock(i,3,1)*avec(i,1)
    39: |              >                     - ablock(i,3,2)*avec(i,2)
    40: |              >                     - ablock(i,3,3)*avec(i,3)
    41: |              >                     - ablock(i,3,4)*avec(i,4)
    42: |              >                     - ablock(i,3,5)*avec(i,5)
    43: |       F          bvec(i,4) = bvec(i,4) - ablock(i,4,1)*avec(i,1)
    44: |              >                     - ablock(i,4,2)*avec(i,2)
    45: |              >                     - ablock(i,4,3)*avec(i,3)
    46: |              >                     - ablock(i,4,4)*avec(i,4)
    47: |              >                     - ablock(i,4,5)*avec(i,5)
    48: |       F          bvec(i,5) = bvec(i,5) - ablock(i,5,1)*avec(i,1)
    49: |              >                     - ablock(i,5,2)*avec(i,2)
    50: |              >                     - ablock(i,5,3)*avec(i,3)
    51: |              >                     - ablock(i,5,4)*avec(i,4)
    52: |              >                     - ablock(i,5,5)*avec(i,5)
    53: V------          end do
    54:           
    55:           
    56:                 return
    57:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: MATMUL_SUB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    81: vec( 101): Vectorized loop.
    82: vec( 128): Fused multiply-add operation applied.
    87: vec( 128): Fused multiply-add operation applied.
    92: vec( 128): Fused multiply-add operation applied.
    97: vec( 128): Fused multiply-add operation applied.
   102: vec( 128): Fused multiply-add operation applied.
   107: vec( 128): Fused multiply-add operation applied.
   112: vec( 128): Fused multiply-add operation applied.
   117: vec( 128): Fused multiply-add operation applied.
   122: vec( 128): Fused multiply-add operation applied.
   127: vec( 128): Fused multiply-add operation applied.
   132: vec( 128): Fused multiply-add operation applied.
   137: vec( 128): Fused multiply-add operation applied.
   142: vec( 128): Fused multiply-add operation applied.
   147: vec( 128): Fused multiply-add operation applied.
   152: vec( 128): Fused multiply-add operation applied.
   157: vec( 128): Fused multiply-add operation applied.
   162: vec( 128): Fused multiply-add operation applied.
   167: vec( 128): Fused multiply-add operation applied.
   172: vec( 128): Fused multiply-add operation applied.
   177: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   187: vec( 128): Fused multiply-add operation applied.
   192: vec( 128): Fused multiply-add operation applied.
   197: vec( 128): Fused multiply-add operation applied.
   202: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: MATMUL_SUB
FORMAT LIST

 LINE   LOOP      STATEMENT

    62:                 subroutine matmul_sub(ablock, bblock, cblock)
    63:           
    64:           c---------------------------------------------------------------------
    65:           c---------------------------------------------------------------------
    66:           
    67:           c---------------------------------------------------------------------
    68:           c     subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
    69:           c---------------------------------------------------------------------
    70:           
    71:                 implicit none
    72:                 include 'blk_par.h'
    73:           
    74:                 double precision ablock, bblock, cblock
    75:                 dimension ablock(blkdim,5,5), bblock(blkdim,5,5),
    76:                >          cblock(blkdim,5,5)
    77:           
    78:                 integer i
    79:           
    80:           !dir$ vector always
    81: V------>        do i = 1, bsize
    82: |       F          cblock(i,1,1) = cblock(i,1,1) - ablock(i,1,1)*bblock(i,1,1)
    83: |              >                             - ablock(i,1,2)*bblock(i,2,1)
    84: |              >                             - ablock(i,1,3)*bblock(i,3,1)
    85: |              >                             - ablock(i,1,4)*bblock(i,4,1)
    86: |              >                             - ablock(i,1,5)*bblock(i,5,1)
    87: |       F          cblock(i,2,1) = cblock(i,2,1) - ablock(i,2,1)*bblock(i,1,1)
    88: |              >                             - ablock(i,2,2)*bblock(i,2,1)
    89: |              >                             - ablock(i,2,3)*bblock(i,3,1)
    90: |              >                             - ablock(i,2,4)*bblock(i,4,1)
    91: |              >                             - ablock(i,2,5)*bblock(i,5,1)
    92: |       F          cblock(i,3,1) = cblock(i,3,1) - ablock(i,3,1)*bblock(i,1,1)
    93: |              >                             - ablock(i,3,2)*bblock(i,2,1)
    94: |              >                             - ablock(i,3,3)*bblock(i,3,1)
    95: |              >                             - ablock(i,3,4)*bblock(i,4,1)
    96: |              >                             - ablock(i,3,5)*bblock(i,5,1)
    97: |       F          cblock(i,4,1) = cblock(i,4,1) - ablock(i,4,1)*bblock(i,1,1)
    98: |              >                             - ablock(i,4,2)*bblock(i,2,1)
    99: |              >                             - ablock(i,4,3)*bblock(i,3,1)
   100: |              >                             - ablock(i,4,4)*bblock(i,4,1)
   101: |              >                             - ablock(i,4,5)*bblock(i,5,1)
   102: |       F          cblock(i,5,1) = cblock(i,5,1) - ablock(i,5,1)*bblock(i,1,1)
   103: |              >                             - ablock(i,5,2)*bblock(i,2,1)
   104: |              >                             - ablock(i,5,3)*bblock(i,3,1)
   105: |              >                             - ablock(i,5,4)*bblock(i,4,1)
   106: |              >                             - ablock(i,5,5)*bblock(i,5,1)
   107: |       F          cblock(i,1,2) = cblock(i,1,2) - ablock(i,1,1)*bblock(i,1,2)
   108: |              >                             - ablock(i,1,2)*bblock(i,2,2)
   109: |              >                             - ablock(i,1,3)*bblock(i,3,2)
   110: |              >                             - ablock(i,1,4)*bblock(i,4,2)
   111: |              >                             - ablock(i,1,5)*bblock(i,5,2)
   112: |       F          cblock(i,2,2) = cblock(i,2,2) - ablock(i,2,1)*bblock(i,1,2)
   113: |              >                             - ablock(i,2,2)*bblock(i,2,2)
   114: |              >                             - ablock(i,2,3)*bblock(i,3,2)
   115: |              >                             - ablock(i,2,4)*bblock(i,4,2)
   116: |              >                             - ablock(i,2,5)*bblock(i,5,2)
   117: |       F          cblock(i,3,2) = cblock(i,3,2) - ablock(i,3,1)*bblock(i,1,2)
   118: |              >                             - ablock(i,3,2)*bblock(i,2,2)
   119: |              >                             - ablock(i,3,3)*bblock(i,3,2)
   120: |              >                             - ablock(i,3,4)*bblock(i,4,2)
   121: |              >                             - ablock(i,3,5)*bblock(i,5,2)
   122: |       F          cblock(i,4,2) = cblock(i,4,2) - ablock(i,4,1)*bblock(i,1,2)
   123: |              >                             - ablock(i,4,2)*bblock(i,2,2)
   124: |              >                             - ablock(i,4,3)*bblock(i,3,2)
   125: |              >                             - ablock(i,4,4)*bblock(i,4,2)
   126: |              >                             - ablock(i,4,5)*bblock(i,5,2)
   127: |       F          cblock(i,5,2) = cblock(i,5,2) - ablock(i,5,1)*bblock(i,1,2)
   128: |              >                             - ablock(i,5,2)*bblock(i,2,2)
   129: |              >                             - ablock(i,5,3)*bblock(i,3,2)
   130: |              >                             - ablock(i,5,4)*bblock(i,4,2)
   131: |              >                             - ablock(i,5,5)*bblock(i,5,2)
   132: |       F          cblock(i,1,3) = cblock(i,1,3) - ablock(i,1,1)*bblock(i,1,3)
   133: |              >                             - ablock(i,1,2)*bblock(i,2,3)
   134: |              >                             - ablock(i,1,3)*bblock(i,3,3)
   135: |              >                             - ablock(i,1,4)*bblock(i,4,3)
   136: |              >                             - ablock(i,1,5)*bblock(i,5,3)
   137: |       F          cblock(i,2,3) = cblock(i,2,3) - ablock(i,2,1)*bblock(i,1,3)
   138: |              >                             - ablock(i,2,2)*bblock(i,2,3)
   139: |              >                             - ablock(i,2,3)*bblock(i,3,3)
   140: |              >                             - ablock(i,2,4)*bblock(i,4,3)
   141: |              >                             - ablock(i,2,5)*bblock(i,5,3)
   142: |       F          cblock(i,3,3) = cblock(i,3,3) - ablock(i,3,1)*bblock(i,1,3)
   143: |              >                             - ablock(i,3,2)*bblock(i,2,3)
   144: |              >                             - ablock(i,3,3)*bblock(i,3,3)
   145: |              >                             - ablock(i,3,4)*bblock(i,4,3)
   146: |              >                             - ablock(i,3,5)*bblock(i,5,3)
   147: |       F          cblock(i,4,3) = cblock(i,4,3) - ablock(i,4,1)*bblock(i,1,3)
   148: |              >                             - ablock(i,4,2)*bblock(i,2,3)
   149: |              >                             - ablock(i,4,3)*bblock(i,3,3)
   150: |              >                             - ablock(i,4,4)*bblock(i,4,3)
   151: |              >                             - ablock(i,4,5)*bblock(i,5,3)
   152: |       F          cblock(i,5,3) = cblock(i,5,3) - ablock(i,5,1)*bblock(i,1,3)
   153: |              >                             - ablock(i,5,2)*bblock(i,2,3)
   154: |              >                             - ablock(i,5,3)*bblock(i,3,3)
   155: |              >                             - ablock(i,5,4)*bblock(i,4,3)
   156: |              >                             - ablock(i,5,5)*bblock(i,5,3)
   157: |       F          cblock(i,1,4) = cblock(i,1,4) - ablock(i,1,1)*bblock(i,1,4)
   158: |              >                             - ablock(i,1,2)*bblock(i,2,4)
   159: |              >                             - ablock(i,1,3)*bblock(i,3,4)
   160: |              >                             - ablock(i,1,4)*bblock(i,4,4)
   161: |              >                             - ablock(i,1,5)*bblock(i,5,4)
   162: |       F          cblock(i,2,4) = cblock(i,2,4) - ablock(i,2,1)*bblock(i,1,4)
   163: |              >                             - ablock(i,2,2)*bblock(i,2,4)
   164: |              >                             - ablock(i,2,3)*bblock(i,3,4)
   165: |              >                             - ablock(i,2,4)*bblock(i,4,4)
   166: |              >                             - ablock(i,2,5)*bblock(i,5,4)
   167: |       F          cblock(i,3,4) = cblock(i,3,4) - ablock(i,3,1)*bblock(i,1,4)
   168: |              >                             - ablock(i,3,2)*bblock(i,2,4)
   169: |              >                             - ablock(i,3,3)*bblock(i,3,4)
   170: |              >                             - ablock(i,3,4)*bblock(i,4,4)
   171: |              >                             - ablock(i,3,5)*bblock(i,5,4)
   172: |       F          cblock(i,4,4) = cblock(i,4,4) - ablock(i,4,1)*bblock(i,1,4)
   173: |              >                             - ablock(i,4,2)*bblock(i,2,4)
   174: |              >                             - ablock(i,4,3)*bblock(i,3,4)
   175: |              >                             - ablock(i,4,4)*bblock(i,4,4)
   176: |              >                             - ablock(i,4,5)*bblock(i,5,4)
   177: |       F          cblock(i,5,4) = cblock(i,5,4) - ablock(i,5,1)*bblock(i,1,4)
   178: |              >                             - ablock(i,5,2)*bblock(i,2,4)
   179: |              >                             - ablock(i,5,3)*bblock(i,3,4)
   180: |              >                             - ablock(i,5,4)*bblock(i,4,4)
   181: |              >                             - ablock(i,5,5)*bblock(i,5,4)
   182: |       F          cblock(i,1,5) = cblock(i,1,5) - ablock(i,1,1)*bblock(i,1,5)
   183: |              >                             - ablock(i,1,2)*bblock(i,2,5)
   184: |              >                             - ablock(i,1,3)*bblock(i,3,5)
   185: |              >                             - ablock(i,1,4)*bblock(i,4,5)
   186: |              >                             - ablock(i,1,5)*bblock(i,5,5)
   187: |       F          cblock(i,2,5) = cblock(i,2,5) - ablock(i,2,1)*bblock(i,1,5)
   188: |              >                             - ablock(i,2,2)*bblock(i,2,5)
   189: |              >                             - ablock(i,2,3)*bblock(i,3,5)
   190: |              >                             - ablock(i,2,4)*bblock(i,4,5)
   191: |              >                             - ablock(i,2,5)*bblock(i,5,5)
   192: |       F          cblock(i,3,5) = cblock(i,3,5) - ablock(i,3,1)*bblock(i,1,5)
   193: |              >                             - ablock(i,3,2)*bblock(i,2,5)
   194: |              >                             - ablock(i,3,3)*bblock(i,3,5)
   195: |              >                             - ablock(i,3,4)*bblock(i,4,5)
   196: |              >                             - ablock(i,3,5)*bblock(i,5,5)
   197: |       F          cblock(i,4,5) = cblock(i,4,5) - ablock(i,4,1)*bblock(i,1,5)
   198: |              >                             - ablock(i,4,2)*bblock(i,2,5)
   199: |              >                             - ablock(i,4,3)*bblock(i,3,5)
   200: |              >                             - ablock(i,4,4)*bblock(i,4,5)
   201: |              >                             - ablock(i,4,5)*bblock(i,5,5)
   202: |       F          cblock(i,5,5) = cblock(i,5,5) - ablock(i,5,1)*bblock(i,1,5)
   203: |              >                             - ablock(i,5,2)*bblock(i,2,5)
   204: |              >                             - ablock(i,5,3)*bblock(i,3,5)
   205: |              >                             - ablock(i,5,4)*bblock(i,4,5)
   206: |              >                             - ablock(i,5,5)*bblock(i,5,5)
   207: V------         end do
   208:            
   209:                 return
   210:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: BINVCRHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   240: vec( 101): Vectorized loop.
   254: vec( 128): Fused multiply-add operation applied.
   255: vec( 128): Fused multiply-add operation applied.
   256: vec( 128): Fused multiply-add operation applied.
   257: vec( 128): Fused multiply-add operation applied.
   258: vec( 128): Fused multiply-add operation applied.
   259: vec( 128): Fused multiply-add operation applied.
   260: vec( 128): Fused multiply-add operation applied.
   261: vec( 128): Fused multiply-add operation applied.
   262: vec( 128): Fused multiply-add operation applied.
   263: vec( 128): Fused multiply-add operation applied.
   266: vec( 128): Fused multiply-add operation applied.
   267: vec( 128): Fused multiply-add operation applied.
   268: vec( 128): Fused multiply-add operation applied.
   269: vec( 128): Fused multiply-add operation applied.
   270: vec( 128): Fused multiply-add operation applied.
   271: vec( 128): Fused multiply-add operation applied.
   272: vec( 128): Fused multiply-add operation applied.
   273: vec( 128): Fused multiply-add operation applied.
   274: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   278: vec( 128): Fused multiply-add operation applied.
   279: vec( 128): Fused multiply-add operation applied.
   280: vec( 128): Fused multiply-add operation applied.
   281: vec( 128): Fused multiply-add operation applied.
   282: vec( 128): Fused multiply-add operation applied.
   283: vec( 128): Fused multiply-add operation applied.
   284: vec( 128): Fused multiply-add operation applied.
   285: vec( 128): Fused multiply-add operation applied.
   286: vec( 128): Fused multiply-add operation applied.
   287: vec( 128): Fused multiply-add operation applied.
   290: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   292: vec( 128): Fused multiply-add operation applied.
   293: vec( 128): Fused multiply-add operation applied.
   294: vec( 128): Fused multiply-add operation applied.
   295: vec( 128): Fused multiply-add operation applied.
   296: vec( 128): Fused multiply-add operation applied.
   297: vec( 128): Fused multiply-add operation applied.
   298: vec( 128): Fused multiply-add operation applied.
   299: vec( 128): Fused multiply-add operation applied.
   325: vec( 128): Fused multiply-add operation applied.
   326: vec( 128): Fused multiply-add operation applied.
   327: vec( 128): Fused multiply-add operation applied.
   328: vec( 128): Fused multiply-add operation applied.
   329: vec( 128): Fused multiply-add operation applied.
   330: vec( 128): Fused multiply-add operation applied.
   331: vec( 128): Fused multiply-add operation applied.
   332: vec( 128): Fused multiply-add operation applied.
   333: vec( 128): Fused multiply-add operation applied.
   336: vec( 128): Fused multiply-add operation applied.
   337: vec( 128): Fused multiply-add operation applied.
   338: vec( 128): Fused multiply-add operation applied.
   339: vec( 128): Fused multiply-add operation applied.
   340: vec( 128): Fused multiply-add operation applied.
   341: vec( 128): Fused multiply-add operation applied.
   342: vec( 128): Fused multiply-add operation applied.
   343: vec( 128): Fused multiply-add operation applied.
   344: vec( 128): Fused multiply-add operation applied.
   347: vec( 128): Fused multiply-add operation applied.
   348: vec( 128): Fused multiply-add operation applied.
   349: vec( 128): Fused multiply-add operation applied.
   350: vec( 128): Fused multiply-add operation applied.
   351: vec( 128): Fused multiply-add operation applied.
   352: vec( 128): Fused multiply-add operation applied.
   353: vec( 128): Fused multiply-add operation applied.
   354: vec( 128): Fused multiply-add operation applied.
   355: vec( 128): Fused multiply-add operation applied.
   369: vec( 128): Fused multiply-add operation applied.
   370: vec( 128): Fused multiply-add operation applied.
   371: vec( 128): Fused multiply-add operation applied.
   372: vec( 128): Fused multiply-add operation applied.
   373: vec( 128): Fused multiply-add operation applied.
   374: vec( 128): Fused multiply-add operation applied.
   375: vec( 128): Fused multiply-add operation applied.
   376: vec( 128): Fused multiply-add operation applied.
   389: vec( 128): Fused multiply-add operation applied.
   390: vec( 128): Fused multiply-add operation applied.
   391: vec( 128): Fused multiply-add operation applied.
   392: vec( 128): Fused multiply-add operation applied.
   393: vec( 128): Fused multiply-add operation applied.
   394: vec( 128): Fused multiply-add operation applied.
   395: vec( 128): Fused multiply-add operation applied.
   396: vec( 128): Fused multiply-add operation applied.
   399: vec( 128): Fused multiply-add operation applied.
   400: vec( 128): Fused multiply-add operation applied.
   401: vec( 128): Fused multiply-add operation applied.
   402: vec( 128): Fused multiply-add operation applied.
   403: vec( 128): Fused multiply-add operation applied.
   404: vec( 128): Fused multiply-add operation applied.
   405: vec( 128): Fused multiply-add operation applied.
   406: vec( 128): Fused multiply-add operation applied.
   419: vec( 128): Fused multiply-add operation applied.
   420: vec( 128): Fused multiply-add operation applied.
   421: vec( 128): Fused multiply-add operation applied.
   422: vec( 128): Fused multiply-add operation applied.
   423: vec( 128): Fused multiply-add operation applied.
   424: vec( 128): Fused multiply-add operation applied.
   425: vec( 128): Fused multiply-add operation applied.
   428: vec( 128): Fused multiply-add operation applied.
   429: vec( 128): Fused multiply-add operation applied.
   430: vec( 128): Fused multiply-add operation applied.
   431: vec( 128): Fused multiply-add operation applied.
   432: vec( 128): Fused multiply-add operation applied.
   433: vec( 128): Fused multiply-add operation applied.
   434: vec( 128): Fused multiply-add operation applied.
   446: vec( 128): Fused multiply-add operation applied.
   447: vec( 128): Fused multiply-add operation applied.
   448: vec( 128): Fused multiply-add operation applied.
   449: vec( 128): Fused multiply-add operation applied.
   450: vec( 128): Fused multiply-add operation applied.
   451: vec( 128): Fused multiply-add operation applied.
   452: vec( 128): Fused multiply-add operation applied.
   464: vec( 128): Fused multiply-add operation applied.
   465: vec( 128): Fused multiply-add operation applied.
   466: vec( 128): Fused multiply-add operation applied.
   467: vec( 128): Fused multiply-add operation applied.
   468: vec( 128): Fused multiply-add operation applied.
   469: vec( 128): Fused multiply-add operation applied.
   472: vec( 128): Fused multiply-add operation applied.
   473: vec( 128): Fused multiply-add operation applied.
   474: vec( 128): Fused multiply-add operation applied.
   475: vec( 128): Fused multiply-add operation applied.
   476: vec( 128): Fused multiply-add operation applied.
   477: vec( 128): Fused multiply-add operation applied.
   480: vec( 128): Fused multiply-add operation applied.
   481: vec( 128): Fused multiply-add operation applied.
   482: vec( 128): Fused multiply-add operation applied.
   483: vec( 128): Fused multiply-add operation applied.
   484: vec( 128): Fused multiply-add operation applied.
   485: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: BINVCRHS
FORMAT LIST

 LINE   LOOP      STATEMENT

   217:                 subroutine binvcrhs( lhs,c,r )
   218:           
   219:           c---------------------------------------------------------------------
   220:           c---------------------------------------------------------------------
   221:           
   222:           c---------------------------------------------------------------------
   223:           c
   224:           c---------------------------------------------------------------------
   225:           
   226:                 implicit none
   227:                 include 'blk_par.h'
   228:           
   229:                 double precision pivot, coeff, lhs
   230:                 dimension lhs(blkdim,5,5)
   231:                 double precision c(blkdim,5,5), r(blkdim,5)
   232:           
   233:                 integer i
   234:           
   235:           c---------------------------------------------------------------------
   236:           c
   237:           c---------------------------------------------------------------------
   238:           
   239:           !dir$ vector always
   240: V------>        do i = 1, bsize
   241: |               pivot = 1.00d0/lhs(i,1,1)
   242: |               lhs(i,1,2) = lhs(i,1,2)*pivot
   243: |               lhs(i,1,3) = lhs(i,1,3)*pivot
   244: |               lhs(i,1,4) = lhs(i,1,4)*pivot
   245: |               lhs(i,1,5) = lhs(i,1,5)*pivot
   246: |               c(i,1,1) = c(i,1,1)*pivot
   247: |               c(i,1,2) = c(i,1,2)*pivot
   248: |               c(i,1,3) = c(i,1,3)*pivot
   249: |               c(i,1,4) = c(i,1,4)*pivot
   250: |               c(i,1,5) = c(i,1,5)*pivot
   251: |               r(i,1)   = r(i,1)  *pivot
   252: |         
   253: |               coeff = lhs(i,2,1)
   254: |       F       lhs(i,2,2)= lhs(i,2,2) - coeff*lhs(i,1,2)
   255: |       F       lhs(i,2,3)= lhs(i,2,3) - coeff*lhs(i,1,3)
   256: |       F       lhs(i,2,4)= lhs(i,2,4) - coeff*lhs(i,1,4)
   257: |       F       lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,1,5)
   258: |       F       c(i,2,1) = c(i,2,1) - coeff*c(i,1,1)
   259: |       F       c(i,2,2) = c(i,2,2) - coeff*c(i,1,2)
   260: |       F       c(i,2,3) = c(i,2,3) - coeff*c(i,1,3)
   261: |       F       c(i,2,4) = c(i,2,4) - coeff*c(i,1,4)
   262: |       F       c(i,2,5) = c(i,2,5) - coeff*c(i,1,5)
   263: |       F       r(i,2)   = r(i,2)   - coeff*r(i,1)
   264: |         
   265: |               coeff = lhs(i,3,1)
   266: |       F       lhs(i,3,2)= lhs(i,3,2) - coeff*lhs(i,1,2)
   267: |       F       lhs(i,3,3)= lhs(i,3,3) - coeff*lhs(i,1,3)
   268: |       F       lhs(i,3,4)= lhs(i,3,4) - coeff*lhs(i,1,4)
   269: |       F       lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,1,5)
   270: |       F       c(i,3,1) = c(i,3,1) - coeff*c(i,1,1)
   271: |       F       c(i,3,2) = c(i,3,2) - coeff*c(i,1,2)
   272: |       F       c(i,3,3) = c(i,3,3) - coeff*c(i,1,3)
   273: |       F       c(i,3,4) = c(i,3,4) - coeff*c(i,1,4)
   274: |       F       c(i,3,5) = c(i,3,5) - coeff*c(i,1,5)
   275: |       F       r(i,3)   = r(i,3)   - coeff*r(i,1)
   276: |         
   277: |               coeff = lhs(i,4,1)
   278: |       F       lhs(i,4,2)= lhs(i,4,2) - coeff*lhs(i,1,2)
   279: |       F       lhs(i,4,3)= lhs(i,4,3) - coeff*lhs(i,1,3)
   280: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,1,4)
   281: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,1,5)
   282: |       F       c(i,4,1) = c(i,4,1) - coeff*c(i,1,1)
   283: |       F       c(i,4,2) = c(i,4,2) - coeff*c(i,1,2)
   284: |       F       c(i,4,3) = c(i,4,3) - coeff*c(i,1,3)
   285: |       F       c(i,4,4) = c(i,4,4) - coeff*c(i,1,4)
   286: |       F       c(i,4,5) = c(i,4,5) - coeff*c(i,1,5)
   287: |       F       r(i,4)   = r(i,4)   - coeff*r(i,1)
   288: |         
   289: |               coeff = lhs(i,5,1)
   290: |       F       lhs(i,5,2)= lhs(i,5,2) - coeff*lhs(i,1,2)
   291: |       F       lhs(i,5,3)= lhs(i,5,3) - coeff*lhs(i,1,3)
   292: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,1,4)
   293: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,1,5)
   294: |       F       c(i,5,1) = c(i,5,1) - coeff*c(i,1,1)
   295: |       F       c(i,5,2) = c(i,5,2) - coeff*c(i,1,2)
   296: |       F       c(i,5,3) = c(i,5,3) - coeff*c(i,1,3)
   297: |       F       c(i,5,4) = c(i,5,4) - coeff*c(i,1,4)
   298: |       F       c(i,5,5) = c(i,5,5) - coeff*c(i,1,5)
   299: |       F       r(i,5)   = r(i,5)   - coeff*r(i,1)
   300: |         
   301: |         
   302: |               pivot = 1.00d0/lhs(i,2,2)
   303: |               lhs(i,2,3) = lhs(i,2,3)*pivot
   304: |               lhs(i,2,4) = lhs(i,2,4)*pivot
   305: |               lhs(i,2,5) = lhs(i,2,5)*pivot
   306: |               c(i,2,1) = c(i,2,1)*pivot
   307: |               c(i,2,2) = c(i,2,2)*pivot
   308: |               c(i,2,3) = c(i,2,3)*pivot
   309: |               c(i,2,4) = c(i,2,4)*pivot
   310: |               c(i,2,5) = c(i,2,5)*pivot
   311: |               r(i,2)   = r(i,2)  *pivot
   312: |         
   313: |               coeff = lhs(i,1,2)
   314: |               lhs(i,1,3)= lhs(i,1,3) - coeff*lhs(i,2,3)
   315: |               lhs(i,1,4)= lhs(i,1,4) - coeff*lhs(i,2,4)
   316: |               lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,2,5)
   317: |               c(i,1,1) = c(i,1,1) - coeff*c(i,2,1)
   318: |               c(i,1,2) = c(i,1,2) - coeff*c(i,2,2)
   319: |               c(i,1,3) = c(i,1,3) - coeff*c(i,2,3)
   320: |               c(i,1,4) = c(i,1,4) - coeff*c(i,2,4)
   321: |               c(i,1,5) = c(i,1,5) - coeff*c(i,2,5)
   322: |               r(i,1)   = r(i,1)   - coeff*r(i,2)
   323: |         
   324: |               coeff = lhs(i,3,2)
   325: |       F       lhs(i,3,3)= lhs(i,3,3) - coeff*lhs(i,2,3)
   326: |       F       lhs(i,3,4)= lhs(i,3,4) - coeff*lhs(i,2,4)
   327: |       F       lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,2,5)
   328: |       F       c(i,3,1) = c(i,3,1) - coeff*c(i,2,1)
   329: |       F       c(i,3,2) = c(i,3,2) - coeff*c(i,2,2)
   330: |       F       c(i,3,3) = c(i,3,3) - coeff*c(i,2,3)
   331: |       F       c(i,3,4) = c(i,3,4) - coeff*c(i,2,4)
   332: |       F       c(i,3,5) = c(i,3,5) - coeff*c(i,2,5)
   333: |       F       r(i,3)   = r(i,3)   - coeff*r(i,2)
   334: |         
   335: |               coeff = lhs(i,4,2)
   336: |       F       lhs(i,4,3)= lhs(i,4,3) - coeff*lhs(i,2,3)
   337: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,2,4)
   338: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,2,5)
   339: |       F       c(i,4,1) = c(i,4,1) - coeff*c(i,2,1)
   340: |       F       c(i,4,2) = c(i,4,2) - coeff*c(i,2,2)
   341: |       F       c(i,4,3) = c(i,4,3) - coeff*c(i,2,3)
   342: |       F       c(i,4,4) = c(i,4,4) - coeff*c(i,2,4)
   343: |       F       c(i,4,5) = c(i,4,5) - coeff*c(i,2,5)
   344: |       F       r(i,4)   = r(i,4)   - coeff*r(i,2)
   345: |         
   346: |               coeff = lhs(i,5,2)
   347: |       F       lhs(i,5,3)= lhs(i,5,3) - coeff*lhs(i,2,3)
   348: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,2,4)
   349: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,2,5)
   350: |       F       c(i,5,1) = c(i,5,1) - coeff*c(i,2,1)
   351: |       F       c(i,5,2) = c(i,5,2) - coeff*c(i,2,2)
   352: |       F       c(i,5,3) = c(i,5,3) - coeff*c(i,2,3)
   353: |       F       c(i,5,4) = c(i,5,4) - coeff*c(i,2,4)
   354: |       F       c(i,5,5) = c(i,5,5) - coeff*c(i,2,5)
   355: |       F       r(i,5)   = r(i,5)   - coeff*r(i,2)
   356: |         
   357: |         
   358: |               pivot = 1.00d0/lhs(i,3,3)
   359: |               lhs(i,3,4) = lhs(i,3,4)*pivot
   360: |               lhs(i,3,5) = lhs(i,3,5)*pivot
   361: |               c(i,3,1) = c(i,3,1)*pivot
   362: |               c(i,3,2) = c(i,3,2)*pivot
   363: |               c(i,3,3) = c(i,3,3)*pivot
   364: |               c(i,3,4) = c(i,3,4)*pivot
   365: |               c(i,3,5) = c(i,3,5)*pivot
   366: |               r(i,3)   = r(i,3)  *pivot
   367: |         
   368: |               coeff = lhs(i,1,3)
   369: |       F       lhs(i,1,4)= lhs(i,1,4) - coeff*lhs(i,3,4)
   370: |       F       lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,3,5)
   371: |       F       c(i,1,1) = c(i,1,1) - coeff*c(i,3,1)
   372: |       F       c(i,1,2) = c(i,1,2) - coeff*c(i,3,2)
   373: |       F       c(i,1,3) = c(i,1,3) - coeff*c(i,3,3)
   374: |       F       c(i,1,4) = c(i,1,4) - coeff*c(i,3,4)
   375: |       F       c(i,1,5) = c(i,1,5) - coeff*c(i,3,5)
   376: |       F       r(i,1)   = r(i,1)   - coeff*r(i,3)
   377: |         
   378: |               coeff = lhs(i,2,3)
   379: |               lhs(i,2,4)= lhs(i,2,4) - coeff*lhs(i,3,4)
   380: |               lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,3,5)
   381: |               c(i,2,1) = c(i,2,1) - coeff*c(i,3,1)
   382: |               c(i,2,2) = c(i,2,2) - coeff*c(i,3,2)
   383: |               c(i,2,3) = c(i,2,3) - coeff*c(i,3,3)
   384: |               c(i,2,4) = c(i,2,4) - coeff*c(i,3,4)
   385: |               c(i,2,5) = c(i,2,5) - coeff*c(i,3,5)
   386: |               r(i,2)   = r(i,2)   - coeff*r(i,3)
   387: |         
   388: |               coeff = lhs(i,4,3)
   389: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,3,4)
   390: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,3,5)
   391: |       F       c(i,4,1) = c(i,4,1) - coeff*c(i,3,1)
   392: |       F       c(i,4,2) = c(i,4,2) - coeff*c(i,3,2)
   393: |       F       c(i,4,3) = c(i,4,3) - coeff*c(i,3,3)
   394: |       F       c(i,4,4) = c(i,4,4) - coeff*c(i,3,4)
   395: |       F       c(i,4,5) = c(i,4,5) - coeff*c(i,3,5)
   396: |       F       r(i,4)   = r(i,4)   - coeff*r(i,3)
   397: |         
   398: |               coeff = lhs(i,5,3)
   399: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,3,4)
   400: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,3,5)
   401: |       F       c(i,5,1) = c(i,5,1) - coeff*c(i,3,1)
   402: |       F       c(i,5,2) = c(i,5,2) - coeff*c(i,3,2)
   403: |       F       c(i,5,3) = c(i,5,3) - coeff*c(i,3,3)
   404: |       F       c(i,5,4) = c(i,5,4) - coeff*c(i,3,4)
   405: |       F       c(i,5,5) = c(i,5,5) - coeff*c(i,3,5)
   406: |       F       r(i,5)   = r(i,5)   - coeff*r(i,3)
   407: |         
   408: |         
   409: |               pivot = 1.00d0/lhs(i,4,4)
   410: |               lhs(i,4,5) = lhs(i,4,5)*pivot
   411: |               c(i,4,1) = c(i,4,1)*pivot
   412: |               c(i,4,2) = c(i,4,2)*pivot
   413: |               c(i,4,3) = c(i,4,3)*pivot
   414: |               c(i,4,4) = c(i,4,4)*pivot
   415: |               c(i,4,5) = c(i,4,5)*pivot
   416: |               r(i,4)   = r(i,4)  *pivot
   417: |         
   418: |               coeff = lhs(i,1,4)
   419: |       F       lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,4,5)
   420: |       F       c(i,1,1) = c(i,1,1) - coeff*c(i,4,1)
   421: |       F       c(i,1,2) = c(i,1,2) - coeff*c(i,4,2)
   422: |       F       c(i,1,3) = c(i,1,3) - coeff*c(i,4,3)
   423: |       F       c(i,1,4) = c(i,1,4) - coeff*c(i,4,4)
   424: |       F       c(i,1,5) = c(i,1,5) - coeff*c(i,4,5)
   425: |       F       r(i,1)   = r(i,1)   - coeff*r(i,4)
   426: |         
   427: |               coeff = lhs(i,2,4)
   428: |       F       lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,4,5)
   429: |       F       c(i,2,1) = c(i,2,1) - coeff*c(i,4,1)
   430: |       F       c(i,2,2) = c(i,2,2) - coeff*c(i,4,2)
   431: |       F       c(i,2,3) = c(i,2,3) - coeff*c(i,4,3)
   432: |       F       c(i,2,4) = c(i,2,4) - coeff*c(i,4,4)
   433: |       F       c(i,2,5) = c(i,2,5) - coeff*c(i,4,5)
   434: |       F       r(i,2)   = r(i,2)   - coeff*r(i,4)
   435: |         
   436: |               coeff = lhs(i,3,4)
   437: |               lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,4,5)
   438: |               c(i,3,1) = c(i,3,1) - coeff*c(i,4,1)
   439: |               c(i,3,2) = c(i,3,2) - coeff*c(i,4,2)
   440: |               c(i,3,3) = c(i,3,3) - coeff*c(i,4,3)
   441: |               c(i,3,4) = c(i,3,4) - coeff*c(i,4,4)
   442: |               c(i,3,5) = c(i,3,5) - coeff*c(i,4,5)
   443: |               r(i,3)   = r(i,3)   - coeff*r(i,4)
   444: |         
   445: |               coeff = lhs(i,5,4)
   446: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,4,5)
   447: |       F       c(i,5,1) = c(i,5,1) - coeff*c(i,4,1)
   448: |       F       c(i,5,2) = c(i,5,2) - coeff*c(i,4,2)
   449: |       F       c(i,5,3) = c(i,5,3) - coeff*c(i,4,3)
   450: |       F       c(i,5,4) = c(i,5,4) - coeff*c(i,4,4)
   451: |       F       c(i,5,5) = c(i,5,5) - coeff*c(i,4,5)
   452: |       F       r(i,5)   = r(i,5)   - coeff*r(i,4)
   453: |         
   454: |         
   455: |               pivot = 1.00d0/lhs(i,5,5)
   456: |               c(i,5,1) = c(i,5,1)*pivot
   457: |               c(i,5,2) = c(i,5,2)*pivot
   458: |               c(i,5,3) = c(i,5,3)*pivot
   459: |               c(i,5,4) = c(i,5,4)*pivot
   460: |               c(i,5,5) = c(i,5,5)*pivot
   461: |               r(i,5)   = r(i,5)  *pivot
   462: |         
   463: |               coeff = lhs(i,1,5)
   464: |       F       c(i,1,1) = c(i,1,1) - coeff*c(i,5,1)
   465: |       F       c(i,1,2) = c(i,1,2) - coeff*c(i,5,2)
   466: |       F       c(i,1,3) = c(i,1,3) - coeff*c(i,5,3)
   467: |       F       c(i,1,4) = c(i,1,4) - coeff*c(i,5,4)
   468: |       F       c(i,1,5) = c(i,1,5) - coeff*c(i,5,5)
   469: |       F       r(i,1)   = r(i,1)   - coeff*r(i,5)
   470: |         
   471: |               coeff = lhs(i,2,5)
   472: |       F       c(i,2,1) = c(i,2,1) - coeff*c(i,5,1)
   473: |       F       c(i,2,2) = c(i,2,2) - coeff*c(i,5,2)
   474: |       F       c(i,2,3) = c(i,2,3) - coeff*c(i,5,3)
   475: |       F       c(i,2,4) = c(i,2,4) - coeff*c(i,5,4)
   476: |       F       c(i,2,5) = c(i,2,5) - coeff*c(i,5,5)
   477: |       F       r(i,2)   = r(i,2)   - coeff*r(i,5)
   478: |         
   479: |               coeff = lhs(i,3,5)
   480: |       F       c(i,3,1) = c(i,3,1) - coeff*c(i,5,1)
   481: |       F       c(i,3,2) = c(i,3,2) - coeff*c(i,5,2)
   482: |       F       c(i,3,3) = c(i,3,3) - coeff*c(i,5,3)
   483: |       F       c(i,3,4) = c(i,3,4) - coeff*c(i,5,4)
   484: |       F       c(i,3,5) = c(i,3,5) - coeff*c(i,5,5)
   485: |       F       r(i,3)   = r(i,3)   - coeff*r(i,5)
   486: |         
   487: |               coeff = lhs(i,4,5)
   488: |               c(i,4,1) = c(i,4,1) - coeff*c(i,5,1)
   489: |               c(i,4,2) = c(i,4,2) - coeff*c(i,5,2)
   490: |               c(i,4,3) = c(i,4,3) - coeff*c(i,5,3)
   491: |               c(i,4,4) = c(i,4,4) - coeff*c(i,5,4)
   492: |               c(i,4,5) = c(i,4,5) - coeff*c(i,5,5)
   493: |               r(i,4)   = r(i,4)   - coeff*r(i,5)
   494: V------         end do
   495:           
   496:           
   497:                 return
   498:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: BINVRHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   528: vec( 101): Vectorized loop.
   538: vec( 128): Fused multiply-add operation applied.
   539: vec( 128): Fused multiply-add operation applied.
   540: vec( 128): Fused multiply-add operation applied.
   541: vec( 128): Fused multiply-add operation applied.
   542: vec( 128): Fused multiply-add operation applied.
   545: vec( 128): Fused multiply-add operation applied.
   546: vec( 128): Fused multiply-add operation applied.
   547: vec( 128): Fused multiply-add operation applied.
   548: vec( 128): Fused multiply-add operation applied.
   549: vec( 128): Fused multiply-add operation applied.
   552: vec( 128): Fused multiply-add operation applied.
   553: vec( 128): Fused multiply-add operation applied.
   554: vec( 128): Fused multiply-add operation applied.
   555: vec( 128): Fused multiply-add operation applied.
   556: vec( 128): Fused multiply-add operation applied.
   559: vec( 128): Fused multiply-add operation applied.
   560: vec( 128): Fused multiply-add operation applied.
   561: vec( 128): Fused multiply-add operation applied.
   562: vec( 128): Fused multiply-add operation applied.
   563: vec( 128): Fused multiply-add operation applied.
   579: vec( 128): Fused multiply-add operation applied.
   580: vec( 128): Fused multiply-add operation applied.
   581: vec( 128): Fused multiply-add operation applied.
   582: vec( 128): Fused multiply-add operation applied.
   585: vec( 128): Fused multiply-add operation applied.
   586: vec( 128): Fused multiply-add operation applied.
   587: vec( 128): Fused multiply-add operation applied.
   588: vec( 128): Fused multiply-add operation applied.
   591: vec( 128): Fused multiply-add operation applied.
   592: vec( 128): Fused multiply-add operation applied.
   593: vec( 128): Fused multiply-add operation applied.
   594: vec( 128): Fused multiply-add operation applied.
   603: vec( 128): Fused multiply-add operation applied.
   604: vec( 128): Fused multiply-add operation applied.
   605: vec( 128): Fused multiply-add operation applied.
   613: vec( 128): Fused multiply-add operation applied.
   614: vec( 128): Fused multiply-add operation applied.
   615: vec( 128): Fused multiply-add operation applied.
   618: vec( 128): Fused multiply-add operation applied.
   619: vec( 128): Fused multiply-add operation applied.
   620: vec( 128): Fused multiply-add operation applied.
   628: vec( 128): Fused multiply-add operation applied.
   629: vec( 128): Fused multiply-add operation applied.
   632: vec( 128): Fused multiply-add operation applied.
   633: vec( 128): Fused multiply-add operation applied.
   640: vec( 128): Fused multiply-add operation applied.
   641: vec( 128): Fused multiply-add operation applied.
   648: vec( 128): Fused multiply-add operation applied.
   651: vec( 128): Fused multiply-add operation applied.
   654: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:57 2020
FILE NAME: solve_subs_blk.f

PROCEDURE NAME: BINVRHS
FORMAT LIST

 LINE   LOOP      STATEMENT

   505:                 subroutine binvrhs( lhs,r )
   506:           
   507:           c---------------------------------------------------------------------
   508:           c---------------------------------------------------------------------
   509:           
   510:           c---------------------------------------------------------------------
   511:           c
   512:           c---------------------------------------------------------------------
   513:           
   514:                 implicit none
   515:                 include 'blk_par.h'
   516:           
   517:                 double precision pivot, coeff, lhs
   518:                 dimension lhs(blkdim,5,5)
   519:                 double precision r(blkdim,5)
   520:           
   521:                 integer i
   522:           
   523:           c---------------------------------------------------------------------
   524:           c
   525:           c---------------------------------------------------------------------
   526:           
   527:           !dir$ vector always
   528: V------>        do i = 1, bsize
   529: |         
   530: |               pivot = 1.00d0/lhs(i,1,1)
   531: |               lhs(i,1,2) = lhs(i,1,2)*pivot
   532: |               lhs(i,1,3) = lhs(i,1,3)*pivot
   533: |               lhs(i,1,4) = lhs(i,1,4)*pivot
   534: |               lhs(i,1,5) = lhs(i,1,5)*pivot
   535: |               r(i,1)   = r(i,1)  *pivot
   536: |         
   537: |               coeff = lhs(i,2,1)
   538: |       F       lhs(i,2,2)= lhs(i,2,2) - coeff*lhs(i,1,2)
   539: |       F       lhs(i,2,3)= lhs(i,2,3) - coeff*lhs(i,1,3)
   540: |       F       lhs(i,2,4)= lhs(i,2,4) - coeff*lhs(i,1,4)
   541: |       F       lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,1,5)
   542: |       F       r(i,2)   = r(i,2)   - coeff*r(i,1)
   543: |         
   544: |               coeff = lhs(i,3,1)
   545: |       F       lhs(i,3,2)= lhs(i,3,2) - coeff*lhs(i,1,2)
   546: |       F       lhs(i,3,3)= lhs(i,3,3) - coeff*lhs(i,1,3)
   547: |       F       lhs(i,3,4)= lhs(i,3,4) - coeff*lhs(i,1,4)
   548: |       F       lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,1,5)
   549: |       F       r(i,3)   = r(i,3)   - coeff*r(i,1)
   550: |         
   551: |               coeff = lhs(i,4,1)
   552: |       F       lhs(i,4,2)= lhs(i,4,2) - coeff*lhs(i,1,2)
   553: |       F       lhs(i,4,3)= lhs(i,4,3) - coeff*lhs(i,1,3)
   554: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,1,4)
   555: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,1,5)
   556: |       F       r(i,4)   = r(i,4)   - coeff*r(i,1)
   557: |         
   558: |               coeff = lhs(i,5,1)
   559: |       F       lhs(i,5,2)= lhs(i,5,2) - coeff*lhs(i,1,2)
   560: |       F       lhs(i,5,3)= lhs(i,5,3) - coeff*lhs(i,1,3)
   561: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,1,4)
   562: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,1,5)
   563: |       F       r(i,5)   = r(i,5)   - coeff*r(i,1)
   564: |         
   565: |         
   566: |               pivot = 1.00d0/lhs(i,2,2)
   567: |               lhs(i,2,3) = lhs(i,2,3)*pivot
   568: |               lhs(i,2,4) = lhs(i,2,4)*pivot
   569: |               lhs(i,2,5) = lhs(i,2,5)*pivot
   570: |               r(i,2)   = r(i,2)  *pivot
   571: |         
   572: |               coeff = lhs(i,1,2)
   573: |               lhs(i,1,3)= lhs(i,1,3) - coeff*lhs(i,2,3)
   574: |               lhs(i,1,4)= lhs(i,1,4) - coeff*lhs(i,2,4)
   575: |               lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,2,5)
   576: |               r(i,1)   = r(i,1)   - coeff*r(i,2)
   577: |         
   578: |               coeff = lhs(i,3,2)
   579: |       F       lhs(i,3,3)= lhs(i,3,3) - coeff*lhs(i,2,3)
   580: |       F       lhs(i,3,4)= lhs(i,3,4) - coeff*lhs(i,2,4)
   581: |       F       lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,2,5)
   582: |       F       r(i,3)   = r(i,3)   - coeff*r(i,2)
   583: |         
   584: |               coeff = lhs(i,4,2)
   585: |       F       lhs(i,4,3)= lhs(i,4,3) - coeff*lhs(i,2,3)
   586: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,2,4)
   587: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,2,5)
   588: |       F       r(i,4)   = r(i,4)   - coeff*r(i,2)
   589: |         
   590: |               coeff = lhs(i,5,2)
   591: |       F       lhs(i,5,3)= lhs(i,5,3) - coeff*lhs(i,2,3)
   592: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,2,4)
   593: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,2,5)
   594: |       F       r(i,5)   = r(i,5)   - coeff*r(i,2)
   595: |         
   596: |         
   597: |               pivot = 1.00d0/lhs(i,3,3)
   598: |               lhs(i,3,4) = lhs(i,3,4)*pivot
   599: |               lhs(i,3,5) = lhs(i,3,5)*pivot
   600: |               r(i,3)   = r(i,3)  *pivot
   601: |         
   602: |               coeff = lhs(i,1,3)
   603: |       F       lhs(i,1,4)= lhs(i,1,4) - coeff*lhs(i,3,4)
   604: |       F       lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,3,5)
   605: |       F       r(i,1)   = r(i,1)   - coeff*r(i,3)
   606: |         
   607: |               coeff = lhs(i,2,3)
   608: |               lhs(i,2,4)= lhs(i,2,4) - coeff*lhs(i,3,4)
   609: |               lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,3,5)
   610: |               r(i,2)   = r(i,2)   - coeff*r(i,3)
   611: |         
   612: |               coeff = lhs(i,4,3)
   613: |       F       lhs(i,4,4)= lhs(i,4,4) - coeff*lhs(i,3,4)
   614: |       F       lhs(i,4,5)= lhs(i,4,5) - coeff*lhs(i,3,5)
   615: |       F       r(i,4)   = r(i,4)   - coeff*r(i,3)
   616: |         
   617: |               coeff = lhs(i,5,3)
   618: |       F       lhs(i,5,4)= lhs(i,5,4) - coeff*lhs(i,3,4)
   619: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,3,5)
   620: |       F       r(i,5)   = r(i,5)   - coeff*r(i,3)
   621: |         
   622: |         
   623: |               pivot = 1.00d0/lhs(i,4,4)
   624: |               lhs(i,4,5) = lhs(i,4,5)*pivot
   625: |               r(i,4)   = r(i,4)  *pivot
   626: |         
   627: |               coeff = lhs(i,1,4)
   628: |       F       lhs(i,1,5)= lhs(i,1,5) - coeff*lhs(i,4,5)
   629: |       F       r(i,1)   = r(i,1)   - coeff*r(i,4)
   630: |         
   631: |               coeff = lhs(i,2,4)
   632: |       F       lhs(i,2,5)= lhs(i,2,5) - coeff*lhs(i,4,5)
   633: |       F       r(i,2)   = r(i,2)   - coeff*r(i,4)
   634: |         
   635: |               coeff = lhs(i,3,4)
   636: |               lhs(i,3,5)= lhs(i,3,5) - coeff*lhs(i,4,5)
   637: |               r(i,3)   = r(i,3)   - coeff*r(i,4)
   638: |         
   639: |               coeff = lhs(i,5,4)
   640: |       F       lhs(i,5,5)= lhs(i,5,5) - coeff*lhs(i,4,5)
   641: |       F       r(i,5)   = r(i,5)   - coeff*r(i,4)
   642: |         
   643: |         
   644: |               pivot = 1.00d0/lhs(i,5,5)
   645: |               r(i,5)   = r(i,5)  *pivot
   646: |         
   647: |               coeff = lhs(i,1,5)
   648: |       F       r(i,1)   = r(i,1)   - coeff*r(i,5)
   649: |         
   650: |               coeff = lhs(i,2,5)
   651: |       F       r(i,2)   = r(i,2)   - coeff*r(i,5)
   652: |         
   653: |               coeff = lhs(i,3,5)
   654: |       F       r(i,3)   = r(i,3)   - coeff*r(i,5)
   655: |         
   656: |               coeff = lhs(i,4,5)
   657: |               r(i,4)   = r(i,4)   - coeff*r(i,5)
   658: V------         end do
   659:           
   660:           
   661:                 return
   662:                 end


