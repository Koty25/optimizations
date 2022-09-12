NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MORTAR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: par(1801): Parallel routine generated.: MORTAR$1
    40: par(1803): Parallelized by "do".
    40: vec( 101): Vectorized loop.
   715: par(1801): Parallel routine generated.: MORTAR$2
   716: par(1803): Parallelized by "do".
   716: vec( 101): Vectorized loop.
   728: par(1801): Parallel routine generated.: MORTAR$3
   729: par(1803): Parallelized by "do".
   734: vec( 110): Vectorization obstructive procedure reference.: MORTAR_VERTEX
   734: vec( 103): Unvectorized loop.
   739: opt(1025): Reference to this procedure inhibits optimization.: MORTAR_VERTEX
   777: par(1801): Parallel routine generated.: MORTAR$4
   780: par(1803): Parallelized by "do".
   780: vec( 102): Partially vectorized loop.
   780: vec( 113): Overhead of loop division is too large.
   844: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
   867: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
   906: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
   945: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
   998: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1033: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1053: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1102: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1137: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1157: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1177: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: NCON_EDGE
  1336: par(1801): Parallel routine generated.: MORTAR$5
  1337: par(1803): Parallelized by "do".
  1337: vec( 101): Vectorized loop.
  1351: par(1801): Parallel routine generated.: MORTAR$6
  1353: par(1803): Parallelized by "do".
  1379: vec( 110): Vectorization obstructive procedure reference.: MOR_NE
  1379: vec( 110): Vectorization obstructive procedure reference.: MOR_EDGE
  1379: vec( 110): Vectorization obstructive procedure reference.: MOR_ASSIGN
  1379: vec( 103): Unvectorized loop.
  1383: opt(1025): Reference to this procedure inhibits optimization.: MOR_ASSIGN
  1442: opt(1418): Constant-length loop is expanded.
  1443: opt(1418): Constant-length loop is expanded.
  1455: opt(1418): Constant-length loop is expanded.
  1456: opt(1418): Constant-length loop is expanded.
  1472: par(1801): Parallel routine generated.: MORTAR$7
  1474: par(1803): Parallelized by "do".
  1475: vec( 110): Vectorization obstructive procedure reference.: EDGECOPY_S
  1487: vec( 103): Unvectorized loop.
  1487: vec( 107): Iteration count is too small.
  1489: opt(1418): Constant-length loop is expanded.
  1490: opt(1418): Constant-length loop is expanded.
  1491: opt(1036): Potential feedback - use directive or compiler option if OK.: IDMO


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MORTAR
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine mortar
     3:           c-----------------------------------------------------------------
     4:           c     generate mortar point index number
     5:           c-----------------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:           
    10:                 integer count, iel, jface, ntemp, i, ii, jj, ntemp1,
    11:                &        iii, jjj, face2, ne, ie, edge_g, ie2,
    12:                &        mor_v(3), cb, cb1, cb2, cb3, cb4, cb5, cb6,
    13:                &        space, sumcb, ij1, ij2, n1, n2, n3, n4, n5
    14:           
    15:                 n1=lx1*lx1*6*4*nelt
    16:                 n2=8*nelt
    17:                 n3=2*64*nelt
    18:                 n4=12*nelt
    19:                 n5=2*12*nelt
    20:           
    21:                 call nr_init_omp(idmo,n1,0)
    22:                 call nr_init_omp(nemo,n2,0)
    23:                 call nr_init_omp(vassign,n2,0)
    24:                 call nr_init_omp(emo,n3,0)
    25:                 call  l_init_omp(if_1_edge,n4,.false.)
    26:                 call nr_init_omp(diagn,n5,0)
    27:           c.....Mortar points indices are generated in two steps: first generate
    28:           c     them for all element vertices (corner points), then for conforming
    29:           c     edge and conforming face interiors. Each time a new mortar index
    30:           c     is generated for a mortar point, it is broadcast to all elements
    31:           c     sharing this mortar point.
    32:           
    33:           c.....VERTICES
    34:                 count=0
    35:           
    36:           c.....assign mortar point indices to element vertices
    37:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,sumcb,ij1,ij2,
    38:           c$OMP& cb,cb1,cb2,ntemp,ntemp1)
    39:           
    40: Y------>        do iel=1,nelt
    41: |         
    42: |         c.......first calculate how many new mortar indices will be generated for
    43: |         c       each element.
    44: |         
    45: |         c.......For each element, at least one vertex (vertex 8) will be new mortar
    46: |         c       point. All possible new mortar points will be on face 2,4 or 6. By
    47: |         c       checking the type of these three faces, we are able to tell
    48: |         c       how many new mortar vertex points will be generated in each element.
    49: |         
    50: |                 cb=cbc(6,iel)
    51: |                 cb1=cbc(4,iel)
    52: |                 cb2=cbc(2,iel)
    53: |         
    54: |         c.......For different combinations of the type of these three faces,
    55: |         c       we group them into 27 configurations.
    56: |         c       For different face types we assign the following integers:
    57: |         c              1 for type 2 or 3
    58: |         c              2 for type 0
    59: |         c              5 for type 1
    60: |         c       By summing these integers for faces 2,4 and 6, sumcb will have
    61: |         c       10 different numbers indicating 10 different combinations.
    62: |         
    63: |                 sumcb=0
    64: |                 if(cb.eq.2.or.cb.eq.3)then
    65: |                   sumcb=sumcb+1
    66: |                 elseif(cb.eq.0)then
    67: |                   sumcb=sumcb+2
    68: |                 elseif(cb.eq.1)then
    69: |                   sumcb=sumcb+5
    70: |                 end if
    71: |                 if(cb1.eq.2.or.cb1.eq.3)then
    72: |                   sumcb=sumcb+1
    73: |                 elseif(cb1.eq.0)then
    74: |                   sumcb=sumcb+2
    75: |                 elseif(cb1.eq.1)then
    76: |                   sumcb=sumcb+5
    77: |                 end if
    78: |                 if(cb2.eq.2.or.cb2.eq.3)then
    79: |                   sumcb=sumcb+1
    80: |                 elseif(cb2.eq.0)then
    81: |                   sumcb=sumcb+2
    82: |                 elseif(cb2.eq.1)then
    83: |                   sumcb=sumcb+5
    84: |                 end if
    85: |         
    86: |         c.......compute newc(iel)
    87: |         c       newc(iel) records how many new mortar indices will be generated
    88: |         c                 for element iel
    89: |         c       vassign(i,iel) records the element vertex of the i'th new mortar
    90: |         c                 vertex point for element iel. e.g. vassign(2,iel)=8 means
    91: |         c                 the 2nd new mortar vertex point generated on element
    92: |         c                 iel is iel's 8th vertex.
    93: |          
    94: |                 if(sumcb.eq.3)then
    95: |         c.......the three face types for face 2,4, and 6 are 2 2 2
    96: |                   newc(iel)=1
    97: |                   vassign(1,iel)=8
    98: |          
    99: |                 elseif(sumcb.eq.4)then
   100: |         c.......the three face types for face 2,4 and 6 are 2 2 0 (not
   101: |         c       necessarily in this order)
   102: |                   newc(iel)=2
   103: |                   if(cb.eq.0)then
   104: |                     vassign(1,iel)=4
   105: |                   elseif(cb1.eq.0)then
   106: |                     vassign(1,iel)=6
   107: |                   elseif(cb2.eq.0)then
   108: |                     vassign(1,iel)=7
   109: |                   end if
   110: |                   vassign(2,iel)=8
   111: |         
   112: |                 elseif(sumcb.eq.7)then
   113: |         c.......the three face types for face 2,4 and 6 are 2 2 1 (not
   114: |         c       necessarily in this order)
   115: |                   if(cb.eq.1)then
   116: |                     ij1=ijel(1,6,iel)
   117: |                     ij2=ijel(2,6,iel)
   118: |                     if(ij1.eq.1.and.ij2.eq.1)then
   119: |                       newc(iel)=2
   120: |                       vassign(1,iel)=4
   121: |                       vassign(2,iel)=8
   122: |                     elseif(ij1.eq.1.and.ij2.eq.2)then
   123: |                       ntemp=sje(1,1,6,iel)
   124: |       G               if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   125: |                         newc(iel)=1
   126: |                         vassign(1,iel)=8
   127: |                       else
   128: |                         newc(iel)=2
   129: |                         vassign(1,iel)=4
   130: |                         vassign(2,iel)=8
   131: |                       end if
   132: |                     elseif(ij1.eq.2.and.ij2.eq.1)then
   133: |                       ntemp=sje(1,1,6,iel)
   134: |       G               if(cbc(3,ntemp).eq.3.and.sje(1,1,3,ntemp).lt.iel)then
   135: |                         newc(iel)=1
   136: |                         vassign(1,iel)=8
   137: |                       else
   138: |                         newc(iel)=2
   139: |                         vassign(1,iel)=4
   140: |                         vassign(2,iel)=8
   141: |                       endif
   142: |                     else
   143: |                       newc(iel)=1
   144: |                       vassign(1,iel)=8
   145: |                     end if
   146: |                   elseif(cb1.eq.1)then
   147: |                     ij1=ijel(1,4,iel)
   148: |                     ij2=ijel(2,4,iel)
   149: |                     if(ij1.eq.1.and.ij2.eq.1)then
   150: |                       newc(iel)=2
   151: |                       vassign(1,iel)=6
   152: |                       vassign(2,iel)=8
   153: |                     elseif(ij1.eq.1.and.ij2.eq.2)then
   154: |                       ntemp=sje(1,1,4,iel)
   155: |       G               if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   156: |                         newc(iel)=1
   157: |                         vassign(1,iel)=8
   158: |                       else
   159: |                         newc(iel)=2
   160: |                         vassign(1,iel)=6
   161: |                         vassign(2,iel)=8
   162: |                       endif
   163: |                     elseif(ij1.eq.2.and.ij2.eq.1)then
   164: |                       ntemp=sje(1,1,4,iel)
   165: |       G               if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   166: |                         newc(iel)=1
   167: |                         vassign(1,iel)=8
   168: |                       else
   169: |                         newc(iel)=2
   170: |                         vassign(1,iel)=6
   171: |                         vassign(2,iel)=8
   172: |                       endif
   173: |                     else
   174: |                       newc(iel)=1
   175: |                       vassign(1,iel)=8
   176: |                     end if
   177: |         
   178: |                   elseif(cb2.eq.1)then
   179: |                     ij1=ijel(1,2,iel)
   180: |                     ij2=ijel(2,2,iel)
   181: |                     if(ij1.eq.1.and.ij2.eq.1)then
   182: |                       newc(iel)=2
   183: |                       vassign(1,iel)=7
   184: |                       vassign(2,iel)=8
   185: |                     elseif(ij1.eq.1.and.ij2.eq.2)then
   186: |                       ntemp=sje(1,1,2,iel)
   187: |       G               if(cbc(3,ntemp).eq.3.and.sje(1,1,3,ntemp).lt.iel)then
   188: |                         newc(iel)=1
   189: |                         vassign(1,iel)=8
   190: |                       else
   191: |                         newc(iel)=2
   192: |                         vassign(1,iel)=7
   193: |                         vassign(2,iel)=8
   194: |                       end if
   195: |         
   196: |                     elseif(ij1.eq.2.and.ij2.eq.1)then
   197: |                       ntemp=sje(1,1,2,iel)
   198: |       G               if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   199: |                         newc(iel)=1
   200: |                         vassign(1,iel)=8
   201: |                       else
   202: |                         newc(iel)=2
   203: |                         vassign(1,iel)=7
   204: |                         vassign(2,iel)=8
   205: |                       end if
   206: |                     else
   207: |                       newc(iel)=1
   208: |                       vassign(1,iel)=8
   209: |                     end if
   210: |                   end if
   211: |         
   212: |                 elseif(sumcb.eq.5)then
   213: |         c.......the three face types for face 2,4 and 6 are 2/3 0 0 (not
   214: |         c       necessarily in this order)
   215: |                   newc(iel)=4
   216: |                   if(cb.eq.2.or.cb.eq.3)then
   217: |                     vassign(1,iel)=5
   218: |                     vassign(2,iel)=6
   219: |                     vassign(3,iel)=7
   220: |                     vassign(4,iel)=8
   221: |                   elseif(cb1.eq.2.or.cb1.eq.3)then
   222: |                     vassign(1,iel)=3
   223: |                     vassign(2,iel)=4
   224: |                     vassign(3,iel)=7
   225: |                     vassign(4,iel)=8
   226: |                   elseif(cb2.eq.2.or.cb2.eq.3)then
   227: |                     vassign(1,iel)=2
   228: |                     vassign(2,iel)=4
   229: |                     vassign(3,iel)=6
   230: |                     vassign(4,iel)=8
   231: |                   end if
   232: |         
   233: |                 elseif(sumcb.eq.8)then
   234: |         c.......the three face types for face 2,4 and 6 are 2 0 1 (not
   235: |         c       necessarily in this order)
   236: |         
   237: |         c.........if face 2 of type 1
   238: |                   if(cb.eq.1)then
   239: |                     if(cb1.eq.2.or.cb1.eq.3)then
   240: |                       ij1=ijel(1,6,iel)
   241: |                       if(ij1.eq.1)then
   242: |                         newc(iel)=4
   243: |                         vassign(1,iel)=3
   244: |                         vassign(2,iel)=4
   245: |                         vassign(3,iel)=7
   246: |                         vassign(4,iel)=8
   247: |                       else
   248: |                         ntemp=sje(1,1,6,iel)
   249: |       G                 if(cbc(3,ntemp).eq.3.and.sje(1,1,3,ntemp).lt.iel)then
   250: |                           newc(iel)=2
   251: |                           vassign(1,iel)=7
   252: |                           vassign(2,iel)=8
   253: |                         else
   254: |                           newc(iel)=3
   255: |                           vassign(1,iel)=4
   256: |                           vassign(2,iel)=7
   257: |                           vassign(3,iel)=8
   258: |                         end if
   259: |                       end if
   260: |         
   261: |                     elseif(cb2.eq.2.or.cb2.eq.3)then
   262: |                       if(ijel(2,6,iel).eq.1)then
   263: |                         newc(iel)=4
   264: |                         vassign(1,iel)=2
   265: |                         vassign(2,iel)=4
   266: |                         vassign(3,iel)=6
   267: |                         vassign(4,iel)=8
   268: |                       else
   269: |                         ntemp=sje(1,1,6,iel)
   270: |       G                 if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   271: |                           newc(iel)=2
   272: |                           vassign(1,iel)=6
   273: |                           vassign(2,iel)=8
   274: |                         else
   275: |                           newc(iel)=3
   276: |                           vassign(1,iel)=4
   277: |                           vassign(2,iel)=6
   278: |                           vassign(3,iel)=8
   279: |                         end if
   280: |                       end if
   281: |                     end if
   282: |         
   283: |         c.........if face 4 of type 1
   284: |                   elseif(cb1.eq.1)then
   285: |                     if(cb.eq.2.or.cb.eq.3)then
   286: |                       ij1=ijel(1,4,iel)
   287: |                       ij2=ijel(2,4,iel)
   288: |         
   289: |                       if(ij1.eq.1.and.ij2.eq.1)then
   290: |                         ntemp=sje(1,1,4,iel)
   291: |       G                 if(cbc(2,ntemp).eq.3.and.sje(1,1,2,ntemp).lt.iel)then
   292: |                           newc(iel)=3
   293: |                           vassign(1,iel)=6
   294: |                           vassign(2,iel)=7
   295: |                           vassign(3,iel)=8
   296: |                         else
   297: |                           newc(iel)=4
   298: |                           vassign(1,iel)=5
   299: |                           vassign(2,iel)=6
   300: |                           vassign(3,iel)=7
   301: |                           vassign(4,iel)=8
   302: |                         end if
   303: |                       elseif(ij1.eq.1.and.ij2.eq.2)then
   304: |                         ntemp=sje(1,1,4,iel)
   305: |       G                 if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   306: |                           newc(iel)=3
   307: |                           vassign(1,iel)=5
   308: |                           vassign(2,iel)=7
   309: |                           vassign(3,iel)=8
   310: |                         else
   311: |                           newc(iel)=4
   312: |                           vassign(1,iel)=5
   313: |                           vassign(2,iel)=6
   314: |                           vassign(3,iel)=7
   315: |                           vassign(4,iel)=8
   316: |                         end if
   317: |                       elseif(ij1.eq.2.and.ij2.eq.1)then
   318: |                         ntemp=sje(1,1,4,iel)
   319: |       G                 if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   320: |                           newc(iel)=2
   321: |                           vassign(1,iel)=7
   322: |                           vassign(2,iel)=8
   323: |                         else
   324: |                           newc(iel)=3
   325: |                           vassign(1,iel)=6
   326: |                           vassign(2,iel)=7
   327: |                           vassign(3,iel)=8
   328: |                         end if
   329: |                       elseif(ij1.eq.2.and.ij2.eq.2)then
   330: |                         ntemp=sje(1,1,4,iel)
   331: |       G                 if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   332: |                           newc(iel)=2
   333: |                           vassign(1,iel)=7
   334: |                           vassign(2,iel)=8
   335: |                         else
   336: |                           newc(iel)=3
   337: |                           vassign(1,iel)=5
   338: |                           vassign(2,iel)=7
   339: |                           vassign(3,iel)=8
   340: |                         end if
   341: |                       end if
   342: |                     else
   343: |                       if(ijel(2,4,iel).eq.1)then
   344: |                         newc(iel)=4
   345: |                         vassign(1,iel)=2
   346: |                         vassign(2,iel)=4
   347: |                         vassign(3,iel)=6
   348: |                         vassign(4,iel)=8
   349: |                       else
   350: |                         ntemp=sje(1,1,4,iel)
   351: |       G                 if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   352: |                           newc(iel)=2
   353: |                           vassign(1,iel)=4
   354: |                           vassign(2,iel)=8
   355: |                         else
   356: |                           newc(iel)=3
   357: |                           vassign(1,iel)=4
   358: |                           vassign(2,iel)=6
   359: |                           vassign(3,iel)=8
   360: |                         end if
   361: |                       end if
   362: |                     endif
   363: |         c.........if face 6 of type 1
   364: |                   elseif(cb2.eq.1)then
   365: |                     if(cb.eq.2.or.cb.eq.3)then
   366: |                       if(ijel(1,2,iel).eq.1)then
   367: |                         newc(iel)=4
   368: |                         vassign(1,iel)=5
   369: |                         vassign(2,iel)=6
   370: |                         vassign(3,iel)=7
   371: |                         vassign(4,iel)=8
   372: |                       else
   373: |                         ntemp=sje(1,1,2,iel)
   374: |       G                 if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   375: |                           newc(iel)=2
   376: |                           vassign(1,iel)=6
   377: |                           vassign(2,iel)=8
   378: |                         else
   379: |                           newc(iel)=3
   380: |                           vassign(1,iel)=6
   381: |                           vassign(2,iel)=7
   382: |                           vassign(3,iel)=8
   383: |                         end if
   384: |                       end if
   385: |                     else
   386: |                       if(ijel(2,2,iel).eq.1)then
   387: |                         newc(iel)=4
   388: |                         vassign(1,iel)=3
   389: |                         vassign(2,iel)=4
   390: |                         vassign(3,iel)=7
   391: |                         vassign(4,iel)=8
   392: |                       else
   393: |                         ntemp=sje(1,1,2,iel)
   394: |       G                 if(cbc(3,ntemp).eq.3.and.sje(1,1,3,ntemp).lt.iel)then
   395: |                           newc(iel)=2
   396: |                           vassign(1,iel)=4
   397: |                           vassign(2,iel)=8
   398: |                         else
   399: |                           newc(iel)=3
   400: |                           vassign(1,iel)=4
   401: |                           vassign(2,iel)=7
   402: |                           vassign(3,iel)=8
   403: |                         end if
   404: |                       end if
   405: |                     end if
   406: |                   end if
   407: |         
   408: |                 elseif(sumcb.eq.11)then
   409: |         c.......the three face type for face 2,4 and 6 are 2 1 1(not
   410: |         c       necessarily in this order)
   411: |                   if(cb.eq.2.or.cb.eq.3)then
   412: |                     if(ijel(1,4,iel).eq.1)then
   413: |                       ntemp=sje(1,1,4,iel)
   414: |       G               if(cbc(2,ntemp).eq.3.and.sje(1,1,2,ntemp).lt.iel)then
   415: |                         newc(iel)=3
   416: |                         vassign(1,iel)=6
   417: |                         vassign(2,iel)=7
   418: |                         vassign(3,iel)=8
   419: |                       else
   420: |                         newc(iel)=4
   421: |                         vassign(1,iel)=5
   422: |                         vassign(2,iel)=6
   423: |                         vassign(3,iel)=7
   424: |                         vassign(4,iel)=8
   425: |                       end if
   426: |         
   427: |         c...........if ijel(1,4,iel)=2
   428: |                     else
   429: |                       ntemp=sje(1,1,2,iel)
   430: |       G               if(cbc(5,ntemp).eq.3.and.sje(1,1,5,ntemp).lt.iel)then
   431: |                         ntemp1=sje(1,1,4,iel)
   432: |       G                 if(cbc(5,ntemp1).eq.3.and.
   433: |              &             sje(1,1,5,ntemp1).lt.iel)then
   434: |                           newc(iel)=1
   435: |                           vassign(1,iel)=8
   436: |                         else
   437: |                           newc(iel)=2
   438: |                           vassign(1,iel)=6
   439: |                           vassign(2,iel)=8
   440: |                         end if
   441: |                       else
   442: |                         ntemp1=sje(1,1,4,iel)
   443: |       G                 if(cbc(5,ntemp1).eq.3.and.
   444: |              &             sje(1,1,5,ntemp1).lt.iel)then
   445: |                           newc(iel)=2
   446: |                           vassign(1,iel)=7
   447: |                           vassign(2,iel)=8
   448: |                         else
   449: |                           newc(iel)=3
   450: |                           vassign(1,iel)=6
   451: |                           vassign(2,iel)=7
   452: |                           vassign(3,iel)=8
   453: |                         end if
   454: |                       end if
   455: |                     end if
   456: |                   elseif(cb1.eq.2.or.cb1.eq.3)then
   457: |                     if(ijel(2,2,iel).eq.1)then
   458: |                       ntemp=sje(1,1,2,iel)
   459: |       G               if(cbc(6,ntemp).eq.3.and.sje(1,1,6,ntemp).lt.iel)then
   460: |                         newc(iel)=3
   461: |                         vassign(1,iel)=4
   462: |                         vassign(2,iel)=7
   463: |                         vassign(3,iel)=8
   464: |                       else
   465: |                         newc(iel)=4
   466: |                         vassign(1,iel)=3
   467: |                         vassign(2,iel)=4
   468: |                         vassign(3,iel)=7
   469: |                         vassign(4,iel)=8
   470: |                       end if
   471: |         c...........if ijel(2,2,iel)=2
   472: |                     else
   473: |                       ntemp=sje(1,1,2,iel)
   474: |       G               if(cbc(3,ntemp).eq.3.and.sje(1,1,3,ntemp).lt.iel)then
   475: |                         ntemp1=sje(1,1,6,iel)
   476: |       G                 if(cbc(3,ntemp1).eq.3.and.
   477: |              &            sje(1,1,3,ntemp1).lt.iel)then
   478: |                           newc(iel)=1
   479: |                           vassign(1,iel)=8
   480: |                         else
   481: |                           newc(iel)=2
   482: |                           vassign(1,iel)=4
   483: |                           vassign(2,iel)=8
   484: |                         end if
   485: |                       else
   486: |                         ntemp1=sje(1,1,6,iel)
   487: |       G                 if(cbc(3,ntemp1).eq.3.and.
   488: |              &            sje(1,1,3,ntemp1).lt.iel)then
   489: |                           newc(iel)=2
   490: |                           vassign(1,iel)=7
   491: |                           vassign(2,iel)=8
   492: |                         else
   493: |                           newc(iel)=3
   494: |                           vassign(1,iel)=4
   495: |                           vassign(2,iel)=7
   496: |                           vassign(3,iel)=8
   497: |                         end if
   498: |                       end if
   499: |                     end if
   500: |                   elseif(cb2.eq.2.or.cb2.eq.3)then
   501: |                     if(ijel(2,6,iel).eq.1)then
   502: |                       ntemp=sje(1,1,4,iel)
   503: |       G               if(cbc(6,ntemp).eq.3.and.sje(1,1,6,ntemp).lt.iel)then
   504: |                         newc(iel)=3
   505: |                         vassign(1,iel)=4
   506: |                         vassign(2,iel)=6
   507: |                         vassign(3,iel)=8
   508: |                       else
   509: |                         newc(iel)=4
   510: |                         vassign(1,iel)=2
   511: |                         vassign(2,iel)=4
   512: |                         vassign(3,iel)=6
   513: |                         vassign(4,iel)=8
   514: |                       end if
   515: |         c...........if ijel(2,6,iel)=2
   516: |                     else
   517: |                       ntemp=sje(1,1,4,iel)
   518: |       G               if(cbc(1,ntemp).eq.3.and.sje(1,1,1,ntemp).lt.iel)then
   519: |                         ntemp1=sje(1,1,6,iel)
   520: |       G                 if(cbc(1,ntemp1).eq.3.and.
   521: |              &            sje(1,1,1,ntemp1).lt.iel)then
   522: |                           newc(iel)=1
   523: |                           vassign(1,iel)=8
   524: |                         else
   525: |                           newc(iel)=2
   526: |                           vassign(1,iel)=4
   527: |                           vassign(2,iel)=8
   528: |                         end if
   529: |                       else
   530: |                         ntemp1=sje(1,1,6,iel)
   531: |       G                 if(cbc(1,ntemp1).eq.3.and.
   532: |              &              sje(1,1,1,ntemp1).lt.iel)then
   533: |                           newc(iel)=2
   534: |                           vassign(1,iel)=6
   535: |                           vassign(2,iel)=8
   536: |                         else
   537: |                           newc(iel)=3
   538: |                           vassign(1,iel)=4
   539: |                           vassign(2,iel)=6
   540: |                           vassign(3,iel)=8
   541: |                         end if
   542: |                       end if
   543: |                     end if
   544: |         
   545: |                   end if
   546: |          
   547: |                 elseif(sumcb.eq.6)then
   548: |         c.......the three face type for face 2,4 and 6 are 0 0 0(not
   549: |         c       necessarily in this order)
   550: |                   newc(iel)=8
   551: |                   vassign(1,iel)=1
   552: |                   vassign(2,iel)=2
   553: |                   vassign(3,iel)=3
   554: |                   vassign(4,iel)=4
   555: |                   vassign(5,iel)=5
   556: |                   vassign(6,iel)=6
   557: |                   vassign(7,iel)=7
   558: |                   vassign(8,iel)=8
   559: |         
   560: |                 elseif(sumcb.eq.9)then
   561: |         c.......the three face type for face 2,4 and 6 are 0 0 1(not
   562: |         c       necessarily in this order)
   563: |                   newc(iel)=7
   564: |                   vassign(1,iel)=2
   565: |                   vassign(2,iel)=3
   566: |                   vassign(3,iel)=4
   567: |                   vassign(4,iel)=5
   568: |                   vassign(5,iel)=6
   569: |                   vassign(6,iel)=7
   570: |                   vassign(7,iel)=8
   571: |         
   572: |                 elseif(sumcb.eq.12)then
   573: |         c.......the three face type for face 2,4 and 6 are 0 1 1(not
   574: |         c       necessarily in this order)
   575: |                   if(cb.eq.0)then
   576: |                     ntemp=sje(1,1,2,iel)
   577: |       G             if(cbc(4,ntemp).eq.3.and.sje(1,1,4,ntemp).lt.iel)then
   578: |                       newc(iel)=6
   579: |                       vassign(1,iel)=2
   580: |                       vassign(2,iel)=3
   581: |                       vassign(3,iel)=4
   582: |                       vassign(4,iel)=6
   583: |                       vassign(5,iel)=7
   584: |                       vassign(6,iel)=8
   585: |                     else
   586: |                       newc(iel)=7
   587: |                       vassign(1,iel)=2
   588: |                       vassign(2,iel)=3
   589: |                       vassign(3,iel)=4
   590: |                       vassign(4,iel)=5
   591: |                       vassign(5,iel)=6
   592: |                       vassign(6,iel)=7
   593: |                       vassign(7,iel)=8
   594: |                     end if
   595: |                   elseif(cb1.eq.0)then
   596: |                     newc(iel)=7
   597: |                     vassign(1,iel)=2
   598: |                     vassign(2,iel)=3
   599: |                     vassign(3,iel)=4
   600: |                     vassign(4,iel)=5
   601: |                     vassign(5,iel)=6
   602: |                     vassign(6,iel)=7
   603: |                     vassign(7,iel)=8
   604: |                   elseif(cb2.eq.0)then
   605: |                     ntemp=sje(1,1,4,iel)
   606: |       G             if(cbc(6,ntemp).eq.3.and.sje(1,1,6,ntemp).lt.iel)then
   607: |                       newc(iel)=6
   608: |                       vassign(1,iel)=3
   609: |                       vassign(2,iel)=4
   610: |                       vassign(3,iel)=5
   611: |                       vassign(4,iel)=6
   612: |                       vassign(5,iel)=7
   613: |                       vassign(6,iel)=8
   614: |                     else
   615: |                       newc(iel)=7
   616: |                       vassign(1,iel)=2
   617: |                       vassign(2,iel)=3
   618: |                       vassign(3,iel)=4
   619: |                       vassign(4,iel)=5
   620: |                       vassign(5,iel)=6
   621: |                       vassign(6,iel)=7
   622: |                       vassign(7,iel)=8
   623: |                     end if
   624: |                   end if
   625: |          
   626: |                 elseif(sumcb.eq.15)then
   627: |         c.......the three face type for face 2,4 and 6 are 1 1 1(not
   628: |         c       necessarily in this order)
   629: |                   ntemp=sje(1,1,4,iel)
   630: |                   ntemp1=sje(1,1,2,iel)
   631: |       G           if(cbc(6,ntemp).eq.3.and.sje(1,1,6,ntemp).lt.iel)then
   632: |       G             if(cbc(2,ntemp).eq.3.and.sje(1,1,2,ntemp).lt.iel)then
   633: |       G               if(cbc(6,ntemp1).eq.3.and.sje(1,1,6,ntemp1).lt.iel)then
   634: |                         newc(iel)=4
   635: |                         vassign(1,iel)=4
   636: |                         vassign(2,iel)=6
   637: |                         vassign(3,iel)=7
   638: |                         vassign(4,iel)=8
   639: |                       else
   640: |                         newc(iel)=5
   641: |                         vassign(1,iel)=3
   642: |                         vassign(2,iel)=4
   643: |                         vassign(3,iel)=6
   644: |                         vassign(4,iel)=7
   645: |                         vassign(5,iel)=8
   646: |                       end if
   647: |                     else
   648: |       G               if(cbc(6,ntemp1).eq.3.and.sje(1,1,6,ntemp1).lt.iel)then
   649: |                         newc(iel)=5
   650: |                         vassign(1,iel)=4
   651: |                         vassign(2,iel)=5
   652: |                         vassign(3,iel)=6
   653: |                         vassign(4,iel)=7
   654: |                         vassign(5,iel)=8
   655: |                       else
   656: |                         newc(iel)=6
   657: |                         vassign(1,iel)=3
   658: |                         vassign(2,iel)=4
   659: |                         vassign(3,iel)=5
   660: |                         vassign(4,iel)=6
   661: |                         vassign(5,iel)=7
   662: |                         vassign(6,iel)=8
   663: |                       end if
   664: |                     end if
   665: |                   else
   666: |       G             if(cbc(2,ntemp).eq.3.and.sje(1,1,2,ntemp).lt.iel)then
   667: |       G               if(cbc(6,ntemp1).eq.3.and.sje(1,1,6,ntemp1).lt.iel)then
   668: |                         newc(iel)=5
   669: |                         vassign(1,iel)=2
   670: |                         vassign(2,iel)=4
   671: |                         vassign(3,iel)=6
   672: |                         vassign(4,iel)=7
   673: |                         vassign(5,iel)=8
   674: |                       else
   675: |                         newc(iel)=6
   676: |                         vassign(1,iel)=2
   677: |                         vassign(2,iel)=3
   678: |                         vassign(3,iel)=4
   679: |                         vassign(4,iel)=6
   680: |                         vassign(5,iel)=7
   681: |                         vassign(6,iel)=8
   682: |                       end if
   683: |                     else
   684: |       G               if(cbc(6,ntemp1).eq.3.and.sje(1,1,6,ntemp1).lt.iel)then
   685: |                         newc(iel)=6
   686: |                         vassign(1,iel)=2
   687: |                         vassign(2,iel)=4
   688: |                         vassign(3,iel)=5
   689: |                         vassign(4,iel)=6
   690: |                         vassign(5,iel)=7
   691: |                         vassign(6,iel)=8
   692: |         
   693: |                       else
   694: |                         newc(iel)=7
   695: |                         vassign(1,iel)=2
   696: |                         vassign(2,iel)=3
   697: |                         vassign(3,iel)=4
   698: |                         vassign(4,iel)=5
   699: |                         vassign(5,iel)=6
   700: |                         vassign(6,iel)=7
   701: |                         vassign(7,iel)=8
   702: |                       end if
   703: |                     end if
   704: |                   end if
   705: |                 end if
   706: Y------         end do
   707:           c$OMP END PARALLEL DO
   708:           c.....end computing how many new mortar vertex points will be generated
   709:           c     on each element.
   710:           
   711:           c.....Compute (potentially in parallel) front(iel), which records how many
   712:           c     new mortar point indices are to be generated from element 1 to iel.
   713:           c     front(iel)=newc(1)+newc(2)+...+newc(iel)
   714:           
   715:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel)
   716: Y------>        do iel=1,nelt
   717: |                 front(iel)=newc(iel)
   718: Y------         end do
   719:           c$OMP END PARALLEL DO
   720:           
   721:                 call parallel_add(front)
   722:           
   723:           c.....On each element, generate new mortar point indices and assign them
   724:           c     to all elements sharing this mortar point. Note, if a mortar point
   725:           c     is shared by several elements, the mortar point index of it will only
   726:           c     be generated on the element with the lowest element index.
   727:           
   728:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,i,count)
   729: P------>        do iel=1,nelt
   730: |         
   731: |         c.......compute the starting vertex mortar point index in element iel
   732: |                 front(iel)=front(iel)-newc(iel)
   733: |         
   734: |+----->          do i=1,newc(iel)
   735: ||        c.........count is the new mortar index number, which will be assigned
   736: ||        c         to a vertex of iel and broadcast to all other elements sharing
   737: ||        c         this vertex point.
   738: ||                  count=front(iel)+i
   739: ||                  call mortar_vertex(vassign(i,iel),iel,count)
   740: |+-----           end do
   741: P------         end do
   742:           c$OMP END PARALLEL DO
   743:           
   744:           c.....nvertex records how many mortar indices are for element vertices.
   745:           c     It is used in the computation of the preconditioner.
   746:                 count=front(nelt)+newc(nelt)
   747:                 nvertex=count
   748:           
   749:           c.....CONFORMING EDGE AND FACE INTERIOR
   750:           
   751:           c.....find out how many new mortar point indices will be assigned to all
   752:           c.....conforming edges and all conforming face interiors on each element
   753:           
   754:                 n1=12*nelt
   755:                 n2=6*nelt
   756:           
   757:           c.....eassign(i,iel)=.true.   indicates that the i'th edge on iel will
   758:           c                             generate new mortar points.
   759:           c     ncon_edge(i,iel)=.true. indicates that the i'th edge on iel is
   760:           c                             nonconforming
   761:                 call l_init_omp(ncon_edge,n1,.false.)
   762:                 call l_init_omp(eassign,n1,.false.)
   763:           c.....fassign(i,iel)=.true. indicates that the i'th face of iel will
   764:           c                           generate new mortar points
   765:                 call l_init_omp(fassign,n2,.false.)
   766:           
   767:           c.....newe records how many new edges are to be assigned
   768:           c     diagn(1,n,iel) records the element index of neighbor element of iel,
   769:           c                    that shares edge n of iel
   770:           c     diagn(2,n,iel) records the neighbor element diagn(1,n,iel) shares which
   771:           c                    part of edge n of iel. diagn(2,n,iel)=1 refers to left
   772:           c                    or bottom half of the edge n, diagn(2,n,iel)=2 refers
   773:           c                    to the right or top part of edge n.
   774:           c     if_1_edge(n,iel)=.true. indicates that the size of iel is smaller than
   775:           c                    that of its neighbor connected, neighbored by edge n only
   776:           
   777:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,cb1,cb2,cb3,cb4,cb5
   778:           c$OMP& ,cb6,ntemp)
   779:           
   780: Y------>        do iel=1,nelt
   781: |                 newc(iel)=0
   782: |                 newe(iel)=0
   783: |                 newi(iel)=0
   784: |                 cb1=cbc(1,iel)
   785: |                 cb2=cbc(2,iel)
   786: |                 cb3=cbc(3,iel)
   787: |                 cb4=cbc(4,iel)
   788: |                 cb5=cbc(5,iel)
   789: |                 cb6=cbc(6,iel)
   790: |         
   791: |         c.......on face 6
   792: |         
   793: |                 if(cb6.eq.0)then
   794: |                   if(cb4.eq.0.or.cb4.eq.1)then
   795: |         c...........if face 6 is of type 0 and face 4 is of type 0 or type 1, the edge
   796: |         c           shared by face 4 and 6 (edge 11) will generate new mortar point
   797: |         c           indices.
   798: |                     newe(iel)=newe(iel)+1
   799: |                     eassign(11,iel)=.true.
   800: |                   end if
   801: |                   if(cb1.ne.3)then
   802: |         c...........if face 1 is of type 3, the edge shared by face 6 and 1 (edge 1)
   803: |         c           will generate new mortar points indices.
   804: |                     newe(iel)=newe(iel)+1
   805: |                     eassign(1,iel)=.true.
   806: |                   end if
   807: |                   if(cb3.ne.3)then
   808: |                     newe(iel)=newe(iel)+1
   809: |                     eassign(9,iel)=.true.
   810: |                   end if
   811: |                   if(cb2.eq.0.or.cb2.eq.1)then
   812: |                     newe(iel)=newe(iel)+1
   813: |                     eassign(5,iel)=.true.
   814: |                   end if
   815: |                 elseif(cb6.eq.1)then
   816: |                   if(cb4.eq.0)then
   817: |                     newe(iel)=newe(iel)+1
   818: |                     eassign(11,iel)=.true.
   819: |                   elseif(cb4.eq.1)then
   820: |         
   821: |         c...........If face 6 and face 4 both are of type 1, ntemp is the neighbor
   822: |         c           element on face 4.
   823: |                     ntemp=sje(1,1,4,iel)
   824: |         
   825: |         c...........if ntemp's face 6 is not noncoforming or the neighbor element
   826: |         c           of ntemp on face 6 has an element index larger than iel, the
   827: |         c           edge shared by face 6 and 4 (edge 11) will generate new mortar
   828: |         c           point indices.
   829: |                     if(cbc(6,ntemp).ne.3.or.sje(1,1,6,ntemp).gt.iel)then
   830: |         
   831: |                       newe(iel)=newe(iel)+1
   832: |                       eassign(11,iel)=.true.
   833: |         c.............if the face 6 of ntemp is of type 2
   834: |                       if(cbc(6,ntemp).eq.2)then
   835: |         c...............The neighbor element of iel, neighbored by edge 11, is
   836: |         c               sje(1,1,6,ntemp) (the neighbor element of ntemp on ntemp's
   837: |         c               face 6).
   838: |                         diagn(1,11,iel)=sje(1,1,6,ntemp)
   839: |         c...............The neighbor element of iel, neighbored by edge 11 shares
   840: |         c               the ijel(2,6,iel) part of edge 11 of iel
   841: |                         diagn(2,11,iel)=ijel(2,6,iel)
   842: |         c...............edge 10 of element sje(1,1,6,ntemp) (the neighbor element of
   843: |         c               ntemp on ntemp's face 6) is a nonconforming edge
   844: |                         ncon_edge(10,sje(1,1,6,ntemp))=.true.
   845: |         c...............if_1_edge(n,iel)=.true. indicates that iel is of a smaller
   846: |         c               size than its neighbor element, neighbored by edge n of iel only.
   847: |                         if_1_edge(11,iel)=.true.
   848: |                       endif
   849: |                       if(cbc(6,ntemp).eq.3.and.
   850: |              &          sje(1,1,6,ntemp).gt.iel)then
   851: |                         diagn(1,11,iel)=sje(2,ijel(2,6,iel),6,ntemp)
   852: |                       endif
   853: |                     end if
   854: |                   endif
   855: |         
   856: |                   if(cb1.eq.0)then
   857: |                     newe(iel)=newe(iel)+1
   858: |                     eassign(1,iel)=.true.
   859: |                   elseif(cb1.eq.1)then
   860: |                     ntemp=sje(1,1,1,iel)
   861: |                     if(cbc(6,ntemp).ne.3.or.sje(1,1,6,ntemp).gt.iel)then
   862: |                       newe(iel)=newe(iel)+1
   863: |                       eassign(1,iel)=.true.
   864: |                       if(cbc(6,ntemp).eq.2)then
   865: |                         diagn(1,1,iel)=sje(1,1,6,ntemp)
   866: |                         diagn(2,1,iel)=ijel(1,6,iel)
   867: |                         ncon_edge(7,sje(1,1,6,ntemp))=.true.
   868: |                         if_1_edge(1,iel)=.true.
   869: |                       endif
   870: |                       if(cbc(6,ntemp).eq.3.and.
   871: |              &          sje(1,1,6,ntemp).gt.iel)then
   872: |                         diagn(1,1,iel)=sje(ijel(1,6,iel),1,6,ntemp)
   873: |                       endif
   874: |                     end if
   875: |                   elseif(cb1.eq.2)then
   876: |                     if(ijel(2,6,iel).eq.2)then
   877: |                       ntemp=sje(1,1,1,iel)
   878: |                       if(cbc(6,ntemp).eq.1)then
   879: |                         newe(iel)=newe(iel)+1
   880: |                         eassign(1,iel)=.true.
   881: |         c.............if cbc(6,ntemp)=2
   882: |                       else
   883: |                         if(sje(1,1,6,ntemp).gt.iel)then
   884: |                           newe(iel)=newe(iel)+1
   885: |                           eassign(1,iel)=.true.
   886: |                           diagn(1,1,iel)=sje(1,1,6,ntemp)
   887: |                         end if
   888: |                       end if
   889: |                     else
   890: |                       newe(iel)=newe(iel)+1
   891: |                       eassign(1,iel)=.true.
   892: |                     end if
   893: |                   end if
   894: |         
   895: |                   if(cb3.eq.0)then
   896: |                     newe(iel)=newe(iel)+1
   897: |                     eassign(9,iel)=.true.
   898: |                   elseif(cb3.eq.1)then
   899: |                     ntemp=sje(1,1,3,iel)
   900: |                     if(cbc(6,ntemp).ne.3.or.sje(1,1,6,ntemp).gt.iel)then
   901: |                       newe(iel)=newe(iel)+1
   902: |                       eassign(9,iel)=.true.
   903: |                       if(cbc(6,ntemp).eq.2)then
   904: |                         diagn(1,9,iel)=sje(1,1,6,ntemp)
   905: |                         diagn(2,9,iel)=ijel(2,6,iel)
   906: |                         ncon_edge(12,sje(1,1,6,ntemp))=.true.
   907: |                         if_1_edge(9,iel)=.true.
   908: |                       endif
   909: |                       if(cbc(6,ntemp).eq.3.and.
   910: |              &           sje(1,1,6,ntemp).gt.iel)then
   911: |                         diagn(1,9,iel)=sje(2,ijel(2,6,iel),6,ntemp)
   912: |                       endif
   913: |                     end if
   914: |                   elseif(cb3.eq.2)then
   915: |                     if(ijel(1,6,iel).eq.2)then
   916: |                       ntemp=sje(1,1,3,iel)
   917: |                       if(cbc(6,ntemp).eq.1)then
   918: |                         newe(iel)=newe(iel)+1
   919: |                         eassign(9,iel)=.true.
   920: |         c.............if cbc(6,ntemp)=2
   921: |                       else
   922: |                         if(sje(1,1,6,ntemp).gt.iel)then
   923: |                           newe(iel)=newe(iel)+1
   924: |                           eassign(9,iel)=.true.
   925: |                           diagn(1,9,iel)=sje(1,1,6,ntemp)
   926: |                         end if
   927: |                       end if
   928: |                     else
   929: |                       newe(iel)=newe(iel)+1
   930: |                       eassign(9,iel)=.true.
   931: |                     end if
   932: |                   end if
   933: |         
   934: |                   if(cb2.eq.0)then
   935: |                     newe(iel)=newe(iel)+1
   936: |                     eassign(5,iel)=.true.
   937: |                   elseif(cb2.eq.1)then
   938: |                     ntemp=sje(1,1,2,iel)
   939: |                     if(cbc(6,ntemp).ne.3.or.sje(1,1,6,ntemp).gt.iel)then
   940: |                       newe(iel)=newe(iel)+1
   941: |                       eassign(5,iel)=.true.
   942: |                       if(cbc(6,ntemp).eq.2)then
   943: |                         diagn(1,5,iel)=sje(1,1,6,ntemp)
   944: |                         diagn(2,5,iel)=ijel(1,6,iel)
   945: |                         ncon_edge(3,sje(1,1,6,ntemp))=.true.
   946: |                         if_1_edge(5,iel)=.true.
   947: |                       endif
   948: |                       if(cbc(6,ntemp).eq.3.and.
   949: |              &          sje(1,1,6,ntemp).gt.iel)then
   950: |                         diagn(1,9,iel)=sje(2,ijel(2,6,iel),6,ntemp)
   951: |                       endif
   952: |                     endif
   953: |                   end if
   954: |                 end if
   955: |         
   956: |         c.......one face 4
   957: |                 if(cb4.eq.0)then
   958: |                   if(cb1.ne.3)then
   959: |                     newe(iel)=newe(iel)+1
   960: |                     eassign(4,iel)=.true.
   961: |                   endif
   962: |                   if(cb5.ne.3)then
   963: |                     newe(iel)=newe(iel)+1
   964: |                     eassign(12,iel)=.true.
   965: |                   endif
   966: |                   if(cb2.eq.0.or.cb2.eq.1)then
   967: |                     newe(iel)=newe(iel)+1
   968: |                     eassign(8,iel)=.true.
   969: |                   end if
   970: |          
   971: |                 elseif(cb4.eq.1)then
   972: |                   if(cb1.eq.2)then
   973: |                     if(ijel(2,4,iel).eq.1)then
   974: |                       newe(iel)=newe(iel)+1
   975: |                       eassign(4,iel)=.true.
   976: |                     else
   977: |                       ntemp=sje(1,1,4,iel)
   978: |                       if(cbc(1,ntemp).ne.3.or.sje(1,1,1,ntemp).gt.iel)then
   979: |                         newe(iel)=newe(iel)+1
   980: |                         eassign(4,iel)=.true.
   981: |                         if(cbc(1,ntemp).eq.3.and.
   982: |              &            sje(1,1,1,ntemp).gt.iel)then
   983: |                           diagn(1,4,iel)=sje(ijel(1,4,iel),2,1,ntemp)
   984: |                         endif
   985: |                       endif
   986: |                     end if
   987: |                   elseif(cb1.eq.0)then
   988: |                     newe(iel)=newe(iel)+1
   989: |                     eassign(4,iel)=.true.
   990: |                   elseif(cb1.eq.1)then
   991: |                     ntemp=sje(1,1,4,iel)
   992: |                     if(cbc(1,ntemp).ne.3.or.sje(1,1,1,ntemp).gt.iel)then
   993: |                       newe(iel)=newe(iel)+1
   994: |                       eassign(4,iel)=.true.
   995: |                       if(cbc(1,ntemp).eq.2)then
   996: |                         diagn(1,4,iel)=sje(1,1,1,ntemp)
   997: |                         diagn(2,4,iel)=ijel(1,4,iel)
   998: |                         ncon_edge(6,sje(1,1,1,ntemp))=.true.
   999: |                         if_1_edge(4,iel)=.true.
  1000: |                       endif
  1001: |                       if(cbc(1,ntemp).eq.3.and.
  1002: |              &          sje(1,1,1,ntemp).gt.iel)then
  1003: |                         diagn(1,4,iel)=sje(ijel(1,4,iel),2,1,ntemp)
  1004: |                       endif
  1005: |                     end if
  1006: |                   end if
  1007: |                   if(cb5.eq.2)then
  1008: |                     if(ijel(1,4,iel).eq.1)then
  1009: |                       newe(iel)=newe(iel)+1
  1010: |                       eassign(12,iel)=.true.
  1011: |                     else
  1012: |                       ntemp=sje(1,1,4,iel)
  1013: |                       if(cbc(5,ntemp).ne.3.or.sje(1,1,5,ntemp).gt.iel)then
  1014: |                         newe(iel)=newe(iel)+1
  1015: |                         eassign(12,iel)=.true.
  1016: |                         if(cbc(5,ntemp).eq.3.and.
  1017: |              &            sje(1,1,5,ntemp).gt.iel)then
  1018: |                           diagn(1,12,iel)=sje(2,ijel(2,4,iel),5,ntemp)
  1019: |                         endif
  1020: |                       endif
  1021: |                     end if
  1022: |                   elseif(cb5.eq.0)then
  1023: |                     newe(iel)=newe(iel)+1
  1024: |                     eassign(12,iel)=.true.
  1025: |                   elseif(cb5.eq.1)then
  1026: |                     ntemp=sje(1,1,4,iel)
  1027: |                     if(cbc(5,ntemp).ne.3.or.sje(1,1,5,ntemp).gt.iel)then
  1028: |                       newe(iel)=newe(iel)+1
  1029: |                       eassign(12,iel)=.true.
  1030: |                       if(cbc(5,ntemp).eq.2)then
  1031: |                         diagn(1,12,iel)=sje(1,1,5,ntemp)
  1032: |                         diagn(2,12,iel)=ijel(2,4,iel)
  1033: |                         ncon_edge(9,sje(1,1,5,ntemp))=.true.
  1034: |                         if_1_edge(12,iel)=.true.
  1035: |                       endif
  1036: |                       if(cbc(5,ntemp).eq.3.and.
  1037: |              &          sje(1,1,5,ntemp).gt.iel)then
  1038: |                         diagn(1,12,iel)=sje(2,ijel(2,4,iel),5,ntemp)
  1039: |                       endif
  1040: |                     end if
  1041: |                   end if
  1042: |                   if(cb2.eq.0)then
  1043: |                     newe(iel)=newe(iel)+1
  1044: |                     eassign(8,iel)=.true.
  1045: |                   elseif(cb2.eq.1)then
  1046: |                     ntemp=sje(1,1,4,iel)
  1047: |                     if(cbc(2,ntemp).ne.3.or.sje(1,1,2,ntemp).gt.iel)then
  1048: |                       newe(iel)=newe(iel)+1
  1049: |                       eassign(8,iel)=.true.
  1050: |                       if(cbc(2,ntemp).eq.2)then
  1051: |                         diagn(1,8,iel)=sje(1,1,2,ntemp)
  1052: |                         diagn(2,8,iel)=ijel(1,4,iel)
  1053: |                         ncon_edge(2,sje(1,1,2,ntemp))=.true.
  1054: |                         if_1_edge(8,iel)=.true.
  1055: |                       endif
  1056: |                       if(cbc(2,ntemp).eq.3.and.
  1057: |              &          sje(1,1,2,ntemp).gt.iel)then
  1058: |                         diagn(1,8,iel)=sje(ijel(1,4,iel),2,3,ntemp)
  1059: |                       endif
  1060: |                     endif
  1061: |                   end if
  1062: |                 end if
  1063: |         
  1064: |         c.......on face 2
  1065: |                 if(cb2.eq.0)then
  1066: |                   if(cb3.ne.3)then
  1067: |                     newe(iel)=newe(iel)+1
  1068: |                     eassign(6,iel)=.true.
  1069: |                   endif
  1070: |                   if(cb5.ne.3)then
  1071: |                     newe(iel)=newe(iel)+1
  1072: |                     eassign(7,iel)=.true.
  1073: |                   endif
  1074: |                 elseif(cb2.eq.1)then
  1075: |                   if(cb3.eq.2)then
  1076: |                     if(ijel(2,2,iel).eq.1)then
  1077: |                       newe(iel)=newe(iel)+1
  1078: |                       eassign(6,iel)=.true.
  1079: |                     else
  1080: |                       ntemp=sje(1,1,2,iel)
  1081: |                       if(cbc(3,ntemp).ne.3.or.
  1082: |              &          sje(1,1,3,ntemp).gt.iel)then
  1083: |                         newe(iel)=newe(iel)+1
  1084: |                         eassign(6,iel)=.true.
  1085: |                         if(cbc(3,ntemp).eq.3.and.
  1086: |              &            sje(1,1,3,ntemp).gt.iel)then
  1087: |                           diagn(1,6,iel)=sje(ijel(1,2,iel),2,3,ntemp)
  1088: |                         endif
  1089: |                       endif
  1090: |                     endif
  1091: |                   elseif(cb3.eq.0)then
  1092: |                     newe(iel)=newe(iel)+1
  1093: |                     eassign(6,iel)=.true.
  1094: |                   elseif(cb3.eq.1)then
  1095: |                     ntemp=sje(1,1,2,iel)
  1096: |                     if(cbc(3,ntemp).ne.3.or.sje(1,1,3,ntemp).gt.iel)then
  1097: |                       newe(iel)=newe(iel)+1
  1098: |                       eassign(6,iel)=.true.
  1099: |                       if(cbc(3,ntemp).eq.2)then
  1100: |                         diagn(1,6,iel)=sje(1,1,3,ntemp)
  1101: |                         diagn(2,6,iel)=ijel(1,2,iel)
  1102: |                         ncon_edge(4,sje(1,1,3,ntemp))=.true.
  1103: |                         if_1_edge(6,iel)=.true.
  1104: |                       endif
  1105: |                       if(cbc(3,ntemp).eq.3.and.
  1106: |              &          sje(1,1,3,ntemp).gt.iel)then
  1107: |                         diagn(1,6,iel)=sje(ijel(1,4,iel),2,3,ntemp)
  1108: |                       endif
  1109: |                     endif
  1110: |                   endif
  1111: |                   if(cb5.eq.2)then
  1112: |                     if(ijel(1,2,iel).eq.1)then
  1113: |                       newe(iel)=newe(iel)+1
  1114: |                       eassign(7,iel)=.true.
  1115: |                     else
  1116: |                       ntemp=sje(1,1,2,iel)
  1117: |                       if(cbc(5,ntemp).ne.3.or.sje(1,1,5,ntemp).gt.iel)then
  1118: |                         newe(iel)=newe(iel)+1
  1119: |                         eassign(7,iel)=.true.
  1120: |                         if(cbc(5,ntemp).eq.3.and.
  1121: |              &            sje(1,1,5,ntemp).gt.iel)then
  1122: |                           diagn(1,7,iel)=sje(ijel(2,2,iel),2,5,ntemp)
  1123: |                         endif
  1124: |                       endif
  1125: |                     endif
  1126: |                   elseif(cb5.eq.0)then
  1127: |                     newe(iel)=newe(iel)+1
  1128: |                     eassign(7,iel)=.true.
  1129: |                   elseif(cb5.eq.1)then
  1130: |                     ntemp=sje(1,1,2,iel)
  1131: |                     if(cbc(5,ntemp).ne.3.or.sje(1,1,5,ntemp).gt.iel)then
  1132: |                       newe(iel)=newe(iel)+1
  1133: |                       eassign(7,iel)=.true.
  1134: |                       if(cbc(5,ntemp).eq.2)then
  1135: |                         diagn(1,7,iel)=sje(1,1,5,ntemp)
  1136: |                         diagn(2,7,iel)=ijel(2,2,iel)
  1137: |                         ncon_edge(1,sje(1,1,5,ntemp))=.true.
  1138: |                         if_1_edge(7,iel)=.true.
  1139: |                       endif
  1140: |                       if(cbc(5,ntemp).eq.3.and.
  1141: |              &          sje(1,1,5,ntemp).gt.iel)then
  1142: |                         diagn(1,7,iel)=sje(2,ijel(2,4,iel),5,ntemp)
  1143: |                       endif
  1144: |                     endif
  1145: |                   endif
  1146: |                 end if
  1147: |         
  1148: |         c.......on face 1
  1149: |                 if(cb1.eq.1)then
  1150: |                   newe(iel)=newe(iel)+2
  1151: |                   eassign(2,iel)=.true.
  1152: |                   if(cb3.eq.1)then
  1153: |                     ntemp=sje(1,1,1,iel)
  1154: |                     if(cbc(3,ntemp).eq.2)then
  1155: |                       diagn(1,2,iel)=sje(1,1,3,ntemp)
  1156: |                       diagn(2,2,iel)=ijel(1,1,iel)
  1157: |                       ncon_edge(8,sje(1,1,3,ntemp))=.true.
  1158: |                       if_1_edge(2,iel)=.true.
  1159: |                     elseif(cbc(3,ntemp).eq.3)then
  1160: |                       diagn(1,2,iel)=sje(ijel(1,1,iel),1,3,ntemp)
  1161: |                     endif
  1162: |                   elseif(cb3.eq.2)then
  1163: |                     ntemp=sje(1,1,3,iel)
  1164: |                     if(ijel(2,1,iel).eq.2)then
  1165: |                       if(cbc(1,ntemp).eq.2)then
  1166: |                         diagn(1,2,iel)=sje(1,1,1,ntemp)
  1167: |                       end if
  1168: |                     endif
  1169: |                   end if
  1170: |         
  1171: |                   eassign(3,iel)=.true.
  1172: |                   if(cb5.eq.1)then
  1173: |                     ntemp=sje(1,1,1,iel)
  1174: |                     if(cbc(5,ntemp).eq.2)then
  1175: |                       diagn(1,3,iel)=sje(1,1,5,ntemp)
  1176: |                       diagn(2,3,iel)=ijel(2,1,iel)
  1177: |                       ncon_edge(5,sje(1,1,5,ntemp))=.true.
  1178: |                       if_1_edge(3,iel)=.true.
  1179: |                     elseif(cbc(5,ntemp).eq.3)then
  1180: |                       diagn(1,3,iel)=sje(ijel(2,1,iel),1,5,ntemp)
  1181: |                     endif
  1182: |                   elseif(cb5.eq.2)then
  1183: |                     ntemp=sje(1,1,5,iel)
  1184: |                     if(ijel(1,1,iel).eq.2)then
  1185: |                       if(cbc(1,ntemp).eq.2)then
  1186: |                         diagn(1,3,iel)=sje(1,1,1,ntemp)
  1187: |                       end if
  1188: |                     endif
  1189: |          
  1190: |                   end if
  1191: |                 elseif(cb1.eq.2)then
  1192: |                   if(cb3.eq.2)then
  1193: |                     ntemp=sje(1,1,1,iel)
  1194: |                     if(cbc(3,ntemp).ne.3)then
  1195: |                       newe(iel)=newe(iel)+1
  1196: |                       eassign(2,iel)=.true.
  1197: |                       if(cbc(3,ntemp).eq.2)then
  1198: |                         diagn(1,2,iel)=sje(1,1,3,ntemp)
  1199: |                       endif
  1200: |                     endif
  1201: |                   elseif(cb3.eq.0.or.cb3.eq.1)then
  1202: |                     newe(iel)=newe(iel)+1
  1203: |                     eassign(2,iel)=.true.
  1204: |                     if(cb3.eq.1)then
  1205: |                       ntemp=sje(1,1,1,iel)
  1206: |                       if(cbc(3,ntemp).eq.2)then
  1207: |                         diagn(1,2,iel)=sje(1,1,3,ntemp)
  1208: |                       endif
  1209: |                     endif
  1210: |                   end if
  1211: |                   if(cb5.eq.2)then
  1212: |                     ntemp=sje(1,1,1,iel)
  1213: |                     if(cbc(5,ntemp).ne.3)then
  1214: |                       newe(iel)=newe(iel)+1
  1215: |                       eassign(3,iel)=.true.
  1216: |                       if(cbc(5,ntemp).eq.2)then
  1217: |                         diagn(1,3,iel)=sje(1,1,5,ntemp)
  1218: |                       endif
  1219: |                     endif
  1220: |                   elseif(cb5.eq.0.or.cb5.eq.1)then
  1221: |                     newe(iel)=newe(iel)+1
  1222: |                     eassign(3,iel)=.true.
  1223: |                     if(cb5.eq.1)then
  1224: |                       ntemp=sje(1,1,1,iel)
  1225: |                       if(cbc(5,ntemp).eq.2)then
  1226: |                         diagn(1,3,iel)=sje(1,1,5,ntemp)
  1227: |                       endif
  1228: |                     endif
  1229: |                   end if
  1230: |                 elseif(cb1.eq.0)then
  1231: |                   if(cb3.ne.3)then
  1232: |                     newe(iel)=newe(iel)+1
  1233: |                     eassign(2,iel)=.true.
  1234: |                   endif
  1235: |                   if(cb5.ne.3)then
  1236: |                     newe(iel)=newe(iel)+1
  1237: |                     eassign(3,iel)=.true.
  1238: |                   endif
  1239: |                 endif
  1240: |         
  1241: |         c.......on face 3
  1242: |                 if(cb3.eq.1)then
  1243: |                   newe(iel)=newe(iel)+1
  1244: |                   eassign(10,iel)=.true.
  1245: |                   if(cb5.eq.1)then
  1246: |                     ntemp=sje(1,1,3,iel)
  1247: |                     if(cbc(5,ntemp).eq.2)then
  1248: |                       diagn(1,10,iel)=sje(1,1,5,ntemp)
  1249: |                       diagn(2,10,iel)=ijel(2,3,iel)
  1250: |                       ncon_edge(11,sje(1,1,5,ntemp))=.true.
  1251: |                       if_1_edge(10,iel)=.true.
  1252: |                     endif
  1253: |                   endif
  1254: |                   if(ijel(1,3,iel).eq.2)then
  1255: |                     ntemp=sje(1,1,3,iel)
  1256: |                     if(cbc(5,ntemp).eq.3)then
  1257: |                       diagn(1,10,iel)=sje(1,ijel(2,3,iel),5,ntemp)
  1258: |                     endif
  1259: |                   endif
  1260: |                 elseif(cb3.eq.2)then
  1261: |                   if(cb5.eq.2)then
  1262: |                     ntemp=sje(1,1,3,iel)
  1263: |                     if(cbc(5,ntemp).ne.3)then
  1264: |                       newe(iel)=newe(iel)+1
  1265: |                       eassign(10,iel)=.true.
  1266: |                       if(cbc(5,ntemp).eq.2)then
  1267: |                         diagn(1,10,iel)=sje(1,1,5,ntemp)
  1268: |                       endif
  1269: |                     endif
  1270: |                   elseif(cb5.eq.0.or.cb5.eq.1)then
  1271: |                     newe(iel)=newe(iel)+1
  1272: |                     eassign(10,iel)=.true.
  1273: |                     if(cb5.eq.1)then
  1274: |                       ntemp=sje(1,1,3,iel)
  1275: |                       if(cbc(5,ntemp).eq.2)then
  1276: |                         diagn(1,10,iel)=sje(1,1,5,ntemp)
  1277: |                       endif
  1278: |                     endif
  1279: |                   end if
  1280: |                 elseif(cb3.eq.0)then
  1281: |                   if(cb5.ne.3)then
  1282: |                     newe(iel)=newe(iel)+1
  1283: |                     eassign(10,iel)=.true.
  1284: |                   endif
  1285: |                 endif
  1286: |         
  1287: |         c       CONFORMING FACE INTERIOR
  1288: |         
  1289: |         c.......find how many new mortar point indices will be assigned
  1290: |         c       to face interiors on all faces on each element
  1291: |         
  1292: |         c.......newi record how many new face interior points will be assigned
  1293: |         
  1294: |         c.......on face 6
  1295: |                 if(cb6.eq.1.or.cb6.eq.0)then
  1296: |                   newi(iel)=newi(iel)+9
  1297: |                   fassign(6,iel)=.true.
  1298: |                 end if
  1299: |         c.......on face 4
  1300: |                 if(cb4.eq.1.or.cb4.eq.0)then
  1301: |                   newi(iel)=newi(iel)+9
  1302: |                   fassign(4,iel)=.true.
  1303: |                 end if
  1304: |         c.......on face 2
  1305: |                 if(cb2.eq.1.or.cb2.eq.0)then
  1306: |                   newi(iel)=newi(iel)+9
  1307: |                   fassign(2,iel)=.true.
  1308: |                 end if
  1309: |         c.......on face 1
  1310: |                 if(cb1.ne.3)then
  1311: |                   newi(iel)=newi(iel)+9
  1312: |                   fassign(1,iel)=.true.
  1313: |                 end if
  1314: |         c.......on face 3
  1315: |                 if(cb3.ne.3)then
  1316: |                   newi(iel)=newi(iel)+9
  1317: |                   fassign(3,iel)=.true.
  1318: |                 endif
  1319: |         c.......on face 5
  1320: |                 if(cb5.ne.3)then
  1321: |                   newi(iel)=newi(iel)+9
  1322: |                   fassign(5,iel)=.true.
  1323: |                 endif
  1324: |         
  1325: |         c.......newc is the total number of new mortar point indices
  1326: |         c       to be assigned to each element.
  1327: |                 newc(iel)=newe(iel)*3+newi(iel)
  1328: Y------         end do
  1329:           c$OMP END PARALLEL DO
  1330:           
  1331:           c.....Compute (potentially in parallel) front(iel), which records how
  1332:           c     many new mortar point indices are to be assigned (to conforming
  1333:           c     edges and conforming face interiors) from element 1 to iel.
  1334:           c     front(iel)=newc(1)+newc(2)+...+newc(iel)
  1335:           
  1336:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel)
  1337: Y------>        do iel=1,nelt
  1338: |                 front(iel)=newc(iel)
  1339: Y------         end do
  1340:           c$OMP END PARALLEL DO
  1341:           
  1342:                 call parallel_add(front)
  1343:           
  1344:           c.....nmor is the total number or mortar points
  1345:                 nmor=nvertex+front(nelt)
  1346:           
  1347:           c.....Generate (potentially in parallel) new mortar point indices on
  1348:           c     each conforming element face. On each face, first visit all
  1349:           c     conforming edges, and then the face interior.
  1350:           
  1351:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,count,i,cb1,ne,
  1352:           c$OMP& space,ie,edge_g,face2,ie2,ntemp,ii,jj,jface,cb,mor_v)
  1353: P------>        do iel=1,nelt
  1354: |                 front(iel)=front(iel)-newc(iel)
  1355: |                 count=nvertex+front(iel)
  1356: |+----->          do i=1,6
  1357: ||                  cb1=cbc(i,iel)
  1358: ||                  if (i.le.2) then
  1359: ||                    ne=4
  1360: ||                    space=1
  1361: ||                  elseif (i.le.4)then
  1362: ||                    ne=3
  1363: ||                    space=2
  1364: ||        
  1365: ||        c.........i loops over faces. Only 4 faces need to be examed for edge visit.
  1366: ||        c         On face 1, edge 1,2,3 and 4 will be visited. On face 2, edge 5,6,7
  1367: ||        c         and 8 will be visited. On face 3, edge 9 and 10 will be visited and
  1368: ||        c         on face 4, edge 11 and 12 will be visited. The 12 edges can be
  1369: ||        c         covered by four faces, there is no need to visit edges on face
  1370: ||        c         5 and 6.  So ne is set to be 0.
  1371: ||        c         However, i still needs to loop over 5 and 6, since the interiors
  1372: ||        c         of face 5 and 6 still need to be visited.
  1373: ||        
  1374: ||                  else
  1375: ||                    ne=0
  1376: ||                    space=1
  1377: ||                  end if
  1378: ||        
  1379: ||+---->            do ie=1,ne,space
  1380: |||                   edge_g=edgenumber(ie,i)
  1381: |||                   if(eassign(edge_g,iel))then
  1382: |||       c.............generate the new mortar points index, mor_v
  1383: |||                     call mor_assign(mor_v,count)
  1384: |||       c.............assign mor_v to local edge ie of face i on element iel
  1385: |||                     call mor_edge(ie,i,iel,mor_v)
  1386: |||       
  1387: |||       c.............Since this edge is shared by another face of element
  1388: |||       c             iel, assign mor_v to the corresponding edge on the other
  1389: |||       c             face also.
  1390: |||       
  1391: |||       c.............find the other face
  1392: |||                     face2=f_e_ef(ie,i)
  1393: |||       c.............find the local edge index of this edge on the other face
  1394: |||                     ie2=localedgenumber(face2,edge_g)
  1395: |||       c.............asssign mor_v  to local edge ie2 of face face2 on element iel
  1396: |||                     call mor_edge(ie2,face2,iel,mor_v)
  1397: |||       
  1398: |||       c.............There are some neighbor elements also sharing this edge. Assign
  1399: |||       c             mor_v to neighbor element, neighbored by face i.
  1400: |||                     if (cbc(i,iel).eq.2)then
  1401: |||                       ntemp=sje(1,1,i,iel)
  1402: |||                       call mor_edge(ie,jjface(i),ntemp,mor_v)
  1403: |||                       call mor_edge(op(ie2),face2,ntemp,mor_v)
  1404: |||                     end if
  1405: |||       
  1406: |||       c.............assign mor_v  to neighbor element neighbored by face face2
  1407: |||                     if (cbc(face2,iel).eq.2)then
  1408: |||                       ntemp=sje(1,1,face2,iel)
  1409: |||                       call mor_edge(ie2,jjface(face2),ntemp,mor_v)
  1410: |||                       call mor_edge(op(ie),i,ntemp,mor_v)
  1411: |||                     end if
  1412: |||       
  1413: |||       c.............assign mor_v to neighbor element sharing this edge
  1414: |||       
  1415: |||       c.............if the neighbor is of the same size of iel
  1416: |||                     if(.not.if_1_edge(edgenumber(ie,i),iel))then
  1417: |||                       if(diagn(1,edgenumber(ie,i),iel).ne.0)then
  1418: |||                         ntemp=diagn(1,edgenumber(ie,i),iel)
  1419: |||                         call mor_edge(op(ie2),jjface(face2),ntemp,mor_v)
  1420: |||                         call mor_edge(op(ie),jjface(i),ntemp,mor_v)
  1421: |||                       endif
  1422: |||       
  1423: |||       c.............if the neighbor has a size larger than iel's
  1424: |||                     else
  1425: |||                       if(diagn(1,edgenumber(ie,i),iel).ne.0)then
  1426: |||                         ntemp=diagn(1,edgenumber(ie,i),iel)
  1427: |||                         call mor_ne(mor_v,diagn(2,edgenumber(ie,i),iel),
  1428: |||            &            ie,i,ie2,face2,iel,ntemp)
  1429: |||                       end if
  1430: |||                     endif
  1431: |||        
  1432: |||                   endif
  1433: ||+----             end do
  1434: ||        
  1435: ||                  if(fassign(i,iel))then
  1436: ||        c...........generate new mortar points index in face interior.
  1437: ||        c           if face i is of type 2 or iel doesn't have a neighbor element,
  1438: ||        c           assign new mortar point indices to interior mortar points
  1439: ||        c           of face i of iel.
  1440: ||                    cb=cbc(i,iel)
  1441: ||                    if (cb.eq.1.or.cb.eq.0) then
  1442: ||*---->                do jj =2,lx1-1
  1443: |||*--->                  do ii=2,lx1-1
  1444: ||||                        count=count+1
  1445: ||||                        idmo(ii,jj,1,1,i,iel)=count
  1446: |||*---                   end do
  1447: ||*----                 end do
  1448: ||        
  1449: ||        c...........if face i is of type 2, assign new mortar point indices
  1450: ||        c           to iel as well as to the neighboring element on face i
  1451: ||                    elseif (cb.eq.2) then
  1452: ||                      if (idmo(2,2,1,1,i,iel).eq.0) then
  1453: ||                        ntemp=sje(1,1,i,iel)
  1454: ||                        jface = jjface(i)
  1455: ||*---->                  do jj =2,lx1-1
  1456: |||*--->                    do ii=2,lx1-1
  1457: ||||                          count=count+1
  1458: ||||                          idmo(ii,jj,1,1,i,iel)=count
  1459: ||||                          idmo(ii,jj,1,1,jface,ntemp)=count
  1460: |||*---                     end do
  1461: ||*----                   end do
  1462: ||                      end if
  1463: ||                    end if
  1464: ||                  end if
  1465: |+-----           end do
  1466: P------         end do
  1467:           c$OMP END  PARALLEL DO
  1468:           
  1469:            
  1470:           c.....for edges on nonconforming faces, copy the mortar points indices
  1471:           c     from neighbors.
  1472:           c$OMP PARALLEL DO DEFAULT(SHARED)
  1473:           c$OMP& PRIVATE(iel,i,cb,jface,iii,jjj,ntemp,ii,jj)
  1474: P------>        do iel=1,nelt
  1475: |+----->          do i=1,6
  1476: ||                  cb=cbc(i,iel)
  1477: ||                  if (cb.eq.3) then
  1478: ||        c...........edges
  1479: ||                    call edgecopy_s(i,iel)
  1480: ||                  end if
  1481: ||        
  1482: ||        c.........face interior
  1483: ||        
  1484: ||                  jface = jjface(i)
  1485: ||                  if (cb.eq.3) then
  1486: ||+---->              do iii=1,2
  1487: |||+--->                do jjj=1,2
  1488: ||||                      ntemp=sje(iii,jjj,i,iel)
  1489: ||||*-->                  do jj =1,lx1
  1490: |||||*->                    do ii=1,lx1
  1491: ||||||                        idmo(ii,jj,iii,jjj,i,iel)=
  1492: ||||||         &                         idmo(ii,jj,1,1,jface,ntemp)
  1493: |||||*-                     end do
  1494: ||||*--                   end do
  1495: ||||                      idmo(1,1,iii,jjj,i,iel)=idmo(1,1,1,1,jface,ntemp)
  1496: ||||                      idmo(lx1,1,iii,jjj,i,iel)=idmo(lx1,1,1,2,jface,ntemp)
  1497: ||||                      idmo(1,lx1,iii,jjj,i,iel)=idmo(1,lx1,2,1,jface,ntemp)
  1498: ||||                      idmo(lx1,lx1,iii,jjj,i,iel)=
  1499: ||||           &                         idmo(lx1,lx1,2,2,jface,ntemp)
  1500: |||+---                 end do
  1501: ||+----               end do
  1502: ||                  end if
  1503: |+-----           end do
  1504: P------         end do
  1505:           c$OMP END PARALLEL DO
  1506:                 return
  1507:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: GET_EMO
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1529: vec( 101): Vectorized loop.
  1530: opt(1084): Branch out of the loop inhibits optimization.
  1530: vec( 126): Idiom detected.: SEARCH


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: GET_EMO
FORMAT LIST

 LINE   LOOP      STATEMENT

  1510:                  subroutine get_emo(ie,n,ng)
  1511:           c-----------------------------------------------------------------
  1512:           c      This subroutine fills array emo.
  1513:           c      emo  records all elements sharing the same mortar point
  1514:           c                 (only applies to element vertices) .
  1515:           c      emo(1,i,n) gives the element ID of the i'th element sharing
  1516:           c                 mortar point n. (emo(1,i,n)=ie), ie is element
  1517:           c                 index.
  1518:           c      emo(2,i,n) gives the vertex index of mortar point n on this
  1519:           c                 element (emo(2,i,n)=ng), ng is the vertex index.
  1520:           c      nemo(n) records the total number of elements sharing mortar
  1521:           c                 point n.
  1522:           c-----------------------------------------------------------------
  1523:            
  1524:                  use ua_data
  1525:                  implicit none
  1526:           
  1527:                  integer ie, n, ntemp, i,ng
  1528:           
  1529: V------>         do i=1,nemo(n)
  1530: |                  if (emo(1,i,n).eq.ie) return
  1531: V------          end do
  1532:           
  1533:           c$     call omp_set_lock(tlock(n))
  1534:                  ntemp=nemo(n)+1
  1535:                  nemo(n)=ntemp
  1536:                  emo(1,ntemp,n)=ie
  1537:                  emo(2,ntemp,n)=ng
  1538:           c$     call omp_unset_lock(tlock(n))
  1539:           
  1540:                  return
  1541:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: IFSAME
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: IFSAME
FORMAT LIST

 LINE   LOOP      STATEMENT

  1544:                 logical function ifsame(ntemp,j,iel,i)
  1545:           c-----------------------------------------------------------------
  1546:           c     Check whether the i's vertex of element iel is at the same
  1547:           c     location as j's vertex of element ntemp.
  1548:           c-----------------------------------------------------------------
  1549:           
  1550:                 use ua_data
  1551:                 implicit none
  1552:           
  1553:                 integer iel, i, ntemp, j
  1554:           
  1555:                 ifsame=.false.
  1556:                 if (ntemp.eq.0 .or. iel.eq.0) return
  1557:                 if (xc(i,iel).eq.xc(j,ntemp).and.
  1558:                &    yc(i,iel).eq.yc(j,ntemp).and.
  1559:                &    zc(i,iel).eq.zc(j,ntemp)) then
  1560:                   ifsame=.true.
  1561:                 end if
  1562:           
  1563:                 return
  1564:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_ASSIGN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1577: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_ASSIGN
FORMAT LIST

 LINE   LOOP      STATEMENT

  1567:                 subroutine mor_assign(mor_v,count)
  1568:           c-----------------------------------------------------------------
  1569:           c     Assign three consecutive numbers for mor_v, which will
  1570:           c     be assigned to the three interior points of an edge as the
  1571:           c     mortar point indices.
  1572:           c-----------------------------------------------------------------
  1573:            
  1574:                 implicit none
  1575:                 integer mor_v(3),count,i
  1576:            
  1577: *------>        do i=1,3
  1578: |                 count=count+1
  1579: |                 mor_v(i)=count
  1580: *------         end do
  1581:           
  1582:                 return
  1583:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_EDGE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1600: opt(1418): Constant-length loop is expanded.
  1605: opt(1418): Constant-length loop is expanded.
  1610: opt(1418): Constant-length loop is expanded.
  1615: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_EDGE
FORMAT LIST

 LINE   LOOP      STATEMENT

  1586:                 subroutine mor_edge(ie,face,iel,mor_v)
  1587:           c-----------------------------------------------------------------
  1588:           c     Copy the mortar points index from mor_v to local
  1589:           c     edge ie of the face'th face on element iel.
  1590:           c     The edge is conforming.
  1591:           c-----------------------------------------------------------------
  1592:           
  1593:                 use ua_data
  1594:                 implicit none
  1595:           
  1596:                 integer ie,i,iel,mor_v(3),j,nn,face
  1597:           
  1598:                 if (ie.eq.1) then
  1599:                   j=1
  1600: *------>          do nn=2,lx1-1
  1601: |                   idmo(nn,j,1,1,face,iel)=mor_v(nn-1)
  1602: *------           end do
  1603:                 elseif (ie.eq.2) then
  1604:                   i=lx1
  1605: *------>          do nn=2,lx1-1
  1606: |                   idmo(i,nn,1,1,face,iel)=mor_v(nn-1)
  1607: *------           end do
  1608:                 elseif (ie.eq.3) then
  1609:                   j=lx1
  1610: *------>          do nn=2,lx1-1
  1611: |                   idmo(nn,j,1,1,face,iel)=mor_v(nn-1)
  1612: *------           end do
  1613:                 elseif (ie.eq.4) then
  1614:                   i=1
  1615: *------>          do nn=2,lx1-1
  1616: |                   idmo(i,nn,1,1,face,iel)=mor_v(nn-1)
  1617: *------           end do
  1618:                 end if
  1619:           
  1620:                 return
  1621:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: EDGECOPY_S
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1646: opt(1418): Constant-length loop is expanded.
  1650: opt(1418): Constant-length loop is expanded.
  1665: opt(1418): Constant-length loop is expanded.
  1671: opt(1418): Constant-length loop is expanded.
  1682: opt(1418): Constant-length loop is expanded.
  1688: opt(1418): Constant-length loop is expanded.
  1699: opt(1418): Constant-length loop is expanded.
  1704: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: EDGECOPY_S
FORMAT LIST

 LINE   LOOP      STATEMENT

  1624:                 subroutine edgecopy_s(face,iel)
  1625:           c------------------------------------------------------------
  1626:           c     Copy mortar points index on edges from neighbor elements
  1627:           c     to an element face of the 3rd type.
  1628:           c------------------------------------------------------------
  1629:           
  1630:                  use ua_data
  1631:                  implicit none
  1632:           
  1633:                  integer face, iel, ntemp1, ntemp2, ntemp3, ntemp4,
  1634:                &         edge_g, edge_l, face2, mor_s_v(4,2), i
  1635:           
  1636:           c......find four neighbors on this face (3rd type)
  1637:                  ntemp1=sje(1,1,face,iel)
  1638:                  ntemp2=sje(1,2,face,iel)
  1639:                  ntemp3=sje(2,1,face,iel)
  1640:                  ntemp4=sje(2,2,face,iel)
  1641:           
  1642:           c......local edge 1
  1643:           
  1644:           c......mor_s_v is the array of mortar indices to  be copied.
  1645:                  call nrzero(mor_s_v,4*2)
  1646: *------>         do i=2,lx1-1
  1647: |                   mor_s_v(i-1,1)=idmo(i,1,1,1,jjface(face),ntemp1)
  1648: *------          end do
  1649:                  mor_s_v(lx1-1,1)=idmo(lx1,1,1,2,jjface(face),ntemp1)
  1650: *------>         do i=1,lx1-1
  1651: |                   mor_s_v(i,2)=idmo(i,1,1,1,jjface(face),ntemp2)
  1652: *------          end do
  1653:           
  1654:           c......copy mor_s_v to local edge 1 on this face
  1655:                  call mor_s_e(1,face,iel,mor_s_v)
  1656:           
  1657:           c......copy mor_s_v to the corresponding edge on the other face sharing
  1658:           c      local edge 1
  1659:                  face2=f_e_ef(1,face)
  1660:                  edge_g=edgenumber(1,face)
  1661:                  edge_l=localedgenumber(face2,edge_g)
  1662:                  call mor_s_e(edge_l,face2,iel,mor_s_v)
  1663:           
  1664:           c......local edge 2
  1665: *------>         do i=2,lx1-1
  1666: |                   mor_s_v(i-1,1)=idmo(lx1,i,1,1,jjface(face),ntemp2)
  1667: *------          end do
  1668:                  mor_s_v(lx1-1,1)=idmo(lx1,lx1,2,2,jjface(face),ntemp2)
  1669:           
  1670:                  mor_s_v(1,2)=idmo(lx1,1,1,2,jjface(face),ntemp4)
  1671: *------>         do i=2,lx1-1
  1672: |                   mor_s_v(i,2)=idmo(lx1,i,1,1,jjface(face),ntemp4)
  1673: *------          end do
  1674:           
  1675:                  call mor_s_e(2,face,iel,mor_s_v)
  1676:                  face2=f_e_ef(2,face)
  1677:                  edge_g=edgenumber(2,face)
  1678:                  edge_l=localedgenumber(face2,edge_g)
  1679:                  call mor_s_e(edge_l,face2,iel,mor_s_v)
  1680:           
  1681:           c......local edge 3
  1682: *------>         do i=2,lx1-1
  1683: |                   mor_s_v(i-1,1)=idmo(i,lx1,1,1,jjface(face),ntemp3)
  1684: *------          end do
  1685:                  mor_s_v(lx1-1,1)=idmo(lx1,lx1,2,2,jjface(face),ntemp3)
  1686:           
  1687:                  mor_s_v(1,2)=idmo(1,lx1,2,1,jjface(face),ntemp4)
  1688: *------>         do i=2,lx1-1
  1689: |                   mor_s_v(i,2)=idmo(i,lx1,1,1,jjface(face),ntemp4)
  1690: *------          end do
  1691:           
  1692:                  call mor_s_e(3,face,iel,mor_s_v)
  1693:                  face2=f_e_ef(3,face)
  1694:                  edge_g=edgenumber(3,face)
  1695:                  edge_l=localedgenumber(face2,edge_g)
  1696:                  call mor_s_e(edge_l,face2,iel,mor_s_v)
  1697:           
  1698:           c......local edge 4
  1699: *------>         do i=2,lx1-1
  1700: |                   mor_s_v(i-1,1)=idmo(1,i,1,1,jjface(face),ntemp1)
  1701: *------          end do
  1702:                  mor_s_v(lx1-1,1)=idmo(1,lx1,2,1,jjface(face),ntemp1)
  1703:           
  1704: *------>         do i=1,lx1-1
  1705: |                   mor_s_v(i,2)=idmo(1,i,1,1,jjface(face),ntemp3)
  1706: *------          end do
  1707:           
  1708:                  call mor_s_e(4,face,iel,mor_s_v)
  1709:                  face2=f_e_ef(4,face)
  1710:                  edge_g=edgenumber(4,face)
  1711:                  edge_l=localedgenumber(face2,edge_g)
  1712:                  call mor_s_e(edge_l,face2,iel,mor_s_v)
  1713:           
  1714:                  return
  1715:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_S_E
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1730: opt(1418): Constant-length loop is expanded.
  1733: opt(1418): Constant-length loop is expanded.
  1737: opt(1418): Constant-length loop is expanded.
  1740: opt(1418): Constant-length loop is expanded.
  1744: opt(1418): Constant-length loop is expanded.
  1747: opt(1418): Constant-length loop is expanded.
  1751: opt(1418): Constant-length loop is expanded.
  1754: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_S_E
FORMAT LIST

 LINE   LOOP      STATEMENT

  1718:                  subroutine mor_s_e(n,face,iel,mor_s_v)
  1719:           c------------------------------------------------------------
  1720:           c      Copy mortar points index from mor_s_v to local edge n
  1721:           c      on face "face" of element iel. The edge is nonconforming.
  1722:           c------------------------------------------------------------
  1723:           
  1724:                  use ua_data
  1725:                  implicit none
  1726:           
  1727:                  integer n,face,iel,mor_s_v(4,2), i
  1728:           
  1729:                  if (n.eq.1) then
  1730: *------>           do i=2,lx1
  1731: |                    idmo(i,1,1,1,face,iel)=mor_s_v(i-1,1)
  1732: *------            end do
  1733: *------>           do i=1,lx1-1
  1734: |                    idmo(i,1,1,2,face,iel)=mor_s_v(i,2)
  1735: *------            end do
  1736:                  else if (n.eq.2) then
  1737: *------>           do i=2,lx1
  1738: |                   idmo(lx1,i,1,2,face,iel)=mor_s_v(i-1,1)
  1739: *------            end do
  1740: *------>           do i=1,lx1-1
  1741: |                   idmo(lx1,i,2,2,face,iel)=mor_s_v(i,2)
  1742: *------            end do
  1743:                  else if (n.eq.3) then
  1744: *------>           do i=2,lx1
  1745: |                    idmo(i,lx1,2,1,face,iel)=mor_s_v(i-1,1)
  1746: *------            end do
  1747: *------>           do i=1,lx1-1
  1748: |                    idmo(i,lx1,2,2,face,iel)=mor_s_v(i,2)
  1749: *------            end do
  1750:                  else if (n.eq.4) then
  1751: *------>           do i=2,lx1
  1752: |                    idmo(1,i,1,1,face,iel)=mor_s_v(i-1,1)
  1753: *------            end do
  1754: *------>           do i=1,lx1-1
  1755: |                    idmo(1,i,2,1,face,iel)=mor_s_v(i,2)
  1756: *------            end do
  1757:                  end if
  1758:                  return
  1759:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_S_E_NN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1777: opt(1418): Constant-length loop is expanded.
  1781: opt(1418): Constant-length loop is expanded.
  1787: opt(1418): Constant-length loop is expanded.
  1791: opt(1418): Constant-length loop is expanded.
  1797: opt(1418): Constant-length loop is expanded.
  1801: opt(1418): Constant-length loop is expanded.
  1807: opt(1418): Constant-length loop is expanded.
  1811: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_S_E_NN
FORMAT LIST

 LINE   LOOP      STATEMENT

  1762:                  subroutine mor_s_e_nn(n,face,iel,mor_s_v,nn)
  1763:           c------------------------------------------------------------
  1764:           c      Copy mortar point indices from mor_s_v to local edge n
  1765:           c      on face "face" of element iel. nn is the edge mortar index,
  1766:           c      which indicates that mor_s_v  corresponds to left/bottom or
  1767:           c      right/top part of the edge.
  1768:           c------------------------------------------------------------
  1769:           
  1770:                  use ua_data
  1771:                  implicit none
  1772:           
  1773:                  integer n,face,iel,mor_s_v(4), i,nn
  1774:           
  1775:                  if (n.eq.1) then
  1776:                    if(nn.eq.1)then
  1777: *------>              do i=2,lx1
  1778: |                       idmo(i,1,1,1,face,iel)=mor_s_v(i-1)
  1779: *------               end do
  1780:                    else
  1781: *------>             do i=1,lx1-1
  1782: |                      idmo(i,1,1,2,face,iel)=mor_s_v(i)
  1783: *------              end do
  1784:                    endif
  1785:                  else if (n.eq.2) then
  1786:                    if(nn.eq.1)then
  1787: *------>             do i=2,lx1
  1788: |                     idmo(lx1,i,1,2,face,iel)=mor_s_v(i-1)
  1789: *------              end do
  1790:                    else
  1791: *------>             do i=1,lx1-1
  1792: |                     idmo(lx1,i,2,2,face,iel)=mor_s_v(i)
  1793: *------              end do
  1794:                    endif
  1795:                  else if (n.eq.3) then
  1796:                    if(nn.eq.1)then
  1797: *------>             do i=2,lx1
  1798: |                      idmo(i,lx1,2,1,face,iel)=mor_s_v(i-1)
  1799: *------              end do
  1800:                    else
  1801: *------>             do i=1,lx1-1
  1802: |                     idmo(i,lx1,2,2,face,iel)=mor_s_v(i)
  1803: *------              end do
  1804:                    endif
  1805:                  else if (n.eq.4) then
  1806:                    if(nn.eq.1)then
  1807: *------>             do i=2,lx1
  1808: |                     idmo(1,i,1,1,face,iel)=mor_s_v(i-1)
  1809: *------              end do
  1810:                    else
  1811: *------>             do i=1,lx1-1
  1812: |                     idmo(1,i,2,1,face,iel)=mor_s_v(i)
  1813: *------              end do
  1814:                    endif
  1815:                  end if
  1816:                  return
  1817:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MORTAR_VERTEX
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1835: vec( 101): Vectorized loop.
  1844: opt(1418): Constant-length loop is expanded.
  1863: opt(1418): Constant-length loop is expanded.
  1932: vec( 110): Vectorization obstructive procedure reference.: IFSAME
  1932: vec( 103): Unvectorized loop.
  1947: opt(1025): Reference to this procedure inhibits optimization.: IFSAME


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MORTAR_VERTEX
FORMAT LIST

 LINE   LOOP      STATEMENT

  1821:                 subroutine mortar_vertex(i,iel,count)
  1822:           c---------------------------------------------------------------
  1823:           c     Assign mortar point index "count" to iel's i'th vertex
  1824:           c     and also to all elements sharing this vertex.
  1825:           c---------------------------------------------------------------
  1826:           
  1827:                 use ua_data
  1828:                 implicit none
  1829:           
  1830:                 integer i,iel,count,ntempx(8),ifntempx(8),lc_a(3),nnb(3),
  1831:                &        face_a(3),itemp,ntemp,ii, jj,j(3),
  1832:                &        iintempx(3),l,nbe, lc, temp
  1833:                 logical ifsame,if_temp
  1834:           
  1835: V------>        do l= 1,8
  1836: |                 ntempx(l)=0
  1837: |                 ifntempx(l)=0
  1838: V------         end do
  1839:           
  1840:           c.....face_a records the three faces sharing this vertex on iel.
  1841:           c     lc_a gives the local corner number of this vertex on each
  1842:           c     face in face_a.
  1843:           
  1844: *------>        do l=1,3
  1845: |                 face_a(l)=f_c(l,i)
  1846: |                 lc_a(l)=local_corner(i,face_a(l))
  1847: *------         end do
  1848:           
  1849:           c.....each vertex is shared by at most 8 elements.
  1850:           c     ntempx(j) gives the element index of a POSSIBLE element with its
  1851:           c               j'th  vertex is iel's i'th vertex
  1852:           c     ifntempx(i)=ntempx(i) means  ntempx(i) exists
  1853:           c     ifntempx(i)=0 means ntempx(i) does not exist.
  1854:           
  1855:                 ntempx(9-i)=iel
  1856:                 ifntempx(9-i)=iel
  1857:           
  1858:           c.....first find all elements sharing this vertex, ifntempx
  1859:           
  1860:           c.....find the three possible neighbors of iel, neighbored by faces
  1861:           c     listed in array face_a
  1862:           
  1863: *------>        do itemp= 1, 3
  1864: |         
  1865: |         c.......j(itemp) is the local corner number of this vertex on the
  1866: |         c       neighbor element on the corresponding face.
  1867: |                 j(itemp)=c_f(lc_a(itemp),jjface(face_a(itemp)))
  1868: |         
  1869: |         c.......iitempx(itemp) records the vertex index of i on the
  1870: |         c       neighbor element, neighborned by face_a(itemp)
  1871: |                 iintempx(itemp)=cal_intempx(lc_a(itemp),face_a(itemp))
  1872: |         
  1873: |         c.......ntemp refers the neighbor element
  1874: |                 ntemp=0
  1875: |         
  1876: |         c.......if the face is nonconforming, find out in which piece of the
  1877: |         c       mortar the vertex is located
  1878: |                 ii=cal_iijj(1,lc_a(itemp))
  1879: |                 jj=cal_iijj(2,lc_a(itemp))
  1880: |                 ntemp=sje(ii,jj,face_a(itemp),iel)
  1881: |         
  1882: |         c.......if the face is conforming
  1883: |                 if(ntemp.eq.0)then
  1884: |                   ntemp=sje(1,1,face_a(itemp),iel)
  1885: |         c.........find the possible neighbor
  1886: |                   ntempx(iintempx(itemp))=ntemp
  1887: |         c.........check whether this possible neighbor is a real neighbor or not
  1888: |                   if(ntemp.ne.0)then
  1889: |                     if(ifsame(ntemp,j(itemp),iel,i))then
  1890: |                       ifntempx(iintempx(itemp))=ntemp
  1891: |                     end if
  1892: |                   end if
  1893: |         
  1894: |         c.......if the face is nonconforming
  1895: |                 else
  1896: |                   if(ntemp.ne.0)then
  1897: |                     if(ifsame(ntemp,j(itemp),iel,i))then
  1898: |                       ifntempx(iintempx(itemp))=ntemp
  1899: |                       ntempx(iintempx(itemp))=ntemp
  1900: |                     end if
  1901: |                   end if
  1902: |                 end if
  1903: *------         end do
  1904:           
  1905:           c.....find the possible three neighbors, neighbored by an edge only
  1906: +------>        do l=1,3
  1907: |         
  1908: |         c.....find first existing neighbor of any of the faces in array face_a
  1909: |                 if_temp=.false.
  1910: |                 if(l.eq.1)then
  1911: |                   if_temp=.true.
  1912: |                 elseif(l.eq.2)then
  1913: |                   if(ifntempx(iintempx(l-1)).eq.0)then
  1914: |                     if_temp=.true.
  1915: |                   end if
  1916: |                 elseif(l.eq.3)then
  1917: |                   if(ifntempx(iintempx(l-1)).eq.0
  1918: |              &       .and.ifntempx(iintempx(l-2)).eq.0) then
  1919: |                     if_temp=.true.
  1920: |                   end if
  1921: |                 end if
  1922: |          
  1923: |                 if(if_temp)then
  1924: |                   if (ifntempx(iintempx(l)).ne.0) then
  1925: |                     nbe=ifntempx(iintempx(l))
  1926: |         c...........if 1st neighor exists, check the neighbor's two neighbors in
  1927: |         c           the other two directions.
  1928: |         c           e.g. if l=1, check directions 2 and 3,i.e. itemp=2,3,1
  1929: |         c           if l=2, itemp=3,1,-2
  1930: |         c           if l=3, itemp=1,2,1
  1931: |         c
  1932: |+----->              do itemp=face_l1(l),face_l2(l),face_ld(l)
  1933: ||        c.............lc is the local corner number of this vertex on face face_a(itemp)
  1934: ||        c             on the neighbor element of iel, neighbored by a face face_a(l)
  1935: ||                      lc=local_corner(j(l),face_a(itemp))
  1936: ||        c.............temp is the vertex index of this vertex on the neighbor element
  1937: ||        c             neighbored by an edge
  1938: ||                      temp=cal_intempx(lc,face_a(itemp))
  1939: ||                      ii=cal_iijj(1,lc)
  1940: ||                      jj=cal_iijj(2,lc)
  1941: ||                      ntemp=sje(ii,jj,face_a(itemp),nbe)
  1942: ||        
  1943: ||        c.............if the face face_a(itemp) is conforming
  1944: ||                      if(ntemp.eq.0)then
  1945: ||                        ntemp=sje(1,1,face_a(itemp),nbe)
  1946: ||                        if(ntemp.ne.0)then
  1947: ||                          if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  1948: ||             &            nbe,j(l)))then
  1949: ||                            ntempx(temp)=ntemp
  1950: ||                            ifntempx(temp)=ntemp
  1951: ||        c...................nnb(itemp) records the neighbor element neighbored by an
  1952: ||        c                   edge only
  1953: ||                            nnb(itemp)=ntemp
  1954: ||                          end if
  1955: ||                        end if
  1956: ||        
  1957: ||        c.............if the face face_a(itemp) is nonconforming
  1958: ||                      else
  1959: ||                        if(ntemp.ne.0)then
  1960: ||                          if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  1961: ||             &              nbe,j(l)))then
  1962: ||                            ntempx(temp)=ntemp
  1963: ||                            ifntempx(temp)=ntemp
  1964: ||                            nnb(itemp)=ntemp
  1965: ||                          end if
  1966: ||                        end if
  1967: ||                      end if
  1968: |+-----               end do
  1969: |         
  1970: |         c...........check the last neighbor element, neighbored by an edge
  1971: |         
  1972: |         c...........ifntempx(iintempx(l)) has been visited in the above, now
  1973: |         c           check another neighbor element(nbe) neighbored by a face
  1974: |         
  1975: |         c...........if the neighbor element is neighbored by face
  1976: |         c           face_a(face_l1(l)) exists
  1977: |                     if(ifntempx(iintempx(face_l1(l))).ne.0)then
  1978: |                       nbe=ifntempx(iintempx(face_l1(l)))
  1979: |         c.............itemp is the last direction other than l and face_l1(l)
  1980: |                       itemp=face_l2(l)
  1981: |                       lc=local_corner(j(face_l1(l)),face_a(itemp))
  1982: |                       temp=cal_intempx(lc,face_a(itemp))
  1983: |                       ii=cal_iijj(1,lc)
  1984: |                       jj=cal_iijj(2,lc)
  1985: |         
  1986: |         c.............ntemp records the last neighbor element neighbored by an edge
  1987: |         c             with element iel
  1988: |                       ntemp=sje(ii,jj,face_a(itemp),nbe)
  1989: |         c.............if conforming
  1990: |                       if(ntemp.eq.0)then
  1991: |                         ntemp=sje(1,1,face_a(itemp),nbe)
  1992: |                         if(ntemp.ne.0)then
  1993: |                           if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  1994: |              &              nbe,j(face_l1(l))))then
  1995: |                             ntempx(temp)=ntemp
  1996: |                             ifntempx(temp)=ntemp
  1997: |                             nnb(l)=ntemp
  1998: |                           end if
  1999: |                         end if
  2000: |         c.............if nonconforming
  2001: |                       else
  2002: |                         if(ntemp.ne.0)then
  2003: |                           if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  2004: |              &            nbe,j(face_l1(l))))then
  2005: |                             ntempx(temp)=ntemp
  2006: |                             ifntempx(temp)=ntemp
  2007: |                             nnb(l)=ntemp
  2008: |                           end if
  2009: |                         end if
  2010: |                       end if
  2011: |         
  2012: |         c...........if the neighbor element neighbored by face face_a(face_l2(l))
  2013: |         c           does not exist
  2014: |                     elseif(ifntempx(iintempx(face_l2(l))).ne.0)then
  2015: |                       nbe=ifntempx(iintempx(face_l2(l)))
  2016: |                       itemp=face_l1(l)
  2017: |                       lc=local_corner(j(face_l2(l)),face_a(itemp))
  2018: |                       temp=cal_intempx(lc,face_a(itemp))
  2019: |                       ii=cal_iijj(1,lc)
  2020: |                       jj=cal_iijj(2,lc)
  2021: |                       ntemp=sje(ii,jj,face_a(itemp),nbe)
  2022: |                       if(ntemp.eq.0)then
  2023: |                         ntemp=sje(1,1,face_a(itemp),nbe)
  2024: |                         if(ntemp.ne.0)then
  2025: |                           if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  2026: |              &            nbe,j(face_l2(l))))then
  2027: |                             ntempx(temp)=ntemp
  2028: |                             ifntempx(temp)=ntemp
  2029: |                             nnb(l)=ntemp
  2030: |                           end if
  2031: |                         end if
  2032: |                       else
  2033: |                         if(ntemp.ne.0.)then
  2034: |                           if(ifsame(ntemp,c_f(lc,jjface(face_a(itemp))),
  2035: |              &            nbe,j(face_l2(l))))then
  2036: |                             ntempx(temp)=ntemp
  2037: |                             ifntempx(temp)=ntemp
  2038: |                             nnb(l)=ntemp
  2039: |                           end if
  2040: |                         end if
  2041: |                       end if
  2042: |                     endif
  2043: |                   endif
  2044: |                 end if
  2045: +------         end do
  2046:           
  2047:           c.....check the neighbor element, neighbored by a vertex only
  2048:           
  2049:           c.....nnb are the three possible neighbor elements neighbored by an edge
  2050:           
  2051:                 nnb(1)=ifntempx(cal_nnb(1,i))
  2052:                 nnb(2)=ifntempx(cal_nnb(2,i))
  2053:                 nnb(3)=ifntempx(cal_nnb(3,i))
  2054:                 ntemp=0
  2055:           
  2056:           c.....the neighbor element neighbored by a vertex must be a neighbor of
  2057:           c     a valid(nonzero) nnb(i), neighbored by a face
  2058:           
  2059:                 if(nnb(1).ne.0)then
  2060:                   lc=oplc(local_corner(i,face_a(3)))
  2061:                   ii=cal_iijj(1,lc)
  2062:                   jj=cal_iijj(2,lc)
  2063:           c.......ntemp records the neighbor of iel, neighbored by vertex i
  2064:                   ntemp=sje(ii,jj,face_a(3),nnb(1))
  2065:           c.......temp is the vertex index of i on ntemp
  2066:                   temp=cal_intempx(lc,face_a(3))
  2067:                   if(ntemp.eq.0)then
  2068:                     ntemp=sje(1,1,face_a(3),nnb(1))
  2069:                     if(ntemp.ne.0)then
  2070:                       if(ifsame(ntemp,c_f(lc,jjface(face_a(3))),
  2071:                &         iel,i))then
  2072:                         ntempx(temp)=ntemp
  2073:                         ifntempx(temp)=ntemp
  2074:                       end if
  2075:                     end if
  2076:                   else
  2077:                     if(ntemp.ne.0)then
  2078:                       if(ifsame(ntemp,c_f(lc,jjface(face_a(3))),
  2079:                &         iel,i))then
  2080:                         ntempx(temp)=ntemp
  2081:                         ifntempx(temp)=ntemp
  2082:                       end if
  2083:                     end if
  2084:                   end if
  2085:                 elseif(nnb(2).ne.0)then
  2086:                   lc=oplc(local_corner(i,face_a(1)))
  2087:                   ii=cal_iijj(1,lc)
  2088:                   jj=cal_iijj(2,lc)
  2089:                   ntemp=sje(ii,jj,face_a(1),nnb(2))
  2090:                   temp=cal_intempx(lc,face_a(1))
  2091:                   if(ntemp.eq.0)then
  2092:                     ntemp=sje(1,1,face_a(1),nnb(2))
  2093:                     if(ntemp.ne.0)then
  2094:                       if(ifsame(ntemp,
  2095:                &        c_f(lc,jjface(face_a(1))),iel,i))then
  2096:                         ntempx(temp)=ntemp
  2097:                         ifntempx(temp)=ntemp
  2098:                       end if
  2099:                     end if
  2100:                   else
  2101:                     if(ntemp.ne.0)then
  2102:                       if(ifsame(ntemp,
  2103:                &      c_f(lc,jjface(face_a(1))),iel,i))then
  2104:                         ntempx(temp)=ntemp
  2105:                         ifntempx(temp)=ntemp
  2106:                       end if
  2107:                     end if
  2108:                   end if
  2109:                 elseif(nnb(3).ne.0)then
  2110:                   lc=oplc(local_corner(i,face_a(2)))
  2111:                   ii=cal_iijj(1,lc)
  2112:                   jj=cal_iijj(2,lc)
  2113:                   ntemp=sje(ii,jj,face_a(2),nnb(3))
  2114:                   temp=cal_intempx(lc, face_a(2))
  2115:                   if(ntemp.eq.0)then
  2116:                     ntemp=sje(1,1,face_a(2),nnb(3))
  2117:                     if(ntemp.ne.0)then
  2118:                       if(ifsame(ntemp,
  2119:                &         c_f(lc,jjface(face_a(2))),iel,i))then
  2120:                         ifntempx(temp)=ntemp
  2121:                         ntempx(temp)=ntemp
  2122:                       end if
  2123:                     end if
  2124:                   else
  2125:                     if(ntemp.ne.0)then
  2126:                       if(ifsame(ntemp,
  2127:                &        c_f(lc,jjface(face_a(2))),iel,i))then
  2128:                         ifntempx(temp)=ntemp
  2129:                         ntempx(temp)=ntemp
  2130:                       end if
  2131:                     end if
  2132:                   end if
  2133:                 end if
  2134:           
  2135:           c.....ifntempx records all elements sharing this vertex, assign count
  2136:           c     to all these elements.
  2137:           
  2138:                 if (ifntempx(1).ne.0) then
  2139:                   idmo(lx1,lx1,2,2,1,ntempx(1))=count
  2140:                   idmo(lx1,lx1,2,2,3,ntempx(1))=count
  2141:                   idmo(lx1,lx1,2,2,5,ntempx(1))=count
  2142:                   call get_emo(ntempx(1),count,8)
  2143:                 end if
  2144:           
  2145:                 if (ifntempx(2).ne.0) then
  2146:                   idmo(lx1,lx1,2,2,2,ntempx(2))=count
  2147:                   idmo(1,lx1,2,1,3,ntempx(2))=count
  2148:                   idmo(1,lx1,2,1,5,ntempx(2))=count
  2149:                   call get_emo(ntempx(2),count,7)
  2150:                 end if
  2151:           
  2152:                 if (ifntempx(3).ne.0) then
  2153:                   idmo(1,lx1,2,1,1,ntempx(3))=count
  2154:                   idmo(lx1,lx1,2,2,4,ntempx(3))=count
  2155:                   idmo(lx1,1,1,2,5,ntempx(3))=count
  2156:                   call get_emo(ntempx(3),count,6)
  2157:                 end if
  2158:                 if (ifntempx(4).ne.0) then
  2159:                   idmo(1,lx1,2,1,2,ntempx(4))=count
  2160:                   idmo(1,lx1,2,1,4,ntempx(4))=count
  2161:                   idmo(1,1,1,1,5,ntempx(4))=count
  2162:                   call get_emo(ntempx(4),count,5)
  2163:                 end if
  2164:           
  2165:                 if (ifntempx(5).ne.0) then
  2166:                   idmo(lx1,1,1,2,1,ntempx(5))=count
  2167:                   idmo(lx1,1,1,2,3,ntempx(5))=count
  2168:                   idmo(lx1,lx1,2,2,6,ntempx(5))=count
  2169:                   call get_emo(ntempx(5),count,4)
  2170:                 end if
  2171:           
  2172:           
  2173:                 if (ifntempx(6).ne.0) then
  2174:                   idmo(lx1,1,1,2,2,ntempx(6))=count
  2175:                   idmo(1,1,1,1,3,ntempx(6))=count
  2176:                   idmo(1,lx1,2,1,6,ntempx(6))=count
  2177:                   call get_emo(ntempx(6),count,3)
  2178:                 end if
  2179:           
  2180:                 if (ifntempx(7).ne.0) then
  2181:                   idmo(1,1,1,1,1,ntempx(7))=count
  2182:                   idmo(lx1,1,1,2,4,ntempx(7))=count
  2183:                   idmo(lx1,1,1,2,6,ntempx(7))=count
  2184:                   call get_emo(ntempx(7),count,2)
  2185:                 end if
  2186:           
  2187:                 if (ifntempx(8).ne.0) then
  2188:                   idmo(1,1,1,1,2,ntempx(8))=count
  2189:                   idmo(1,1,1,1,4,ntempx(8))=count
  2190:                   idmo(1,1,1,1,6,ntempx(8))=count
  2191:                   call get_emo(ntempx(8),count,1)
  2192:                 end if
  2193:           
  2194:                 return
  2195:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_NE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  2222: opt(1418): Constant-length loop is expanded.
  2228: opt(1418): Constant-length loop is expanded.
  2235: opt(1418): Constant-length loop is expanded.
  2241: opt(1418): Constant-length loop is expanded.
  2248: opt(1418): Constant-length loop is expanded.
  2254: opt(1418): Constant-length loop is expanded.
  2261: opt(1418): Constant-length loop is expanded.
  2267: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: mason.f

PROCEDURE NAME: MOR_NE
FORMAT LIST

 LINE   LOOP      STATEMENT

  2199:                 subroutine mor_ne(mor_v,nn,edge,face,edge2,face2,ntemp,iel)
  2200:           c---------------------------------------------------------------
  2201:           c     Copy the mortar points index  (mor_v + vertex mortar point) from
  2202:           c     edge'th local edge on face'th face on element ntemp to iel.
  2203:           c     ntemp is iel's neighbor, neighbored by this edge only.
  2204:           c     This subroutine is for the situation that iel is of larger
  2205:           c     size than ntemp.
  2206:           c     face, face2 are face indices
  2207:           c     edge and edge2 are local edge numbers of this edge on face and face2
  2208:           c     nn is edge motar index, which indicate whether this edge
  2209:           c     corresponds to the left/bottom or right/top part of the edge
  2210:           c     on iel.
  2211:           c---------------------------------------------------------------
  2212:           
  2213:                 use ua_data
  2214:                 implicit none
  2215:           
  2216:                 integer mor_v(3),nn,edge,face,edge2,face2,ntemp,iel, i,
  2217:                &mor_s_v(4)
  2218:           
  2219:           c.....get mor_s_v which is the mor_v + vertex mortar
  2220:                 if (edge.eq.3) then
  2221:                   if(nn.eq.1)then
  2222: *------>            do i=2,lx1-1
  2223: |                     mor_s_v(i-1)=mor_v(i-1)
  2224: *------             end do
  2225:                     mor_s_v(4)=idmo(lx1,lx1,2,2,face,ntemp)
  2226:                   else
  2227:                     mor_s_v(1)=idmo(1,lx1,2,1,face,ntemp)
  2228: *------>            do i=2,lx1-1
  2229: |                     mor_s_v(i)=mor_v(i-1)
  2230: *------             end do
  2231:                   endif
  2232:            
  2233:                 elseif (edge.eq.4) then
  2234:                   if(nn.eq.1)then
  2235: *------>            do i=2,lx1-1
  2236: |                     mor_s_v(i-1)=mor_v(i-1)
  2237: *------             end do
  2238:                     mor_s_v(4)=idmo(1,lx1,2,1,face,ntemp)
  2239:                   else
  2240:                     mor_s_v(1)=idmo(1,1,1,1,face,ntemp)
  2241: *------>            do i=2,lx1-1
  2242: |                     mor_s_v(i)=mor_v(i-1)
  2243: *------             end do
  2244:                   endif
  2245:           
  2246:                 elseif (edge.eq.1) then
  2247:                   if(nn.eq.1)then
  2248: *------>            do i=2,lx1-1
  2249: |                     mor_s_v(i-1)=mor_v(i-1)
  2250: *------             end do
  2251:                     mor_s_v(4)=idmo(lx1,1,1,2,face,ntemp)
  2252:                   else
  2253:                     mor_s_v(1)=idmo(1,1,1,1,face,ntemp)
  2254: *------>            do i=2,lx1-1
  2255: |                     mor_s_v(i)=mor_v(i-1)
  2256: *------             end do
  2257:                    endif
  2258:           
  2259:                 else if (edge.eq.2) then
  2260:                   if(nn.eq.1)then
  2261: *------>            do i=2,lx1-1
  2262: |                      mor_s_v(i-1)=mor_v(i-1)
  2263: *------             end do
  2264:                     mor_s_v(4)=idmo(lx1,lx1,2,2,face,ntemp)
  2265:                   else
  2266:                     mor_s_v(1)=idmo(lx1,1,1,2,face,ntemp)
  2267: *------>            do i=2,lx1-1
  2268: |                      mor_s_v(i)=mor_v(i-1)
  2269: *------             end do
  2270:                   endif
  2271:                 end if
  2272:           
  2273:           c.....copy mor_s_v to iel's local edge(op(edge)), on face jjface(face)
  2274:                 call mor_s_e_nn(op(edge),jjface(face),iel,mor_s_v,nn)
  2275:           c.....copy mor_s_v to iel's local edge(op(edge2)), on face jjface(face2)
  2276:           c     since this edge is shared by two faces on iel
  2277:                 call mor_s_e_nn(op(edge2),jjface(face2),iel,mor_s_v,nn)
  2278:           
  2279:                 return
  2280:                 end


