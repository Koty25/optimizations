NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:33 2020
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    20: par(1801): Parallel routine generated.: EXACT_RHS$1
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    26: vec( 101): Vectorized loop.
    28: vec( 101): Vectorized loop.
    29: opt(1418): Constant-length loop is expanded.
    40: par(1803): Parallelized by "do".
    40: opt(2495): Loop nest collapsed by "collapse" clause.
    45: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    45: vec( 103): Unvectorized loop.
    48: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    49: opt(1418): Constant-length loop is expanded.
    55: opt(1418): Constant-length loop is expanded.
    67: vec( 124): Iteration count is assumed. Iteration count.: 103
    67: vec( 101): Vectorized loop.
    71: vec( 128): Fused multiply-add operation applied.
    75: vec( 128): Fused multiply-add operation applied.
    81: vec( 128): Fused multiply-add operation applied.
    86: vec( 128): Fused multiply-add operation applied.
    91: vec( 128): Fused multiply-add operation applied.
   105: opt(1418): Constant-length loop is expanded.
   115: vec( 101): Vectorized loop.
   116: opt(1418): Constant-length loop is expanded.
   117: vec( 128): Fused multiply-add operation applied.
   123: opt(1418): Constant-length loop is expanded.
   141: par(1803): Parallelized by "do".
   141: opt(2495): Loop nest collapsed by "collapse" clause.
   146: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   146: vec( 103): Unvectorized loop.
   149: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   150: opt(1418): Constant-length loop is expanded.
   156: opt(1418): Constant-length loop is expanded.
   167: vec( 124): Iteration count is assumed. Iteration count.: 103
   167: vec( 101): Vectorized loop.
   171: vec( 128): Fused multiply-add operation applied.
   175: vec( 128): Fused multiply-add operation applied.
   180: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   191: vec( 128): Fused multiply-add operation applied.
   204: opt(1418): Constant-length loop is expanded.
   214: vec( 101): Vectorized loop.
   215: opt(1418): Constant-length loop is expanded.
   216: vec( 128): Fused multiply-add operation applied.
   222: opt(1418): Constant-length loop is expanded.
   241: par(1803): Parallelized by "do".
   241: opt(2495): Loop nest collapsed by "collapse" clause.
   246: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   246: vec( 103): Unvectorized loop.
   249: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
   250: opt(1418): Constant-length loop is expanded.
   256: opt(1418): Constant-length loop is expanded.
   267: vec( 124): Iteration count is assumed. Iteration count.: 103
   267: vec( 101): Vectorized loop.
   271: vec( 128): Fused multiply-add operation applied.
   275: vec( 128): Fused multiply-add operation applied.
   280: vec( 128): Fused multiply-add operation applied.
   285: vec( 128): Fused multiply-add operation applied.
   291: vec( 128): Fused multiply-add operation applied.
   304: opt(1418): Constant-length loop is expanded.
   314: vec( 101): Vectorized loop.
   315: opt(1418): Constant-length loop is expanded.
   316: vec( 128): Fused multiply-add operation applied.
   322: opt(1418): Constant-length loop is expanded.
   340: par(1803): Parallelized by "do".
   340: opt(2495): Loop nest collapsed by "collapse" clause.
   340: vec( 101): Vectorized loop.
   342: vec( 101): Vectorized loop.
   343: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:33 2020
FILE NAME: exact_rhs.f

PROCEDURE NAME: EXACT_RHS
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine exact_rhs
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c     compute the right hand side based on exact solution
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use bt_data
    15:                 implicit none
    16:           
    17:                 double precision dtemp(5), xi, eta, zeta, dtpp
    18:                 integer m, i, j, k, ip1, im1, jp1, jm1, km1, kp1
    19:           
    20:           !$omp parallel default(shared) private(i,j,k,m,zeta,eta,xi,
    21:           !$omp&  dtpp,im1,ip1,jm1,jp1,km1,kp1,dtemp)
    22:           c---------------------------------------------------------------------
    23:           c     initialize
    24:           c---------------------------------------------------------------------
    25:           !$omp do schedule(static) collapse(2)
    26: Y------>        do k= 0, grid_points(3)-1
    27: |+----->           do j = 0, grid_points(2)-1
    28: ||V---->              do i = 0, grid_points(1)-1
    29: |||*--->                 do m = 1, 5
    30: ||||                        forcing(m,i,j,k) = 0.0d0
    31: |||*---                  enddo
    32: ||V----               enddo
    33: |+-----            enddo
    34: Y------         enddo
    35:           
    36:           c---------------------------------------------------------------------
    37:           c     xi-direction flux differences
    38:           c---------------------------------------------------------------------
    39:           !$omp do schedule(static) collapse(2)
    40: P------>        do k = 1, grid_points(3)-2
    41: |+----->           do j = 1, grid_points(2)-2
    42: ||                    zeta = dble(k) * dnzm1
    43: ||                    eta = dble(j) * dnym1
    44: ||        
    45: ||+---->              do i=0, grid_points(1)-1
    46: |||                      xi = dble(i) * dnxm1
    47: |||       
    48: |||                      call exact_solution(xi, eta, zeta, dtemp)
    49: |||*--->                 do m = 1, 5
    50: ||||                        ue(i,m) = dtemp(m)
    51: |||*---                  enddo
    52: |||       
    53: |||                      dtpp = 1.0d0 / dtemp(1)
    54: |||       
    55: |||*--->                 do m = 2, 5
    56: ||||                        buf(i,m) = dtpp * dtemp(m)
    57: |||*---                  enddo
    58: |||       
    59: |||                      cuf(i)   = buf(i,2) * buf(i,2)
    60: |||                      buf(i,1) = cuf(i) + buf(i,3) * buf(i,3) +
    61: |||            >                 buf(i,4) * buf(i,4)
    62: |||                      q(i) = 0.5d0*(buf(i,2)*ue(i,2) + buf(i,3)*ue(i,3) +
    63: |||            >                 buf(i,4)*ue(i,4))
    64: |||       
    65: ||+----               enddo
    66: ||         
    67: ||V---->              do i = 1, grid_points(1)-2
    68: |||                      im1 = i-1
    69: |||                      ip1 = i+1
    70: |||       
    71: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
    72: |||            >                 tx2*( ue(ip1,2)-ue(im1,2) )+
    73: |||            >                 dx1tx1*(ue(ip1,1)-2.0d0*ue(i,1)+ue(im1,1))
    74: |||       
    75: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - tx2 * (
    76: |||            >                 (ue(ip1,2)*buf(ip1,2)+c2*(ue(ip1,5)-q(ip1)))-
    77: |||            >                 (ue(im1,2)*buf(im1,2)+c2*(ue(im1,5)-q(im1))))+
    78: |||            >                 xxcon1*(buf(ip1,2)-2.0d0*buf(i,2)+buf(im1,2))+
    79: |||            >                 dx2tx1*( ue(ip1,2)-2.0d0* ue(i,2)+ue(im1,2))
    80: |||       
    81: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - tx2 * (
    82: |||            >                 ue(ip1,3)*buf(ip1,2)-ue(im1,3)*buf(im1,2))+
    83: |||            >                 xxcon2*(buf(ip1,3)-2.0d0*buf(i,3)+buf(im1,3))+
    84: |||            >                 dx3tx1*( ue(ip1,3)-2.0d0*ue(i,3) +ue(im1,3))
    85: |||        
    86: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - tx2*(
    87: |||            >                 ue(ip1,4)*buf(ip1,2)-ue(im1,4)*buf(im1,2))+
    88: |||            >                 xxcon2*(buf(ip1,4)-2.0d0*buf(i,4)+buf(im1,4))+
    89: |||            >                 dx4tx1*( ue(ip1,4)-2.0d0* ue(i,4)+ ue(im1,4))
    90: |||       
    91: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - tx2*(
    92: |||            >                 buf(ip1,2)*(c1*ue(ip1,5)-c2*q(ip1))-
    93: |||            >                 buf(im1,2)*(c1*ue(im1,5)-c2*q(im1)))+
    94: |||            >                 0.5d0*xxcon3*(buf(ip1,1)-2.0d0*buf(i,1)+
    95: |||            >                 buf(im1,1))+
    96: |||            >                 xxcon4*(cuf(ip1)-2.0d0*cuf(i)+cuf(im1))+
    97: |||            >                 xxcon5*(buf(ip1,5)-2.0d0*buf(i,5)+buf(im1,5))+
    98: |||            >                 dx5tx1*( ue(ip1,5)-2.0d0* ue(i,5)+ ue(im1,5))
    99: ||V----               enddo
   100: ||        
   101: ||        c---------------------------------------------------------------------
   102: ||        c     Fourth-order dissipation
   103: ||        c---------------------------------------------------------------------
   104: ||        
   105: ||*---->              do m = 1, 5
   106: |||                      i = 1
   107: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   108: |||            >                    (5.0d0*ue(i,m) - 4.0d0*ue(i+1,m) +ue(i+2,m))
   109: |||                      i = 2
   110: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   111: |||            >                    (-4.0d0*ue(i-1,m) + 6.0d0*ue(i,m) -
   112: |||            >                    4.0d0*ue(i+1,m) +       ue(i+2,m))
   113: ||*----               enddo
   114: ||        
   115: ||V---->              do i = 3, grid_points(1)-4
   116: |||*--->                 do m = 1, 5
   117: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   118: ||||           >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   119: ||||           >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m) + ue(i+2,m))
   120: |||*---                  enddo
   121: ||V----               enddo
   122: ||        
   123: ||*---->              do m = 1, 5
   124: |||                      i = grid_points(1)-3
   125: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   126: |||            >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
   127: |||            >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m))
   128: |||                      i = grid_points(1)-2
   129: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   130: |||            >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) + 5.0d0*ue(i,m))
   131: ||*----               enddo
   132: ||        
   133: |+-----            enddo
   134: P------         enddo
   135:           !$omp end do
   136:           
   137:           c---------------------------------------------------------------------
   138:           c     eta-direction flux differences
   139:           c---------------------------------------------------------------------
   140:           !$omp do schedule(static) collapse(2)
   141: P------>        do k = 1, grid_points(3)-2
   142: |+----->           do i=1, grid_points(1)-2
   143: ||                    zeta = dble(k) * dnzm1
   144: ||                    xi = dble(i) * dnxm1
   145: ||        
   146: ||+---->              do j=0, grid_points(2)-1
   147: |||                      eta = dble(j) * dnym1
   148: |||       
   149: |||                      call exact_solution(xi, eta, zeta, dtemp)
   150: |||*--->                 do m = 1, 5
   151: ||||                        ue(j,m) = dtemp(m)
   152: |||*---                  enddo
   153: |||        
   154: |||                      dtpp = 1.0d0/dtemp(1)
   155: |||       
   156: |||*--->                 do m = 2, 5
   157: ||||                        buf(j,m) = dtpp * dtemp(m)
   158: |||*---                  enddo
   159: |||       
   160: |||                      cuf(j)   = buf(j,3) * buf(j,3)
   161: |||                      buf(j,1) = cuf(j) + buf(j,2) * buf(j,2) +
   162: |||            >                 buf(j,4) * buf(j,4)
   163: |||                      q(j) = 0.5d0*(buf(j,2)*ue(j,2) + buf(j,3)*ue(j,3) +
   164: |||            >                 buf(j,4)*ue(j,4))
   165: ||+----               enddo
   166: ||        
   167: ||V---->              do j = 1, grid_points(2)-2
   168: |||                      jm1 = j-1
   169: |||                      jp1 = j+1
   170: |||        
   171: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
   172: |||            >                 ty2*( ue(jp1,3)-ue(jm1,3) )+
   173: |||            >                 dy1ty1*(ue(jp1,1)-2.0d0*ue(j,1)+ue(jm1,1))
   174: |||       
   175: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - ty2*(
   176: |||            >                 ue(jp1,2)*buf(jp1,3)-ue(jm1,2)*buf(jm1,3))+
   177: |||            >                 yycon2*(buf(jp1,2)-2.0d0*buf(j,2)+buf(jm1,2))+
   178: |||            >                 dy2ty1*( ue(jp1,2)-2.0* ue(j,2)+ ue(jm1,2))
   179: |||       
   180: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - ty2*(
   181: |||            >                 (ue(jp1,3)*buf(jp1,3)+c2*(ue(jp1,5)-q(jp1)))-
   182: |||            >                 (ue(jm1,3)*buf(jm1,3)+c2*(ue(jm1,5)-q(jm1))))+
   183: |||            >                 yycon1*(buf(jp1,3)-2.0d0*buf(j,3)+buf(jm1,3))+
   184: |||            >                 dy3ty1*( ue(jp1,3)-2.0d0*ue(j,3) +ue(jm1,3))
   185: |||       
   186: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - ty2*(
   187: |||            >                 ue(jp1,4)*buf(jp1,3)-ue(jm1,4)*buf(jm1,3))+
   188: |||            >                 yycon2*(buf(jp1,4)-2.0d0*buf(j,4)+buf(jm1,4))+
   189: |||            >                 dy4ty1*( ue(jp1,4)-2.0d0*ue(j,4)+ ue(jm1,4))
   190: |||       
   191: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - ty2*(
   192: |||            >                 buf(jp1,3)*(c1*ue(jp1,5)-c2*q(jp1))-
   193: |||            >                 buf(jm1,3)*(c1*ue(jm1,5)-c2*q(jm1)))+
   194: |||            >                 0.5d0*yycon3*(buf(jp1,1)-2.0d0*buf(j,1)+
   195: |||            >                 buf(jm1,1))+
   196: |||            >                 yycon4*(cuf(jp1)-2.0d0*cuf(j)+cuf(jm1))+
   197: |||            >                 yycon5*(buf(jp1,5)-2.0d0*buf(j,5)+buf(jm1,5))+
   198: |||            >                 dy5ty1*(ue(jp1,5)-2.0d0*ue(j,5)+ue(jm1,5))
   199: ||V----               enddo
   200: ||        
   201: ||        c---------------------------------------------------------------------
   202: ||        c     Fourth-order dissipation
   203: ||        c---------------------------------------------------------------------
   204: ||*---->              do m = 1, 5
   205: |||                      j = 1
   206: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   207: |||            >                    (5.0d0*ue(j,m) - 4.0d0*ue(j+1,m) +ue(j+2,m))
   208: |||                      j = 2
   209: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   210: |||            >                    (-4.0d0*ue(j-1,m) + 6.0d0*ue(j,m) -
   211: |||            >                    4.0d0*ue(j+1,m) +       ue(j+2,m))
   212: ||*----               enddo
   213: ||        
   214: ||V---->              do j = 3, grid_points(2)-4
   215: |||*--->                 do m = 1, 5
   216: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   217: ||||           >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   218: ||||           >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m) + ue(j+2,m))
   219: |||*---                  enddo
   220: ||V----               enddo
   221: ||        
   222: ||*---->              do m = 1, 5
   223: |||                      j = grid_points(2)-3
   224: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   225: |||            >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
   226: |||            >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m))
   227: |||                      j = grid_points(2)-2
   228: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   229: |||            >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) + 5.0d0*ue(j,m))
   230: |||       
   231: ||*----               enddo
   232: ||        
   233: |+-----            enddo
   234: P------         enddo
   235:           !$omp end do
   236:           
   237:           c---------------------------------------------------------------------
   238:           c     zeta-direction flux differences
   239:           c---------------------------------------------------------------------
   240:           !$omp do schedule(static) collapse(2)
   241: P------>        do j=1, grid_points(2)-2
   242: |+----->           do i = 1, grid_points(1)-2
   243: ||                    eta = dble(j) * dnym1
   244: ||                    xi = dble(i) * dnxm1
   245: ||        
   246: ||+---->              do k=0, grid_points(3)-1
   247: |||                      zeta = dble(k) * dnzm1
   248: |||       
   249: |||                      call exact_solution(xi, eta, zeta, dtemp)
   250: |||*--->                 do m = 1, 5
   251: ||||                        ue(k,m) = dtemp(m)
   252: |||*---                  enddo
   253: |||       
   254: |||                      dtpp = 1.0d0/dtemp(1)
   255: |||       
   256: |||*--->                 do m = 2, 5
   257: ||||                        buf(k,m) = dtpp * dtemp(m)
   258: |||*---                  enddo
   259: |||       
   260: |||                      cuf(k)   = buf(k,4) * buf(k,4)
   261: |||                      buf(k,1) = cuf(k) + buf(k,2) * buf(k,2) +
   262: |||            >                 buf(k,3) * buf(k,3)
   263: |||                      q(k) = 0.5d0*(buf(k,2)*ue(k,2) + buf(k,3)*ue(k,3) +
   264: |||            >                 buf(k,4)*ue(k,4))
   265: ||+----               enddo
   266: ||        
   267: ||V---->              do k=1, grid_points(3)-2
   268: |||                      km1 = k-1
   269: |||                      kp1 = k+1
   270: |||        
   271: |||     F                forcing(1,i,j,k) = forcing(1,i,j,k) -
   272: |||            >                 tz2*( ue(kp1,4)-ue(km1,4) )+
   273: |||            >                 dz1tz1*(ue(kp1,1)-2.0d0*ue(k,1)+ue(km1,1))
   274: |||       
   275: |||     F                forcing(2,i,j,k) = forcing(2,i,j,k) - tz2 * (
   276: |||            >                 ue(kp1,2)*buf(kp1,4)-ue(km1,2)*buf(km1,4))+
   277: |||            >                 zzcon2*(buf(kp1,2)-2.0d0*buf(k,2)+buf(km1,2))+
   278: |||            >                 dz2tz1*( ue(kp1,2)-2.0d0* ue(k,2)+ ue(km1,2))
   279: |||       
   280: |||     F                forcing(3,i,j,k) = forcing(3,i,j,k) - tz2 * (
   281: |||            >                 ue(kp1,3)*buf(kp1,4)-ue(km1,3)*buf(km1,4))+
   282: |||            >                 zzcon2*(buf(kp1,3)-2.0d0*buf(k,3)+buf(km1,3))+
   283: |||            >                 dz3tz1*(ue(kp1,3)-2.0d0*ue(k,3)+ue(km1,3))
   284: |||       
   285: |||     F                forcing(4,i,j,k) = forcing(4,i,j,k) - tz2 * (
   286: |||            >                 (ue(kp1,4)*buf(kp1,4)+c2*(ue(kp1,5)-q(kp1)))-
   287: |||            >                 (ue(km1,4)*buf(km1,4)+c2*(ue(km1,5)-q(km1))))+
   288: |||            >                 zzcon1*(buf(kp1,4)-2.0d0*buf(k,4)+buf(km1,4))+
   289: |||            >                 dz4tz1*( ue(kp1,4)-2.0d0*ue(k,4) +ue(km1,4))
   290: |||       
   291: |||     F                forcing(5,i,j,k) = forcing(5,i,j,k) - tz2 * (
   292: |||            >                 buf(kp1,4)*(c1*ue(kp1,5)-c2*q(kp1))-
   293: |||            >                 buf(km1,4)*(c1*ue(km1,5)-c2*q(km1)))+
   294: |||            >                 0.5d0*zzcon3*(buf(kp1,1)-2.0d0*buf(k,1)
   295: |||            >                 +buf(km1,1))+
   296: |||            >                 zzcon4*(cuf(kp1)-2.0d0*cuf(k)+cuf(km1))+
   297: |||            >                 zzcon5*(buf(kp1,5)-2.0d0*buf(k,5)+buf(km1,5))+
   298: |||            >                 dz5tz1*( ue(kp1,5)-2.0d0*ue(k,5)+ ue(km1,5))
   299: ||V----               enddo
   300: ||        
   301: ||        c---------------------------------------------------------------------
   302: ||        c     Fourth-order dissipation
   303: ||        c---------------------------------------------------------------------
   304: ||*---->              do m = 1, 5
   305: |||                      k = 1
   306: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   307: |||            >                    (5.0d0*ue(k,m) - 4.0d0*ue(k+1,m) +ue(k+2,m))
   308: |||                      k = 2
   309: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   310: |||            >                    (-4.0d0*ue(k-1,m) + 6.0d0*ue(k,m) -
   311: |||            >                    4.0d0*ue(k+1,m) +       ue(k+2,m))
   312: ||*----               enddo
   313: ||        
   314: ||V---->              do k = 3, grid_points(3)-4
   315: |||*--->                 do m = 1, 5
   316: ||||    F                   forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
   317: ||||           >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   318: ||||           >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m) + ue(k+2,m))
   319: |||*---                  enddo
   320: ||V----               enddo
   321: ||        
   322: ||*---->              do m = 1, 5
   323: |||                      k = grid_points(3)-3
   324: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   325: |||            >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
   326: |||            >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m))
   327: |||                      k = grid_points(3)-2
   328: |||                      forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
   329: |||            >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) + 5.0d0*ue(k,m))
   330: ||*----               enddo
   331: ||        
   332: |+-----            enddo
   333: P------         enddo
   334:           !$omp end do
   335:           
   336:           c---------------------------------------------------------------------
   337:           c     now change the sign of the forcing function,
   338:           c---------------------------------------------------------------------
   339:           !$omp do schedule(static) collapse(2)
   340: Y------>        do k = 1, grid_points(3)-2
   341: |+----->           do j = 1, grid_points(2)-2
   342: ||V---->              do i = 1, grid_points(1)-2
   343: |||*--->                 do m = 1, 5
   344: ||||                        forcing(m,i,j,k) = -1.d0 * forcing(m,i,j,k)
   345: |||*---                  enddo
   346: ||V----               enddo
   347: |+-----            enddo
   348: Y------         enddo
   349:           !$omp end do nowait
   350:           !$omp end parallel
   351:           
   352:                 return
   353:                 end


