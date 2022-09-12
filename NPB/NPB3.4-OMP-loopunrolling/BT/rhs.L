NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:34 2020
FILE NAME: rhs.f

PROCEDURE NAME: COMPUTE_RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    19: par(1801): Parallel routine generated.: COMPUTE_RHS$1
    27: par(1803): Parallelized by "do".
    27: opt(2495): Loop nest collapsed by "collapse" clause.
    27: vec( 101): Vectorized loop.
    29: vec( 101): Vectorized loop.
    35: vec( 128): Fused multiply-add operation applied.
    52: par(1803): Parallelized by "do".
    52: opt(2495): Loop nest collapsed by "collapse" clause.
    52: vec( 101): Vectorized loop.
    54: vec( 101): Vectorized loop.
    55: opt(1418): Constant-length loop is expanded.
    64: par(1805): Master-construct generated.
    71: par(1803): Parallelized by "do".
    71: opt(2495): Loop nest collapsed by "collapse" clause.
    71: vec( 102): Partially vectorized loop.
    73: vec( 101): Vectorized loop.
    78: vec( 128): Fused multiply-add operation applied.
    83: vec( 128): Fused multiply-add operation applied.
    93: vec( 128): Fused multiply-add operation applied.
   101: vec( 128): Fused multiply-add operation applied.
   109: vec( 128): Fused multiply-add operation applied.
   129: opt(1418): Constant-length loop is expanded.
   136: opt(1418): Constant-length loop is expanded.
   142: vec( 101): Vectorized loop.
   143: opt(1418): Constant-length loop is expanded.
   144: vec( 128): Fused multiply-add operation applied.
   152: opt(1418): Constant-length loop is expanded.
   159: opt(1418): Constant-length loop is expanded.
   167: par(1805): Master-construct generated.
   176: par(1803): Parallelized by "do".
   176: opt(2495): Loop nest collapsed by "collapse" clause.
   176: vec( 102): Partially vectorized loop.
   178: vec( 101): Vectorized loop.
   182: vec( 128): Fused multiply-add operation applied.
   186: vec( 128): Fused multiply-add operation applied.
   193: vec( 128): Fused multiply-add operation applied.
   202: vec( 128): Fused multiply-add operation applied.
   209: vec( 128): Fused multiply-add operation applied.
   229: vec( 101): Vectorized loop.
   230: opt(1418): Constant-length loop is expanded.
   231: vec( 128): Fused multiply-add operation applied.
   238: vec( 101): Vectorized loop.
   239: opt(1418): Constant-length loop is expanded.
   240: vec( 128): Fused multiply-add operation applied.
   247: vec( 101): Vectorized loop.
   248: opt(1418): Constant-length loop is expanded.
   249: vec( 128): Fused multiply-add operation applied.
   256: vec( 101): Vectorized loop.
   257: opt(1418): Constant-length loop is expanded.
   258: vec( 128): Fused multiply-add operation applied.
   265: vec( 101): Vectorized loop.
   266: opt(1418): Constant-length loop is expanded.
   267: vec( 128): Fused multiply-add operation applied.
   277: par(1805): Master-construct generated.
   286: par(1803): Parallelized by "do".
   286: opt(2495): Loop nest collapsed by "collapse" clause.
   286: vec( 102): Partially vectorized loop.
   288: vec( 101): Vectorized loop.
   293: vec( 128): Fused multiply-add operation applied.
   297: vec( 128): Fused multiply-add operation applied.
   304: vec( 128): Fused multiply-add operation applied.
   311: vec( 128): Fused multiply-add operation applied.
   320: vec( 128): Fused multiply-add operation applied.
   340: vec( 101): Vectorized loop.
   341: opt(1418): Constant-length loop is expanded.
   342: vec( 128): Fused multiply-add operation applied.
   349: vec( 101): Vectorized loop.
   350: opt(1418): Constant-length loop is expanded.
   351: vec( 128): Fused multiply-add operation applied.
   358: vec( 101): Vectorized loop.
   359: opt(1418): Constant-length loop is expanded.
   360: vec( 128): Fused multiply-add operation applied.
   367: vec( 101): Vectorized loop.
   368: opt(1418): Constant-length loop is expanded.
   369: vec( 128): Fused multiply-add operation applied.
   376: vec( 101): Vectorized loop.
   377: opt(1418): Constant-length loop is expanded.
   378: vec( 128): Fused multiply-add operation applied.
   388: par(1805): Master-construct generated.
   393: par(1803): Parallelized by "do".
   393: opt(2495): Loop nest collapsed by "collapse" clause.
   393: vec( 101): Vectorized loop.
   395: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:34 2020
FILE NAME: rhs.f

PROCEDURE NAME: COMPUTE_RHS
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine compute_rhs
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:                 use bt_data
    10:                 implicit none
    11:           
    12:                 integer i, j, k, m
    13:                 double precision rho_inv, uijk, up1, um1, vijk, vp1, vm1,
    14:                >     wijk, wp1, wm1
    15:           
    16:           
    17:                 if (timeron) call timer_start(t_rhs)
    18:           
    19:           !$omp parallel default(shared) private(i,j,k,m,rho_inv,uijk,up1,um1,
    20:           !$omp&   vijk,vp1,vm1,wijk,wp1,wm1)
    21:           
    22:           c---------------------------------------------------------------------
    23:           c     compute the reciprocal of density, and the kinetic energy,
    24:           c     and the speed of sound.
    25:           c---------------------------------------------------------------------
    26:           !$omp do schedule(static) collapse(2)
    27: Y------>        do k = 0, grid_points(3)-1
    28: |+----->           do j = 0, grid_points(2)-1
    29: ||V---->              do i = 0, grid_points(1)-1
    30: |||                      rho_inv = 1.0d0/u(1,i,j,k)
    31: |||                      rho_i(i,j,k) = rho_inv
    32: |||                      us(i,j,k) = u(2,i,j,k) * rho_inv
    33: |||                      vs(i,j,k) = u(3,i,j,k) * rho_inv
    34: |||                      ws(i,j,k) = u(4,i,j,k) * rho_inv
    35: |||     F                square(i,j,k)     = 0.5d0* (
    36: |||            >                 u(2,i,j,k)*u(2,i,j,k) +
    37: |||            >                 u(3,i,j,k)*u(3,i,j,k) +
    38: |||            >                 u(4,i,j,k)*u(4,i,j,k) ) * rho_inv
    39: |||                      qs(i,j,k) = square(i,j,k) * rho_inv
    40: ||V----               enddo
    41: |+-----            enddo
    42: Y------         enddo
    43:           !$omp end do nowait
    44:           
    45:           c---------------------------------------------------------------------
    46:           c copy the exact forcing term to the right hand side;  because
    47:           c this forcing term is known, we can store it on the whole grid
    48:           c including the boundary
    49:           c---------------------------------------------------------------------
    50:           
    51:           !$omp do schedule(static) collapse(2)
    52: Y------>        do k = 0, grid_points(3)-1
    53: |+----->           do j = 0, grid_points(2)-1
    54: ||V---->              do i = 0, grid_points(1)-1
    55: |||*--->                 do m = 1, 5
    56: ||||                        rhs(m,i,j,k) = forcing(m,i,j,k)
    57: |||*---                  enddo
    58: ||V----               enddo
    59: |+-----            enddo
    60: Y------         enddo
    61:           !$omp end do
    62:           
    63:           
    64:           !$omp master
    65:                 if (timeron) call timer_start(t_rhsx)
    66:           !$omp end master
    67:           c---------------------------------------------------------------------
    68:           c     compute xi-direction fluxes
    69:           c---------------------------------------------------------------------
    70:           !$omp do schedule(static) collapse(2)
    71: Y------>        do k = 1, grid_points(3)-2
    72: |+----->           do j = 1, grid_points(2)-2
    73: ||V---->              do i = 1, grid_points(1)-2
    74: |||                      uijk = us(i,j,k)
    75: |||                      up1  = us(i+1,j,k)
    76: |||                      um1  = us(i-1,j,k)
    77: |||       
    78: |||     F                rhs(1,i,j,k) = rhs(1,i,j,k) + dx1tx1 *
    79: |||            >                 (u(1,i+1,j,k) - 2.0d0*u(1,i,j,k) +
    80: |||            >                 u(1,i-1,j,k)) -
    81: |||            >                 tx2 * (u(2,i+1,j,k) - u(2,i-1,j,k))
    82: |||       
    83: |||     F                rhs(2,i,j,k) = rhs(2,i,j,k) + dx2tx1 *
    84: |||            >                 (u(2,i+1,j,k) - 2.0d0*u(2,i,j,k) +
    85: |||            >                 u(2,i-1,j,k)) +
    86: |||            >                 xxcon2*con43 * (up1 - 2.0d0*uijk + um1) -
    87: |||            >                 tx2 * (u(2,i+1,j,k)*up1 -
    88: |||            >                 u(2,i-1,j,k)*um1 +
    89: |||            >                 (u(5,i+1,j,k)- square(i+1,j,k)-
    90: |||            >                 u(5,i-1,j,k)+ square(i-1,j,k))*
    91: |||            >                 c2)
    92: |||       
    93: |||     F                rhs(3,i,j,k) = rhs(3,i,j,k) + dx3tx1 *
    94: |||            >                 (u(3,i+1,j,k) - 2.0d0*u(3,i,j,k) +
    95: |||            >                 u(3,i-1,j,k)) +
    96: |||            >                 xxcon2 * (vs(i+1,j,k) - 2.0d0*vs(i,j,k) +
    97: |||            >                 vs(i-1,j,k)) -
    98: |||            >                 tx2 * (u(3,i+1,j,k)*up1 -
    99: |||            >                 u(3,i-1,j,k)*um1)
   100: |||       
   101: |||     F                rhs(4,i,j,k) = rhs(4,i,j,k) + dx4tx1 *
   102: |||            >                 (u(4,i+1,j,k) - 2.0d0*u(4,i,j,k) +
   103: |||            >                 u(4,i-1,j,k)) +
   104: |||            >                 xxcon2 * (ws(i+1,j,k) - 2.0d0*ws(i,j,k) +
   105: |||            >                 ws(i-1,j,k)) -
   106: |||            >                 tx2 * (u(4,i+1,j,k)*up1 -
   107: |||            >                 u(4,i-1,j,k)*um1)
   108: |||       
   109: |||     F                rhs(5,i,j,k) = rhs(5,i,j,k) + dx5tx1 *
   110: |||            >                 (u(5,i+1,j,k) - 2.0d0*u(5,i,j,k) +
   111: |||            >                 u(5,i-1,j,k)) +
   112: |||            >                 xxcon3 * (qs(i+1,j,k) - 2.0d0*qs(i,j,k) +
   113: |||            >                 qs(i-1,j,k)) +
   114: |||            >                 xxcon4 * (up1*up1 -       2.0d0*uijk*uijk +
   115: |||            >                 um1*um1) +
   116: |||            >                 xxcon5 * (u(5,i+1,j,k)*rho_i(i+1,j,k) -
   117: |||            >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   118: |||            >                 u(5,i-1,j,k)*rho_i(i-1,j,k)) -
   119: |||            >                 tx2 * ( (c1*u(5,i+1,j,k) -
   120: |||            >                 c2*square(i+1,j,k))*up1 -
   121: |||            >                 (c1*u(5,i-1,j,k) -
   122: |||            >                 c2*square(i-1,j,k))*um1 )
   123: ||V----               enddo
   124: ||        
   125: ||        c---------------------------------------------------------------------
   126: ||        c     add fourth order xi-direction dissipation
   127: ||        c---------------------------------------------------------------------
   128: ||                    i = 1
   129: ||*---->              do m = 1, 5
   130: |||                      rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   131: |||            >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
   132: |||            >                    u(m,i+2,j,k))
   133: ||*----               enddo
   134: ||        
   135: ||                    i = 2
   136: ||*---->              do m = 1, 5
   137: |||                      rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   138: |||            >                    (-4.0d0*u(m,i-1,j,k) + 6.0d0*u(m,i,j,k) -
   139: |||            >                    4.0d0*u(m,i+1,j,k) + u(m,i+2,j,k))
   140: ||*----               enddo
   141: ||        
   142: ||V---->              do i = 3,grid_points(1)-4
   143: |||*--->                 do m = 1, 5
   144: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   145: ||||           >                    (  u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) +
   146: ||||           >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
   147: ||||           >                    u(m,i+2,j,k) )
   148: |||*---                  enddo
   149: ||V----               enddo
   150: ||         
   151: ||                    i = grid_points(1)-3
   152: ||*---->              do m = 1, 5
   153: |||                      rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   154: |||            >                    ( u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) +
   155: |||            >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) )
   156: ||*----               enddo
   157: ||        
   158: ||                    i = grid_points(1)-2
   159: ||*---->              do m = 1, 5
   160: |||                      rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   161: |||            >                    ( u(m,i-2,j,k) - 4.d0*u(m,i-1,j,k) +
   162: |||            >                    5.d0*u(m,i,j,k) )
   163: ||*----               enddo
   164: |+-----            enddo
   165: Y------         enddo
   166:           !$omp end do nowait
   167:           !$omp master
   168:                 if (timeron) call timer_stop(t_rhsx)
   169:           
   170:                 if (timeron) call timer_start(t_rhsy)
   171:           !$omp end master
   172:           c---------------------------------------------------------------------
   173:           c     compute eta-direction fluxes
   174:           c---------------------------------------------------------------------
   175:           !$omp do schedule(static) collapse(2)
   176: Y------>        do k = 1, grid_points(3)-2
   177: |+----->           do j = 1, grid_points(2)-2
   178: ||V---->              do i = 1, grid_points(1)-2
   179: |||                      vijk = vs(i,j,k)
   180: |||                      vp1  = vs(i,j+1,k)
   181: |||                      vm1  = vs(i,j-1,k)
   182: |||     F                rhs(1,i,j,k) = rhs(1,i,j,k) + dy1ty1 *
   183: |||            >                 (u(1,i,j+1,k) - 2.0d0*u(1,i,j,k) +
   184: |||            >                 u(1,i,j-1,k)) -
   185: |||            >                 ty2 * (u(3,i,j+1,k) - u(3,i,j-1,k))
   186: |||     F                rhs(2,i,j,k) = rhs(2,i,j,k) + dy2ty1 *
   187: |||            >                 (u(2,i,j+1,k) - 2.0d0*u(2,i,j,k) +
   188: |||            >                 u(2,i,j-1,k)) +
   189: |||            >                 yycon2 * (us(i,j+1,k) - 2.0d0*us(i,j,k) +
   190: |||            >                 us(i,j-1,k)) -
   191: |||            >                 ty2 * (u(2,i,j+1,k)*vp1 -
   192: |||            >                 u(2,i,j-1,k)*vm1)
   193: |||     F                rhs(3,i,j,k) = rhs(3,i,j,k) + dy3ty1 *
   194: |||            >                 (u(3,i,j+1,k) - 2.0d0*u(3,i,j,k) +
   195: |||            >                 u(3,i,j-1,k)) +
   196: |||            >                 yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
   197: |||            >                 ty2 * (u(3,i,j+1,k)*vp1 -
   198: |||            >                 u(3,i,j-1,k)*vm1 +
   199: |||            >                 (u(5,i,j+1,k) - square(i,j+1,k) -
   200: |||            >                 u(5,i,j-1,k) + square(i,j-1,k))
   201: |||            >                 *c2)
   202: |||     F                rhs(4,i,j,k) = rhs(4,i,j,k) + dy4ty1 *
   203: |||            >                 (u(4,i,j+1,k) - 2.0d0*u(4,i,j,k) +
   204: |||            >                 u(4,i,j-1,k)) +
   205: |||            >                 yycon2 * (ws(i,j+1,k) - 2.0d0*ws(i,j,k) +
   206: |||            >                 ws(i,j-1,k)) -
   207: |||            >                 ty2 * (u(4,i,j+1,k)*vp1 -
   208: |||            >                 u(4,i,j-1,k)*vm1)
   209: |||     F                rhs(5,i,j,k) = rhs(5,i,j,k) + dy5ty1 *
   210: |||            >                 (u(5,i,j+1,k) - 2.0d0*u(5,i,j,k) +
   211: |||            >                 u(5,i,j-1,k)) +
   212: |||            >                 yycon3 * (qs(i,j+1,k) - 2.0d0*qs(i,j,k) +
   213: |||            >                 qs(i,j-1,k)) +
   214: |||            >                 yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk +
   215: |||            >                 vm1*vm1) +
   216: |||            >                 yycon5 * (u(5,i,j+1,k)*rho_i(i,j+1,k) -
   217: |||            >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   218: |||            >                 u(5,i,j-1,k)*rho_i(i,j-1,k)) -
   219: |||            >                 ty2 * ((c1*u(5,i,j+1,k) -
   220: |||            >                 c2*square(i,j+1,k)) * vp1 -
   221: |||            >                 (c1*u(5,i,j-1,k) -
   222: |||            >                 c2*square(i,j-1,k)) * vm1)
   223: ||V----               enddo
   224: ||        
   225: ||        c---------------------------------------------------------------------
   226: ||        c     add fourth order eta-direction dissipation
   227: ||        c---------------------------------------------------------------------
   228: ||                    if (j .eq. 1) then
   229: ||V---->                 do i = 1, grid_points(1)-2
   230: |||*--->                 do m = 1, 5
   231: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   232: ||||           >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
   233: ||||           >                    u(m,i,j+2,k))
   234: |||*---                  enddo
   235: ||V----                  enddo
   236: ||        
   237: ||                    else if (j .eq. 2) then
   238: ||V---->                 do i = 1, grid_points(1)-2
   239: |||*--->                 do m = 1, 5
   240: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   241: ||||           >                    (-4.0d0*u(m,i,j-1,k) + 6.0d0*u(m,i,j,k) -
   242: ||||           >                    4.0d0*u(m,i,j+1,k) + u(m,i,j+2,k))
   243: |||*---                  enddo
   244: ||V----                  enddo
   245: ||         
   246: ||                    else if (j .eq. grid_points(2)-3) then
   247: ||V---->                 do i = 1, grid_points(1)-2
   248: |||*--->                 do m = 1, 5
   249: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   250: ||||           >                    ( u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) +
   251: ||||           >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) )
   252: |||*---                  enddo
   253: ||V----                  enddo
   254: ||        
   255: ||                    else if (j .eq. grid_points(2)-2) then
   256: ||V---->                 do i = 1, grid_points(1)-2
   257: |||*--->                 do m = 1, 5
   258: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   259: ||||           >                    ( u(m,i,j-2,k) - 4.d0*u(m,i,j-1,k) +
   260: ||||           >                    5.d0*u(m,i,j,k) )
   261: |||*---                  enddo
   262: ||V----                  enddo
   263: ||        
   264: ||                    else
   265: ||V---->                 do i = 1,grid_points(1)-2
   266: |||*--->                 do m = 1, 5
   267: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   268: ||||           >                    (  u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) +
   269: ||||           >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
   270: ||||           >                    u(m,i,j+2,k) )
   271: |||*---                  enddo
   272: ||V----                  enddo
   273: ||                    endif
   274: |+-----            enddo
   275: Y------         enddo
   276:           !$omp end do nowait
   277:           !$omp master
   278:                 if (timeron) call timer_stop(t_rhsy)
   279:           
   280:                 if (timeron) call timer_start(t_rhsz)
   281:           !$omp end master
   282:           c---------------------------------------------------------------------
   283:           c     compute zeta-direction fluxes
   284:           c---------------------------------------------------------------------
   285:           !$omp do schedule(static) collapse(2)
   286: Y------>        do k = 1, grid_points(3)-2
   287: |+----->           do j = 1, grid_points(2)-2
   288: ||V---->              do i = 1, grid_points(1)-2
   289: |||                      wijk = ws(i,j,k)
   290: |||                      wp1  = ws(i,j,k+1)
   291: |||                      wm1  = ws(i,j,k-1)
   292: |||       
   293: |||     F                rhs(1,i,j,k) = rhs(1,i,j,k) + dz1tz1 *
   294: |||            >                 (u(1,i,j,k+1) - 2.0d0*u(1,i,j,k) +
   295: |||            >                 u(1,i,j,k-1)) -
   296: |||            >                 tz2 * (u(4,i,j,k+1) - u(4,i,j,k-1))
   297: |||     F                rhs(2,i,j,k) = rhs(2,i,j,k) + dz2tz1 *
   298: |||            >                 (u(2,i,j,k+1) - 2.0d0*u(2,i,j,k) +
   299: |||            >                 u(2,i,j,k-1)) +
   300: |||            >                 zzcon2 * (us(i,j,k+1) - 2.0d0*us(i,j,k) +
   301: |||            >                 us(i,j,k-1)) -
   302: |||            >                 tz2 * (u(2,i,j,k+1)*wp1 -
   303: |||            >                 u(2,i,j,k-1)*wm1)
   304: |||     F                rhs(3,i,j,k) = rhs(3,i,j,k) + dz3tz1 *
   305: |||            >                 (u(3,i,j,k+1) - 2.0d0*u(3,i,j,k) +
   306: |||            >                 u(3,i,j,k-1)) +
   307: |||            >                 zzcon2 * (vs(i,j,k+1) - 2.0d0*vs(i,j,k) +
   308: |||            >                 vs(i,j,k-1)) -
   309: |||            >                 tz2 * (u(3,i,j,k+1)*wp1 -
   310: |||            >                 u(3,i,j,k-1)*wm1)
   311: |||     F                rhs(4,i,j,k) = rhs(4,i,j,k) + dz4tz1 *
   312: |||            >                 (u(4,i,j,k+1) - 2.0d0*u(4,i,j,k) +
   313: |||            >                 u(4,i,j,k-1)) +
   314: |||            >                 zzcon2*con43 * (wp1 - 2.0d0*wijk + wm1) -
   315: |||            >                 tz2 * (u(4,i,j,k+1)*wp1 -
   316: |||            >                 u(4,i,j,k-1)*wm1 +
   317: |||            >                 (u(5,i,j,k+1) - square(i,j,k+1) -
   318: |||            >                 u(5,i,j,k-1) + square(i,j,k-1))
   319: |||            >                 *c2)
   320: |||     F                rhs(5,i,j,k) = rhs(5,i,j,k) + dz5tz1 *
   321: |||            >                 (u(5,i,j,k+1) - 2.0d0*u(5,i,j,k) +
   322: |||            >                 u(5,i,j,k-1)) +
   323: |||            >                 zzcon3 * (qs(i,j,k+1) - 2.0d0*qs(i,j,k) +
   324: |||            >                 qs(i,j,k-1)) +
   325: |||            >                 zzcon4 * (wp1*wp1 - 2.0d0*wijk*wijk +
   326: |||            >                 wm1*wm1) +
   327: |||            >                 zzcon5 * (u(5,i,j,k+1)*rho_i(i,j,k+1) -
   328: |||            >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   329: |||            >                 u(5,i,j,k-1)*rho_i(i,j,k-1)) -
   330: |||            >                 tz2 * ( (c1*u(5,i,j,k+1) -
   331: |||            >                 c2*square(i,j,k+1))*wp1 -
   332: |||            >                 (c1*u(5,i,j,k-1) -
   333: |||            >                 c2*square(i,j,k-1))*wm1)
   334: ||V----               enddo
   335: ||        
   336: ||        c---------------------------------------------------------------------
   337: ||        c     add fourth order zeta-direction dissipation
   338: ||        c---------------------------------------------------------------------
   339: ||                    if (k.eq.1) then
   340: ||V---->                 do i = 1, grid_points(1)-2
   341: |||*--->                 do m = 1, 5
   342: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   343: ||||           >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
   344: ||||           >                    u(m,i,j,k+2))
   345: |||*---                  enddo
   346: ||V----                  enddo
   347: ||        
   348: ||                    else if (k.eq.2) then
   349: ||V---->                 do i = 1, grid_points(1)-2
   350: |||*--->                 do m = 1, 5
   351: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   352: ||||           >                    (-4.0d0*u(m,i,j,k-1) + 6.0d0*u(m,i,j,k) -
   353: ||||           >                    4.0d0*u(m,i,j,k+1) + u(m,i,j,k+2))
   354: |||*---                  enddo
   355: ||V----                  enddo
   356: ||        
   357: ||                    else if (k.eq.grid_points(3)-3) then
   358: ||V---->                 do i = 1, grid_points(1)-2
   359: |||*--->                 do m = 1, 5
   360: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   361: ||||           >                    ( u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) +
   362: ||||           >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) )
   363: |||*---                  enddo
   364: ||V----                  enddo
   365: ||        
   366: ||                    else if (k.eq.grid_points(3)-2) then
   367: ||V---->                 do i = 1, grid_points(1)-2
   368: |||*--->                 do m = 1, 5
   369: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   370: ||||           >                    ( u(m,i,j,k-2) - 4.d0*u(m,i,j,k-1) +
   371: ||||           >                    5.d0*u(m,i,j,k) )
   372: |||*---                  enddo
   373: ||V----                  enddo
   374: ||        
   375: ||                    else
   376: ||V---->                 do i = 1,grid_points(1)-2
   377: |||*--->                 do m = 1, 5
   378: ||||    F                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   379: ||||           >                    (  u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) +
   380: ||||           >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
   381: ||||           >                    u(m,i,j,k+2) )
   382: |||*---                  enddo
   383: ||V----                  enddo
   384: ||                    endif
   385: |+-----            enddo
   386: Y------         enddo
   387:           !$omp end do nowait
   388:           !$omp master
   389:                 if (timeron) call timer_stop(t_rhsz)
   390:           !$omp end master
   391:           
   392:           !$omp do schedule(static) collapse(2)
   393: Y------>        do k = 1, grid_points(3)-2
   394: |+----->           do j = 1, grid_points(2)-2
   395: ||V---->              do i = 1, grid_points(1)-2
   396: |||                      rhs(1,i,j,k) = rhs(1,i,j,k) * dt
   397: |||                      rhs(2,i,j,k) = rhs(2,i,j,k) * dt
   398: |||                      rhs(3,i,j,k) = rhs(3,i,j,k) * dt
   399: |||                      rhs(4,i,j,k) = rhs(4,i,j,k) * dt
   400: |||                      rhs(5,i,j,k) = rhs(5,i,j,k) * dt
   401: ||V----               enddo
   402: |+-----            enddo
   403: Y------         enddo
   404:           !$omp end do nowait
   405:           
   406:           !$omp end parallel
   407:           
   408:                 if (timeron) call timer_stop(t_rhs)
   409:           
   410:                 return
   411:                 end


