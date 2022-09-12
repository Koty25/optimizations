NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: rhs.f

PROCEDURE NAME: COMPUTE_RHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    17: inl(1212): Source for routine not found.: TIMER_START
    18: par(1801): Parallel routine generated.: COMPUTE_RHS$1
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    26: vec( 101): Vectorized loop.
    28: vec( 101): Vectorized loop.
    34: vec( 128): Fused multiply-add operation applied.
    56: par(1803): Parallelized by "do".
    56: opt(2495): Loop nest collapsed by "collapse" clause.
    56: vec( 101): Vectorized loop.
    58: vec( 101): Vectorized loop.
    59: opt(1418): Constant-length loop is expanded.
    70: par(1805): Master-construct generated.
    74: par(1803): Parallelized by "do".
    74: opt(2495): Loop nest collapsed by "collapse" clause.
    74: vec( 102): Partially vectorized loop.
    76: vec( 101): Vectorized loop.
    81: vec( 128): Fused multiply-add operation applied.
    86: vec( 128): Fused multiply-add operation applied.
    96: vec( 128): Fused multiply-add operation applied.
   104: vec( 128): Fused multiply-add operation applied.
   112: vec( 128): Fused multiply-add operation applied.
   132: opt(1418): Constant-length loop is expanded.
   139: opt(1418): Constant-length loop is expanded.
   145: vec( 101): Vectorized loop.
   146: opt(1418): Constant-length loop is expanded.
   147: vec( 128): Fused multiply-add operation applied.
   155: opt(1418): Constant-length loop is expanded.
   162: opt(1418): Constant-length loop is expanded.
   170: par(1805): Master-construct generated.
   171: inl(1212): Source for routine not found.: TIMER_STOP
   179: par(1803): Parallelized by "do".
   179: opt(2495): Loop nest collapsed by "collapse" clause.
   179: vec( 102): Partially vectorized loop.
   181: vec( 101): Vectorized loop.
   185: vec( 128): Fused multiply-add operation applied.
   189: vec( 128): Fused multiply-add operation applied.
   196: vec( 128): Fused multiply-add operation applied.
   205: vec( 128): Fused multiply-add operation applied.
   212: vec( 128): Fused multiply-add operation applied.
   234: vec( 101): Vectorized loop.
   235: opt(1418): Constant-length loop is expanded.
   236: vec( 128): Fused multiply-add operation applied.
   243: vec( 101): Vectorized loop.
   244: opt(1418): Constant-length loop is expanded.
   245: vec( 128): Fused multiply-add operation applied.
   252: vec( 101): Vectorized loop.
   253: opt(1418): Constant-length loop is expanded.
   254: vec( 128): Fused multiply-add operation applied.
   261: vec( 101): Vectorized loop.
   262: opt(1418): Constant-length loop is expanded.
   263: vec( 128): Fused multiply-add operation applied.
   270: vec( 101): Vectorized loop.
   271: opt(1418): Constant-length loop is expanded.
   272: vec( 128): Fused multiply-add operation applied.
   282: par(1805): Master-construct generated.
   291: par(1803): Parallelized by "do".
   291: opt(2495): Loop nest collapsed by "collapse" clause.
   291: vec( 102): Partially vectorized loop.
   293: vec( 101): Vectorized loop.
   298: vec( 128): Fused multiply-add operation applied.
   302: vec( 128): Fused multiply-add operation applied.
   309: vec( 128): Fused multiply-add operation applied.
   316: vec( 128): Fused multiply-add operation applied.
   325: vec( 128): Fused multiply-add operation applied.
   346: vec( 101): Vectorized loop.
   347: opt(1418): Constant-length loop is expanded.
   348: vec( 128): Fused multiply-add operation applied.
   355: vec( 101): Vectorized loop.
   356: opt(1418): Constant-length loop is expanded.
   357: vec( 128): Fused multiply-add operation applied.
   364: vec( 101): Vectorized loop.
   365: opt(1418): Constant-length loop is expanded.
   366: vec( 128): Fused multiply-add operation applied.
   373: vec( 101): Vectorized loop.
   374: opt(1418): Constant-length loop is expanded.
   375: vec( 128): Fused multiply-add operation applied.
   382: vec( 101): Vectorized loop.
   383: opt(1418): Constant-length loop is expanded.
   384: vec( 128): Fused multiply-add operation applied.
   394: par(1805): Master-construct generated.
   399: par(1803): Parallelized by "do".
   399: opt(2495): Loop nest collapsed by "collapse" clause.
   399: vec( 101): Vectorized loop.
   401: vec( 101): Vectorized loop.
   402: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: rhs.f

PROCEDURE NAME: COMPUTE_RHS
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                  subroutine compute_rhs
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:                  use sp_data
    10:                  implicit none
    11:           
    12:                  integer i, j, k, m
    13:                  double precision aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,
    14:                >                  wijk, wp1, wm1
    15:           
    16:           
    17:                  if (timeron) call timer_start(t_rhs)
    18:           !$omp parallel default(shared) private(i,j,k,m,rho_inv,aux,uijk,up1,um1,
    19:           !$omp&   vijk,vp1,vm1,wijk,wp1,wm1)
    20:           c---------------------------------------------------------------------
    21:           c      compute the reciprocal of density, and the kinetic energy,
    22:           c      and the speed of sound.
    23:           c---------------------------------------------------------------------
    24:           
    25:           !$omp do schedule(static) collapse(2)
    26: Y------>         do    k = 0, grid_points(3)-1
    27: |+----->            do    j = 0, grid_points(2)-1
    28: ||V---->               do    i = 0, grid_points(1)-1
    29: |||                       rho_inv = 1.0d0/u(1,i,j,k)
    30: |||                       rho_i(i,j,k) = rho_inv
    31: |||                       us(i,j,k) = u(2,i,j,k) * rho_inv
    32: |||                       vs(i,j,k) = u(3,i,j,k) * rho_inv
    33: |||                       ws(i,j,k) = u(4,i,j,k) * rho_inv
    34: |||     F                 square(i,j,k)     = 0.5d0* (
    35: |||            >                        u(2,i,j,k)*u(2,i,j,k) +
    36: |||            >                        u(3,i,j,k)*u(3,i,j,k) +
    37: |||            >                        u(4,i,j,k)*u(4,i,j,k) ) * rho_inv
    38: |||                       qs(i,j,k) = square(i,j,k) * rho_inv
    39: |||       c---------------------------------------------------------------------
    40: |||       c               (don't need speed and ainx until the lhs computation)
    41: |||       c---------------------------------------------------------------------
    42: |||                       aux = c1c2*rho_inv* (u(5,i,j,k) - square(i,j,k))
    43: |||                       speed(i,j,k) = dsqrt(aux)
    44: ||V----                end do
    45: |+-----             end do
    46: Y------          end do
    47:           !$omp end do nowait
    48:           
    49:           c---------------------------------------------------------------------
    50:           c copy the exact forcing term to the right hand side;  because
    51:           c this forcing term is known, we can store it on the whole grid
    52:           c including the boundary
    53:           c---------------------------------------------------------------------
    54:           
    55:           !$omp do schedule(static) collapse(2)
    56: Y------>         do    k = 0, nz2+1
    57: |+----->            do    j = 0, ny2+1
    58: ||V---->               do    i = 0, nx2+1
    59: |||*--->                  do    m = 1, 5
    60: ||||                         rhs(m,i,j,k) = forcing(m,i,j,k)
    61: |||*---                   end do
    62: ||V----                end do
    63: |+-----             end do
    64: Y------          end do
    65:           !$omp end do
    66:           
    67:           c---------------------------------------------------------------------
    68:           c      compute xi-direction fluxes
    69:           c---------------------------------------------------------------------
    70:           !$omp master
    71:                  if (timeron) call timer_start(t_rhsx)
    72:           !$omp end master
    73:           !$omp do schedule(static) collapse(2)
    74: Y------>         do    k = 1, nz2
    75: |+----->            do    j = 1, ny2
    76: ||V---->               do    i = 1, nx2
    77: |||                       uijk = us(i,j,k)
    78: |||                       up1  = us(i+1,j,k)
    79: |||                       um1  = us(i-1,j,k)
    80: |||       
    81: |||     F                 rhs(1,i,j,k) = rhs(1,i,j,k) + dx1tx1 *
    82: |||            >                    (u(1,i+1,j,k) - 2.0d0*u(1,i,j,k) +
    83: |||            >                     u(1,i-1,j,k)) -
    84: |||            >                    tx2 * (u(2,i+1,j,k) - u(2,i-1,j,k))
    85: |||       
    86: |||     F                 rhs(2,i,j,k) = rhs(2,i,j,k) + dx2tx1 *
    87: |||            >                    (u(2,i+1,j,k) - 2.0d0*u(2,i,j,k) +
    88: |||            >                     u(2,i-1,j,k)) +
    89: |||            >                    xxcon2*con43 * (up1 - 2.0d0*uijk + um1) -
    90: |||            >                    tx2 * (u(2,i+1,j,k)*up1 -
    91: |||            >                           u(2,i-1,j,k)*um1 +
    92: |||            >                           (u(5,i+1,j,k)- square(i+1,j,k)-
    93: |||            >                            u(5,i-1,j,k)+ square(i-1,j,k))*
    94: |||            >                            c2)
    95: |||       
    96: |||     F                 rhs(3,i,j,k) = rhs(3,i,j,k) + dx3tx1 *
    97: |||            >                    (u(3,i+1,j,k) - 2.0d0*u(3,i,j,k) +
    98: |||            >                     u(3,i-1,j,k)) +
    99: |||            >                    xxcon2 * (vs(i+1,j,k) - 2.0d0*vs(i,j,k) +
   100: |||            >                              vs(i-1,j,k)) -
   101: |||            >                    tx2 * (u(3,i+1,j,k)*up1 -
   102: |||            >                           u(3,i-1,j,k)*um1)
   103: |||       
   104: |||     F                 rhs(4,i,j,k) = rhs(4,i,j,k) + dx4tx1 *
   105: |||            >                    (u(4,i+1,j,k) - 2.0d0*u(4,i,j,k) +
   106: |||            >                     u(4,i-1,j,k)) +
   107: |||            >                    xxcon2 * (ws(i+1,j,k) - 2.0d0*ws(i,j,k) +
   108: |||            >                              ws(i-1,j,k)) -
   109: |||            >                    tx2 * (u(4,i+1,j,k)*up1 -
   110: |||            >                           u(4,i-1,j,k)*um1)
   111: |||       
   112: |||     F                 rhs(5,i,j,k) = rhs(5,i,j,k) + dx5tx1 *
   113: |||            >                    (u(5,i+1,j,k) - 2.0d0*u(5,i,j,k) +
   114: |||            >                     u(5,i-1,j,k)) +
   115: |||            >                    xxcon3 * (qs(i+1,j,k) - 2.0d0*qs(i,j,k) +
   116: |||            >                              qs(i-1,j,k)) +
   117: |||            >                    xxcon4 * (up1*up1 -       2.0d0*uijk*uijk +
   118: |||            >                              um1*um1) +
   119: |||            >                    xxcon5 * (u(5,i+1,j,k)*rho_i(i+1,j,k) -
   120: |||            >                              2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   121: |||            >                              u(5,i-1,j,k)*rho_i(i-1,j,k)) -
   122: |||            >                    tx2 * ( (c1*u(5,i+1,j,k) -
   123: |||            >                             c2*square(i+1,j,k))*up1 -
   124: |||            >                            (c1*u(5,i-1,j,k) -
   125: |||            >                             c2*square(i-1,j,k))*um1 )
   126: ||V----                end do
   127: ||        
   128: ||        c---------------------------------------------------------------------
   129: ||        c      add fourth order xi-direction dissipation
   130: ||        c---------------------------------------------------------------------
   131: ||                     i = 1
   132: ||*---->               do    m = 1, 5
   133: |||                       rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   134: |||            >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
   135: |||            >                            u(m,i+2,j,k))
   136: ||*----                end do
   137: ||        
   138: ||                     i = 2
   139: ||*---->               do    m = 1, 5
   140: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   141: |||            >                    (-4.0d0*u(m,i-1,j,k) + 6.0d0*u(m,i,j,k) -
   142: |||            >                      4.0d0*u(m,i+1,j,k) + u(m,i+2,j,k))
   143: ||*----                end do
   144: ||        
   145: ||V---->               do  i = 3, nx2-2
   146: |||*--->                  do     m = 1, 5
   147: ||||    F                    rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   148: ||||           >                    (  u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) +
   149: ||||           >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
   150: ||||           >                         u(m,i+2,j,k) )
   151: |||*---                   end do
   152: ||V----                end do
   153: ||        
   154: ||                     i = nx2-1
   155: ||*---->               do     m = 1, 5
   156: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   157: |||            >                    ( u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) +
   158: |||            >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) )
   159: ||*----                end do
   160: ||        
   161: ||                     i = nx2
   162: ||*---->               do     m = 1, 5
   163: |||                       rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   164: |||            >                    ( u(m,i-2,j,k) - 4.d0*u(m,i-1,j,k) +
   165: |||            >                      5.d0*u(m,i,j,k) )
   166: ||*----                end do
   167: |+-----             end do
   168: Y------          end do
   169:           !$omp end do nowait
   170:           !$omp master
   171:                  if (timeron) call timer_stop(t_rhsx)
   172:           
   173:           c---------------------------------------------------------------------
   174:           c      compute eta-direction fluxes
   175:           c---------------------------------------------------------------------
   176:                  if (timeron) call timer_start(t_rhsy)
   177:           !$omp end master
   178:           !$omp do schedule(static) collapse(2)
   179: Y------>         do     k = 1, nz2
   180: |+----->            do     j = 1, ny2
   181: ||V---->               do     i = 1, nx2
   182: |||                       vijk = vs(i,j,k)
   183: |||                       vp1  = vs(i,j+1,k)
   184: |||                       vm1  = vs(i,j-1,k)
   185: |||     F                 rhs(1,i,j,k) = rhs(1,i,j,k) + dy1ty1 *
   186: |||            >                   (u(1,i,j+1,k) - 2.0d0*u(1,i,j,k) +
   187: |||            >                    u(1,i,j-1,k)) -
   188: |||            >                   ty2 * (u(3,i,j+1,k) - u(3,i,j-1,k))
   189: |||     F                 rhs(2,i,j,k) = rhs(2,i,j,k) + dy2ty1 *
   190: |||            >                   (u(2,i,j+1,k) - 2.0d0*u(2,i,j,k) +
   191: |||            >                    u(2,i,j-1,k)) +
   192: |||            >                   yycon2 * (us(i,j+1,k) - 2.0d0*us(i,j,k) +
   193: |||            >                             us(i,j-1,k)) -
   194: |||            >                   ty2 * (u(2,i,j+1,k)*vp1 -
   195: |||            >                          u(2,i,j-1,k)*vm1)
   196: |||     F                 rhs(3,i,j,k) = rhs(3,i,j,k) + dy3ty1 *
   197: |||            >                   (u(3,i,j+1,k) - 2.0d0*u(3,i,j,k) +
   198: |||            >                    u(3,i,j-1,k)) +
   199: |||            >                   yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
   200: |||            >                   ty2 * (u(3,i,j+1,k)*vp1 -
   201: |||            >                          u(3,i,j-1,k)*vm1 +
   202: |||            >                          (u(5,i,j+1,k) - square(i,j+1,k) -
   203: |||            >                           u(5,i,j-1,k) + square(i,j-1,k))
   204: |||            >                          *c2)
   205: |||     F                 rhs(4,i,j,k) = rhs(4,i,j,k) + dy4ty1 *
   206: |||            >                   (u(4,i,j+1,k) - 2.0d0*u(4,i,j,k) +
   207: |||            >                    u(4,i,j-1,k)) +
   208: |||            >                   yycon2 * (ws(i,j+1,k) - 2.0d0*ws(i,j,k) +
   209: |||            >                             ws(i,j-1,k)) -
   210: |||            >                   ty2 * (u(4,i,j+1,k)*vp1 -
   211: |||            >                          u(4,i,j-1,k)*vm1)
   212: |||     F                 rhs(5,i,j,k) = rhs(5,i,j,k) + dy5ty1 *
   213: |||            >                   (u(5,i,j+1,k) - 2.0d0*u(5,i,j,k) +
   214: |||            >                    u(5,i,j-1,k)) +
   215: |||            >                   yycon3 * (qs(i,j+1,k) - 2.0d0*qs(i,j,k) +
   216: |||            >                             qs(i,j-1,k)) +
   217: |||            >                   yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk +
   218: |||            >                             vm1*vm1) +
   219: |||            >                   yycon5 * (u(5,i,j+1,k)*rho_i(i,j+1,k) -
   220: |||            >                             2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   221: |||            >                             u(5,i,j-1,k)*rho_i(i,j-1,k)) -
   222: |||            >                   ty2 * ((c1*u(5,i,j+1,k) -
   223: |||            >                           c2*square(i,j+1,k)) * vp1 -
   224: |||            >                          (c1*u(5,i,j-1,k) -
   225: |||            >                           c2*square(i,j-1,k)) * vm1)
   226: ||V----                end do
   227: ||        
   228: ||        
   229: ||        c---------------------------------------------------------------------
   230: ||        c      add fourth order eta-direction dissipation
   231: ||        c---------------------------------------------------------------------
   232: ||        
   233: ||                  if (j .eq. 1) then
   234: ||V---->            do     i = 1, nx2
   235: |||*--->               do     m = 1, 5
   236: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   237: ||||           >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
   238: ||||           >                            u(m,i,j+2,k))
   239: |||*---                end do
   240: ||V----             end do
   241: ||        
   242: ||                  else if (j .eq. 2) then
   243: ||V---->            do     i = 1, nx2
   244: |||*--->               do     m = 1, 5
   245: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   246: ||||           >                    (-4.0d0*u(m,i,j-1,k) + 6.0d0*u(m,i,j,k) -
   247: ||||           >                      4.0d0*u(m,i,j+1,k) + u(m,i,j+2,k))
   248: |||*---                end do
   249: ||V----             end do
   250: ||         
   251: ||                  else if (j .eq. ny2-1) then
   252: ||V---->            do     i = 1, nx2
   253: |||*--->               do     m = 1, 5
   254: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   255: ||||           >                    ( u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) +
   256: ||||           >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) )
   257: |||*---                end do
   258: ||V----             end do
   259: ||        
   260: ||                  else if (j .eq. ny2) then
   261: ||V---->            do     i = 1, nx2
   262: |||*--->               do     m = 1, 5
   263: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   264: ||||           >                    ( u(m,i,j-2,k) - 4.d0*u(m,i,j-1,k) +
   265: ||||           >                      5.d0*u(m,i,j,k) )
   266: |||*---                end do
   267: ||V----             end do
   268: ||        
   269: ||                  else  !do    j = 3, ny2-2
   270: ||V---->               do  i = 1,nx2
   271: |||*--->                  do     m = 1, 5
   272: ||||    F                    rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   273: ||||           >                    (  u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) +
   274: ||||           >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
   275: ||||           >                         u(m,i,j+2,k) )
   276: |||*---                   end do
   277: ||V----                end do
   278: ||                  endif
   279: |+-----             end do
   280: Y------          end do
   281:           !$omp end do nowait
   282:           !$omp master
   283:                  if (timeron) call timer_stop(t_rhsy)
   284:           
   285:           c---------------------------------------------------------------------
   286:           c      compute zeta-direction fluxes
   287:           c---------------------------------------------------------------------
   288:                  if (timeron) call timer_start(t_rhsz)
   289:           !$omp end master
   290:           !$omp do schedule(static) collapse(2)
   291: Y------>         do    k = 1, grid_points(3)-2
   292: |+----->            do     j = 1, grid_points(2)-2
   293: ||V---->               do     i = 1, grid_points(1)-2
   294: |||                       wijk = ws(i,j,k)
   295: |||                       wp1  = ws(i,j,k+1)
   296: |||                       wm1  = ws(i,j,k-1)
   297: |||       
   298: |||     F                 rhs(1,i,j,k) = rhs(1,i,j,k) + dz1tz1 *
   299: |||            >                   (u(1,i,j,k+1) - 2.0d0*u(1,i,j,k) +
   300: |||            >                    u(1,i,j,k-1)) -
   301: |||            >                   tz2 * (u(4,i,j,k+1) - u(4,i,j,k-1))
   302: |||     F                 rhs(2,i,j,k) = rhs(2,i,j,k) + dz2tz1 *
   303: |||            >                   (u(2,i,j,k+1) - 2.0d0*u(2,i,j,k) +
   304: |||            >                    u(2,i,j,k-1)) +
   305: |||            >                   zzcon2 * (us(i,j,k+1) - 2.0d0*us(i,j,k) +
   306: |||            >                             us(i,j,k-1)) -
   307: |||            >                   tz2 * (u(2,i,j,k+1)*wp1 -
   308: |||            >                          u(2,i,j,k-1)*wm1)
   309: |||     F                 rhs(3,i,j,k) = rhs(3,i,j,k) + dz3tz1 *
   310: |||            >                   (u(3,i,j,k+1) - 2.0d0*u(3,i,j,k) +
   311: |||            >                    u(3,i,j,k-1)) +
   312: |||            >                   zzcon2 * (vs(i,j,k+1) - 2.0d0*vs(i,j,k) +
   313: |||            >                             vs(i,j,k-1)) -
   314: |||            >                   tz2 * (u(3,i,j,k+1)*wp1 -
   315: |||            >                          u(3,i,j,k-1)*wm1)
   316: |||     F                 rhs(4,i,j,k) = rhs(4,i,j,k) + dz4tz1 *
   317: |||            >                   (u(4,i,j,k+1) - 2.0d0*u(4,i,j,k) +
   318: |||            >                    u(4,i,j,k-1)) +
   319: |||            >                   zzcon2*con43 * (wp1 - 2.0d0*wijk + wm1) -
   320: |||            >                   tz2 * (u(4,i,j,k+1)*wp1 -
   321: |||            >                          u(4,i,j,k-1)*wm1 +
   322: |||            >                          (u(5,i,j,k+1) - square(i,j,k+1) -
   323: |||            >                           u(5,i,j,k-1) + square(i,j,k-1))
   324: |||            >                          *c2)
   325: |||     F                 rhs(5,i,j,k) = rhs(5,i,j,k) + dz5tz1 *
   326: |||            >                   (u(5,i,j,k+1) - 2.0d0*u(5,i,j,k) +
   327: |||            >                    u(5,i,j,k-1)) +
   328: |||            >                   zzcon3 * (qs(i,j,k+1) - 2.0d0*qs(i,j,k) +
   329: |||            >                             qs(i,j,k-1)) +
   330: |||            >                   zzcon4 * (wp1*wp1 - 2.0d0*wijk*wijk +
   331: |||            >                             wm1*wm1) +
   332: |||            >                   zzcon5 * (u(5,i,j,k+1)*rho_i(i,j,k+1) -
   333: |||            >                             2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
   334: |||            >                             u(5,i,j,k-1)*rho_i(i,j,k-1)) -
   335: |||            >                   tz2 * ( (c1*u(5,i,j,k+1) -
   336: |||            >                            c2*square(i,j,k+1))*wp1 -
   337: |||            >                           (c1*u(5,i,j,k-1) -
   338: |||            >                            c2*square(i,j,k-1))*wm1)
   339: ||V----                end do
   340: ||        
   341: ||        c---------------------------------------------------------------------
   342: ||        c      add fourth order zeta-direction dissipation
   343: ||        c---------------------------------------------------------------------
   344: ||        
   345: ||                  if (k .eq. 1) then
   346: ||V---->            do     i = 1, grid_points(1)-2
   347: |||*--->               do     m = 1, 5
   348: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
   349: ||||           >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
   350: ||||           >                            u(m,i,j,k+2))
   351: |||*---                end do
   352: ||V----             end do
   353: ||        
   354: ||                  else if (k .eq. 2) then
   355: ||V---->            do     i = 1, grid_points(1)-2
   356: |||*--->               do     m = 1, 5
   357: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   358: ||||           >                    (-4.0d0*u(m,i,j,k-1) + 6.0d0*u(m,i,j,k) -
   359: ||||           >                      4.0d0*u(m,i,j,k+1) + u(m,i,j,k+2))
   360: |||*---                end do
   361: ||V----             end do
   362: ||         
   363: ||                  else if (k .eq. grid_points(3)-3) then
   364: ||V---->            do     i = 1, grid_points(1)-2
   365: |||*--->               do     m = 1, 5
   366: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   367: ||||           >                    ( u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) +
   368: ||||           >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) )
   369: |||*---                end do
   370: ||V----             end do
   371: ||        
   372: ||                  else if (k .eq. grid_points(3)-2) then
   373: ||V---->            do     i = 1, grid_points(1)-2
   374: |||*--->               do     m = 1, 5
   375: ||||    F                 rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   376: ||||           >                    ( u(m,i,j,k-2) - 4.d0*u(m,i,j,k-1) +
   377: ||||           >                      5.d0*u(m,i,j,k) )
   378: |||*---                end do
   379: ||V----             end do
   380: ||        
   381: ||                  else !do     k = 3, grid_points(3)-4
   382: ||V---->               do     i = 1,grid_points(1)-2
   383: |||*--->                  do     m = 1, 5
   384: ||||    F                    rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
   385: ||||           >                    (  u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) +
   386: ||||           >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
   387: ||||           >                         u(m,i,j,k+2) )
   388: |||*---                   end do
   389: ||V----                end do
   390: ||                  endif
   391: |+-----             end do
   392: Y------          end do
   393:           !$omp end do nowait
   394:           !$omp master
   395:                  if (timeron) call timer_stop(t_rhsz)
   396:           !$omp end master
   397:           
   398:           !$omp do schedule(static) collapse(2)
   399: Y------>         do    k = 1, nz2
   400: |+----->            do    j = 1, ny2
   401: ||V---->               do    i = 1, nx2
   402: |||*--->                  do    m = 1, 5
   403: ||||                         rhs(m,i,j,k) = rhs(m,i,j,k) * dt
   404: |||*---                   end do
   405: ||V----                end do
   406: |+-----             end do
   407: Y------          end do
   408:           !$omp end do nowait
   409:           !$omp end parallel
   410:                   if (timeron) call timer_stop(t_rhs)
   411:            
   412:                  return
   413:                  end


