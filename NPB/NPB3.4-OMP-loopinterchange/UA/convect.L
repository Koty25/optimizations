NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:45 2020
FILE NAME: convect.f

PROCEDURE NAME: CONVECT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    33: opt(1418): Constant-length loop is expanded.
    39: par(1801): Parallel routine generated.: CONVECT$1
    43: par(1803): Parallelized by "do".
    48: opt(1418): Constant-length loop is expanded.
    51: opt(1418): Constant-length loop is expanded.
    54: opt(1418): Constant-length loop is expanded.
    58: vec( 102): Partially vectorized loop.
    59: vec( 101): Vectorized loop.
    60: vec( 101): Vectorized loop.
    61: vec( 128): Fused multiply-add operation applied.
    69: opt(1418): Constant-length loop is expanded.
    70: vec( 128): Fused multiply-add operation applied.
    74: opt(1418): Constant-length loop is expanded.
    75: vec( 128): Fused multiply-add operation applied.
    79: opt(1418): Constant-length loop is expanded.
    80: vec( 128): Fused multiply-add operation applied.
    84: vec( 128): Fused multiply-add operation applied.
    85: vec( 128): Fused multiply-add operation applied.
    91: vec( 102): Partially vectorized loop.
    92: vec( 101): Vectorized loop.
    93: vec( 101): Vectorized loop.
    94: vec( 128): Fused multiply-add operation applied.
   102: opt(1418): Constant-length loop is expanded.
   103: vec( 128): Fused multiply-add operation applied.
   107: opt(1418): Constant-length loop is expanded.
   108: vec( 128): Fused multiply-add operation applied.
   112: opt(1418): Constant-length loop is expanded.
   113: vec( 128): Fused multiply-add operation applied.
   117: vec( 128): Fused multiply-add operation applied.
   118: vec( 128): Fused multiply-add operation applied.
   123: vec( 102): Partially vectorized loop.
   124: vec( 101): Vectorized loop.
   125: vec( 101): Vectorized loop.
   126: vec( 128): Fused multiply-add operation applied.
   134: opt(1418): Constant-length loop is expanded.
   135: vec( 128): Fused multiply-add operation applied.
   139: opt(1418): Constant-length loop is expanded.
   140: vec( 128): Fused multiply-add operation applied.
   144: opt(1418): Constant-length loop is expanded.
   145: vec( 128): Fused multiply-add operation applied.
   149: vec( 128): Fused multiply-add operation applied.
   150: vec( 128): Fused multiply-add operation applied.
   155: vec( 102): Partially vectorized loop.
   156: vec( 101): Vectorized loop.
   157: vec( 101): Vectorized loop.
   158: vec( 128): Fused multiply-add operation applied.
   166: opt(1418): Constant-length loop is expanded.
   167: vec( 128): Fused multiply-add operation applied.
   171: opt(1418): Constant-length loop is expanded.
   172: vec( 128): Fused multiply-add operation applied.
   176: opt(1418): Constant-length loop is expanded.
   177: vec( 128): Fused multiply-add operation applied.
   181: vec( 128): Fused multiply-add operation applied.
   182: vec( 128): Fused multiply-add operation applied.
   189: vec( 103): Unvectorized loop.
   189: vec( 110): Vectorization obstructive procedure reference.: FACEV
   196: vec( 101): Vectorized loop.
   197: opt(1418): Constant-length loop is expanded.
   198: vec( 128): Fused multiply-add operation applied.
   200: vec( 128): Fused multiply-add operation applied.
   218: par(1801): Parallel routine generated.: CONVECT$2
   219: par(1803): Parallelized by "do".
   219: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:45 2020
FILE NAME: convect.f

PROCEDURE NAME: CONVECT
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine convect(ifmortar)
     3:           c---------------------------------------------------------
     4:           c     Advance the convection term using 4th order RK
     5:           c     1.ta1 is solution from last time step
     6:           c     2.the heat source is considered part of d/dx
     7:           c     3.trhs is right hand side for the diffusion equation
     8:           c     4.tmor is solution on mortar points, which will be used
     9:           c       as the initial guess when advancing the diffusion term
    10:           c---------------------------------------------------------
    11:           
    12:                 use ua_data
    13:                 implicit none
    14:           
    15:                 double precision alpha2, tempa(lx1,lx1,lx1),
    16:                &       rdtime, pidivalpha, sixth,
    17:                &       dtx1, dtx2, dtx3, src, rk1(lx1,lx1,lx1), rk2(lx1,lx1,lx1),
    18:                &       rk3(lx1,lx1,lx1), rk4(lx1,lx1,lx1), temp(lx1,lx1,lx1),
    19:                &       subtime(3), xx0(3), yy0(3), zz0(3), dtime2, r2, sum,
    20:                &       xloc(lx1), yloc(lx1), zloc(lx1)
    21:                 integer k,iel,i,j,iside,isize, substep, ip
    22:                 logical ifmortar
    23:                 parameter (sixth=1.d0/6.d0)
    24:           
    25:                 if (timeron) call timer_start(t_convect)
    26:                 pidivalpha = dacos(-1.d0)/alpha
    27:                 alpha2     = alpha*alpha
    28:                 dtime2     = dtime/2.d0
    29:                 rdtime     = 1.d0/dtime
    30:                 subtime(1) = time
    31:                 subtime(2) = time+dtime2
    32:                 subtime(3) = time+dtime
    33: *------>        do substep = 1, 3
    34: |                 xx0(substep) = x00+velx*subtime(substep)
    35: |                 yy0(substep) = y00+vely*subtime(substep)
    36: |                 zz0(substep) = z00+velz*subtime(substep)
    37: *------         end do
    38:           
    39:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rk4,rk3,rk2,temp,rk1,dtx3,
    40:           c$OMP$ dtx2,dtx1,iside,ip,sum,src,r2,i,j,k,isize,iel,tempa,
    41:           c$OMP$ xloc,yloc,zloc)
    42:           
    43: P------>        do iel = 1, nelt
    44: |                 isize=size_e(iel)
    45: |         c.......xloc(i) is the location of i'th collocation in x direction in an element.
    46: |         c       yloc(i) is the location of j'th collocation in y direction in an element.
    47: |         c       zloc(i) is the location of k'th collocation in z direction in an element.
    48: |*----->          do i = 1, lx1
    49: ||                  xloc(i) = xfrac(i)*(xc(2,iel)-xc(1,iel))+xc(1,iel)
    50: |*-----           end do
    51: |*----->          do j = 1, lx1
    52: ||                  yloc(j) = xfrac(j)*(yc(4,iel)-yc(1,iel))+yc(1,iel)
    53: |*-----           end do
    54: |*----->          do k = 1, lx1
    55: ||                  zloc(k) = xfrac(k)*(zc(5,iel)-zc(1,iel))+zc(1,iel)
    56: |*-----           end do
    57: |         
    58: |S----->          do k = 1, lx1
    59: ||V---->            do j = 1, lx1
    60: |||V--->              do i = 1, lx1
    61: ||||    F               r2 = (xloc(i)-xx0(1))**2+(yloc(j)-yy0(1))**2+
    62: ||||           &             (zloc(k)-zz0(1))**2
    63: ||||                    if (r2.le.alpha2) then
    64: ||||                      src = dcos(dsqrt(r2)*pidivalpha)+1.d0
    65: ||||                    else
    66: ||||                      src = 0.d0
    67: ||||                    endif
    68: ||||                    sum = 0.d0
    69: ||||*-->                do ip = 1, lx1
    70: |||||   F                 sum = sum + dxm1(i,ip) * ta1(ip,j,k,iel)
    71: ||||*--                 end do
    72: ||||                    dtx1 = -velx*sum*xrm1_s(i,j,k,isize)
    73: ||||                    sum  = 0.d0
    74: ||||*-->                do ip = 1, lx1
    75: |||||   F                 sum = sum + dxm1(j,ip) * ta1(i,ip,k,iel)
    76: ||||*--                 end do
    77: ||||                    dtx2=-vely*sum*xrm1_s(i,j,k,isize)
    78: ||||                    sum = 0.d0
    79: ||||*-->                do ip = 1, lx1
    80: |||||   F                 sum = sum + dxm1(k,ip) * ta1(i,j,ip,iel)
    81: ||||*--                 end do
    82: ||||                    dtx3=-velz*sum*xrm1_s(i,j,k,isize)
    83: ||||      
    84: ||||    F               rk1(i,j,k)= dtx1 + dtx2 + dtx3 + src
    85: ||||    F               temp(i,j,k)=ta1(i,j,k,iel)+dtime2*rk1(i,j,k)
    86: ||||      
    87: |||V---               end do
    88: ||V----             end do
    89: |S-----           end do
    90: |         
    91: |S----->          do k = 1, lx1
    92: ||V---->            do j = 1, lx1
    93: |||V--->              do i = 1, lx1
    94: ||||    F               r2 = (xloc(i)-xx0(2))**2 + (yloc(j)-yy0(2))**2 +
    95: ||||           &             (zloc(k)-zz0(2))**2
    96: ||||                    if (r2.le.alpha2) then
    97: ||||                      src = dcos(dsqrt(r2)*pidivalpha)+1.d0
    98: ||||                    else
    99: ||||                      src = 0.d0
   100: ||||                    endif
   101: ||||                    sum = 0.d0
   102: ||||*-->                do ip = 1, lx1
   103: |||||   F                 sum = sum + dxm1(i,ip) * temp(ip,j,k)
   104: ||||*--                 end do
   105: ||||                    dtx1 =-velx*sum*xrm1_s(i,j,k,isize)
   106: ||||                    sum = 0.d0
   107: ||||*-->                do ip = 1, lx1
   108: |||||   F                 sum = sum + dxm1(j,ip) * temp(i,ip,k)
   109: ||||*--                 end do
   110: ||||                    dtx2 =-vely*sum*xrm1_s(i,j,k,isize)
   111: ||||                    sum = 0.d0
   112: ||||*-->                do ip = 1, lx1
   113: |||||   F                 sum = sum + dxm1(k,ip) * temp(i,j,ip)
   114: ||||*--                 end do
   115: ||||                    dtx3 =-velz*sum*xrm1_s(i,j,k,isize)
   116: ||||      
   117: ||||    F               rk2(i,j,k)= dtx1 + dtx2 + dtx3 + src
   118: ||||    F               tempa(i,j,k)=ta1(i,j,k,iel)+dtime2*rk2(i,j,k)
   119: |||V---               end do
   120: ||V----             end do
   121: |S-----           end do
   122: |         
   123: |S----->          do k = 1, lx1
   124: ||V---->            do j = 1, lx1
   125: |||V--->              do i = 1, lx1
   126: ||||    F               r2 = (xloc(i)-xx0(2))**2 + (yloc(j)-yy0(2))**2 +
   127: ||||           &             (zloc(k)-zz0(2))**2
   128: ||||                    if (r2.le.alpha2) then
   129: ||||                      src = dcos(dsqrt(r2)*pidivalpha)+1.d0
   130: ||||                    else
   131: ||||                      src = 0.d0
   132: ||||                    endif
   133: ||||                    sum = 0.d0
   134: ||||*-->                do ip = 1, lx1
   135: |||||   F                 sum = sum + dxm1(i,ip) * tempa(ip,j,k)
   136: ||||*--                 end do
   137: ||||                    dtx1 =-velx*sum*xrm1_s(i,j,k,isize)
   138: ||||                    sum = 0.d0
   139: ||||*-->                do ip = 1, lx1
   140: |||||   F                 sum = sum + dxm1(j,ip) * tempa(i,ip,k)
   141: ||||*--                 end do
   142: ||||                    dtx2 =-vely*sum*xrm1_s(i,j,k,isize)
   143: ||||                    sum = 0.d0
   144: ||||*-->                do ip = 1, lx1
   145: |||||   F                 sum = sum + dxm1(k,ip) * tempa(i,j,ip)
   146: ||||*--                 end do
   147: ||||                    dtx3 =-velz*sum*xrm1_s(i,j,k,isize)
   148: ||||      
   149: ||||    F               rk3(i,j,k)= dtx1 + dtx2 + dtx3 + src
   150: ||||    F               temp(i,j,k)=ta1(i,j,k,iel)+dtime*rk3(i,j,k)
   151: |||V---               end do
   152: ||V----             end do
   153: |S-----           end do
   154: |         
   155: |S----->          do k = 1, lx1
   156: ||V---->            do j = 1, lx1
   157: |||V--->              do i = 1, lx1
   158: ||||    F               r2 = (xloc(i)-xx0(3))**2 + (yloc(j)-yy0(3))**2 +
   159: ||||           &             (zloc(k)-zz0(3))**2
   160: ||||                    if (r2.le.alpha2) then
   161: ||||                      src = dcos(dsqrt(r2)*pidivalpha)+1.d0
   162: ||||                    else
   163: ||||                      src = 0.d0
   164: ||||                    endif
   165: ||||                    sum = 0.d0
   166: ||||*-->                do ip = 1, lx1
   167: |||||   F                 sum = sum + dxm1(i,ip) * temp(ip,j,k)
   168: ||||*--                 end do
   169: ||||                    dtx1 =-velx*sum*xrm1_s(i,j,k,isize)
   170: ||||                    sum = 0.d0
   171: ||||*-->                do ip = 1, lx1
   172: |||||   F                 sum = sum + dxm1(j,ip) * temp(i,ip,k)
   173: ||||*--                 end do
   174: ||||                    dtx2 =-vely*sum*xrm1_s(i,j,k,isize)
   175: ||||                    sum = 0.d0
   176: ||||*-->                do ip = 1, lx1
   177: |||||   F                 sum = sum + dxm1(k,ip) * temp(i,j,ip)
   178: ||||*--                 end do
   179: ||||                    dtx3 =-velz*sum*xrm1_s(i,j,k,isize)
   180: ||||      
   181: ||||    F               rk4(i,j,k)= dtx1 + dtx2 + dtx3 + src
   182: ||||    F               tempa(i,j,k)=sixth*(rk1(i,j,k)+2.d0*
   183: ||||           &                   rk2(i,j,k)+2.d0*rk3(i,j,k)+rk4(i,j,k))
   184: |||V---               end do
   185: ||V----             end do
   186: |S-----           end do
   187: |         
   188: |         c.......apply boundary condition
   189: |+----->          do iside=1,nsides
   190: ||                  if(cbc(iside,iel).eq.0)then
   191: ||                    call facev(tempa,iside,0.0d0)
   192: ||                  end if
   193: |+-----           end do
   194: |          
   195: |+----->          do k=1,lx1
   196: ||V---->            do j=1,lx1
   197: |||*--->              do i=1,lx1
   198: ||||    F               trhs(i,j,k,iel)=bm1_s(i,j,k,isize)*(ta1(i,j,k,iel)*rdtime+
   199: ||||           &                        tempa(i,j,k))
   200: ||||    F               ta1(i,j,k,iel)=ta1(i,j,k,iel)+tempa(i,j,k)*dtime
   201: |||*---               end do
   202: ||V----             end do
   203: |+-----           end do
   204: |         
   205: P------         end do
   206:           c$OMP END PARALLEL DO
   207:           
   208:           c.....get mortar for intial guess for CG
   209:           
   210:                 if (timeron) call timer_start(t_transfb_c)
   211:                 if(ifmortar)then
   212:                   call transfb_c_2(ta1)
   213:                 else
   214:                   call transfb_c(ta1)
   215:                 end if
   216:                 if (timeron) call timer_stop(t_transfb_c)
   217:           
   218:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   219: Y------>        do i=1,nmor
   220: |                tmort(i)=tmort(i)/mormult(i)
   221: Y------         end do
   222:           c$OMP END PARALLEL DO
   223:                 if (timeron) call timer_stop(t_convect)
   224:           
   225:                 return
   226:                 end


