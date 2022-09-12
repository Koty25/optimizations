NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:59 2021
FILE NAME: diffuse.f

PROCEDURE NAME: DIFFUSION
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    25: par(1801): Parallel routine generated.: DIFFUSION$1
    26: par(1803): Parallelized by "do".
    26: vec( 101): Vectorized loop.
    38: par(1801): Parallel routine generated.: DIFFUSION$2
    40: par(1803): Parallelized by "do".
    42: vec( 101): Vectorized loop.
    43: opt(1418): Constant-length loop is expanded.
    45: vec( 128): Fused multiply-add operation applied.
    45: vec( 126): Idiom detected.: SUM
    54: par(1803): Parallelized by "do".
    54: vec( 101): Vectorized loop.
    56: vec( 128): Fused multiply-add operation applied.
    56: vec( 126): Idiom detected.: SUM
    59: par(1807): Critical-construct generated.
    59: par(1809): Barrier synchronization.
    65: vec( 110): Vectorization obstructive procedure reference.: ADDS2M1
    65: vec( 110): Vectorization obstructive procedure reference.: ADDS1M1
    65: vec( 110): Vectorization obstructive procedure reference.: TRANSFB
    65: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
    65: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
    65: vec( 110): Vectorization obstructive procedure reference.: TRANSF
    65: vec( 109): Vectorization obstructive statement.
    65: vec( 103): Unvectorized loop.
    68: par(1801): Parallel routine generated.: DIFFUSION$3
    72: par(1803): Parallelized by "do".
    74: vec( 101): Vectorized loop.
    75: opt(1418): Constant-length loop is expanded.
    77: vec( 128): Fused multiply-add operation applied.
    77: vec( 126): Idiom detected.: SUM
    85: par(1803): Parallelized by "do".
    85: vec( 101): Vectorized loop.
    87: vec( 128): Fused multiply-add operation applied.
    87: vec( 126): Idiom detected.: SUM
    90: par(1807): Critical-construct generated.
    90: par(1809): Barrier synchronization.
   108: par(1801): Parallel routine generated.: DIFFUSION$4
   109: vec( 110): Vectorization obstructive procedure reference.: LAPLACIAN
   109: par(1803): Parallelized by "do".
   109: vec( 103): Unvectorized loop.
   121: par(1801): Parallel routine generated.: DIFFUSION$5
   122: par(1803): Parallelized by "do".
   123: vec( 103): Unvectorized loop.
   123: vec( 110): Vectorization obstructive procedure reference.: FACEV
   133: par(1801): Parallel routine generated.: DIFFUSION$6
   135: par(1803): Parallelized by "do".
   135: vec( 103): Unvectorized loop.
   135: vec( 113): Overhead of loop division is too large.
   136: vec( 101): Vectorized loop.
   137: opt(1418): Constant-length loop is expanded.
   138: opt(1418): Constant-length loop is expanded.
   139: vec( 128): Fused multiply-add operation applied.
   139: vec( 126): Idiom detected.: SUM
   147: par(1803): Parallelized by "do".
   147: vec( 101): Vectorized loop.
   149: vec( 128): Fused multiply-add operation applied.
   149: vec( 126): Idiom detected.: SUM
   152: par(1807): Critical-construct generated.
   152: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:59 2021
FILE NAME: diffuse.f

PROCEDURE NAME: DIFFUSION
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine diffusion(ifmortar)
     3:           c---------------------------------------------------------------------
     4:           c     advance the diffusion term using CG iterations
     5:           c---------------------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:           
    10:                 double precision  rho_aux, rho1, rho2, beta, cona
    11:                 logical ifmortar
    12:                 integer iter,ie, im,iside,i,j,k
    13:           
    14:                 if (timeron) call timer_start(t_diffusion)
    15:           c.....set up diagonal preconditioner
    16:                 if (ifmortar) then
    17:                   call setuppc
    18:                   call setpcmo
    19:                 end if
    20:           
    21:           c.....arrays t and umor are accumlators of (am pm) in the CG algorithm
    22:           c     (see the specification)
    23:           
    24:                 call r_init_omp(t,ntot,0.d0)
    25:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    26: Y------>        do i=1,nmor
    27: |                 umor(i)=0.d0
    28: Y------         end do
    29:           c$OMP END PARALLEL DO
    30:           
    31:           c.....calculate initial am (see specification) in CG algorithm
    32:           
    33:           c.....trhs and rmor are combined to generate r0 in CG algorithm.
    34:           c     pdiff and pmorx are combined to generate q0 in the CG algorithm.
    35:           c     rho1 is  (qm,rm) in the CG algorithm.
    36:           
    37:                 rho1 = 0.d0
    38:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(im,ie,i,j,k) REDUCTION(+:rho1)
    39:           c$OMP DO
    40: P------>         do ie=1,nelt
    41: |+----->           do k=1,lx1
    42: ||V---->             do j=1,lx1
    43: |||*--->               do i=1,lx1
    44: ||||                     pdiff(i,j,k,ie) = dpcelm(i,j,k,ie)*trhs(i,j,k,ie)
    45: ||||    F                rho1            = rho1 + trhs(i,j,k,ie)*pdiff(i,j,k,ie)*
    46: ||||           &                                          tmult(i,j,k,ie)
    47: |||*---                end do
    48: ||V----              end do
    49: |+-----            end do
    50: P------          end do
    51:           c$OMP END DO nowait
    52:           
    53:           c$OMP DO
    54: Y------>        do im = 1, nmor
    55: |                 pmorx(im) = dpcmor(im)*rmor(im)
    56: |       F         rho1      = rho1 + rmor(im)*pmorx(im)
    57: Y------         end do
    58:           c$OMP END DO nowait
    59:           c$OMP END PARALLEL
    60:           
    61:           c.................................................................
    62:           c     commence conjugate gradient iteration
    63:           c.................................................................
    64:           
    65: +------>        do iter=1, nmxh
    66: |                 if(iter.gt.1) then
    67: |                   rho_aux = 0.d0
    68: |         c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(im,ie,i,j,k) REDUCTION(+:rho_aux)
    69: |         c$OMP DO
    70: |         c.........pdiffp and ppmor are combined to generate q_m+1 in the specification
    71: |         c         rho_aux is (q_m+1,r_m+1)
    72: |P----->            do ie = 1, nelt
    73: ||+---->              do k=1,lx1
    74: |||V--->                do j=1,lx1
    75: ||||*-->                  do i=1,lx1
    76: |||||                       pdiffp(i,j,k,ie) = dpcelm(i,j,k,ie)*trhs(i,j,k,ie)
    77: |||||   F                   rho_aux =rho_aux+trhs(i,j,k,ie)*pdiffp(i,j,k,ie)*
    78: |||||          &                                            tmult(i,j,k,ie)
    79: ||||*--                   end do
    80: |||V---                 end do
    81: ||+----               end do
    82: |P-----             end do
    83: |         c$OMP END DO nowait
    84: |         c$OMP DO
    85: |Y----->            do im = 1, nmor
    86: ||                    ppmor(im) = dpcmor(im)*rmor(im)
    87: ||      F             rho_aux = rho_aux + rmor(im)*ppmor(im)
    88: |Y-----             end do
    89: |         c$OMP END DO nowait
    90: |         c$OMP END PARALLEL
    91: |         
    92: |         c.........compute bm (beta) in the specification
    93: |                   rho2 = rho1
    94: |                   rho1 = rho_aux
    95: |                   beta = rho1/rho2
    96: |         c.........update p_m+1 in the specification
    97: |                   call adds1m1(pdiff, pdiffp, beta,ntot)
    98: |                   call adds1m1(pmorx, ppmor,  beta, nmor)
    99: |                 end if
   100: |          
   101: |         c.......compute matrix vector product: (theta pm) in the specification
   102: |         
   103: |                 if (timeron) call timer_start(t_transf)
   104: |                 call transf(pmorx,pdiff)
   105: |                 if (timeron) call timer_stop(t_transf)
   106: |         
   107: |         c.......compute pdiffp which is (A theta pm) in the specification
   108: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ie)
   109: |P----->          do ie=1, nelt
   110: ||                  call laplacian(pdiffp(1,1,1,ie),pdiff(1,1,1,ie),size_e(ie))
   111: |P-----           end do
   112: |         c$OMP END PARALLEL DO
   113: |         
   114: |         c.......compute ppmor which will be used to compute (thetaT A theta pm)
   115: |         c       in the specification
   116: |                 if (timeron) call timer_start(t_transfb)
   117: |                 call transfb(ppmor,pdiffp)
   118: |                 if (timeron) call timer_stop(t_transfb)
   119: |          
   120: |         c.......apply boundary condition
   121: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ie,iside)
   122: |P----->          do ie=1,nelt
   123: ||+---->            do iside=1,nsides
   124: |||                   if(cbc(iside,ie).eq.0)then
   125: |||                     call facev(pdiffp(1,1,1,ie),iside,0.d0)
   126: |||                   end if
   127: ||+----             end do
   128: |P-----           end do
   129: |         c$OMP END PARALLEL DO
   130: |         
   131: |         c.......compute cona which is (pm,theta T A theta pm)
   132: |                 cona = 0.d0
   133: |         c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(im,ie,i,j,k) REDUCTION(+:cona)
   134: |         c$OMP DO
   135: |P----->          do ie = 1, nelt
   136: ||V---->            do k=1,lx1
   137: |||*--->              do j=1,lx1
   138: ||||*-->                do i=1,lx1
   139: |||||   F                 cona = cona +
   140: |||||          &          pdiff(i,j,k,ie)*pdiffp(i,j,k,ie)*tmult(i,j,k,ie)
   141: ||||*--                 end do
   142: |||*---                end do
   143: ||V----             end do
   144: |P-----           end do
   145: |         c$OMP END DO nowait
   146: |         c$OMP DO
   147: |Y----->          do im = 1, nmor
   148: ||                  ppmor(im) = ppmor(im)*tmmor(im)
   149: ||      F           cona = cona + pmorx(im)*ppmor(im)
   150: |Y-----           end do
   151: |         c$OMP END DO nowait
   152: |         c$OMP END PARALLEL
   153: |         
   154: |         c.......compute am
   155: |                 cona = rho1/cona
   156: |         c.......compute (am pm)
   157: |                 call adds2m1(t,    pdiff,   cona, ntot)
   158: |                 call adds2m1(umor, pmorx,   cona, nmor)
   159: |         c.......compute r_m+1
   160: |                 call adds2m1(trhs, pdiffp, -cona, ntot)
   161: |                 call adds2m1(rmor, ppmor,  -cona, nmor)
   162: |          
   163: +------         end do
   164:           
   165:                 if (timeron) call timer_start(t_transf)
   166:                 call transf(umor,t)
   167:                 if (timeron) call timer_stop(t_transf)
   168:                 if (timeron) call timer_stop(t_diffusion)
   169:           
   170:                 return
   171:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:59 2021
FILE NAME: diffuse.f

PROCEDURE NAME: LAPLACIAN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   191: vec( 103): Unvectorized loop.
   191: vec( 113): Overhead of loop division is too large.
   192: vec( 101): Vectorized loop.
   193: opt(1418): Constant-length loop is expanded.
   194: opt(1418): Constant-length loop is expanded.
   195: vec( 126): Idiom detected.: SUM
   202: vec( 103): Unvectorized loop.
   202: vec( 113): Overhead of loop division is too large.
   203: vec( 101): Vectorized loop.
   204: opt(1418): Constant-length loop is expanded.
   205: opt(1418): Constant-length loop is expanded.
   206: vec( 126): Idiom detected.: SUM
   214: vec( 101): Vectorized loop.
   215: opt(1418): Constant-length loop is expanded.
   216: opt(1418): Constant-length loop is expanded.
   217: vec( 128): Fused multiply-add operation applied.
   231: vec( 101): Vectorized loop.
   232: opt(1418): Constant-length loop is expanded.
   233: opt(1418): Constant-length loop is expanded.
   234: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:59 2021
FILE NAME: diffuse.f

PROCEDURE NAME: LAPLACIAN
FORMAT LIST

 LINE   LOOP      STATEMENT

   175:                 subroutine laplacian(r,u,sizei)
   176:           c------------------------------------------------------------------
   177:           c     compute  r = visc*[A]x +[B]x on a given element.
   178:           c------------------------------------------------------------------
   179:           
   180:                 use ua_data
   181:                 implicit none
   182:           
   183:                 double precision r(lx1,lx1,lx1), u(lx1,lx1,lx1), rdtime
   184:                 integer i,j,k, ix,iz, sizei
   185:           
   186:                 double precision tm1(lx1,lx1,lx1),tm2(lx1,lx1,lx1)
   187:           
   188:                 rdtime = 1.d0/dtime
   189:           
   190:                 call r_init(tm1,nxyz,0.d0)
   191: +------>        do iz=1,lx1
   192: |V----->          do k = 1, lx1
   193: ||*---->            do j = 1, lx1
   194: |||*--->              do i = 1, lx1
   195: ||||                    tm1(i,j,iz) = tm1(i,j,iz)+wdtdr(i,k)*u(k,j,iz)
   196: |||*---               end do
   197: ||*----             end do
   198: |V-----           end do
   199: +------         end do
   200:            
   201:                 call r_init(tm2,nxyz,0.d0)
   202: +------>        do iz=1,lx1
   203: |V----->          do k = 1, lx1
   204: ||*---->            do j = 1, lx1
   205: |||*--->              do i = 1, lx1
   206: ||||                    tm2(i,j,iz) = tm2(i,j,iz)+u(i,k,iz)*wdtdr(k,j)
   207: |||*---               end do
   208: ||*----             end do
   209: |V-----           end do
   210: +------         end do
   211:            
   212:                 call r_init(r,nxyz,0.d0)
   213: +------>        do k = 1, lx1
   214: |V----->          do iz=1, lx1
   215: ||*---->            do j = 1, lx1
   216: |||*--->              do i = 1, lx1
   217: ||||    F               r(i,j,iz) = r(i,j,iz)+u(i,j,k)*wdtdr(k,iz)
   218: |||*---               end do
   219: ||*----             end do
   220: |V-----           end do
   221: +------         end do
   222:           
   223:           c.....collocate with remaining weights and sum to complete factorization.
   224:            
   225:           c      do ix=1,nxyz
   226:           c         r(ix,1,1)=visc*(tm1(ix,1,1)*g4m1_s(ix,1,1,sizei)+
   227:           c     &                   tm2(ix,1,1)*g5m1_s(ix,1,1,sizei)+
   228:           c     &                     r(ix,1,1)*g6m1_s(ix,1,1,sizei))+
   229:           c     &               bm1_s(ix,1,1,sizei)*rdtime*u(ix,1,1)
   230:           c      end do
   231: V------>        do k=1,lx1
   232: |*----->          do j=1,lx1
   233: ||*---->            do i=1,lx1
   234: |||     F             r(i,j,k)=visc*(tm1(i,j,k)*g4m1_s(i,j,k,sizei)+
   235: |||            &                   tm2(i,j,k)*g5m1_s(i,j,k,sizei)+
   236: |||            &                    r(i,j,k)*g6m1_s(i,j,k,sizei))+
   237: |||            &               bm1_s(i,j,k,sizei)*rdtime*u(i,j,k)
   238: ||*----             end do
   239: |*-----           end do
   240: V------         end do
   241:           
   242:                 return
   243:                 end


