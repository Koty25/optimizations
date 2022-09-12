NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: CG
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    77: vec( 103): Unvectorized loop.
    77: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   167: par(1801): Parallel routine generated.: MAIN_$1
   178: par(1809): Barrier synchronization.
   190: par(1803): Parallelized by "do".
   190: vec( 101): Vectorized loop.
   191: vec( 101): Vectorized loop.
   201: par(1803): Parallelized by "do".
   201: vec( 101): Vectorized loop.
   206: par(1803): Parallelized by "do".
   206: vec( 101): Vectorized loop.
   222: vec( 103): Unvectorized loop.
   222: vec( 109): Vectorization obstructive statement.
   222: vec( 110): Vectorization obstructive procedure reference.: CONJ_GRAD
   237: par(1801): Parallel routine generated.: MAIN_$2
   239: par(1807): Critical-construct generated.
   239: par(1809): Barrier synchronization.
   239: par(1803): Parallelized by "do".
   239: vec( 101): Vectorized loop.
   240: vec( 128): Fused multiply-add operation applied.
   240: vec( 126): Idiom detected.: SUM
   241: vec( 128): Fused multiply-add operation applied.
   241: vec( 126): Idiom detected.: SUM
   252: par(1803): Parallelized by "do".
   252: vec( 101): Vectorized loop.
   268: par(1801): Parallel routine generated.: MAIN_$3
   269: par(1803): Parallelized by "do".
   269: vec( 101): Vectorized loop.
   288: vec( 103): Unvectorized loop.
   288: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
   288: vec( 109): Vectorization obstructive statement.
   288: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
   288: vec( 110): Vectorization obstructive procedure reference.: CONJ_GRAD
   306: par(1801): Parallel routine generated.: MAIN_$4
   308: par(1807): Critical-construct generated.
   308: par(1809): Barrier synchronization.
   308: par(1803): Parallelized by "do".
   308: vec( 101): Vectorized loop.
   309: vec( 128): Fused multiply-add operation applied.
   309: vec( 126): Idiom detected.: SUM
   310: vec( 128): Fused multiply-add operation applied.
   310: vec( 126): Idiom detected.: SUM
   318: par(1805): Master-construct generated.
   331: par(1803): Parallelized by "do".
   331: vec( 101): Vectorized loop.
   415: vec( 103): Unvectorized loop.
   415: vec( 180): I/O statement obstructs vectorization.
   415: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: CG
FORMAT LIST

 LINE   LOOP      STATEMENT

    50:                 program cg
    51:           c---------------------------------------------------------------------
    52:           c---------------------------------------------------------------------
    53:           
    54:                 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    55:           
    56:                 use cg_data
    57:           
    58:                 implicit none
    59:           
    60:           
    61:                 integer            i, j, k, it
    62:           
    63:                 double precision   zeta, randlc
    64:                 external           randlc
    65:                 double precision   rnorm
    66:                 double precision   norm_temp1,norm_temp2,norm_temp3
    67:           
    68:                 double precision   t, mflops, tmax
    69:                 character          class
    70:                 logical            verified
    71:                 double precision   zeta_verify_value, epsilon, err
    72:           
    73:                 character t_names(t_last)*8
    74:           !$    integer   omp_get_max_threads
    75:           !$    external  omp_get_max_threads
    76:           
    77: +------>        do i = 1, T_last
    78: |                  call timer_clear( i )
    79: +------         end do
    80:           
    81:                 call check_timer_flag( timeron )
    82:                 if (timeron) then
    83:                    t_names(t_init) = 'init'
    84:                    t_names(t_bench) = 'total'
    85:                    t_names(t_conj_grad) = 'conjgd'
    86:                 endif
    87:           
    88:                 call timer_start( T_init )
    89:           
    90:                 firstrow = 1
    91:                 lastrow  = na
    92:                 firstcol = 1
    93:                 lastcol  = na
    94:           
    95:           
    96:                 if( na .eq. 1400 .and.
    97:                &    nonzer .eq. 7 .and.
    98:                &    niter .eq. 15 .and.
    99:                &    shift .eq. 10.d0 ) then
   100:                    class = 'S'
   101:                    zeta_verify_value = 8.5971775078648d0
   102:                 else if( na .eq. 7000 .and.
   103:                &         nonzer .eq. 8 .and.
   104:                &         niter .eq. 15 .and.
   105:                &         shift .eq. 12.d0 ) then
   106:                    class = 'W'
   107:                    zeta_verify_value = 10.362595087124d0
   108:                 else if( na .eq. 14000 .and.
   109:                &         nonzer .eq. 11 .and.
   110:                &         niter .eq. 15 .and.
   111:                &         shift .eq. 20.d0 ) then
   112:                    class = 'A'
   113:                    zeta_verify_value = 17.130235054029d0
   114:                 else if( na .eq. 75000 .and.
   115:                &         nonzer .eq. 13 .and.
   116:                &         niter .eq. 75 .and.
   117:                &         shift .eq. 60.d0 ) then
   118:                    class = 'B'
   119:                    zeta_verify_value = 22.712745482631d0
   120:                 else if( na .eq. 150000 .and.
   121:                &         nonzer .eq. 15 .and.
   122:                &         niter .eq. 75 .and.
   123:                &         shift .eq. 110.d0 ) then
   124:                    class = 'C'
   125:                    zeta_verify_value = 28.973605592845d0
   126:                 else if( na .eq. 1500000 .and.
   127:                &         nonzer .eq. 21 .and.
   128:                &         niter .eq. 100 .and.
   129:                &         shift .eq. 500.d0 ) then
   130:                    class = 'D'
   131:                    zeta_verify_value = 52.514532105794d0
   132:                 else if( na .eq. 9000000 .and.
   133:                &         nonzer .eq. 26 .and.
   134:                &         niter .eq. 100 .and.
   135:                &         shift .eq. 1.5d3 ) then
   136:                    class = 'E'
   137:                    zeta_verify_value = 77.522164599383d0
   138:                 else if( na .eq. 54000000 .and.
   139:                &         nonzer .eq. 31 .and.
   140:                &         niter .eq. 100 .and.
   141:                &         shift .eq. 5.0d3 ) then
   142:                    class = 'F'
   143:                    zeta_verify_value = 107.3070826433d0
   144:                 else
   145:                    class = 'U'
   146:                 endif
   147:           
   148:                 write( *,1000 )
   149:                 write( *,1001 ) na
   150:                 write( *,1002 ) niter
   151:           !$    write( *,1003 ) omp_get_max_threads()
   152:                 write( *,* )
   153:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
   154:                >          ' - CG Benchmark', /)
   155:            1001 format(' Size: ', i11 )
   156:            1002 format(' Iterations:                  ', i5 )
   157:            1003 format(' Number of available threads: ', i5)
   158:           
   159:                 naa = na
   160:                 nzz = nz
   161:           
   162:                 call alloc_space
   163:           
   164:           c---------------------------------------------------------------------
   165:           c  Inialize random number generator
   166:           c---------------------------------------------------------------------
   167:           !$omp parallel default(shared) private(i,j,k,zeta)
   168:                 tran    = 314159265.0D0
   169:                 amult   = 1220703125.0D0
   170:                 zeta    = randlc( tran, amult )
   171:           
   172:           c---------------------------------------------------------------------
   173:           c
   174:           c---------------------------------------------------------------------
   175:                 call makea(naa, nzz, a, colidx, rowstr,
   176:                >           firstrow, lastrow, firstcol, lastcol,
   177:                >           arow, acol, aelt, v, iv)
   178:           !$omp barrier
   179:           
   180:           
   181:           c---------------------------------------------------------------------
   182:           c  Note: as a result of the above call to makea:
   183:           c        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
   184:           c        values of colidx which are col indexes go from firstcol --> lastcol
   185:           c        So:
   186:           c        Shift the col index vals from actual (firstcol --> lastcol )
   187:           c        to local, i.e., (1 --> lastcol-firstcol+1)
   188:           c---------------------------------------------------------------------
   189:           !$omp do
   190: Y------>        do j=1,lastrow-firstrow+1
   191: |V----->           do k=rowstr(j),rowstr(j+1)-1
   192: ||                    colidx(k) = colidx(k) - firstcol + 1
   193: |V-----            enddo
   194: Y------         enddo
   195:           !$omp end do nowait
   196:           
   197:           c---------------------------------------------------------------------
   198:           c  set starting vector to (1, 1, .... 1)
   199:           c---------------------------------------------------------------------
   200:           !$omp do
   201: Y------>        do i = 1, na+1
   202: |                  x(i) = 1.0D0
   203: Y------         enddo
   204:           !$omp end do nowait
   205:           !$omp do
   206: Y------>        do j=1, lastcol-firstcol+1
   207: |                  q(j) = 0.0d0
   208: |                  z(j) = 0.0d0
   209: |                  r(j) = 0.0d0
   210: |                  p(j) = 0.0d0
   211: Y------         enddo
   212:           !$omp end do nowait
   213:           !$omp end parallel
   214:           
   215:                 zeta  = 0.0d0
   216:           
   217:           c---------------------------------------------------------------------
   218:           c---->
   219:           c  Do one iteration untimed to init all code and data page tables
   220:           c---->                    (then reinit, start timing, to niter its)
   221:           c---------------------------------------------------------------------
   222: +------>        do it = 1, 1
   223: |         
   224: |         c---------------------------------------------------------------------
   225: |         c  The call to the conjugate gradient routine:
   226: |         c---------------------------------------------------------------------
   227: |                  call conj_grad ( rnorm )
   228: |         
   229: |         c---------------------------------------------------------------------
   230: |         c  zeta = shift + 1/(x.z)
   231: |         c  So, first: (x.z)
   232: |         c  Also, find norm of z
   233: |         c  So, first: (z.z)
   234: |         c---------------------------------------------------------------------
   235: |                  norm_temp1 = 0.0d0
   236: |                  norm_temp2 = 0.0d0
   237: |         !$omp parallel default(shared) private(j,norm_temp3)
   238: |         !$omp do reduction(+:norm_temp1,norm_temp2)
   239: |Y----->           do j=1, lastcol-firstcol+1
   240: ||      F             norm_temp1 = norm_temp1 + x(j)*z(j)
   241: ||      F             norm_temp2 = norm_temp2 + z(j)*z(j)
   242: |Y-----            enddo
   243: |         !$omp end do
   244: |         
   245: |                  norm_temp3 = 1.0d0 / sqrt( norm_temp2 )
   246: |         
   247: |         
   248: |         c---------------------------------------------------------------------
   249: |         c  Normalize z to obtain x
   250: |         c---------------------------------------------------------------------
   251: |         !$omp do
   252: |Y----->           do j=1, lastcol-firstcol+1
   253: ||                    x(j) = norm_temp3*z(j)
   254: |Y-----            enddo
   255: |         !$omp end do nowait
   256: |         !$omp end parallel
   257: |         
   258: |         
   259: +------         enddo                              ! end of do one iteration untimed
   260:           
   261:           
   262:           c---------------------------------------------------------------------
   263:           c  set starting vector to (1, 1, .... 1)
   264:           c---------------------------------------------------------------------
   265:           c
   266:           c
   267:           c
   268:           !$omp parallel do default(shared) private(i)
   269: Y------>        do i = 1, na+1
   270: |                  x(i) = 1.0D0
   271: Y------         enddo
   272:           !$omp end parallel do
   273:           
   274:                 zeta  = 0.0d0
   275:           
   276:                 call timer_stop( T_init )
   277:           
   278:                 write (*, 2000) timer_read(T_init)
   279:            2000 format(' Initialization time = ',f15.3,' seconds')
   280:           
   281:                 call timer_start( T_bench )
   282:           
   283:           c---------------------------------------------------------------------
   284:           c---->
   285:           c  Main Iteration for inverse power method
   286:           c---->
   287:           c---------------------------------------------------------------------
   288: +------>        do it = 1, niter
   289: |         
   290: |         c---------------------------------------------------------------------
   291: |         c  The call to the conjugate gradient routine:
   292: |         c---------------------------------------------------------------------
   293: |                  if ( timeron ) call timer_start( T_conj_grad )
   294: |                  call conj_grad ( rnorm )
   295: |                  if ( timeron ) call timer_stop( T_conj_grad )
   296: |         
   297: |         
   298: |         c---------------------------------------------------------------------
   299: |         c  zeta = shift + 1/(x.z)
   300: |         c  So, first: (x.z)
   301: |         c  Also, find norm of z
   302: |         c  So, first: (z.z)
   303: |         c---------------------------------------------------------------------
   304: |                  norm_temp1 = 0.0d0
   305: |                  norm_temp2 = 0.0d0
   306: |         !$omp parallel default(shared) private(j,norm_temp3)
   307: |         !$omp do reduction(+:norm_temp1,norm_temp2)
   308: |Y----->           do j=1, lastcol-firstcol+1
   309: ||      F             norm_temp1 = norm_temp1 + x(j)*z(j)
   310: ||      F             norm_temp2 = norm_temp2 + z(j)*z(j)
   311: |Y-----            enddo
   312: |         !$omp end do
   313: |         
   314: |         
   315: |                  norm_temp3 = 1.0d0 / sqrt( norm_temp2 )
   316: |         
   317: |         
   318: |         !$omp master
   319: |                  zeta = shift + 1.0d0 / norm_temp1
   320: |                  if( it .eq. 1 ) write( *,9000 )
   321: |                  write( *,9001 ) it, rnorm, zeta
   322: |         !$omp end master
   323: |         
   324: |          9000    format( /,'   iteration           ||r||                 zeta' )
   325: |          9001    format( 4x, i5, 6x, e21.14, f20.13 )
   326: |         
   327: |         c---------------------------------------------------------------------
   328: |         c  Normalize z to obtain x
   329: |         c---------------------------------------------------------------------
   330: |         !$omp do
   331: |Y----->           do j=1, lastcol-firstcol+1
   332: ||                    x(j) = norm_temp3*z(j)
   333: |Y-----            enddo
   334: |         !$omp end do nowait
   335: |         !$omp end parallel
   336: |         
   337: |         
   338: +------         enddo                              ! end of main iter inv pow meth
   339:           
   340:                 call timer_stop( T_bench )
   341:           
   342:           c---------------------------------------------------------------------
   343:           c  End of timed section
   344:           c---------------------------------------------------------------------
   345:           
   346:                 t = timer_read( T_bench )
   347:           
   348:           
   349:                 write(*,100)
   350:            100  format(' Benchmark completed ')
   351:           
   352:                 epsilon = 1.d-10
   353:                 if (class .ne. 'U') then
   354:           
   355:           c         err = abs( zeta - zeta_verify_value)
   356:                    err = abs( zeta - zeta_verify_value )/zeta_verify_value
   357:                    if( (.not.ieee_is_nan(err)) .and. (err .le. epsilon) ) then
   358:                       verified = .TRUE.
   359:                       write(*, 200)
   360:                       write(*, 201) zeta
   361:                       write(*, 202) err
   362:            200        format(' VERIFICATION SUCCESSFUL ')
   363:            201        format(' Zeta is    ', E20.13)
   364:            202        format(' Error is   ', E20.13)
   365:                    else
   366:                       verified = .FALSE.
   367:                       write(*, 300)
   368:                       write(*, 301) zeta
   369:                       write(*, 302) zeta_verify_value
   370:            300        format(' VERIFICATION FAILED')
   371:            301        format(' Zeta                ', E20.13)
   372:            302        format(' The correct zeta is ', E20.13)
   373:                    endif
   374:                 else
   375:                    verified = .FALSE.
   376:                    write (*, 400)
   377:                    write (*, 401)
   378:                    write (*, 201) zeta
   379:            400     format(' Problem size unknown')
   380:            401     format(' NO VERIFICATION PERFORMED')
   381:                 endif
   382:           
   383:           
   384:                 if( t .ne. 0. ) then
   385:                    mflops = 1.0d-6 * 2*niter*dble( na )
   386:                &               * ( 3.+nonzer*dble(nonzer+1)
   387:                &                 + 25.*(5.+nonzer*dble(nonzer+1))
   388:                &                 + 3. ) / t
   389:                 else
   390:                    mflops = 0.d0
   391:                 endif
   392:           
   393:           
   394:                    call print_results('CG', class, na, 0, 0,
   395:                >                      niter, t,
   396:                >                      mflops, '          floating point',
   397:                >                      verified, npbversion, compiletime,
   398:                >                      cs1, cs2, cs3, cs4, cs5, cs6, cs7)
   399:           
   400:           
   401:           
   402:            600  format( i4, 2e19.12)
   403:           
   404:           
   405:           c---------------------------------------------------------------------
   406:           c      More timers
   407:           c---------------------------------------------------------------------
   408:                 if (.not.timeron) goto 999
   409:           
   410:                 tmax = timer_read(T_bench)
   411:                 if (tmax .eq. 0.0) tmax = 1.0
   412:           
   413:                 write(*,800)
   414:            800  format('  SECTION   Time (secs)')
   415: +------>        do i=1, t_last
   416: |                  t = timer_read(i)
   417: |                  if (i.eq.t_init) then
   418: |                     write(*,810) t_names(i), t
   419: |                  else
   420: |                     write(*,810) t_names(i), t, t*100./tmax
   421: |                     if (i.eq.t_conj_grad) then
   422: |                        t = tmax - t
   423: |                        write(*,820) 'rest ', t, t*100./tmax
   424: |                     endif
   425: |                  endif
   426: |          810     format(2x,a8,':',f9.3:'  ( ',f6.2,' %)')
   427: |          820     format(2x,a8,':',f9.3,'  ( ',f6.2,' %)')
   428: +------         end do
   429:           
   430:            999  continue
   431:           
   432:           
   433:                 end                              ! end main


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: CONJ_GRAD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   462: par(1801): Parallel routine generated.: CONJ_GRAD$1
   469: par(1803): Parallelized by "do".
   469: vec( 101): Vectorized loop.
   483: par(1807): Critical-construct generated.
   483: par(1809): Barrier synchronization.
   483: par(1803): Parallelized by "do".
   483: vec( 101): Vectorized loop.
   484: vec( 128): Fused multiply-add operation applied.
   484: vec( 126): Idiom detected.: SUM
   493: vec( 109): Vectorization obstructive statement.
   495: par(1805): Master-construct generated.
   503: par(1809): Barrier synchronization.
   518: par(1803): Parallelized by "do".
   518: vec( 102): Partially vectorized loop.
   520: vec( 101): Vectorized loop.
   521: vec( 128): Fused multiply-add operation applied.
   521: vec( 126): Idiom detected.: SUM
   567: par(1807): Critical-construct generated.
   567: par(1809): Barrier synchronization.
   567: par(1803): Parallelized by "do".
   567: vec( 101): Vectorized loop.
   568: vec( 128): Fused multiply-add operation applied.
   568: vec( 126): Idiom detected.: SUM
   583: par(1807): Critical-construct generated.
   583: par(1809): Barrier synchronization.
   583: par(1803): Parallelized by "do".
   583: vec( 101): Vectorized loop.
   584: vec( 128): Fused multiply-add operation applied.
   585: vec( 128): Fused multiply-add operation applied.
   593: vec( 128): Fused multiply-add operation applied.
   593: vec( 126): Idiom detected.: SUM
   606: par(1803): Parallelized by "do".
   606: vec( 101): Vectorized loop.
   607: vec( 128): Fused multiply-add operation applied.
   621: par(1803): Parallelized by "do".
   621: vec( 102): Partially vectorized loop.
   623: vec( 101): Vectorized loop.
   624: vec( 128): Fused multiply-add operation applied.
   624: vec( 126): Idiom detected.: SUM
   635: par(1807): Critical-construct generated.
   635: par(1809): Barrier synchronization.
   635: par(1803): Parallelized by "do".
   635: vec( 101): Vectorized loop.
   637: vec( 128): Fused multiply-add operation applied.
   637: vec( 126): Idiom detected.: SUM


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: CONJ_GRAD
FORMAT LIST

 LINE   LOOP      STATEMENT

   439:                 subroutine conj_grad ( rnorm )
   440:           c---------------------------------------------------------------------
   441:           c---------------------------------------------------------------------
   442:           
   443:           c---------------------------------------------------------------------
   444:           c  Floaging point arrays here are named as in NPB1 spec discussion of
   445:           c  CG algorithm
   446:           c---------------------------------------------------------------------
   447:            
   448:                 use cg_data
   449:                 implicit none
   450:           
   451:                 integer   j, k
   452:                 integer   cgit, cgitmax
   453:           
   454:                 double precision   d, sum, rho, rho0, alpha, beta, rnorm, suml
   455:           
   456:                 data      cgitmax / 25 /
   457:           
   458:           
   459:                 rho = 0.0d0
   460:                 sum = 0.0d0
   461:           
   462:           !$omp parallel default(shared) private(j,k,cgit,suml,alpha,beta)
   463:           !$omp&  shared(d,rho0,rho,sum)
   464:           
   465:           c---------------------------------------------------------------------
   466:           c  Initialize the CG algorithm:
   467:           c---------------------------------------------------------------------
   468:           !$omp do
   469: Y------>        do j=1,naa+1
   470: |                  q(j) = 0.0d0
   471: |                  z(j) = 0.0d0
   472: |                  r(j) = x(j)
   473: |                  p(j) = r(j)
   474: Y------         enddo
   475:           !$omp end do
   476:           
   477:           
   478:           c---------------------------------------------------------------------
   479:           c  rho = r.r
   480:           c  Now, obtain the norm of r: First, sum squares of r elements locally...
   481:           c---------------------------------------------------------------------
   482:           !$omp do reduction(+:rho)
   483: Y------>        do j=1, lastcol-firstcol+1
   484: |       F          rho = rho + r(j)*r(j)
   485: Y------         enddo
   486:           !$omp end do
   487:           
   488:           c---------------------------------------------------------------------
   489:           c---->
   490:           c  The conj grad iteration loop
   491:           c---->
   492:           c---------------------------------------------------------------------
   493: +------>        do cgit = 1, cgitmax
   494: |         
   495: |         !$omp master
   496: |         c---------------------------------------------------------------------
   497: |         c  Save a temporary of rho and initialize reduction variables
   498: |         c---------------------------------------------------------------------
   499: |                  rho0 = rho
   500: |                  d = 0.d0
   501: |                  rho = 0.d0
   502: |         !$omp end master
   503: |         !$omp barrier
   504: |         
   505: |         c---------------------------------------------------------------------
   506: |         c  q = A.p
   507: |         c  The partition submatrix-vector multiply: use workspace w
   508: |         c---------------------------------------------------------------------
   509: |         C
   510: |         C  NOTE: this version of the multiply is actually (slightly: maybe %5)
   511: |         C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version
   512: |         C        below.   On the Cray t3d, the reverse is true, i.e., the
   513: |         C        unrolled-by-two version is some 10% faster.
   514: |         C        The unrolled-by-8 version below is significantly faster
   515: |         C        on the Cray t3d - overall speed of code is 1.5 times faster.
   516: |         C
   517: |         !$omp do
   518: |Y----->           do j=1,lastrow-firstrow+1
   519: ||                    suml = 0.d0
   520: ||V---->              do k=rowstr(j),rowstr(j+1)-1
   521: |||     G                suml = suml + a(k)*p(colidx(k))
   522: ||V----               enddo
   523: ||                    q(j) = suml
   524: |Y-----            enddo
   525: |         !$omp end do
   526: |         
   527: |         CC          do j=1,lastrow-firstrow+1
   528: |         CC             i = rowstr(j)
   529: |         CC             iresidue = mod( rowstr(j+1)-i, 2 )
   530: |         CC             sum1 = 0.d0
   531: |         CC             sum2 = 0.d0
   532: |         CC             if( iresidue .eq. 1 )
   533: |         CC      &          sum1 = sum1 + a(i)*p(colidx(i))
   534: |         CC             do k=i+iresidue, rowstr(j+1)-2, 2
   535: |         CC                sum1 = sum1 + a(k)  *p(colidx(k))
   536: |         CC                sum2 = sum2 + a(k+1)*p(colidx(k+1))
   537: |         CC             enddo
   538: |         CC             q(j) = sum1 + sum2
   539: |         CC          enddo
   540: |         
   541: |         CC          do j=1,lastrow-firstrow+1
   542: |         CC             i = rowstr(j)
   543: |         CC             iresidue = mod( rowstr(j+1)-i, 8 )
   544: |         CC             suml = 0.d0
   545: |         CC             do k=i,i+iresidue-1
   546: |         CC                suml = suml +  a(k)*p(colidx(k))
   547: |         CC             enddo
   548: |         CC             do k=i+iresidue, rowstr(j+1)-8, 8
   549: |         CC                suml = suml + a(k  )*p(colidx(k  ))
   550: |         CC      &                   + a(k+1)*p(colidx(k+1))
   551: |         CC      &                   + a(k+2)*p(colidx(k+2))
   552: |         CC      &                   + a(k+3)*p(colidx(k+3))
   553: |         CC      &                   + a(k+4)*p(colidx(k+4))
   554: |         CC      &                   + a(k+5)*p(colidx(k+5))
   555: |         CC      &                   + a(k+6)*p(colidx(k+6))
   556: |         CC      &                   + a(k+7)*p(colidx(k+7))
   557: |         CC             enddo
   558: |         CC             q(j) = suml
   559: |         CC          enddo
   560: |          
   561: |         
   562: |         
   563: |         c---------------------------------------------------------------------
   564: |         c  Obtain p.q
   565: |         c---------------------------------------------------------------------
   566: |         !$omp do reduction(+:d)
   567: |Y----->           do j=1, lastcol-firstcol+1
   568: ||      F             d = d + p(j)*q(j)
   569: |Y-----            enddo
   570: |         !$omp end do
   571: |         
   572: |         
   573: |         c---------------------------------------------------------------------
   574: |         c  Obtain alpha = rho / (p.q)
   575: |         c---------------------------------------------------------------------
   576: |                  alpha = rho0 / d
   577: |         
   578: |         c---------------------------------------------------------------------
   579: |         c  Obtain z = z + alpha*p
   580: |         c  and    r = r - alpha*q
   581: |         c---------------------------------------------------------------------
   582: |         !$omp do reduction(+:rho)
   583: |Y----->           do j=1, lastcol-firstcol+1
   584: ||      F             z(j) = z(j) + alpha*p(j)
   585: ||      F             r(j) = r(j) - alpha*q(j)
   586: ||        c         enddo
   587: ||         
   588: ||        c---------------------------------------------------------------------
   589: ||        c  rho = r.r
   590: ||        c  Now, obtain the norm of r: First, sum squares of r elements locally...
   591: ||        c---------------------------------------------------------------------
   592: ||        c         do j=1, lastcol-firstcol+1
   593: ||      F             rho = rho + r(j)*r(j)
   594: |Y-----            enddo
   595: |         !$omp end do
   596: |         
   597: |         c---------------------------------------------------------------------
   598: |         c  Obtain beta:
   599: |         c---------------------------------------------------------------------
   600: |                  beta = rho / rho0
   601: |         
   602: |         c---------------------------------------------------------------------
   603: |         c  p = r + beta*p
   604: |         c---------------------------------------------------------------------
   605: |         !$omp do
   606: |Y----->           do j=1, lastcol-firstcol+1
   607: ||      F             p(j) = r(j) + beta*p(j)
   608: |Y-----            enddo
   609: |         !$omp end do
   610: |         
   611: |         
   612: +------         enddo                             ! end of do cgit=1,cgitmax
   613:           
   614:           
   615:           c---------------------------------------------------------------------
   616:           c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
   617:           c  First, form A.z
   618:           c  The partition submatrix-vector multiply
   619:           c---------------------------------------------------------------------
   620:           !$omp do
   621: Y------>        do j=1,lastrow-firstrow+1
   622: |                  suml = 0.d0
   623: |V----->           do k=rowstr(j),rowstr(j+1)-1
   624: ||      G             suml = suml + a(k)*z(colidx(k))
   625: |V-----            enddo
   626: |                  r(j) = suml
   627: Y------         enddo
   628:           !$omp end do
   629:           
   630:           
   631:           c---------------------------------------------------------------------
   632:           c  At this point, r contains A.z
   633:           c---------------------------------------------------------------------
   634:           !$omp do reduction(+:sum)
   635: Y------>        do j=1, lastcol-firstcol+1
   636: |                  suml = x(j) - r(j)
   637: |       F          sum  = sum + suml*suml
   638: Y------         enddo
   639:           !$omp end do nowait
   640:           !$omp end parallel
   641:           
   642:                 rnorm = sqrt( sum )
   643:           
   644:           
   645:           
   646:                 return
   647:                 end                               ! end of routine conj_grad


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: MAKEA
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   717: vec( 103): Unvectorized loop.
   717: vec( 108): Unvectorizable loop structure.
   738: vec( 110): Vectorization obstructive procedure reference.: VECSET
   738: vec( 110): Vectorization obstructive procedure reference.: SPRNVC
   740: opt(1025): Reference to this procedure inhibits optimization.: SPRNVC
   744: vec( 101): Vectorized loop.
   750: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: MAKEA
FORMAT LIST

 LINE   LOOP      STATEMENT

   653:                 subroutine makea( n, nz, a, colidx, rowstr,
   654:                >                  firstrow, lastrow, firstcol, lastcol,
   655:                >                  arow, acol, aelt, v, iv )
   656:           c---------------------------------------------------------------------
   657:           c---------------------------------------------------------------------
   658:           
   659:                 use tinfo
   660:                 use cg_data, only : nonzer, rcond, shift
   661:           
   662:                 implicit none
   663:           
   664:                 integer             n, nz
   665:                 integer             firstrow, lastrow, firstcol, lastcol
   666:                 integer             colidx(nz), rowstr(n+1)
   667:                 integer             iv(n+nz), arow(n), acol(nonzer+1,n)
   668:                 double precision    aelt(nonzer+1,n), v(nz)
   669:                 double precision    a(nz)
   670:           
   671:           c---------------------------------------------------------------------
   672:           c       generate the test problem for benchmark 6
   673:           c       makea generates a sparse matrix with a
   674:           c       prescribed sparsity distribution
   675:           c
   676:           c       parameter    type        usage
   677:           c
   678:           c       input
   679:           c
   680:           c       n            i           number of cols/rows of matrix
   681:           c       nz           i           nonzeros as declared array size
   682:           c       rcond        r*8         condition number
   683:           c       shift        r*8         main diagonal shift
   684:           c
   685:           c       output
   686:           c
   687:           c       a            r*8         array for nonzeros
   688:           c       colidx       i           col indices
   689:           c       rowstr       i           row pointers
   690:           c
   691:           c       workspace
   692:           c
   693:           c       iv, arow, acol i
   694:           c       v, aelt        r*8
   695:           c---------------------------------------------------------------------
   696:           
   697:                 integer          i, iouter, ivelt, nzv, nn1
   698:                 integer          ivc(nonzer+1)
   699:                 double precision vc(nonzer+1)
   700:           
   701:           c---------------------------------------------------------------------
   702:           c      nonzer is approximately  (int(sqrt(nnza /n)));
   703:           c---------------------------------------------------------------------
   704:           
   705:                 external          sparse, sprnvc, vecset
   706:           !$    integer           omp_get_num_threads, omp_get_thread_num
   707:           !$    external          omp_get_num_threads, omp_get_thread_num
   708:                 integer           work
   709:           
   710:           
   711:           c---------------------------------------------------------------------
   712:           c    nn1 is the smallest power of two not less than n
   713:           c---------------------------------------------------------------------
   714:           
   715:                 nn1 = 1
   716:            50   continue
   717:                   nn1 = 2 * nn1
   718:                   if (nn1 .lt. n) goto 50
   719:           
   720:           c---------------------------------------------------------------------
   721:           c  Generate nonzero positions and save for the use in sparse.
   722:           c---------------------------------------------------------------------
   723:                 num_threads = 1
   724:           !$    num_threads = omp_get_num_threads()
   725:                 myid = 0
   726:           !$    myid  = omp_get_thread_num()
   727:                 if (num_threads .gt. max_threads) then
   728:                    if (myid .eq. 0) write(*,100) num_threads, max_threads
   729:           100      format(' Warning: num_threads',i6,
   730:                &          ' exceeded an internal limit',i6)
   731:                    num_threads = max_threads
   732:                 endif
   733:                 work  = (n + num_threads - 1)/num_threads
   734:                 ilow  = work * myid + 1
   735:                 ihigh = ilow + work - 1
   736:                 if (ihigh .gt. n) ihigh = n
   737:           
   738: +------>        do iouter = 1, ihigh
   739: |                  nzv = nonzer
   740: |                  call sprnvc( n, nzv, nn1, vc, ivc )
   741: |                  if ( iouter .ge. ilow ) then
   742: |                     call vecset( n, vc, ivc, nzv, iouter, .5D0 )
   743: |                     arow(iouter) = nzv
   744: |V----->              do ivelt = 1, nzv
   745: ||                       acol(ivelt, iouter) = ivc(ivelt)
   746: ||                       aelt(ivelt, iouter) = vc(ivelt)
   747: |V-----               enddo
   748: |                  endif
   749: +------         enddo
   750:           !$omp barrier
   751:           
   752:           c---------------------------------------------------------------------
   753:           c       ... make the sparse matrix from list of elements with duplicates
   754:           c           (v and iv are used as  workspace)
   755:           c---------------------------------------------------------------------
   756:                 call sparse( a, colidx, rowstr, n, nz, nonzer, arow, acol,
   757:                >             aelt, firstrow, lastrow,
   758:                >             v, iv(1), iv(nz+1), rcond, shift )
   759:                 return
   760:           
   761:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: SPARSE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   805: vec( 101): Vectorized loop.
   809: vec( 103): Unvectorized loop.
   809: vec( 113): Overhead of loop division is too large.
   810: vec( 102): Partially vectorized loop.
   814: opt(1036): Potential feedback - use directive or compiler option if OK.
   814: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: ROWSTR
   823: vec( 101): Vectorized loop.
   824: opt(1037): Feedback of array elements.
   824: vec( 126): Idiom detected.: ITERATION
   827: par(1809): Barrier synchronization.
   831: vec( 101): Vectorized loop.
   832: vec( 126): Idiom detected.: SUM
   836: vec( 101): Vectorized loop.
   840: par(1809): Barrier synchronization.
   860: vec( 102): Partially vectorized loop.
   861: vec( 101): Vectorized loop.
   882: vec( 109): Vectorization obstructive statement.
   882: vec( 180): I/O statement obstructs vectorization.
   894: opt(1135): Outer loop conditionally executes inner loop.
   894: vec( 103): Unvectorized loop.
   894: vec( 108): Unvectorizable loop structure.
   899: vec( 101): Vectorized loop.
   928: par(1809): Barrier synchronization.
   934: vec( 101): Vectorized loop.
   935: opt(1037): Feedback of array elements.
   935: vec( 126): Idiom detected.: ITERATION
   938: par(1809): Barrier synchronization.
   942: vec( 101): Vectorized loop.
   943: vec( 126): Idiom detected.: SUM
   947: vec( 101): Vectorized loop.
   951: par(1809): Barrier synchronization.
   954: par(1803): Parallelized by "do".
   954: vec( 101): Vectorized loop.
   962: vec( 101): Vectorized loop.
   970: par(1803): Parallelized by "do".
   970: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: SPARSE
FORMAT LIST

 LINE   LOOP      STATEMENT

   766:                 subroutine sparse( a, colidx, rowstr, n, nz, nonzer, arow, acol,
   767:                >                   aelt, firstrow, lastrow,
   768:                >                   v, iv, nzloc, rcond, shift )
   769:           c---------------------------------------------------------------------
   770:           c---------------------------------------------------------------------
   771:           
   772:                 use tinfo
   773:           
   774:                 implicit           none
   775:           
   776:                 integer            colidx(*), rowstr(*), iv(*)
   777:                 integer            firstrow, lastrow
   778:                 integer            n, nz, nonzer, arow(*), acol(nonzer+1,*)
   779:                 double precision   a(*), aelt(nonzer+1,*), v(*), rcond, shift
   780:           
   781:           c---------------------------------------------------------------------
   782:           c       rows range from firstrow to lastrow
   783:           c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
   784:           c---------------------------------------------------------------------
   785:                 integer            nzloc(n), nrows
   786:           
   787:           c---------------------------------------------------
   788:           c       generate a sparse matrix from a list of
   789:           c       [col, row, element] tri
   790:           c---------------------------------------------------
   791:           
   792:                 integer            i, j, j1, j2, nza, k, kk, nzrow, jcol
   793:                 double precision   xi, size, scale, ratio, va
   794:           
   795:           c---------------------------------------------------------------------
   796:           c    how many rows of result
   797:           c---------------------------------------------------------------------
   798:                 nrows = lastrow - firstrow + 1
   799:                 j1 = ilow + 1
   800:                 j2 = ihigh + 1
   801:           
   802:           c---------------------------------------------------------------------
   803:           c     ...count the number of triples in each row
   804:           c---------------------------------------------------------------------
   805: V------>        do j = j1, j2
   806: |                  rowstr(j) = 0
   807: V------         enddo
   808:           
   809: +------>        do i = 1, n
   810: |S----->           do nza = 1, arow(i)
   811: ||                    j = acol(nza, i)
   812: ||                    if (j.ge.ilow .and. j.le.ihigh) then
   813: ||                       j = j + 1
   814: ||                       rowstr(j) = rowstr(j) + arow(i)
   815: ||                    endif
   816: |S-----            end do
   817: +------         end do
   818:           
   819:                 if (myid .eq. 0) then
   820:                    rowstr(1) = 1
   821:                    j1 = 1
   822:                 endif
   823: V------>        do j = j1+1, j2
   824: |                  rowstr(j) = rowstr(j) + rowstr(j-1)
   825: V------         enddo
   826:                 if (myid .lt. num_threads) last_n(myid) = rowstr(j2)
   827:           !$omp barrier
   828:           
   829:                 nzrow = 0
   830:                 if (myid .lt. num_threads) then
   831: V------>           do i = 0, myid-1
   832: |                     nzrow = nzrow + last_n(i)
   833: V------            end do
   834:                 endif
   835:                 if (nzrow .gt. 0) then
   836: V------>           do j = j1, j2
   837: |                     rowstr(j) = rowstr(j) + nzrow
   838: V------            enddo
   839:                 endif
   840:           !$omp barrier
   841:                 nza = rowstr(nrows+1) - 1
   842:           
   843:           c---------------------------------------------------------------------
   844:           c     ... rowstr(j) now is the location of the first nonzero
   845:           c           of row j of a
   846:           c---------------------------------------------------------------------
   847:           
   848:                 if (nza .gt. nz) then
   849:           !$omp master
   850:                    write(*,*) 'Space for matrix elements exceeded in sparse'
   851:                    write(*,*) 'nza, nzmax = ',nza, nz
   852:           !$omp end master
   853:                    stop
   854:                 endif
   855:           
   856:           
   857:           c---------------------------------------------------------------------
   858:           c     ... preload data pages
   859:           c---------------------------------------------------------------------
   860: S------>        do j = ilow, ihigh
   861: |V----->           do k = rowstr(j), rowstr(j+1)-1
   862: ||                     v(k) = 0.d0
   863: ||                     iv(k) = 0
   864: |V-----            enddo
   865: |                  nzloc(j) = 0
   866: S------         enddo
   867:           
   868:           c---------------------------------------------------------------------
   869:           c     ... generate actual values by summing duplicates
   870:           c---------------------------------------------------------------------
   871:           
   872:                 size = 1.0D0
   873:                 ratio = rcond ** (1.0D0 / dfloat(n))
   874:           
   875: +------>        do i = 1, n
   876: |+----->           do nza = 1, arow(i)
   877: ||                    j = acol(nza, i)
   878: ||        
   879: ||                    if (j .lt. ilow .or. j .gt. ihigh) goto 60
   880: ||        
   881: ||                    scale = size * aelt(nza, i)
   882: ||+---->              do nzrow = 1, arow(i)
   883: |||                      jcol = acol(nzrow, i)
   884: |||                      va = aelt(nzrow, i) * scale
   885: |||       
   886: |||       c---------------------------------------------------------------------
   887: |||       c       ... add the identity * rcond to the generated matrix to bound
   888: |||       c           the smallest eigenvalue from below by rcond
   889: |||       c---------------------------------------------------------------------
   890: |||                      if (jcol .eq. j .and. j .eq. i) then
   891: |||                         va = va + rcond - shift
   892: |||                      endif
   893: |||       
   894: |||+--->                 do k = rowstr(j), rowstr(j+1)-1
   895: ||||                        if (iv(k) .gt. jcol) then
   896: ||||      c---------------------------------------------------------------------
   897: ||||      c       ... insert colidx here orderly
   898: ||||      c---------------------------------------------------------------------
   899: ||||V-->                       do kk = rowstr(j+1)-2, k, -1
   900: |||||                             if (iv(kk) .gt. 0) then
   901: |||||                                v(kk+1)  = v(kk)
   902: |||||                                iv(kk+1) = iv(kk)
   903: |||||                             endif
   904: ||||V--                        enddo
   905: ||||                           iv(k) = jcol
   906: ||||                           v(k)  = 0.d0
   907: ||||                           goto 40
   908: ||||                        else if (iv(k) .eq. 0) then
   909: ||||                           iv(k) = jcol
   910: ||||                           goto 40
   911: ||||                        else if (iv(k) .eq. jcol) then
   912: ||||      c---------------------------------------------------------------------
   913: ||||      c       ... mark the duplicated entry
   914: ||||      c---------------------------------------------------------------------
   915: ||||                           nzloc(j) = nzloc(j) + 1
   916: ||||                           goto 40
   917: ||||                        endif
   918: |||+---                  enddo
   919: |||                      print *,'internal error in sparse: i=',i
   920: |||                      stop
   921: |||          40          continue
   922: |||                      v(k) = v(k) + va
   923: ||+----               enddo
   924: ||           60       continue
   925: |+-----            enddo
   926: |                  size = size * ratio
   927: +------         enddo
   928:           !$omp barrier
   929:           
   930:           
   931:           c---------------------------------------------------------------------
   932:           c       ... remove empty entries and generate final results
   933:           c---------------------------------------------------------------------
   934: V------>        do j = ilow+1, ihigh
   935: |                  nzloc(j) = nzloc(j) + nzloc(j-1)
   936: V------         enddo
   937:                 if (myid .lt. num_threads) last_n(myid) = nzloc(ihigh)
   938:           !$omp barrier
   939:           
   940:                 nzrow = 0
   941:                 if (myid .lt. num_threads) then
   942: V------>           do i = 0, myid-1
   943: |                     nzrow = nzrow + last_n(i)
   944: V------            end do
   945:                 endif
   946:                 if (nzrow .gt. 0) then
   947: V------>           do j = ilow, ihigh
   948: |                     nzloc(j) = nzloc(j) + nzrow
   949: V------            enddo
   950:                 endif
   951:           !$omp barrier
   952:           
   953:           !$omp do
   954: Y------>        do j = 1, nrows
   955: |                  if (j .gt. 1) then
   956: |                     j1 = rowstr(j) - nzloc(j-1)
   957: |                  else
   958: |                     j1 = 1
   959: |                  endif
   960: |                  j2 = rowstr(j+1) - nzloc(j) - 1
   961: |                  nza = rowstr(j)
   962: |V----->           do k = j1, j2
   963: ||                    a(k) = v(nza)
   964: ||                    colidx(k) = iv(nza)
   965: ||                    nza = nza + 1
   966: |V-----            enddo
   967: Y------         enddo
   968:           !$omp end do
   969:           !$omp do
   970: Y------>        do j = 2, nrows+1
   971: |                  rowstr(j) = rowstr(j) - nzloc(j-1)
   972: Y------         enddo
   973:           !$omp end do
   974:                 nza = rowstr(nrows+1) - 1
   975:           
   976:           
   977:           CC       write (*, 11000) nza
   978:                 return
   979:           11000   format ( //,'final nonzero count in sparse ',
   980:                1            /,'number of nonzeros       = ', i16 )
   981:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: SPRNVC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1018: vec( 108): Unvectorizable loop structure.
  1018: vec( 110): Vectorization obstructive procedure reference.: ICNVRT
  1018: vec( 110): Vectorization obstructive procedure reference.: RANDLC
  1032: vec( 101): Vectorized loop.
  1033: vec( 126): Idiom detected.: SEARCH


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: SPRNVC
FORMAT LIST

 LINE   LOOP      STATEMENT

   987:                 subroutine sprnvc( n, nz, nn1, v, iv )
   988:           c---------------------------------------------------------------------
   989:           c---------------------------------------------------------------------
   990:           
   991:                 use cg_data, only : amult, tran
   992:           
   993:                 implicit           none
   994:           
   995:                 double precision   v(*)
   996:                 integer            n, nz, nn1, iv(*)
   997:           
   998:           
   999:           c---------------------------------------------------------------------
  1000:           c       generate a sparse n-vector (v, iv)
  1001:           c       having nzv nonzeros
  1002:           c
  1003:           c       mark(i) is set to 1 if position i is nonzero.
  1004:           c       mark is all zero on entry and is reset to all zero before exit
  1005:           c       this corrects a performance bug found by John G. Lewis, caused by
  1006:           c       reinitialization of mark on every one of the n calls to sprnvc
  1007:           c---------------------------------------------------------------------
  1008:           
  1009:                   integer            nzv, ii, i, icnvrt
  1010:           
  1011:                   external           randlc, icnvrt
  1012:                   double precision   randlc, vecelt, vecloc
  1013:           
  1014:           
  1015:                   nzv = 0
  1016:           
  1017:           100     continue
  1018:                   if (nzv .ge. nz) goto 110
  1019:           
  1020:                    vecelt = randlc( tran, amult )
  1021:           
  1022:           c---------------------------------------------------------------------
  1023:           c   generate an integer between 1 and n in a portable manner
  1024:           c---------------------------------------------------------------------
  1025:                    vecloc = randlc(tran, amult)
  1026:                    i = icnvrt(vecloc, nn1) + 1
  1027:                    if (i .gt. n) goto 100
  1028:           
  1029:           c---------------------------------------------------------------------
  1030:           c  was this integer generated already?
  1031:           c---------------------------------------------------------------------
  1032: V------>           do ii = 1, nzv
  1033: |                     if (iv(ii) .eq. i) goto 100
  1034: V------            enddo
  1035:                    nzv = nzv + 1
  1036:                    v(nzv) = vecelt
  1037:                    iv(nzv) = i
  1038:                    goto 100
  1039:           110     continue
  1040:           
  1041:                 return
  1042:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: ICNVRT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: ICNVRT
FORMAT LIST

 LINE   LOOP      STATEMENT

  1048:                 function icnvrt(x, ipwr2)
  1049:           c---------------------------------------------------------------------
  1050:           c---------------------------------------------------------------------
  1051:           
  1052:                 implicit           none
  1053:           
  1054:                 double precision   x
  1055:                 integer            ipwr2, icnvrt
  1056:           
  1057:           c---------------------------------------------------------------------
  1058:           c    scale a double precision number x in (0,1) by a power of 2 and chop it
  1059:           c---------------------------------------------------------------------
  1060:                 icnvrt = int(ipwr2 * x)
  1061:           
  1062:                 return
  1063:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: VECSET
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1086: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:36 2020
FILE NAME: cg.f

PROCEDURE NAME: VECSET
FORMAT LIST

 LINE   LOOP      STATEMENT

  1069:                 subroutine vecset(n, v, iv, nzv, i, val)
  1070:           c---------------------------------------------------------------------
  1071:           c---------------------------------------------------------------------
  1072:           
  1073:                 implicit           none
  1074:           
  1075:                 integer            n, iv(*), nzv, i, k
  1076:                 double precision   v(*), val
  1077:           
  1078:           c---------------------------------------------------------------------
  1079:           c       set ith element of sparse vector (v, iv) with
  1080:           c       nzv nonzeros to val
  1081:           c---------------------------------------------------------------------
  1082:           
  1083:                 logical set
  1084:           
  1085:                 set = .false.
  1086: V------>        do k = 1, nzv
  1087: |                  if (iv(k) .eq. i) then
  1088: |                     v(k) = val
  1089: |                     set  = .true.
  1090: |                  endif
  1091: V------         enddo
  1092:                 if (.not. set) then
  1093:                    nzv     = nzv + 1
  1094:                    v(nzv)  = val
  1095:                    iv(nzv) = i
  1096:                 endif
  1097:                 return
  1098:                 end


