NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:37 2020
FILE NAME: ep.f

PROCEDURE NAME: EMBAR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   121: par(1801): Parallel routine generated.: MAIN_$1
   122: vec( 101): Vectorized loop.
   140: vec( 103): Unvectorized loop.
   140: vec( 110): Vectorization obstructive procedure reference.: RANDLC
   141: opt(1025): Reference to this procedure inhibits optimization.: RANDLC
   150: vec( 101): Vectorized loop.
   160: par(1801): Parallel routine generated.: MAIN_$2
   162: vec( 101): Vectorized loop.
   167: par(1803): Parallelized by "do".
   175: vec( 103): Unvectorized loop.
   175: vec( 110): Vectorization obstructive procedure reference.: RANDLC
   177: opt(1025): Reference to this procedure inhibits optimization.: RANDLC
   195: vec( 103): Unvectorized loop.
   204: opt(1036): Potential feedback - use directive or compiler option if OK.: QQ
   204: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: QQ
   215: vec( 103): Unvectorized loop.
   215: vec( 109): Vectorization obstructive statement.
   216: par(3007): Atomic-construct generated.
   219: par(1807): Critical-construct generated.
   219: par(1809): Barrier synchronization.
   221: vec( 101): Vectorized loop.
   222: vec( 126): Idiom detected.: SUM
   268: opt(1118): This I/O statement inhibits optimization of loop.
   268: vec( 180): I/O statement obstructs vectorization.
   268: vec( 103): Unvectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:37 2020
FILE NAME: ep.f

PROCEDURE NAME: EMBAR
FORMAT LIST

 LINE   LOOP      STATEMENT

    48:                 program EMBAR
    49:           c---------------------------------------------------------------------
    50:           
    51:           c   This is the OpenMP version of the APP Benchmark 1,
    52:           c   the "embarassingly parallel" benchmark.
    53:           
    54:                 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    55:           
    56:                 use ep_data
    57:           
    58:                 implicit none
    59:           
    60:                 double precision Mops, t1, t2, t3, t4, x1, x2,
    61:                >                 sx, sy, tm, an, tt, gc, dum(3)
    62:                 double precision sx_verify_value, sy_verify_value, sx_err, sy_err
    63:           
    64:                 integer          i, ik, kk, l, k, nit,
    65:                >                 np, k_offset, j
    66:           
    67:                 logical          verified, timers_enabled
    68:           
    69:                 external         randlc, timer_read
    70:                 double precision randlc, timer_read
    71:           
    72:                 character*15     size
    73:           
    74:                 double precision epsilon
    75:                 parameter       (epsilon=1.d-8)
    76:           
    77:           !$    integer  omp_get_max_threads
    78:           !$    external omp_get_max_threads
    79:           
    80:                 data             dum /1.d0, 1.d0, 1.d0/
    81:           
    82:           
    83:                 call check_timer_flag( timers_enabled )
    84:           
    85:           c   Because the size of the problem is too large to store in a 32-bit
    86:           c   integer for some classes, we put it into a string (for printing).
    87:           c   Have to strip off the decimal point put in there by the floating
    88:           c   point print statement (internal file)
    89:           
    90:                 write(*, 1000)
    91:                 write(size, '(f15.0)' ) 2.d0**(m+1)
    92:                 j = 15
    93:                 if (size(j:j) .eq. '.') j = j - 1
    94:                 write (*,1001) size(1:j)
    95:           !$    write (*,1003) omp_get_max_threads()
    96:                 write (*,*)
    97:           
    98:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
    99:                >          ' - EP Benchmark', /)
   100:            1001 format(' Number of random numbers generated: ', a15)
   101:            1003 format(' Number of available threads:        ', 2x,i13)
   102:           
   103:                 verified = .false.
   104:           
   105:           c   Compute the number of "batches" of random number pairs generated
   106:           c   per processor. Adjust if the number of processors does not evenly
   107:           c   divide the total number
   108:           
   109:                 np = nn
   110:           
   111:           
   112:           c   Call the random number generator functions and initialize
   113:           c   the x-array to reduce the effects of paging on the timings.
   114:           c   Also, call all mathematical functions that are used. Make
   115:           c   sure these initializations cannot be eliminated as dead code.
   116:           
   117:                 call vranlc(0, dum(1), dum(2), dum(3))
   118:                 dum(1) = randlc(dum(2), dum(3))
   119:                 Mops = log(sqrt(abs(max(1.d0,1.d0))))
   120:           
   121:           !$omp parallel default(shared) private(i)
   122: V------>        do 5    i = 1, 2*nk
   123: |                  x(i) = -1.d99
   124: V------    5    continue
   125:           
   126:                 call timer_clear(1)
   127:                 if (timers_enabled) call timer_clear(2)
   128:                 if (timers_enabled) call timer_clear(3)
   129:           !$omp end parallel
   130:           
   131:                 call timer_start(1)
   132:           
   133:                 t1 = a
   134:                 call vranlc(0, t1, a, x)
   135:           
   136:           c   Compute AN = A ^ (2 * NK) (mod 2^46).
   137:           
   138:                 t1 = a
   139:           
   140: +------>        do 100 i = 1, mk + 1
   141: |                  t2 = randlc(t1, t1)
   142: +------    100  continue
   143:           
   144:                 an = t1
   145:                 tt = s
   146:                 gc = 0.d0
   147:                 sx = 0.d0
   148:                 sy = 0.d0
   149:           
   150: V------>        do 110 i = 0, nq - 1
   151: |                  q(i) = 0.d0
   152: V------    110  continue
   153:           
   154:           c   Each instance of this loop may be performed independently. We compute
   155:           c   the k offsets separately to take into account the fact that some nodes
   156:           c   have more numbers to generate than others
   157:           
   158:                 k_offset = -1
   159:           
   160:           !$omp parallel default(shared) reduction(+:sx,sy)
   161:           !$omp&  private(k,kk,t1,t2,t3,t4,i,ik,x1,x2,l)
   162: V------>        do 115 i = 0, nq - 1
   163: |                  qq(i) = 0.d0
   164: V------    115  continue
   165:           
   166:           !$omp do schedule(static)
   167: P------>        do 150 k = 1, np
   168: |                  kk = k_offset + k
   169: |                  t1 = s
   170: |                  t2 = an
   171: |         
   172: |         c        Find starting seed t1 for this kk.
   173: |         
   174: |                  if (timers_enabled) call timer_start(3)
   175: |+----->           do 120 i = 1, 100
   176: ||                    ik = kk / 2
   177: ||                    if (2 * ik .ne. kk) t3 = randlc(t1, t2)
   178: ||                    if (ik .eq. 0) goto 130
   179: ||                    t3 = randlc(t2, t2)
   180: ||                    kk = ik
   181: |+-----    120     continue
   182: |         
   183: |         c        Compute uniform pseudorandom numbers.
   184: |          130     continue
   185: |         
   186: |                  call vranlc(2 * nk, t1, a, x)
   187: |                  if (timers_enabled) call timer_stop(3)
   188: |         
   189: |         c        Compute Gaussian deviates by acceptance-rejection method and
   190: |         c        tally counts in concentric square annuli.  This loop is not
   191: |         c        vectorizable.
   192: |         
   193: |                  if (timers_enabled) call timer_start(2)
   194: |         
   195: |+----->           do 140 i = 1, nk
   196: ||                    x1 = 2.d0 * x(2*i-1) - 1.d0
   197: ||                    x2 = 2.d0 * x(2*i) - 1.d0
   198: ||                    t1 = x1 ** 2 + x2 ** 2
   199: ||                    if (t1 .le. 1.d0) then
   200: ||                       t2   = sqrt(-2.d0 * log(t1) / t1)
   201: ||                       t3   = (x1 * t2)
   202: ||                       t4   = (x2 * t2)
   203: ||                       l    = max(abs(t3), abs(t4))
   204: ||                       qq(l) = qq(l) + 1.d0
   205: ||                       sx   = sx + t3
   206: ||                       sy   = sy + t4
   207: ||                    endif
   208: |+-----    140     continue
   209: |         
   210: |                  if (timers_enabled) call timer_stop(2)
   211: |         
   212: P------    150  continue
   213:           !$omp end do nowait
   214:           
   215: +------>        do 155 i = 0, nq - 1
   216: |         !$omp atomic
   217: |                  q(i) = q(i) + qq(i)
   218: +------    155  continue
   219:           !$omp end parallel
   220:           
   221: V------>        do 160 i = 0, nq - 1
   222: |                  gc = gc + q(i)
   223: V------    160  continue
   224:           
   225:                 call timer_stop(1)
   226:                 tm  = timer_read(1)
   227:           
   228:                 nit=0
   229:                 verified = .true.
   230:                 if (m.eq.24) then
   231:                    sx_verify_value = -3.247834652034740D+3
   232:                    sy_verify_value = -6.958407078382297D+3
   233:                 elseif (m.eq.25) then
   234:                    sx_verify_value = -2.863319731645753D+3
   235:                    sy_verify_value = -6.320053679109499D+3
   236:                 elseif (m.eq.28) then
   237:                    sx_verify_value = -4.295875165629892D+3
   238:                    sy_verify_value = -1.580732573678431D+4
   239:                 elseif (m.eq.30) then
   240:                    sx_verify_value =  4.033815542441498D+4
   241:                    sy_verify_value = -2.660669192809235D+4
   242:                 elseif (m.eq.32) then
   243:                    sx_verify_value =  4.764367927995374D+4
   244:                    sy_verify_value = -8.084072988043731D+4
   245:                 elseif (m.eq.36) then
   246:                    sx_verify_value =  1.982481200946593D+5
   247:                    sy_verify_value = -1.020596636361769D+5
   248:                 elseif (m.eq.40) then
   249:                    sx_verify_value = -5.319717441530D+05
   250:                    sy_verify_value = -3.688834557731D+05
   251:                 elseif (m.eq.44) then
   252:                    sx_verify_value = -1.986768071074039D+1
   253:                    sy_verify_value = -2.164441190520301D+1
   254:                 else
   255:                    verified = .false.
   256:                 endif
   257:                 if (verified) then
   258:                    sx_err = abs((sx - sx_verify_value)/sx_verify_value)
   259:                    sy_err = abs((sy - sy_verify_value)/sy_verify_value)
   260:                    if (ieee_is_nan(sx_err) .or. ieee_is_nan(sy_err)) then
   261:                       verified = .false.
   262:                    else
   263:                       verified = ((sx_err.le.epsilon) .and. (sy_err.le.epsilon))
   264:                    endif
   265:                 endif
   266:                 Mops = 2.d0**(m+1)/tm/1000000.d0
   267:           
   268: +------>        write (6,11) tm, m, gc, sx, sy, (i, q(i), i = 0, nq - 1)
   269:            11   format ('EP Benchmark Results:'//'CPU Time =',f10.4/'N = 2^',
   270:                >        i5/'No. Gaussian Pairs =',f15.0/'Sums = ',1p,2d25.15/
   271:                >        'Counts:'/(i3,0p,f15.0))
   272:           
   273:                 call print_results('EP', class, m+1, 0, 0, nit,
   274:                >                   tm, Mops,
   275:                >                   'Random numbers generated',
   276:                >                   verified, npbversion, compiletime, cs1,
   277:                >                   cs2, cs3, cs4, cs5, cs6, cs7)
   278:           
   279:           
   280:                 if (timers_enabled) then
   281:                    if (tm .le. 0.d0) tm = 1.0
   282:                  	 write(*,800)
   283:           800      format('  SECTION   Time (secs)')
   284:                    tt = timer_read(1)
   285:                    print 810, 'Total ', tt, tt*100./tm
   286:                    tt = timer_read(2)
   287:                    print 810, 'Gausd_pairs ', tt, tt*100./tm
   288:                    tt = timer_read(3)
   289:                    print 810, 'Random_n ', tt, tt*100./tm
   290:           810      format(2x,a12,':',f9.3,'  ( ',f6.2,' %)')
   291:                 endif
   292:           
   293:           
   294:                 end


