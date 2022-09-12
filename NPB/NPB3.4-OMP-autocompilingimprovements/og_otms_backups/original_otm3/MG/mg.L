NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME : mg.f

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fassume-contiguous                      : disable
        -fcopyin-intent-out                      : enable
        -fcse-after-vectorization                : disable
        -ffast-formatted-io                      : enable
        -ffast-math                              : enable
        -fignore-asynchronous                    : disable
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
        -fivdep                                  : disable
        -floop-collapse                          : enable
        -floop-count                             : 5000
        -floop-fusion                            : enable
        -floop-interchange                       : enable
        -floop-normalize                         : enable
        -floop-split                             : disable
        -floop-strip-mine                        : enable
        -floop-unroll                            : enable
        -floop-unroll-complete                   : 4
        -floop-unroll-complete-nest              : 3
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : disable
        -fmove-nested-loop-invariants-outer      : enable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fnamed-noalias-aggressive               : disable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freorder-logical-expression             : enable
        -freplace-loop-equation                  : disable
        -marray-io                               : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mstack-arrays                           : enable
        -mno-stack-arrays                        : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : enable
        -mvector-loop-count-test                 : disable
        -mvector-low-precise-divide-function     : disable
        -mvector-merge-conditional               : enable
        -mvector-packed                          : disable
        -mvector-power-to-explog                 : disable
        -mvector-power-to-sqrt                   : enable
        -mvector-reduction                       : enable
        -mvector-shortloop-reduction             : disable
        -mvector-sqrt-instruction                : disable
        -mvector-threshold                       : 5
        -mwork-vector-kind=none                  : disable

      Parallelization Options : 
        -fopenmp                                 : enable
        -mparallel                               : disable
        -mparallel-innerloop                     : disable
        -mparallel-omp-routine                   : enable
        -mparallel-outerloop-strip-mine          : disable
        -mparallel-sections                      : disable
        -mparallel-threshold                     : 2000
        -mschedule-dynamic                       : disable
        -mschedule-runtime                       : disable
        -mschedule-static                        : disable
        -mschedule-chunk-size                    : 0
        -pthread                                 : disable

      Inlining Options : 
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
        -finline-max-depth                       : 2
        -finline-max-function-size               : 50
        -finline-max-times                       : 6

      Code Generation Options : 
        -finstrument-functions                   : disable
        -fpic | -fPIC                            : disable
        -ftrace                                  : disable
        -p                                       : disable

      Debugging Options : 
        -fbounds-check                           : disable
        -fcheck=all                              : disable
        -fcheck=alias                            : disable
        -fcheck=bits                             : disable
        -fcheck=bounds                           : disable
        -fcheck=dangling                         : disable
        -fcheck=do                               : disable
        -fcheck=iovf                             : disable
        -fcheck=pointer                          : disable
        -fcheck=present                          : disable
        -fcheck=recursion                        : disable
        -g                                       : disable
        -minit-stack                             : not specified
        -mmemory-trace                           : disable
        -mmemory-trace-full                      : disable
        -traceback                               : disable

      Language Options : 
        -bss                                     : disable
        -fdefault-integer                        : 4
        -fdefault-double                         : 8
        -fdefault-real                           : 4
        -fextend-source                          : disable
        -ffree-form                              : disable
        -ffixed-form                             : enable
        -ff90-sign                               : disable
        -fmax-continuation-lines                 : 1023
        -frealloc-lhs                            : enable
        -frealloc-lhs-array                      : enable
        -frealloc-lhs-scalar                     : enable
        -save                                    : disable
        -std                                     : Fortran 2008

      Message Options : 
        -Wall                                    : disable
        -Werror                                  : disable
        -Wextension                              : enable
        -Wobsolescent                            : enable
        -Woverflow                               : disable
        -Woverflow-errors                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -pedantic-errors                         : disable
        -w                                       : disable

      List Output Options : 
        -report-file                             : not specified
        -report-append-mode                      : disable
        -report-all                              : enable
        -report-cg                               : enable
        -report-diagnostics                      : enable
        -report-format                           : enable
        -report-inline                           : enable
        -report-option                           : enable
        -report-vector                           : enable

      Directory Options : 
        -fintrinsic-modules-path                 : /opt/nec/ve/nfort/3.2.0/include
        -module                                  : not specified

      Miscellaneous Options : 
        -v                                       : disable


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    83: vec( 103): Unvectorized loop.
    83: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   115: opt(1118): This I/O statement inhibits optimization of loop.
   115: vec( 180): I/O statement obstructs vectorization.
   115: vec( 103): Unvectorized loop.
   124: vec( 101): Vectorized loop.
   234: vec( 103): Unvectorized loop.
   234: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   247: opt(3008): Reference within a conditional branch moved outside loop.
   247: vec( 110): Vectorization obstructive procedure reference.: MG3P
   247: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
   247: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
   247: vec( 110): Vectorization obstructive procedure reference.: RESID
   247: vec( 180): I/O statement obstructs vectorization.
   247: vec( 103): Unvectorized loop.
   249: opt(1118): This I/O statement inhibits optimization of loop.
   347: opt(3008): Reference within a conditional branch moved outside loop.
   347: vec( 103): Unvectorized loop.
   347: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   347: vec( 180): I/O statement obstructs vectorization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

    51:                 program mg
    52:           c---------------------------------------------------------------------
    53:           
    54:                 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    55:           
    56:                 use mg_data
    57:                 use mg_fields
    58:           
    59:                 implicit none
    60:           
    61:           c---------------------------------------------------------------------------c
    62:           c k is the current level. It is passed down through subroutine args
    63:           c and is NOT global. it is the current iteration
    64:           c---------------------------------------------------------------------------c
    65:           
    66:                 integer k, it
    67:            
    68:                 external timer_read
    69:                 double precision t, tinit, mflops, timer_read
    70:           
    71:                 double precision rnm2, rnmu, old2, oldu, epsilon
    72:                 integer n1, n2, n3, nit
    73:                 double precision nn, verify_value, err
    74:                 logical verified
    75:           
    76:                 integer i, fstatus
    77:                 character t_names(t_last)*8
    78:                 double precision tmax
    79:           !$    integer  omp_get_max_threads
    80:           !$    external omp_get_max_threads
    81:           
    82:           
    83: +------>        do i = T_init, T_last
    84: |                  call timer_clear(i)
    85: +------         end do
    86:           
    87:                 call timer_start(T_init)
    88:           
    89:           c---------------------------------------------------------------------
    90:           c Read in input data
    91:           c---------------------------------------------------------------------
    92:           
    93:                 call check_timer_flag( timeron )
    94:                 if (timeron) then
    95:                    t_names(t_init) = 'init'
    96:                    t_names(t_bench) = 'benchmk'
    97:                    t_names(t_mg3P) = 'mg3P'
    98:                    t_names(t_psinv) = 'psinv'
    99:                    t_names(t_resid) = 'resid'
   100:                    t_names(t_rprj3) = 'rprj3'
   101:                    t_names(t_interp) = 'interp'
   102:                    t_names(t_norm2) = 'norm2'
   103:                    t_names(t_comm3) = 'comm3'
   104:                 endif
   105:           
   106:                 write (*, 1000)
   107:           
   108:                 open(unit=7,file="mg.input", status="old", iostat=fstatus)
   109:                 if (fstatus .eq. 0) then
   110:                    write(*,50)
   111:            50      format(' Reading from input file mg.input')
   112:                    read(7,*) lt
   113:                    read(7,*) nx(lt), ny(lt), nz(lt)
   114:                    read(7,*) nit
   115: +------>           read(7,*) (debug_vec(i),i=0,7)
   116:                 else
   117:                    write(*,51)
   118:            51      format(' No input file. Using compiled defaults ')
   119:                    lt = lt_default
   120:                    nit = nit_default
   121:                    nx(lt) = nx_default
   122:                    ny(lt) = ny_default
   123:                    nz(lt) = nz_default
   124: V------>           do i = 0,7
   125: |                     debug_vec(i) = debug_default
   126: V------            end do
   127:                 endif
   128:           
   129:           
   130:                 if ( (nx(lt) .ne. ny(lt)) .or. (nx(lt) .ne. nz(lt)) ) then
   131:                    Class = 'U'
   132:                 else if( nx(lt) .eq. 32 .and. nit .eq. 4 ) then
   133:                    Class = 'S'
   134:                 else if( nx(lt) .eq. 128 .and. nit .eq. 4 ) then
   135:                    Class = 'W'
   136:                 else if( nx(lt) .eq. 256 .and. nit .eq. 4 ) then
   137:                    Class = 'A'
   138:                 else if( nx(lt) .eq. 256 .and. nit .eq. 20 ) then
   139:                    Class = 'B'
   140:                 else if( nx(lt) .eq. 512 .and. nit .eq. 20 ) then
   141:                    Class = 'C'
   142:                 else if( nx(lt) .eq. 1024 .and. nit .eq. 50 ) then
   143:                    Class = 'D'
   144:                 else if( nx(lt) .eq. 2048 .and. nit .eq. 50 ) then
   145:                    Class = 'E'
   146:                 else if( nx(lt) .eq. 4096 .and. nit .eq. 50 ) then
   147:                    Class = 'F'
   148:                 else
   149:                    Class = 'U'
   150:                 endif
   151:           
   152:           c---------------------------------------------------------------------
   153:           c  Use these for debug info:
   154:           c---------------------------------------------------------------------
   155:           c     debug_vec(0) = 1 !=> report all norms
   156:           c     debug_vec(1) = 1 !=> some setup information
   157:           c     debug_vec(1) = 2 !=> more setup information
   158:           c     debug_vec(2) = k => at level k or below, show result of resid
   159:           c     debug_vec(3) = k => at level k or below, show result of psinv
   160:           c     debug_vec(4) = k => at level k or below, show result of rprj
   161:           c     debug_vec(5) = k => at level k or below, show result of interp
   162:           c     debug_vec(6) = 1 => (unused)
   163:           c     debug_vec(7) = 1 => (unused)
   164:           c---------------------------------------------------------------------
   165:                 a(0) = -8.0D0/3.0D0
   166:                 a(1) =  0.0D0
   167:                 a(2) =  1.0D0/6.0D0
   168:                 a(3) =  1.0D0/12.0D0
   169:            
   170:                 if(Class .eq. 'A' .or. Class .eq. 'S'.or. Class .eq.'W') then
   171:           c---------------------------------------------------------------------
   172:           c     Coefficients for the S(a) smoother
   173:           c---------------------------------------------------------------------
   174:                    c(0) =  -3.0D0/8.0D0
   175:                    c(1) =  +1.0D0/32.0D0
   176:                    c(2) =  -1.0D0/64.0D0
   177:                    c(3) =   0.0D0
   178:                 else
   179:           c---------------------------------------------------------------------
   180:           c     Coefficients for the S(b) smoother
   181:           c---------------------------------------------------------------------
   182:                    c(0) =  -3.0D0/17.0D0
   183:                    c(1) =  +1.0D0/33.0D0
   184:                    c(2) =  -1.0D0/61.0D0
   185:                    c(3) =   0.0D0
   186:                 endif
   187:                 lb = 1
   188:                 k  = lt
   189:           
   190:                 call alloc_space
   191:           
   192:                 call setup(n1,n2,n3,k)
   193:                 call zero3(u,n1,n2,n3)
   194:                 call zran3(v,n1,n2,n3,nx(lt),ny(lt),k)
   195:           
   196:                 call norm2u3(v,n1,n2,n3,rnm2,rnmu,nx(lt),ny(lt),nz(lt))
   197:           c     write(*,*)
   198:           c     write(*,*)' norms of random v are'
   199:           c     write(*,600) 0, rnm2, rnmu
   200:           c     write(*,*)' about to evaluate resid, k=',k
   201:           
   202:                 write (*, 1001) nx(lt),ny(lt),nz(lt), Class
   203:                 write (*, 1002) nit
   204:           !$    write (*, 1003) omp_get_max_threads()
   205:                 write (*, *)
   206:           
   207:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
   208:                >          ' - MG Benchmark', /)
   209:            1001 format(' Size: ', i4, 'x', i4, 'x', i4, '  (class ', A, ')' )
   210:            1002 format(' Iterations:                  ', i5)
   211:            1003 format(' Number of available threads: ', i5)
   212:           
   213:           
   214:                 call resid(u,v,r,n1,n2,n3,a,k)
   215:                 call norm2u3(r,n1,n2,n3,rnm2,rnmu,nx(lt),ny(lt),nz(lt))
   216:                 old2 = rnm2
   217:                 oldu = rnmu
   218:           
   219:           c---------------------------------------------------------------------
   220:           c     One iteration for startup
   221:           c---------------------------------------------------------------------
   222:                 call mg3P(u,v,r,a,c,n1,n2,n3,k)
   223:                 call resid(u,v,r,n1,n2,n3,a,k)
   224:                 call setup(n1,n2,n3,k)
   225:                 call zero3(u,n1,n2,n3)
   226:                 call zran3(v,n1,n2,n3,nx(lt),ny(lt),k)
   227:           
   228:                 call timer_stop(T_init)
   229:                 tinit = timer_read(T_init)
   230:           
   231:                 write( *,'(A,F15.3,A/)' )
   232:                >     ' Initialization time: ',tinit, ' seconds'
   233:           
   234: +------>        do i = T_bench, T_last
   235: |                  call timer_clear(i)
   236: +------         end do
   237:           
   238:                 call timer_start(T_bench)
   239:           
   240:                 if (timeron) call timer_start(T_resid2)
   241:                 call resid(u,v,r,n1,n2,n3,a,k)
   242:                 if (timeron) call timer_stop(T_resid2)
   243:                 call norm2u3(r,n1,n2,n3,rnm2,rnmu,nx(lt),ny(lt),nz(lt))
   244:                 old2 = rnm2
   245:                 oldu = rnmu
   246:           
   247: +------>        do  it=1,nit
   248: |                  if (it.eq.1 .or. it.eq.nit .or. mod(it,5).eq.0) then
   249: |                     write(*,80) it
   250: |            80       format('  iter ',i3)
   251: |                  endif
   252: |                  if (timeron) call timer_start(T_mg3P)
   253: |                  call mg3P(u,v,r,a,c,n1,n2,n3,k)
   254: |                  if (timeron) call timer_stop(T_mg3P)
   255: |                  if (timeron) call timer_start(T_resid2)
   256: |                  call resid(u,v,r,n1,n2,n3,a,k)
   257: |                  if (timeron) call timer_stop(T_resid2)
   258: +------         enddo
   259:           
   260:           
   261:                 call norm2u3(r,n1,n2,n3,rnm2,rnmu,nx(lt),ny(lt),nz(lt))
   262:           
   263:                 call timer_stop(T_bench)
   264:           
   265:                 t = timer_read(T_bench)
   266:           
   267:                 verified = .FALSE.
   268:                 verify_value = 0.0
   269:           
   270:                 write(*,100)
   271:            100  format(/' Benchmark completed ')
   272:           
   273:                 epsilon = 1.d-8
   274:                 if (Class .ne. 'U') then
   275:                    if(Class.eq.'S') then
   276:                       verify_value = 0.5307707005734d-04
   277:                    elseif(Class.eq.'W') then
   278:                       verify_value = 0.6467329375339d-05
   279:                    elseif(Class.eq.'A') then
   280:                       verify_value = 0.2433365309069d-05
   281:                    elseif(Class.eq.'B') then
   282:                       verify_value = 0.1800564401355d-05
   283:                    elseif(Class.eq.'C') then
   284:                       verify_value = 0.5706732285740d-06
   285:                    elseif(Class.eq.'D') then
   286:                       verify_value = 0.1583275060440d-09
   287:                    elseif(Class.eq.'E') then
   288:                       verify_value = 0.5630442584711d-10
   289:                    elseif(Class.eq.'F') then
   290:                       verify_value = 0.1889225697989d-10
   291:                    endif
   292:           
   293:                    err = abs( rnm2 - verify_value ) / verify_value
   294:                    if( (.not.ieee_is_nan(err)) .and. (err .le. epsilon) ) then
   295:                       verified = .TRUE.
   296:                       write(*, 200)
   297:                       write(*, 201) rnm2
   298:                       write(*, 202) err
   299:            200        format(' VERIFICATION SUCCESSFUL ')
   300:            201        format(' L2 Norm is ', E20.13)
   301:            202        format(' Error is   ', E20.13)
   302:                    else
   303:                       verified = .FALSE.
   304:                       write(*, 300)
   305:                       write(*, 301) rnm2
   306:                       write(*, 302) verify_value
   307:            300        format(' VERIFICATION FAILED')
   308:            301        format(' L2 Norm is             ', E20.13)
   309:            302        format(' The correct L2 Norm is ', E20.13)
   310:                    endif
   311:                 else
   312:                    verified = .FALSE.
   313:                    write (*, 400)
   314:                    write (*, 401)
   315:                    write (*, 201) rnm2
   316:            400     format(' Problem size unknown')
   317:            401     format(' NO VERIFICATION PERFORMED')
   318:                 endif
   319:           
   320:                 nn = 1.0d0*nx(lt)*ny(lt)*nz(lt)
   321:           
   322:                 if( t .ne. 0. ) then
   323:                    mflops = 58.*nit*nn*1.0D-6 /t
   324:                 else
   325:                    mflops = 0.0
   326:                 endif
   327:           
   328:                 call print_results('MG', class, nx(lt), ny(lt), nz(lt),
   329:                >                   nit, t,
   330:                >                   mflops, '          floating point',
   331:                >                   verified, npbversion, compiletime,
   332:                >                   cs1, cs2, cs3, cs4, cs5, cs6, cs7)
   333:           
   334:           
   335:            600  format( i4, 2e19.12)
   336:           
   337:           c---------------------------------------------------------------------
   338:           c      More timers
   339:           c---------------------------------------------------------------------
   340:                 if (.not.timeron) goto 999
   341:           
   342:                 tmax = timer_read(t_bench)
   343:                 if (tmax .eq. 0.0) tmax = 1.0
   344:           
   345:                 write(*,800)
   346:            800  format('  SECTION   Time (secs)')
   347: +------>        do i=t_bench, t_last
   348: |                  t = timer_read(i)
   349: |                  if (i.eq.t_resid2) then
   350: |                     t = timer_read(T_resid) - t
   351: |                     write(*,820) 'mg-resid ', t, t*100./tmax
   352: |                  else
   353: |                     write(*,810) t_names(i), t, t*100./tmax
   354: |                  endif
   355: |          810     format(2x,a15,':',f9.3,'  ( ',f6.2,' %)')
   356: |          820     format(2x,a15,':',f9.3,'  ( ',f6.2,' %)')
   357: +------         end do
   358:           
   359:            999  continue
   360:           
   361:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:83)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (mg.f:83)
  LOOP END

  LOOP BEGIN: (mg.f:124)
    <Vectorized loop.>
    **  Runtime-tested short loop. (mg.f:124)
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:124)
    *** The number of VLOAD, VSTORE. :  0,  1. (mg.f:124)
  LOOP END

  LOOP BEGIN: (mg.f:115)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (mg.f:115)
  LOOP END

  LOOP BEGIN: (mg.f:234)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (mg.f:234)
  LOOP END

  LOOP BEGIN: (mg.f:247)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (mg.f:247)
    *** Vectorization obstructive procedure reference. : MG3P (mg.f:247)
    *** Vectorization obstructive procedure reference. : TIMER_START (mg.f:247)
    *** Vectorization obstructive procedure reference. : TIMER_STOP (mg.f:247)
    *** Vectorization obstructive procedure reference. : RESID (mg.f:247)
  LOOP END

  LOOP BEGIN: (mg.f:347)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_READ (mg.f:347)
    *** I/O statement obstructs vectorization. (mg.f:347)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 58 [s0-s12 s15-s16 s18-s44 s48-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 2512 bytes
      Register spill area      :  424 bytes
      Parameter area           :  248 bytes
      Register save area       :  176 bytes
      User data area           : 1656 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:83)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (mg.f:115)
    *** Estimated execution cycle                       : 9
  LOOP END

  LOOP BEGIN: (mg.f:234)
    *** Estimated execution cycle                       : 10
  LOOP END

  LOOP BEGIN: (mg.f:247)
    *** Estimated execution cycle                       : 56
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 9
            Across calls                                : 7
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 37
            Across calls                                :  8
            Over basic blocks                           : 29
    *** The number of SCALAR REGISTER TRANSFER          : 30
  LOOP END

  LOOP BEGIN: (mg.f:347)
    *** Estimated execution cycle                       : 85
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Across calls                                : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 8
            Across calls                                : 7
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 22
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SETUP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   384: opt(1418): Constant-length loop is expanded.
   385: vec( 103): Unvectorized loop.
   386: opt(1037): Feedback of array elements.
   386: vec( 120): Unvectorizable dependency.: NG
   390: vec( 101): Vectorized loop.
   396: vec( 101): Vectorized loop.
   397: opt(1418): Constant-length loop is expanded.
   420: vec( 101): Vectorized loop.
   421: opt(1037): Feedback of array elements.
   421: vec( 126): Idiom detected.: ITERATION


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SETUP
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   366:                 subroutine setup(n1,n2,n3,k)
   367:           
   368:           c---------------------------------------------------------------------
   369:           c---------------------------------------------------------------------
   370:           
   371:                 use mg_data
   372:                 implicit none
   373:           
   374:                 integer n1,n2,n3,k
   375:                 integer j
   376:           
   377:                 integer ax, mi(3,maxlevel)
   378:                 integer ng(3,maxlevel)
   379:           
   380:           
   381:                 ng(1,lt) = nx(lt)
   382:                 ng(2,lt) = ny(lt)
   383:                 ng(3,lt) = nz(lt)
   384: *------>        do  ax=1,3
   385: |+----->           do  k=lt-1,1,-1
   386: ||                    ng(ax,k) = ng(ax,k+1)/2
   387: |+-----            enddo
   388: *------         enddo
   389:            61   format(10i4)
   390: V------>        do  k=lt,1,-1
   391: |                  nx(k) = ng(1,k)
   392: |                  ny(k) = ng(2,k)
   393: |                  nz(k) = ng(3,k)
   394: V------         enddo
   395:           
   396: V------>        do  k = lt,1,-1
   397: |*----->           do  ax = 1,3
   398: ||                    mi(ax,k) = 2 + ng(ax,k)
   399: |*-----            enddo
   400: |         
   401: |                  m1(k) = mi(1,k)
   402: |                  m2(k) = mi(2,k)
   403: |                  m3(k) = mi(3,k)
   404: |         
   405: V------         enddo
   406:           
   407:                 k = lt
   408:                 is1 = 2 + ng(1,k) - ng(1,lt)
   409:                 ie1 = 1 + ng(1,k)
   410:                 n1 = 3 + ie1 - is1
   411:                 is2 = 2 + ng(2,k) - ng(2,lt)
   412:                 ie2 = 1 + ng(2,k)
   413:                 n2 = 3 + ie2 - is2
   414:                 is3 = 2 + ng(3,k) - ng(3,lt)
   415:                 ie3 = 1 + ng(3,k)
   416:                 n3 = 3 + ie3 - is3
   417:           
   418:           
   419:                 ir(lt)=1
   420: V------>        do  j = lt-1, 1, -1
   421: |                  ir(j)=ir(j+1)+one*m1(j+1)*m2(j+1)*m3(j+1)
   422: V------         enddo
   423:           
   424:           
   425:                 if( debug_vec(1) .ge. 1 )then
   426:                    write(*,*)' in setup, '
   427:                    write(*,*)' k  lt  nx  ny  nz ',
   428:                >        ' n1  n2  n3 is1 is2 is3 ie1 ie2 ie3'
   429:                    write(*,9) k,lt,ng(1,k),ng(2,k),ng(3,k),
   430:                >              n1,n2,n3,is1,is2,is3,ie1,ie2,ie3
   431:            9       format(15i4)
   432:                 endif
   433:           
   434:                 k = lt
   435:           
   436:                 return
   437:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SETUP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:385)
    <Unvectorized loop.>
    *** Unvectorizable dependency. : NG (mg.f:386)
  LOOP END

  LOOP BEGIN: (mg.f:385)
    <Unvectorized loop.>
    *** Unvectorizable dependency. : NG (mg.f:386)
  LOOP END

  LOOP BEGIN: (mg.f:385)
    <Unvectorized loop.>
    *** Unvectorizable dependency. : NG (mg.f:386)
  LOOP END

  LOOP BEGIN: (mg.f:390)
    <Vectorized loop.>
    **  Runtime-tested short loop. (mg.f:390)
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:390)
    *** The number of VLOAD, VSTORE. :  3,  3. (mg.f:390)
  LOOP END

  LOOP BEGIN: (mg.f:396)
    <Vectorized loop.>
    **  Runtime-tested short loop. (mg.f:396)
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:396)
    *** The number of VLOAD, VSTORE. :  3,  6. (mg.f:396)
  LOOP END

  LOOP BEGIN: (mg.f:420)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:420)
    *** The number of VLOAD, VSTORE. :  3,  1. (mg.f:420)
    *** Idiom detected. : ITERATION (mg.f:421)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SETUP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 51 [s0-s12 s15-s16 s18-s20 s23-s34 s43-s63]
      Vector registers         : 17 [v47-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1584 bytes
      Register spill area      :   32 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1296 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:385)
    *** Estimated execution cycle                       : 17
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END

  LOOP BEGIN: (mg.f:385)
    *** Estimated execution cycle                       : 17
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END

  LOOP BEGIN: (mg.f:385)
    *** Estimated execution cycle                       : 17
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END

  LOOP BEGIN: (mg.f:420)
    *** Estimated execution cycle                       : 508
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG3P
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   465: vec( 110): Vectorization obstructive procedure reference.: RPRJ3
   465: vec( 103): Unvectorized loop.
   467: opt(1025): Reference to this procedure inhibits optimization.: RPRJ3
   478: vec( 110): Vectorization obstructive procedure reference.: ZERO3
   478: vec( 110): Vectorization obstructive procedure reference.: PSINV
   478: vec( 110): Vectorization obstructive procedure reference.: RESID
   478: vec( 110): Vectorization obstructive procedure reference.: INTERP
   478: vec( 103): Unvectorized loop.
   483: opt(1025): Reference to this procedure inhibits optimization.: ZERO3


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG3P
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   442:                 subroutine mg3P(u,v,r,a,c,n1,n2,n3,k)
   443:           
   444:           c---------------------------------------------------------------------
   445:           c---------------------------------------------------------------------
   446:           
   447:           c---------------------------------------------------------------------
   448:           c     multigrid V-cycle routine
   449:           c---------------------------------------------------------------------
   450:           
   451:                 use mg_data
   452:                 implicit none
   453:           
   454:                 integer n1, n2, n3, k
   455:                 double precision u(nr),v(nv),r(nr)
   456:                 double precision a(0:3),c(0:3)
   457:           
   458:                 integer j
   459:           
   460:           c---------------------------------------------------------------------
   461:           c     down cycle.
   462:           c     restrict the residual from the find grid to the coarse
   463:           c---------------------------------------------------------------------
   464:           
   465: +------>        do  k= lt, lb+1 , -1
   466: |                  j = k-1
   467: |                  call rprj3(r(ir(k)),m1(k),m2(k),m3(k),
   468: |              >        r(ir(j)),m1(j),m2(j),m3(j),k)
   469: +------         enddo
   470:           
   471:                 k = lb
   472:           c---------------------------------------------------------------------
   473:           c     compute an approximate solution on the coarsest grid
   474:           c---------------------------------------------------------------------
   475:                 call zero3(u(ir(k)),m1(k),m2(k),m3(k))
   476:                 call psinv(r(ir(k)),u(ir(k)),m1(k),m2(k),m3(k),c,k)
   477:           
   478: +------>        do  k = lb+1, lt-1
   479: |                   j = k-1
   480: |         c---------------------------------------------------------------------
   481: |         c        prolongate from level k-1  to k
   482: |         c---------------------------------------------------------------------
   483: |                  call zero3(u(ir(k)),m1(k),m2(k),m3(k))
   484: |                  call interp(u(ir(j)),m1(j),m2(j),m3(j),
   485: |              >               u(ir(k)),m1(k),m2(k),m3(k),k)
   486: |         c---------------------------------------------------------------------
   487: |         c        compute residual for level k
   488: |         c---------------------------------------------------------------------
   489: |                  call resid(u(ir(k)),r(ir(k)),r(ir(k)),m1(k),m2(k),m3(k),a,k)
   490: |         c---------------------------------------------------------------------
   491: |         c        apply smoother
   492: |         c---------------------------------------------------------------------
   493: |                  call psinv(r(ir(k)),u(ir(k)),m1(k),m2(k),m3(k),c,k)
   494: +------         enddo
   495:            200  continue
   496:                 j = lt - 1
   497:                 k = lt
   498:                 call interp(u(ir(j)),m1(j),m2(j),m3(j),u,n1,n2,n3,k)
   499:                 call resid(u,v,r,n1,n2,n3,a,k)
   500:                 call psinv(r,u,n1,n2,n3,c,k)
   501:           
   502:                 return
   503:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG3P
VECTORIZATION LIST

  LOOP BEGIN: (mg.f:465)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : RPRJ3 (mg.f:465)
  LOOP END

  LOOP BEGIN: (mg.f:478)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : ZERO3 (mg.f:478)
    *** Vectorization obstructive procedure reference. : PSINV (mg.f:478)
    *** Vectorization obstructive procedure reference. : RESID (mg.f:478)
    *** Vectorization obstructive procedure reference. : INTERP (mg.f:478)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: MG3P
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 52 [s0-s12 s15-s16 s18-s39 s49-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 608 bytes
      Register spill area      : 176 bytes
      Parameter area           :  72 bytes
      Register save area       : 176 bytes
      User data area           : 176 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:465)
    *** Estimated execution cycle                       : 30
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 14
            Over basic blocks                           : 14
    *** The number of SCALAR REGISTER TRANSFER          : 11
  LOOP END

  LOOP BEGIN: (mg.f:478)
    *** Estimated execution cycle                       : 86
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 5
            Across calls                                : 5
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 10
            Across calls                                :  5
            Over basic blocks                           :  5
    *** The number of SCALAR REGISTER TRANSFER          : 24
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: PSINV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   535: par(1801): Parallel routine generated.: PSINV$1
   536: par(1803): Parallelized by "do".
   536: opt(2495): Loop nest collapsed by "collapse" clause.
   536: vec( 101): Vectorized loop.
   538: vec( 101): Vectorized loop.
   544: vec( 101): Vectorized loop.
   545: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: PSINV
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   508:                 subroutine psinv( r,u,n1,n2,n3,c,k)
   509:           
   510:           c---------------------------------------------------------------------
   511:           c---------------------------------------------------------------------
   512:           
   513:           c---------------------------------------------------------------------
   514:           c     psinv applies an approximate inverse as smoother:  u = u + Cr
   515:           c
   516:           c     This  implementation costs  15A + 4M per result, where
   517:           c     A and M denote the costs of Addition and Multiplication.
   518:           c     Presuming coefficient c(3) is zero (the NPB assumes this,
   519:           c     but it is thus not a general case), 2A + 1M may be eliminated,
   520:           c     resulting in 13A + 3M.
   521:           c     Note that this vectorizes, and is also fine for cache
   522:           c     based machines.
   523:           c---------------------------------------------------------------------
   524:           
   525:                 use mg_data
   526:                 implicit none
   527:           
   528:                 integer n1,n2,n3,k
   529:                 double precision u(n1,n2,n3),r(n1,n2,n3),c(0:3)
   530:                 integer i3, i2, i1
   531:           
   532:                 double precision r1(m), r2(m)
   533:           
   534:                 if (timeron) call timer_start(T_psinv)
   535:           !$omp parallel do default(shared) private(i1,i2,i3,r1,r2) collapse(2)
   536: Y------>        do i3=2,n3-1
   537: |+----->           do i2=2,n2-1
   538: ||V---->              do i1=1,n1
   539: |||                      r1(i1) = r(i1,i2-1,i3) + r(i1,i2+1,i3)
   540: |||            >                + r(i1,i2,i3-1) + r(i1,i2,i3+1)
   541: |||                      r2(i1) = r(i1,i2-1,i3-1) + r(i1,i2+1,i3-1)
   542: |||            >                + r(i1,i2-1,i3+1) + r(i1,i2+1,i3+1)
   543: ||V----               enddo
   544: ||V---->              do i1=2,n1-1
   545: |||     F                u(i1,i2,i3) = u(i1,i2,i3)
   546: |||            >                     + c(0) * r(i1,i2,i3)
   547: |||            >                     + c(1) * ( r(i1-1,i2,i3) + r(i1+1,i2,i3)
   548: |||            >                              + r1(i1) )
   549: |||            >                     + c(2) * ( r2(i1) + r1(i1-1) + r1(i1+1) )
   550: |||       c---------------------------------------------------------------------
   551: |||       c  Assume c(3) = 0    (Enable line below if c(3) not= 0)
   552: |||       c---------------------------------------------------------------------
   553: |||       c    >                     + c(3) * ( r2(i1-1) + r2(i1+1) )
   554: |||       c---------------------------------------------------------------------
   555: ||V----               enddo
   556: |+-----            enddo
   557: Y------         enddo
   558:                 if (timeron) call timer_stop(T_psinv)
   559:           
   560:           c---------------------------------------------------------------------
   561:           c     exchange boundary points
   562:           c---------------------------------------------------------------------
   563:                 call comm3(u,n1,n2,n3,k)
   564:           
   565:                 if( debug_vec(0) .ge. 1 )then
   566:                    call rep_nrm(u,n1,n2,n3,'   psinv',k)
   567:                 endif
   568:           
   569:                 if( debug_vec(3) .ge. k )then
   570:                    call showall(u,n1,n2,n3)
   571:                 endif
   572:           
   573:                 return
   574:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: PSINV
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:536)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:536)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:536)

    LOOP BEGIN: (mg.f:538)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:538)
      *** The number of VLOAD, VSTORE. :  8,  2. (mg.f:538)
    LOOP END

    LOOP BEGIN: (mg.f:544)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:544)
      *** The number of VLOAD, VSTORE. :  8,  1. (mg.f:544)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: PSINV
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 55 [s0-s4 s7-s12 s15-s16 s18 s23-s63]
      Vector registers         : 42 [v22-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 8528 bytes
      Register spill area      :    8 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 8328 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:536)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:536)
    *** Estimated execution cycle                       : 50
    *** The number of SCALAR REGISTER TRANSFER          : 23

    LOOP BEGIN: (mg.f:538)
      *** Estimated execution cycle                     : 459
      *** The number of SCALAR REGISTER TRANSFER        : 7
    LOOP END

    LOOP BEGIN: (mg.f:544)
      *** Estimated execution cycle                     : 490
      *** The number of SCALAR REGISTER TRANSFER        : 8
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RESID
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   606: par(1801): Parallel routine generated.: RESID$1
   607: par(1803): Parallelized by "do".
   607: opt(2495): Loop nest collapsed by "collapse" clause.
   607: vec( 101): Vectorized loop.
   609: vec( 101): Vectorized loop.
   615: vec( 101): Vectorized loop.
   616: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RESID
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   579:                 subroutine resid( u,v,r,n1,n2,n3,a,k )
   580:           
   581:           c---------------------------------------------------------------------
   582:           c---------------------------------------------------------------------
   583:           
   584:           c---------------------------------------------------------------------
   585:           c     resid computes the residual:  r = v - Au
   586:           c
   587:           c     This  implementation costs  15A + 4M per result, where
   588:           c     A and M denote the costs of Addition (or Subtraction) and
   589:           c     Multiplication, respectively.
   590:           c     Presuming coefficient a(1) is zero (the NPB assumes this,
   591:           c     but it is thus not a general case), 3A + 1M may be eliminated,
   592:           c     resulting in 12A + 3M.
   593:           c     Note that this vectorizes, and is also fine for cache
   594:           c     based machines.
   595:           c---------------------------------------------------------------------
   596:           
   597:                 use mg_data
   598:                 implicit none
   599:           
   600:                 integer n1,n2,n3,k
   601:                 double precision u(n1,n2,n3),v(n1,n2,n3),r(n1,n2,n3),a(0:3)
   602:                 integer i3, i2, i1
   603:                 double precision u1(m), u2(m)
   604:           
   605:                 if (timeron) call timer_start(T_resid)
   606:           !$omp parallel do default(shared) private(i1,i2,i3,u1,u2) collapse(2)
   607: Y------>        do i3=2,n3-1
   608: |+----->           do i2=2,n2-1
   609: ||V---->              do i1=1,n1
   610: |||                      u1(i1) = u(i1,i2-1,i3) + u(i1,i2+1,i3)
   611: |||            >                + u(i1,i2,i3-1) + u(i1,i2,i3+1)
   612: |||                      u2(i1) = u(i1,i2-1,i3-1) + u(i1,i2+1,i3-1)
   613: |||            >                + u(i1,i2-1,i3+1) + u(i1,i2+1,i3+1)
   614: ||V----               enddo
   615: ||V---->              do i1=2,n1-1
   616: |||     F                r(i1,i2,i3) = v(i1,i2,i3)
   617: |||            >                     - a(0) * u(i1,i2,i3)
   618: |||       c---------------------------------------------------------------------
   619: |||       c  Assume a(1) = 0      (Enable 2 lines below if a(1) not= 0)
   620: |||       c---------------------------------------------------------------------
   621: |||       c    >                     - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
   622: |||       c    >                              + u1(i1) )
   623: |||       c---------------------------------------------------------------------
   624: |||            >                     - a(2) * ( u2(i1) + u1(i1-1) + u1(i1+1) )
   625: |||            >                     - a(3) * ( u2(i1-1) + u2(i1+1) )
   626: ||V----               enddo
   627: |+-----            enddo
   628: Y------         enddo
   629:                 if (timeron) call timer_stop(T_resid)
   630:           
   631:           c---------------------------------------------------------------------
   632:           c     exchange boundary data
   633:           c---------------------------------------------------------------------
   634:                 call comm3(r,n1,n2,n3,k)
   635:           
   636:                 if( debug_vec(0) .ge. 1 )then
   637:                    call rep_nrm(r,n1,n2,n3,'   resid',k)
   638:                 endif
   639:           
   640:                 if( debug_vec(2) .ge. k )then
   641:                    call showall(r,n1,n2,n3)
   642:                 endif
   643:           
   644:                 return
   645:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RESID
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:607)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:607)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:607)

    LOOP BEGIN: (mg.f:609)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:609)
      *** The number of VLOAD, VSTORE. :  8,  2. (mg.f:609)
    LOOP END

    LOOP BEGIN: (mg.f:615)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:615)
      *** The number of VLOAD, VSTORE. :  7,  1. (mg.f:615)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RESID
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 57 [s0-s12 s15-s16 s18 s23-s63]
      Vector registers         : 40 [v24-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 8560 bytes
      Register spill area      :    8 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 8360 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:607)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:607)
    *** Estimated execution cycle                       : 50
    *** The number of SCALAR REGISTER TRANSFER          : 26

    LOOP BEGIN: (mg.f:609)
      *** Estimated execution cycle                     : 459
      *** The number of SCALAR REGISTER TRANSFER        : 7
    LOOP END

    LOOP BEGIN: (mg.f:615)
      *** Estimated execution cycle                     : 428
      *** The number of SCALAR REGISTER TRANSFER        : 8
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RPRJ3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   693: par(1801): Parallel routine generated.: RPRJ3$1
   695: par(1803): Parallelized by "do".
   695: opt(2495): Loop nest collapsed by "collapse" clause.
   695: vec( 101): Vectorized loop.
   700: vec( 101): Vectorized loop.
   708: vec( 101): Vectorized loop.
   714: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RPRJ3
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   650:                 subroutine rprj3( r,m1k,m2k,m3k,s,m1j,m2j,m3j,k )
   651:           
   652:           c---------------------------------------------------------------------
   653:           c---------------------------------------------------------------------
   654:           
   655:           c---------------------------------------------------------------------
   656:           c     rprj3 projects onto the next coarser grid,
   657:           c     using a trilinear Finite Element projection:  s = r' = P r
   658:           c
   659:           c     This  implementation costs  20A + 4M per result, where
   660:           c     A and M denote the costs of Addition and Multiplication.
   661:           c     Note that this vectorizes, and is also fine for cache
   662:           c     based machines.
   663:           c---------------------------------------------------------------------
   664:           
   665:                 use mg_data
   666:                 implicit none
   667:           
   668:                 integer m1k, m2k, m3k, m1j, m2j, m3j,k
   669:                 double precision r(m1k,m2k,m3k), s(m1j,m2j,m3j)
   670:                 integer j3, j2, j1, i3, i2, i1, d1, d2, d3, j
   671:           
   672:                 double precision x1(m), y1(m), x2,y2
   673:           
   674:                 if (timeron) call timer_start(T_rprj3)
   675:                 if(m1k.eq.3)then
   676:                   d1 = 2
   677:                 else
   678:                   d1 = 1
   679:                 endif
   680:           
   681:                 if(m2k.eq.3)then
   682:                   d2 = 2
   683:                 else
   684:                   d2 = 1
   685:                 endif
   686:           
   687:                 if(m3k.eq.3)then
   688:                   d3 = 2
   689:                 else
   690:                   d3 = 1
   691:                 endif
   692:           
   693:           !$omp parallel do default(shared) collapse(2)
   694:           !$omp& private(j1,j2,j3,i1,i2,i3,x1,y1,x2,y2)
   695: Y------>        do  j3=2,m3j-1
   696: |+----->           do  j2=2,m2j-1
   697: ||                    i3 = 2*j3-d3
   698: ||                    i2 = 2*j2-d2
   699: ||        
   700: ||V---->              do j1=2,m1j
   701: |||                     i1 = 2*j1-d1
   702: |||                     x1(i1-1) = r(i1-1,i2-1,i3  ) + r(i1-1,i2+1,i3  )
   703: |||            >                 + r(i1-1,i2,  i3-1) + r(i1-1,i2,  i3+1)
   704: |||                     y1(i1-1) = r(i1-1,i2-1,i3-1) + r(i1-1,i2-1,i3+1)
   705: |||            >                 + r(i1-1,i2+1,i3-1) + r(i1-1,i2+1,i3+1)
   706: ||V----               enddo
   707: ||        
   708: ||V---->              do  j1=2,m1j-1
   709: |||                     i1 = 2*j1-d1
   710: |||                     y2 = r(i1,  i2-1,i3-1) + r(i1,  i2-1,i3+1)
   711: |||            >           + r(i1,  i2+1,i3-1) + r(i1,  i2+1,i3+1)
   712: |||                     x2 = r(i1,  i2-1,i3  ) + r(i1,  i2+1,i3  )
   713: |||            >           + r(i1,  i2,  i3-1) + r(i1,  i2,  i3+1)
   714: |||     F               s(j1,j2,j3) =
   715: |||            >               0.5D0 * r(i1,i2,i3)
   716: |||            >             + 0.25D0 * ( r(i1-1,i2,i3) + r(i1+1,i2,i3) + x2)
   717: |||            >             + 0.125D0 * ( x1(i1-1) + x1(i1+1) + y2)
   718: |||            >             + 0.0625D0 * ( y1(i1-1) + y1(i1+1) )
   719: ||V----               enddo
   720: ||        
   721: |+-----            enddo
   722: Y------         enddo
   723:                 if (timeron) call timer_stop(T_rprj3)
   724:           
   725:           
   726:                 j = k-1
   727:                 call comm3(s,m1j,m2j,m3j,j)
   728:           
   729:                 if( debug_vec(0) .ge. 1 )then
   730:                    call rep_nrm(s,m1j,m2j,m3j,'   rprj3',k-1)
   731:                 endif
   732:           
   733:                 if( debug_vec(4) .ge. k )then
   734:                    call showall(s,m1j,m2j,m3j)
   735:                 endif
   736:           
   737:                 return
   738:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RPRJ3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:695)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:695)
    *** The number of VLOAD, VSTORE. :  0,  4. (mg.f:695)

    LOOP BEGIN: (mg.f:700)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:700)
      *** The number of VLOAD, VSTORE. :  8,  2. (mg.f:700)
    LOOP END

    LOOP BEGIN: (mg.f:708)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:708)
      *** The number of VLOAD, VSTORE. : 15,  1. (mg.f:708)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: RPRJ3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 61 [v3-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 8592 bytes
      Register spill area      :   32 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 8368 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:695)
    *** Estimated execution cycle                       : 703
  LOOP END

  LOOP BEGIN: (mg.f:695)
    *** Estimated execution cycle                       : 66
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER TRANSFER          : 35

    LOOP BEGIN: (mg.f:700)
      *** Estimated execution cycle                     : 459
      *** The number of SCALAR REGISTER TRANSFER        : 11
    LOOP END

    LOOP BEGIN: (mg.f:708)
      *** Estimated execution cycle                     : 972
      *** The number of SCALAR REGISTER TRANSFER        : 7
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: INTERP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   775: par(1801): Parallel routine generated.: INTERP$1
   776: par(1803): Parallelized by "do".
   776: opt(2495): Loop nest collapsed by "collapse" clause.
   776: vec( 101): Vectorized loop.
   779: vec( 101): Vectorized loop.
   785: opt(1112): Loop fused with previous loop.
   785: vec( 101): Vectorized loop.
   788: vec( 128): Fused multiply-add operation applied.
   791: opt(1112): Loop fused with previous loop.: I1
   792: vec( 128): Fused multiply-add operation applied.
   794: vec( 128): Fused multiply-add operation applied.
   797: opt(1112): Loop fused with previous loop.: I1
   798: vec( 128): Fused multiply-add operation applied.
   800: vec( 128): Fused multiply-add operation applied.
   803: opt(1112): Loop fused with previous loop.: I1
   804: vec( 128): Fused multiply-add operation applied.
   806: vec( 128): Fused multiply-add operation applied.
   838: par(1801): Parallel routine generated.: INTERP$2
   840: par(1803): Parallelized by "do".
   840: opt(2495): Loop nest collapsed by "collapse" clause.
   840: vec( 101): Vectorized loop.
   842: vec( 101): Vectorized loop.
   846: vec( 101): Vectorized loop.
   847: vec( 128): Fused multiply-add operation applied.
   853: par(1803): Parallelized by "do".
   853: opt(2495): Loop nest collapsed by "collapse" clause.
   853: vec( 101): Vectorized loop.
   855: vec( 101): Vectorized loop.
   856: vec( 128): Fused multiply-add operation applied.
   859: vec( 101): Vectorized loop.
   860: vec( 128): Fused multiply-add operation applied.
   868: par(1803): Parallelized by "do".
   868: opt(2495): Loop nest collapsed by "collapse" clause.
   868: vec( 101): Vectorized loop.
   870: vec( 101): Vectorized loop.
   871: vec( 128): Fused multiply-add operation applied.
   874: vec( 101): Vectorized loop.
   875: vec( 128): Fused multiply-add operation applied.
   882: par(1803): Parallelized by "do".
   882: opt(2495): Loop nest collapsed by "collapse" clause.
   882: vec( 101): Vectorized loop.
   884: vec( 101): Vectorized loop.
   885: vec( 128): Fused multiply-add operation applied.
   889: vec( 101): Vectorized loop.
   890: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: INTERP
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   743:                 subroutine interp( z,mm1,mm2,mm3,u,n1,n2,n3,k )
   744:           
   745:           c---------------------------------------------------------------------
   746:           c---------------------------------------------------------------------
   747:           
   748:           c---------------------------------------------------------------------
   749:           c     interp adds the trilinear interpolation of the correction
   750:           c     from the coarser grid to the current approximation:  u = u + Qu'
   751:           c
   752:           c     Observe that this  implementation costs  16A + 4M, where
   753:           c     A and M denote the costs of Addition and Multiplication.
   754:           c     Note that this vectorizes, and is also fine for cache
   755:           c     based machines.  Vector machines may get slightly better
   756:           c     performance however, with 8 separate "do i1" loops, rather than 4.
   757:           c---------------------------------------------------------------------
   758:           
   759:                 use mg_data
   760:                 implicit none
   761:           
   762:                 integer mm1, mm2, mm3, n1, n2, n3,k
   763:                 double precision z(mm1,mm2,mm3),u(n1,n2,n3)
   764:                 integer i3, i2, i1, d1, d2, d3, t1, t2, t3
   765:           
   766:           c note that m = 1037 in globals.h but for this only need to be
   767:           c 535 to handle up to 1024^3
   768:           c      integer m
   769:           c      parameter( m=535 )
   770:                 double precision z1(m),z2(m),z3(m)
   771:           
   772:                 if (timeron) call timer_start(T_interp)
   773:                 if( n1 .ne. 3 .and. n2 .ne. 3 .and. n3 .ne. 3 ) then
   774:           
   775:           !$omp parallel do default(shared) private(i1,i2,i3,z1,z2,z3) collapse(2)
   776: Y------>           do  i3=1,mm3-1
   777: |+----->              do  i2=1,mm2-1
   778: ||        
   779: ||V---->                 do i1=1,mm1
   780: |||                         z1(i1) = z(i1,i2+1,i3) + z(i1,i2,i3)
   781: |||                         z2(i1) = z(i1,i2,i3+1) + z(i1,i2,i3)
   782: |||                         z3(i1) = z(i1,i2+1,i3+1) + z(i1,i2,i3+1) + z1(i1)
   783: ||V----                  enddo
   784: ||        
   785: ||V---->                 do  i1=1,mm1-1
   786: |||                         u(2*i1-1,2*i2-1,2*i3-1)=u(2*i1-1,2*i2-1,2*i3-1)
   787: |||            >                 +z(i1,i2,i3)
   788: |||     F                   u(2*i1,2*i2-1,2*i3-1)=u(2*i1,2*i2-1,2*i3-1)
   789: |||            >                 +0.5d0*(z(i1+1,i2,i3)+z(i1,i2,i3))
   790: |||                      enddo
   791: |||                      do i1=1,mm1-1
   792: |||     F                   u(2*i1-1,2*i2,2*i3-1)=u(2*i1-1,2*i2,2*i3-1)
   793: |||            >                 +0.5d0 * z1(i1)
   794: |||     F                   u(2*i1,2*i2,2*i3-1)=u(2*i1,2*i2,2*i3-1)
   795: |||            >                 +0.25d0*( z1(i1) + z1(i1+1) )
   796: |||                      enddo
   797: |||                      do i1=1,mm1-1
   798: |||     F                   u(2*i1-1,2*i2-1,2*i3)=u(2*i1-1,2*i2-1,2*i3)
   799: |||            >                 +0.5d0 * z2(i1)
   800: |||     F                   u(2*i1,2*i2-1,2*i3)=u(2*i1,2*i2-1,2*i3)
   801: |||            >                 +0.25d0*( z2(i1) + z2(i1+1) )
   802: |||                      enddo
   803: |||                      do i1=1,mm1-1
   804: |||     F                   u(2*i1-1,2*i2,2*i3)=u(2*i1-1,2*i2,2*i3)
   805: |||            >                 +0.25d0* z3(i1)
   806: |||     F                   u(2*i1,2*i2,2*i3)=u(2*i1,2*i2,2*i3)
   807: |||            >                 +0.125d0*( z3(i1) + z3(i1+1) )
   808: ||V----                  enddo
   809: |+-----               enddo
   810: Y------            enddo
   811:           
   812:                 else
   813:           
   814:                    if(n1.eq.3)then
   815:                       d1 = 2
   816:                       t1 = 1
   817:                    else
   818:                       d1 = 1
   819:                       t1 = 0
   820:                    endif
   821:            
   822:                    if(n2.eq.3)then
   823:                       d2 = 2
   824:                       t2 = 1
   825:                    else
   826:                       d2 = 1
   827:                       t2 = 0
   828:                    endif
   829:            
   830:                    if(n3.eq.3)then
   831:                       d3 = 2
   832:                       t3 = 1
   833:                    else
   834:                       d3 = 1
   835:                       t3 = 0
   836:                    endif
   837:            
   838:           !$omp parallel default(shared) private(i1,i2,i3)
   839:           !$omp do collapse(2)
   840: Y------>           do  i3=d3,mm3-1
   841: |+----->              do  i2=d2,mm2-1
   842: ||V---->                 do  i1=d1,mm1-1
   843: |||                         u(2*i1-d1,2*i2-d2,2*i3-d3)=u(2*i1-d1,2*i2-d2,2*i3-d3)
   844: |||            >                 +z(i1,i2,i3)
   845: ||V----                  enddo
   846: ||V---->                 do  i1=1,mm1-1
   847: |||     F                   u(2*i1-t1,2*i2-d2,2*i3-d3)=u(2*i1-t1,2*i2-d2,2*i3-d3)
   848: |||            >                 +0.5D0*(z(i1+1,i2,i3)+z(i1,i2,i3))
   849: ||V----                  enddo
   850: |+-----               enddo
   851: Y------            enddo
   852:           !$omp do collapse(2)
   853: Y------>           do  i3=d3,mm3-1
   854: |+----->              do  i2=1,mm2-1
   855: ||V---->                 do  i1=d1,mm1-1
   856: |||     F                   u(2*i1-d1,2*i2-t2,2*i3-d3)=u(2*i1-d1,2*i2-t2,2*i3-d3)
   857: |||            >                 +0.5D0*(z(i1,i2+1,i3)+z(i1,i2,i3))
   858: ||V----                  enddo
   859: ||V---->                 do  i1=1,mm1-1
   860: |||     F                   u(2*i1-t1,2*i2-t2,2*i3-d3)=u(2*i1-t1,2*i2-t2,2*i3-d3)
   861: |||            >                 +0.25D0*(z(i1+1,i2+1,i3)+z(i1+1,i2,i3)
   862: |||            >                 +z(i1,  i2+1,i3)+z(i1,  i2,i3))
   863: ||V----                  enddo
   864: |+-----               enddo
   865: Y------            enddo
   866:           
   867:           !$omp do collapse(2)
   868: Y------>           do  i3=1,mm3-1
   869: |+----->              do  i2=d2,mm2-1
   870: ||V---->                 do  i1=d1,mm1-1
   871: |||     F                   u(2*i1-d1,2*i2-d2,2*i3-t3)=u(2*i1-d1,2*i2-d2,2*i3-t3)
   872: |||            >                 +0.5D0*(z(i1,i2,i3+1)+z(i1,i2,i3))
   873: ||V----                  enddo
   874: ||V---->                 do  i1=1,mm1-1
   875: |||     F                   u(2*i1-t1,2*i2-d2,2*i3-t3)=u(2*i1-t1,2*i2-d2,2*i3-t3)
   876: |||            >                 +0.25D0*(z(i1+1,i2,i3+1)+z(i1,i2,i3+1)
   877: |||            >                 +z(i1+1,i2,i3  )+z(i1,i2,i3  ))
   878: ||V----                  enddo
   879: |+-----               enddo
   880: Y------            enddo
   881:           !$omp do collapse(2)
   882: Y------>           do  i3=1,mm3-1
   883: |+----->              do  i2=1,mm2-1
   884: ||V---->                 do  i1=d1,mm1-1
   885: |||     F                   u(2*i1-d1,2*i2-t2,2*i3-t3)=u(2*i1-d1,2*i2-t2,2*i3-t3)
   886: |||            >                 +0.25D0*(z(i1,i2+1,i3+1)+z(i1,i2,i3+1)
   887: |||            >                 +z(i1,i2+1,i3  )+z(i1,i2,i3  ))
   888: ||V----                  enddo
   889: ||V---->                 do  i1=1,mm1-1
   890: |||     F                   u(2*i1-t1,2*i2-t2,2*i3-t3)=u(2*i1-t1,2*i2-t2,2*i3-t3)
   891: |||            >                 +0.125D0*(z(i1+1,i2+1,i3+1)+z(i1+1,i2,i3+1)
   892: |||            >                 +z(i1  ,i2+1,i3+1)+z(i1  ,i2,i3+1)
   893: |||            >                 +z(i1+1,i2+1,i3  )+z(i1+1,i2,i3  )
   894: |||            >                 +z(i1  ,i2+1,i3  )+z(i1  ,i2,i3  ))
   895: ||V----                  enddo
   896: |+-----               enddo
   897: Y------            enddo
   898:           !$omp end do nowait
   899:           !$omp end parallel
   900:           
   901:                 endif
   902:                 if (timeron) call timer_stop(T_interp)
   903:           
   904:                 if( debug_vec(0) .ge. 1 )then
   905:                    call rep_nrm(z,mm1,mm2,mm3,'z: inter',k-1)
   906:                    call rep_nrm(u,n1,n2,n3,'u: inter',k)
   907:                 endif
   908:           
   909:                 if( debug_vec(5) .ge. k )then
   910:                    call showall(z,mm1,mm2,mm3)
   911:                    call showall(u,n1,n2,n3)
   912:                 endif
   913:           
   914:                 return
   915:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: INTERP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:840)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:840)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:840)

    LOOP BEGIN: (mg.f:842)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:842)
      *** The number of VLOAD, VSTORE. :  2,  1. (mg.f:842)
    LOOP END

    LOOP BEGIN: (mg.f:846)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:846)
      *** The number of VLOAD, VSTORE. :  3,  1. (mg.f:846)
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:853)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:853)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:853)

    LOOP BEGIN: (mg.f:855)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:855)
      *** The number of VLOAD, VSTORE. :  3,  1. (mg.f:855)
    LOOP END

    LOOP BEGIN: (mg.f:859)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:859)
      *** The number of VLOAD, VSTORE. :  5,  1. (mg.f:859)
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:868)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:868)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:868)

    LOOP BEGIN: (mg.f:870)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:870)
      *** The number of VLOAD, VSTORE. :  3,  1. (mg.f:870)
    LOOP END

    LOOP BEGIN: (mg.f:874)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:874)
      *** The number of VLOAD, VSTORE. :  5,  1. (mg.f:874)
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:882)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:882)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:882)

    LOOP BEGIN: (mg.f:884)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:884)
      *** The number of VLOAD, VSTORE. :  5,  1. (mg.f:884)
    LOOP END

    LOOP BEGIN: (mg.f:889)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:889)
      *** The number of VLOAD, VSTORE. :  9,  1. (mg.f:889)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: INTERP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 54 [s0-s2 s5-s12 s15-s16 s18-s19 s23-s29 s32-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 384 bytes
      Register spill area      :  32 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           : 160 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:840)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:840)
    *** Estimated execution cycle                       : 39
    *** The number of SCALAR REGISTER TRANSFER          : 15

    LOOP BEGIN: (mg.f:842)
      *** Estimated execution cycle                     : 104
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (mg.f:846)
      *** Estimated execution cycle                     : 169
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:853)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:853)
    *** Estimated execution cycle                       : 42
    *** The number of SCALAR REGISTER TRANSFER          : 15

    LOOP BEGIN: (mg.f:855)
      *** Estimated execution cycle                     : 170
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (mg.f:859)
      *** Estimated execution cycle                     : 298
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:868)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:868)
    *** Estimated execution cycle                       : 42
    *** The number of SCALAR REGISTER TRANSFER          : 15

    LOOP BEGIN: (mg.f:870)
      *** Estimated execution cycle                     : 170
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (mg.f:874)
      *** Estimated execution cycle                     : 297
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:882)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:882)
    *** Estimated execution cycle                       : 49
    *** The number of SCALAR REGISTER TRANSFER          : 15

    LOOP BEGIN: (mg.f:884)
      *** Estimated execution cycle                     : 298
      *** The number of SCALAR REGISTER TRANSFER        : 6
    LOOP END

    LOOP BEGIN: (mg.f:889)
      *** Estimated execution cycle                     : 554
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: NORM2U3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   951: par(1801): Parallel routine generated.: NORM2U3$1
   953: par(1803): Parallelized by "do".
   953: opt(2495): Loop nest collapsed by "collapse" clause.
   953: vec( 101): Vectorized loop.
   955: vec( 101): Vectorized loop.
   956: vec( 126): Idiom detected.: SUM
   956: vec( 128): Fused multiply-add operation applied.
   958: vec( 126): Idiom detected.: MAX/MIN
   961: par(1807): Critical-construct generated.
   961: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: NORM2U3
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   920:                 subroutine norm2u3(r,n1,n2,n3,rnm2,rnmu,nx,ny,nz)
   921:           
   922:           c---------------------------------------------------------------------
   923:           c---------------------------------------------------------------------
   924:           
   925:           c---------------------------------------------------------------------
   926:           c     norm2u3 evaluates approximations to the L2 norm and the
   927:           c     uniform (or L-infinity or Chebyshev) norm, under the
   928:           c     assumption that the boundaries are periodic or zero.  Add the
   929:           c     boundaries in with half weight (quarter weight on the edges
   930:           c     and eighth weight at the corners) for inhomogeneous boundaries.
   931:           c---------------------------------------------------------------------
   932:                 use mg_data, only : timeron
   933:           
   934:                 implicit none
   935:           
   936:                 integer n1, n2, n3, nx, ny, nz
   937:                 double precision rnm2, rnmu, r(n1,n2,n3)
   938:                 double precision s, a
   939:                 integer i3, i2, i1
   940:           
   941:                 double precision dn
   942:           
   943:                 integer T_norm2
   944:                 parameter (T_norm2=9)
   945:           
   946:                 if (timeron) call timer_start(T_norm2)
   947:                 dn = 1.0d0*nx*ny*nz
   948:           
   949:                 s=0.0D0
   950:                 rnmu = 0.0D0
   951:           !$omp parallel do default(shared) private(i1,i2,i3,a) collapse(2)
   952:           !$omp& reduction(+:s) reduction(max:rnmu)
   953: Y------>        do  i3=2,n3-1
   954: |+----->           do  i2=2,n2-1
   955: ||V---->              do  i1=2,n1-1
   956: |||     F                s=s+r(i1,i2,i3)**2
   957: |||                      a=abs(r(i1,i2,i3))
   958: |||                      rnmu=dmax1(rnmu,a)
   959: ||V----               enddo
   960: |+-----            enddo
   961: Y------         enddo
   962:           
   963:                 rnm2=sqrt( s / dn )
   964:                 if (timeron) call timer_stop(T_norm2)
   965:           
   966:                 return
   967:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: NORM2U3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:953)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:953)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:953)

    LOOP BEGIN: (mg.f:955)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:955)
      *** The number of VLOAD, VSTORE. :  1,  0. (mg.f:955)
      *** Idiom detected. : SUM (mg.f:956)
      *** Idiom detected. : MAX/MIN (mg.f:958)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: NORM2U3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 40 [s0-s2 s8-s12 s15-s16 s18 s23-s30 s43-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 4368 bytes
      Register spill area      : 2056 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 2120 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:953)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:953)
    *** Estimated execution cycle                       : 127
    *** The number of SCALAR REGISTER TRANSFER          : 7

    LOOP BEGIN: (mg.f:955)
      *** Estimated execution cycle                     : 133
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: REP_NRM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: REP_NRM
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   972:                 subroutine rep_nrm(u,n1,n2,n3,title,kk)
   973:           
   974:           c---------------------------------------------------------------------
   975:           c---------------------------------------------------------------------
   976:           
   977:           c---------------------------------------------------------------------
   978:           c     report on norm
   979:           c---------------------------------------------------------------------
   980:           
   981:                 use mg_data
   982:                 implicit none
   983:           
   984:                 integer n1, n2, n3, kk
   985:                 double precision u(n1,n2,n3)
   986:                 character*8 title
   987:           
   988:                 double precision rnm2, rnmu
   989:           
   990:           
   991:                 call norm2u3(u,n1,n2,n3,rnm2,rnmu,nx(kk),ny(kk),nz(kk))
   992:                 write(*,7)kk,title,rnm2,rnmu
   993:            7    format(' Level',i2,' in ',a8,': norms =',D21.14,D21.14)
   994:           
   995:                 return
   996:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: REP_NRM
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: REP_NRM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 23 [s0-s12 s15-s16 s23-s27 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1296 bytes
      Register spill area      :    0 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1040 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: COMM3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1018: par(1801): Parallel routine generated.: COMM3$1
  1020: par(1803): Parallelized by "do".
  1021: vec( 101): Vectorized loop.
  1028: vec( 101): Vectorized loop.
  1035: par(1803): Parallelized by "do".
  1036: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: COMM3
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1001:                 subroutine comm3(u,n1,n2,n3,kk)
  1002:           
  1003:           c---------------------------------------------------------------------
  1004:           c---------------------------------------------------------------------
  1005:           
  1006:           c---------------------------------------------------------------------
  1007:           c     comm3 organizes the communication on all borders
  1008:           c---------------------------------------------------------------------
  1009:           
  1010:                 use mg_data
  1011:                 implicit none
  1012:           
  1013:                 integer n1, n2, n3, kk
  1014:                 double precision u(n1,n2,n3)
  1015:                 integer i1, i2, i3
  1016:           
  1017:                 if (timeron) call timer_start(T_comm3)
  1018:           !$omp parallel default(shared) private(i1,i2,i3)
  1019:           !$omp do
  1020: P------>        do  i3=2,n3-1
  1021: |V----->           do  i2=2,n2-1
  1022: ||                    u( 1,i2,i3) = u(n1-1,i2,i3)
  1023: ||                    u(n1,i2,i3) = u(   2,i2,i3)
  1024: |V-----            enddo
  1025: |         c      enddo
  1026: |         
  1027: |         c      do  i3=2,n3-1
  1028: |V----->           do  i1=1,n1
  1029: ||                    u(i1, 1,i3) = u(i1,n2-1,i3)
  1030: ||                    u(i1,n2,i3) = u(i1,   2,i3)
  1031: |V-----            enddo
  1032: P------         enddo
  1033:           
  1034:           !$omp do
  1035: P------>        do  i2=1,n2
  1036: |V----->           do  i1=1,n1
  1037: ||                    u(i1,i2, 1) = u(i1,i2,n3-1)
  1038: ||                    u(i1,i2,n3) = u(i1,i2,   2)
  1039: |V-----            enddo
  1040: P------         enddo
  1041:           !$omp end do nowait
  1042:           !$omp end parallel
  1043:                 if (timeron) call timer_stop(T_comm3)
  1044:           
  1045:                 return
  1046:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: COMM3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:1020)
    <Unvectorized loop.>

    LOOP BEGIN: (mg.f:1021)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:1021)
      *** The number of VLOAD, VSTORE. :  2,  2. (mg.f:1021)
    LOOP END

    LOOP BEGIN: (mg.f:1028)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:1028)
      *** The number of VLOAD, VSTORE. :  2,  2. (mg.f:1028)
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:1035)
    <Unvectorized loop.>

    LOOP BEGIN: (mg.f:1036)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:1036)
      *** The number of VLOAD, VSTORE. :  2,  2. (mg.f:1036)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: COMM3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 42 [s0-s2 s8-s12 s15-s16 s18 s23-s25 s36-s63]
      Vector registers         :  6 [v58-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 240 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  48 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1020)
    *** Estimated execution cycle                       : 27
    *** The number of SCALAR REGISTER TRANSFER          : 19

    LOOP BEGIN: (mg.f:1021)
      *** Estimated execution cycle                     : 73
      *** The number of SCALAR REGISTER TRANSFER        : 8
    LOOP END

    LOOP BEGIN: (mg.f:1028)
      *** Estimated execution cycle                     : 73
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (mg.f:1035)
    *** Estimated execution cycle                       : 12
    *** The number of SCALAR REGISTER TRANSFER          : 3

    LOOP BEGIN: (mg.f:1036)
      *** Estimated execution cycle                     : 73
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZRAN3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1101: vec( 110): Vectorization obstructive procedure reference.: RANDLC
  1101: vec( 103): Unvectorized loop.
  1103: opt(1025): Reference to this procedure inhibits optimization.: RANDLC
  1109: par(1801): Parallel routine generated.: ZRAN3$1
  1111: par(1803): Parallelized by "do".
  1113: vec( 110): Vectorization obstructive procedure reference.: VRANLC
  1113: vec( 110): Vectorization obstructive procedure reference.: RANDLC
  1113: vec( 103): Unvectorized loop.
  1115: opt(1025): Reference to this procedure inhibits optimization.: VRANLC
  1129: par(1801): Parallel routine generated.: ZRAN3$2
  1131: vec( 101): Vectorized loop.
  1143: par(1803): Parallelized by "do".
  1143: opt(2495): Loop nest collapsed by "collapse" clause.
  1145: vec( 110): Vectorization obstructive procedure reference.: BUBBLE
  1145: vec( 103): Unvectorized loop.
  1151: opt(1025): Reference to this procedure inhibits optimization.: BUBBLE
  1180: par(1805): Master-construct generated.
  1186: vec( 109): Vectorization obstructive statement.
  1186: par(1803): Parallelized by "do".
  1186: vec( 103): Unvectorized loop.
  1238: par(1801): Parallel routine generated.: ZRAN3$3
  1239: par(1803): Parallelized by "do".
  1239: opt(2495): Loop nest collapsed by "collapse" clause.
  1239: vec( 101): Vectorized loop.
  1241: vec( 101): Vectorized loop.
  1248: vec( 101): Vectorized loop.
  1251: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZRAN3
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1052:                 subroutine zran3(z,n1,n2,n3,nx1,ny1,k)
  1053:           
  1054:           c---------------------------------------------------------------------
  1055:           c---------------------------------------------------------------------
  1056:           
  1057:           c---------------------------------------------------------------------
  1058:           c     zran3  loads +1 at ten randomly chosen points,
  1059:           c     loads -1 at a different ten random points,
  1060:           c     and zero elsewhere.
  1061:           c---------------------------------------------------------------------
  1062:           
  1063:                 use mg_data
  1064:                 implicit none
  1065:           
  1066:                 integer n1, n2, n3, k, nx1, ny1, i0, mm0, mm1
  1067:                 double precision z(n1,n2,n3)
  1068:           
  1069:                 integer mm, i1, i2, i3, d1, e1, e2, e3
  1070:                 double precision x, a
  1071:                 double precision xx, x0, x1, a1, a2, ai, power
  1072:                 parameter( mm = 10,  a = 5.D0 ** 13, x = 314159265.D0)
  1073:                 double precision ten( mm, 0:1 ), best0, best1
  1074:                 integer i, j1( mm, 0:1 ), j2( mm, 0:1 ), j3( mm, 0:1 )
  1075:                 integer jg( 0:3, mm, 0:1 )
  1076:           
  1077:                 external randlc
  1078:                 double precision randlc, rdummy
  1079:           !$    integer  omp_get_thread_num, omp_get_num_threads
  1080:           !$    external omp_get_thread_num, omp_get_num_threads
  1081:                 integer myid, num_threads
  1082:           
  1083:                 a1 = power( a, nx1 )
  1084:                 a2 = power( a, nx1*ny1 )
  1085:           
  1086:                 call zero3(z,n1,n2,n3)
  1087:           
  1088:                 i = is1-2+nx1*(is2-2+ny1*(is3-2))
  1089:           
  1090:                 ai = power( a, i )
  1091:                 d1 = ie1 - is1 + 1
  1092:                 e1 = ie1 - is1 + 2
  1093:                 e2 = ie2 - is2 + 2
  1094:                 e3 = ie3 - is3 + 2
  1095:                 x0 = x
  1096:                 rdummy = randlc( x0, ai )
  1097:           
  1098:           c---------------------------------------------------------------------
  1099:           c     save the starting seeds for the following loop
  1100:           c---------------------------------------------------------------------
  1101: +------>        do  i3 = 2, e3
  1102: |                  starts(i3) = x0
  1103: |                  rdummy = randlc( x0, a2 )
  1104: +------         end do
  1105:           
  1106:           c---------------------------------------------------------------------
  1107:           c     fill array
  1108:           c---------------------------------------------------------------------
  1109:           !$omp parallel do default(shared) private(i2,i3,x1,xx,rdummy)
  1110:           !$omp&  shared(e2,e3,d1,a1)
  1111: P------>        do  i3 = 2, e3
  1112: |                  x1 = starts(i3)
  1113: |+----->           do  i2 = 2, e2
  1114: ||                    xx = x1
  1115: ||                    call vranlc( d1, xx, a, z( 2, i2, i3 ))
  1116: ||                    rdummy = randlc( x1, a1 )
  1117: |+-----            enddo
  1118: P------         enddo
  1119:           !$omp end parallel do
  1120:           
  1121:           c---------------------------------------------------------------------
  1122:           c       call comm3(z,n1,n2,n3)
  1123:           c       call showall(z,n1,n2,n3)
  1124:           c---------------------------------------------------------------------
  1125:           
  1126:           c---------------------------------------------------------------------
  1127:           c     each thread looks for twenty candidates
  1128:           c---------------------------------------------------------------------
  1129:           !$omp parallel default(shared) private(i,i0,i1,i2,i3,j1,j2,j3,ten,
  1130:           !$omp&  myid,num_threads) shared(best0,best1,n1,n2,n3)
  1131: V------>        do  i=1,mm
  1132: |                  ten( i, 1 ) = 0.0D0
  1133: |                  j1( i, 1 ) = 0
  1134: |                  j2( i, 1 ) = 0
  1135: |                  j3( i, 1 ) = 0
  1136: |                  ten( i, 0 ) = 1.0D0
  1137: |                  j1( i, 0 ) = 0
  1138: |                  j2( i, 0 ) = 0
  1139: |                  j3( i, 0 ) = 0
  1140: V------         enddo
  1141:           
  1142:           !$omp do collapse(2)
  1143: P------>        do  i3=2,n3-1
  1144: |+----->           do  i2=2,n2-1
  1145: ||+---->              do  i1=2,n1-1
  1146: |||                      if( z(i1,i2,i3) .gt. ten( 1, 1 ) )then
  1147: |||                         ten(1,1) = z(i1,i2,i3)
  1148: |||                         j1(1,1) = i1
  1149: |||                         j2(1,1) = i2
  1150: |||                         j3(1,1) = i3
  1151: |||                         call bubble( ten, j1, j2, j3, mm, 1 )
  1152: |||                      endif
  1153: |||                      if( z(i1,i2,i3) .lt. ten( 1, 0 ) )then
  1154: |||                         ten(1,0) = z(i1,i2,i3)
  1155: |||                         j1(1,0) = i1
  1156: |||                         j2(1,0) = i2
  1157: |||                         j3(1,0) = i3
  1158: |||                         call bubble( ten, j1, j2, j3, mm, 0 )
  1159: |||                      endif
  1160: ||+----               enddo
  1161: |+-----            enddo
  1162: P------         enddo
  1163:           !$omp end do
  1164:           
  1165:           
  1166:           c---------------------------------------------------------------------
  1167:           c     Now which of these are globally best?
  1168:           c---------------------------------------------------------------------
  1169:                 i1 = mm
  1170:                 i0 = mm
  1171:                 myid = 0
  1172:           !$    myid = omp_get_thread_num()
  1173:           !$    num_threads = omp_get_num_threads()
  1174: +------>        do  i=mm,1,-1
  1175: |         
  1176: |         c ... ORDERED access is required here for sequential consistency
  1177: |         c ... in case that two values are identical.
  1178: |         c ... Since an "ORDERED" section is only defined in OpenMP 2,
  1179: |         c ... we use a dummy loop to emulate ordered access in OpenMP 1.x.
  1180: |         !$omp master
  1181: |                  best1 = 0.0D0
  1182: |                  best0 = 1.0D0
  1183: |         !$omp end master
  1184: |         
  1185: |         !$omp do ordered schedule(static)
  1186: |P----->  !$       do i2=1,num_threads
  1187: ||        !$omp ordered
  1188: ||                 if (ten(i1,1) .gt. best1) then
  1189: ||                    best1 = ten(i1,1)
  1190: ||                    jg( 0, i, 1 ) = myid
  1191: ||                 endif
  1192: ||                 if (ten(i0,0) .lt. best0) then
  1193: ||                    best0 = ten(i0,0)
  1194: ||                    jg( 0, i, 0 ) = myid
  1195: ||                 endif
  1196: ||        !$omp end ordered
  1197: |P-----   !$       end do
  1198: |         
  1199: |                  if (myid .eq. jg( 0, i, 1 )) then
  1200: |                     jg( 1, i, 1 ) = j1( i1, 1 )
  1201: |                     jg( 2, i, 1 ) = j2( i1, 1 )
  1202: |                     jg( 3, i, 1 ) = j3( i1, 1 )
  1203: |                     i1 = i1-1
  1204: |                  endif
  1205: |         
  1206: |                  if (myid .eq. jg( 0, i, 0 )) then
  1207: |                     jg( 1, i, 0 ) = j1( i0, 0 )
  1208: |                     jg( 2, i, 0 ) = j2( i0, 0 )
  1209: |                     jg( 3, i, 0 ) = j3( i0, 0 )
  1210: |                     i0 = i0-1
  1211: |                  endif
  1212: |         
  1213: +------         enddo
  1214:           !$omp end parallel
  1215:           
  1216:           c      mm1 = i1+1
  1217:           c      mm0 = i0+1
  1218:                 mm1 = 1
  1219:                 mm0 = 1
  1220:           
  1221:           c     write(*,*)' '
  1222:           c     write(*,*)' negative charges at'
  1223:           c     write(*,9)(jg(1,i,0),jg(2,i,0),jg(3,i,0),i=1,mm)
  1224:           c     write(*,*)' positive charges at'
  1225:           c     write(*,9)(jg(1,i,1),jg(2,i,1),jg(3,i,1),i=1,mm)
  1226:           c     write(*,*)' small random numbers were'
  1227:           c     write(*,8)(ten( i,0),i=mm,1,-1)
  1228:           c     write(*,*)' and they were found on processor number'
  1229:           c     write(*,7)(jg(0,i,0),i=mm,1,-1)
  1230:           c     write(*,*)' large random numbers were'
  1231:           c     write(*,8)(ten( i,1),i=mm,1,-1)
  1232:           c     write(*,*)' and they were found on processor number'
  1233:           c     write(*,7)(jg(0,i,1),i=mm,1,-1)
  1234:           c 9    format(5(' (',i3,2(',',i3),')'))
  1235:           c 8    format(5D15.8)
  1236:           c 7    format(10i4)
  1237:           
  1238:           !$omp parallel do default(shared) private(i1,i2,i3) collapse(2)
  1239: Y------>        do  i3=1,n3
  1240: |+----->           do  i2=1,n2
  1241: ||V---->              do  i1=1,n1
  1242: |||                      z(i1,i2,i3) = 0.0D0
  1243: ||V----               enddo
  1244: |+-----            enddo
  1245: Y------         enddo
  1246:           !$omp end parallel do
  1247:           
  1248: V------>        do  i=mm,mm0,-1
  1249: |       C          z( jg(1,i,0), jg(2,i,0), jg(3,i,0) ) = -1.0D0
  1250: V------         enddo
  1251: V------>        do  i=mm,mm1,-1
  1252: |       C          z( jg(1,i,1), jg(2,i,1), jg(3,i,1) ) = +1.0D0
  1253: V------         enddo
  1254:           
  1255:                 call comm3(z,n1,n2,n3,k)
  1256:           
  1257:           c---------------------------------------------------------------------
  1258:           c          call showall(z,n1,n2,n3)
  1259:           c---------------------------------------------------------------------
  1260:           
  1261:                 return
  1262:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZRAN3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:1239)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:1239)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:1239)

    LOOP BEGIN: (mg.f:1241)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:1241)
      *** The number of VLOAD, VSTORE. :  0,  1. (mg.f:1241)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZRAN3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 37 [s0-s2 s8-s12 s15-s16 s18 s23-s28 s44-s63]
      Vector registers         : 14 [v50-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 256 bytes
      Register spill area      :   8 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  56 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1239)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:1239)
    *** Estimated execution cycle                       : 22
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (mg.f:1241)
      *** Estimated execution cycle                     : 39
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SHOWALL
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1287: opt(1118): This I/O statement inhibits optimization of loop.
  1287: vec( 180): I/O statement obstructs vectorization.
  1287: vec( 103): Unvectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SHOWALL
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1268:                 subroutine showall(z,n1,n2,n3)
  1269:           
  1270:           c---------------------------------------------------------------------
  1271:           c---------------------------------------------------------------------
  1272:           
  1273:                 implicit none
  1274:           
  1275:           
  1276:                 integer n1,n2,n3,i1,i2,i3
  1277:                 double precision z(n1,n2,n3)
  1278:                 integer m1, m2, m3
  1279:           
  1280:                 m1 = min(n1,18)
  1281:                 m2 = min(n2,14)
  1282:                 m3 = min(n3,18)
  1283:           
  1284:                 write(*,*)'  '
  1285: +------>        do  i3=1,m3
  1286: |+----->           do  i1=1,m1
  1287: ||+---->              write(*,6)(z(i1,i2,i3),i2=1,m2)
  1288: |+-----            enddo
  1289: |                  write(*,*)' - - - - - - - '
  1290: +------         enddo
  1291:                 write(*,*)'  '
  1292:            6    format(15f6.3)
  1293:           
  1294:                 return
  1295:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SHOWALL
VECTORIZATION LIST

  LOOP BEGIN: (mg.f:1285)
    <Unvectorized loop.>

    LOOP BEGIN: (mg.f:1286)
      <Unvectorized loop.>

      LOOP BEGIN: (mg.f:1287)
        <Unvectorized loop.>
        *** I/O statement obstructs vectorization. (mg.f:1287)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: SHOWALL
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 39 [s0-s12 s15-s16 s18-s38 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1424 bytes
      Register spill area      :   72 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1096 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1285)
    *** Estimated execution cycle                       : 49
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 12
            Across calls                                :  4
            Over basic blocks                           :  8
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 6
            Across calls                                : 4
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER TRANSFER          : 17

    LOOP BEGIN: (mg.f:1286)
      *** Estimated execution cycle                     : 20
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 8
              Over basic blocks                         : 8
      *** The number of SCALAR REGISTER TRANSFER        : 17

      LOOP BEGIN: (mg.f:1287)
        *** Estimated execution cycle                   : 8
        *** The number of SCALAR REGISTER TRANSFER      : 6
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: POWER
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1301: vec( 103): Unvectorized loop.
  1301: vec( 110): Vectorization obstructive procedure reference.: RANDLC


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: POWER
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1301:                 double precision function power( a, n )
  1302:           
  1303:           c---------------------------------------------------------------------
  1304:           c---------------------------------------------------------------------
  1305:           
  1306:           c---------------------------------------------------------------------
  1307:           c     power  raises an integer, disguised as a double
  1308:           c     precision real, to an integer power
  1309:           c---------------------------------------------------------------------
  1310:                 implicit none
  1311:           
  1312:                 double precision a, aj
  1313:                 integer n, nj
  1314:                 external randlc
  1315:                 double precision randlc, rdummy
  1316:           
  1317:                 power = 1.0D0
  1318:                 nj = n
  1319:                 aj = a
  1320:            100  continue
  1321:           
  1322:                 if( nj .eq. 0 ) goto 200
  1323:                 if( mod(nj,2) .eq. 1 ) rdummy =  randlc( power, aj )
  1324:                 rdummy = randlc( aj, aj )
  1325:                 nj = nj/2
  1326:                 go to 100
  1327:           
  1328:            200  continue
  1329:                 return
  1330:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: POWER
VECTORIZATION LIST

  LOOP BEGIN: (mg.f:1301)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : RANDLC (mg.f:1301)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: POWER
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 15 [s0-s1 s8-s12 s15-s16 s23-s25 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1301)
    *** Estimated execution cycle                       : 17
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: BUBBLE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1353: vec( 103): Unvectorized loop.
  1357: vec( 120): Unvectorizable dependency.: TEN
  1361: vec( 120): Unvectorizable dependency.: J1
  1365: vec( 120): Unvectorizable dependency.: J2
  1369: vec( 120): Unvectorizable dependency.: J3
  1373: opt(1084): Branch out of the loop inhibits optimization.
  1379: vec( 103): Unvectorized loop.
  1383: vec( 120): Unvectorizable dependency.: TEN
  1387: vec( 120): Unvectorizable dependency.: J1
  1391: vec( 120): Unvectorizable dependency.: J2
  1395: vec( 120): Unvectorizable dependency.: J3
  1399: opt(1084): Branch out of the loop inhibits optimization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: BUBBLE
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1335:                 subroutine bubble( ten, j1, j2, j3, m, ind )
  1336:           
  1337:           c---------------------------------------------------------------------
  1338:           c---------------------------------------------------------------------
  1339:           
  1340:           c---------------------------------------------------------------------
  1341:           c     bubble        does a bubble sort in direction dir
  1342:           c---------------------------------------------------------------------
  1343:                 implicit none
  1344:           
  1345:           
  1346:                 integer m, ind, j1( m, 0:1 ), j2( m, 0:1 ), j3( m, 0:1 )
  1347:                 double precision ten( m, 0:1 )
  1348:                 double precision temp
  1349:                 integer i, j_temp
  1350:           
  1351:                 if( ind .eq. 1 )then
  1352:           
  1353: +------>           do  i=1,m-1
  1354: |                     if( ten(i,ind) .gt. ten(i+1,ind) )then
  1355: |         
  1356: |                        temp = ten( i+1, ind )
  1357: |                        ten( i+1, ind ) = ten( i, ind )
  1358: |                        ten( i, ind ) = temp
  1359: |         
  1360: |                        j_temp           = j1( i+1, ind )
  1361: |                        j1( i+1, ind ) = j1( i,   ind )
  1362: |                        j1( i,   ind ) = j_temp
  1363: |         
  1364: |                        j_temp           = j2( i+1, ind )
  1365: |                        j2( i+1, ind ) = j2( i,   ind )
  1366: |                        j2( i,   ind ) = j_temp
  1367: |         
  1368: |                        j_temp           = j3( i+1, ind )
  1369: |                        j3( i+1, ind ) = j3( i,   ind )
  1370: |                        j3( i,   ind ) = j_temp
  1371: |         
  1372: |                     else
  1373: |                        return
  1374: |                     endif
  1375: +------            enddo
  1376:           
  1377:                 else
  1378:           
  1379: +------>           do  i=1,m-1
  1380: |                     if( ten(i,ind) .lt. ten(i+1,ind) )then
  1381: |         
  1382: |                        temp = ten( i+1, ind )
  1383: |                        ten( i+1, ind ) = ten( i, ind )
  1384: |                        ten( i, ind ) = temp
  1385: |         
  1386: |                        j_temp           = j1( i+1, ind )
  1387: |                        j1( i+1, ind ) = j1( i,   ind )
  1388: |                        j1( i,   ind ) = j_temp
  1389: |         
  1390: |                        j_temp           = j2( i+1, ind )
  1391: |                        j2( i+1, ind ) = j2( i,   ind )
  1392: |                        j2( i,   ind ) = j_temp
  1393: |         
  1394: |                        j_temp           = j3( i+1, ind )
  1395: |                        j3( i+1, ind ) = j3( i,   ind )
  1396: |                        j3( i,   ind ) = j_temp
  1397: |         
  1398: |                     else
  1399: |                        return
  1400: |                     endif
  1401: +------            enddo
  1402:           
  1403:                 endif
  1404:           
  1405:                 return
  1406:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: BUBBLE
VECTORIZATION LIST

  LOOP BEGIN: (mg.f:1379)
    <Unvectorized loop.>
    *** Unvectorizable dependency. : TEN (mg.f:1383)
    *** Unvectorizable dependency. : J1 (mg.f:1387)
    *** Unvectorizable dependency. : J2 (mg.f:1391)
    *** Unvectorizable dependency. : J3 (mg.f:1395)
  LOOP END

  LOOP BEGIN: (mg.f:1353)
    <Unvectorized loop.>
    *** Unvectorizable dependency. : TEN (mg.f:1357)
    *** Unvectorizable dependency. : J1 (mg.f:1361)
    *** Unvectorizable dependency. : J2 (mg.f:1365)
    *** Unvectorizable dependency. : J3 (mg.f:1369)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: BUBBLE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 20 [s0-s5 s8-s11 s18-s19 s56-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1379)
    *** Estimated execution cycle                       : 30
  LOOP END

  LOOP BEGIN: (mg.f:1353)
    *** Estimated execution cycle                       : 30
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZERO3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1423: par(1801): Parallel routine generated.: ZERO3$1
  1424: par(1803): Parallelized by "do".
  1424: opt(2495): Loop nest collapsed by "collapse" clause.
  1424: vec( 101): Vectorized loop.
  1426: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZERO3
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

  1411:                 subroutine zero3(z,n1,n2,n3)
  1412:           
  1413:           c---------------------------------------------------------------------
  1414:           c---------------------------------------------------------------------
  1415:           
  1416:                 implicit none
  1417:           
  1418:           
  1419:                 integer n1, n2, n3
  1420:                 double precision z(n1,n2,n3)
  1421:                 integer i1, i2, i3
  1422:           
  1423:           !$omp parallel do default(shared) private(i1,i2,i3) collapse(2)
  1424: Y------>        do  i3=1,n3
  1425: |+----->           do  i2=1,n2
  1426: ||V---->              do  i1=1,n1
  1427: |||                      z(i1,i2,i3)=0.0D0
  1428: ||V----               enddo
  1429: |+-----            enddo
  1430: Y------         enddo
  1431:           
  1432:                 return
  1433:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZERO3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (mg.f:1424)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (mg.f:1424)
    *** The number of VLOAD, VSTORE. :  0,  2. (mg.f:1424)

    LOOP BEGIN: (mg.f:1426)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (mg.f:1426)
      *** The number of VLOAD, VSTORE. :  0,  1. (mg.f:1426)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: mg.f

PROCEDURE NAME: ZERO3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 37 [s0-s2 s8-s12 s15-s16 s18 s23-s28 s44-s63]
      Vector registers         : 14 [v50-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 256 bytes
      Register spill area      :   8 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  56 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (mg.f:1424)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (mg.f:1424)
    *** Estimated execution cycle                       : 22
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (mg.f:1426)
      *** Estimated execution cycle                     : 39
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


