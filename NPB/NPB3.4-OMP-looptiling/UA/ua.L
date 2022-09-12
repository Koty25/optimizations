NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:27 2020
FILE NAME: ua.f

PROCEDURE NAME: UA
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    65: inl(1212): Source for routine not found.: CHECK_TIMER_FLAG
   122: vec( 103): Unvectorized loop.
   122: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   123: inl(1212): Source for routine not found.: TIMER_CLEAR
   125: inl(1212): Source for routine not found.: TIMER_START
   127: inl(1212): Source for routine not found.: ALLOC_SPACE
   130: inl(1212): Source for routine not found.: CREATE_INITIAL_GRID
   132: inl(1212): Source for routine not found.: R_INIT_OMP
   133: inl(1212): Source for routine not found.: NR_INIT_OMP
   135: inl(1212): Source for routine not found.: INIT_LOCKS
   138: inl(1212): Source for routine not found.: COEF
   139: inl(1212): Source for routine not found.: GEOM1
   142: inl(1212): Source for routine not found.: SETDEF
   145: inl(1212): Source for routine not found.: SETPCMO_PRE
   149: inl(1212): Source for routine not found.: MORTAR
   150: inl(1212): Source for routine not found.: PREPWORK
   151: inl(1212): Source for routine not found.: ADAPTATION
   152: inl(1212): Source for routine not found.: TIMER_STOP
   162: inl(1212): Source for routine not found.: R_INIT
   166: vec( 103): Unvectorized loop.
   166: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   173: inl(1212): Source for routine not found.: CONVECT
   177: inl(1212): Source for routine not found.: TRANSF
   182: par(1801): Parallel routine generated.: MAIN_$1
   184: vec( 110): Vectorization obstructive procedure reference.: LAPLACIAN
   184: par(1803): Parallelized by "do".
   184: vec( 103): Unvectorized loop.
   185: inl(1212): Source for routine not found.: LAPLACIAN
   190: par(1803): Parallelized by "do".
   190: vec( 103): Unvectorized loop.
   190: vec( 113): Overhead of loop division is too large.
   191: vec( 101): Vectorized loop.
   192: opt(1418): Constant-length loop is expanded.
   193: opt(1418): Constant-length loop is expanded.
   202: inl(1212): Source for routine not found.: TRANSFB
   208: par(1801): Parallel routine generated.: MAIN_$2
   209: par(1803): Parallelized by "do".
   210: vec( 103): Unvectorized loop.
   210: vec( 110): Vectorization obstructive procedure reference.: FACEV
   212: inl(1212): Source for routine not found.: FACEV
   218: inl(1212): Source for routine not found.: COL2
   221: inl(1212): Source for routine not found.: DIFFUSION
   225: inl(1212): Source for routine not found.: ADD2
   242: inl(1212): Source for routine not found.: TIMER_READ
   244: inl(1212): Source for routine not found.: VERIFY
   250: inl(1212): Source for routine not found.: PRINT_RESULTS
   260: vec( 103): Unvectorized loop.
   260: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   267: vec( 103): Unvectorized loop.
   267: vec( 180): I/O statement obstructs vectorization.
   268: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:27 2020
FILE NAME: ua.f

PROCEDURE NAME: UA
FORMAT LIST

 LINE   LOOP      STATEMENT

    44:                 program ua
    45:           
    46:                 use ua_data
    47:                 implicit none
    48:           
    49:                 integer          step, ie,iside,i,j, fstatus,k
    50:                 external         timer_read
    51:                 double precision timer_read, mflops, tmax, nelt_tot
    52:                 character        class
    53:                 logical          ifmortar, verified
    54:           !$    integer          omp_get_max_threads
    55:           !$    external         omp_get_max_threads
    56:           
    57:                 double precision t2, trecs(t_last)
    58:                 character t_names(t_last)*10
    59:           
    60:           c---------------------------------------------------------------------
    61:           c     Read input file (if it exists), else take
    62:           c     defaults from parameters
    63:           c---------------------------------------------------------------------
    64:            
    65:                 call check_timer_flag( timeron )
    66:                 if (timeron) then
    67:                    t_names(t_total) = 'total'
    68:                    t_names(t_init) = 'init'
    69:                    t_names(t_convect) = 'convect'
    70:                    t_names(t_transfb_c) = 'transfb_c'
    71:                    t_names(t_diffusion) = 'diffusion'
    72:                    t_names(t_transf) = 'transf'
    73:                    t_names(t_transfb) = 'transfb'
    74: +======>           t_names(t_adaptation) = 'adaptation'
    75:                    t_names(t_transf2) = 'transf+b'
    76:                    t_names(t_add2) = 'add2'
    77:                 endif
    78:           
    79:                 write (*,1000)
    80:                 open (unit=2,file='inputua.data',status='old', iostat=fstatus)
    81:           
    82:                 if (fstatus .eq. 0) then
    83:                   write(*,233)
    84:            233    format(' Reading from input file inputua.data')
    85:                   read (2,*) fre
    86:                   read (2,*) niter
    87:                   read (2,*) nmxh
    88:                   read (2,*) alpha
    89:                   class = 'U'
    90:                   close(2)
    91:                 else
    92:                   write(*,234)
    93:                   fre        = fre_default
    94:                   niter      = niter_default
    95:                   nmxh       = nmxh_default
    96:                   alpha      = alpha_default
    97:                   class      = class_default
    98:                 endif
    99:            234  format(' No input file inputua.data. Using compiled defaults')
   100:           
   101:                 dlmin = 0.5d0**refine_max
   102:                 dtime = 0.04d0*dlmin
   103:           
   104:                 write (*,1001) refine_max
   105:                 write (*,1002) fre
   106:                 write (*,1003) niter, dtime
   107:                 write (*,1004) nmxh
   108:                 write (*,1005) alpha
   109:           !$    write (*,1006) omp_get_max_threads()
   110:                 write (*,*)
   111:           
   112:           
   113:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
   114:                >          ' - UA Benchmark', /)
   115:            1001 format(' Levels of refinement:        ', i8)
   116:            1002 format(' Adaptation frequency:        ', i8)
   117:            1003 format(' Time steps:                  ', i8, '    dt: ', g15.6)
   118:            1004 format(' CG iterations:               ', i8)
   119:            1005 format(' Heat source radius:          ', f8.4)
   120:            1006 format(' Number of available threads: ', i8)
   121:           
   122: +------>        do i = 1, t_last
   123: |                  call timer_clear(i)
   124: +------         end do
   125:                 if (timeron) call timer_start(t_init)
   126:           
   127:                 call alloc_space
   128:           
   129:           c.....set up initial mesh (single element) and solution (all zero)
   130:                 call create_initial_grid
   131:           
   132:                 call r_init_omp(ta1,ntot,0.d0)
   133:                 call nr_init_omp(sje,4*6*nelt,0)
   134:           
   135:                 call init_locks
   136:           
   137:           c.....compute tables of coefficients and weights
   138:                 call coef
   139:                 call geom1
   140:           
   141:           c.....compute the discrete laplacian operators
   142:                 call setdef
   143:           
   144:           c.....prepare for the preconditioner
   145:                 call setpcmo_pre
   146:           
   147:           c.....refine initial mesh and do some preliminary work
   148:                 time = 0.d0
   149:                 call mortar
   150:                 call prepwork
   151:                 call adaptation(ifmortar,0)
   152:                 if (timeron) call timer_stop(t_init)
   153:           
   154:                 call timer_clear(1)
   155:           
   156:                 time = 0.d0
   157: +------>        do step= 0, niter
   158: |         
   159: |                 if (step .eq. 1) then
   160: |         c.........reset the solution and start the timer, keep track of total no elms
   161: |         
   162: |                   call r_init(ta1,ntot,0.d0)
   163: |         
   164: |                   time = 0.d0
   165: |                   nelt_tot = 0.d0
   166: |+----->            do i = 1, t_last
   167: ||                     if (i.ne.t_init) call timer_clear(i)
   168: |+-----             end do
   169: |                   call timer_start(1)
   170: |                 endif
   171: |         
   172: |         c.......advance the convection step
   173: |                 call convect(ifmortar)
   174: |         
   175: |                 if (timeron) call timer_start(t_transf2)
   176: |         c.......prepare the intital guess for cg
   177: |                 call transf(tmort,ta1)
   178: |         
   179: |         c.......compute residual for diffusion term based on intital guess
   180: |         
   181: |         c.......compute the left hand side of equation, lapacian t
   182: |         c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie,k,j,i)
   183: |         c$OMP DO
   184: |P----->          do ie = 1,nelt
   185: ||                  call laplacian(ta2(1,1,1,ie),ta1(1,1,1,ie),size_e(ie))
   186: |P-----           end do
   187: |         c$OMP END DO
   188: |         c.......compute the residual
   189: |         c$OMP DO
   190: |P----->          do ie = 1, nelt
   191: ||V---->            do k=1,lx1
   192: |||*--->              do j=1,lx1
   193: ||||*-->                do i=1,lx1
   194: |||||                     trhs(i,j,k,ie) = trhs(i,j,k,ie) - ta2(i,j,k,ie)
   195: ||||*--                 end do
   196: |||*---               end do
   197: ||V----             end do
   198: |P-----           end do
   199: |         c$OMP END DO
   200: |         c$OMP END PARALLEL
   201: |         c.......get the residual on mortar
   202: |                 call transfb(rmor,trhs)
   203: |                 if (timeron) call timer_stop(t_transf2)
   204: |         
   205: |         c.......apply boundary condition: zero out the residual on domain boundaries
   206: |         
   207: |         c.......apply boundary conidtion to trhs
   208: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ie,iside)
   209: |P----->          do ie=1,nelt
   210: ||+---->            do iside=1,nsides
   211: |||                   if (cbc(iside,ie).eq.0) then
   212: |||                     call facev(trhs(1,1,1,ie),iside,0.d0)
   213: |||                   end if
   214: ||+----             end do
   215: |P-----           end do
   216: |         c$OMP END PARALLELDO
   217: |         c.......apply boundary condition to rmor
   218: |                 call col2(rmor,tmmor,nmor)
   219: |         
   220: |         c.......call the conjugate gradient iterative solver
   221: |                 call diffusion(ifmortar)
   222: |         
   223: |         c.......add convection and diffusion
   224: |                 if (timeron) call timer_start(t_add2)
   225: |                 call add2(ta1,t,ntot)
   226: |                 if (timeron) call timer_stop(t_add2)
   227: |         
   228: |          
   229: |         c.......perform mesh adaptation
   230: |                 time=time+dtime
   231: |                 if ((step.ne.0).and.(step/fre*fre .eq. step)) then
   232: |                    if (step .ne. niter) then
   233: |                      call adaptation(ifmortar,step)
   234: |                    end if
   235: |                 else
   236: |                   ifmortar = .false.
   237: |                 end if
   238: |                 nelt_tot = nelt_tot + dble(nelt)
   239: +------         end do
   240:           
   241:                 call timer_stop(1)
   242:                 tmax = timer_read(1)
   243:            
   244:                 call verify(class, verified)
   245:           
   246:           c.....compute millions of collocation points advanced per second.
   247:           c.....diffusion: nmxh advancements, convection: 1 advancement
   248:                 mflops = nelt_tot*dble(lx1*lx1*lx1*(nmxh+1))/(tmax*1.d6)
   249:           
   250:                 call print_results('UA', class, refine_max, 0, 0, niter,
   251:                &     tmax, mflops, '    coll. point advanced',
   252:                &     verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5,
   253:                &     cs6, '(none)')
   254:           
   255:           c---------------------------------------------------------------------
   256:           c      More timers
   257:           c---------------------------------------------------------------------
   258:                 if (.not.timeron) goto 999
   259:           
   260: +------>        do i=1, t_last
   261: |                  trecs(i) = timer_read(i)
   262: +------         end do
   263:                 if (tmax .eq. 0.0) tmax = 1.0
   264:           
   265:                 write(*,800)
   266:            800  format('  SECTION     Time (secs)')
   267: +------>        do i=1, t_last
   268: |                  write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
   269: |                  if (i.eq.t_transfb_c) then
   270: |                     t2 = trecs(t_convect) - trecs(t_transfb_c)
   271: |                     write(*,820) 'sub-convect ', t2, t2*100./tmax
   272: |                  else if (i.eq.t_transfb) then
   273: |                     t2 = trecs(t_diffusion) - trecs(t_transf) - trecs(t_transfb)
   274: |                     write(*,820) 'sub-diffuse ', t2, t2*100./tmax
   275: |                  endif
   276: |          810     format(2x,a15,' :',f9.3,'  ( ',f6.2,' %)')
   277: |          820     format(2x,a15,':',f9.3,'  ( ',f6.2,' %)')
   278: +------         end do
   279:           
   280:            999  continue
   281:           
   282:                 end


