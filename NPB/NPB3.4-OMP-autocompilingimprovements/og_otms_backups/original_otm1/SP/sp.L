NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME : sp.f

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
FILE NAME: sp.f

PROCEDURE NAME: SP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   131: vec( 103): Unvectorized loop.
   131: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   149: vec( 103): Unvectorized loop.
   149: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   154: opt(3008): Reference within a conditional branch moved outside loop.
   154: vec( 110): Vectorization obstructive procedure reference.: ADI
   154: vec( 180): I/O statement obstructs vectorization.
   154: vec( 103): Unvectorized loop.
   157: opt(1118): This I/O statement inhibits optimization of loop.
   192: vec( 103): Unvectorized loop.
   192: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   200: opt(3008): Reference within a conditional branch moved outside loop.
   200: vec( 103): Unvectorized loop.
   200: vec( 180): I/O statement obstructs vectorization.
   201: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: sp.f

PROCEDURE NAME: SP
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

    46:                  program SP
    47:           c---------------------------------------------------------------------
    48:           
    49:                  use sp_data
    50:                  implicit none
    51:           
    52:                  include 'blk_par.h'
    53:           
    54:                  integer          i, niter, step, fstatus
    55:                  external         timer_read
    56:                  double precision mflops, n3, t, tmax, timer_read, trecs(t_last)
    57:                  logical          verified
    58:                  character        class
    59:                  character        t_names(t_last)*8
    60:           !$     integer  omp_get_max_threads
    61:           !$     external omp_get_max_threads
    62:           
    63:           c---------------------------------------------------------------------
    64:           c      Read input file (if it exists), else take
    65:           c      defaults from parameters
    66:           c---------------------------------------------------------------------
    67:            
    68:                  call check_timer_flag( timeron )
    69:                  if (timeron) then
    70:                    t_names(t_total) = 'total'
    71:                    t_names(t_rhsx) = 'rhsx'
    72:                    t_names(t_rhsy) = 'rhsy'
    73:                    t_names(t_rhsz) = 'rhsz'
    74:                    t_names(t_rhs) = 'rhs'
    75:                    t_names(t_xsolve) = 'xsolve'
    76:                    t_names(t_ysolve) = 'ysolve'
    77:                    t_names(t_zsolve) = 'zsolve'
    78:                    t_names(t_rdis1) = 'redist1'
    79:                    t_names(t_rdis2) = 'redist2'
    80:                    t_names(t_tzetar) = 'tzetar'
    81:                    t_names(t_ninvr) = 'ninvr'
    82:                    t_names(t_pinvr) = 'pinvr'
    83:                    t_names(t_txinvr) = 'txinvr'
    84:                    t_names(t_add) = 'add'
    85:                  endif
    86:           
    87:                  write(*, 1000)
    88:                  open (unit=2,file='inputsp.data',status='old', iostat=fstatus)
    89:           
    90:                  if (fstatus .eq. 0) then
    91:                    write(*,233)
    92:            233     format(' Reading from input file inputsp.data')
    93:                    read (2,*) niter
    94:                    read (2,*) dt
    95:                    read (2,*) grid_points(1), grid_points(2), grid_points(3)
    96:                    close(2)
    97:                  else
    98:                    write(*,234)
    99:                    niter = niter_default
   100:                    dt    = dt_default
   101:                    grid_points(1) = problem_size
   102:                    grid_points(2) = problem_size
   103:                    grid_points(3) = problem_size
   104:                  endif
   105:            234   format(' No input file inputsp.data. Using compiled defaults')
   106:           
   107:                  write(*, 1001) grid_points(1), grid_points(2), grid_points(3)
   108:                  write(*, 1002) niter, dt
   109:                  if (blkdim .gt. 0) write(*, 1004) blkdim
   110:           !$     write(*, 1003) omp_get_max_threads()
   111:                  write(*, *)
   112:           
   113:            1000  format(//, ' NAS Parallel Benchmarks (NPB3.4-OMP)',
   114:                >            ' - SP Benchmark', /)
   115:            1001  format(' Size: ', i4, 'x', i4, 'x', i4)
   116:            1002  format(' Iterations: ', i4, '    dt:  ', f11.7)
   117:            1003  format(' Number of available threads: ', i5)
   118:            1004  format(' Dimension blocking size: ', i5)
   119:           
   120:                  if ( (grid_points(1) .gt. IMAX) .or.
   121:                >      (grid_points(2) .gt. JMAX) .or.
   122:                >      (grid_points(3) .gt. KMAX) ) then
   123: +------>               print *, (grid_points(i),i=1,3)
   124:                        print *,' Problem size too big for compiled array sizes'
   125:                        goto 999
   126:                  endif
   127:                  nx2 = grid_points(1) - 2
   128:                  ny2 = grid_points(2) - 2
   129:                  nz2 = grid_points(3) - 2
   130:           
   131: +------>         do i = 1, t_last
   132: |                   call timer_clear(i)
   133: +------          end do
   134:           
   135:                  call alloc_space
   136:           
   137:                  call set_constants
   138:           
   139:                  call exact_rhs
   140:           
   141:                  call initialize
   142:           
   143:           c---------------------------------------------------------------------
   144:           c      do one time step to touch all code, and reinitialize
   145:           c---------------------------------------------------------------------
   146:                  call adi
   147:                  call initialize
   148:           
   149: +------>         do i = 1, t_last
   150: |                   call timer_clear(i)
   151: +------          end do
   152:                  call timer_start(1)
   153:           
   154: +------>         do  step = 1, niter
   155: |         
   156: |                   if (mod(step, 20) .eq. 0 .or. step .eq. 1) then
   157: |                      write(*, 200) step
   158: |          200         format(' Time step ', i4)
   159: |                   endif
   160: |         
   161: |                   call adi
   162: |         
   163: +------          end do
   164:           
   165:                  call timer_stop(1)
   166:                  tmax = timer_read(1)
   167:            
   168:                  call verify(niter, class, verified)
   169:           
   170:                  if( tmax .ne. 0. ) then
   171:                     n3 = dble(grid_points(1))*grid_points(2)*grid_points(3)
   172:                     t = (grid_points(1)+grid_points(2)+grid_points(3))/3.d0
   173:                     mflops = 1.0d-6*dble( niter )*(881.174 * n3
   174:                >             -4683.91 * t**2
   175:                >             +11484.5 * t
   176:                >             -19272.4) / tmax
   177:                  else
   178:                     mflops = 0.d0
   179:                  endif
   180:           
   181:                 call print_results('SP', class, grid_points(1),
   182:                >     grid_points(2), grid_points(3), niter,
   183:                >     tmax, mflops, '          floating point',
   184:                >     verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5,
   185:                >     cs6, '(none)')
   186:           
   187:           c---------------------------------------------------------------------
   188:           c      More timers
   189:           c---------------------------------------------------------------------
   190:                  if (.not.timeron) goto 999
   191:           
   192: +------>         do i=1, t_last
   193: |                   trecs(i) = timer_read(i)
   194: +------          end do
   195:                  if (tmax .eq. 0.0) tmax = 1.0
   196:           
   197:                  write(*,800)
   198:            800   format('  SECTION   Time (secs)')
   199:           
   200: +------>         do i=1, t_last
   201: |                   write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
   202: |                   if (i.eq.t_rhs) then
   203: |                      t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
   204: |                      write(*,820) 'sub-rhs ', t, t*100./tmax
   205: |                      t = trecs(t_rhs) - t
   206: |                      write(*,820) 'rest-rhs ', t, t*100./tmax
   207: |                   elseif (i.eq.t_zsolve) then
   208: |                      t = trecs(t_zsolve) - trecs(t_rdis1) - trecs(t_rdis2)
   209: |                      write(*,820) 'sub-zsol ', t, t*100./tmax
   210: |                   elseif (i.eq.t_rdis2) then
   211: |                      t = trecs(t_rdis1) + trecs(t_rdis2)
   212: |                      write(*,820) 'redist ', t, t*100./tmax
   213: |                   endif
   214: |          810      format(2x,a12,':',f9.3,'  ( ',f6.2,' %)')
   215: |          820      format(2x,a12,':',f9.3,'  ( ',f6.2,' %)')
   216: +------          end do
   217:           
   218:            999   continue
   219:           
   220:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: sp.f

PROCEDURE NAME: SP
VECTORIZATION LIST

  LOOP BEGIN: (sp.f:131)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (sp.f:131)
  LOOP END

  LOOP BEGIN: (sp.f:149)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (sp.f:149)
  LOOP END

  LOOP BEGIN: (sp.f:154)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : ADI (sp.f:154)
    *** I/O statement obstructs vectorization. (sp.f:154)
  LOOP END

  LOOP BEGIN: (sp.f:192)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_READ (sp.f:192)
  LOOP END

  LOOP BEGIN: (sp.f:200)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (sp.f:200)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:03 2021
FILE NAME: sp.f

PROCEDURE NAME: SP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 51 [s0-s12 s15-s16 s18-s37 s48-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 2032 bytes
      Register spill area      :  144 bytes
      Parameter area           :  248 bytes
      Register save area       :  176 bytes
      User data area           : 1456 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (sp.f:131)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (sp.f:149)
    *** Estimated execution cycle                       : 10
  LOOP END

  LOOP BEGIN: (sp.f:154)
    *** Estimated execution cycle                       : 35
    *** The number of SCALAR REGISTER TRANSFER          : 20
  LOOP END

  LOOP BEGIN: (sp.f:192)
    *** Estimated execution cycle                       : 13
  LOOP END

  LOOP BEGIN: (sp.f:200)
    *** Estimated execution cycle                       : 147
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 8
            Across calls                                : 7
            Others                                      : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 26
            Across calls                                : 12
            Over basic blocks                           : 14
    *** The number of SCALAR REGISTER TRANSFER          : 45
  LOOP END


