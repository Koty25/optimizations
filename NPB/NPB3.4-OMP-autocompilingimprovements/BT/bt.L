NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME : bt.f

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
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
        -fignore-induction-variable-overflow     : enable
        -fignore-volatile                        : enable
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
        -fmove-loop-invariants-unsafe            : enable
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
        -freplace-loop-equation                  : enable
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
        -finline-copy-arguments                  : disable
        -finline-functions                       : enable
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME: bt.f

PROCEDURE NAME: BT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    72: inl(1212): Source for routine not found.: CHECK_TIMER_FLAG
   128: vec( 103): Unvectorized loop.
   128: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   129: inl(1212): Source for routine not found.: TIMER_CLEAR
   132: inl(1212): Source for routine not found.: ALLOC_SPACE
   134: inl(1212): Source for routine not found.: SET_CONSTANTS
   136: inl(1212): Source for routine not found.: INITIALIZE
   138: inl(1212): Source for routine not found.: EXACT_RHS
   143: inl(1212): Source for routine not found.: ADI
   146: vec( 103): Unvectorized loop.
   146: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   149: inl(1212): Source for routine not found.: TIMER_START
   151: opt(3008): Reference within a conditional branch moved outside loop.
   151: vec( 110): Vectorization obstructive procedure reference.: ADI
   151: vec( 180): I/O statement obstructs vectorization.
   151: vec( 103): Unvectorized loop.
   155: opt(1118): This I/O statement inhibits optimization of loop.
   163: inl(1212): Source for routine not found.: TIMER_STOP
   164: inl(1212): Source for routine not found.: TIMER_READ
   166: inl(1212): Source for routine not found.: VERIFY
   177: inl(1212): Source for routine not found.: PRINT_RESULTS
   188: vec( 103): Unvectorized loop.
   188: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   195: opt(3008): Reference within a conditional branch moved outside loop.
   195: vec( 103): Unvectorized loop.
   195: vec( 180): I/O statement obstructs vectorization.
   196: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME: bt.f

PROCEDURE NAME: BT
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

    48:                  program BT
    49:           c---------------------------------------------------------------------
    50:           
    51:                  use bt_data
    52:                  implicit none
    53:           
    54:                  include 'blk_par.h'
    55:           
    56:                  integer i, niter, step, fstatus
    57:                  double precision navg, mflops, n3
    58:           
    59:                  external timer_read
    60:                  double precision tmax, timer_read, t, trecs(t_last)
    61:                  logical verified
    62:                  character class
    63:                  character t_names(t_last)*8
    64:           !$     integer  omp_get_max_threads
    65:           !$     external omp_get_max_threads
    66:           
    67:           c---------------------------------------------------------------------
    68:           c      Root node reads input file (if it exists) else takes
    69:           c      defaults from parameters
    70:           c---------------------------------------------------------------------
    71:           
    72:                  call check_timer_flag( timeron )
    73:                  if (timeron) then
    74:                    t_names(t_total) = 'total'
    75:                    t_names(t_rhsx) = 'rhsx'
    76:                    t_names(t_rhsy) = 'rhsy'
    77:                    t_names(t_rhsz) = 'rhsz'
    78:                    t_names(t_rhs) = 'rhs'
    79:                    t_names(t_xsolve) = 'xsolve'
    80:                    t_names(t_ysolve) = 'ysolve'
    81:                    t_names(t_zsolve) = 'zsolve'
    82:                    t_names(t_rdis1) = 'redist1'
    83:                    t_names(t_solsub) = 'solsubs'
    84:                    t_names(t_add) = 'add'
    85:                  endif
    86:           
    87:                  write(*, 1000)
    88:                  open (unit=2,file='inputbt.data',status='old', iostat=fstatus)
    89:           
    90:                  if (fstatus .eq. 0) then
    91:                    write(*,233)
    92:            233     format(' Reading from input file inputbt.data')
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
   105:            234   format(' No input file inputbt.data. Using compiled defaults')
   106:           
   107:                  write(*, 1001) grid_points(1), grid_points(2), grid_points(3)
   108:                  write(*, 1002) niter, dt
   109:                  if (blkdim .gt. 0) write(*, 1004) blkdim
   110:           !$     write(*, 1003) omp_get_max_threads()
   111:                  write(*, *)
   112:           
   113:            1000  format(//, ' NAS Parallel Benchmarks (NPB3.4-OMP)',
   114:                >            ' - BT Benchmark', /)
   115:            1001  format(' Size: ', i4, 'x', i4, 'x', i4)
   116:            1002  format(' Iterations: ', i4, '       dt: ', f11.7)
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
   127:           
   128: +------>         do i = 1, t_last
   129: |                   call timer_clear(i)
   130: +------          end do
   131:           
   132:                  call alloc_space
   133:           
   134:                  call set_constants
   135:           
   136:                  call initialize
   137:           
   138:                  call exact_rhs
   139:           
   140:           c---------------------------------------------------------------------
   141:           c      do one time step to touch all code, and reinitialize
   142:           c---------------------------------------------------------------------
   143:                  call adi
   144:                  call initialize
   145:           
   146: +------>         do i = 1, t_last
   147: |                   call timer_clear(i)
   148: +------          end do
   149:                  call timer_start(1)
   150:           
   151: +------>         do  step = 1, niter
   152: |         
   153: |                   if (mod(step, 20) .eq. 0 .or.
   154: |              >        step .eq. 1) then
   155: |                      write(*, 200) step
   156: |          200         format(' Time step ', i4)
   157: |                   endif
   158: |         
   159: |                   call adi
   160: |         
   161: +------          end do
   162:           
   163:                  call timer_stop(1)
   164:                  tmax = timer_read(1)
   165:            
   166:                  call verify(niter, class, verified)
   167:           
   168:                  n3 = dble(grid_points(1))*grid_points(2)*grid_points(3)
   169:                  navg = (grid_points(1)+grid_points(2)+grid_points(3))/3.d0
   170:                  if( tmax .ne. 0. ) then
   171:                     mflops = 1.0d-6*dble(niter)*
   172:                >            (3478.8*n3-17655.7*navg**2+28023.7*navg)
   173:                >            / tmax
   174:                  else
   175:                     mflops = 0.d0
   176:                  endif
   177:                  call print_results('BT', class, grid_points(1),
   178:                >  grid_points(2), grid_points(3), niter,
   179:                >  tmax, mflops, '          floating point',
   180:                >  verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5,
   181:                >  cs6, '(none)')
   182:           
   183:           c---------------------------------------------------------------------
   184:           c      More timers
   185:           c---------------------------------------------------------------------
   186:                  if (.not.timeron) goto 999
   187:           
   188: +------>         do i=1, t_last
   189: |                   trecs(i) = timer_read(i)
   190: +------          end do
   191:                  if (tmax .eq. 0.0) tmax = 1.0
   192:           
   193:                  write(*,800)
   194:            800   format('  SECTION   Time (secs)')
   195: +------>         do i=1, t_last
   196: |                   write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
   197: |                   if (i.eq.t_rhs) then
   198: |                      t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
   199: |                      write(*,820) 'sub-rhs ', t, t*100./tmax
   200: |                      t = trecs(t_rhs) - t
   201: |                      write(*,820) 'rest-rhs ', t, t*100./tmax
   202: |                   elseif (i.eq.t_solsub) then
   203: |                      t = trecs(t_xsolve) + trecs(t_ysolve) + trecs(t_zsolve)
   204: |              >           - trecs(t_rdis1) - trecs(t_solsub)
   205: |                      write(*,820) 'rest-sol ', t, t*100./tmax
   206: |                   endif
   207: |          810      format(2x,a15,':',f9.3,'  ( ',f6.2,' %)')
   208: |          820      format(2x,a15,':',f9.3,'  ( ',f6.2,' %)')
   209: +------          end do
   210:           
   211:            999   continue
   212:           
   213:                  end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME: bt.f

PROCEDURE NAME: BT
INLINE LIST

  ROOT: BT (bt.f:48)
  -> NOINLINE: CHECK_TIMER_FLAG (bt.f:72)
     *** Source for routine not found.
  -> NOINLINE: TIMER_CLEAR (bt.f:129)
     *** Source for routine not found.
  -> NOINLINE: ALLOC_SPACE (bt.f:132)
     *** Source for routine not found.
  -> NOINLINE: SET_CONSTANTS (bt.f:134)
     *** Source for routine not found.
  -> NOINLINE: INITIALIZE (bt.f:136)
     *** Source for routine not found.
  -> NOINLINE: EXACT_RHS (bt.f:138)
     *** Source for routine not found.
  -> NOINLINE: ADI (bt.f:143)
     *** Source for routine not found.
  -> NOINLINE: INITIALIZE (bt.f:144)
     *** Source for routine not found.
  -> NOINLINE: TIMER_CLEAR (bt.f:147)
     *** Source for routine not found.
  -> NOINLINE: TIMER_START (bt.f:149)
     *** Source for routine not found.
  -> NOINLINE: ADI (bt.f:159)
     *** Source for routine not found.
  -> NOINLINE: TIMER_STOP (bt.f:163)
     *** Source for routine not found.
  -> NOINLINE: TIMER_READ (bt.f:164)
     *** Source for routine not found.
  -> NOINLINE: VERIFY (bt.f:166)
     *** Source for routine not found.
  -> NOINLINE: PRINT_RESULTS (bt.f:177)
     *** Source for routine not found.
  -> NOINLINE: TIMER_READ (bt.f:189)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME: bt.f

PROCEDURE NAME: BT
VECTORIZATION LIST

  LOOP BEGIN: (bt.f:128)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (bt.f:128)
  LOOP END

  LOOP BEGIN: (bt.f:146)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (bt.f:146)
  LOOP END

  LOOP BEGIN: (bt.f:151)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : ADI (bt.f:151)
    *** I/O statement obstructs vectorization. (bt.f:151)
  LOOP END

  LOOP BEGIN: (bt.f:188)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_READ (bt.f:188)
  LOOP END

  LOOP BEGIN: (bt.f:195)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (bt.f:195)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:44 2021
FILE NAME: bt.f

PROCEDURE NAME: BT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 48 [s0-s12 s15-s16 s18-s37 s51-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1904 bytes
      Register spill area      :  112 bytes
      Parameter area           :  248 bytes
      Register save area       :  176 bytes
      User data area           : 1360 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (bt.f:128)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (bt.f:146)
    *** Estimated execution cycle                       : 10
  LOOP END

  LOOP BEGIN: (bt.f:151)
    *** Estimated execution cycle                       : 35
    *** The number of SCALAR REGISTER TRANSFER          : 20
  LOOP END

  LOOP BEGIN: (bt.f:188)
    *** Estimated execution cycle                       : 13
  LOOP END

  LOOP BEGIN: (bt.f:195)
    *** Estimated execution cycle                       : 118
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 7
            Across calls                                : 6
            Others                                      : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 18
            Across calls                                : 10
            Over basic blocks                           :  8
    *** The number of SCALAR REGISTER TRANSFER          : 38
  LOOP END


