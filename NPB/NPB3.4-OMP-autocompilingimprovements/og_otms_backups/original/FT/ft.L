NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME : ft.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    82: vec( 103): Unvectorized loop.
    82: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
    99: vec( 103): Unvectorized loop.
    99: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
   117: vec( 110): Vectorization obstructive procedure reference.: EVOLVE
   117: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
   117: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
   117: vec( 110): Vectorization obstructive procedure reference.: FFT
   117: vec( 110): Vectorization obstructive procedure reference.: CHECKSUM
   117: vec( 103): Unvectorized loop.
   118: opt(1025): Reference to this procedure inhibits optimization.: TIMER_START


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FT
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

    55:                 program ft
    56:           
    57:           c---------------------------------------------------------------------
    58:           c---------------------------------------------------------------------
    59:           
    60:           c---------------------------------------------------------------------
    61:           c Module ft_fields defines main arrays (u0, u1, u2) in the problem
    62:           c---------------------------------------------------------------------
    63:           
    64:                 use ft_data
    65:                 use ft_fields
    66:           
    67:                 implicit none
    68:           
    69:                 integer i
    70:           
    71:                 integer iter
    72:                 double precision total_time, mflops
    73:                 logical verified
    74:                 character class
    75:           
    76:           
    77:           c---------------------------------------------------------------------
    78:           c Run the entire problem once to make sure all data is touched.
    79:           c This reduces variable startup costs, which is important for such a
    80:           c short benchmark. The other NPB 2 implementations are similar.
    81:           c---------------------------------------------------------------------
    82: +------>        do i = 1, t_max
    83: |                  call timer_clear(i)
    84: +------         end do
    85:           
    86:                 call alloc_space
    87:           
    88:                 call setup()
    89:                 call init_ui(u0, u1, twiddle, dims(1), dims(2), dims(3))
    90:                 call compute_indexmap(twiddle, dims(1), dims(2), dims(3))
    91:                 call compute_initial_conditions(u1, dims(1), dims(2), dims(3))
    92:                 call fft_init (dims(1))
    93:                 call fft(1, u1, u0)
    94:           
    95:           c---------------------------------------------------------------------
    96:           c Start over from the beginning. Note that all operations must
    97:           c be timed, in contrast to other benchmarks.
    98:           c---------------------------------------------------------------------
    99: +------>        do i = 1, t_max
   100: |                  call timer_clear(i)
   101: +------         end do
   102:           
   103:                 call timer_start(T_total)
   104:                 if (timers_enabled) call timer_start(T_setup)
   105:           
   106:                 call compute_indexmap(twiddle, dims(1), dims(2), dims(3))
   107:           
   108:                 call compute_initial_conditions(u1, dims(1), dims(2), dims(3))
   109:           
   110:                 call fft_init (dims(1))
   111:           
   112:                 if (timers_enabled) call timer_stop(T_setup)
   113:                 if (timers_enabled) call timer_start(T_fft)
   114:                 call fft(1, u1, u0)
   115:                 if (timers_enabled) call timer_stop(T_fft)
   116:           
   117: +------>        do iter = 1, niter
   118: |                  if (timers_enabled) call timer_start(T_evolve)
   119: |                  call evolve(u0, u1, twiddle, dims(1), dims(2), dims(3))
   120: |                  if (timers_enabled) call timer_stop(T_evolve)
   121: |                  if (timers_enabled) call timer_start(T_fft)
   122: |         c         call fft(-1, u1, u2)
   123: |                  call fft(-1, u1, u1)
   124: |                  if (timers_enabled) call timer_stop(T_fft)
   125: |                  if (timers_enabled) call timer_start(T_checksum)
   126: |         c         call checksum(iter, u2, dims(1), dims(2), dims(3))
   127: |                  call checksum(iter, u1, dims(1), dims(2), dims(3))
   128: |                  if (timers_enabled) call timer_stop(T_checksum)
   129: +------         end do
   130:           
   131:                 call verify(nx, ny, nz, niter, verified, class)
   132:           
   133:                 call timer_stop(t_total)
   134:                 total_time = timer_read(t_total)
   135:           
   136:                 if( total_time .ne. 0. ) then
   137:                    mflops = 1.0d-6*ntotal_f *
   138:                >             (14.8157+7.19641*log(ntotal_f)
   139:                >          +  (5.23518+7.21113*log(ntotal_f))*niter)
   140:                >                 /total_time
   141:                 else
   142:                    mflops = 0.0
   143:                 endif
   144:                 call print_results('FT', class, nx, ny, nz, niter,
   145:                >  total_time, mflops, '          floating point', verified,
   146:                >  npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
   147:                 if (timers_enabled) call print_timers()
   148:           
   149:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FT
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:82)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (ft.f:82)
  LOOP END

  LOOP BEGIN: (ft.f:99)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_CLEAR (ft.f:99)
  LOOP END

  LOOP BEGIN: (ft.f:117)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : EVOLVE (ft.f:117)
    *** Vectorization obstructive procedure reference. : TIMER_START (ft.f:117)
    *** Vectorization obstructive procedure reference. : TIMER_STOP (ft.f:117)
    *** Vectorization obstructive procedure reference. : FFT (ft.f:117)
    *** Vectorization obstructive procedure reference. : CHECKSUM (ft.f:117)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 44 [s0-s12 s15-s16 s18-s34 s52-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 656 bytes
      Register spill area      :  48 bytes
      Parameter area           : 248 bytes
      Register save area       : 176 bytes
      User data area           : 176 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:82)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (ft.f:99)
    *** Estimated execution cycle                       : 10
  LOOP END

  LOOP BEGIN: (ft.f:117)
    *** Estimated execution cycle                       : 15
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 5
            Over basic blocks                           : 5
    *** The number of SCALAR REGISTER TRANSFER          : 17
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: INIT_UI
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   170: par(1801): Parallel routine generated.: INIT_UI$1
   171: par(1803): Parallelized by "do".
   171: opt(2495): Loop nest collapsed by "collapse" clause.
   171: vec( 101): Vectorized loop.
   173: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: INIT_UI
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

   154:                 subroutine init_ui(u0, u1, twiddle, d1, d2, d3)
   155:           
   156:           c---------------------------------------------------------------------
   157:           c---------------------------------------------------------------------
   158:           
   159:           c---------------------------------------------------------------------
   160:           c touch all the big data
   161:           c---------------------------------------------------------------------
   162:           
   163:                 implicit none
   164:                 integer d1, d2, d3
   165:                 double complex   u0(d1+1,d2,d3)
   166:                 double complex   u1(d1+1,d2,d3)
   167:                 double precision twiddle(d1+1,d2,d3)
   168:                 integer i, j, k
   169:           
   170:           !$omp parallel do default(shared) private(i,j,k) collapse(2)
   171: Y------>        do k = 1, d3
   172: |+----->           do j = 1, d2
   173: ||V---->              do i = 1, d1
   174: |||                      u0(i,j,k) = 0.d0
   175: |||                      u1(i,j,k) = 0.d0
   176: |||                      twiddle(i,j,k) = 0.d0
   177: ||V----               end do
   178: |+-----            end do
   179: Y------         end do
   180:           
   181:                 return
   182:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: INIT_UI
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:171)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ft.f:171)
    *** The number of VLOAD, VSTORE. :  0,  2. (ft.f:171)

    LOOP BEGIN: (ft.f:173)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:173)
      *** The number of VLOAD, VSTORE. :  0,  5. (ft.f:173)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: INIT_UI
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 47 [s0-s2 s8-s12 s15-s16 s18 s23-s28 s34-s63]
      Vector registers         : 18 [v46-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 304 bytes
      Register spill area      :   8 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           : 104 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:171)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (ft.f:171)
    *** Estimated execution cycle                       : 32
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (ft.f:173)
      *** Estimated execution cycle                     : 169
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: EVOLVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   205: par(1801): Parallel routine generated.: EVOLVE$1
   206: par(1803): Parallelized by "do".
   206: opt(2495): Loop nest collapsed by "collapse" clause.
   206: vec( 101): Vectorized loop.
   208: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: EVOLVE
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

   187:                 subroutine evolve(u0, u1, twiddle, d1, d2, d3)
   188:           
   189:           c---------------------------------------------------------------------
   190:           c---------------------------------------------------------------------
   191:           
   192:           c---------------------------------------------------------------------
   193:           c evolve u0 -> u1 (t time steps) in fourier space
   194:           c---------------------------------------------------------------------
   195:           
   196:                 use ft_data
   197:                 implicit none
   198:           
   199:                 integer d1, d2, d3
   200:                 double complex   u0(d1+1,d2,d3)
   201:                 double complex   u1(d1+1,d2,d3)
   202:                 double precision twiddle(d1+1,d2,d3)
   203:                 integer i, j, k
   204:           
   205:           !$omp parallel do default(shared) private(i,j,k) collapse(2)
   206: Y------>        do k = 1, d3
   207: |+----->           do j = 1, d2
   208: ||V---->              do i = 1, d1
   209: |||                      u0(i,j,k) = u0(i,j,k) * twiddle(i,j,k)
   210: |||                      u1(i,j,k) = u0(i,j,k)
   211: ||V----               end do
   212: |+-----            end do
   213: Y------         end do
   214:           
   215:                 return
   216:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: EVOLVE
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:206)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ft.f:206)
    *** The number of VLOAD, VSTORE. :  0,  2. (ft.f:206)

    LOOP BEGIN: (ft.f:208)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:208)
      *** The number of VLOAD, VSTORE. :  3,  4. (ft.f:208)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: EVOLVE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 46 [s0-s2 s8-s12 s15-s16 s18 s23-s28 s35-s63]
      Vector registers         : 18 [v46-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 304 bytes
      Register spill area      :   8 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           : 104 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:206)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (ft.f:206)
    *** Estimated execution cycle                       : 32
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (ft.f:208)
      *** Estimated execution cycle                     : 234
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INITIAL_CONDITIONS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   250: vec( 110): Vectorization obstructive procedure reference.: RANDLC
   250: vec( 103): Unvectorized loop.
   251: opt(1025): Reference to this procedure inhibits optimization.: RANDLC
   258: par(1801): Parallel routine generated.: COMPUTE_INITIAL_CONDITIONS$1
   259: par(1803): Parallelized by "do".
   261: vec( 110): Vectorization obstructive procedure reference.: VRANLC
   261: vec( 103): Unvectorized loop.
   262: opt(1025): Reference to this procedure inhibits optimization.: VRANLC


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INITIAL_CONDITIONS
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

   222:                 subroutine compute_initial_conditions(u0, d1, d2, d3)
   223:           
   224:           c---------------------------------------------------------------------
   225:           c---------------------------------------------------------------------
   226:           
   227:           c---------------------------------------------------------------------
   228:           c Fill in array u0 with initial conditions from
   229:           c random number generator
   230:           c---------------------------------------------------------------------
   231:           
   232:                 use ft_data
   233:                 implicit none
   234:           
   235:                 integer d1, d2, d3
   236:                 double complex u0(d1+1, d2, d3)
   237:                 integer k, j
   238:                 double precision x0, start, an, dummy, starts(nz)
   239:            
   240:           
   241:                 start = seed
   242:           c---------------------------------------------------------------------
   243:           c Jump to the starting element for our first plane.
   244:           c---------------------------------------------------------------------
   245:                 call ipow46(a, 0, an)
   246:                 dummy = randlc(start, an)
   247:                 call ipow46(a, 2*nx*ny, an)
   248:           
   249:                 starts(1) = start
   250: +------>        do k = 2, dims(3)
   251: |                  dummy = randlc(start, an)
   252: |                  starts(k) = start
   253: +------         end do
   254:            
   255:           c---------------------------------------------------------------------
   256:           c Go through by z planes filling in one square at a time.
   257:           c---------------------------------------------------------------------
   258:           !$omp parallel do default(shared) private(k,j,x0)
   259: P------>        do k = 1, dims(3)
   260: |                  x0 = starts(k)
   261: |+----->           do j = 1, dims(2)
   262: ||                    call vranlc(2*nx, x0, a, u0(1, j, k))
   263: |+-----            end do
   264: P------         end do
   265:           
   266:                 return
   267:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INITIAL_CONDITIONS
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:259)
    <Unvectorized loop.>

    LOOP BEGIN: (ft.f:261)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : VRANLC (ft.f:261)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INITIAL_CONDITIONS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 35 [s0-s3 s8-s12 s15-s16 s18-s30 s53-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 272 bytes
      Register spill area      :   0 bytes
      Parameter area           :  32 bytes
      Register save area       : 176 bytes
      User data area           :  64 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:259)
    *** Estimated execution cycle                       : 10
    *** The number of SCALAR REGISTER TRANSFER          : 7

    LOOP BEGIN: (ft.f:261)
      *** Estimated execution cycle                     : 20
      *** The number of SCALAR REGISTER TRANSFER        : 14
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: IPOW46
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   299: vec( 103): Unvectorized loop.
   299: vec( 110): Vectorization obstructive procedure reference.: RANDLC


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: IPOW46
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

   273:                 subroutine ipow46(a, exponent, result)
   274:           
   275:           c---------------------------------------------------------------------
   276:           c---------------------------------------------------------------------
   277:           
   278:           c---------------------------------------------------------------------
   279:           c compute a^exponent mod 2^46
   280:           c---------------------------------------------------------------------
   281:           
   282:                 implicit none
   283:                 double precision a, result, dummy, q, r
   284:                 integer exponent, n, n2
   285:                 external randlc
   286:                 double precision randlc
   287:           c---------------------------------------------------------------------
   288:           c Use
   289:           c   a^n = a^(n/2)*a^(n/2) if n even else
   290:           c   a^n = a*a^(n-1)       if n odd
   291:           c---------------------------------------------------------------------
   292:                 result = 1
   293:                 if (exponent .eq. 0) return
   294:                 q = a
   295:                 r = 1
   296:                 n = exponent
   297:           
   298:           
   299: +------>        do while (n .gt. 1)
   300: |                  n2 = n/2
   301: |                  if (n2 * 2 .eq. n) then
   302: |                     dummy = randlc(q, q)
   303: |                     n = n2
   304: |                  else
   305: |                     dummy = randlc(r, q)
   306: |                     n = n-1
   307: |                  endif
   308: +------         end do
   309:                 dummy = randlc(r, q)
   310:                 result = r
   311:                 return
   312:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: IPOW46
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:299)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : RANDLC (ft.f:299)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: IPOW46
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 19 [s0-s2 s8-s12 s15-s16 s18-s20 s23-s25 s61-s63]
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

  LOOP BEGIN: (ft.f:299)
    *** Estimated execution cycle                       : 13
    *** The number of SCALAR REGISTER TRANSFER          : 6
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: SETUP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   371: opt(1317): Logical clause simplified.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

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

   318:                 subroutine setup
   319:           
   320:           c---------------------------------------------------------------------
   321:           c---------------------------------------------------------------------
   322:           
   323:                 use ft_data
   324:                 implicit none
   325:           
   326:           !$    integer  omp_get_max_threads
   327:           !$    external omp_get_max_threads
   328:                 debug = .FALSE.
   329:           
   330:                 call check_timer_flag( timers_enabled )
   331:           
   332:                 write(*, 1000)
   333:           
   334:                 niter = niter_default
   335:           
   336:                 write(*, 1001) nx, ny, nz
   337:                 write(*, 1002) niter
   338:           !$    write(*, 1003) omp_get_max_threads()
   339:                 write(*, *)
   340:           
   341:           
   342:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
   343:                >          ' - FT Benchmark', /)
   344:            1001 format(' Size                : ', i4, 'x', i4, 'x', i4)
   345:            1002 format(' Iterations                  :', i7)
   346:            1003 format(' Number of available threads :', i7)
   347:           
   348:                 dims(1) = nx
   349:                 dims(2) = ny
   350:                 dims(3) = nz
   351:           
   352:           
   353:           c---------------------------------------------------------------------
   354:           c Set up info for blocking of ffts and transposes.  This improves
   355:           c performance on cache-based systems. Blocking involves
   356:           c working on a chunk of the problem at a time, taking chunks
   357:           c along the first, second, or third dimension.
   358:           c
   359:           c - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
   360:           c - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)
   361:           
   362:           c Since 1st dim is always in processor, we'll assume it's long enough
   363:           c (default blocking factor is 16 so min size for 1st dim is 16)
   364:           c The only case we have to worry about is cffts1 in a 2d decomposition.
   365:           c so the blocking factor should not be larger than the 2nd dimension.
   366:           c---------------------------------------------------------------------
   367:           
   368:                 fftblock = fftblock_default
   369:                 fftblockpad = fftblockpad_default
   370:           
   371:                 if (fftblock .ne. fftblock_default) fftblockpad = fftblock+3
   372:           
   373:                 return
   374:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: SETUP
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: SETUP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 26 [s0-s12 s15-s16 s23-s30 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1280 bytes
      Register spill area      :    0 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1024 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INDEXMAP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   409: par(1801): Parallel routine generated.: COMPUTE_INDEXMAP$1
   411: par(1803): Parallelized by "do".
   411: opt(2495): Loop nest collapsed by "collapse" clause.
   411: vec( 101): Vectorized loop.
   417: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INDEXMAP
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

   380:                 subroutine compute_indexmap(twiddle, d1, d2, d3)
   381:           
   382:           c---------------------------------------------------------------------
   383:           c---------------------------------------------------------------------
   384:           
   385:           c---------------------------------------------------------------------
   386:           c compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2
   387:           c for time evolution exponent.
   388:           c---------------------------------------------------------------------
   389:           
   390:                 use ft_data
   391:                 implicit none
   392:           
   393:                 integer d1, d2, d3
   394:                 double precision twiddle(d1+1, d2, d3)
   395:                 integer i, j, k, kk, kk2, jj, kj2, ii
   396:                 double precision ap
   397:           
   398:           c---------------------------------------------------------------------
   399:           c basically we want to convert the fortran indices
   400:           c   1 2 3 4 5 6 7 8
   401:           c to
   402:           c   0 1 2 3 -4 -3 -2 -1
   403:           c The following magic formula does the trick:
   404:           c mod(i-1+n/2, n) - n/2
   405:           c---------------------------------------------------------------------
   406:           
   407:                 ap = - 4.d0 * alpha * pi *pi
   408:           
   409:           !$omp parallel do default(shared) private(i,j,k,kk,kk2,jj,kj2,ii)
   410:           !$omp&  collapse(2)
   411: Y------>        do k = 1, dims(3)
   412: |+----->           do j = 1, dims(2)
   413: ||                    kk =  mod(k-1+nz/2, nz) - nz/2
   414: ||                    kk2 = kk*kk
   415: ||                    jj = mod(j-1+ny/2, ny) - ny/2
   416: ||                    kj2 = jj*jj+kk2
   417: ||V---->              do i = 1, dims(1)
   418: |||                      ii = mod(i-1+nx/2, nx) - nx/2
   419: |||                      twiddle(i,j,k) = dexp(ap*dble(ii*ii+kj2))
   420: ||V----               end do
   421: |+-----            end do
   422: Y------         end do
   423:           
   424:                 return
   425:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INDEXMAP
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:411)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ft.f:411)
    *** The number of VLOAD, VSTORE. :  0,  3. (ft.f:411)

    LOOP BEGIN: (ft.f:417)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:417)
      *** The number of VLOAD, VSTORE. :  0,  1. (ft.f:417)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: COMPUTE_INDEXMAP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 51 [s0-s2 s8-s12 s15-s16 s18-s33 s39-s63]
      Vector registers         : 37 [v0 v28-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 4416 bytes
      Register spill area      : 2088 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 2136 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:411)
    *** Estimated execution cycle                       : 1111
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (ft.f:411)
    *** Estimated execution cycle                       : 51
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 4
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (ft.f:417)
      *** Estimated execution cycle                     : 438
      *** The number of SCALAR REGISTER TRANSFER        : 5
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: PRINT_TIMERS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   456: vec( 103): Unvectorized loop.
   456: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   456: vec( 180): I/O statement obstructs vectorization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: PRINT_TIMERS
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

   432:                 subroutine print_timers()
   433:           
   434:           c---------------------------------------------------------------------
   435:           c---------------------------------------------------------------------
   436:           
   437:                 use ft_data
   438:                 implicit none
   439:           
   440:                 integer i
   441:                 double precision t, t_m
   442:                 character*25 tstrings(T_max)
   443:                 data tstrings / '          total ',
   444:                >                '          setup ',
   445:                >                '            fft ',
   446:                >                '         evolve ',
   447:                >                '       checksum ',
   448:                >                '           fftx ',
   449:                >                '           ffty ',
   450:                >                '           fftz ' /
   451:           
   452:                 t_m = timer_read(T_total)
   453:                 if (t_m .le. 0.0d0) t_m = 1.0d0
   454:                  write(*,800)
   455:            800   format('  SECTION   Time (secs)')
   456: +------>        do i = 1, t_max
   457: |                  t = timer_read(i)
   458: |                  write(*, 100) tstrings(i), t, t*100.0/t_m
   459: +------         end do
   460:            100  format(2x,A16,':',f9.3,'  ( ',f6.2,' %)')
   461:                 return
   462:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: PRINT_TIMERS
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:456)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_READ (ft.f:456)
    *** I/O statement obstructs vectorization. (ft.f:456)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: PRINT_TIMERS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 31 [s0-s12 s15-s16 s23-s35 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1328 bytes
      Register spill area      :   16 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1056 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:456)
    *** Estimated execution cycle                       : 65
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Across calls                                : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 1
            Across calls                                : 1
    *** The number of SCALAR REGISTER TRANSFER          : 10
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT
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

   469:                 subroutine fft(dir, x1, x2)
   470:           
   471:           c---------------------------------------------------------------------
   472:           c---------------------------------------------------------------------
   473:           
   474:                 use ft_data
   475:                 implicit none
   476:           
   477:                 integer dir
   478:                 double complex x1(ntotalp), x2(ntotalp)
   479:           
   480:                 double complex y1(fftblockpad_default*maxdim),
   481:                >               y2(fftblockpad_default*maxdim)
   482:           
   483:           c---------------------------------------------------------------------
   484:           c note: args x1, x2 must be different arrays
   485:           c note: args for cfftsx are (direction, layout, xin, xout, scratch)
   486:           c       xin/xout may be the same and it can be somewhat faster
   487:           c       if they are
   488:           c---------------------------------------------------------------------
   489:           
   490:                 if (dir .eq. 1) then
   491:                    call cffts1(1, dims(1), dims(2), dims(3), x1, x1, y1, y2)
   492:                    call cffts2(1, dims(1), dims(2), dims(3), x1, x1, y1, y2)
   493:                    call cffts3(1, dims(1), dims(2), dims(3), x1, x2, y1, y2)
   494:                 else
   495:                    call cffts3(-1, dims(1), dims(2), dims(3), x1, x1, y1, y2)
   496:                    call cffts2(-1, dims(1), dims(2), dims(3), x1, x1, y1, y2)
   497:                    call cffts1(-1, dims(1), dims(2), dims(3), x1, x2, y1, y2)
   498:                 endif
   499:                 return
   500:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 26 [s0-s12 s15-s16 s23-s30 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 540944 bytes
      Register spill area      :      0 bytes
      Parameter area           :     64 bytes
      Register save area       :    176 bytes
      User data area           : 540704 bytes
      Others                   :      0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   524: par(1801): Parallel routine generated.: CFFTS1$1
   526: opt(3008): Reference within a conditional branch moved outside loop.
   526: vec( 110): Vectorization obstructive procedure reference.: CFFTZ
   526: par(1803): Parallelized by "do".
   526: opt(2495): Loop nest collapsed by "collapse" clause.
   530: opt(1592): Outer loop unrolled inside inner loop.: J
   531: vec( 101): Vectorized loop.
   539: opt(1592): Outer loop unrolled inside inner loop.: J
   540: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS1
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

   507:                 subroutine cffts1(is, d1, d2, d3, x, xout, y1, y2)
   508:           
   509:           c---------------------------------------------------------------------
   510:           c---------------------------------------------------------------------
   511:           
   512:                 use ft_data
   513:                 implicit none
   514:           
   515:                 integer is, d1, d2, d3, logd1
   516:                 double complex x(d1+1,d2,d3)
   517:                 double complex xout(d1+1,d2,d3)
   518:                 double complex y1(fftblockpad, d1), y2(fftblockpad, d1)
   519:                 integer i, j, k, jj, jn
   520:           
   521:                 logd1 = ilog2(d1)
   522:           
   523:                 if (timers_enabled) call timer_start(T_fftx)
   524:           !$omp parallel do default(shared) private(i,j,k,jj,y1,y2,jn)
   525:           !$omp&  shared(is,logd1,d1) collapse(2)
   526: P------>        do k = 1, d3
   527: |+----->           do jn = 0, d2/fftblock - 1
   528: ||        !         do jj = 0, d2 - fftblock, fftblock
   529: ||                    jj = jn*fftblock
   530: ||U---->              do j = 1, fftblock
   531: |||V--->                 do i = 1, d1
   532: ||||                        y1(j,i) = x(i,j+jj,k)
   533: |||V---                  enddo
   534: ||U----               enddo
   535: ||         
   536: ||                    call cfftz (is, logd1, d1, y1, y2)
   537: ||        
   538: ||        
   539: ||U---->              do j = 1, fftblock
   540: |||V--->                 do i = 1, d1
   541: ||||                        xout(i,j+jj,k) = y1(j,i)
   542: |||V---                  enddo
   543: ||U----               enddo
   544: |+-----            enddo
   545: P------         enddo
   546:                 if (timers_enabled) call timer_stop(T_fftx)
   547:           
   548:                 return
   549:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS1
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:526)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : CFFTZ (ft.f:526)

    LOOP BEGIN: (ft.f:530)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:530)

      LOOP BEGIN: (ft.f:531)
        <Vectorized loop.>
        **  Copied loop. (ft.f:531)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:531)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:531)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:530)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:530)

      LOOP BEGIN: (ft.f:531)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:531)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:531)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:539)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:539)

      LOOP BEGIN: (ft.f:540)
        <Vectorized loop.>
        **  Copied loop. (ft.f:540)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:540)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:540)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:539)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:539)

      LOOP BEGIN: (ft.f:540)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:540)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:540)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS1
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 640 bytes
      Register spill area      : 152 bytes
      Parameter area           :  40 bytes
      Register save area       : 176 bytes
      User data area           : 264 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:526)
    *** Estimated execution cycle                       : 59
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 24
            Over basic blocks                           : 24
    *** The number of SCALAR REGISTER TRANSFER          : 37

    LOOP BEGIN: (ft.f:530)
      *** Estimated execution cycle                     : 16
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 14

      LOOP BEGIN: (ft.f:531)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 4
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:530)
      *** Estimated execution cycle                     : 26
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 23

      LOOP BEGIN: (ft.f:531)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 1
                Over basic blocks                       : 1
        *** The number of SCALAR REGISTER TRANSFER      : 7
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:539)
      *** Estimated execution cycle                     : 17
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 9

      LOOP BEGIN: (ft.f:540)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 4
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:539)
      *** Estimated execution cycle                     : 25
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 13

      LOOP BEGIN: (ft.f:540)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 2
                Over basic blocks                       : 2
        *** The number of SCALAR REGISTER TRANSFER      : 6
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   572: par(1801): Parallel routine generated.: CFFTS2$1
   574: opt(3008): Reference within a conditional branch moved outside loop.
   574: vec( 110): Vectorization obstructive procedure reference.: CFFTZ
   574: par(1803): Parallelized by "do".
   574: opt(2495): Loop nest collapsed by "collapse" clause.
   578: opt(1592): Outer loop unrolled inside inner loop.: J
   579: vec( 101): Vectorized loop.
   586: opt(1592): Outer loop unrolled inside inner loop.: J
   587: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS2
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

   555:                 subroutine cffts2(is, d1, d2, d3, x, xout, y1, y2)
   556:           
   557:           c---------------------------------------------------------------------
   558:           c---------------------------------------------------------------------
   559:           
   560:                 use ft_data
   561:                 implicit none
   562:           
   563:                 integer is, d1, d2, d3, logd2
   564:                 double complex x(d1+1,d2,d3)
   565:                 double complex xout(d1+1,d2,d3)
   566:                 double complex y1(fftblockpad, d2), y2(fftblockpad, d2)
   567:                 integer i, j, k, ii, in
   568:           
   569:                 logd2 = ilog2(d2)
   570:           
   571:                 if (timers_enabled) call timer_start(T_ffty)
   572:           !$omp parallel do default(shared) private(i,j,k,ii,y1,y2,in)
   573:           !$omp&  shared(is,logd2,d2) collapse(2)
   574: P------>        do k = 1, d3
   575: |+----->          do in = 0, d1/fftblock - 1
   576: ||        !        do ii = 0, d1 - fftblock, fftblock
   577: ||                   ii = in*fftblock
   578: ||U---->             do j = 1, d2
   579: |||V--->                do i = 1, fftblock
   580: ||||                       y1(i,j) = x(i+ii,j,k)
   581: |||V---                 enddo
   582: ||U----              enddo
   583: ||        
   584: ||                   call cfftz (is, logd2, d2, y1, y2)
   585: ||         
   586: ||U---->             do j = 1, d2
   587: |||V--->                do i = 1, fftblock
   588: ||||                       xout(i+ii,j,k) = y1(i,j)
   589: |||V---                 enddo
   590: ||U----              enddo
   591: |+-----           enddo
   592: P------         enddo
   593:                 if (timers_enabled) call timer_stop(T_ffty)
   594:           
   595:                 return
   596:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:574)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : CFFTZ (ft.f:574)

    LOOP BEGIN: (ft.f:578)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:578)

      LOOP BEGIN: (ft.f:579)
        <Vectorized loop.>
        **  Copied loop. (ft.f:579)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:579)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:579)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:578)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:578)

      LOOP BEGIN: (ft.f:579)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:579)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:579)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:586)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:586)

      LOOP BEGIN: (ft.f:587)
        <Vectorized loop.>
        **  Copied loop. (ft.f:587)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:587)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:587)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:586)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:586)

      LOOP BEGIN: (ft.f:587)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:587)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:587)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 608 bytes
      Register spill area      : 136 bytes
      Parameter area           :  40 bytes
      Register save area       : 176 bytes
      User data area           : 248 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:574)
    *** Estimated execution cycle                       : 65
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 24
            Over basic blocks                           : 24
    *** The number of SCALAR REGISTER TRANSFER          : 45

    LOOP BEGIN: (ft.f:578)
      *** Estimated execution cycle                     : 18
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 20

      LOOP BEGIN: (ft.f:579)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 5
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:578)
      *** Estimated execution cycle                     : 36
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 32

      LOOP BEGIN: (ft.f:579)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:586)
      *** Estimated execution cycle                     : 17
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 11

      LOOP BEGIN: (ft.f:587)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 5
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:586)
      *** Estimated execution cycle                     : 34
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 24

      LOOP BEGIN: (ft.f:587)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS3
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   619: par(1801): Parallel routine generated.: CFFTS3$1
   621: opt(3008): Reference within a conditional branch moved outside loop.
   621: vec( 110): Vectorization obstructive procedure reference.: CFFTZ
   621: par(1803): Parallelized by "do".
   621: opt(2495): Loop nest collapsed by "collapse" clause.
   625: opt(1592): Outer loop unrolled inside inner loop.: K
   626: vec( 101): Vectorized loop.
   633: opt(1592): Outer loop unrolled inside inner loop.: K
   634: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS3
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

   602:                 subroutine cffts3(is, d1, d2, d3, x, xout, y1, y2)
   603:           
   604:           c---------------------------------------------------------------------
   605:           c---------------------------------------------------------------------
   606:           
   607:                 use ft_data
   608:                 implicit none
   609:           
   610:                 integer is, d1, d2, d3, logd3
   611:                 double complex x(d1+1,d2,d3)
   612:                 double complex xout(d1+1,d2,d3)
   613:                 double complex y1(fftblockpad, d3), y2(fftblockpad, d3)
   614:                 integer i, j, k, ii, in
   615:           
   616:                 logd3 = ilog2(d3)
   617:           
   618:                 if (timers_enabled) call timer_start(T_fftz)
   619:           !$omp parallel do default(shared) private(i,j,k,ii,y1,y2,in)
   620:           !$omp&  shared(is) collapse(2)
   621: P------>        do j = 1, d2
   622: |+----->          do in = 0, d1/fftblock - 1
   623: ||        !        do ii = 0, d1 - fftblock, fftblock
   624: ||                   ii = in*fftblock
   625: ||U---->             do k = 1, d3
   626: |||V--->                do i = 1, fftblock
   627: ||||                       y1(i,k) = x(i+ii,j,k)
   628: |||V---                 enddo
   629: ||U----              enddo
   630: ||        
   631: ||                   call cfftz (is, logd3, d3, y1, y2)
   632: ||        
   633: ||U---->             do k = 1, d3
   634: |||V--->                do i = 1, fftblock
   635: ||||                       xout(i+ii,j,k) = y1(i,k)
   636: |||V---                 enddo
   637: ||U----              enddo
   638: |+-----           enddo
   639: P------         enddo
   640:                 if (timers_enabled) call timer_stop(T_fftz)
   641:           
   642:                 return
   643:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS3
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:621)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : CFFTZ (ft.f:621)

    LOOP BEGIN: (ft.f:625)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:625)

      LOOP BEGIN: (ft.f:626)
        <Vectorized loop.>
        **  Copied loop. (ft.f:626)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:626)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:626)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:625)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:625)

      LOOP BEGIN: (ft.f:626)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:626)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:626)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:633)
      <Unvectorized loop.>
      **  Copied loop. (ft.f:633)

      LOOP BEGIN: (ft.f:634)
        <Vectorized loop.>
        **  Copied loop. (ft.f:634)
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:634)
        *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:634)
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:633)
      <Unvectorized loop.>
      **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:633)

      LOOP BEGIN: (ft.f:634)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:634)
        *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:634)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTS3
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 20 [v44-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 608 bytes
      Register spill area      : 136 bytes
      Parameter area           :  40 bytes
      Register save area       : 176 bytes
      User data area           : 248 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:621)
    *** Estimated execution cycle                       : 65
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 24
            Over basic blocks                           : 24
    *** The number of SCALAR REGISTER TRANSFER          : 45

    LOOP BEGIN: (ft.f:625)
      *** Estimated execution cycle                     : 18
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 20

      LOOP BEGIN: (ft.f:626)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 5
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:625)
      *** Estimated execution cycle                     : 36
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 1
              Over basic blocks                         : 1
      *** The number of SCALAR REGISTER TRANSFER        : 32

      LOOP BEGIN: (ft.f:626)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:633)
      *** Estimated execution cycle                     : 17
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 11

      LOOP BEGIN: (ft.f:634)
        *** Estimated execution cycle                   : 104
        *** The number of SCALAR REGISTER TRANSFER      : 5
      LOOP END
    LOOP END

    LOOP BEGIN: (ft.f:633)
      *** Estimated execution cycle                     : 34
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 2
              Over basic blocks                         : 2
      *** The number of SCALAR REGISTER TRANSFER        : 24

      LOOP BEGIN: (ft.f:634)
        *** Estimated execution cycle                   : 311
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   676: opt(1019): Feedback of scalar value from one loop pass to another.: LN
   678: vec( 101): Vectorized loop.
   680: opt(1019): Feedback of scalar value from one loop pass to another.: KU


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT_INIT
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

   649:                 subroutine fft_init (n)
   650:           
   651:           c---------------------------------------------------------------------
   652:           c---------------------------------------------------------------------
   653:           
   654:           c---------------------------------------------------------------------
   655:           c compute the roots-of-unity array that will be used for subsequent FFTs.
   656:           c---------------------------------------------------------------------
   657:           
   658:                 use ft_data
   659:                 implicit none
   660:           
   661:                 integer m,n,nu,ku,i,j,ln
   662:                 double precision t, ti
   663:           
   664:           
   665:           c---------------------------------------------------------------------
   666:           c   Initialize the U array with sines and cosines in a manner that permits
   667:           c   stride one access at each FFT iteration.
   668:           c---------------------------------------------------------------------
   669:                 nu = n
   670:                 m = ilog2(n)
   671:                 u(1) = m
   672:                 ku = 2
   673:                 ln = 1
   674:           
   675: +------>        do j = 1, m
   676: |                  t = pi / ln
   677: |          
   678: |V----->           do i = 0, ln - 1
   679: ||                    ti = i * t
   680: ||                    u(i+ku) = dcmplx (cos (ti), sin(ti))
   681: |V-----            enddo
   682: |          
   683: |                  ku = ku + ln
   684: |                  ln = 2 * ln
   685: +------         enddo
   686:            
   687:                 return
   688:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT_INIT
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:675)
    <Unvectorized loop.>

    LOOP BEGIN: (ft.f:678)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:678)
      *** The number of VLOAD, VSTORE. :  0,  2. (ft.f:678)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFT_INIT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 26 [s0-s2 s8-s12 s15-s16 s18 s23-s33 s60-s63]
      Vector registers         :  5 [v0-v1 v61-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 4288 bytes
      Register spill area      : 2048 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 2048 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:675)
    *** Estimated execution cycle                       : 43
    *** The number of SCALAR REGISTER TRANSFER          : 7

    LOOP BEGIN: (ft.f:678)
      *** Estimated execution cycle                     : 172
      *** The number of SCALAR REGISTER TRANSFER        : 3
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTZ
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   730: opt(3008): Reference within a conditional branch moved outside loop.
   730: vec( 110): Vectorization obstructive procedure reference.: FFTZ2
   730: vec( 103): Unvectorized loop.
   741: opt(1592): Outer loop unrolled inside inner loop.: J
   742: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTZ
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

   693:                 subroutine cfftz (is, m, n, x, y)
   694:           
   695:           c---------------------------------------------------------------------
   696:           c---------------------------------------------------------------------
   697:           
   698:           c---------------------------------------------------------------------
   699:           c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
   700:           c   to Swarztrauber.  X is both the input and the output array, while Y is a
   701:           c   scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to
   702:           c   perform FFTs, the array U must be initialized by calling CFFTZ with IS
   703:           c   set to 0 and M set to MX, where MX is the maximum value of M for any
   704:           c   subsequent call.
   705:           c---------------------------------------------------------------------
   706:           
   707:                 use ft_data
   708:                 implicit none
   709:           
   710:                 integer is,m,n,i,j,l,mx
   711:                 double complex x, y
   712:           
   713:                 dimension x(fftblockpad,n), y(fftblockpad,n)
   714:           
   715:           c---------------------------------------------------------------------
   716:           c   Check if input parameters are invalid.
   717:           c---------------------------------------------------------------------
   718:                 mx = u(1)
   719:                 if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 1 .or. m .gt. mx)
   720:                >  then
   721:                   write (*, 1)  is, m, mx
   722:            1      format ('CFFTZ: Either U has not been initialized, or else'/
   723:                >    'one of the input parameters is invalid', 3I5)
   724:                   stop
   725:                 endif
   726:           
   727:           c---------------------------------------------------------------------
   728:           c   Perform one variant of the Stockham FFT.
   729:           c---------------------------------------------------------------------
   730: +------>        do l = 1, m, 2
   731: |                 call fftz2 (is, l, m, n, fftblock, fftblockpad, u, x, y)
   732: |                 if (l .eq. m) goto 160
   733: |                 call fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, y, x)
   734: +------         enddo
   735:           
   736:                 goto 180
   737:           
   738:           c---------------------------------------------------------------------
   739:           c   Copy Y to X.
   740:           c---------------------------------------------------------------------
   741: U------>   160  do j = 1, n
   742: |V----->          do i = 1, fftblock
   743: ||                  x(i,j) = y(i,j)
   744: |V-----           enddo
   745: U------         enddo
   746:           
   747:            180  continue
   748:           
   749:                 return
   750:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTZ
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:730)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : FFTZ2 (ft.f:730)
  LOOP END

  LOOP BEGIN: (ft.f:741)
    <Unvectorized loop.>
    **  Copied loop. (ft.f:741)

    LOOP BEGIN: (ft.f:742)
      <Vectorized loop.>
      **  Copied loop. (ft.f:742)
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:742)
      *** The number of VLOAD, VSTORE. :  2,  2. (ft.f:742)
    LOOP END
  LOOP END

  LOOP BEGIN: (ft.f:741)
    <Unvectorized loop.>
    **  Outer loop unrolled loop. : THE UNROLL TIME IS 4. (ft.f:741)

    LOOP BEGIN: (ft.f:742)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ft.f:742)
      *** The number of VLOAD, VSTORE. :  8,  8. (ft.f:742)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CFFTZ
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 60 [s0-s12 s15-s16 s18-s34 s36-s63]
      Vector registers         : 10 [v54-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1360 bytes
      Register spill area      :   32 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1072 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:730)
    *** Estimated execution cycle                       : 22
    *** The number of SCALAR REGISTER TRANSFER          : 23
  LOOP END

  LOOP BEGIN: (ft.f:741)
    *** Estimated execution cycle                       : 10
    *** The number of SCALAR REGISTER TRANSFER          : 4

    LOOP BEGIN: (ft.f:742)
      *** Estimated execution cycle                     : 104
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (ft.f:741)
    *** Estimated execution cycle                       : 10
    *** The number of SCALAR REGISTER TRANSFER          : 3

    LOOP BEGIN: (ft.f:742)
      *** Estimated execution cycle                     : 311
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFTZ2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   781: vec( 102): Partially vectorized loop.
   796: vec( 101): Vectorized loop.
   800: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFTZ2
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

   755:                 subroutine fftz2 (is, l, m, n, ny, ny1, u, x, y)
   756:           
   757:           c---------------------------------------------------------------------
   758:           c---------------------------------------------------------------------
   759:           
   760:           c---------------------------------------------------------------------
   761:           c   Performs the L-th iteration of the second variant of the Stockham FFT.
   762:           c---------------------------------------------------------------------
   763:           
   764:                 implicit none
   765:           
   766:                 integer is,k,l,m,n,ny,ny1,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22
   767:                 double complex u,x,y,u1,x11,x21
   768:                 dimension u(n), x(ny1,n), y(ny1,n)
   769:           
   770:           
   771:           c---------------------------------------------------------------------
   772:           c   Set initial parameters.
   773:           c---------------------------------------------------------------------
   774:           
   775:                 n1 = n / 2
   776:                 lk = 2 ** (l - 1)
   777:                 li = 2 ** (m - l)
   778:                 lj = 2 * lk
   779:                 ku = li + 1
   780:           
   781: S------>        do i = 0, li - 1
   782: |                 i11 = i * lk + 1
   783: |                 i12 = i11 + n1
   784: |                 i21 = i * lj + 1
   785: |                 i22 = i21 + lk
   786: |                 if (is .ge. 1) then
   787: |                   u1 = u(ku+i)
   788: |                 else
   789: |                   u1 = dconjg (u(ku+i))
   790: |                 endif
   791: |         
   792: |         c---------------------------------------------------------------------
   793: |         c   This loop is vectorizable.
   794: |         c---------------------------------------------------------------------
   795: |+----->          do k = 0, lk - 1
   796: ||V---->            do j = 1, ny
   797: |||                   x11 = x(j,i11+k)
   798: |||                   x21 = x(j,i12+k)
   799: |||                   y(j,i21+k) = x11 + x21
   800: |||     F             y(j,i22+k) = u1 * (x11 - x21)
   801: ||V----             enddo
   802: |+-----           enddo
   803: S------         enddo
   804:           
   805:                 return
   806:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFTZ2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:781)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ft.f:781)
    *** The number of VLOAD, VSTORE. :  4,  4. (ft.f:781)

    LOOP BEGIN: (ft.f:795)
      <Unvectorized loop.>

      LOOP BEGIN: (ft.f:796)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (ft.f:796)
        *** The number of VLOAD, VSTORE. :  4,  4. (ft.f:796)
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: FFTZ2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 48 [s0-s12 s15-s16 s18-s29 s43-s63]
      Vector registers         : 17 [v47-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 4320 bytes
      Register spill area      : 2064 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           : 2064 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:781)
    *** Estimated execution cycle                       : 108
  LOOP END

  LOOP BEGIN: (ft.f:795)
    *** Estimated execution cycle                       : 7

    LOOP BEGIN: (ft.f:795)
      *** Estimated execution cycle                     : 12
      *** The number of SCALAR REGISTER TRANSFER        : 13

      LOOP BEGIN: (ft.f:796)
        *** Estimated execution cycle                   : 462
        *** The number of SCALAR REGISTER TRANSFER      : 2
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: ILOG2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   827: vec( 103): Unvectorized loop.
   827: vec( 108): Unvectorizable loop structure.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: ILOG2
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

   814:                 integer function ilog2(n)
   815:           
   816:           c---------------------------------------------------------------------
   817:           c---------------------------------------------------------------------
   818:           
   819:                 implicit none
   820:                 integer n, nn, lg
   821:                 if (n .eq. 1) then
   822:                    ilog2=0
   823:                    return
   824:                 endif
   825:                 lg = 1
   826:                 nn = 2
   827: +------>        do while (nn .lt. n)
   828: |                  nn = nn*2
   829: |                  lg = lg+1
   830: +------         end do
   831:                 ilog2 = lg
   832:                 return
   833:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: ILOG2
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:827)
    <Unvectorized loop.>
    *** Unvectorizable loop structure. (ft.f:827)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: ILOG2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         :  8 [s0 s8-s11 s18 s62-s63]
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

  LOOP BEGIN: (ft.f:827)
    *** Estimated execution cycle                       : 3
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CHECKSUM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   853: par(1801): Parallel routine generated.: CHECKSUM$1
   854: par(1803): Parallelized by "do".
   854: vec( 101): Vectorized loop.
   858: vec( 126): Idiom detected.: SUM
   859: par(1807): Critical-construct generated.
   859: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CHECKSUM
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

   839:                 subroutine checksum(i, u1, d1, d2, d3)
   840:           
   841:           c---------------------------------------------------------------------
   842:           c---------------------------------------------------------------------
   843:           
   844:                 use ft_data
   845:                 implicit none
   846:           
   847:                 integer i, d1, d2, d3
   848:                 double complex u1(d1+1,d2,d3)
   849:                 integer j, q,r,s
   850:                 double complex chk
   851:                 chk = (0.0,0.0)
   852:           
   853:           !$omp parallel do default(shared) private(i,q,r,s) reduction(+:chk)
   854: Y------>        do j=1,1024
   855: |                  q = mod(j, nx)+1
   856: |                  r = mod(3*j,ny)+1
   857: |                  s = mod(5*j,nz)+1
   858: |       G          chk=chk+u1(q,r,s)
   859: Y------         end do
   860:           
   861:                 chk = chk/ntotal_f
   862:            
   863:                 write (*, 30) i, chk
   864:            30   format (' T =',I5,5X,'Checksum =',1P2D22.12)
   865:                 sums(i) = chk
   866:                 return
   867:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CHECKSUM
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ft.f:854)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  2,  0. (ft.f:854)
    *** The number of VLOAD, VSTORE. :  0,  0. (ft.f:854)
    *** VGT generated (ft.f:858)
    *** Idiom detected. : SUM (ft.f:858)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: CHECKSUM
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 32 [s0-s2 s8-s12 s15-s16 s23-s27 s47-s63]
      Vector registers         : 37 [v27-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 20704 bytes
      Register spill area      : 10240 bytes
      Parameter area           :     8 bytes
      Register save area       :   176 bytes
      User data area           : 10272 bytes
      Others                   :     8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:854)
    *** Estimated execution cycle                       : 1306
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1111: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
  1111: vec( 103): Unvectorized loop.
  1113: opt(1025): Reference to this procedure inhibits optimization.: IEEE_ARITHMETIC::IEEE_IS_NAN


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: VERIFY
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

   873:                 subroutine verify (d1, d2, d3, nt, verified, class)
   874:           
   875:           c---------------------------------------------------------------------
   876:           c---------------------------------------------------------------------
   877:           
   878:                 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   879:           
   880:                 use ft_data
   881:           
   882:                 implicit none
   883:           
   884:                 integer d1, d2, d3, nt
   885:                 character class
   886:                 logical verified
   887:                 integer i
   888:                 double precision err, epsilon
   889:           
   890:           c---------------------------------------------------------------------
   891:           c   Reference checksums
   892:           c---------------------------------------------------------------------
   893:                 double complex csum_ref(25)
   894:           
   895:           
   896:                 class = 'U'
   897:           
   898:                 epsilon = 1.0d-12
   899:                 verified = .FALSE.
   900:           
   901:                 if (d1 .eq. 64 .and.
   902:                >    d2 .eq. 64 .and.
   903:                >    d3 .eq. 64 .and.
   904:                >    nt .eq. 6) then
   905:           c---------------------------------------------------------------------
   906:           c   Sample size reference checksums
   907:           c---------------------------------------------------------------------
   908:                    class = 'S'
   909:                    csum_ref(1) = dcmplx(5.546087004964D+02, 4.845363331978D+02)
   910:                    csum_ref(2) = dcmplx(5.546385409189D+02, 4.865304269511D+02)
   911:                    csum_ref(3) = dcmplx(5.546148406171D+02, 4.883910722336D+02)
   912:                    csum_ref(4) = dcmplx(5.545423607415D+02, 4.901273169046D+02)
   913:                    csum_ref(5) = dcmplx(5.544255039624D+02, 4.917475857993D+02)
   914:                    csum_ref(6) = dcmplx(5.542683411902D+02, 4.932597244941D+02)
   915:           
   916:                 else if (d1 .eq. 128 .and.
   917:                >    d2 .eq. 128 .and.
   918:                >    d3 .eq. 32 .and.
   919:                >    nt .eq. 6) then
   920:           c---------------------------------------------------------------------
   921:           c   Class W size reference checksums
   922:           c---------------------------------------------------------------------
   923:                    class = 'W'
   924:                    csum_ref(1) = dcmplx(5.673612178944D+02, 5.293246849175D+02)
   925:                    csum_ref(2) = dcmplx(5.631436885271D+02, 5.282149986629D+02)
   926:                    csum_ref(3) = dcmplx(5.594024089970D+02, 5.270996558037D+02)
   927:                    csum_ref(4) = dcmplx(5.560698047020D+02, 5.260027904925D+02)
   928:                    csum_ref(5) = dcmplx(5.530898991250D+02, 5.249400845633D+02)
   929:                    csum_ref(6) = dcmplx(5.504159734538D+02, 5.239212247086D+02)
   930:           
   931:                 else if (d1 .eq. 256 .and.
   932:                >    d2 .eq. 256 .and.
   933:                >    d3 .eq. 128 .and.
   934:                >    nt .eq. 6) then
   935:           c---------------------------------------------------------------------
   936:           c   Class A size reference checksums
   937:           c---------------------------------------------------------------------
   938:                    class = 'A'
   939:                    csum_ref(1) = dcmplx(5.046735008193D+02, 5.114047905510D+02)
   940:                    csum_ref(2) = dcmplx(5.059412319734D+02, 5.098809666433D+02)
   941:                    csum_ref(3) = dcmplx(5.069376896287D+02, 5.098144042213D+02)
   942:                    csum_ref(4) = dcmplx(5.077892868474D+02, 5.101336130759D+02)
   943:                    csum_ref(5) = dcmplx(5.085233095391D+02, 5.104914655194D+02)
   944:                    csum_ref(6) = dcmplx(5.091487099959D+02, 5.107917842803D+02)
   945:            
   946:                 else if (d1 .eq. 512 .and.
   947:                >    d2 .eq. 256 .and.
   948:                >    d3 .eq. 256 .and.
   949:                >    nt .eq. 20) then
   950:           c---------------------------------------------------------------------
   951:           c   Class B size reference checksums
   952:           c---------------------------------------------------------------------
   953:                    class = 'B'
   954:                    csum_ref(1)  = dcmplx(5.177643571579D+02, 5.077803458597D+02)
   955:                    csum_ref(2)  = dcmplx(5.154521291263D+02, 5.088249431599D+02)
   956:                    csum_ref(3)  = dcmplx(5.146409228649D+02, 5.096208912659D+02)
   957:                    csum_ref(4)  = dcmplx(5.142378756213D+02, 5.101023387619D+02)
   958:                    csum_ref(5)  = dcmplx(5.139626667737D+02, 5.103976610617D+02)
   959:                    csum_ref(6)  = dcmplx(5.137423460082D+02, 5.105948019802D+02)
   960:                    csum_ref(7)  = dcmplx(5.135547056878D+02, 5.107404165783D+02)
   961:                    csum_ref(8)  = dcmplx(5.133910925466D+02, 5.108576573661D+02)
   962:                    csum_ref(9)  = dcmplx(5.132470705390D+02, 5.109577278523D+02)
   963:                    csum_ref(10) = dcmplx(5.131197729984D+02, 5.110460304483D+02)
   964:                    csum_ref(11) = dcmplx(5.130070319283D+02, 5.111252433800D+02)
   965:                    csum_ref(12) = dcmplx(5.129070537032D+02, 5.111968077718D+02)
   966:                    csum_ref(13) = dcmplx(5.128182883502D+02, 5.112616233064D+02)
   967:                    csum_ref(14) = dcmplx(5.127393733383D+02, 5.113203605551D+02)
   968:                    csum_ref(15) = dcmplx(5.126691062020D+02, 5.113735928093D+02)
   969:                    csum_ref(16) = dcmplx(5.126064276004D+02, 5.114218460548D+02)
   970:                    csum_ref(17) = dcmplx(5.125504076570D+02, 5.114656139760D+02)
   971:                    csum_ref(18) = dcmplx(5.125002331720D+02, 5.115053595966D+02)
   972:                    csum_ref(19) = dcmplx(5.124551951846D+02, 5.115415130407D+02)
   973:                    csum_ref(20) = dcmplx(5.124146770029D+02, 5.115744692211D+02)
   974:           
   975:                 else if (d1 .eq. 512 .and.
   976:                >    d2 .eq. 512 .and.
   977:                >    d3 .eq. 512 .and.
   978:                >    nt .eq. 20) then
   979:           c---------------------------------------------------------------------
   980:           c   Class C size reference checksums
   981:           c---------------------------------------------------------------------
   982:                    class = 'C'
   983:                    csum_ref(1)  = dcmplx(5.195078707457D+02, 5.149019699238D+02)
   984:                    csum_ref(2)  = dcmplx(5.155422171134D+02, 5.127578201997D+02)
   985:                    csum_ref(3)  = dcmplx(5.144678022222D+02, 5.122251847514D+02)
   986:                    csum_ref(4)  = dcmplx(5.140150594328D+02, 5.121090289018D+02)
   987:                    csum_ref(5)  = dcmplx(5.137550426810D+02, 5.121143685824D+02)
   988:                    csum_ref(6)  = dcmplx(5.135811056728D+02, 5.121496764568D+02)
   989:                    csum_ref(7)  = dcmplx(5.134569343165D+02, 5.121870921893D+02)
   990:                    csum_ref(8)  = dcmplx(5.133651975661D+02, 5.122193250322D+02)
   991:                    csum_ref(9)  = dcmplx(5.132955192805D+02, 5.122454735794D+02)
   992:                    csum_ref(10) = dcmplx(5.132410471738D+02, 5.122663649603D+02)
   993:                    csum_ref(11) = dcmplx(5.131971141679D+02, 5.122830879827D+02)
   994:                    csum_ref(12) = dcmplx(5.131605205716D+02, 5.122965869718D+02)
   995:                    csum_ref(13) = dcmplx(5.131290734194D+02, 5.123075927445D+02)
   996:                    csum_ref(14) = dcmplx(5.131012720314D+02, 5.123166486553D+02)
   997:                    csum_ref(15) = dcmplx(5.130760908195D+02, 5.123241541685D+02)
   998:                    csum_ref(16) = dcmplx(5.130528295923D+02, 5.123304037599D+02)
   999:                    csum_ref(17) = dcmplx(5.130310107773D+02, 5.123356167976D+02)
  1000:                    csum_ref(18) = dcmplx(5.130103090133D+02, 5.123399592211D+02)
  1001:                    csum_ref(19) = dcmplx(5.129905029333D+02, 5.123435588985D+02)
  1002:                    csum_ref(20) = dcmplx(5.129714421109D+02, 5.123465164008D+02)
  1003:           
  1004:                 else if (d1 .eq. 2048 .and.
  1005:                >    d2 .eq. 1024 .and.
  1006:                >    d3 .eq. 1024 .and.
  1007:                >    nt .eq. 25) then
  1008:           c---------------------------------------------------------------------
  1009:           c   Class D size reference checksums
  1010:           c---------------------------------------------------------------------
  1011:                    class = 'D'
  1012:                    csum_ref(1)  = dcmplx(5.122230065252D+02, 5.118534037109D+02)
  1013:                    csum_ref(2)  = dcmplx(5.120463975765D+02, 5.117061181082D+02)
  1014:                    csum_ref(3)  = dcmplx(5.119865766760D+02, 5.117096364601D+02)
  1015:                    csum_ref(4)  = dcmplx(5.119518799488D+02, 5.117373863950D+02)
  1016:                    csum_ref(5)  = dcmplx(5.119269088223D+02, 5.117680347632D+02)
  1017:                    csum_ref(6)  = dcmplx(5.119082416858D+02, 5.117967875532D+02)
  1018:                    csum_ref(7)  = dcmplx(5.118943814638D+02, 5.118225281841D+02)
  1019:                    csum_ref(8)  = dcmplx(5.118842385057D+02, 5.118451629348D+02)
  1020:                    csum_ref(9)  = dcmplx(5.118769435632D+02, 5.118649119387D+02)
  1021:                    csum_ref(10) = dcmplx(5.118718203448D+02, 5.118820803844D+02)
  1022:                    csum_ref(11) = dcmplx(5.118683569061D+02, 5.118969781011D+02)
  1023:                    csum_ref(12) = dcmplx(5.118661708593D+02, 5.119098918835D+02)
  1024:                    csum_ref(13) = dcmplx(5.118649768950D+02, 5.119210777066D+02)
  1025:                    csum_ref(14) = dcmplx(5.118645605626D+02, 5.119307604484D+02)
  1026:                    csum_ref(15) = dcmplx(5.118647586618D+02, 5.119391362671D+02)
  1027:                    csum_ref(16) = dcmplx(5.118654451572D+02, 5.119463757241D+02)
  1028:                    csum_ref(17) = dcmplx(5.118665212451D+02, 5.119526269238D+02)
  1029:                    csum_ref(18) = dcmplx(5.118679083821D+02, 5.119580184108D+02)
  1030:                    csum_ref(19) = dcmplx(5.118695433664D+02, 5.119626617538D+02)
  1031:                    csum_ref(20) = dcmplx(5.118713748264D+02, 5.119666538138D+02)
  1032:                    csum_ref(21) = dcmplx(5.118733606701D+02, 5.119700787219D+02)
  1033:                    csum_ref(22) = dcmplx(5.118754661974D+02, 5.119730095953D+02)
  1034:                    csum_ref(23) = dcmplx(5.118776626738D+02, 5.119755100241D+02)
  1035:                    csum_ref(24) = dcmplx(5.118799262314D+02, 5.119776353561D+02)
  1036:                    csum_ref(25) = dcmplx(5.118822370068D+02, 5.119794338060D+02)
  1037:           
  1038:                 else if (d1 .eq. 4096 .and.
  1039:                >    d2 .eq. 2048 .and.
  1040:                >    d3 .eq. 2048 .and.
  1041:                >    nt .eq. 25) then
  1042:           c---------------------------------------------------------------------
  1043:           c   Class E size reference checksums
  1044:           c---------------------------------------------------------------------
  1045:                    class = 'E'
  1046:                    csum_ref(1)  = dcmplx(5.121601045346D+02, 5.117395998266D+02)
  1047:                    csum_ref(2)  = dcmplx(5.120905403678D+02, 5.118614716182D+02)
  1048:                    csum_ref(3)  = dcmplx(5.120623229306D+02, 5.119074203747D+02)
  1049:                    csum_ref(4)  = dcmplx(5.120438418997D+02, 5.119345900733D+02)
  1050:                    csum_ref(5)  = dcmplx(5.120311521872D+02, 5.119551325550D+02)
  1051:                    csum_ref(6)  = dcmplx(5.120226088809D+02, 5.119720179919D+02)
  1052:                    csum_ref(7)  = dcmplx(5.120169296534D+02, 5.119861371665D+02)
  1053:                    csum_ref(8)  = dcmplx(5.120131225172D+02, 5.119979364402D+02)
  1054:                    csum_ref(9)  = dcmplx(5.120104767108D+02, 5.120077674092D+02)
  1055:                    csum_ref(10) = dcmplx(5.120085127969D+02, 5.120159443121D+02)
  1056:                    csum_ref(11) = dcmplx(5.120069224127D+02, 5.120227453670D+02)
  1057:                    csum_ref(12) = dcmplx(5.120055158164D+02, 5.120284096041D+02)
  1058:                    csum_ref(13) = dcmplx(5.120041820159D+02, 5.120331373793D+02)
  1059:                    csum_ref(14) = dcmplx(5.120028605402D+02, 5.120370938679D+02)
  1060:                    csum_ref(15) = dcmplx(5.120015223011D+02, 5.120404138831D+02)
  1061:                    csum_ref(16) = dcmplx(5.120001570022D+02, 5.120432068837D+02)
  1062:                    csum_ref(17) = dcmplx(5.119987650555D+02, 5.120455615860D+02)
  1063:                    csum_ref(18) = dcmplx(5.119973525091D+02, 5.120475499442D+02)
  1064:                    csum_ref(19) = dcmplx(5.119959279472D+02, 5.120492304629D+02)
  1065:                    csum_ref(20) = dcmplx(5.119945006558D+02, 5.120506508902D+02)
  1066:                    csum_ref(21) = dcmplx(5.119930795911D+02, 5.120518503782D+02)
  1067:                    csum_ref(22) = dcmplx(5.119916728462D+02, 5.120528612016D+02)
  1068:                    csum_ref(23) = dcmplx(5.119902874185D+02, 5.120537101195D+02)
  1069:                    csum_ref(24) = dcmplx(5.119889291565D+02, 5.120544194514D+02)
  1070:                    csum_ref(25) = dcmplx(5.119876028049D+02, 5.120550079284D+02)
  1071:           
  1072:                 else if (d1 .eq. 8192 .and.
  1073:                >    d2 .eq. 4096 .and.
  1074:                >    d3 .eq. 4096 .and.
  1075:                >    nt .eq. 25) then
  1076:           c---------------------------------------------------------------------
  1077:           c   Class F size reference checksums
  1078:           c---------------------------------------------------------------------
  1079:                    class = 'F'
  1080:                    csum_ref( 1) = dcmplx(5.119892866928D+02, 5.121457822747D+02)
  1081:                    csum_ref( 2) = dcmplx(5.119560157487D+02, 5.121009044434D+02)
  1082:                    csum_ref( 3) = dcmplx(5.119437960123D+02, 5.120761074285D+02)
  1083:                    csum_ref( 4) = dcmplx(5.119395628845D+02, 5.120614320496D+02)
  1084:                    csum_ref( 5) = dcmplx(5.119390371879D+02, 5.120514085624D+02)
  1085:                    csum_ref( 6) = dcmplx(5.119405091840D+02, 5.120438117102D+02)
  1086:                    csum_ref( 7) = dcmplx(5.119430444528D+02, 5.120376348915D+02)
  1087:                    csum_ref( 8) = dcmplx(5.119460702242D+02, 5.120323831062D+02)
  1088:                    csum_ref( 9) = dcmplx(5.119492377036D+02, 5.120277980818D+02)
  1089:                    csum_ref(10) = dcmplx(5.119523446268D+02, 5.120237368268D+02)
  1090:                    csum_ref(11) = dcmplx(5.119552825361D+02, 5.120201137845D+02)
  1091:                    csum_ref(12) = dcmplx(5.119580008777D+02, 5.120168723492D+02)
  1092:                    csum_ref(13) = dcmplx(5.119604834177D+02, 5.120139707209D+02)
  1093:                    csum_ref(14) = dcmplx(5.119627332821D+02, 5.120113749334D+02)
  1094:                    csum_ref(15) = dcmplx(5.119647637538D+02, 5.120090554887D+02)
  1095:                    csum_ref(16) = dcmplx(5.119665927740D+02, 5.120069857863D+02)
  1096:                    csum_ref(17) = dcmplx(5.119682397643D+02, 5.120051414260D+02)
  1097:                    csum_ref(18) = dcmplx(5.119697238718D+02, 5.120034999132D+02)
  1098:                    csum_ref(19) = dcmplx(5.119710630664D+02, 5.120020405355D+02)
  1099:                    csum_ref(20) = dcmplx(5.119722737384D+02, 5.120007442976D+02)
  1100:                    csum_ref(21) = dcmplx(5.119733705802D+02, 5.119995938652D+02)
  1101:                    csum_ref(22) = dcmplx(5.119743666226D+02, 5.119985735001D+02)
  1102:                    csum_ref(23) = dcmplx(5.119752733481D+02, 5.119976689792D+02)
  1103:                    csum_ref(24) = dcmplx(5.119761008382D+02, 5.119968675026D+02)
  1104:                    csum_ref(25) = dcmplx(5.119768579280D+02, 5.119961575929D+02)
  1105:           
  1106:                 endif
  1107:           
  1108:           
  1109:                 if (class .ne. 'U') then
  1110:           
  1111: +------>           do i = 1, nt
  1112: |                     err = abs( (sums(i) - csum_ref(i)) / csum_ref(i) )
  1113: |                     if (ieee_is_nan(err) .or. (err .gt. epsilon)) goto 100
  1114: +------            end do
  1115:                    verified = .TRUE.
  1116:            100     continue
  1117:           
  1118:                 endif
  1119:           
  1120:            
  1121:                 if (class .ne. 'U') then
  1122:                    if (verified) then
  1123:                       write(*,2000)
  1124:            2000       format(' Result verification successful')
  1125:                    else
  1126:                       write(*,2001)
  1127:            2001       format(' Result verification failed')
  1128:                    endif
  1129:                 endif
  1130:                 print *, 'class = ', class
  1131:           
  1132:                 return
  1133:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: VERIFY
VECTORIZATION LIST

  LOOP BEGIN: (ft.f:1111)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : IEEE_ARITHMETIC::IEEE_IS_NAN (ft.f:1111)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:58 2021
FILE NAME: ft.f

PROCEDURE NAME: VERIFY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 31 [s0-s12 s15-s16 s18-s19 s23-s28 s33 s57-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1696 bytes
      Register spill area      :    0 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1440 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ft.f:1111)
    *** Estimated execution cycle                       : 35
    *** The number of SCALAR REGISTER TRANSFER          : 5
  LOOP END


