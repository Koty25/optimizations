NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME : timers.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
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

    17:                 subroutine timer_clear(n)
    18:           
    19:           c---------------------------------------------------------------------
    20:           c---------------------------------------------------------------------
    21:           
    22:                 use timers
    23:                 implicit none
    24:           
    25:                 integer n
    26:           
    27:                 elapsed(n) = 0.0
    28:                 return
    29:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         :  9 [s0 s8-s11 s23 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    45: inl(1222): Inlined: ELAPSED_TIME


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
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

    35:                 subroutine timer_start(n)
    36:           
    37:           c---------------------------------------------------------------------
    38:           c---------------------------------------------------------------------
    39:           
    40:                 use timers
    41:                 implicit none
    42:           
    43:                 integer n
    44:           
    45:         I       start(n) = elapsed_time()
    46:           
    47:                 return
    48:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
INLINE LIST

  ROOT: TIMER_START (timers.f:35)
  -> INLINE: ELAPSED_TIME (timers.f:45)
   -> NOINLINE: WTIME (timers.f:118)
      *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0 s8-s12 s15-s16 s23-s24 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: inl(1222): Inlined: ELAPSED_TIME


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
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

    54:                 subroutine timer_stop(n)
    55:           
    56:           c---------------------------------------------------------------------
    57:           c---------------------------------------------------------------------
    58:           
    59:                 use timers
    60:                 implicit none
    61:           
    62:                 integer n
    63:           
    64:                 double precision t, now
    65:           
    66:         I       now = elapsed_time()
    67:                 t = now - start(n)
    68:                 elapsed(n) = elapsed(n) + t
    69:           
    70:                 return
    71:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
INLINE LIST

  ROOT: TIMER_STOP (timers.f:54)
  -> INLINE: ELAPSED_TIME (timers.f:66)
   -> NOINLINE: WTIME (timers.f:118)
      *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0 s8-s12 s15-s16 s23-s24 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
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

    77:                 double precision function timer_read(n)
    78:           
    79:           c---------------------------------------------------------------------
    80:           c---------------------------------------------------------------------
    81:           
    82:                 use timers
    83:                 implicit none
    84:           
    85:                 integer n
    86:            
    87:                 timer_read = elapsed(n)
    88:           
    89:                 return
    90:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         :  9 [s0 s8-s11 s23 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   113: opt(1317): Logical clause simplified.
   118: inl(1212): Source for routine not found.: WTIME


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
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

    96:                 double precision function elapsed_time()
    97:           
    98:           c---------------------------------------------------------------------
    99:           c---------------------------------------------------------------------
   100:           
   101:                 implicit none
   102:           c$    external         omp_get_wtime
   103:           c$    double precision omp_get_wtime
   104:           
   105:                 double precision t
   106:                 logical          mp
   107:           
   108:           c ... Use the OpenMP timer if we can (via C$ conditional compilation)
   109:                 mp = .false.
   110:           c$    mp = .true.
   111:           c$    t = omp_get_wtime()
   112:           
   113:                 if (.not.mp) then
   114:           c This function must measure wall clock time, not CPU time.
   115:           c Since there is no portable timer in Fortran (77)
   116:           c we call a routine compiled in C (though the C source may have
   117:           c to be tweaked).
   118:                    call wtime(t)
   119:           c The following is not ok for "official" results because it reports
   120:           c CPU time not wall clock time. It may be useful for developing/testing
   121:           c on timeshared Crays, though.
   122:           c        call second(t)
   123:                 endif
   124:           
   125:                 elapsed_time = t
   126:           
   127:                 return
   128:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
INLINE LIST

  ROOT: ELAPSED_TIME (timers.f:96)
  -> NOINLINE: WTIME (timers.f:118)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 11 [s0 s8-s12 s15-s16 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
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

   133:                 subroutine check_timer_flag( timeron )
   134:           
   135:           c---------------------------------------------------------------------
   136:           c---------------------------------------------------------------------
   137:           
   138:                 implicit none
   139:                 logical timeron
   140:           
   141:                 integer nc, ios
   142:                 character(len=20) val
   143:           
   144:                 timeron = .false.
   145:           
   146:           c ... Check environment variable "NPB_TIMER_FLAG"
   147:                 call get_environment_variable('NPB_TIMER_FLAG', val, nc, ios)
   148:                 if (ios .eq. 0) then
   149:                    if (nc .le. 0) then
   150:                       timeron = .true.
   151:                    else if (val(1:1) .ge. '1' .and. val(1:1) .le. '9') then
   152:                       timeron = .true.
   153:                    else if (val .eq. 'on' .or. val .eq. 'ON' .or.
   154:                &            val .eq. 'yes' .or. val .eq. 'YES' .or.
   155:                &            val .eq. 'true' .or. val .eq. 'TRUE') then
   156:                       timeron = .true.
   157:                    endif
   158:           
   159:                 else
   160:           
   161:           c ... Check if the "timer.flag" file exists
   162:                    open (unit=2, file='timer.flag', status='old', iostat=ios)
   163:                    if (ios .eq. 0) then
   164:                       close(2)
   165:                       timeron = .true.
   166:                    endif
   167:           
   168:                 endif
   169:           
   170:                 return
   171:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 26 [s0-s12 s15-s16 s18-s19 s23-s28 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 304 bytes
      Register spill area      :   0 bytes
      Parameter area           :  96 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


