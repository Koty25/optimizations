NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME : c_timers.c

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fcse-after-vectorization                : disable
        -ffast-math                              : enable
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
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : enable
        -fmove-nested-loop-invariants-outer      : enable
        -fnaked-ivdep                            : disable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freplace-loop-equation                  : enable
        -fstrict-aliasing                        : enable
        -fthis-pointer-alias                     : disable
        -fthis-pointer-noalias                   : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : disable
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
        -finline                                 : enable
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
        -g                                       : disable
        -minit-stack                             : not specified
        -traceback                               : disable

      Language Options : 
        -fallow-keyword-macros                   : enable
        -fgnu89-inline                           : disable
        -frestrict                               : enable
        -fsigned-char                            : enable
        -std                                     : gnu11

      Message Options : 
        -Wall                                    : disable
        -Wcomment                                : disable
        -Werror                                  : disable
        -Wno-div-by-zero                         : disable
        -Wunknown-pragma                         : disable
        -Wunused                                 : disable
        -Wunused-but-set-parameter               : disable
        -Wunused-but-set-variable                : disable
        -Wunused-parameter                       : disable
        -Wunused-value                           : disable
        -Wunused-variable                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -fdiag-system-header                     : disable
        -pedantic                                : disable
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

      Miscellaneous Options : 
        -v                                       : disable


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
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

    16:           double elapsed_time( void )
    17:           {
    18:               double t;
    19:           
    20:           #if defined(_OPENMP) && (_OPENMP > 200010)
    21:           /*  Use the OpenMP timer if we can */
    22:               t = omp_get_wtime();
    23:           #else
    24:               wtime( &t );
    25:           #endif
    26:               return( t );
    27:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
INLINE LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
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
    Total size                 : 176 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
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

    38:           void timer_clear( int n )
    39:           {
    40:               elapsed[n] = 0.0;
    41:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
INLINE LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    49: inl(1222): Inlined: elapsed_time


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
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

    47:           void timer_start( int n )
    48:           {
    49:         I     start[n] = elapsed_time();
    50:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
INLINE LIST

  ROOT: timer_start (c_timers.c:47)
  -> INLINE: elapsed_time (c_timers.c:49)


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    60: inl(1222): Inlined: elapsed_time


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
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

    56:           void timer_stop( int n )
    57:           {
    58:               double t, now;
    59:           
    60:         I     now = elapsed_time();
    61:               t = now - start[n];
    62:               elapsed[n] += t;
    63:           
    64:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
INLINE LIST

  ROOT: timer_stop (c_timers.c:56)
  -> INLINE: elapsed_time (c_timers.c:60)


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
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

    70:           double timer_read( int n )
    71:           {
    72:               return( elapsed[n] );
    73:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
INLINE LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
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

    79:           int check_timer_flag( void )
    80:           {
    81:               int timer_on = 0;
    82:               char *ev = getenv("NPB_TIMER_FLAG");
    83:           
    84:               if (ev) {
    85:                   if (*ev == '\0')
    86:                       timer_on = 1;
    87:                   else if (*ev >= '1' && *ev <= '9')
    88:                       timer_on = 1;
    89:                   else if (strcmp(ev, "on") == 0 || strcmp(ev, "ON") == 0 ||
    90:                            strcmp(ev, "yes") == 0 || strcmp(ev, "YES") == 0 ||
    91:                            strcmp(ev, "true") == 0 || strcmp(ev, "TRUE") == 0)
    92:                       timer_on = 1;
    93:               }
    94:               else {
    95:                   FILE *fp = fopen("timer.flag", "r");
    96:                   if (fp != NULL) {
    97:                       fclose(fp);
    98:                       timer_on = 1;
    99:                   }
   100:               }
   101:           
   102:               return timer_on;
   103:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
INLINE LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
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
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


