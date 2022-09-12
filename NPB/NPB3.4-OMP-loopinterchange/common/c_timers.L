NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: elapsed_time
FORMAT LIST

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


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_clear
FORMAT LIST

 LINE   LOOP      STATEMENT

    38:           void timer_clear( int n )
    39:           {
    40:               elapsed[n] = 0.0;
    41:           }


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_start
FORMAT LIST

 LINE   LOOP      STATEMENT

    47:           void timer_start( int n )
    48:           {
    49:               start[n] = elapsed_time();
    50:           }


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_stop
FORMAT LIST

 LINE   LOOP      STATEMENT

    56:           void timer_stop( int n )
    57:           {
    58:               double t, now;
    59:           
    60:               now = elapsed_time();
    61:               t = now - start[n];
    62:               elapsed[n] += t;
    63:           
    64:           }


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: timer_read
FORMAT LIST

 LINE   LOOP      STATEMENT

    70:           double timer_read( int n )
    71:           {
    72:               return( elapsed[n] );
    73:           }


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Sat Jul 18 20:59:37 2020
FILE NAME: c_timers.c

FUNCTION NAME: check_timer_flag
FORMAT LIST

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


