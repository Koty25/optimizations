NEC C/C++ Compiler (3.0.4) for Vector Engine     Thu Feb 11 10:50:52 2021
FILE NAME: ../common/wtime.c

FUNCTION NAME: wtime_
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Thu Feb 11 10:50:52 2021
FILE NAME: ../common/wtime.c

FUNCTION NAME: wtime_
FORMAT LIST

 LINE   LOOP      STATEMENT

     7:           void wtime(double *t)
     8:           {
     9:              /* a generic timer */
    10:              static int sec = -1;
    11:              struct timeval tv;
    12:              gettimeofday(&tv, (void *)0);
    13:              if (sec < 0) sec = tv.tv_sec;
    14:              *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
    15:           }


