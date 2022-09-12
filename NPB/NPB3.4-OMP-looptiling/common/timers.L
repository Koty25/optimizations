NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_CLEAR
FORMAT LIST

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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    45: inl(1222): Inlined: ELAPSED_TIME


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_START
FORMAT LIST

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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: inl(1222): Inlined: ELAPSED_TIME


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_STOP
FORMAT LIST

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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: TIMER_READ
FORMAT LIST

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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   113: opt(1317): Logical clause simplified.
   118: inl(1212): Source for routine not found.: WTIME


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: ELAPSED_TIME
FORMAT LIST

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


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: timers.f

PROCEDURE NAME: CHECK_TIMER_FLAG
FORMAT LIST

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


