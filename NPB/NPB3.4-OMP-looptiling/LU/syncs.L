NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    28:                 subroutine sync_init( jdim )
    29:           
    30:           c---------------------------------------------------------------------
    31:           c   Initialize sync-related variables
    32:           c---------------------------------------------------------------------
    33:           
    34:                 use syncs
    35:                 implicit none
    36:           
    37:                 integer jdim
    38:           
    39:           !$    integer, external :: omp_get_num_threads, omp_get_thread_num
    40:           
    41:                 mthreadnum = 0
    42:           !$    mthreadnum = omp_get_num_threads() - 1
    43:                 if (mthreadnum .gt. jdim) mthreadnum = jdim
    44:                 iam = 0
    45:           !$    iam = omp_get_thread_num()
    46:                 if (iam .le. mthreadnum) isync(1,iam) = 0
    47:           
    48:                 return
    49:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    79: vec( 103): Unvectorized loop.
    79: vec( 109): Vectorization obstructive statement.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
FORMAT LIST

 LINE   LOOP      STATEMENT

    57:                 subroutine sync_left( ldmx, ldmy, ldmz, v )
    58:           
    59:           c---------------------------------------------------------------------
    60:           c   Thread synchronization for pipeline operation
    61:           c---------------------------------------------------------------------
    62:           
    63:                 use syncs
    64:                 implicit none
    65:           
    66:                 integer ldmx, ldmy, ldmz
    67:                 double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)
    68:           
    69:           c---------------------------------------------------------------------
    70:           c---------------------------------------------------------------------
    71:           
    72:                 integer neigh, iv
    73:           
    74:           
    75:                 if (iam .gt. 0 .and. iam .le. mthreadnum) then
    76:                    neigh = iam - 1
    77:           !$omp atomic read
    78:                    iv = isync(1,neigh)
    79: +------>           do while (iv .eq. 0)
    80: |         !$omp atomic read
    81: |                     iv = isync(1,neigh)
    82: +------            end do
    83:           !$omp atomic write
    84:                    isync(1,neigh) = 0
    85:                 endif
    86:           !$omp flush(isync,v)
    87:           
    88:           
    89:                 return
    90:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   117: vec( 103): Unvectorized loop.
   117: vec( 109): Vectorization obstructive statement.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
FORMAT LIST

 LINE   LOOP      STATEMENT

    95:                 subroutine sync_right( ldmx, ldmy, ldmz, v )
    96:           
    97:           c---------------------------------------------------------------------
    98:           c   Thread synchronization for pipeline operation
    99:           c---------------------------------------------------------------------
   100:           
   101:                 use syncs
   102:                 implicit none
   103:           
   104:                 integer ldmx, ldmy, ldmz
   105:                 double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)
   106:           
   107:           c---------------------------------------------------------------------
   108:           c---------------------------------------------------------------------
   109:           
   110:                 integer iv
   111:           
   112:           
   113:           !$omp flush(isync,v)
   114:                 if (iam .lt. mthreadnum) then
   115:           !$omp atomic read
   116:                    iv = isync(1,iam)
   117: +------>           do while (iv .eq. 1)
   118: |         !$omp atomic read
   119: |                     iv = isync(1,iam)
   120: +------            end do
   121:           !$omp atomic write
   122:                    isync(1,iam) = 1
   123:                 endif
   124:           
   125:           
   126:                 return
   127:                 end


