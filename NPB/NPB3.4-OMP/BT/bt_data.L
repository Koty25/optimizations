NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:47 2021
FILE NAME: bt_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:47 2021
FILE NAME: bt_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

    92:                 subroutine alloc_space
    93:           
    94:           c---------------------------------------------------------------------
    95:           c---------------------------------------------------------------------
    96:           
    97:           c---------------------------------------------------------------------
    98:           c allocate space dynamically for data arrays
    99:           c---------------------------------------------------------------------
   100:           
   101:                 use bt_data
   102:                 implicit none
   103:           
   104:                 integer ios
   105:           
   106:                 integer IMAXP, JMAXP
   107:                 parameter (IMAXP=IMAX/2*2,JMAXP=JMAX/2*2)
   108:           
   109:           c
   110:           c   to improve cache performance, grid dimensions padded by 1
   111:           c   for even number sizes only.
   112:           c
   113:                 allocate (
   114:                >   us      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   115:                >   vs      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   116:                >   ws      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   117:                >   qs      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   118:                >   rho_i   (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   119:                >   square  (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   120:                >   forcing (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   121:                >   u       (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   122:                >   rhs     (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   123:                >         stat = ios)
   124:           
   125:                 if (ios .ne. 0) then
   126:                    write(*,*) 'Error encountered in allocating space'
   127:                    stop
   128:                 endif
   129:           
   130:                 return
   131:                 end


