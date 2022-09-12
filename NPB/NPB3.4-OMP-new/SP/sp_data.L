NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: sp_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: sp_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

    93:                 subroutine alloc_space
    94:           
    95:           c---------------------------------------------------------------------
    96:           c---------------------------------------------------------------------
    97:           
    98:           c---------------------------------------------------------------------
    99:           c allocate space dynamically for data arrays
   100:           c---------------------------------------------------------------------
   101:           
   102:                 use sp_data
   103:                 implicit none
   104:           
   105:                 integer ios
   106:           
   107:                 integer IMAXP, JMAXP
   108:                 parameter (IMAXP=IMAX/2*2,JMAXP=JMAX/2*2)
   109:           
   110:           c
   111:           c   To improve cache performance, first two dimensions padded by 1
   112:           c   for even number sizes only
   113:           c
   114:                 allocate (
   115:                >   u       (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   116:                >   us      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   117:                >   vs      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   118:                >   ws      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   119:                >   qs      (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   120:                >   rho_i   (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   121:                >   speed   (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   122:                >   square  (   0:IMAXP, 0:JMAXP, 0:KMAX-1),
   123:                >   rhs     (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   124:                >   forcing (5, 0:IMAXP, 0:JMAXP, 0:KMAX-1),
   125:                >         stat = ios)
   126:           
   127:                 if (ios .ne. 0) then
   128:                    write(*,*) 'Error encountered in allocating space'
   129:                    stop
   130:                 endif
   131:           
   132:                 return
   133:                 end


