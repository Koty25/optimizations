NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: mg_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: mg_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

    97:                 subroutine alloc_space
    98:           
    99:           c---------------------------------------------------------------------
   100:           c---------------------------------------------------------------------
   101:           
   102:           c---------------------------------------------------------------------
   103:           c allocate space dynamically for field arrays
   104:           c---------------------------------------------------------------------
   105:           
   106:                 use mg_data, only : nr, nv
   107:                 use mg_fields
   108:           
   109:                 implicit none
   110:           
   111:                 integer ios
   112:           
   113:                 allocate( u(nr), v(nv), r(nr),
   114:                >          stat = ios )
   115:           
   116:                 if (ios .ne. 0) then
   117:                    write(*,*) 'Error encountered in allocating space'
   118:                    stop
   119:                 endif
   120:           
   121:                 return
   122:                 end


