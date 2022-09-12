NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: cg_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: cg_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

    82:                 subroutine alloc_space
    83:           
    84:           c---------------------------------------------------------------------
    85:           c---------------------------------------------------------------------
    86:           
    87:           c---------------------------------------------------------------------
    88:           c allocate space dynamically for data arrays
    89:           c---------------------------------------------------------------------
    90:           
    91:                 use cg_data
    92:                 implicit none
    93:           
    94:                 integer ios
    95:           
    96:           
    97:                 allocate (
    98:                >          colidx(nz), rowstr(na+1),
    99:                >          iv(nz+na),  arow(na), acol(naz),
   100:                >          v(nz), aelt(naz), a(nz),
   101:                >          x(na+2),
   102:                >          z(na+2),
   103:                >          p(na+2),
   104:                >          q(na+2),
   105:                >          r(na+2),
   106:                >          stat = ios)
   107:           
   108:                 if (ios .ne. 0) then
   109:                    write(*,*) 'Error encountered in allocating space'
   110:                    stop
   111:                 endif
   112:           
   113:                 return
   114:                 end


