NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:37 2020
FILE NAME: ft_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:37 2020
FILE NAME: ft_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

   157:                 subroutine alloc_space
   158:           
   159:           c---------------------------------------------------------------------
   160:           c---------------------------------------------------------------------
   161:           
   162:           c---------------------------------------------------------------------
   163:           c allocate space dynamically for data arrays
   164:           c---------------------------------------------------------------------
   165:           
   166:                 use ft_data, only : ntotalp
   167:                 use ft_fields
   168:           
   169:                 implicit none
   170:           
   171:                 integer ios
   172:           
   173:           
   174:                 allocate (
   175:                >          u0(ntotalp), pad1(3),
   176:                >          u1(ntotalp), pad2(3),
   177:           c     >          u2(ntotalp),
   178:                >          twiddle(ntotalp),
   179:                >          stat = ios)
   180:           
   181:                 if (ios .ne. 0) then
   182:                    write(*,*) 'Error encountered in allocating space'
   183:                    stop
   184:                 endif
   185:           
   186:                 return
   187:                 end


