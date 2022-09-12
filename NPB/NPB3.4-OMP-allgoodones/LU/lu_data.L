NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: lu_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: lu_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

   132:                 subroutine alloc_space
   133:           
   134:           c---------------------------------------------------------------------
   135:           c---------------------------------------------------------------------
   136:           
   137:           c---------------------------------------------------------------------
   138:           c allocate space dynamically for data arrays
   139:           c---------------------------------------------------------------------
   140:           
   141:                 use lu_data
   142:                 implicit none
   143:           
   144:                 integer ios
   145:           
   146:           c---------------------------------------------------------------------
   147:           c   to improve cache performance, second two dimensions padded by 1
   148:           c   for even number sizes only.
   149:           c   Note: corresponding array (called "v") in routines blts, buts,
   150:           c   and l2norm are similarly padded
   151:           c---------------------------------------------------------------------
   152:           
   153:                 allocate (
   154:                >          u   (5,isiz1/2*2+1,isiz2/2*2+1,isiz3),
   155:                >          rsd (5,isiz1/2*2+1,isiz2/2*2+1,isiz3),
   156:                >          frct(5,isiz1/2*2+1,isiz2/2*2+1,isiz3),
   157:                >          qs    (isiz1/2*2+1,isiz2/2*2+1,isiz3),
   158:                >          rho_i (isiz1/2*2+1,isiz2/2*2+1,isiz3),
   159:                >          phi1  (0:isiz2+1,0:isiz3+1),
   160:                >          phi2  (0:isiz2+1,0:isiz3+1),
   161:                >          stat = ios)
   162:           
   163:                 if (ios .ne. 0) then
   164:                    write(*,*) 'Error encountered in allocating space'
   165:                    stop
   166:                 endif
   167:           
   168:                 return
   169:                 end


