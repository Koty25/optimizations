NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:34 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:34 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine  adi
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:                 call compute_rhs
    10:           
    11:                 call x_solve
    12:           
    13:                 call y_solve
    14:           
    15:                 call z_solve
    16:           
    17:                 call add
    18:           
    19:                 return
    20:                 end


