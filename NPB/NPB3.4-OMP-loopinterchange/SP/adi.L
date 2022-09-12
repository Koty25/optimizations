NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:42 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:42 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  adi
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:                  call compute_rhs
    11:           
    12:                  call txinvr
    13:           
    14:                  call x_solve
    15:           
    16:                  call y_solve
    17:           
    18:                  call z_solve
    19:           
    20:                  call add
    21:           
    22:                  return
    23:                  end


