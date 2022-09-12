NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:25 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    10: inl(1212): Source for routine not found.: COMPUTE_RHS
    12: inl(1212): Source for routine not found.: TXINVR
    14: inl(1212): Source for routine not found.: X_SOLVE
    16: inl(1212): Source for routine not found.: Y_SOLVE
    18: inl(1212): Source for routine not found.: Z_SOLVE
    20: inl(1212): Source for routine not found.: ADD


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:25 2020
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


