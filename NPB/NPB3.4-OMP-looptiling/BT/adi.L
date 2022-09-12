NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:16 2020
FILE NAME: adi.f

PROCEDURE NAME: ADI
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

     9: inl(1212): Source for routine not found.: COMPUTE_RHS
    11: inl(1212): Source for routine not found.: X_SOLVE
    13: inl(1212): Source for routine not found.: Y_SOLVE
    15: inl(1212): Source for routine not found.: Z_SOLVE
    17: inl(1212): Source for routine not found.: ADD


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:16 2020
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


