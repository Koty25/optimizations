NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:15 2020
FILE NAME: exact_solution.f

PROCEDURE NAME: EXACT_SOLUTION
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    19: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:15 2020
FILE NAME: exact_solution.f

PROCEDURE NAME: EXACT_SOLUTION
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine exact_solution(xi,eta,zeta,dtemp)
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     this function returns the exact solution at point xi, eta, zeta
    11:           c---------------------------------------------------------------------
    12:           
    13:                 use bt_data
    14:                 implicit none
    15:           
    16:                 double precision  xi, eta, zeta, dtemp(5)
    17:                 integer m
    18:           
    19: *------>        do m = 1, 5
    20: |                  dtemp(m) =  ce(m,1) +
    21: |              >     xi*(ce(m,2) + xi*(ce(m,5) + xi*(ce(m,8) + xi*ce(m,11)))) +
    22: |              >     eta*(ce(m,3) + eta*(ce(m,6) + eta*(ce(m,9) + eta*ce(m,12))))+
    23: |              >     zeta*(ce(m,4) + zeta*(ce(m,7) + zeta*(ce(m,10) +
    24: |              >     zeta*ce(m,13))))
    25: *------         enddo
    26:           
    27:                 return
    28:                 end


