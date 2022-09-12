NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:42 2020
FILE NAME: exact_solution.f

PROCEDURE NAME: EXACT_SOLUTION
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    21: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:42 2020
FILE NAME: exact_solution.f

PROCEDURE NAME: EXACT_SOLUTION
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine exact_solution(xi,eta,zeta,dtemp)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function returns the exact solution at point xi, eta, zeta
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  double precision  xi, eta, zeta
    18:                  double precision  dtemp(5)
    19:                  integer m
    20:           
    21: *------>         do  m = 1, 5
    22: |                   dtemp(m) =  ce(m,1) +
    23: |              >    xi*(ce(m,2) + xi*(ce(m,5) + xi*(ce(m,8) + xi*ce(m,11)))) +
    24: |              >    eta*(ce(m,3) + eta*(ce(m,6) + eta*(ce(m,9) + eta*ce(m,12))))+
    25: |              >    zeta*(ce(m,4) + zeta*(ce(m,7) + zeta*(ce(m,10) +
    26: |              >    zeta*ce(m,13))))
    27: *------          end do
    28:           
    29:                  return
    30:                  end


