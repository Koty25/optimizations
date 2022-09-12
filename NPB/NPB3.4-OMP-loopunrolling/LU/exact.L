NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: exact.f

PROCEDURE NAME: EXACT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    35: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: exact.f

PROCEDURE NAME: EXACT
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine exact( i, j, k, u000ijk )
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c
    11:           c   compute the exact solution at (i,j,k)
    12:           c
    13:           c---------------------------------------------------------------------
    14:           
    15:                 use lu_data
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  input parameters
    20:           c---------------------------------------------------------------------
    21:                 integer i, j, k
    22:                 double precision u000ijk(*)
    23:           
    24:           c---------------------------------------------------------------------
    25:           c  local variables
    26:           c---------------------------------------------------------------------
    27:                 integer m
    28:                 double precision xi, eta, zeta
    29:           
    30:                 xi  = ( dble ( i - 1 ) ) / ( nx0 - 1 )
    31:                 eta  = ( dble ( j - 1 ) ) / ( ny0 - 1 )
    32:                 zeta = ( dble ( k - 1 ) ) / ( nz - 1 )
    33:           
    34:           
    35: *------>        do m = 1, 5
    36: |                  u000ijk(m) =  ce(m,1)
    37: |              >        + (ce(m,2)
    38: |              >        + (ce(m,5)
    39: |              >        + (ce(m,8)
    40: |              >        +  ce(m,11) * xi) * xi) * xi) * xi
    41: |              >        + (ce(m,3)
    42: |              >        + (ce(m,6)
    43: |              >        + (ce(m,9)
    44: |              >        +  ce(m,12) * eta) * eta) * eta) * eta
    45: |              >        + (ce(m,4)
    46: |              >        + (ce(m,7)
    47: |              >        + (ce(m,10)
    48: |              >        +  ce(m,13) * zeta) * zeta) * zeta) * zeta
    49: *------         end do
    50:           
    51:                 return
    52:                 end


