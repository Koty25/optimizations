NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:20 2020
FILE NAME: domain.f

PROCEDURE NAME: DOMAIN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:20 2020
FILE NAME: domain.f

PROCEDURE NAME: DOMAIN
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine domain
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:                 use lu_data
    11:                 implicit none
    12:           
    13:           c---------------------------------------------------------------------
    14:           c  local variables
    15:           c---------------------------------------------------------------------
    16:           
    17:           
    18:                 nx = nx0
    19:                 ny = ny0
    20:                 nz = nz0
    21:           
    22:           c---------------------------------------------------------------------
    23:           c   check the sub-domain size
    24:           c---------------------------------------------------------------------
    25:                 if ( ( nx .lt. 4 ) .or.
    26:                >     ( ny .lt. 4 ) .or.
    27:                >     ( nz .lt. 4 ) ) then
    28:                    write (*,2001) nx, ny, nz
    29:            2001    format (5x,'SUBDOMAIN SIZE IS TOO SMALL - ',
    30:                >        /5x,'ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS',
    31:                >        /5x,'SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL',
    32:                >        /5x,'TO 4 THEY ARE CURRENTLY', 3I3)
    33:                    stop
    34:                 end if
    35:           
    36:                 if ( ( nx .gt. isiz1 ) .or.
    37:                >     ( ny .gt. isiz2 ) .or.
    38:                >     ( nz .gt. isiz3 ) ) then
    39:                    write (*,2002) nx, ny, nz
    40:            2002    format (5x,'SUBDOMAIN SIZE IS TOO LARGE - ',
    41:                >        /5x,'ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS',
    42:                >        /5x,'SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO ',
    43:                >        /5x,'ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE',
    44:                >        /5x,'CURRENTLY', 3I4)
    45:                    stop
    46:                 end if
    47:           
    48:           c---------------------------------------------------------------------
    49:           c   set up the start and end in i and j extents for all processors
    50:           c---------------------------------------------------------------------
    51:                 ist = 2
    52:                 iend = nx - 1
    53:           
    54:                 jst = 2
    55:                 jend = ny - 1
    56:           
    57:                 ii1 = 2
    58:                 ii2 = nx0 - 1
    59:                 ji1 = 2
    60:                 ji2 = ny0 - 2
    61:                 ki1 = 3
    62:                 ki2 = nz0 - 1
    63:           
    64:                 return
    65:                 end


