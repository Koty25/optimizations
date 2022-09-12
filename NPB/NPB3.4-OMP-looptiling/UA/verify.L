NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:33 2020
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    18: inl(1212): Source for routine not found.: CALC_NORM
    67: inl(1212): Source for routine not found.: IEEE_ARITHMETIC::IEEE_IS_NAN


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:33 2020
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
FORMAT LIST

 LINE   LOOP      STATEMENT

     1:                 subroutine verify(class, verified)
     2:           
     3:                 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
     4:           
     5:                 use ua_data
     6:           
     7:                 implicit none
     8:           
     9:                 double precision norm, calc_norm, epsilon, norm_dif, norm_ref
    10:                 external         calc_norm
    11:                 character        class
    12:                 logical          verified
    13:            
    14:           c.....tolerance level
    15:                 epsilon = 1.0d-08
    16:           
    17:           c.....compute the temperature integral over the whole domain
    18:                 norm = calc_norm()
    19:           
    20:                 verified = .true.
    21:                 if     ( class .eq. 'S' ) then
    22:                   norm_ref = 0.1890013110962D-02
    23:                 elseif ( class .eq. 'W' ) then
    24:                   norm_ref = 0.2569794837076D-04
    25:                 elseif ( class .eq. 'A' ) then
    26:                   norm_ref = 0.8939996281443D-04
    27:                 elseif ( class .eq. 'B' ) then
    28:                   norm_ref = 0.4507561922901D-04
    29:                 elseif ( class .eq. 'C' ) then
    30:                   norm_ref = 0.1544736587100D-04
    31:                 elseif ( class .eq. 'D' ) then
    32:                   norm_ref = 0.1577586272355D-05
    33:                 else
    34:                   class = 'U'
    35:                   norm_ref = 1.d0
    36:                   verified = .false.
    37:                 endif
    38:           
    39:                 norm_dif = dabs((norm - norm_ref)/norm_ref)
    40:           
    41:           c---------------------------------------------------------------------
    42:           c    Output the comparison of computed results to known cases.
    43:           c---------------------------------------------------------------------
    44:           
    45:                 print *
    46:           
    47:                 if (class .ne. 'U') then
    48:                    write(*, 1990) class
    49:            1990    format(' Verification being performed for class ', a)
    50:                    write (*,2000) epsilon
    51:            2000    format(' accuracy setting for epsilon = ', E20.13)
    52:                 else
    53:                    write(*, 1995)
    54:            1995    format(' Unknown class')
    55:                 endif
    56:           
    57:                 if (class .ne. 'U') then
    58:                    write (*,2001)
    59:                 else
    60:                    write (*, 2005)
    61:                 endif
    62:           
    63:            2001 format(' Comparison of temperature integrals')
    64:            2005 format(' Temperature integral')
    65:                 if (class .eq. 'U') then
    66:                    write(*, 2015) norm
    67:                 else if ((.not.ieee_is_nan(norm_dif)) .and.
    68:                &         norm_dif .le. epsilon) then
    69:                    write (*,2011) norm, norm_ref, norm_dif
    70:                 else
    71:                    verified = .false.
    72:                    write (*,2010) norm, norm_ref, norm_dif
    73:                 endif
    74:           
    75:            2010 format(' FAILURE: ', E20.13, E20.13, E20.13)
    76:            2011 format('          ', E20.13, E20.13, E20.13)
    77:            2015 format('          ', E20.13)
    78:            
    79:                 if (class .eq. 'U') then
    80:                   write(*, 2022)
    81:                   write(*, 2023)
    82:            2022   format(' No reference values provided')
    83:            2023   format(' No verification performed')
    84:                 else if (verified) then
    85:                   write(*, 2020)
    86:            2020   format(' Verification Successful')
    87:                 else
    88:                   write(*, 2021)
    89:            2021   format(' Verification failed')
    90:                 endif
    91:           
    92:                 return
    93:                 end


