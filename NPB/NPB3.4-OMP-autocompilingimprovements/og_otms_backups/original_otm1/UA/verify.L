NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME : verify.f

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fassume-contiguous                      : disable
        -fcopyin-intent-out                      : enable
        -fcse-after-vectorization                : disable
        -ffast-formatted-io                      : enable
        -ffast-math                              : enable
        -fignore-asynchronous                    : disable
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
        -fivdep                                  : disable
        -floop-collapse                          : enable
        -floop-count                             : 5000
        -floop-fusion                            : enable
        -floop-interchange                       : enable
        -floop-normalize                         : enable
        -floop-split                             : disable
        -floop-strip-mine                        : enable
        -floop-unroll                            : enable
        -floop-unroll-complete                   : 4
        -floop-unroll-complete-nest              : 3
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : disable
        -fmove-nested-loop-invariants-outer      : enable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fnamed-noalias-aggressive               : disable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freorder-logical-expression             : enable
        -freplace-loop-equation                  : disable
        -marray-io                               : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mstack-arrays                           : enable
        -mno-stack-arrays                        : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : enable
        -mvector-loop-count-test                 : disable
        -mvector-low-precise-divide-function     : disable
        -mvector-merge-conditional               : enable
        -mvector-packed                          : disable
        -mvector-power-to-explog                 : disable
        -mvector-power-to-sqrt                   : enable
        -mvector-reduction                       : enable
        -mvector-shortloop-reduction             : disable
        -mvector-sqrt-instruction                : disable
        -mvector-threshold                       : 5
        -mwork-vector-kind=none                  : disable

      Parallelization Options : 
        -fopenmp                                 : enable
        -mparallel                               : disable
        -mparallel-innerloop                     : disable
        -mparallel-omp-routine                   : enable
        -mparallel-outerloop-strip-mine          : disable
        -mparallel-sections                      : disable
        -mparallel-threshold                     : 2000
        -mschedule-dynamic                       : disable
        -mschedule-runtime                       : disable
        -mschedule-static                        : disable
        -mschedule-chunk-size                    : 0
        -pthread                                 : disable

      Inlining Options : 
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
        -finline-max-depth                       : 2
        -finline-max-function-size               : 50
        -finline-max-times                       : 6

      Code Generation Options : 
        -finstrument-functions                   : disable
        -fpic | -fPIC                            : disable
        -ftrace                                  : disable
        -p                                       : disable

      Debugging Options : 
        -fbounds-check                           : disable
        -fcheck=all                              : disable
        -fcheck=alias                            : disable
        -fcheck=bits                             : disable
        -fcheck=bounds                           : disable
        -fcheck=dangling                         : disable
        -fcheck=do                               : disable
        -fcheck=iovf                             : disable
        -fcheck=pointer                          : disable
        -fcheck=present                          : disable
        -fcheck=recursion                        : disable
        -g                                       : disable
        -minit-stack                             : not specified
        -mmemory-trace                           : disable
        -mmemory-trace-full                      : disable
        -traceback                               : disable

      Language Options : 
        -bss                                     : disable
        -fdefault-integer                        : 4
        -fdefault-double                         : 8
        -fdefault-real                           : 4
        -fextend-source                          : disable
        -ffree-form                              : disable
        -ffixed-form                             : enable
        -ff90-sign                               : disable
        -fmax-continuation-lines                 : 1023
        -frealloc-lhs                            : enable
        -frealloc-lhs-array                      : enable
        -frealloc-lhs-scalar                     : enable
        -save                                    : disable
        -std                                     : Fortran 2008

      Message Options : 
        -Wall                                    : disable
        -Werror                                  : disable
        -Wextension                              : enable
        -Wobsolescent                            : enable
        -Woverflow                               : disable
        -Woverflow-errors                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -pedantic-errors                         : disable
        -w                                       : disable

      List Output Options : 
        -report-file                             : not specified
        -report-append-mode                      : disable
        -report-all                              : enable
        -report-cg                               : enable
        -report-diagnostics                      : enable
        -report-format                           : enable
        -report-inline                           : enable
        -report-option                           : enable
        -report-vector                           : enable

      Directory Options : 
        -fintrinsic-modules-path                 : /opt/nec/ve/nfort/3.2.0/include
        -module                                  : not specified

      Miscellaneous Options : 
        -v                                       : disable


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:11 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 28 [s0-s12 s15-s16 s18-s19 s23-s30 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1296 bytes
      Register spill area      :    0 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1040 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


