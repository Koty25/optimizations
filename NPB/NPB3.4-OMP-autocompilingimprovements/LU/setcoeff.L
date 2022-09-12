NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME : setcoeff.f

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
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
        -fignore-induction-variable-overflow     : enable
        -fignore-volatile                        : enable
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
        -fmove-loop-invariants-unsafe            : enable
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
        -freplace-loop-equation                  : enable
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
        -finline-copy-arguments                  : disable
        -finline-functions                       : enable
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME: setcoeff.f

PROCEDURE NAME: SETCOEFF
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME: setcoeff.f

PROCEDURE NAME: SETCOEFF
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

     5:                 subroutine setcoeff
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
    18:           c---------------------------------------------------------------------
    19:           c   set up coefficients
    20:           c---------------------------------------------------------------------
    21:                 dxi = 1.0d+00 / ( nx0 - 1 )
    22:                 deta = 1.0d+00 / ( ny0 - 1 )
    23:                 dzeta = 1.0d+00 / ( nz0 - 1 )
    24:           
    25:                 tx1 = 1.0d+00 / ( dxi * dxi )
    26:                 tx2 = 1.0d+00 / ( 2.0d+00 * dxi )
    27:                 tx3 = 1.0d+00 / dxi
    28:           
    29:                 ty1 = 1.0d+00 / ( deta * deta )
    30:                 ty2 = 1.0d+00 / ( 2.0d+00 * deta )
    31:                 ty3 = 1.0d+00 / deta
    32:           
    33:                 tz1 = 1.0d+00 / ( dzeta * dzeta )
    34:                 tz2 = 1.0d+00 / ( 2.0d+00 * dzeta )
    35:                 tz3 = 1.0d+00 / dzeta
    36:           
    37:           c---------------------------------------------------------------------
    38:           c   diffusion coefficients
    39:           c---------------------------------------------------------------------
    40:                 dx1 = 0.75d+00
    41:                 dx2 = dx1
    42:                 dx3 = dx1
    43:                 dx4 = dx1
    44:                 dx5 = dx1
    45:           
    46:                 dy1 = 0.75d+00
    47:                 dy2 = dy1
    48:                 dy3 = dy1
    49:                 dy4 = dy1
    50:                 dy5 = dy1
    51:           
    52:                 dz1 = 1.00d+00
    53:                 dz2 = dz1
    54:                 dz3 = dz1
    55:                 dz4 = dz1
    56:                 dz5 = dz1
    57:           
    58:           c---------------------------------------------------------------------
    59:           c   fourth difference dissipation
    60:           c---------------------------------------------------------------------
    61:                 dssp = ( max (dx1, dy1, dz1 ) ) / 4.0d+00
    62:           
    63:           c---------------------------------------------------------------------
    64:           c   coefficients of the exact solution to the first pde
    65:           c---------------------------------------------------------------------
    66:                 ce(1,1) = 2.0d+00
    67:                 ce(1,2) = 0.0d+00
    68:                 ce(1,3) = 0.0d+00
    69:                 ce(1,4) = 4.0d+00
    70:                 ce(1,5) = 5.0d+00
    71:                 ce(1,6) = 3.0d+00
    72:                 ce(1,7) = 5.0d-01
    73:                 ce(1,8) = 2.0d-02
    74:                 ce(1,9) = 1.0d-02
    75:                 ce(1,10) = 3.0d-02
    76:                 ce(1,11) = 5.0d-01
    77:                 ce(1,12) = 4.0d-01
    78:                 ce(1,13) = 3.0d-01
    79:           
    80:           c---------------------------------------------------------------------
    81:           c   coefficients of the exact solution to the second pde
    82:           c---------------------------------------------------------------------
    83:                 ce(2,1) = 1.0d+00
    84:                 ce(2,2) = 0.0d+00
    85:                 ce(2,3) = 0.0d+00
    86:                 ce(2,4) = 0.0d+00
    87:                 ce(2,5) = 1.0d+00
    88:                 ce(2,6) = 2.0d+00
    89:                 ce(2,7) = 3.0d+00
    90:                 ce(2,8) = 1.0d-02
    91:                 ce(2,9) = 3.0d-02
    92:                 ce(2,10) = 2.0d-02
    93:                 ce(2,11) = 4.0d-01
    94:                 ce(2,12) = 3.0d-01
    95:                 ce(2,13) = 5.0d-01
    96:           
    97:           c---------------------------------------------------------------------
    98:           c   coefficients of the exact solution to the third pde
    99:           c---------------------------------------------------------------------
   100:                 ce(3,1) = 2.0d+00
   101:                 ce(3,2) = 2.0d+00
   102:                 ce(3,3) = 0.0d+00
   103:                 ce(3,4) = 0.0d+00
   104:                 ce(3,5) = 0.0d+00
   105:                 ce(3,6) = 2.0d+00
   106:                 ce(3,7) = 3.0d+00
   107:                 ce(3,8) = 4.0d-02
   108:                 ce(3,9) = 3.0d-02
   109:                 ce(3,10) = 5.0d-02
   110:                 ce(3,11) = 3.0d-01
   111:                 ce(3,12) = 5.0d-01
   112:                 ce(3,13) = 4.0d-01
   113:           
   114:           c---------------------------------------------------------------------
   115:           c   coefficients of the exact solution to the fourth pde
   116:           c---------------------------------------------------------------------
   117:                 ce(4,1) = 2.0d+00
   118:                 ce(4,2) = 2.0d+00
   119:                 ce(4,3) = 0.0d+00
   120:                 ce(4,4) = 0.0d+00
   121:                 ce(4,5) = 0.0d+00
   122:                 ce(4,6) = 2.0d+00
   123:                 ce(4,7) = 3.0d+00
   124:                 ce(4,8) = 3.0d-02
   125:                 ce(4,9) = 5.0d-02
   126:                 ce(4,10) = 4.0d-02
   127:                 ce(4,11) = 2.0d-01
   128:                 ce(4,12) = 1.0d-01
   129:                 ce(4,13) = 3.0d-01
   130:           
   131:           c---------------------------------------------------------------------
   132:           c   coefficients of the exact solution to the fifth pde
   133:           c---------------------------------------------------------------------
   134:                 ce(5,1) = 5.0d+00
   135:                 ce(5,2) = 4.0d+00
   136:                 ce(5,3) = 3.0d+00
   137:                 ce(5,4) = 2.0d+00
   138:                 ce(5,5) = 1.0d-01
   139:                 ce(5,6) = 4.0d-01
   140:                 ce(5,7) = 3.0d-01
   141:                 ce(5,8) = 5.0d-02
   142:                 ce(5,9) = 4.0d-02
   143:                 ce(5,10) = 3.0d-02
   144:                 ce(5,11) = 1.0d-01
   145:                 ce(5,12) = 3.0d-01
   146:                 ce(5,13) = 2.0d-01
   147:           
   148:                 return
   149:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME: setcoeff.f

PROCEDURE NAME: SETCOEFF
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME: setcoeff.f

PROCEDURE NAME: SETCOEFF
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:51 2021
FILE NAME: setcoeff.f

PROCEDURE NAME: SETCOEFF
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 25 [s8-s11 s43-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


