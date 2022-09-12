NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME : read_input.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
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

     5:                 subroutine read_input
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:                 use lu_data
    11:                 implicit none
    12:           
    13:                 integer  fstatus
    14:           !$    integer  omp_get_max_threads
    15:           !$    external omp_get_max_threads
    16:           
    17:           
    18:           c---------------------------------------------------------------------
    19:           c    if input file does not exist, it uses defaults
    20:           c       ipr = 1 for detailed progress output
    21:           c       inorm = how often the norm is printed (once every inorm iterations)
    22:           c       itmax = number of pseudo time steps
    23:           c       dt = time step
    24:           c       omega 1 over-relaxation factor for SSOR
    25:           c       tolrsd = steady state residual tolerance levels
    26:           c       nx, ny, nz = number of grid points in x, y, z directions
    27:           c---------------------------------------------------------------------
    28:           
    29:                    write(*, 1000)
    30:           
    31:                    open (unit=3,file='inputlu.data',status='old',
    32:                >         access='sequential',form='formatted', iostat=fstatus)
    33:                    if (fstatus .eq. 0) then
    34:           
    35:                       write(*, *) 'Reading from input file inputlu.data'
    36:           
    37:                       read (3,*)
    38:                       read (3,*)
    39:                       read (3,*) ipr, inorm
    40:                       read (3,*)
    41:                       read (3,*)
    42:                       read (3,*) itmax
    43:                       read (3,*)
    44:                       read (3,*)
    45:                       read (3,*) dt
    46:                       read (3,*)
    47:                       read (3,*)
    48:                       read (3,*) omega
    49:                       read (3,*)
    50:                       read (3,*)
    51:                       read (3,*) tolrsd(1),tolrsd(2),tolrsd(3),tolrsd(4),tolrsd(5)
    52:                       read (3,*)
    53:                       read (3,*)
    54:                       read (3,*) nx0, ny0, nz0
    55:                       close(3)
    56:                    else
    57:                       ipr = ipr_default
    58:                       inorm = inorm_default
    59:                       itmax = itmax_default
    60:                       dt = dt_default
    61:                       omega = omega_default
    62:                       tolrsd(1) = tolrsd1_def
    63:                       tolrsd(2) = tolrsd2_def
    64:                       tolrsd(3) = tolrsd3_def
    65:                       tolrsd(4) = tolrsd4_def
    66:                       tolrsd(5) = tolrsd5_def
    67:                       nx0 = isiz1
    68:                       ny0 = isiz2
    69:                       nz0 = isiz3
    70:                    endif
    71:           
    72:           c---------------------------------------------------------------------
    73:           c   check problem size
    74:           c---------------------------------------------------------------------
    75:           
    76:                    if ( ( nx0 .lt. 4 ) .or.
    77:                >        ( ny0 .lt. 4 ) .or.
    78:                >        ( nz0 .lt. 4 ) ) then
    79:           
    80:                       write (*,2001)
    81:            2001       format (5x,'PROBLEM SIZE IS TOO SMALL - ',
    82:                >           /5x,'SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5')
    83:                       stop
    84:           
    85:                    end if
    86:           
    87:                    if ( ( nx0 .gt. isiz1 ) .or.
    88:                >        ( ny0 .gt. isiz2 ) .or.
    89:                >        ( nz0 .gt. isiz3 ) ) then
    90:           
    91:                       write (*,2002)
    92:            2002       format (5x,'PROBLEM SIZE IS TOO LARGE - ',
    93:                >           /5x,'NX, NY AND NZ SHOULD BE EQUAL TO ',
    94:                >           /5x,'ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY')
    95:                       stop
    96:           
    97:                    end if
    98:           
    99:           
   100:                    write(*, 1001) nx0, ny0, nz0
   101:                    write(*, 1002) itmax
   102:           !$       write(*, 1003) omp_get_max_threads()
   103:                    write(*, *)
   104:           
   105:           
   106:            1000 format(//,' NAS Parallel Benchmarks (NPB3.4-OMP)',
   107:                >          ' - LU Benchmark', /)
   108:            1001    format(' Size: ', i4, 'x', i4, 'x', i4)
   109:            1002    format(' Iterations:                  ', i5)
   110:            1003    format(' Number of available threads: ', i5)
   111:            
   112:           
   113:           
   114:                 return
   115:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:50 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 32 [s0-s12 s15-s16 s23-s33 s58-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1296 bytes
      Register spill area      :    0 bytes
      Parameter area           :   96 bytes
      Register save area       :  176 bytes
      User data area           : 1024 bytes
      Others                   :    0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


