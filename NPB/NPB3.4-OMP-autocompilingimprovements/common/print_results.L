NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME : print_results.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    25: par(1801): Parallel routine generated.: PRINT_RESULTS$1
    26: par(1805): Master-construct generated.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
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

     2:                 subroutine print_results(name, class, n1, n2, n3, niter,
     3:                >               t, mops, optype, verified, npbversion,
     4:                >               compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
     5:            
     6:                 implicit none
     7:                 character name*(*)
     8:                 character class*1
     9:                 integer   n1, n2, n3, niter, j
    10:                 double precision t, mops
    11:                 character optype*24, size*15
    12:                 logical   verified
    13:                 character*(*) npbversion, compiletime,
    14:                >              cs1, cs2, cs3, cs4, cs5, cs6, cs7
    15:                 integer   num_threads, max_threads, i
    16:           c$    integer omp_get_num_threads, omp_get_max_threads
    17:           c$    external omp_get_num_threads, omp_get_max_threads
    18:           
    19:           
    20:                 max_threads = 0
    21:           c$    max_threads = omp_get_max_threads()
    22:           
    23:           c     figure out number of threads used
    24:                 num_threads = 0
    25:           c$omp parallel shared(num_threads)
    26:           c$omp master
    27:           c$    num_threads = omp_get_num_threads()
    28:           c$omp end master
    29:           c$omp end parallel
    30:           
    31:           
    32:                    write (*, 2) name
    33:            2       format(//, ' ', A, ' Benchmark Completed.')
    34:           
    35:                    write (*, 3) Class
    36:            3       format(' Class           = ', 12x, a12)
    37:           
    38:           c   If this is not a grid-based problem (EP, FT, CG), then
    39:           c   we only print n1, which contains some measure of the
    40:           c   problem size. In that case, n2 and n3 are both zero.
    41:           c   Otherwise, we print the grid size n1xn2xn3
    42:           
    43:                    if ((n2 .eq. 0) .and. (n3 .eq. 0)) then
    44:                       if (name(1:2) .eq. 'EP') then
    45:                          write(size, '(f15.0)' ) 2.d0**n1
    46:                          j = 15
    47:                          if (size(j:j) .eq. '.') j = j - 1
    48:                          write (*,42) size(1:j)
    49:            42            format(' Size            = ',9x, a15)
    50:                       else
    51:                          write (*,44) n1
    52:            44            format(' Size            = ',12x, i12)
    53:                       endif
    54:                    else
    55:                       write (*, 4) n1,n2,n3
    56:            4          format(' Size            =  ',9x, i4,'x',i4,'x',i4)
    57:                    endif
    58:           
    59:                    write (*, 5) niter
    60:            5       format(' Iterations      = ', 12x, i12)
    61:            
    62:                    write (*, 6) t
    63:            6       format(' Time in seconds = ',12x, f12.2)
    64:           
    65:                    if (num_threads .gt. 0) write (*,7) num_threads
    66:            7       format(' Total threads   = ', 12x, i12)
    67:            
    68:                    if (max_threads .gt. 0) write (*,8) max_threads
    69:            8       format(' Avail threads   = ', 12x, i12)
    70:           
    71:                    if (num_threads .ne. max_threads) write (*,88)
    72:            88      format(' Warning: Threads used differ from threads available')
    73:           
    74:                    write (*,9) mops
    75:            9       format(' Mop/s total     = ',12x, f12.2)
    76:           
    77:                    if (num_threads .gt. 0) write (*,10) mops/float( num_threads )
    78:            10      format(' Mop/s/thread    = ', 12x, f12.2)
    79:           
    80:                    write(*, 11) optype
    81:            11      format(' Operation type  = ', a24)
    82:           
    83:                    if (verified) then
    84:                       write(*,12) '  SUCCESSFUL'
    85:                    else
    86:                       write(*,12) 'UNSUCCESSFUL'
    87:                    endif
    88:            12      format(' Verification    = ', 12x, a)
    89:           
    90:                    write(*,13) npbversion
    91:            13      format(' Version         = ', 12x, a12)
    92:           
    93:                    write(*,14) compiletime
    94:            14      format(' Compile date    = ', 12x, a12)
    95:           
    96:           
    97:                    write (*,121) cs1
    98:            121     format(/, ' Compile options:', /,
    99:                >          '    F77          = ', A)
   100:           
   101:                    write (*,122) cs2
   102:            122     format('    FLINK        = ', A)
   103:           
   104:                    write (*,123) cs3
   105:            123     format('    F_LIB        = ', A)
   106:           
   107:                    write (*,124) cs4
   108:            124     format('    F_INC        = ', A)
   109:           
   110:                    write (*,125) cs5
   111:            125     format('    FFLAGS       = ', A)
   112:           
   113:                    write (*,126) cs6
   114:            126     format('    FLINKFLAGS   = ', A)
   115:           
   116:                    write(*, 127) cs7
   117:            127     format('    RAND         = ', A)
   118:            
   119:                    write (*,130)
   120:            130     format(//' Please send all errors/feedbacks to:'//
   121:                >            ' NPB Development Team'/
   122:                >            ' npb@nas.nasa.gov'//)
   123:           c 130     format(//' Please send the results of this run to:'//
   124:           c     >            ' NPB Development Team '/
   125:           c     >            ' Internet: npb@nas.nasa.gov'/
   126:           c     >            ' '/
   127:           c     >            ' If email is not available, send this to:'//
   128:           c     >            ' MS T27A-1'/
   129:           c     >            ' NASA Ames Research Center'/
   130:           c     >            ' Moffett Field, CA  94035-1000'//
   131:           c     >            ' Fax: 650-604-3957'//)
   132:           
   133:           
   134:                    return
   135:                    end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0 s2 s8-s12 s15-s16 s23 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


