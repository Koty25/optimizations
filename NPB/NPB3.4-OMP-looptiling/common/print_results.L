NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    25: par(1801): Parallel routine generated.: PRINT_RESULTS$1
    26: par(1805): Master-construct generated.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:18 2020
FILE NAME: print_results.f

PROCEDURE NAME: PRINT_RESULTS
FORMAT LIST

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


