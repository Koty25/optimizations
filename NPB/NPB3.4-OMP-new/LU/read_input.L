NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:52 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:52 2021
FILE NAME: read_input.f

PROCEDURE NAME: READ_INPUT
FORMAT LIST

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


