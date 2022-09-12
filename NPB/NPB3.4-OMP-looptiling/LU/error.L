NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: error.f

PROCEDURE NAME: ERROR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    26: opt(1418): Constant-length loop is expanded.
    30: par(1801): Parallel routine generated.: ERROR$1
    30: opt(1418): Constant-length loop is expanded.
    32: par(1803): Parallelized by "do".
    32: opt(2495): Loop nest collapsed by "collapse" clause.
    34: vec( 110): Vectorization obstructive procedure reference.: EXACT
    34: vec( 103): Unvectorized loop.
    35: inl(1212): Source for routine not found.: EXACT
    36: opt(1418): Constant-length loop is expanded.
    43: par(1807): Critical-construct generated.
    43: par(1809): Barrier synchronization.
    43: opt(1418): Constant-length loop is expanded.
    45: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: error.f

PROCEDURE NAME: ERROR
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine error
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c
    11:           c   compute the solution error
    12:           c
    13:           c---------------------------------------------------------------------
    14:           
    15:                 use lu_data
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  local variables
    20:           c---------------------------------------------------------------------
    21:                 integer i, j, k, m
    22:                 double precision  tmp
    23:                 double precision  u000ijk(5)
    24:           
    25:           
    26: *------>        do m = 1, 5
    27: |                  errnm(m) = 0.0d+00
    28: *------         end do
    29:           
    30:           !$omp parallel do schedule(static) collapse(2) default(shared)
    31:           !$omp&  private(i,j,k,m,tmp,u000ijk) reduction(+: errnm)
    32: P------>        do k = 2, nz-1
    33: |+----->           do j = jst, jend
    34: ||+---->              do i = ist, iend
    35: |||                      call exact( i, j, k, u000ijk )
    36: |||*--->                 do m = 1, 5
    37: ||||                        tmp = ( u000ijk(m) - u(m,i,j,k) )
    38: ||||                        errnm(m) = errnm(m) + tmp * tmp
    39: |||*---                  end do
    40: ||+----               end do
    41: |+-----            end do
    42: P------         end do
    43:           !$omp end parallel do
    44:           
    45: *------>        do m = 1, 5
    46: |                  errnm(m) = sqrt ( errnm(m) / ( dble(nx0-2)*(ny0-2)*(nz0-2) ) )
    47: *------         end do
    48:           
    49:           c        write (*,1002) ( errnm(m), m = 1, 5 )
    50:           
    51:            1002 format (1x/1x,'RMS-norm of error in soln. to ',
    52:                > 'first pde  = ',1pe12.5/,
    53:                > 1x,'RMS-norm of error in soln. to ',
    54:                > 'second pde = ',1pe12.5/,
    55:                > 1x,'RMS-norm of error in soln. to ',
    56:                > 'third pde  = ',1pe12.5/,
    57:                > 1x,'RMS-norm of error in soln. to ',
    58:                > 'fourth pde = ',1pe12.5/,
    59:                > 1x,'RMS-norm of error in soln. to ',
    60:                > 'fifth pde  = ',1pe12.5)
    61:           
    62:                 return
    63:                 end


