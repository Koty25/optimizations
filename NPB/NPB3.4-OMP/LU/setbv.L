NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:52 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    26: par(1801): Parallel routine generated.: SETBV$1
    29: vec( 110): Vectorization obstructive procedure reference.: EXACT
    29: par(1803): Parallelized by "do".
    29: opt(2495): Loop nest collapsed by "collapse" clause.
    29: vec( 103): Unvectorized loop.
    33: opt(1418): Constant-length loop is expanded.
    45: vec( 110): Vectorization obstructive procedure reference.: EXACT
    45: par(1803): Parallelized by "do".
    45: opt(2495): Loop nest collapsed by "collapse" clause.
    45: vec( 103): Unvectorized loop.
    49: opt(1418): Constant-length loop is expanded.
    61: vec( 110): Vectorization obstructive procedure reference.: EXACT
    61: par(1803): Parallelized by "do".
    61: opt(2495): Loop nest collapsed by "collapse" clause.
    61: vec( 103): Unvectorized loop.
    65: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:52 2021
FILE NAME: setbv.f

PROCEDURE NAME: SETBV
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine setbv
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   set the boundary values of dependent variables
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use lu_data
    15:                 implicit none
    16:           
    17:           c---------------------------------------------------------------------
    18:           c   local variables
    19:           c---------------------------------------------------------------------
    20:                 integer i, j, k, m
    21:                 double precision temp1(5), temp2(5)
    22:           
    23:           c---------------------------------------------------------------------
    24:           c   set the dependent variable values along the top and bottom faces
    25:           c---------------------------------------------------------------------
    26:           !$omp parallel default(shared) private(i,j,k,m,temp1,temp2)
    27:           !$omp& shared(nx,ny,nz)
    28:           !$omp do schedule(static) collapse(2)
    29: P------>        do j = 1, ny
    30: |+----->           do i = 1, nx
    31: ||                    call exact( i, j, 1, temp1 )
    32: ||                    call exact( i, j, nz, temp2 )
    33: ||*---->              do m = 1, 5
    34: |||                      u( m, i, j, 1 ) = temp1(m)
    35: |||                      u( m, i, j, nz ) = temp2(m)
    36: ||*----               end do
    37: |+-----            end do
    38: P------         end do
    39:           !$omp end do
    40:           
    41:           c---------------------------------------------------------------------
    42:           c   set the dependent variable values along north and south faces
    43:           c---------------------------------------------------------------------
    44:           !$omp do schedule(static) collapse(2)
    45: P------>        do k = 1, nz
    46: |+----->           do i = 1, nx
    47: ||                    call exact( i, 1, k, temp1 )
    48: ||                    call exact( i, ny, k, temp2 )
    49: ||*---->              do m = 1, 5
    50: |||                      u( m, i, 1, k ) = temp1(m)
    51: |||                      u( m, i, ny, k ) = temp2(m)
    52: ||*----               end do
    53: |+-----            end do
    54: P------         end do
    55:           !$omp end do
    56:           
    57:           c---------------------------------------------------------------------
    58:           c   set the dependent variable values along east and west faces
    59:           c---------------------------------------------------------------------
    60:           !$omp do schedule(static) collapse(2)
    61: P------>        do k = 1, nz
    62: |+----->           do j = 1, ny
    63: ||                    call exact( 1, j, k, temp1 )
    64: ||                    call exact( nx, j, k, temp2 )
    65: ||*---->              do m = 1, 5
    66: |||                      u( m, 1, j, k ) = temp1(m)
    67: |||                      u( m, nx, j, k ) = temp2(m)
    68: ||*----               end do
    69: |+-----            end do
    70: P------         end do
    71:           !$omp end do nowait
    72:           !$omp end parallel
    73:           
    74:                 return
    75:                 end


