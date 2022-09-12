NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:27 2020
FILE NAME: add.f

PROCEDURE NAME: ADD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    19: inl(1212): Source for routine not found.: TIMER_START
    20: par(1801): Parallel routine generated.: ADD$1
    21: par(1803): Parallelized by "do".
    21: opt(2495): Loop nest collapsed by "collapse" clause.
    21: vec( 101): Vectorized loop.
    23: vec( 101): Vectorized loop.
    24: opt(1418): Constant-length loop is expanded.
    30: inl(1212): Source for routine not found.: TIMER_STOP


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:27 2020
FILE NAME: add.f

PROCEDURE NAME: ADD
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  add
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c addition of update to the vector u
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  integer i,j,k,m
    18:           
    19:                  if (timeron) call timer_start(t_add)
    20:           !$omp parallel do default(shared) private(i,j,k,m) collapse(2)
    21: Y------>         do k = 1, nz2
    22: |+----->            do j = 1, ny2
    23: ||V---->               do i = 1, nx2
    24: |||*--->                  do m = 1, 5
    25: ||||                         u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k)
    26: |||*---                   end do
    27: ||V----                end do
    28: |+-----             end do
    29: Y------          end do
    30:                  if (timeron) call timer_stop(t_add)
    31:           
    32:                  return
    33:                  end


