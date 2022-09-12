NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: add.f

PROCEDURE NAME: ADD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    19: par(1801): Parallel routine generated.: ADD$1
    20: par(1803): Parallelized by "do".
    20: opt(2495): Loop nest collapsed by "collapse" clause.
    20: vec( 101): Vectorized loop.
    22: vec( 101): Vectorized loop.
    23: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: add.f

PROCEDURE NAME: ADD
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine  add
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     addition of update to the vector u
    11:           c---------------------------------------------------------------------
    12:           
    13:                 use bt_data
    14:                 implicit none
    15:           
    16:                 integer i, j, k, m
    17:           
    18:                 if (timeron) call timer_start(t_add)
    19:           !$omp parallel do default(shared) private(i,j,k,m) collapse(2)
    20: Y------>        do     k = 1, grid_points(3)-2
    21: |+----->           do     j = 1, grid_points(2)-2
    22: ||V---->              do     i = 1, grid_points(1)-2
    23: |||*--->                 do    m = 1, 5
    24: ||||                        u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k)
    25: |||*---                  enddo
    26: ||V----               enddo
    27: |+-----            enddo
    28: Y------         enddo
    29:                 if (timeron) call timer_stop(t_add)
    30:           
    31:                 return
    32:                 end


