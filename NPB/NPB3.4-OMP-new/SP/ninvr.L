NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:57 2021
FILE NAME: ninvr.f

PROCEDURE NAME: NINVR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    21: par(1801): Parallel routine generated.: NINVR$1
    23: par(1803): Parallelized by "do".
    23: opt(2495): Loop nest collapsed by "collapse" clause.
    23: vec( 101): Vectorized loop.
    25: vec( 101): Vectorized loop.
    39: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:57 2021
FILE NAME: ninvr.f

PROCEDURE NAME: NINVR
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  ninvr
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c   block-diagonal matrix-vector multiplication
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  integer  i, j, k
    18:                  double precision r1, r2, r3, r4, r5, t1, t2
    19:           
    20:                  if (timeron) call timer_start(t_ninvr)
    21:           !$omp parallel do default(shared) private(i,j,k,r1,r2,r3,r4,r5,t1,t2)
    22:           !$omp&  collapse(2)
    23: Y------>         do k = 1, nz2
    24: |+----->            do j = 1, ny2
    25: ||V---->               do i = 1, nx2
    26: |||       
    27: |||                       r1 = rhs(1,i,j,k)
    28: |||                       r2 = rhs(2,i,j,k)
    29: |||                       r3 = rhs(3,i,j,k)
    30: |||                       r4 = rhs(4,i,j,k)
    31: |||                       r5 = rhs(5,i,j,k)
    32: |||        
    33: |||                       t1 = bt * r3
    34: |||                       t2 = 0.5d0 * ( r4 + r5 )
    35: |||       
    36: |||                       rhs(1,i,j,k) = -r2
    37: |||                       rhs(2,i,j,k) =  r1
    38: |||                       rhs(3,i,j,k) = bt * ( r4 - r5 )
    39: |||     F                 rhs(4,i,j,k) = -t1 + t2
    40: |||                       rhs(5,i,j,k) =  t1 + t2
    41: ||V----                enddo
    42: |+-----             enddo
    43: Y------          enddo
    44:                  if (timeron) call timer_stop(t_ninvr)
    45:           
    46:                  return
    47:                  end


