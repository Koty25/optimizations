NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:58 2021
FILE NAME: tzetar.f

PROCEDURE NAME: TZETAR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    23: par(1801): Parallel routine generated.: TZETAR$1
    27: par(1803): Parallelized by "do".
    27: opt(2495): Loop nest collapsed by "collapse" clause.
    27: vec( 101): Vectorized loop.
    29: vec( 101): Vectorized loop.
    52: vec( 128): Fused multiply-add operation applied.
    53: vec( 128): Fused multiply-add operation applied.
    54: vec( 128): Fused multiply-add operation applied.
    55: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:58 2021
FILE NAME: tzetar.f

PROCEDURE NAME: TZETAR
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  tzetar
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
    17:                  integer i, j, k
    18:                  double precision  t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3,
    19:                >                   r4, r5, btuz, ac2u, uzik1
    20:           
    21:           
    22:                  if (timeron) call timer_start(t_tzetar)
    23:           !$omp parallel do default(shared)
    24:           !$omp& private(i,j,k,t1,t2,t3,ac,xvel,yvel,zvel,r1,r2,r3,
    25:           !$omp&              r4,r5,btuz,ac2u,uzik1)
    26:           !$omp&  collapse(2)
    27: Y------>         do    k = 1, nz2
    28: |+----->            do    j = 1, ny2
    29: ||V---->               do    i = 1, nx2
    30: |||       
    31: |||                       xvel = us(i,j,k)
    32: |||                       yvel = vs(i,j,k)
    33: |||                       zvel = ws(i,j,k)
    34: |||                       ac   = speed(i,j,k)
    35: |||       
    36: |||                       ac2u = ac*ac
    37: |||       
    38: |||                       r1 = rhs(1,i,j,k)
    39: |||                       r2 = rhs(2,i,j,k)
    40: |||                       r3 = rhs(3,i,j,k)
    41: |||                       r4 = rhs(4,i,j,k)
    42: |||                       r5 = rhs(5,i,j,k)
    43: |||       
    44: |||                       uzik1 = u(1,i,j,k)
    45: |||                       btuz  = bt * uzik1
    46: |||       
    47: |||                       t1 = btuz/ac * (r4 + r5)
    48: |||                       t2 = r3 + t1
    49: |||                       t3 = btuz * (r4 - r5)
    50: |||       
    51: |||                       rhs(1,i,j,k) = t2
    52: |||     F                 rhs(2,i,j,k) = -uzik1*r2 + xvel*t2
    53: |||     F                 rhs(3,i,j,k) =  uzik1*r1 + yvel*t2
    54: |||     F                 rhs(4,i,j,k) =  zvel*t2  + t3
    55: |||     F                 rhs(5,i,j,k) =  uzik1*(-xvel*r2 + yvel*r1) +
    56: |||            >                    qs(i,j,k)*t2 + c2iv*ac2u*t1 + zvel*t3
    57: |||       
    58: ||V----                end do
    59: |+-----             end do
    60: Y------          end do
    61:                  if (timeron) call timer_stop(t_tzetar)
    62:           
    63:                  return
    64:                  end


