NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    23: par(1801): Parallel routine generated.: TXINVR$1
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    26: vec( 101): Vectorized loop.
    28: vec( 101): Vectorized loop.
    43: vec( 128): Fused multiply-add operation applied.
    45: vec( 128): Fused multiply-add operation applied.
    49: vec( 128): Fused multiply-add operation applied.
    50: vec( 128): Fused multiply-add operation applied.
    51: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: txinvr.f

PROCEDURE NAME: TXINVR
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  txinvr
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c block-diagonal matrix-vector multiplication
    12:           c---------------------------------------------------------------------
    13:           
    14:                  use sp_data
    15:                  implicit none
    16:           
    17:                  integer i, j, k
    18:                  double precision t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3,
    19:                >                  r4, r5, ac2inv
    20:           
    21:           
    22:                  if (timeron) call timer_start(t_txinvr)
    23:           !$omp parallel do default(shared)
    24:           !$omp& private(i,j,k,t1,t2,t3,ac,ru1,uu,vv,ww,r1,r2,r3,r4,r5,ac2inv)
    25:           !$omp&  collapse(2)
    26: Y------>         do    k = 1, nz2
    27: |+----->            do    j = 1, ny2
    28: ||V---->               do    i = 1, nx2
    29: |||       
    30: |||                       ru1 = rho_i(i,j,k)
    31: |||                       uu = us(i,j,k)
    32: |||                       vv = vs(i,j,k)
    33: |||                       ww = ws(i,j,k)
    34: |||                       ac = speed(i,j,k)
    35: |||                       ac2inv = ac*ac
    36: |||       
    37: |||                       r1 = rhs(1,i,j,k)
    38: |||                       r2 = rhs(2,i,j,k)
    39: |||                       r3 = rhs(3,i,j,k)
    40: |||                       r4 = rhs(4,i,j,k)
    41: |||                       r5 = rhs(5,i,j,k)
    42: |||       
    43: |||     F                 t1 = c2 / ac2inv * ( qs(i,j,k)*r1 - uu*r2  -
    44: |||            >                  vv*r3 - ww*r4 + r5 )
    45: |||     F                 t2 = bt * ru1 * ( uu * r1 - r2 )
    46: |||                       t3 = ( bt * ru1 * ac ) * t1
    47: |||       
    48: |||                       rhs(1,i,j,k) = r1 - t1
    49: |||     F                 rhs(2,i,j,k) = - ru1 * ( ww*r1 - r4 )
    50: |||     F                 rhs(3,i,j,k) =   ru1 * ( vv*r1 - r3 )
    51: |||     F                 rhs(4,i,j,k) = - t2 + t3
    52: |||                       rhs(5,i,j,k) =   t2 + t3
    53: |||       
    54: ||V----                end do
    55: |+-----             end do
    56: Y------          end do
    57:                  if (timeron) call timer_stop(t_txinvr)
    58:           
    59:                  return
    60:                  end


