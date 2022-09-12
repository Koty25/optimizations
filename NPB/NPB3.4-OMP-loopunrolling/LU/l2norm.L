NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:40 2020
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: opt(1418): Constant-length loop is expanded.
    41: par(1801): Parallel routine generated.: L2NORM$1
    41: opt(1418): Constant-length loop is expanded.
    43: par(1803): Parallelized by "do".
    43: opt(2495): Loop nest collapsed by "collapse" clause.
    43: vec( 108): Unvectorizable loop structure.
    45: vec( 101): Vectorized loop.
    46: opt(1418): Constant-length loop is expanded.
    47: vec( 128): Fused multiply-add operation applied.
    47: vec( 126): Idiom detected.: SUM
    52: par(1807): Critical-construct generated.
    52: par(1809): Barrier synchronization.
    52: opt(1418): Constant-length loop is expanded.
    54: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:40 2020
FILE NAME: l2norm.f

PROCEDURE NAME: L2NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine l2norm ( ldx, ldy, ldz,
     5:                >                    nx0, ny0, nz0,
     6:                >                    ist, iend,
     7:                >                    jst, jend,
     8:                >                    v, sum )
     9:           c---------------------------------------------------------------------
    10:           c---------------------------------------------------------------------
    11:           
    12:           c---------------------------------------------------------------------
    13:           c   to compute the l2-norm of vector v.
    14:           c---------------------------------------------------------------------
    15:           
    16:                 implicit none
    17:           
    18:           c---------------------------------------------------------------------
    19:           c  input parameters
    20:           c---------------------------------------------------------------------
    21:                 integer ldx, ldy, ldz
    22:                 integer nx0, ny0, nz0
    23:                 integer ist, iend
    24:                 integer jst, jend
    25:           c---------------------------------------------------------------------
    26:           c   To improve cache performance, second two dimensions padded by 1
    27:           c   for even number sizes only.  Only needed in v.
    28:           c---------------------------------------------------------------------
    29:                 double precision  v(5,ldx/2*2+1,ldy/2*2+1,*), sum(5)
    30:           
    31:           c---------------------------------------------------------------------
    32:           c  local variables
    33:           c---------------------------------------------------------------------
    34:                 integer i, j, k, m
    35:           
    36:           
    37: *------>        do m = 1, 5
    38: |                  sum(m) = 0.0d+00
    39: *------         end do
    40:           
    41:           !$omp parallel do schedule(static) collapse(2) default(shared)
    42:           !$omp&  private(i,j,k,m) reduction(+: sum)
    43: P------>        do k = 2, nz0-1
    44: |+----->           do j = jst, jend
    45: ||V---->              do i = ist, iend
    46: |||*--->                 do m = 1, 5
    47: ||||    F                   sum(m) = sum(m) + v(m,i,j,k)*v(m,i,j,k)
    48: |||*---                  end do
    49: ||V----               end do
    50: |+-----            end do
    51: P------         end do
    52:           !$omp end parallel do
    53:           
    54: *------>        do m = 1, 5
    55: |                  sum(m) = sqrt ( sum(m) / ( dble(nx0-2)*(ny0-2)*(nz0-2) ) )
    56: *------         end do
    57:           
    58:                 return
    59:                 end


