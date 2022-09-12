NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    21: opt(1418): Constant-length loop is expanded.
    25: par(1801): Parallel routine generated.: ERROR_NORM$1
    25: opt(1418): Constant-length loop is expanded.
    29: par(1803): Parallelized by "do".
    29: opt(2495): Loop nest collapsed by "collapse" clause.
    33: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    33: vec( 103): Unvectorized loop.
    35: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    37: opt(1418): Constant-length loop is expanded.
    44: par(1807): Critical-construct generated.
    44: par(1809): Barrier synchronization.
    44: opt(1418): Constant-length loop is expanded.
    46: opt(1418): Constant-length loop is expanded.
    47: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine error_norm(rms)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c this function computes the norm of the difference between the
    12:           c computed solution and the exact solution
    13:           c---------------------------------------------------------------------
    14:           
    15:                  use sp_data
    16:                  implicit none
    17:           
    18:                  integer i, j, k, m, d
    19:                  double precision xi, eta, zeta, u_exact(5), rms(5), add
    20:           
    21: *------>         do m = 1, 5
    22: |                   rms(m) = 0.0d0
    23: *------          enddo
    24:           
    25:           !$omp parallel do default(shared)
    26:           !$omp&   private(i,j,k,m,zeta,eta,xi,add,u_exact)
    27:           !$omp&   reduction(+:rms)
    28:           !$omp&   schedule(static) collapse(2)
    29: P------>         do   k = 0, grid_points(3)-1
    30: |+----->            do   j = 0, grid_points(2)-1
    31: ||                     zeta = dble(k) * dnzm1
    32: ||                     eta = dble(j) * dnym1
    33: ||+---->               do   i = 0, grid_points(1)-1
    34: |||                       xi = dble(i) * dnxm1
    35: |||                       call exact_solution(xi, eta, zeta, u_exact)
    36: |||       
    37: |||*--->                  do   m = 1, 5
    38: ||||                         add = u(m,i,j,k)-u_exact(m)
    39: ||||                         rms(m) = rms(m) + add*add
    40: |||*---                   end do
    41: ||+----                end do
    42: |+-----             end do
    43: P------          end do
    44:           !$omp end parallel do
    45:           
    46: *------>         do    m = 1, 5
    47: |*----->            do    d = 1, 3
    48: ||                     rms(m) = rms(m) / dble(grid_points(d)-2)
    49: |*-----             end do
    50: |                   rms(m) = dsqrt(rms(m))
    51: *------          end do
    52:           
    53:                  return
    54:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: opt(1418): Constant-length loop is expanded.
    70: par(1801): Parallel routine generated.: RHS_NORM$1
    70: opt(1418): Constant-length loop is expanded.
    73: par(1803): Parallelized by "do".
    73: opt(2495): Loop nest collapsed by "collapse" clause.
    73: vec( 108): Unvectorizable loop structure.
    75: vec( 101): Vectorized loop.
    76: opt(1418): Constant-length loop is expanded.
    78: vec( 128): Fused multiply-add operation applied.
    78: vec( 126): Idiom detected.: SUM
    83: par(1807): Critical-construct generated.
    83: par(1809): Barrier synchronization.
    83: opt(1418): Constant-length loop is expanded.
    85: opt(1418): Constant-length loop is expanded.
    86: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:44 2020
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

    58:                  subroutine rhs_norm(rms)
    59:           
    60:                  use sp_data
    61:                  implicit none
    62:           
    63:                  integer i, j, k, d, m
    64:                  double precision rms(5), add
    65:           
    66: *------>         do m = 1, 5
    67: |                   rms(m) = 0.0d0
    68: *------          enddo
    69:           
    70:           !$omp parallel do default(shared) private(i,j,k,m,add)
    71:           !$omp&   reduction(+:rms)
    72:           !$omp&   schedule(static) collapse(2)
    73: P------>         do k = 1, nz2
    74: |+----->            do j = 1, ny2
    75: ||V---->               do i = 1, nx2
    76: |||*--->                  do m = 1, 5
    77: ||||                         add = rhs(m,i,j,k)
    78: ||||    F                    rms(m) = rms(m) + add*add
    79: |||*---                   end do
    80: ||V----                end do
    81: |+-----             end do
    82: P------          end do
    83:           !$omp end parallel do
    84:           
    85: *------>         do   m = 1, 5
    86: |*----->            do   d = 1, 3
    87: ||                     rms(m) = rms(m) / dble(grid_points(d)-2)
    88: |*-----             end do
    89: |                   rms(m) = dsqrt(rms(m))
    90: *------          end do
    91:           
    92:                  return
    93:                  end


