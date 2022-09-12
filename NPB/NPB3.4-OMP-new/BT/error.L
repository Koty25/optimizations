NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    20: opt(1418): Constant-length loop is expanded.
    24: par(1801): Parallel routine generated.: ERROR_NORM$1
    24: opt(1418): Constant-length loop is expanded.
    26: par(1803): Parallelized by "do".
    26: opt(2495): Loop nest collapsed by "collapse" clause.
    30: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    30: vec( 103): Unvectorized loop.
    32: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    34: opt(1418): Constant-length loop is expanded.
    41: par(1807): Critical-construct generated.
    41: par(1809): Barrier synchronization.
    41: opt(1418): Constant-length loop is expanded.
    43: opt(1418): Constant-length loop is expanded.
    44: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: error.f

PROCEDURE NAME: ERROR_NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine error_norm(rms)
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     this function computes the norm of the difference between the
    11:           c     computed solution and the exact solution
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use bt_data
    15:                 implicit none
    16:           
    17:                 integer i, j, k, m, d
    18:                 double precision xi, eta, zeta, u_exact(5), rms(5), add
    19:           
    20: *------>        do m = 1, 5
    21: |                  rms(m) = 0.0d0
    22: *------         enddo
    23:           
    24:           !$omp parallel do schedule(static) collapse(2) default(shared)
    25:           !$omp& private(i,j,k,m,zeta,eta,xi,add,u_exact) reduction(+: rms)
    26: P------>        do k = 0, grid_points(3)-1
    27: |+----->           do j = 0, grid_points(2)-1
    28: ||                    zeta = dble(k) * dnzm1
    29: ||                    eta = dble(j) * dnym1
    30: ||+---->              do i = 0, grid_points(1)-1
    31: |||                      xi = dble(i) * dnxm1
    32: |||                      call exact_solution(xi, eta, zeta, u_exact)
    33: |||       
    34: |||*--->                 do m = 1, 5
    35: ||||                        add = u(m,i,j,k)-u_exact(m)
    36: ||||                        rms(m) = rms(m) + add*add
    37: |||*---                  enddo
    38: ||+----               enddo
    39: |+-----            enddo
    40: P------         enddo
    41:           !$omp end parallel do
    42:           
    43: *------>        do m = 1, 5
    44: |*----->           do d = 1, 3
    45: ||                    rms(m) = rms(m) / dble(grid_points(d)-2)
    46: |*-----            enddo
    47: |                  rms(m) = dsqrt(rms(m))
    48: *------         enddo
    49:           
    50:                 return
    51:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    68: opt(1418): Constant-length loop is expanded.
    72: par(1801): Parallel routine generated.: RHS_NORM$1
    72: opt(1418): Constant-length loop is expanded.
    74: par(1803): Parallelized by "do".
    74: opt(2495): Loop nest collapsed by "collapse" clause.
    74: vec( 108): Unvectorizable loop structure.
    76: vec( 101): Vectorized loop.
    77: opt(1418): Constant-length loop is expanded.
    79: vec( 128): Fused multiply-add operation applied.
    79: vec( 126): Idiom detected.: SUM
    84: par(1807): Critical-construct generated.
    84: par(1809): Barrier synchronization.
    84: opt(1418): Constant-length loop is expanded.
    86: opt(1418): Constant-length loop is expanded.
    87: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: error.f

PROCEDURE NAME: RHS_NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

    57:                 subroutine rhs_norm(rms)
    58:           
    59:           c---------------------------------------------------------------------
    60:           c---------------------------------------------------------------------
    61:           
    62:                 use bt_data
    63:                 implicit none
    64:           
    65:                 integer i, j, k, d, m
    66:                 double precision rms(5), add
    67:           
    68: *------>        do m = 1, 5
    69: |                  rms(m) = 0.0d0
    70: *------         enddo
    71:           
    72:           !$omp parallel do schedule(static) collapse(2)
    73:           !$omp&  default(shared) private(i,j,k,m,add) reduction(+: rms)
    74: P------>        do k = 1, grid_points(3)-2
    75: |+----->           do j = 1, grid_points(2)-2
    76: ||V---->              do i = 1, grid_points(1)-2
    77: |||*--->                 do m = 1, 5
    78: ||||                        add = rhs(m,i,j,k)
    79: ||||    F                   rms(m) = rms(m) + add*add
    80: |||*---                  enddo
    81: ||V----               enddo
    82: |+-----            enddo
    83: P------         enddo
    84:           !$omp end parallel do
    85:           
    86: *------>        do m = 1, 5
    87: |*----->           do d = 1, 3
    88: ||                    rms(m) = rms(m) / dble(grid_points(d)-2)
    89: |*-----            enddo
    90: |                  rms(m) = dsqrt(rms(m))
    91: *------         enddo
    92:           
    93:                 return
    94:                 end


