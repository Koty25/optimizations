NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:33 2020
FILE NAME: initialize.f

PROCEDURE NAME: INITIALIZE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    22: par(1801): Parallel routine generated.: INITIALIZE$1
    31: par(1803): Parallelized by "do".
    31: opt(2495): Loop nest collapsed by "collapse" clause.
    31: vec( 101): Vectorized loop.
    33: vec( 101): Vectorized loop.
    34: opt(1418): Constant-length loop is expanded.
    49: par(1803): Parallelized by "do".
    49: opt(2495): Loop nest collapsed by "collapse" clause.
    53: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    53: vec( 103): Unvectorized loop.
    56: opt(1418): Constant-length loop is expanded.
    57: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    61: opt(1418): Constant-length loop is expanded.
    66: opt(1418): Constant-length loop is expanded.
    71: opt(1418): Constant-length loop is expanded.
    99: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    99: par(1803): Parallelized by "do".
    99: opt(2495): Loop nest collapsed by "collapse" clause.
    99: vec( 103): Unvectorized loop.
   104: opt(1418): Constant-length loop is expanded.
   118: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   118: par(1803): Parallelized by "do".
   118: opt(2495): Loop nest collapsed by "collapse" clause.
   118: vec( 103): Unvectorized loop.
   123: opt(1418): Constant-length loop is expanded.
   136: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   136: par(1803): Parallelized by "do".
   136: opt(2495): Loop nest collapsed by "collapse" clause.
   136: vec( 103): Unvectorized loop.
   141: opt(1418): Constant-length loop is expanded.
   155: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   155: par(1803): Parallelized by "do".
   155: opt(2495): Loop nest collapsed by "collapse" clause.
   155: vec( 103): Unvectorized loop.
   160: opt(1418): Constant-length loop is expanded.
   173: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   173: par(1803): Parallelized by "do".
   173: opt(2495): Loop nest collapsed by "collapse" clause.
   173: vec( 103): Unvectorized loop.
   178: opt(1418): Constant-length loop is expanded.
   191: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   191: par(1803): Parallelized by "do".
   191: opt(2495): Loop nest collapsed by "collapse" clause.
   191: vec( 103): Unvectorized loop.
   196: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:33 2020
FILE NAME: initialize.f

PROCEDURE NAME: INITIALIZE
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine  initialize
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c     This subroutine initializes the field variable u using
    11:           c     tri-linear transfinite interpolation of the boundary values
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use bt_data
    15:                 implicit none
    16:            
    17:                 integer i, j, k, m, ix, iy, iz
    18:                 double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta,
    19:                >     Pzeta, temp(5)
    20:           
    21:           
    22:           !$omp parallel default(shared)
    23:           !$omp& private(i,j,k,m,zeta,eta,xi,ix,iy,iz,Pface,Pxi,Peta,Pzeta,temp)
    24:           c---------------------------------------------------------------------
    25:           c  Later (in compute_rhs) we compute 1/u for every element. A few of
    26:           c  the corner elements are not used, but it convenient (and faster)
    27:           c  to compute the whole thing with a simple loop. Make sure those
    28:           c  values are nonzero by initializing the whole thing here.
    29:           c---------------------------------------------------------------------
    30:           !$omp do schedule(static) collapse(2)
    31: Y------>        do k = 0, grid_points(3)-1
    32: |+----->           do j = 0, grid_points(2)-1
    33: ||V---->              do i = 0, grid_points(1)-1
    34: |||*--->                 do m = 1, 5
    35: ||||                        u(m,i,j,k) = 1.0
    36: |||*---                  end do
    37: ||V----               end do
    38: |+-----            end do
    39: Y------         end do
    40:           !$omp end do nowait
    41:           c---------------------------------------------------------------------
    42:           
    43:           
    44:           c---------------------------------------------------------------------
    45:           c     first store the "interpolated" values everywhere on the grid
    46:           c---------------------------------------------------------------------
    47:           
    48:           !$omp do schedule(static) collapse(2)
    49: P------>        do k = 0, grid_points(3)-1
    50: |+----->           do j = 0, grid_points(2)-1
    51: ||                    zeta = dble(k) * dnzm1
    52: ||                    eta = dble(j) * dnym1
    53: ||+---->              do i = 0, grid_points(1)-1
    54: |||                      xi = dble(i) * dnxm1
    55: |||        
    56: |||*--->                 do ix = 1, 2
    57: ||||                        call exact_solution(dble(ix-1), eta, zeta,
    58: ||||           >                    Pface(1,1,ix))
    59: |||*---                  enddo
    60: |||       
    61: |||*--->                 do iy = 1, 2
    62: ||||                        call exact_solution(xi, dble(iy-1) , zeta,
    63: ||||           >                    Pface(1,2,iy))
    64: |||*---                  enddo
    65: |||       
    66: |||*--->                 do iz = 1, 2
    67: ||||                        call exact_solution(xi, eta, dble(iz-1),
    68: ||||           >                    Pface(1,3,iz))
    69: |||*---                  enddo
    70: |||       
    71: |||*--->                 do m = 1, 5
    72: ||||                        Pxi   = xi   * Pface(m,1,2) +
    73: ||||           >                    (1.0d0-xi)   * Pface(m,1,1)
    74: ||||                        Peta  = eta  * Pface(m,2,2) +
    75: ||||           >                    (1.0d0-eta)  * Pface(m,2,1)
    76: ||||                        Pzeta = zeta * Pface(m,3,2) +
    77: ||||           >                    (1.0d0-zeta) * Pface(m,3,1)
    78: ||||       
    79: ||||                        u(m,i,j,k) = Pxi + Peta + Pzeta -
    80: ||||           >                    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
    81: ||||           >                    Pxi*Peta*Pzeta
    82: ||||      
    83: |||*---                  enddo
    84: ||+----               enddo
    85: |+-----            enddo
    86: P------         enddo
    87:           !$omp end do nowait
    88:           
    89:           c---------------------------------------------------------------------
    90:           c     now store the exact values on the boundaries
    91:           c---------------------------------------------------------------------
    92:           
    93:           c---------------------------------------------------------------------
    94:           c     west face
    95:           c---------------------------------------------------------------------
    96:                 i = 0
    97:                 xi = 0.0d0
    98:           !$omp do schedule(static) collapse(2)
    99: P------>        do k = 0, grid_points(3)-1
   100: |+----->           do j = 0, grid_points(2)-1
   101: ||                    zeta = dble(k) * dnzm1
   102: ||                    eta = dble(j) * dnym1
   103: ||                    call exact_solution(xi, eta, zeta, temp)
   104: ||*---->              do m = 1, 5
   105: |||                      u(m,i,j,k) = temp(m)
   106: ||*----               enddo
   107: |+-----            enddo
   108: P------         enddo
   109:           !$omp end do nowait
   110:           
   111:           c---------------------------------------------------------------------
   112:           c     east face
   113:           c---------------------------------------------------------------------
   114:           
   115:                 i = grid_points(1)-1
   116:                 xi = 1.0d0
   117:           !$omp do schedule(static) collapse(2)
   118: P------>        do k = 0, grid_points(3)-1
   119: |+----->           do j = 0, grid_points(2)-1
   120: ||                    zeta = dble(k) * dnzm1
   121: ||                    eta = dble(j) * dnym1
   122: ||                    call exact_solution(xi, eta, zeta, temp)
   123: ||*---->              do m = 1, 5
   124: |||                      u(m,i,j,k) = temp(m)
   125: ||*----               enddo
   126: |+-----            enddo
   127: P------         enddo
   128:           !$omp end do
   129:           
   130:           c---------------------------------------------------------------------
   131:           c     south face
   132:           c---------------------------------------------------------------------
   133:                 j = 0
   134:                 eta = 0.0d0
   135:           !$omp do schedule(static) collapse(2)
   136: P------>        do k = 0, grid_points(3)-1
   137: |+----->           do i = 0, grid_points(1)-1
   138: ||                    zeta = dble(k) * dnzm1
   139: ||                    xi = dble(i) * dnxm1
   140: ||                    call exact_solution(xi, eta, zeta, temp)
   141: ||*---->              do m = 1, 5
   142: |||                      u(m,i,j,k) = temp(m)
   143: ||*----               enddo
   144: |+-----            enddo
   145: P------         enddo
   146:           !$omp end do nowait
   147:           
   148:           
   149:           c---------------------------------------------------------------------
   150:           c     north face
   151:           c---------------------------------------------------------------------
   152:                 j = grid_points(2)-1
   153:                 eta = 1.0d0
   154:           !$omp do schedule(static) collapse(2)
   155: P------>        do k = 0, grid_points(3)-1
   156: |+----->           do i = 0, grid_points(1)-1
   157: ||                    zeta = dble(k) * dnzm1
   158: ||                    xi = dble(i) * dnxm1
   159: ||                    call exact_solution(xi, eta, zeta, temp)
   160: ||*---->              do m = 1, 5
   161: |||                      u(m,i,j,k) = temp(m)
   162: ||*----               enddo
   163: |+-----            enddo
   164: P------         enddo
   165:           !$omp end do
   166:           
   167:           c---------------------------------------------------------------------
   168:           c     bottom face
   169:           c---------------------------------------------------------------------
   170:                 k = 0
   171:                 zeta = 0.0d0
   172:           !$omp do schedule(static) collapse(2)
   173: P------>        do j = 0, grid_points(2)-1
   174: |+----->           do i =0, grid_points(1)-1
   175: ||                    eta = dble(j) * dnym1
   176: ||                    xi = dble(i) *dnxm1
   177: ||                    call exact_solution(xi, eta, zeta, temp)
   178: ||*---->              do m = 1, 5
   179: |||                      u(m,i,j,k) = temp(m)
   180: ||*----               enddo
   181: |+-----            enddo
   182: P------         enddo
   183:           !$omp end do nowait
   184:           
   185:           c---------------------------------------------------------------------
   186:           c     top face
   187:           c---------------------------------------------------------------------
   188:                 k = grid_points(3)-1
   189:                 zeta = 1.0d0
   190:           !$omp do schedule(static) collapse(2)
   191: P------>        do j = 0, grid_points(2)-1
   192: |+----->           do i =0, grid_points(1)-1
   193: ||                    eta = dble(j) * dnym1
   194: ||                    xi = dble(i) * dnxm1
   195: ||                    call exact_solution(xi, eta, zeta, temp)
   196: ||*---->              do m = 1, 5
   197: |||                      u(m,i,j,k) = temp(m)
   198: ||*----               enddo
   199: |+-----            enddo
   200: P------         enddo
   201:           !$omp end do nowait
   202:           !$omp end parallel
   203:           
   204:                 return
   205:                 end


