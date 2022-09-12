NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: initialize.f

PROCEDURE NAME: INITIALIZE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    23: par(1801): Parallel routine generated.: INITIALIZE$1
    32: par(1803): Parallelized by "do".
    32: opt(2495): Loop nest collapsed by "collapse" clause.
    32: vec( 101): Vectorized loop.
    34: vec( 101): Vectorized loop.
    49: par(1803): Parallelized by "do".
    49: opt(2495): Loop nest collapsed by "collapse" clause.
    53: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
    53: vec( 103): Unvectorized loop.
    56: opt(1418): Constant-length loop is expanded.
    58: opt(1025): Reference to this procedure inhibits optimization.: EXACT_SOLUTION
    62: opt(1418): Constant-length loop is expanded.
    68: opt(1418): Constant-length loop is expanded.
    74: opt(1418): Constant-length loop is expanded.
   103: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   103: par(1803): Parallelized by "do".
   103: opt(2495): Loop nest collapsed by "collapse" clause.
   103: vec( 103): Unvectorized loop.
   108: opt(1418): Constant-length loop is expanded.
   122: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   122: par(1803): Parallelized by "do".
   122: opt(2495): Loop nest collapsed by "collapse" clause.
   122: vec( 103): Unvectorized loop.
   127: opt(1418): Constant-length loop is expanded.
   141: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   141: par(1803): Parallelized by "do".
   141: opt(2495): Loop nest collapsed by "collapse" clause.
   141: vec( 103): Unvectorized loop.
   146: opt(1418): Constant-length loop is expanded.
   161: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   161: par(1803): Parallelized by "do".
   161: opt(2495): Loop nest collapsed by "collapse" clause.
   161: vec( 103): Unvectorized loop.
   166: opt(1418): Constant-length loop is expanded.
   180: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   180: par(1803): Parallelized by "do".
   180: opt(2495): Loop nest collapsed by "collapse" clause.
   180: vec( 103): Unvectorized loop.
   185: opt(1418): Constant-length loop is expanded.
   199: vec( 110): Vectorization obstructive procedure reference.: EXACT_SOLUTION
   199: par(1803): Parallelized by "do".
   199: opt(2495): Loop nest collapsed by "collapse" clause.
   199: vec( 103): Unvectorized loop.
   204: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: initialize.f

PROCEDURE NAME: INITIALIZE
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  initialize
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c This subroutine initializes the field variable u using
    12:           c tri-linear transfinite interpolation of the boundary values
    13:           c---------------------------------------------------------------------
    14:           
    15:                  use sp_data
    16:                  implicit none
    17:            
    18:                  integer i, j, k, m, ix, iy, iz
    19:                  double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta,
    20:                >                   Pzeta, temp(5)
    21:            
    22:           
    23:           !$omp parallel default(shared)
    24:           !$omp& private(i,j,k,m,zeta,eta,xi,ix,iy,iz,Pxi,Peta,Pzeta,Pface,temp)
    25:           c---------------------------------------------------------------------
    26:           c  Later (in compute_rhs) we compute 1/u for every element. A few of
    27:           c  the corner elements are not used, but it convenient (and faster)
    28:           c  to compute the whole thing with a simple loop. Make sure those
    29:           c  values are nonzero by initializing the whole thing here.
    30:           c---------------------------------------------------------------------
    31:           !$omp do schedule(static) collapse(2)
    32: Y------>        do k = 0, grid_points(3)-1
    33: |+----->           do j = 0, grid_points(2)-1
    34: ||V---->              do i = 0, grid_points(1)-1
    35: |||                      u(1,i,j,k) = 1.0
    36: |||                      u(2,i,j,k) = 0.0
    37: |||                      u(3,i,j,k) = 0.0
    38: |||                      u(4,i,j,k) = 0.0
    39: |||                      u(5,i,j,k) = 1.0
    40: ||V----               end do
    41: |+-----            end do
    42: Y------         end do
    43:           !$omp end do
    44:           
    45:           c---------------------------------------------------------------------
    46:           c first store the "interpolated" values everywhere on the grid
    47:           c---------------------------------------------------------------------
    48:           !$omp do schedule(static) collapse(2)
    49: P------>            do  k = 0, grid_points(3)-1
    50: |+----->               do  j = 0, grid_points(2)-1
    51: ||                     zeta = dble(k) * dnzm1
    52: ||                        eta = dble(j) * dnym1
    53: ||+---->                  do   i = 0, grid_points(1)-1
    54: |||                          xi = dble(i) * dnxm1
    55: |||        
    56: |||*--->                     do ix = 1, 2
    57: ||||                            Pxi = dble(ix-1)
    58: ||||                            call exact_solution(Pxi, eta, zeta,
    59: ||||           >                                    Pface(1,1,ix))
    60: |||*---                      end do
    61: |||       
    62: |||*--->                     do    iy = 1, 2
    63: ||||                            Peta = dble(iy-1)
    64: ||||                            call exact_solution(xi, Peta, zeta,
    65: ||||           >                                    Pface(1,2,iy))
    66: |||*---                      end do
    67: |||       
    68: |||*--->                     do    iz = 1, 2
    69: ||||                            Pzeta = dble(iz-1)
    70: ||||                            call exact_solution(xi, eta, Pzeta,
    71: ||||           >                                    Pface(1,3,iz))
    72: |||*---                      end do
    73: |||       
    74: |||*--->                     do   m = 1, 5
    75: ||||                            Pxi   = xi   * Pface(m,1,2) +
    76: ||||           >                        (1.0d0-xi)   * Pface(m,1,1)
    77: ||||                            Peta  = eta  * Pface(m,2,2) +
    78: ||||           >                        (1.0d0-eta)  * Pface(m,2,1)
    79: ||||                            Pzeta = zeta * Pface(m,3,2) +
    80: ||||           >                        (1.0d0-zeta) * Pface(m,3,1)
    81: ||||       
    82: ||||                            u(m,i,j,k) = Pxi + Peta + Pzeta -
    83: ||||           >                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
    84: ||||           >                          Pxi*Peta*Pzeta
    85: ||||      
    86: |||*---                      end do
    87: ||+----                   end do
    88: |+-----                end do
    89: P------             end do
    90:           !$omp end do nowait
    91:           
    92:           c---------------------------------------------------------------------
    93:           c now store the exact values on the boundaries
    94:           c---------------------------------------------------------------------
    95:           
    96:           c---------------------------------------------------------------------
    97:           c west face
    98:           c---------------------------------------------------------------------
    99:           
   100:                  xi = 0.0d0
   101:                  i  = 0
   102:           !$omp do schedule(static) collapse(2)
   103: P------>         do  k = 0, grid_points(3)-1
   104: |+----->            do   j = 0, grid_points(2)-1
   105: ||                  zeta = dble(k) * dnzm1
   106: ||                     eta = dble(j) * dnym1
   107: ||                     call exact_solution(xi, eta, zeta, temp)
   108: ||*---->               do   m = 1, 5
   109: |||                       u(m,i,j,k) = temp(m)
   110: ||*----                end do
   111: |+-----             end do
   112: P------          end do
   113:           !$omp end do nowait
   114:           
   115:           c---------------------------------------------------------------------
   116:           c east face
   117:           c---------------------------------------------------------------------
   118:           
   119:                  xi = 1.0d0
   120:                  i  = grid_points(1)-1
   121:           !$omp do schedule(static) collapse(2)
   122: P------>         do   k = 0, grid_points(3)-1
   123: |+----->            do   j = 0, grid_points(2)-1
   124: ||                  zeta = dble(k) * dnzm1
   125: ||                     eta = dble(j) * dnym1
   126: ||                     call exact_solution(xi, eta, zeta, temp)
   127: ||*---->               do   m = 1, 5
   128: |||                       u(m,i,j,k) = temp(m)
   129: ||*----                end do
   130: |+-----             end do
   131: P------          end do
   132:           !$omp end do
   133:           
   134:           c---------------------------------------------------------------------
   135:           c south face
   136:           c---------------------------------------------------------------------
   137:           
   138:                  eta = 0.0d0
   139:                  j   = 0
   140:           !$omp do schedule(static) collapse(2)
   141: P------>         do  k = 0, grid_points(3)-1
   142: |+----->            do   i = 0, grid_points(1)-1
   143: ||                  zeta = dble(k) * dnzm1
   144: ||                     xi = dble(i) * dnxm1
   145: ||                     call exact_solution(xi, eta, zeta, temp)
   146: ||*---->               do   m = 1, 5
   147: |||                       u(m,i,j,k) = temp(m)
   148: ||*----                end do
   149: |+-----             end do
   150: P------          end do
   151:           !$omp end do nowait
   152:           
   153:           
   154:           c---------------------------------------------------------------------
   155:           c north face
   156:           c---------------------------------------------------------------------
   157:           
   158:                  eta = 1.0d0
   159:                  j   = grid_points(2)-1
   160:           !$omp do schedule(static) collapse(2)
   161: P------>         do   k = 0, grid_points(3)-1
   162: |+----->            do   i = 0, grid_points(1)-1
   163: ||                  zeta = dble(k) * dnzm1
   164: ||                     xi = dble(i) * dnxm1
   165: ||                     call exact_solution(xi, eta, zeta, temp)
   166: ||*---->               do   m = 1, 5
   167: |||                       u(m,i,j,k) = temp(m)
   168: ||*----                end do
   169: |+-----             end do
   170: P------          end do
   171:           !$omp end do
   172:           
   173:           c---------------------------------------------------------------------
   174:           c bottom face
   175:           c---------------------------------------------------------------------
   176:           
   177:                  zeta = 0.0d0
   178:                  k    = 0
   179:           !$omp do schedule(static) collapse(2)
   180: P------>         do   j = 0, grid_points(2)-1
   181: |+----->            do   i =0, grid_points(1)-1
   182: ||                  eta = dble(j) * dnym1
   183: ||                     xi = dble(i) *dnxm1
   184: ||                     call exact_solution(xi, eta, zeta, temp)
   185: ||*---->               do   m = 1, 5
   186: |||                       u(m,i,j,k) = temp(m)
   187: ||*----                end do
   188: |+-----             end do
   189: P------          end do
   190:           !$omp end do nowait
   191:           
   192:           c---------------------------------------------------------------------
   193:           c top face
   194:           c---------------------------------------------------------------------
   195:           
   196:                  zeta = 1.0d0
   197:                  k    = grid_points(3)-1
   198:           !$omp do schedule(static) collapse(2)
   199: P------>         do   j = 0, grid_points(2)-1
   200: |+----->            do   i =0, grid_points(1)-1
   201: ||                  eta = dble(j) * dnym1
   202: ||                     xi = dble(i) * dnxm1
   203: ||                     call exact_solution(xi, eta, zeta, temp)
   204: ||*---->               do   m = 1, 5
   205: |||                       u(m,i,j,k) = temp(m)
   206: ||*----                end do
   207: |+-----             end do
   208: P------          end do
   209:           !$omp end do nowait
   210:           !$omp end parallel
   211:           
   212:                  return
   213:                  end


