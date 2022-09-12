NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: setiv.f

PROCEDURE NAME: SETIV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    29: par(1801): Parallel routine generated.: SETIV$1
    33: par(1803): Parallelized by "do".
    33: opt(2495): Loop nest collapsed by "collapse" clause.
    37: vec( 110): Vectorization obstructive procedure reference.: EXACT
    37: vec( 103): Unvectorized loop.
    39: opt(1025): Reference to this procedure inhibits optimization.: EXACT
    45: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: setiv.f

PROCEDURE NAME: SETIV
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine setiv
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c
    11:           c   set the initial values of independent variables based on tri-linear
    12:           c   interpolation of boundary values in the computational space.
    13:           c
    14:           c---------------------------------------------------------------------
    15:           
    16:                 use lu_data
    17:                 implicit none
    18:           
    19:           c---------------------------------------------------------------------
    20:           c  local variables
    21:           c---------------------------------------------------------------------
    22:                 integer i, j, k, m
    23:                 double precision  xi, eta, zeta
    24:                 double precision  pxi, peta, pzeta
    25:                 double precision  ue_1jk(5),ue_nx0jk(5),ue_i1k(5),
    26:                >        ue_iny0k(5),ue_ij1(5),ue_ijnz(5)
    27:           
    28:           
    29:           !$omp parallel default(shared) private(i,j,k,m,pxi,peta,pzeta,
    30:           !$omp& xi,eta,zeta,ue_ijnz,ue_ij1,ue_iny0k,ue_i1k,ue_nx0jk,ue_1jk)
    31:           !$omp& shared(nx0,ny0,nz)
    32:           !$omp do schedule(static) collapse(2)
    33: P------>        do k = 2, nz - 1
    34: |+----->           do j = 2, ny - 1
    35: ||                    zeta = ( dble (k-1) ) / (nz-1)
    36: ||                    eta = ( dble (j-1) ) / (ny0-1)
    37: ||+---->              do i = 2, nx - 1
    38: |||                      xi = ( dble (i-1) ) / (nx0-1)
    39: |||                      call exact (1,j,k,ue_1jk)
    40: |||                      call exact (nx0,j,k,ue_nx0jk)
    41: |||                      call exact (i,1,k,ue_i1k)
    42: |||                      call exact (i,ny0,k,ue_iny0k)
    43: |||                      call exact (i,j,1,ue_ij1)
    44: |||                      call exact (i,j,nz,ue_ijnz)
    45: |||*--->                 do m = 1, 5
    46: ||||                        pxi =   ( 1.0d+00 - xi ) * ue_1jk(m)
    47: ||||           >                              + xi   * ue_nx0jk(m)
    48: ||||                        peta =  ( 1.0d+00 - eta ) * ue_i1k(m)
    49: ||||           >                              + eta   * ue_iny0k(m)
    50: ||||                        pzeta = ( 1.0d+00 - zeta ) * ue_ij1(m)
    51: ||||           >                              + zeta   * ue_ijnz(m)
    52: ||||      
    53: ||||                        u( m, i, j, k ) = pxi + peta + pzeta
    54: ||||           >                 - pxi * peta - peta * pzeta - pzeta * pxi
    55: ||||           >                 + pxi * peta * pzeta
    56: ||||      
    57: |||*---                  end do
    58: ||+----               end do
    59: |+-----            end do
    60: P------         end do
    61:           !$omp end do nowait
    62:           !$omp end parallel
    63:           
    64:                 return
    65:                 end


