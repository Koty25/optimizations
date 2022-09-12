NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:49 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    37: vec( 101): Vectorized loop.
    38: opt(1418): Constant-length loop is expanded.
    46: vec( 101): Vectorized loop.
    47: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:49 2021
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    24:                 subroutine lhsinit(lhs, ni)
    25:           
    26:           c---------------------------------------------------------------------
    27:           c---------------------------------------------------------------------
    28:            
    29:                 integer i, m, n, ni
    30:                 double precision lhs(5,5,3,0:ni)
    31:           
    32:           c---------------------------------------------------------------------
    33:           c     zero the whole left hand side for starters
    34:           c     set all diagonal values to 1. This is overkill, but convenient
    35:           c---------------------------------------------------------------------
    36:                 i = 0
    37: V------>        do m = 1, 5
    38: |*----->           do n = 1, 5
    39: ||                    lhs(m,n,1,i) = 0.0d0
    40: ||                    lhs(m,n,2,i) = 0.0d0
    41: ||                    lhs(m,n,3,i) = 0.0d0
    42: |*-----            end do
    43: |                  lhs(m,m,2,i) = 1.0d0
    44: V------         end do
    45:                 i = ni
    46: V------>        do m = 1, 5
    47: |*----->           do n = 1, 5
    48: ||                    lhs(m,n,1,i) = 0.0d0
    49: ||                    lhs(m,n,2,i) = 0.0d0
    50: ||                    lhs(m,n,3,i) = 0.0d0
    51: |*-----            end do
    52: |                  lhs(m,m,2,i) = 1.0d0
    53: V------         end do
    54:           
    55:                 return
    56:                 end


