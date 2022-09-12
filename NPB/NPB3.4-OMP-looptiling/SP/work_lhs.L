NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    45: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:26 2020
FILE NAME: work_lhs.f

PROCEDURE NAME: LHSINIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    32:                  subroutine lhsinit(ni, lhs, lhsp, lhsm)
    33:           
    34:                  implicit none
    35:           
    36:                  integer ni
    37:                  double precision lhs(5,0:*), lhsp(5,0:*), lhsm(5,0:*)
    38:           
    39:                  integer m
    40:           
    41:           c---------------------------------------------------------------------
    42:           c     zap the whole left hand side for starters
    43:           c     set all diagonal values to 1. This is overkill, but convenient
    44:           c---------------------------------------------------------------------
    45: *------>         do   m = 1, 5
    46: |                   lhs (m,0) = 0.0d0
    47: |                   lhsp(m,0) = 0.0d0
    48: |                   lhsm(m,0) = 0.0d0
    49: |                   lhs (m,ni) = 0.0d0
    50: |                   lhsp(m,ni) = 0.0d0
    51: |                   lhsm(m,ni) = 0.0d0
    52: *------          end do
    53:                  lhs (3,0) = 1.0d0
    54:                  lhsp(3,0) = 1.0d0
    55:                  lhsm(3,0) = 1.0d0
    56:                  lhs (3,ni) = 1.0d0
    57:                  lhsp(3,ni) = 1.0d0
    58:                  lhsm(3,ni) = 1.0d0
    59:            
    60:                  return
    61:                  end


