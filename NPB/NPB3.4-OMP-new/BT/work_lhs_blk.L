NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:56 2020
FILE NAME: work_lhs_blk.f

PROCEDURE NAME: LHSINIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    48: opt(1589): Outer loop moved inside inner loop(s).: JB
    48: vec( 101): Vectorized loop.
    52: opt(1418): Constant-length loop is expanded.
    54: opt(1034): Multiple store conflict.: LHSB
    54: opt(1033): Potential multiple store conflict -- use directive if OK.: LHSB
    66: vec( 101): Vectorized loop.
    68: opt(1418): Constant-length loop is expanded.
    70: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:26:56 2020
FILE NAME: work_lhs_blk.f

PROCEDURE NAME: LHSINIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    29:                 subroutine lhsinit(ni)
    30:           
    31:           c---------------------------------------------------------------------
    32:           c---------------------------------------------------------------------
    33:                 use work_lhs
    34:                 implicit none
    35:           
    36:                 integer ni
    37:           
    38:                 integer i, m, n, jb
    39:           
    40:           c---------------------------------------------------------------------
    41:           c     zero the whole left hand side for starters
    42:           c     set all diagonal values to 1. This is overkill, but convenient
    43:           c---------------------------------------------------------------------
    44:                 if (ni .gt. 0) goto 20
    45:           
    46: +------>        do i = 0, 2, 2
    47: |         !dir$ vector always
    48: |X----->           do jb = 1, bsize
    49: ||        !dir$ unroll
    50: ||*---->              do m = 1, 5
    51: |||       !dir$ unroll
    52: |||*--->                 do n = 1, 5
    53: ||||                        lhsa(jb,m,n,i) = 0.0d0
    54: ||||                        lhsb(jb,m,n,i) = 0.0d0
    55: |||*---                  end do
    56: |||                      lhsb(jb,m,m,i) = 1.0d0
    57: ||*----               end do
    58: |X-----            end do
    59: +------         end do
    60:           
    61:                 return
    62:           
    63:             20  continue
    64: +------>        do i = 0, ni, ni
    65: |         !dir$ vector always
    66: |V----->           do jb = 1, bsize
    67: ||        !dir$ unroll
    68: ||*---->              do m = 1, 5
    69: |||       !dir$ unroll
    70: |||*--->                 do n = 1, 5
    71: ||||                        lhsc(jb,m,n,i) = 0.0d0
    72: |||*---                  end do
    73: ||*----               end do
    74: |V-----            end do
    75: +------         end do
    76:           
    77:                 return
    78:                 end


