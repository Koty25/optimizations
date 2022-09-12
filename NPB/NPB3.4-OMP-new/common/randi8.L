NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: randi8.f

PROCEDURE NAME: RANDLC
FORMAT LIST

 LINE   LOOP      STATEMENT

     1:                 double precision function randlc(x, a)
     2:           
     3:           c---------------------------------------------------------------------
     4:           c
     5:           c   This routine returns a uniform pseudorandom double precision number in the
     6:           c   range (0, 1) by using the linear congruential generator
     7:           c
     8:           c   x_{k+1} = a x_k  (mod 2^46)
     9:           c
    10:           c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
    11:           c   before repeating.  The argument A is the same as 'a' in the above formula,
    12:           c   and X is the same as x_0.  A and X must be odd double precision integers
    13:           c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
    14:           c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
    15:           c   the new seed x_1, so that subsequent calls to RANDLC using the same
    16:           c   arguments will generate a continuous sequence.
    17:           
    18:                 implicit none
    19:                 double precision x, a
    20:                 integer*8 i246m1, Lx, La
    21:                 double precision d2m46
    22:           
    23:                 parameter(d2m46=0.5d0**46)
    24:           
    25:                 save i246m1
    26:                 data i246m1/X'00003FFFFFFFFFFF'/
    27:           
    28:                 Lx = X
    29:                 La = A
    30:           
    31:                 Lx   = iand(Lx*La,i246m1)
    32:                 randlc = d2m46*dble(Lx)
    33:                 x    = dble(Lx)
    34:                 return
    35:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    71: vec( 102): Partially vectorized loop.
    72: opt(1019): Feedback of scalar value from one loop pass to another.: LX
    72: vec( 121): Unvectorizable dependency.
    74: opt(1512): Loop unrolled.: I


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:50 2021
FILE NAME: randi8.f

PROCEDURE NAME: VRANLC
FORMAT LIST

 LINE   LOOP      STATEMENT

    42:                 SUBROUTINE VRANLC (N, X, A, Y)
    43:                 implicit none
    44:                 integer n, i
    45:                 double precision x, a, y(*)
    46:                 integer*8 i246m1, Lx, La
    47:                 double precision d2m46
    48:           
    49:           c This doesn't work, because the compiler does the calculation in 32
    50:           c bits and overflows. No standard way (without f90 stuff) to specify
    51:           c that the rhs should be done in 64 bit arithmetic.
    52:           c      parameter(i246m1=2**46-1)
    53:           
    54:                 parameter(d2m46=0.5d0**46)
    55:           
    56:                 save i246m1
    57:                 data i246m1/X'00003FFFFFFFFFFF'/
    58:           
    59:           c Note that the v6 compiler on an R8000 does something stupid with
    60:           c the above. Using the following instead (or various other things)
    61:           c makes the calculation run almost 10 times as fast.
    62:           c
    63:           c      save d2m46
    64:           c      data d2m46/0.0d0/
    65:           c      if (d2m46 .eq. 0.0d0) then
    66:           c         d2m46 = 0.5d0**46
    67:           c      endif
    68:           
    69:                 Lx = X
    70:                 La = A
    71: S------>        do i = 1, N
    72: |                  Lx   = iand(Lx*La,i246m1)
    73: |                  y(i) = d2m46*dble(Lx)
    74: S------         end do
    75:                 x    = dble(Lx)
    76:           
    77:                 return
    78:                 end


