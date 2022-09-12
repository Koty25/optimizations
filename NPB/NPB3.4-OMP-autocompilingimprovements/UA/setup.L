NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME : setup.f

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fassume-contiguous                      : disable
        -fcopyin-intent-out                      : enable
        -fcse-after-vectorization                : disable
        -ffast-formatted-io                      : enable
        -ffast-math                              : enable
        -fignore-asynchronous                    : disable
        -fignore-induction-variable-overflow     : enable
        -fignore-volatile                        : enable
        -fivdep                                  : disable
        -floop-collapse                          : enable
        -floop-count                             : 5000
        -floop-fusion                            : enable
        -floop-interchange                       : enable
        -floop-normalize                         : enable
        -floop-split                             : disable
        -floop-strip-mine                        : enable
        -floop-unroll                            : enable
        -floop-unroll-complete                   : 4
        -floop-unroll-complete-nest              : 3
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : enable
        -fmove-nested-loop-invariants-outer      : enable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fnamed-noalias-aggressive               : disable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freorder-logical-expression             : enable
        -freplace-loop-equation                  : enable
        -marray-io                               : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mstack-arrays                           : enable
        -mno-stack-arrays                        : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : enable
        -mvector-loop-count-test                 : disable
        -mvector-low-precise-divide-function     : disable
        -mvector-merge-conditional               : enable
        -mvector-packed                          : disable
        -mvector-power-to-explog                 : disable
        -mvector-power-to-sqrt                   : enable
        -mvector-reduction                       : enable
        -mvector-shortloop-reduction             : disable
        -mvector-sqrt-instruction                : disable
        -mvector-threshold                       : 5
        -mwork-vector-kind=none                  : disable

      Parallelization Options : 
        -fopenmp                                 : enable
        -mparallel                               : disable
        -mparallel-innerloop                     : disable
        -mparallel-omp-routine                   : enable
        -mparallel-outerloop-strip-mine          : disable
        -mparallel-sections                      : disable
        -mparallel-threshold                     : 2000
        -mschedule-dynamic                       : disable
        -mschedule-runtime                       : disable
        -mschedule-static                        : disable
        -mschedule-chunk-size                    : 0
        -pthread                                 : disable

      Inlining Options : 
        -finline-copy-arguments                  : disable
        -finline-functions                       : enable
        -finline-max-depth                       : 2
        -finline-max-function-size               : 50
        -finline-max-times                       : 6

      Code Generation Options : 
        -finstrument-functions                   : disable
        -fpic | -fPIC                            : disable
        -ftrace                                  : disable
        -p                                       : disable

      Debugging Options : 
        -fbounds-check                           : disable
        -fcheck=all                              : disable
        -fcheck=alias                            : disable
        -fcheck=bits                             : disable
        -fcheck=bounds                           : disable
        -fcheck=dangling                         : disable
        -fcheck=do                               : disable
        -fcheck=iovf                             : disable
        -fcheck=pointer                          : disable
        -fcheck=present                          : disable
        -fcheck=recursion                        : disable
        -g                                       : disable
        -minit-stack                             : not specified
        -mmemory-trace                           : disable
        -mmemory-trace-full                      : disable
        -traceback                               : disable

      Language Options : 
        -bss                                     : disable
        -fdefault-integer                        : 4
        -fdefault-double                         : 8
        -fdefault-real                           : 4
        -fextend-source                          : disable
        -ffree-form                              : disable
        -ffixed-form                             : enable
        -ff90-sign                               : disable
        -fmax-continuation-lines                 : 1023
        -frealloc-lhs                            : enable
        -frealloc-lhs-array                      : enable
        -frealloc-lhs-scalar                     : enable
        -save                                    : disable
        -std                                     : Fortran 2008

      Message Options : 
        -Wall                                    : disable
        -Werror                                  : disable
        -Wextension                              : enable
        -Wobsolescent                            : enable
        -Woverflow                               : disable
        -Woverflow-errors                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -pedantic-errors                         : disable
        -w                                       : disable

      List Output Options : 
        -report-file                             : not specified
        -report-append-mode                      : disable
        -report-all                              : enable
        -report-cg                               : enable
        -report-diagnostics                      : enable
        -report-format                           : enable
        -report-inline                           : enable
        -report-option                           : enable
        -report-vector                           : enable

      Directory Options : 
        -fintrinsic-modules-path                 : /opt/nec/ve/nfort/3.2.0/include
        -module                                  : not specified

      Miscellaneous Options : 
        -v                                       : disable


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: CREATE_INITIAL_GRID
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    14: opt(1418): Constant-length loop is expanded.
    19: opt(1418): Constant-length loop is expanded.
    26: opt(1418): Constant-length loop is expanded.
    31: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: CREATE_INITIAL_GRID
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

     2:                 subroutine create_initial_grid
     3:           c------------------------------------------------------------------
     4:            
     5:                 use ua_data
     6:                 implicit none
     7:           
     8:                 integer i
     9:           
    10:                 nelt=1
    11:                 ntot=nelt*lx1*lx1*lx1
    12:                 tree(1)=1
    13:                 mt_to_id(1)=1
    14: *------>        do i=1,7,2
    15: |                 xc(i,1)=0.d0
    16: |                 xc(i+1,1)=1.d0
    17: *------         end do
    18:           
    19: *------>        do i=1,2
    20: |                 yc(i,1)=0.d0
    21: |                 yc(2+i,1)=1.d0
    22: |                 yc(4+i,1)=0.d0
    23: |                 yc(6+i,1)=1.d0
    24: *------         end do
    25:            
    26: *------>        do i=1,4
    27: |                 zc(i,1)=0.d0
    28: |                 zc(4+i,1)=1.d0
    29: *------         end do
    30:            
    31: V------>        do i=1,6
    32: |                 cbc(i,1)=0
    33: V------         end do
    34:           
    35:                 return
    36:           
    37:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: CREATE_INITIAL_GRID
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: CREATE_INITIAL_GRID
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (setup.f:31)
    <Vectorized loop.>
    **  Runtime-tested short loop. (setup.f:31)
    *** The number of VGT,   VSC.    :  0,  0. (setup.f:31)
    *** The number of VLOAD, VSTORE. :  0,  1. (setup.f:31)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: CREATE_INITIAL_GRID
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 17 [s8-s11 s51-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: COEF
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    76: vec( 101): Vectorized loop.
    77: opt(1418): Constant-length loop is expanded.
    78: opt(1418): Constant-length loop is expanded.
   101: opt(1418): Constant-length loop is expanded.
   102: opt(1418): Constant-length loop is expanded.
   106: opt(1418): Constant-length loop is expanded.
   107: opt(1418): Constant-length loop is expanded.
   130: opt(1418): Constant-length loop is expanded.
   131: opt(1418): Constant-length loop is expanded.
   163: opt(1418): Constant-length loop is expanded.
   164: opt(1418): Constant-length loop is expanded.
   169: opt(1418): Constant-length loop is expanded.
   170: opt(1418): Constant-length loop is expanded.
   175: opt(1418): Constant-length loop is expanded.
   176: opt(1418): Constant-length loop is expanded.
   189: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: COEF
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

    40:                 subroutine coef
    41:           c-----------------------------------------------------------------
    42:           c
    43:           c     generate
    44:           c
    45:           c            - collocation points
    46:           c            - weights
    47:           c            - derivative matrices
    48:           c            - projection matrices
    49:           c            - interpolation matrices
    50:           c
    51:           c     associated with the
    52:           c
    53:           c            - gauss-legendre lobatto mesh (suffix m1)
    54:           c
    55:           c----------------------------------------------------------------
    56:           
    57:                 use ua_data
    58:                 implicit none
    59:           
    60:                 integer i,j,k
    61:           
    62:           c.....for gauss-legendre lobatto mesh (suffix m1)
    63:           c.....generate collocation points and weights
    64:           
    65:                 zgm1(1)=-1.d0
    66:                 zgm1(2)=-0.6546536707079771d0
    67:                 zgm1(3)=0.d0
    68:                 zgm1(4)= 0.6546536707079771d0
    69:                 zgm1(5)=1.d0
    70:                 wxm1(1)=0.1d0
    71:                 wxm1(2)=49.d0/90.d0
    72:                 wxm1(3)=32.d0/45.d0
    73:                 wxm1(4)=wxm1(2)
    74:                 wxm1(5)=0.1d0
    75:           
    76: V------>        do k=1,lx1
    77: |*----->          do j=1,lx1
    78: ||*---->            do i=1,lx1
    79: |||                   w3m1(i,j,k)=wxm1(i)*wxm1(j)*wxm1(k)
    80: ||*----             end do
    81: |*-----           end do
    82: V------         end do
    83:           
    84:           c.....generate derivative matrices
    85:           
    86:                 dxm1(1,1)=-5.0d0
    87:                 dxm1(2,1)=-1.240990253030982d0
    88:                 dxm1(3,1)= 0.375d0
    89:                 dxm1(4,1)=-0.2590097469690172d0
    90:                 dxm1(5,1)= 0.5d0
    91:                 dxm1(1,2)= 6.756502488724238d0
    92:                 dxm1(2,2)= 0.d0
    93:                 dxm1(3,2)=-1.336584577695453d0
    94:                 dxm1(4,2)= 0.7637626158259734d0
    95:                 dxm1(5,2)=-1.410164177942427d0
    96:                 dxm1(1,3)=-2.666666666666667d0
    97:                 dxm1(2,3)= 1.745743121887939d0
    98:                 dxm1(3,3)= 0.d0
    99:                 dxm1(4,3)=-dxm1(2,3)
   100:                 dxm1(5,3)=-dxm1(1,3)
   101: *------>        do j=4,lx1
   102: |*----->          do i=1,lx1
   103: ||                  dxm1(i,j)=-dxm1(lx1+1-i,lx1+1-j)
   104: |*-----           end do
   105: *------         end do
   106: *------>        do j=1,lx1
   107: |*----->          do i=1,lx1
   108: ||                  dxtm1(i,j)=dxm1(j,i)
   109: |*-----           end do
   110: *------         end do
   111:           
   112:           c.....generate projection (mapping) matrices
   113:           
   114:                 qbnew(1,1,1)=-0.1772843218615690d0
   115:                 qbnew(2,1,1)=9.375d-02
   116:                 qbnew(3,1,1)=-3.700139242414530d-02
   117:                 qbnew(1,2,1)= 0.7152146412463197d0
   118:                 qbnew(2,2,1)=-0.2285757930375471d0
   119:                 qbnew(3,2,1)= 8.333333333333333d-02
   120:                 qbnew(1,3,1)= 0.4398680650316104d0
   121:                 qbnew(2,3,1)= 0.2083333333333333d0
   122:                 qbnew(3,3,1)=-5.891568407922938d-02
   123:                 qbnew(1,4,1)= 8.333333333333333d-02
   124:                 qbnew(2,4,1)= 0.3561799597042137d0
   125:                 qbnew(3,4,1)=-4.854797457965334d-02
   126:                 qbnew(1,5,1)= 0.d0
   127:                 qbnew(2,5,1)=7.03125d-02
   128:                 qbnew(3,5,1)=0.d0
   129:            
   130: *------>        do j=1,lx1
   131: |*----->          do i=1,3
   132: ||                  qbnew(i,j,2)=qbnew(4-i,lx1+1-j,1)
   133: |*-----           end do
   134: *------         end do
   135:           
   136:           c.....generate interpolation matrices for mesh refinement
   137:           
   138:                 ixtmc1(1,1)=1.d0
   139:                 ixtmc1(2,1)=0.d0
   140:                 ixtmc1(3,1)=0.d0
   141:                 ixtmc1(4,1)=0.d0
   142:                 ixtmc1(5,1)=0.d0
   143:                 ixtmc1(1,2)= 0.3385078435248143d0
   144:                 ixtmc1(2,2)= 0.7898516348912331d0
   145:                 ixtmc1(3,2)=-0.1884018684471238d0
   146:                 ixtmc1(4,2)= 9.202967302175333d-02
   147:                 ixtmc1(5,2)=-3.198728299067715d-02
   148:                 ixtmc1(1,3)=-0.1171875d0
   149:                 ixtmc1(2,3)= 0.8840317166357952d0
   150:                 ixtmc1(3,3)= 0.3125d0
   151:                 ixtmc1(4,3)=-0.118406716635795d0
   152:                 ixtmc1(5,3)= 0.0390625d0
   153:                 ixtmc1(1,4)=-7.065070066767144d-02
   154:                 ixtmc1(2,4)= 0.2829703269782467d0
   155:                 ixtmc1(3,4)= 0.902687582732838d0
   156:                 ixtmc1(4,4)=-0.1648516348912333d0
   157:                 ixtmc1(5,4)= 4.984442584781999d-02
   158:                 ixtmc1(1,5)=0.d0
   159:                 ixtmc1(2,5)=0.d0
   160:                 ixtmc1(3,5)=1.d0
   161:                 ixtmc1(4,5)=0.d0
   162:                 ixtmc1(5,5)=0.d0
   163: *------>        do j=1,lx1
   164: |*----->          do i=1,lx1
   165: ||                  ixmc1(i,j)=ixtmc1(j,i)
   166: |*-----           end do
   167: *------         end do
   168:           
   169: *------>        do j=1,lx1
   170: |*----->          do i=1,lx1
   171: ||                  ixtmc2(i,j)=ixtmc1(lx1+1-i,lx1+1-j)
   172: |*-----           end do
   173: *------         end do
   174:           
   175: *------>        do j=1,lx1
   176: |*----->          do i=1,lx1
   177: ||                  ixmc2(i,j)=ixtmc2(j,i)
   178: |*-----           end do
   179: *------         end do
   180:           
   181:           c.....solution interpolation matrix for mesh coarsening
   182:           
   183:                 map2(1)=-0.1179652785083428d0
   184:                 map2(2)= 0.5505046330389332d0
   185:                 map2(3)= 0.7024534364259963d0
   186:                 map2(4)=-0.1972224518285866d0
   187:                 map2(5)= 6.222966087199998d-02
   188:           
   189: *------>        do i=1,lx1
   190: |                 map4(i)=map2(lx1+1-i)
   191: *------         end do
   192:           
   193:                 return
   194:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: COEF
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: COEF
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (setup.f:76)
    <Vectorized loop.>
    **  Runtime-tested short loop. (setup.f:76)
    *** The number of VGT,   VSC.    :  0,  0. (setup.f:76)
    *** The number of VLOAD, VSTORE. :  1, 25. (setup.f:76)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:05 2021
FILE NAME: setup.f

PROCEDURE NAME: COEF
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 47 [s0-s11 s29-s63]
      Vector registers         : 16 [v48-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 0 bytes
      Register spill area      : 0 bytes
      Parameter area           : 0 bytes
      Register save area       : 0 bytes
      User data area           : 0 bytes
      Others                   : 0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: GEOM1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   221: opt(1418): Constant-length loop is expanded.
   225: par(1801): Parallel routine generated.: GEOM1$1
   227: par(1803): Parallelized by "do".
   227: vec( 103): Unvectorized loop.
   227: vec( 113): Overhead of loop division is too large.
   232: vec( 103): Unvectorized loop.
   232: vec( 113): Overhead of loop division is too large.
   233: vec( 103): Unvectorized loop.
   233: vec( 113): Overhead of loop division is too large.
   234: vec( 101): Vectorized loop.
   249: par(1801): Parallel routine generated.: GEOM1$2
   250: par(1803): Parallelized by "do".
   252: vec( 101): Vectorized loop.
   253: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: GEOM1
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   197:                 subroutine geom1
   198:           c-------------------------------------------------------------------
   199:           c
   200:           c     routine to generate elemental geometry information on mesh m1,
   201:           c     (gauss-legendre lobatto mesh).
   202:           c
   203:           c         xrm1_s   -   dx/dr, dy/dr, dz/dr
   204:           c         rxm1_s   -   dr/dx, dr/dy, dr/dz
   205:           c         g1m1_s  geometric factors used in preconditioner computation
   206:           c         g4m1_s  g5m1_s  g6m1_s :
   207:           c         geometric factors used in lapacian opertor
   208:           c         jacm1    -   jacobian
   209:           c         bm1      -   mass matrix
   210:           c         xfrac    -   will be used in prepwork for calculating collocation
   211:           c                      coordinates
   212:           c         idel     -   collocation points index on element boundaries
   213:           c------------------------------------------------------------------
   214:           
   215:                 use ua_data
   216:                 implicit none
   217:           
   218:                 double precision temp,temp1,temp2,dtemp
   219:                 integer isize,i,j,k,ntemp,iel
   220:            
   221: *------>        do i=1,lx1
   222: |                 xfrac(i)=zgm1(i)*0.5d0 + 0.5d0
   223: *------         end do
   224:           
   225:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ISIZE,TEMP,TEMP1,TEMP2,
   226:           c$OMP&  K,J,I,dtemp)
   227: P------>        do isize=1,refine_max
   228: |                 temp=2.d0**(-isize-1)
   229: |                 dtemp=1.d0/temp
   230: |                 temp1=temp**3
   231: |                 temp2=temp**2
   232: |+----->          do k=1,lx1
   233: ||+---->            do j=1,lx1
   234: |||V--->              do i=1,lx1
   235: ||||                    xrm1_s(i,j,k,isize)=dtemp
   236: ||||                    jacm1_s(i,j,k,isize)=temp1
   237: ||||                    rxm1_s(i,j,k,isize)=temp2
   238: ||||                    g1m1_s(i,j,k,isize)=w3m1(i,j,k)*temp
   239: ||||                    bm1_s(i,j,k,isize)=w3m1(i,j,k)*temp1
   240: ||||                    g4m1_s(i,j,k,isize)=g1m1_s(i,j,k,isize)/wxm1(i)
   241: ||||                    g5m1_s(i,j,k,isize)=g1m1_s(i,j,k,isize)/wxm1(j)
   242: ||||                    g6m1_s(i,j,k,isize)=g1m1_s(i,j,k,isize)/wxm1(k)
   243: |||V---               end do
   244: ||+----             end do
   245: |+-----           end do
   246: P------         end do
   247:           c$OMP END PARALLEL DO
   248:           
   249:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ntemp,i,j,iel)
   250: P------>        do iel = 1, lelt
   251: |                 ntemp=lx1*lx1*lx1*(iel-1)
   252: |V----->          do j = 1, lx1
   253: ||*---->            do i = 1, lx1
   254: |||                   idel(i,j,1,iel)=ntemp+(i-1)*lx1 + (j-1)*lx1*lx1+lx1
   255: |||                   idel(i,j,2,iel)=ntemp+(i-1)*lx1 + (j-1)*lx1*lx1+1
   256: |||                   idel(i,j,3,iel)=ntemp+(i-1)*1 + (j-1)*lx1*lx1+lx1*(lx1-1)+1
   257: |||                   idel(i,j,4,iel)=ntemp+(i-1)*1 + (j-1)*lx1*lx1+1
   258: |||                   idel(i,j,5,iel)=ntemp+(i-1)*1 + (j-1)*lx1+lx1*lx1*(lx1-1)+1
   259: |||                   idel(i,j,6,iel)=ntemp+(i-1)*1 + (j-1)*lx1+1
   260: ||*----             end do
   261: |V-----           end do
   262: P------         end do
   263:           c$OMP END PARALLEL DO
   264:           
   265:                 return
   266:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: GEOM1
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: GEOM1
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (setup.f:250)
    <Unvectorized loop.>

    LOOP BEGIN: (setup.f:252)
      <Vectorized loop.>
      **  Runtime-tested short loop. (setup.f:252)
      *** The number of VGT,   VSC.    :  0,  0. (setup.f:252)
      *** The number of VLOAD, VSTORE. :  0, 30. (setup.f:252)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: GEOM1
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 57 [v7-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 272 bytes
      Register spill area      :  40 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  40 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (setup.f:250)
    *** Estimated execution cycle                       : 1444
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 5
            Over basic blocks                           : 5
    *** The number of SCALAR REGISTER TRANSFER          : 2
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: SETDEF
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   279: inl(1212): Source for routine not found.: R_INIT
   281: vec( 101): Vectorized loop.
   282: opt(1418): Constant-length loop is expanded.
   283: opt(1418): Constant-length loop is expanded.
   284: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: SETDEF
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   269:                 subroutine setdef
   270:           c------------------------------------------------------------------
   271:           c     compute the discrete laplacian operators
   272:           c------------------------------------------------------------------
   273:           
   274:                 use ua_data
   275:                 implicit none
   276:           
   277:                 integer i,j,ip
   278:            
   279:                 call r_init(wdtdr(1,1),lx1*lx1,0.d0)
   280:           
   281: V------>        do i=1,lx1
   282: |*----->          do j=1,lx1
   283: ||*---->            do ip=1,lx1
   284: |||     F             wdtdr(i,j) = wdtdr(i,j) + wxm1(ip)*dxm1(ip,i)*dxm1(ip,j)
   285: ||*----             end do
   286: |*-----           end do
   287: V------         end do
   288:           
   289:                 return
   290:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: SETDEF
INLINE LIST

  ROOT: SETDEF (setup.f:269)
  -> NOINLINE: R_INIT (setup.f:279)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: SETDEF
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (setup.f:281)
    <Vectorized loop.>
    **  Runtime-tested short loop. (setup.f:281)
    *** The number of VGT,   VSC.    :  0,  0. (setup.f:281)
    *** The number of VLOAD, VSTORE. : 10,  5. (setup.f:281)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: SETDEF
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 48 [s0-s12 s15-s16 s23-s25 s34-s63]
      Vector registers         : 40 [v24-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 41200 bytes
      Register spill area      : 20488 bytes
      Parameter area           :    24 bytes
      Register save area       :   176 bytes
      User data area           : 20504 bytes
      Others                   :     8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: PREPWORK
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   310: par(1801): Parallel routine generated.: PREPWORK$1
   315: par(1803): Parallelized by "do".
   315: vec( 101): Vectorized loop.
   316: vec( 128): Fused multiply-add operation applied.
   323: par(1803): Parallelized by "do".
   324: inl(1212): Source for routine not found.: R_INIT
   325: vec( 103): Unvectorized loop.
   325: vec( 110): Vectorization obstructive procedure reference.: FACEV
   326: inl(1212): Source for routine not found.: FACEV
   334: par(1803): Parallelized by "do".
   334: vec( 101): Vectorized loop.
   340: opt(3008): Reference within a conditional branch moved outside loop.
   340: par(1803): Parallelized by "do".
   340: vec( 103): Unvectorized loop.
   340: vec( 113): Overhead of loop division is too large.
   341: vec( 101): Vectorized loop.
   344: opt(1418): Constant-length loop is expanded.
   345: opt(1418): Constant-length loop is expanded.
   351: opt(1418): Constant-length loop is expanded.
   359: opt(1418): Constant-length loop is expanded.
   366: opt(1418): Constant-length loop is expanded.
   371: opt(1418): Constant-length loop is expanded.
   374: opt(1418): Constant-length loop is expanded.
   382: opt(1418): Constant-length loop is expanded.
   386: opt(1418): Constant-length loop is expanded.
   389: opt(1418): Constant-length loop is expanded.
   395: opt(1418): Constant-length loop is expanded.
   400: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: PREPWORK
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

 LINE   LOOP      STATEMENT

   294:                 subroutine prepwork
   295:           c------------------------------------------------------------------
   296:           c     mesh information preparations: calculate refinement levels of
   297:           c     each element, mask matrix for domain boundary and element
   298:           c     boundaries
   299:           c------------------------------------------------------------------
   300:           
   301:                 use ua_data
   302:                 implicit none
   303:           
   304:                 integer i, j, iel, iface, cb
   305:                 double precision rdlog2
   306:           
   307:                 ntot = nelt*nxyz
   308:                 rdlog2 = 1.d0/dlog(2.d0)
   309:           
   310:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,IEL,IFACE,CB)
   311:           
   312:           c.....calculate the refinement levels of each element
   313:           
   314:           c$OMP DO
   315: Y------>        do iel = 1, nelt
   316: |       F         size_e(iel)=-dlog(xc(2,iel)-xc(1,iel))*rdlog2+1.d-8
   317: Y------         end do
   318:           c$OMP END DO nowait
   319:           
   320:           c.....mask matrix for element boundary
   321:           
   322:           c$OMP DO
   323: P------>        do iel = 1, nelt
   324: |                 call r_init(tmult(1,1,1,iel),nxyz,1.d0)
   325: |+----->          do iface=1,nsides
   326: ||                  call facev(tmult(1,1,1,iel),iface,0.0d0)
   327: |+-----           end do
   328: P------         end do
   329:           c$OMP END DO nowait
   330:           
   331:           c.....masks for domain boundary at mortar
   332:           
   333:           c$OMP DO
   334: Y------>        do iel=1,nmor
   335: |                 tmmor(iel)=1.d0
   336: Y------         end do
   337:           c$OMP END DO
   338:           
   339:           c$OMP DO
   340: P------>        do iel = 1, nelt
   341: |V----->          do iface = 1,nsides
   342: ||                  cb=cbc(iface,iel)
   343: ||                  if(cb.eq.0) then
   344: ||*---->              do j=2,lx1-1
   345: |||*--->                do i=2,lx1-1
   346: ||||    C                tmmor(idmo(i,j,1,1,iface,iel))=0.d0
   347: |||*---                 end do
   348: ||*----               end do
   349: ||        
   350: ||                    j=1
   351: ||*---->              do i = 1, lx1-1
   352: |||     C                tmmor(idmo(i,j,1,1,iface,iel))=0.d0
   353: ||*----               end do
   354: ||        
   355: ||                    if(idmo(lx1,1,1,1,iface,iel).eq.0)then
   356: ||      C               tmmor(idmo(lx1,1,1,2,iface,iel))=0.d0
   357: ||                    else
   358: ||      C               tmmor(idmo(lx1,1,1,1,iface,iel))=0.d0
   359: ||*---->                do i=1,lx1
   360: |||     C                 tmmor(idmo(i,j,1,2,iface,iel))=0.d0
   361: ||*----                 end do
   362: ||                    end if
   363: ||        
   364: ||                    i=lx1
   365: ||                    if(idmo(lx1,2,1,2,iface,iel).eq.0)then
   366: ||*---->                do j=2,lx1-1
   367: |||     C                 tmmor(idmo(i,j,1,1,iface,iel))=0.d0
   368: ||*----                 end do
   369: ||      C               tmmor(idmo(lx1,lx1,2,2,iface,iel))=0.d0
   370: ||                    else
   371: ||*---->                do j=2,lx1
   372: |||     C                 tmmor(idmo(i,j,1,2,iface,iel))=0.d0
   373: ||*----                 end do
   374: ||*---->                do j=1,lx1
   375: |||     C                 tmmor(idmo(i,j,2,2,iface,iel))=0.d0
   376: ||*----                 end do
   377: ||                    end if
   378: ||         
   379: ||                    j=lx1
   380: ||      C             tmmor(idmo(1,lx1,2,1,iface,iel))=0.d0
   381: ||                    if(idmo(2,lx1,2,1,iface,iel).eq.0)then
   382: ||*---->                do i=2,lx1-1
   383: |||     C                 tmmor(idmo(i,j,1,1,iface,iel))=0.d0
   384: ||*----                 end do
   385: ||                    else
   386: ||*---->                do i=2,lx1
   387: |||     C                 tmmor(idmo(i,j,2,1,iface,iel))=0.d0
   388: ||*----                 end do
   389: ||*---->                do i=1,lx1-1
   390: |||     C                 tmmor(idmo(i,j,2,2,iface,iel))=0.d0
   391: ||*----                 end do
   392: ||                    end if
   393: ||        
   394: ||                    i=1
   395: ||*---->              do j=2,lx1-1
   396: |||     C              tmmor(idmo(i,j,1,1,iface,iel))=0.d0
   397: ||*----               end do
   398: ||                    if(idmo(1,lx1,1,1,iface,iel).ne.0)then
   399: ||      C               tmmor(idmo(i,lx1,1,1,iface,iel))=0.d0
   400: ||*---->                do j=1,lx1-1
   401: |||     C                tmmor(idmo(i,j,2,1,iface,iel))=0.d0
   402: ||*----                 end do
   403: ||                    end if
   404: ||        
   405: ||                  endif
   406: |V-----           end do
   407: P------          end do
   408:           c$OMP END DO nowait
   409:            
   410:           c$OMP END PARALLEL
   411:                 return
   412:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: PREPWORK
INLINE LIST

  ROOT: PREPWORK (setup.f:294)
  -> NOINLINE: R_INIT (setup.f:324)
     *** Source for routine not found.
  -> NOINLINE: FACEV (setup.f:326)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: PREPWORK
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (setup.f:315)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (setup.f:315)
    *** The number of VLOAD, VSTORE. :  2,  1. (setup.f:315)
  LOOP END

  LOOP BEGIN: (setup.f:323)
    <Unvectorized loop.>

    LOOP BEGIN: (setup.f:325)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : FACEV (setup.f:325)
    LOOP END
  LOOP END

  LOOP BEGIN: (setup.f:334)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (setup.f:334)
    *** The number of VLOAD, VSTORE. :  0,  1. (setup.f:334)
  LOOP END

  LOOP BEGIN: (setup.f:340)
    <Unvectorized loop.>
    *** Overhead of loop division is too large. (setup.f:340)

    LOOP BEGIN: (setup.f:341)
      <Vectorized loop.>
      **  Runtime-tested short loop. (setup.f:341)
      *** The number of VGT,   VSC.    :  0, 53. (setup.f:341)
      *** The number of VLOAD, VSTORE. : 54,  0. (setup.f:341)
      *** VSC generated (setup.f:346)
      *** VSC generated (setup.f:352)
      *** VSC generated (setup.f:358)
      *** VSC generated (setup.f:360)
      *** VSC generated (setup.f:356)
      *** VSC generated (setup.f:372)
      *** VSC generated (setup.f:375)
      *** VSC generated (setup.f:367)
      *** VSC generated (setup.f:369)
      *** VSC generated (setup.f:380)
      *** VSC generated (setup.f:387)
      *** VSC generated (setup.f:390)
      *** VSC generated (setup.f:383)
      *** VSC generated (setup.f:396)
      *** VSC generated (setup.f:399)
      *** VSC generated (setup.f:401)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:06 2021
FILE NAME: setup.f

PROCEDURE NAME: PREPWORK
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  3 [vm13-vm15]
      VREG directive           :  0

  Routine stack
    Total size                 : 459632 bytes
      Register spill area      : 229696 bytes
      Parameter area           :     24 bytes
      Register save area       :    176 bytes
      User data area           : 229728 bytes
      Others                   :      8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (setup.f:315)
    *** Estimated execution cycle                       : 242
    *** The number of SCALAR REGISTER TRANSFER          : 4
  LOOP END

  LOOP BEGIN: (setup.f:323)
    *** Estimated execution cycle                       : 36
    *** The number of SCALAR REGISTER TRANSFER          : 15

    LOOP BEGIN: (setup.f:325)
      *** Estimated execution cycle                     : 31
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (setup.f:334)
    *** Estimated execution cycle                       : 39
  LOOP END

  LOOP BEGIN: (setup.f:340)
    *** Estimated execution cycle                       : 2851
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 34
            Not enough registers                        : 34
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 40
            Over basic blocks                           : 40
  LOOP END


