NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME : transfer.f

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
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
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
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
        -fmove-loop-invariants-unsafe            : disable
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
        -freplace-loop-equation                  : disable
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
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: INIT_LOCKS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    13: par(1801): Parallel routine generated.: INIT_LOCKS$1
    14: par(1803): Parallelized by "do".
    14: vec( 103): Unvectorized loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: INIT_LOCKS
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

     2:                 subroutine init_locks
     3:           c------------------------------------------------------------------
     4:           c     Initialize locks to be used for atomic updates
     5:           c------------------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:           
    10:                 integer i
    11:           
    12:           c.....initialize locks in parallel
    13:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    14: P------>  c$    do i=1,lmor
    15: |         c$      call omp_init_lock(tlock(i))
    16: P------   c$    end do
    17:           
    18:                 return
    19:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: INIT_LOCKS
VECTORIZATION LIST

  LOOP BEGIN: (transfer.f:14)
    <Unvectorized loop.>
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: INIT_LOCKS
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 17 [s0-s2 s8-s12 s15-s16 s18 s23-s25 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:14)
    *** Estimated execution cycle                       : 10
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    39: par(1801): Parallel routine generated.: TRANSF$1
    41: opt(3008): Reference within a conditional branch moved outside loop.
    41: par(1803): Parallelized by "do".
    42: opt(3008): Reference within a conditional branch moved outside loop.
    42: vec( 110): Vectorization obstructive procedure reference.: R_INIT
    75: opt(1025): Reference to this procedure inhibits optimization.: R_INIT
    78: vec( 103): Unvectorized loop.
    78: vec( 113): Overhead of loop division is too large.
    79: vec( 101): Vectorized loop.
    89: opt(1418): Constant-length loop is expanded.
    91: opt(1418): Constant-length loop is expanded.
    93: opt(1019): Feedback of scalar value from one loop pass to another.: TMP
    93: vec( 128): Fused multiply-add operation applied.
   105: vec( 108): Unvectorizable loop structure.
   113: opt(1418): Constant-length loop is expanded.
   115: opt(1418): Constant-length loop is expanded.
   116: opt(1019): Feedback of scalar value from one loop pass to another.
   122: vec( 103): Unvectorized loop.
   122: vec( 107): Iteration count is too small.
   132: opt(1418): Constant-length loop is expanded.
   134: opt(1418): Constant-length loop is expanded.
   142: opt(1418): Constant-length loop is expanded.
   144: opt(1418): Constant-length loop is expanded.
   145: opt(1019): Feedback of scalar value from one loop pass to another.
   155: opt(1418): Constant-length loop is expanded.
   156: opt(1418): Constant-length loop is expanded.
   168: vec( 103): Unvectorized loop.
   168: vec( 107): Iteration count is too small.
   170: opt(1418): Constant-length loop is expanded.
   171: opt(1418): Constant-length loop is expanded.
   180: opt(1418): Constant-length loop is expanded.
   189: vec( 103): Unvectorized loop.
   189: vec( 107): Iteration count is too small.
   191: opt(1418): Constant-length loop is expanded.
   192: opt(1418): Constant-length loop is expanded.
   201: opt(1418): Constant-length loop is expanded.
   210: vec( 103): Unvectorized loop.
   210: vec( 107): Iteration count is too small.
   212: opt(1418): Constant-length loop is expanded.
   213: opt(1418): Constant-length loop is expanded.
   222: opt(1418): Constant-length loop is expanded.
   231: vec( 103): Unvectorized loop.
   231: vec( 107): Iteration count is too small.
   233: opt(1418): Constant-length loop is expanded.
   234: opt(1418): Constant-length loop is expanded.
   242: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF
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

    23:                 subroutine transf(tmor,tx)
    24:           c------------------------------------------------------------------
    25:           c     Map values from mortar(tmor) to element(tx)
    26:           c------------------------------------------------------------------
    27:           
    28:                 use ua_data
    29:                 implicit none
    30:           
    31:                 double precision tmor(*),tx(*), tmp(lx1,lx1,2)
    32:                 integer ig1,ig2,ig3,ig4,ie,iface,il1,il2,il3,il4,
    33:                &        nnje,ije1,ije2,col,i,j,ig,il
    34:           
    35:           
    36:           c.....zero out tx on element boundaries
    37:                 call col2(tx,tmult,ntot)
    38:           
    39:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(il,j,ig,i,col,ije2,ije1,ig4,
    40:           c$OMP& ig3,ig2,ig1,nnje,il4,il3,il2,il1,iface,ie,tmp)
    41: P------>        do ie=1,nelt
    42: |+----->          do iface=1,nsides
    43: ||        
    44: ||        c.........get the collocation point index of the four local corners on the
    45: ||        c         face iface of element ie
    46: ||                  il1=idel(1,1,iface,ie)
    47: ||                  il2=idel(lx1,1,iface,ie)
    48: ||                  il3=idel(1,lx1,iface,ie)
    49: ||                  il4=idel(lx1,lx1,iface,ie)
    50: ||        
    51: ||        c.........get the mortar indices of the four local corners
    52: ||                  ig1= idmo(1,  1  ,1,1,iface,ie)
    53: ||                  ig2= idmo(lx1,1  ,1,2,iface,ie)
    54: ||                  ig3= idmo(1,  lx1,2,1,iface,ie)
    55: ||                  ig4= idmo(lx1,lx1,2,2,iface,ie)
    56: ||         
    57: ||        c.........copy the value from tmor to tx for these four local corners
    58: ||                  tx(il1) = tmor(ig1)
    59: ||                  tx(il2) = tmor(ig2)
    60: ||                  tx(il3) = tmor(ig3)
    61: ||                  tx(il4) = tmor(ig4)
    62: ||         
    63: ||        c.........nnje=1 for conforming faces, nnje=2 for nonconforming faces
    64: ||                  if(cbc(iface,ie).eq.3) then
    65: ||                    nnje=2
    66: ||                  else
    67: ||                    nnje=1
    68: ||                  end if
    69: ||        
    70: ||        c.........for nonconforming faces
    71: ||                  if(nnje.eq.2) then
    72: ||        
    73: ||        c...........nonconforming faces have four pieces of mortar, first map them to
    74: ||        c           two intermediate mortars, stored in tmp
    75: ||                    call r_init(tmp,lx1*lx1*2,0.d0)
    76: ||         
    77: ||+---->              do ije1=1,nnje
    78: |||+--->                do ije2=1,nnje
    79: ||||V-->                  do col=1,lx1
    80: |||||     
    81: |||||     c.................in each row col, when coloumn i=1 or lx1, the value
    82: |||||     c                 in tmor is copied to tmp
    83: |||||                       i = v_end(ije2)
    84: |||||                       ig=idmo(i,col,ije1,ije2,iface,ie)
    85: |||||   G                   tmp(i,col,ije1)=tmor(ig)
    86: |||||     
    87: |||||     c.................in each row col, value in the interior three collocation
    88: |||||     c                 points is computed by apply mapping matrix qbnew to tmor
    89: |||||*->                    do i=2,lx1-1
    90: ||||||                        il= idel(i,col,iface,ie)
    91: ||||||*>                      do j=1,lx1
    92: |||||||                         ig=idmo(j,col,ije1,ije2,iface,ie)
    93: ||||||| G                       tmp(i,col,ije1) = tmp(i,col,ije1) +
    94: |||||||        &                qbnew(i-1,j,ije2)*tmor(ig)
    95: ||||||*                       end do
    96: |||||*-                     end do
    97: |||||     
    98: ||||V--                   end do
    99: |||+---                 end do
   100: ||+----               end do
   101: ||         
   102: ||        c...........mapping from two pieces of intermediate mortar tmp to element
   103: ||        c           face tx
   104: ||        
   105: ||+---->              do ije1=1, nnje
   106: |||       
   107: |||       c.............the first column, col=1, is an edge of face iface.
   108: |||       c             the value on the three interior collocation points, tx, is
   109: |||       c             computed by applying mapping matrices qbnew to tmp.
   110: |||       c             the mapping result is divided by 2, because there will be
   111: |||       c             duplicated contribution from another face sharing this edge.
   112: |||                     col=1
   113: |||*--->                do i=2,lx1-1
   114: ||||                      il= idel(col,i,iface,ie)
   115: ||||*-->                  do j=1,lx1
   116: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*
   117: |||||          &                       tmp(col,j,ije1)*0.5d0
   118: ||||*--                   end do
   119: |||*---                 end do
   120: |||       
   121: |||       c.............for column 2 ~ lx-1
   122: |||+--->                do col=2,lx1-1
   123: ||||      
   124: ||||      c...............when i=1 or lx1, the collocation points are also on an edge of
   125: ||||      c               the face, so the mapping result also needs to be divided by 2
   126: ||||                      i = v_end(ije1)
   127: ||||                      il= idel(col,i,iface,ie)
   128: ||||                      tx(il)=tx(il)+tmp(col,i,ije1)*0.5d0
   129: ||||      
   130: ||||      c...............compute the value at interior collocation points in
   131: ||||      c               columns 2 ~ lx1
   132: ||||*-->                  do i=2,lx1-1
   133: |||||                       il= idel(col,i,iface,ie)
   134: |||||*->                    do j=1,lx1
   135: ||||||                        tx(il) = tx(il) + qbnew(i-1,j,ije1)* tmp(col,j,ije1)
   136: |||||*-                     end do
   137: ||||*--                   end do
   138: |||+---                 end do
   139: |||       
   140: |||       c.............same as col=1
   141: |||                     col=lx1
   142: |||*--->                do  i=2,lx1-1
   143: ||||                      il= idel(col,i,iface,ie)
   144: ||||*-->                  do j=1,lx1
   145: |||||                       tx(il) = tx(il) + qbnew(i-1,j,ije1)*
   146: |||||          &                     tmp(col,j,ije1)*0.5d0
   147: ||||*--                   end do
   148: |||*---                 end do
   149: ||+----               end do
   150: ||        
   151: ||        c.........for conforming faces
   152: ||                  else
   153: ||        
   154: ||        c.........face interior
   155: ||*---->              do col=2,lx1-1
   156: |||*--->                do i=2,lx1-1
   157: ||||                      il= idel(i,col,iface,ie)
   158: ||||                      ig= idmo(i,col,1,1,iface,ie)
   159: ||||                      tx(il)=tmor(ig)
   160: |||*---                 end do
   161: ||*----               end do
   162: ||        
   163: ||         
   164: ||        c...........edges of conforming faces
   165: ||        
   166: ||        c...........if local edge 1 is a nonconforming edge
   167: ||                    if(idmo(lx1,1,1,1,iface,ie).ne.0)then
   168: ||+---->                do i=2,lx1-1
   169: |||                       il= idel(i,1,iface,ie)
   170: |||*--->                  do ije1=1,2
   171: ||||*-->                    do j=1,lx1
   172: |||||                         ig=idmo(j,1,1,ije1,iface,ie)
   173: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   174: ||||*--                     end do
   175: |||*---                   end do
   176: ||+----                 end do
   177: ||        
   178: ||        c...........if local edge 1 is a conforming edge
   179: ||                    else
   180: ||*---->                do i=2,lx1-1
   181: |||                       il= idel(i,1,iface,ie)
   182: |||                       ig= idmo(i,1,1,1,iface,ie)
   183: |||                       tx(il)=tmor(ig)
   184: ||*----                 end do
   185: ||                    end if
   186: ||        
   187: ||        c...........if local edge 2 is a nonconforming edge
   188: ||                    if(idmo(lx1,2,1,2,iface,ie).ne.0)then
   189: ||+---->                do i=2,lx1-1
   190: |||                       il= idel(lx1,i,iface,ie)
   191: |||*--->                  do ije1=1,2
   192: ||||*-->                    do j=1,lx1
   193: |||||                         ig=idmo(lx1,j,ije1,2,iface,ie)
   194: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   195: ||||*--                     end do
   196: |||*---                   end do
   197: ||+----                 end do
   198: ||        
   199: ||        c...........if local edge 2 is a conforming edge
   200: ||                    else
   201: ||*---->                do i=2,lx1-1
   202: |||                       il= idel(lx1,i,iface,ie)
   203: |||                       ig= idmo(lx1,i,1,1,iface,ie)
   204: |||                       tx(il)=tmor(ig)
   205: ||*----                 end do
   206: ||                    end if
   207: ||        
   208: ||        c...........if local edge 3 is a nonconforming edge
   209: ||                    if(idmo(2,lx1,2,1,iface,ie).ne.0)then
   210: ||+---->                do  i=2,lx1-1
   211: |||                       il= idel(i,lx1,iface,ie)
   212: |||*--->                  do ije1=1,2
   213: ||||*-->                    do j=1,lx1
   214: |||||                         ig=idmo(j,lx1,2,ije1,iface,ie)
   215: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   216: ||||*--                     end do
   217: |||*---                   end do
   218: ||+----                 end do
   219: ||        
   220: ||        c...........if local edge 3 is a conforming edge
   221: ||                    else
   222: ||*---->                do i=2,lx1-1
   223: |||                       il= idel(i,lx1,iface,ie)
   224: |||                       ig= idmo(i,lx1,1,1,iface,ie)
   225: |||                       tx(il)=tmor(ig)
   226: ||*----                 end do
   227: ||                    end if
   228: ||        
   229: ||        c...........if local edge 4 is a nonconforming edge
   230: ||                    if(idmo(1,lx1,1,1,iface,ie).ne.0)then
   231: ||+---->                do i=2,lx1-1
   232: |||                       il= idel(1,i,iface,ie)
   233: |||*--->                  do ije1=1,2
   234: ||||*-->                    do j=1,lx1
   235: |||||                         ig=idmo(1,j,ije1,1,iface,ie)
   236: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   237: ||||*--                     end do
   238: |||*---                   end do
   239: ||+----                 end do
   240: ||        c...........if local edge 4 is a conforming edge
   241: ||                    else
   242: ||*---->                do i=2,lx1-1
   243: |||                       il= idel(1,i,iface,ie)
   244: |||                       ig= idmo(1,i,1,1,iface,ie)
   245: |||                       tx(il)=tmor(ig)
   246: ||*----                 end do
   247: ||                    end if
   248: ||                  end if
   249: ||         
   250: |+-----           end do
   251: P------         end do
   252:           c$OMP END PARALLEL DO
   253:           
   254:                 return
   255:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:41)
    <Unvectorized loop.>

    LOOP BEGIN: (transfer.f:42)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : R_INIT (transfer.f:42)

      LOOP BEGIN: (transfer.f:168)
        <Unvectorized loop.>
        *** Iteration count is too small. (transfer.f:168)
      LOOP END

      LOOP BEGIN: (transfer.f:189)
        <Unvectorized loop.>
        *** Iteration count is too small. (transfer.f:189)
      LOOP END

      LOOP BEGIN: (transfer.f:210)
        <Unvectorized loop.>
        *** Iteration count is too small. (transfer.f:210)
      LOOP END

      LOOP BEGIN: (transfer.f:231)
        <Unvectorized loop.>
        *** Iteration count is too small. (transfer.f:231)
      LOOP END

      LOOP BEGIN: (transfer.f:77)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:78)
          <Unvectorized loop.>
          *** Overhead of loop division is too large. (transfer.f:78)

          LOOP BEGIN: (transfer.f:79)
            <Vectorized loop.>
            **  Runtime-tested short loop. (transfer.f:79)
            *** The number of VGT,   VSC.    :  6,  0. (transfer.f:79)
            *** The number of VLOAD, VSTORE. :  9,  4. (transfer.f:79)
            *** VGT generated (transfer.f:85)
            *** VGT generated (transfer.f:93)
          LOOP END
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:105)
        <Unvectorized loop.>
        *** Unvectorizable loop structure. (transfer.f:105)

        LOOP BEGIN: (transfer.f:122)
          <Unvectorized loop.>
          *** Iteration count is too small. (transfer.f:122)
        LOOP END
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 36 [v28-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1424 bytes
      Register spill area      :  384 bytes
      Parameter area           :   24 bytes
      Register save area       :  176 bytes
      User data area           :  832 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:41)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 4
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER TRANSFER          : 2

    LOOP BEGIN: (transfer.f:42)
      *** Estimated execution cycle                     : 144
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 25
              Over basic blocks                         : 25
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 27
              Over basic blocks                         : 27
      *** The number of SCALAR REGISTER TRANSFER        : 62

      LOOP BEGIN: (transfer.f:168)
        *** Estimated execution cycle                   : 78
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END

      LOOP BEGIN: (transfer.f:189)
        *** Estimated execution cycle                   : 78
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 3
                Over basic blocks                       : 3
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 4
                Over basic blocks                       : 4
        *** The number of SCALAR REGISTER TRANSFER      : 11
      LOOP END

      LOOP BEGIN: (transfer.f:210)
        *** Estimated execution cycle                   : 78
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 4
                Over basic blocks                       : 4
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 3
                Over basic blocks                       : 3
        *** The number of SCALAR REGISTER TRANSFER      : 21
      LOOP END

      LOOP BEGIN: (transfer.f:231)
        *** Estimated execution cycle                   : 78
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 2
                Over basic blocks                       : 2
        *** The number of SCALAR REGISTER TRANSFER      : 10
      LOOP END

      LOOP BEGIN: (transfer.f:77)
        *** Estimated execution cycle                   : 30
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 5
                Over basic blocks                       : 5
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 9
                Over basic blocks                       : 9
        *** The number of SCALAR REGISTER TRANSFER      : 9

        LOOP BEGIN: (transfer.f:78)
          *** Estimated execution cycle                 : 1319
          *** The number of SCALAR REGISTER TRANSFER    : 4
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:105)
        *** Estimated execution cycle                   : 126
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 15
                Over basic blocks                       : 15
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 21
                Over basic blocks                       : 21
        *** The number of SCALAR REGISTER TRANSFER      : 6

        LOOP BEGIN: (transfer.f:122)
          *** Estimated execution cycle                 : 80
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 6
                  Over basic blocks                     : 6
        LOOP END
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   277: par(1801): Parallel routine generated.: TRANSFB$1
   282: par(1803): Parallelized by "do".
   282: vec( 101): Vectorized loop.
   288: opt(3008): Reference within a conditional branch moved outside loop.
   288: par(1803): Parallelized by "do".
   289: opt(3008): Reference within a conditional branch moved outside loop.
   336: vec( 101): Vectorized loop.
   336: vec( 124): Iteration count is assumed. Iteration count.: 5
   338: vec( 101): Vectorized loop.
   348: opt(1418): Constant-length loop is expanded.
   350: vec( 128): Fused multiply-add operation applied.
   357: opt(1418): Constant-length loop is expanded.
   359: opt(1418): Constant-length loop is expanded.
   361: vec( 128): Fused multiply-add operation applied.
   391: vec( 103): Unvectorized loop.
   393: opt(1418): Constant-length loop is expanded.
   398: opt(1097): This statement prevents loop optimization.
   415: vec( 103): Unvectorized loop.
   418: opt(1418): Constant-length loop is expanded.
   424: opt(1097): This statement prevents loop optimization.
   435: vec( 103): Unvectorized loop.
   436: opt(1418): Constant-length loop is expanded.
   440: opt(1097): This statement prevents loop optimization.
   451: vec( 103): Unvectorized loop.
   453: opt(1418): Constant-length loop is expanded.
   459: opt(1097): This statement prevents loop optimization.
   467: opt(1418): Constant-length loop is expanded.
   480: vec( 103): Unvectorized loop.
   482: opt(1418): Constant-length loop is expanded.
   488: opt(1097): This statement prevents loop optimization.
   496: opt(1418): Constant-length loop is expanded.
   509: vec( 103): Unvectorized loop.
   511: opt(1418): Constant-length loop is expanded.
   517: opt(1097): This statement prevents loop optimization.
   525: opt(1418): Constant-length loop is expanded.
   538: vec( 103): Unvectorized loop.
   540: opt(1418): Constant-length loop is expanded.
   546: opt(1097): This statement prevents loop optimization.
   554: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB
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

   259:                 subroutine transfb(tmor,tx)
   260:           c------------------------------------------------------------------
   261:           c     Map from element(tx) to mortar(tmor).
   262:           c     tmor sums contributions from all elements.
   263:           c------------------------------------------------------------------
   264:           
   265:                 use ua_data
   266:                 implicit none
   267:           
   268:                 double precision third
   269:                 parameter (third=1.d0/3.d0)
   270:                 integer shift
   271:           
   272:                 double precision tmp,tmp1,tx(*),tmor(*),temp(lx1,lx1,2),
   273:                &                 top(lx1,2)
   274:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,nnje,
   275:                &        ije1,ije2,col,i,j,ije,ig,il
   276:           
   277:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(il,j,ig,i,col,ije2,ije1,ig4,
   278:           c$OMP& ig3,ig2,ig1,nnje,il4,il3,il2,il1,iface,ie,ije,
   279:           c$OMP& tmp,shift,temp,top,tmp1)
   280:           
   281:           c$OMP DO
   282: Y------>        do ie=1,nmor
   283: |                 tmor(ie)=0.d0
   284: Y------         end do
   285:           c$OMP END DO
   286:           
   287:           c$OMP DO
   288: P------>        do ie=1,nelt
   289: |+----->          do iface=1,nsides
   290: ||        c.........nnje=1 for conforming faces, nnje=2 for nonconforming faces
   291: ||                  if(cbc(iface,ie).eq.3) then
   292: ||                    nnje=2
   293: ||                  else
   294: ||                    nnje=1
   295: ||                  end if
   296: ||        
   297: ||        c.........get collocation point index of four local corners on the face
   298: ||                  il1 = idel(1,  1,  iface,ie)
   299: ||                  il2 = idel(lx1,1,  iface,ie)
   300: ||                  il3 = idel(1,  lx1,iface,ie)
   301: ||                  il4 = idel(lx1,lx1,iface,ie)
   302: ||        
   303: ||        c.........get the mortar indices of the four local corners
   304: ||                  ig1 = idmo(1,  1,  1,1,iface,ie)
   305: ||                  ig2 = idmo(lx1,1,  1,2,iface,ie)
   306: ||                  ig3 = idmo(1,  lx1,2,1,iface,ie )
   307: ||                  ig4 = idmo(lx1,lx1,2,2,iface,ie)
   308: ||        
   309: ||        c.........sum the values from tx to tmor for these four local corners
   310: ||        c         only 1/3 of the value is summed, since there will be two duplicated
   311: ||        c         contributions from the other two faces sharing this vertex
   312: ||        c
   313: ||        c$        call omp_set_lock(tlock(ig1))
   314: ||                  tmor(ig1) = tmor(ig1)+tx(il1)*third
   315: ||        c$        call omp_unset_lock(tlock(ig1))
   316: ||        c
   317: ||        c$        call omp_set_lock(tlock(ig2))
   318: ||                  tmor(ig2) = tmor(ig2)+tx(il2)*third
   319: ||        c$        call omp_unset_lock(tlock(ig2))
   320: ||        c
   321: ||        c$        call omp_set_lock(tlock(ig3))
   322: ||                  tmor(ig3) = tmor(ig3)+tx(il3)*third
   323: ||        c$        call omp_unset_lock(tlock(ig3))
   324: ||        c
   325: ||        c$        call omp_set_lock(tlock(ig4))
   326: ||                  tmor(ig4) = tmor(ig4)+tx(il4)*third
   327: ||        c$        call omp_unset_lock(tlock(ig4))
   328: ||        
   329: ||        c.........for nonconforming faces
   330: ||                  if(nnje.eq.2) then
   331: ||                    call r_init(temp,lx1*lx1*2,0.d0)
   332: ||        
   333: ||        c...........nonconforming faces have four pieces of mortar, first map tx to
   334: ||        c           two intermediate mortars stored in temp
   335: ||        
   336: ||V---->              do ije2 = 1, nnje
   337: |||                     shift = ije2-1
   338: |||V--->                do col=1,lx1
   339: ||||      c...............For mortar points on face edge (top and bottom), copy the
   340: ||||      c               value from tx to temp
   341: ||||                      il=idel(col,v_end(ije2),iface,ie)
   342: ||||    G                 temp(col,v_end(ije2),ije2)=tx(il)
   343: ||||      
   344: ||||      c...............For mortar points on face edge (top and bottom), calculate
   345: ||||      c               the interior points' contribution to them, i.e. top()
   346: ||||                      j = v_end(ije2)
   347: ||||                      tmp=0.d0
   348: ||||*-->                  do i=2,lx1-1
   349: |||||                       il=idel(col,i,iface,ie)
   350: |||||   G                   tmp = tmp + qbnew(i-1,j,ije2)*tx(il)
   351: ||||*--                   end do
   352: ||||      
   353: ||||                      top(col,ije2)=tmp
   354: ||||      
   355: ||||      c...............Use mapping matrices qbnew to map the value from tx to temp
   356: ||||      c               for mortar points not on the top bottom face edge.
   357: ||||*-->                  do j=2-shift,lx1-shift
   358: |||||                       tmp=0.d0
   359: |||||*->                    do i=2,lx1-1
   360: ||||||                        il=idel(col,i,iface,ie)
   361: ||||||  F                     tmp = tmp + qbnew(i-1,j,ije2)*tx(il)
   362: |||||*-                     end do
   363: |||||                       temp(col,j,ije2) = tmp + temp(col,j,ije2)
   364: ||||*--                   end do
   365: |||V---                 end do
   366: ||V----               end do
   367: ||        
   368: ||        c...........mapping from temp to tmor
   369: ||        
   370: ||+---->              do ije1=1, nnje
   371: |||                     shift = ije1-1
   372: |||+--->                do ije2=1,nnje
   373: ||||      
   374: ||||      c...............for each column of collocation points on a piece of mortar
   375: ||||+-->                  do col=2-shift,lx1-shift
   376: |||||     
   377: |||||     c.................For the end point, which is on an edge (local edge 2,4),
   378: |||||     c                 the contribution is halved since there will be duplicated
   379: |||||     c                 contribution from another face sharing this edge.
   380: |||||     
   381: |||||                       ig=idmo(v_end(ije2),col,ije1,ije2,iface,ie)
   382: |||||     c
   383: |||||     c$                call omp_set_lock(tlock(ig))
   384: |||||                       tmor(ig)=tmor(ig)+temp(v_end(ije2),col,ije1)*0.5d0
   385: |||||     c$                call omp_unset_lock(tlock(ig))
   386: |||||     
   387: |||||     c.................In each row of collocation points on a piece of mortar,
   388: |||||     c                 sum the contributions from interior collocation points
   389: |||||     c                 (i=2,lx1-1)
   390: |||||     
   391: |||||+->                    do  j=1,lx1
   392: ||||||                        tmp=0.d0
   393: ||||||*>                      do i=2,lx1-1
   394: |||||||                         tmp = tmp + qbnew(i-1,j,ije2) * temp(i,col,ije1)
   395: ||||||*                       end do
   396: ||||||                        ig=idmo(j,col,ije1,ije2,iface,ie)
   397: ||||||    c
   398: ||||||    c$                  call omp_set_lock(tlock(ig))
   399: ||||||                        tmor(ig)=tmor(ig)+tmp
   400: ||||||    c$                  call omp_unset_lock(tlock(ig))
   401: |||||+-                     end do
   402: ||||+--                   end do
   403: ||||      
   404: ||||      c...............For tmor on local edge 1 and 3, tmp is the contribution from
   405: ||||      c               an edge, so it is halved because of duplicated contribution
   406: ||||      c               from another face sharing this edge. tmp1 is contribution
   407: ||||      c               from face interior.
   408: ||||      
   409: ||||                      col = v_end(ije1)
   410: ||||                      ig=idmo(v_end(ije2),col,ije1,ije2,iface,ie)
   411: ||||      c
   412: ||||      c$              call omp_set_lock(tlock(ig))
   413: ||||                      tmor(ig)=tmor(ig)+top(v_end(ije2),ije1)*0.5d0
   414: ||||      c$              call omp_unset_lock(tlock(ig))
   415: ||||+-->                  do  j=1,lx1
   416: |||||                       tmp=0.d0
   417: |||||                       tmp1=0.d0
   418: |||||*->                    do i=2,lx1-1
   419: ||||||                        tmp  = tmp  + qbnew(i-1,j,ije2) * temp(i,col,ije1)
   420: ||||||                        tmp1 = tmp1 + qbnew(i-1,j,ije2) * top(i,ije1)
   421: |||||*-                     end do
   422: |||||                       ig=idmo(j,col,ije1,ije2,iface,ie)
   423: |||||     c
   424: |||||     c$                call omp_set_lock(tlock(ig))
   425: |||||                       tmor(ig)=tmor(ig)+tmp*0.5d0+tmp1
   426: |||||     c$                call omp_unset_lock(tlock(ig))
   427: ||||+--                   end do
   428: |||+---                 end do
   429: ||+----               end do
   430: ||        
   431: ||        c.........for conforming faces
   432: ||                  else
   433: ||        
   434: ||        c.........face interior
   435: ||+---->              do col=2,lx1-1
   436: |||*--->                do j=2,lx1-1
   437: ||||                      il=idel(j,col,iface,ie)
   438: ||||                      ig=idmo(j,col,1,1,iface,ie)
   439: ||||      c
   440: ||||      c$              call omp_set_lock(tlock(ig))
   441: ||||                      tmor(ig)=tmor(ig)+tx(il)
   442: ||||      c$              call omp_unset_lock(tlock(ig))
   443: |||*---                 end do
   444: ||+----               end do
   445: ||        
   446: ||        c...........edges of conforming faces
   447: ||        
   448: ||        c...........if local edge 1 is a nonconforming edge
   449: ||                    if(idmo(lx1,1,1,1,iface,ie).ne.0)then
   450: ||+---->                do ije=1,2
   451: |||+--->                  do j=1,lx1
   452: ||||                        tmp=0.d0
   453: ||||*-->                    do i=2,lx1-1
   454: |||||                         il=idel(i,1,iface,ie)
   455: |||||                         tmp= tmp + qbnew(i-1,j,ije)*tx(il)
   456: ||||*--                     end do
   457: ||||                        ig=idmo(j,1,1,ije,iface,ie)
   458: ||||      c
   459: ||||      c$                call omp_set_lock(tlock(ig))
   460: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   461: ||||      c$                call omp_unset_lock(tlock(ig))
   462: |||+---                   end do
   463: ||+----                 end do
   464: ||        
   465: ||        c...........if local edge 1 is a conforming edge
   466: ||                    else
   467: ||*---->                do j=2,lx1-1
   468: |||                       il=idel(j,1,iface,ie)
   469: |||                       ig=idmo(j,1,1,1,iface,ie)
   470: |||       c
   471: |||       c$              call omp_set_lock(tlock(ig))
   472: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   473: |||       c$              call omp_unset_lock(tlock(ig))
   474: ||*----                 end do
   475: ||                    end if
   476: ||        
   477: ||        c...........if local edge 2 is a nonconforming edge
   478: ||                    if(idmo(lx1,2,1,2,iface,ie).ne.0)then
   479: ||+---->                do ije=1,2
   480: |||+--->                  do j=1,lx1
   481: ||||                        tmp=0.d0
   482: ||||*-->                    do i=2,lx1-1
   483: |||||                         il=idel(lx1,i,iface,ie)
   484: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   485: ||||*--                     end do
   486: ||||                        ig=idmo(lx1,j,ije,2,iface,ie)
   487: ||||      c
   488: ||||      c$                call omp_set_lock(tlock(ig))
   489: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   490: ||||      c$                call omp_unset_lock(tlock(ig))
   491: |||+---                   end do
   492: ||+----                 end do
   493: ||        
   494: ||        c...........if local edge 2 is a conforming edge
   495: ||                    else
   496: ||*---->                do j=2,lx1-1
   497: |||                       il=idel(lx1,j,iface,ie)
   498: |||                       ig=idmo(lx1,j,1,1,iface,ie)
   499: |||       c
   500: |||       c$              call omp_set_lock(tlock(ig))
   501: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   502: |||       c$              call omp_unset_lock(tlock(ig))
   503: ||*----                 end do
   504: ||                    end if
   505: ||        
   506: ||        c...........if local edge 3 is a nonconforming edge
   507: ||                    if(idmo(2,lx1,2,1,iface,ie).ne.0)then
   508: ||+---->                do ije=1,2
   509: |||+--->                  do j=1,lx1
   510: ||||                        tmp=0.d0
   511: ||||*-->                    do i=2,lx1-1
   512: |||||                         il=idel(i,lx1,iface,ie)
   513: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   514: ||||*--                     end do
   515: ||||                        ig=idmo(j,lx1,2,ije,iface,ie)
   516: ||||      c
   517: ||||      c$                call omp_set_lock(tlock(ig))
   518: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   519: ||||      c$                call omp_unset_lock(tlock(ig))
   520: |||+---                   end do
   521: ||+----                 end do
   522: ||        
   523: ||        c...........if local edge 3 is a conforming edge
   524: ||                    else
   525: ||*---->                do j=2,lx1-1
   526: |||                       il=idel(j,lx1,iface,ie)
   527: |||                       ig=idmo(j,lx1,1,1,iface,ie)
   528: |||       c
   529: |||       c$              call omp_set_lock(tlock(ig))
   530: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   531: |||       c$              call omp_unset_lock(tlock(ig))
   532: ||*----                 end do
   533: ||                    end if
   534: ||        
   535: ||        c...........if local edge 4 is a nonconforming edge
   536: ||                    if(idmo(1,lx1,1,1,iface,ie).ne.0)then
   537: ||+---->                do ije=1,2
   538: |||+--->                  do j=1,lx1
   539: ||||                        tmp=0.d0
   540: ||||*-->                    do i=2,lx1-1
   541: |||||                         il=idel(1,i,iface,ie)
   542: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   543: ||||*--                     end do
   544: ||||                        ig=idmo(1,j,ije,1,iface,ie)
   545: ||||      c
   546: ||||      c$                call omp_set_lock(tlock(ig))
   547: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   548: ||||      c$                call omp_unset_lock(tlock(ig))
   549: |||+---                   end do
   550: ||+----                 end do
   551: ||        
   552: ||        c...........if local edge 4 is a conforming edge
   553: ||                    else
   554: ||*---->                do j=2,lx1-1
   555: |||                       il=idel(1,j,iface,ie)
   556: |||                       ig=idmo(1,j,1,1,iface,ie)
   557: |||       c
   558: |||       c$              call omp_set_lock(tlock(ig))
   559: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   560: |||       c$              call omp_unset_lock(tlock(ig))
   561: ||*----                 end do
   562: ||                    end if
   563: ||                  end if
   564: |+-----           end do
   565: P------         end do
   566:           c$OMP END DO NOWAIT
   567:           c$OMP END PARALLEL
   568:           
   569:                 return
   570:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:282)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:282)
    *** The number of VLOAD, VSTORE. :  0,  1. (transfer.f:282)
  LOOP END

  LOOP BEGIN: (transfer.f:288)
    <Unvectorized loop.>

    LOOP BEGIN: (transfer.f:289)
      <Unvectorized loop.>

      LOOP BEGIN: (transfer.f:435)
        <Unvectorized loop.>
      LOOP END

      LOOP BEGIN: (transfer.f:450)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:451)
          <Unvectorized loop.>
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:479)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:480)
          <Unvectorized loop.>
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:508)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:509)
          <Unvectorized loop.>
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:537)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:538)
          <Unvectorized loop.>
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:336)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (transfer.f:336)
        *** The number of VLOAD, VSTORE. :  1, 29. (transfer.f:336)
        *** Iteration count is assumed. Iteration count. : 5 (transfer.f:336)

        LOOP BEGIN: (transfer.f:338)
          <Vectorized loop.>
          **  Runtime-tested short loop. (transfer.f:338)
          *** The number of VGT,   VSC.    :  4,  0. (transfer.f:338)
          *** The number of VLOAD, VSTORE. :  8,  6. (transfer.f:338)
          *** VGT generated (transfer.f:342)
          *** VGT generated (transfer.f:350)
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:370)
        <Unvectorized loop.>

        LOOP BEGIN: (transfer.f:372)
          <Unvectorized loop.>

          LOOP BEGIN: (transfer.f:375)
            <Unvectorized loop.>

            LOOP BEGIN: (transfer.f:391)
              <Unvectorized loop.>
            LOOP END
          LOOP END

          LOOP BEGIN: (transfer.f:415)
            <Unvectorized loop.>
          LOOP END
        LOOP END
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 64 [v0-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 326976 bytes
      Register spill area      : 163112 bytes
      Parameter area           :     24 bytes
      Register save area       :    176 bytes
      User data area           : 163656 bytes
      Others                   :      8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:282)
    *** Estimated execution cycle                       : 39
  LOOP END

  LOOP BEGIN: (transfer.f:288)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 5
            Over basic blocks                           : 5
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 3
            Over basic blocks                           : 3
    *** The number of SCALAR REGISTER TRANSFER          : 1

    LOOP BEGIN: (transfer.f:289)
      *** Estimated execution cycle                     : 1352
      *** The number of VECTOR REGISTER SPILL
            Total                                       : 3
              Over basic blocks                         : 3
      *** The number of VECTOR REGISTER RESTORE
            Total                                       : 3
              Over basic blocks                         : 3
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 115
              Across calls                              :  34
              Not enough registers                      :  17
              Over basic blocks                         :  41
              Others                                    :  23
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 194
              Across calls                              :  51
              Over basic blocks                         : 129
              Others                                    :  14
      *** The number of SCALAR REGISTER TRANSFER        : 110

      LOOP BEGIN: (transfer.f:435)
        *** Estimated execution cycle                   : 197
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 14
                Across calls                            : 10
                Over basic blocks                       :  1
                Others                                  :  3
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 26
                Across calls                            : 15
                Over basic blocks                       :  8
                Others                                  :  3
        *** The number of SCALAR REGISTER TRANSFER      : 7
      LOOP END

      LOOP BEGIN: (transfer.f:450)
        *** Estimated execution cycle                   : 31
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 3
                Over basic blocks                       : 3
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 24
                Over basic blocks                       : 24
        *** The number of SCALAR REGISTER TRANSFER      : 3

        LOOP BEGIN: (transfer.f:451)
          *** Estimated execution cycle                 : 116
          *** The number of SCALAR REGISTER SPILL
                Total                                   : 6
                  Across calls                          : 5
                  Over basic blocks                     : 1
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 5
                  Across calls                          : 5
          *** The number of SCALAR REGISTER TRANSFER    : 4
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:479)
        *** Estimated execution cycle                   : 27
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 9
                Over basic blocks                       : 9
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 17
                Over basic blocks                       : 17
        *** The number of SCALAR REGISTER TRANSFER      : 3

        LOOP BEGIN: (transfer.f:480)
          *** Estimated execution cycle                 : 128
          *** The number of SCALAR REGISTER SPILL
                Total                                   : 8
                  Across calls                          : 6
                  Over basic blocks                     : 2
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 6
                  Across calls                          : 6
          *** The number of SCALAR REGISTER TRANSFER    : 6
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:508)
        *** Estimated execution cycle                   : 28
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 7
                Over basic blocks                       : 7
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 19
                Over basic blocks                       : 19
        *** The number of SCALAR REGISTER TRANSFER      : 4

        LOOP BEGIN: (transfer.f:509)
          *** Estimated execution cycle                 : 132
          *** The number of SCALAR REGISTER SPILL
                Total                                   : 10
                  Across calls                          :  7
                  Over basic blocks                     :  3
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 7
                  Across calls                          : 7
          *** The number of SCALAR REGISTER TRANSFER    : 8
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:537)
        *** Estimated execution cycle                   : 19
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 4
                Over basic blocks                       : 4
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 18
                Over basic blocks                       : 18
        *** The number of SCALAR REGISTER TRANSFER      : 7

        LOOP BEGIN: (transfer.f:538)
          *** Estimated execution cycle                 : 131
          *** The number of SCALAR REGISTER SPILL
                Total                                   : 11
                  Across calls                          :  7
                  Over basic blocks                     :  4
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 7
                  Across calls                          : 7
          *** The number of SCALAR REGISTER TRANSFER    : 9
        LOOP END
      LOOP END

      LOOP BEGIN: (transfer.f:336)
        *** Estimated execution cycle                   : 1140
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 1
                Over basic blocks                       : 1
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 1
                Over basic blocks                       : 1
        *** The number of SCALAR REGISTER TRANSFER      : 1
      LOOP END

      LOOP BEGIN: (transfer.f:370)
        *** Estimated execution cycle                   : 46
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 21
                Over basic blocks                       : 21
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 14
                Over basic blocks                       : 14
        *** The number of SCALAR REGISTER TRANSFER      : 14

        LOOP BEGIN: (transfer.f:372)
          *** Estimated execution cycle                 : 102
          *** The number of SCALAR REGISTER SPILL
                Total                                   : 15
                  Across calls                          :  5
                  Over basic blocks                     : 10
          *** The number of SCALAR REGISTER RESTORE
                Total                                   : 48
                  Across calls                          :  7
                  Over basic blocks                     : 41
          *** The number of SCALAR REGISTER TRANSFER    : 21

          LOOP BEGIN: (transfer.f:375)
            *** Estimated execution cycle               : 82
            *** The number of SCALAR REGISTER SPILL
                  Total                                 : 7
                    Across calls                        : 4
                    Over basic blocks                   : 3
            *** The number of SCALAR REGISTER RESTORE
                  Total                                 : 25
                    Across calls                        :  6
                    Over basic blocks                   : 19
            *** The number of SCALAR REGISTER TRANSFER  : 11

            LOOP BEGIN: (transfer.f:391)
              *** Estimated execution cycle             : 87
              *** The number of SCALAR REGISTER SPILL
                    Total                               : 5
                      Across calls                      : 4
                      Over basic blocks                 : 1
              *** The number of SCALAR REGISTER RESTORE
                    Total                               : 4
                      Across calls                      : 4
              *** The number of SCALAR REGISTER TRANSFER : 4
            LOOP END
          LOOP END

          LOOP BEGIN: (transfer.f:415)
            *** Estimated execution cycle               : 109
            *** The number of SCALAR REGISTER SPILL
                  Total                                 : 12
                    Across calls                        :  5
                    Over basic blocks                   :  7
            *** The number of SCALAR REGISTER RESTORE
                  Total                                 : 5
                    Across calls                        : 5
            *** The number of SCALAR REGISTER TRANSFER  : 10
          LOOP END
        LOOP END
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_E
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   593: opt(1418): Constant-length loop is expanded.
   598: opt(1418): Constant-length loop is expanded.
   604: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_E
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

   574:                 subroutine transfb_cor_e(n,tmor,tx)
   575:           c--------------------------------------------------------------
   576:           c     This subroutine performs the edge to mortar mapping and
   577:           c     calculates the mapping result on the mortar point at a vertex
   578:           c     under situation 1,2, or 3.
   579:           c     n refers to the configuration of three edges sharing a vertex,
   580:           c     n = 1: only one edge is nonconforming
   581:           c     n = 2: two edges are nonconforming
   582:           c     n = 3: three edges are nonconforming
   583:           c-------------------------------------------------------------------
   584:           
   585:                 use ua_data
   586:                 implicit none
   587:           
   588:                 double precision tmor,tx(lx1,lx1,lx1),tmp
   589:                 integer i,n
   590:           
   591:                 tmor=tx(1,1,1)
   592:           
   593: *------>        do i=2,lx1-1
   594: |                 tmor= tmor + qbnew(i-1,1,1)*tx(i,1,1)
   595: *------         end do
   596:           
   597:                 if(n.gt.1)then
   598: *------>          do i=2,lx1-1
   599: |                   tmor= tmor + qbnew(i-1,1,1)*tx(1,i,1)
   600: *------           end do
   601:                 end if
   602:           
   603:                 if(n.eq.3)then
   604: *------>          do i=2,lx1-1
   605: |                   tmor= tmor + qbnew(i-1,1,1)*tx(1,1,i)
   606: *------           end do
   607:                 end if
   608:           
   609:                 return
   610:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_E
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_E
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 15 [s0-s2 s8-s11 s18 s57-s63]
      Vector registers         :  0
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_F
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   632: opt(1418): Constant-length loop is expanded.
   634: opt(1418): Constant-length loop is expanded.
   640: opt(1418): Constant-length loop is expanded.
   645: opt(1418): Constant-length loop is expanded.
   652: opt(1418): Constant-length loop is expanded.
   653: opt(1418): Constant-length loop is expanded.
   658: opt(1418): Constant-length loop is expanded.
   665: opt(1418): Constant-length loop is expanded.
   666: opt(1418): Constant-length loop is expanded.
   670: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_F
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

   613:                 subroutine transfb_cor_f(n,tmor,tx)
   614:           c--------------------------------------------------------------
   615:           c     This subroutine performs the mapping from face to mortar.
   616:           c     Output tmor is the mapping result on a mortar vertex
   617:           c     of situations of three edges and three faces sharing a vertex:
   618:           c     n=4: only one face is nonconforming
   619:           c     n=5: one face and one edge are nonconforming
   620:           c     n=6: two faces are nonconforming
   621:           c     n=7: three faces are nonconforming
   622:           c--------------------------------------------------------------
   623:           
   624:                 use ua_data
   625:                 implicit none
   626:           
   627:                 double precision tx(lx1,lx1,lx1),tmor,temp(lx1)
   628:                 integer col,i,n
   629:           
   630:                 call r_init(temp,lx1,0.d0)
   631:           
   632: *------>        do col=1,lx1
   633: |                 temp(col)=tx(col,1,1)
   634: |*----->          do i=2,lx1-1
   635: ||                  temp(col) = temp(col) + qbnew(i-1,1,1)*tx(col,i,1)
   636: |*-----           end do
   637: *------         end do
   638:                 tmor=temp(1)
   639:           
   640: *------>        do i=2,lx1-1
   641: |                 tmor = tmor + qbnew(i-1,1,1) *temp(i)
   642: *------         end do
   643:           
   644:                 if(n.eq.5)then
   645: *------>          do i=2,lx1-1
   646: |                   tmor = tmor + qbnew(i-1,1,1) *tx(1,1,i)
   647: *------           end do
   648:                 end if
   649:            
   650:                 if(n.ge.6)then
   651:                   call r_init(temp,lx1,0.d0)
   652: *------>          do col=1,lx1
   653: |*----->            do i=2,lx1-1
   654: ||                    temp(col) = temp(col) + qbnew(i-1,1,1)*tx(col,1,i)
   655: |*-----             end do
   656: *------           end do
   657:                   tmor=tmor+temp(1)
   658: *------>          do i=2,lx1-1
   659: |                   tmor = tmor +qbnew(i-1,1,1) *temp(i)
   660: *------           end do
   661:                 end if
   662:            
   663:                 if(n.eq.7)then
   664:                   call r_init(temp,lx1,0.d0)
   665: *------>          do col=2,lx1-1
   666: |*----->            do i=2,lx1-1
   667: ||                    temp(col) = temp(col) + qbnew(i-1,1,1)*tx(1,col,i)
   668: |*-----             end do
   669: *------           end do
   670: *------>          do i=2,lx1-1
   671: |                   tmor = tmor + qbnew(i-1,1,1) *temp(i)
   672: *------           end do
   673:                 end if
   674:           
   675:                 return
   676:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_F
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_COR_F
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 26 [s0-s2 s8-s12 s15-s16 s18 s23-s28 s55-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 288 bytes
      Register spill area      :   0 bytes
      Parameter area           :  24 bytes
      Register save area       : 176 bytes
      User data area           :  80 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF_NC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   696: vec( 101): Vectorized loop.
   699: opt(1418): Constant-length loop is expanded.
   700: opt(1418): Constant-length loop is expanded.
   701: vec( 128): Fused multiply-add operation applied.
   706: vec( 101): Vectorized loop.
   709: opt(1418): Constant-length loop is expanded.
   710: opt(1418): Constant-length loop is expanded.
   711: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF_NC
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

   680:                 subroutine transf_nc(tmor,tx)
   681:           c------------------------------------------------------------------------
   682:           c     Perform mortar to element mapping on a nonconforming face.
   683:           c     This subroutin is used when all entries in tmor are zero except
   684:           c     one tmor(i,j)=1. So this routine is simplified. Only one piece of
   685:           c     mortar  (tmor only has two indices) and one piece of intermediate
   686:           c     mortar (tmp) are involved.
   687:           c------------------------------------------------------------------------
   688:           
   689:                 use ua_data
   690:                 implicit none
   691:           
   692:                 double precision tmor(lx1,lx1), tx(lx1,lx1), tmp(lx1,lx1)
   693:                 integer col,i,j
   694:           
   695:                 call r_init(tmp,lx1*lx1,0.d0)
   696: V------>        do col=1,lx1
   697: |                 i = 1
   698: |                 tmp(i,col)=tmor(i,col)
   699: |*----->          do i=2,lx1-1
   700: ||*---->            do j=1,lx1
   701: |||     F             tmp(i,col) = tmp(i,col) + qbnew(i-1,j,1)*tmor(j,col)
   702: ||*----             end do
   703: |*-----           end do
   704: V------         end do
   705:           
   706: V------>        do col=1,lx1
   707: |                 i = 1
   708: |                 tx(col,i)   = tx(col,i)   + tmp(col,i)
   709: |*----->          do i=2,lx1-1
   710: ||*---->            do j=1,lx1
   711: |||     F             tx(col,i) = tx(col,i) + qbnew(i-1,j,1)*tmp(col,j)
   712: ||*----             end do
   713: |*-----           end do
   714: V------         end do
   715:           
   716:                 return
   717:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF_NC
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:696)
    <Vectorized loop.>
    **  Runtime-tested short loop. (transfer.f:696)
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:696)
    *** The number of VLOAD, VSTORE. :  8,  4. (transfer.f:696)
  LOOP END

  LOOP BEGIN: (transfer.f:706)
    <Vectorized loop.>
    **  Runtime-tested short loop. (transfer.f:706)
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:706)
    *** The number of VLOAD, VSTORE. :  9,  4. (transfer.f:706)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSF_NC
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 42 [s0-s2 s8-s12 s15-s16 s23-s26 s36-s63]
      Vector registers         : 48 [v16-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 65968 bytes
      Register spill area      : 32768 bytes
      Parameter area           :    24 bytes
      Register save area       :   176 bytes
      User data area           : 32992 bytes
      Others                   :     8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC0
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   733: opt(1418): Constant-length loop is expanded.
   734: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC0
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

   720:                 subroutine transfb_nc0(tmor,tx)
   721:           c------------------------------------------------------------------------
   722:           c     Performs mapping from element to mortar when the nonconforming
   723:           c     edges are shared by two conforming faces of an element.
   724:           c------------------------------------------------------------------------
   725:           
   726:                 use ua_data
   727:                 implicit none
   728:           
   729:                 double precision tmor(lx1,lx1),tx(lx1,lx1,lx1)
   730:                 integer i,j
   731:           
   732:                 call r_init(tmor,lx1*lx1,0.d0)
   733: *------>        do j=1,lx1
   734: |*----->          do i=2,lx1-1
   735: ||                  tmor(j,1)= tmor(j,1) + qbnew(i-1,j  ,1)*tx(i,1,1)
   736: |*-----           end do
   737: *------         end do
   738:           
   739:                 return
   740:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC0
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC0
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 19 [s0-s2 s8-s12 s15-s16 s23-s25 s58-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 224 bytes
      Register spill area      :   0 bytes
      Parameter area           :  24 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   764: vec( 101): Vectorized loop.
   768: opt(1418): Constant-length loop is expanded.
   769: vec( 128): Fused multiply-add operation applied.
   772: opt(1418): Constant-length loop is expanded.
   773: opt(1418): Constant-length loop is expanded.
   774: vec( 128): Fused multiply-add operation applied.
   784: opt(1418): Constant-length loop is expanded.
   785: opt(1418): Constant-length loop is expanded.
   791: vec( 103): Unvectorized loop.
   791: vec( 107): Iteration count is too small.
   793: opt(1418): Constant-length loop is expanded.
   794: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC2
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

   743:                 subroutine transfb_nc2(tmor,tx)
   744:           c------------------------------------------------------------------------
   745:           c     Maps values from element to mortar when the nonconforming edges are
   746:           c     shared by two nonconforming faces of an element.
   747:           c     Although each face shall have four pieces of mortar, only value in
   748:           c     one piece (location (1,1)) is used in the calling routine so only
   749:           c     the value in the first mortar is calculated in this subroutine.
   750:           c------------------------------------------------------------------------
   751:           
   752:                 use ua_data
   753:                 implicit none
   754:           
   755:                 double precision tx(lx1,lx1),tmor(lx1,lx1),bottom(lx1),
   756:                &                 temp(lx1,lx1)
   757:                 integer col,j,i
   758:           
   759:                 call r_init(tmor,lx1*lx1,0.d0)
   760:                 call r_init(temp,lx1*lx1,0.d0)
   761:                 tmor(1,1)=tx(1,1)
   762:           
   763:           c.....mapping from tx to intermediate mortar temp + bottom
   764: V------>        do col=1,lx1
   765: |                 temp(col,1)=tx(col,1)
   766: |                 j=1
   767: |                 bottom(col)= 0.d0
   768: |*----->          do i=2,lx1-1
   769: ||      F           bottom(col) = bottom(col) + qbnew(i-1,j,1)*tx(col,i)
   770: |*-----           end do
   771: |         
   772: |*----->          do j=2,lx1
   773: ||*---->            do i=2,lx1-1
   774: |||     F             temp(col,j) = temp(col,j) + qbnew(i-1,j,1)*tx(col,i)
   775: ||*----             end do
   776: |*-----           end do
   777: V------         end do
   778:           
   779:           c.....from intermediate mortar to mortar
   780:           
   781:           c.....On the nonconforming edge, temp is divided by 2 as there will be
   782:           c     a duplicate contribution from another face sharing this edge
   783:                 col=1
   784: *------>        do j=1,lx1
   785: |*----->          do i=2,lx1-1
   786: ||                  tmor(j,col)=tmor(j,col)+ qbnew(i-1,j,1) * bottom(i) +
   787: ||             &                             qbnew(i-1,j,1) * temp(i,col) * 0.5d0
   788: |*-----           end do
   789: *------         end do
   790:           
   791: +------>        do col=2,lx1
   792: |                 tmor(1,col)=tmor(1,col)+temp(1,col)
   793: |*----->          do j=1,lx1
   794: ||*---->            do i=2,lx1-1
   795: |||                   tmor(j,col) = tmor(j,col) + qbnew(i-1,j,1) *temp(i,col)
   796: ||*----             end do
   797: |*-----           end do
   798: +------         end do
   799:           
   800:                 return
   801:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:764)
    <Vectorized loop.>
    **  Runtime-tested short loop. (transfer.f:764)
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:764)
    *** The number of VLOAD, VSTORE. :  8,  6. (transfer.f:764)
  LOOP END

  LOOP BEGIN: (transfer.f:791)
    <Unvectorized loop.>
    *** Iteration count is too small. (transfer.f:791)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 56 [s0-s12 s15-s16 s23-s63]
      Vector registers         : 23 [v41-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 480 bytes
      Register spill area      :   0 bytes
      Parameter area           :  24 bytes
      Register save area       : 176 bytes
      User data area           : 272 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:791)
    *** Estimated execution cycle                       : 61
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   826: vec( 101): Vectorized loop.
   830: opt(1418): Constant-length loop is expanded.
   831: vec( 128): Fused multiply-add operation applied.
   834: opt(1418): Constant-length loop is expanded.
   835: opt(1418): Constant-length loop is expanded.
   836: vec( 128): Fused multiply-add operation applied.
   844: opt(1418): Constant-length loop is expanded.
   845: opt(1418): Constant-length loop is expanded.
   855: vec( 103): Unvectorized loop.
   855: vec( 107): Iteration count is too small.
   857: opt(1418): Constant-length loop is expanded.
   858: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC1
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

   805:                 subroutine transfb_nc1(tmor,tx)
   806:           c------------------------------------------------------------------------
   807:           c     Maps values from element to mortar when the nonconforming edges are
   808:           c     shared by a nonconforming face and a conforming face of an element
   809:           c------------------------------------------------------------------------
   810:           
   811:                 use ua_data
   812:                 implicit none
   813:           
   814:                 double precision tx(lx1,lx1),tmor(lx1,lx1),bottom(lx1),
   815:                &                 temp(lx1,lx1)
   816:                 integer col,j,i
   817:           
   818:                 call r_init(tmor,lx1*lx1,0.d0)
   819:                 call r_init(temp,lx1*lx1,0.d0)
   820:           
   821:                 tmor(1,1)=tx(1,1)
   822:           c.....Contribution from the nonconforming faces
   823:           c     Since the calling subroutine is only interested in the value on the
   824:           c     mortar (location (1,1)), only this piece of mortar is calculated.
   825:           
   826: V------>        do col=1,lx1
   827: |                 temp(col,1)=tx(col,1)
   828: |                 j = 1
   829: |                 bottom(col)= 0.d0
   830: |*----->          do i=2,lx1-1
   831: ||      F           bottom(col)=bottom(col) + qbnew(i-1,j,1)*tx(col,i)
   832: |*-----           end do
   833: |         
   834: |*----->          do j=2,lx1
   835: ||*---->            do i=2,lx1-1
   836: |||     F             temp(col,j) = temp(col,j) + qbnew(i-1,j,1)*tx(col,i)
   837: ||*----             end do
   838: ||        
   839: |*-----           end do
   840: V------         end do
   841:           
   842:                 col=1
   843:                 tmor(1,col)=tmor(1,col)+bottom(1)
   844: *------>        do j=1,lx1
   845: |*----->          do i=2,lx1-1
   846: ||        
   847: ||        c.........temp is not divided by 2 here. It includes the contribution
   848: ||        c         from the other conforming face.
   849: ||        
   850: ||                  tmor(j,col)=tmor(j,col) + qbnew(i-1,j,1) *bottom(i) +
   851: ||             &                              qbnew(i-1,j,1) *temp(i,col)
   852: |*-----           end do
   853: *------         end do
   854:           
   855: +------>        do col=2,lx1
   856: |                 tmor(1,col)=tmor(1,col)+temp(1,col)
   857: |*----->          do j=1,lx1
   858: ||*---->            do i=2,lx1-1
   859: |||                   tmor(j,col) = tmor(j,col) + qbnew(i-1,j,1) *temp(i,col)
   860: ||*----             end do
   861: |*-----           end do
   862: +------         end do
   863:           
   864:                 return
   865:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC1
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:826)
    <Vectorized loop.>
    **  Runtime-tested short loop. (transfer.f:826)
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:826)
    *** The number of VLOAD, VSTORE. :  8,  6. (transfer.f:826)
  LOOP END

  LOOP BEGIN: (transfer.f:855)
    <Unvectorized loop.>
    *** Iteration count is too small. (transfer.f:855)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:12 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_NC1
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 54 [s0-s2 s4-s12 s15-s16 s23-s32 s34-s63]
      Vector registers         : 23 [v41-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 480 bytes
      Register spill area      :   0 bytes
      Parameter area           :  24 bytes
      Register save area       : 176 bytes
      User data area           : 272 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:855)
    *** Estimated execution cycle                       : 61
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   883: par(1801): Parallel routine generated.: TRANSFB_C$1
   887: par(1803): Parallelized by "do".
   887: vec( 101): Vectorized loop.
   893: opt(3008): Reference within a conditional branch moved outside loop.
   893: par(1803): Parallelized by "do".
   894: opt(3008): Reference within a conditional branch moved outside loop.
   921: vec( 103): Unvectorized loop.
   922: opt(1418): Constant-length loop is expanded.
   926: opt(1097): This statement prevents loop optimization.
   933: opt(1418): Constant-length loop is expanded.
   944: opt(1418): Constant-length loop is expanded.
   955: opt(1418): Constant-length loop is expanded.
   966: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C
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

   869:                 subroutine transfb_c(tx)
   870:           c-------------------------------------------------------------------
   871:           c     Prepare initial guess for cg. All values from conforming
   872:           c     boundary are copied and summed on tmor.
   873:           c-------------------------------------------------------------------
   874:           
   875:                 use ua_data
   876:                 implicit none
   877:           
   878:                 double precision third
   879:                 parameter (third = 1.d0/3.d0)
   880:                 double precision tx(*)
   881:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,col,j,ig,il
   882:           
   883:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IE,IFACE,IL1,IL2,
   884:           c$OMP& IL3,IL4,IG1,IG2,IG3,IG4,COL,J,IG,IL)
   885:           
   886:           c$OMP DO
   887: Y------>        do j=1,nmor
   888: |                 tmort(j)=0.d0
   889: Y------         end do
   890:           c$OMP END DO
   891:           
   892:           c$OMP DO
   893: P------>        do ie=1,nelt
   894: |+----->          do iface=1,nsides
   895: ||                  if(cbc(iface,ie).ne.3)then
   896: ||                    il1 = idel(1,1,iface,ie)
   897: ||                    il2 = idel(lx1,1,iface,ie)
   898: ||                    il3 = idel(1,lx1,iface,ie)
   899: ||                    il4 = idel(lx1,lx1,iface,ie)
   900: ||                    ig1 = idmo(1,  1,  1,1,iface,ie)
   901: ||                    ig2 = idmo(lx1,1,  1,2,iface,ie)
   902: ||                    ig3 = idmo(1,  lx1,2,1,iface,ie)
   903: ||                    ig4 = idmo(lx1,lx1,2,2,iface,ie)
   904: ||        c
   905: ||        c$          call omp_set_lock(tlock(ig1))
   906: ||                    tmort(ig1) = tmort(ig1)+tx(il1)*third
   907: ||        c$          call omp_unset_lock(tlock(ig1))
   908: ||        c
   909: ||        c$          call omp_set_lock(tlock(ig2))
   910: ||                    tmort(ig2) = tmort(ig2)+tx(il2)*third
   911: ||        c$          call omp_unset_lock(tlock(ig2))
   912: ||        c
   913: ||        c$          call omp_set_lock(tlock(ig3))
   914: ||                    tmort(ig3) = tmort(ig3)+tx(il3)*third
   915: ||        c$          call omp_unset_lock(tlock(ig3))
   916: ||        c
   917: ||        c$          call omp_set_lock(tlock(ig4))
   918: ||                    tmort(ig4) = tmort(ig4)+tx(il4)*third
   919: ||        c$          call omp_unset_lock(tlock(ig4))
   920: ||        
   921: ||+---->              do  col=2,lx1-1
   922: |||*--->                do j=2,lx1-1
   923: ||||                      il=idel(j,col,iface,ie)
   924: ||||                      ig=idmo(j,col,1,1,iface,ie)
   925: ||||      c
   926: ||||      c$              call omp_set_lock(tlock(ig))
   927: ||||                      tmort(ig)=tmort(ig)+tx(il)
   928: ||||      c$              call omp_unset_lock(tlock(ig))
   929: |||*---                 end do
   930: ||+----               end do
   931: ||        
   932: ||                    if(idmo(lx1,1,1,1,iface,ie).eq.0)then
   933: ||*---->                do j=2,lx1-1
   934: |||                       il=idel(j,1,iface,ie)
   935: |||                       ig=idmo(j,1,1,1,iface,ie)
   936: |||       c
   937: |||       c$              call omp_set_lock(tlock(ig))
   938: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   939: |||       c$              call omp_unset_lock(tlock(ig))
   940: ||*----                 end do
   941: ||                    end if
   942: ||        
   943: ||                    if(idmo(lx1,2,1,2,iface,ie).eq.0)then
   944: ||*---->                do j=2,lx1-1
   945: |||                       il=idel(lx1,j,iface,ie)
   946: |||                       ig=idmo(lx1,j,1,1,iface,ie)
   947: |||       c
   948: |||       c$              call omp_set_lock(tlock(ig))
   949: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   950: |||       c$              call omp_unset_lock(tlock(ig))
   951: ||*----                 end do
   952: ||                    end if
   953: ||        
   954: ||                    if(idmo(2,lx1,2,1,iface,ie).eq.0)then
   955: ||*---->                do j=2,lx1-1
   956: |||                       il=idel(j,lx1,iface,ie)
   957: |||                       ig=idmo(j,lx1,1,1,iface,ie)
   958: |||       c
   959: |||       c$              call omp_set_lock(tlock(ig))
   960: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   961: |||       c$              call omp_unset_lock(tlock(ig))
   962: ||*----                 end do
   963: ||                    end if
   964: ||        
   965: ||                    if(idmo(1,lx1,1,1,iface,ie).eq.0)then
   966: ||*---->                do j=2,lx1-1
   967: |||                       il=idel(1,j,iface,ie)
   968: |||                       ig=idmo(1,j,1,1,iface,ie)
   969: |||       c
   970: |||       c$              call omp_set_lock(tlock(ig))
   971: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   972: |||       c$              call omp_unset_lock(tlock(ig))
   973: ||*----                 end do
   974: ||                    end if
   975: ||                  end if!
   976: |+-----           end do
   977: P------         end do
   978:           c$OMP END DO NOWAIT
   979:           c$OMP END PARALLEL
   980:                 return
   981:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:887)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:887)
    *** The number of VLOAD, VSTORE. :  0,  1. (transfer.f:887)
  LOOP END

  LOOP BEGIN: (transfer.f:893)
    <Unvectorized loop.>

    LOOP BEGIN: (transfer.f:894)
      <Unvectorized loop.>

      LOOP BEGIN: (transfer.f:921)
        <Unvectorized loop.>
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 49 [s0-s2 s8-s12 s15-s16 s18-s44 s52-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1072 bytes
      Register spill area      :  432 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           :  448 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:887)
    *** Estimated execution cycle                       : 39
  LOOP END

  LOOP BEGIN: (transfer.f:893)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 3
            Over basic blocks                           : 3
    *** The number of SCALAR REGISTER TRANSFER          : 1

    LOOP BEGIN: (transfer.f:894)
      *** Estimated execution cycle                     : 195
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 32
              Across calls                              : 18
              Over basic blocks                         :  4
              Others                                    : 10
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 71
              Across calls                              : 21
              Over basic blocks                         : 40
              Others                                    : 10
      *** The number of SCALAR REGISTER TRANSFER        : 34

      LOOP BEGIN: (transfer.f:921)
        *** Estimated execution cycle                   : 160
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 7
                Across calls                            : 3
                Others                                  : 4
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 9
                Across calls                            : 3
                Over basic blocks                       : 3
                Others                                  : 3
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C_2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   999: par(1801): Parallel routine generated.: TRANSFB_C_2$1
  1003: par(1803): Parallelized by "do".
  1003: vec( 101): Vectorized loop.
  1008: par(1803): Parallelized by "do".
  1008: vec( 101): Vectorized loop.
  1014: opt(3008): Reference within a conditional branch moved outside loop.
  1014: par(1803): Parallelized by "do".
  1015: opt(3008): Reference within a conditional branch moved outside loop.
  1047: vec( 103): Unvectorized loop.
  1048: opt(1418): Constant-length loop is expanded.
  1052: opt(1097): This statement prevents loop optimization.
  1060: opt(1418): Constant-length loop is expanded.
  1072: opt(1418): Constant-length loop is expanded.
  1084: opt(1418): Constant-length loop is expanded.
  1096: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C_2
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

   984:                 subroutine transfb_c_2(tx)
   985:           c-------------------------------------------------------------------
   986:           c     Prepare initial guess for CG. All values from conforming
   987:           c     boundary are copied and summed in tmort.
   988:           c     mormult is multiplicity, which is used to average tmort.
   989:           c-------------------------------------------------------------------
   990:           
   991:                 use ua_data
   992:                 implicit none
   993:           
   994:                 double precision third
   995:                 parameter (third = 1.d0/3.d0)
   996:                 double precision tx(*)
   997:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,col,j,ig,il
   998:           
   999:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IE,IFACE,IL1,IL2,
  1000:           c$OMP& IL3,IL4,IG1,IG2,IG3,IG4,COL,J,IG,IL)
  1001:           
  1002:           c$OMP DO
  1003: Y------>        do j=1,nmor
  1004: |                 tmort(j)=0.d0
  1005: Y------         end do
  1006:           c$OMP END DO nowait
  1007:           c$OMP DO
  1008: Y------>        do j=1,nmor
  1009: |                 mormult(j)=0.d0
  1010: Y------         end do
  1011:           c$OMP END DO
  1012:           
  1013:           c$OMP DO
  1014: P------>        do ie=1,nelt
  1015: |+----->          do iface=1,nsides
  1016: ||         
  1017: ||                  if(cbc(iface,ie).ne.3)then
  1018: ||                    il1 = idel(1,  1,  iface,ie)
  1019: ||                    il2 = idel(lx1,1,  iface,ie)
  1020: ||                    il3 = idel(1,  lx1,iface,ie)
  1021: ||                    il4 = idel(lx1,lx1,iface,ie)
  1022: ||                    ig1 = idmo(1,  1,  1,1,iface,ie)
  1023: ||                    ig2 = idmo(lx1,1,  1,2,iface,ie)
  1024: ||                    ig3 = idmo(1,  lx1,2,1,iface,ie)
  1025: ||                    ig4 = idmo(lx1,lx1,2,2,iface,ie)
  1026: ||        c
  1027: ||        c$          call omp_set_lock(tlock(ig1))
  1028: ||                    tmort(ig1) = tmort(ig1)+tx(il1)*third
  1029: ||                    mormult(ig1) = mormult(ig1)+third
  1030: ||        c$          call omp_unset_lock(tlock(ig1))
  1031: ||        c
  1032: ||        c$          call omp_set_lock(tlock(ig2))
  1033: ||                    tmort(ig2) = tmort(ig2)+tx(il2)*third
  1034: ||                    mormult(ig2) = mormult(ig2)+third
  1035: ||        c$          call omp_unset_lock(tlock(ig2))
  1036: ||        c
  1037: ||        c$          call omp_set_lock(tlock(ig3))
  1038: ||                    tmort(ig3) = tmort(ig3)+tx(il3)*third
  1039: ||                    mormult(ig3) = mormult(ig3)+third
  1040: ||        c$          call omp_unset_lock(tlock(ig3))
  1041: ||        c
  1042: ||        c$          call omp_set_lock(tlock(ig4))
  1043: ||                    tmort(ig4) = tmort(ig4)+tx(il4)*third
  1044: ||                    mormult(ig4) = mormult(ig4)+third
  1045: ||        c$          call omp_unset_lock(tlock(ig4))
  1046: ||        
  1047: ||+---->              do  col=2,lx1-1
  1048: |||*--->                do j=2,lx1-1
  1049: ||||                      il=idel(j,col,iface,ie)
  1050: ||||                      ig=idmo(j,col,1,1,iface,ie)
  1051: ||||      c
  1052: ||||      c$              call omp_set_lock(tlock(ig))
  1053: ||||                      tmort(ig)=tmort(ig)+tx(il)
  1054: ||||                      mormult(ig)=mormult(ig)+1.d0
  1055: ||||      c$              call omp_unset_lock(tlock(ig))
  1056: |||*---                 end do
  1057: ||+----               end do
  1058: ||        
  1059: ||                    if(idmo(lx1,1,1,1,iface,ie).eq.0)then
  1060: ||*---->                do j=2,lx1-1
  1061: |||                       il=idel(j,1,iface,ie)
  1062: |||                       ig=idmo(j,1,1,1,iface,ie)
  1063: |||       c
  1064: |||       c$              call omp_set_lock(tlock(ig))
  1065: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1066: |||                       mormult(ig)=mormult(ig)+0.5d0
  1067: |||       c$              call omp_unset_lock(tlock(ig))
  1068: ||*----                 end do
  1069: ||                    end if
  1070: ||        
  1071: ||                    if(idmo(lx1,2,1,2,iface,ie).eq.0)then
  1072: ||*---->                do j=2,lx1-1
  1073: |||                       il=idel(lx1,j,iface,ie)
  1074: |||                       ig=idmo(lx1,j,1,1,iface,ie)
  1075: |||       c
  1076: |||       c$              call omp_set_lock(tlock(ig))
  1077: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1078: |||                       mormult(ig)=mormult(ig)+0.5d0
  1079: |||       c$              call omp_unset_lock(tlock(ig))
  1080: ||*----                 end do
  1081: ||                    end if
  1082: ||        
  1083: ||                    if(idmo(2,lx1,2,1,iface,ie).eq.0)then
  1084: ||*---->                do j=2,lx1-1
  1085: |||                       il=idel(j,lx1,iface,ie)
  1086: |||                       ig=idmo(j,lx1,1,1,iface,ie)
  1087: |||       c
  1088: |||       c$              call omp_set_lock(tlock(ig))
  1089: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1090: |||                       mormult(ig)=mormult(ig)+0.5d0
  1091: |||       c$              call omp_unset_lock(tlock(ig))
  1092: ||*----                  end do
  1093: ||                    end if
  1094: ||        
  1095: ||                    if(idmo(1,lx1,1,1,iface,ie).eq.0)then
  1096: ||*---->                do j=2,lx1-1
  1097: |||                       il=idel(1,j,iface,ie)
  1098: |||                       ig=idmo(1,j,1,1,iface,ie)
  1099: |||       c
  1100: |||       c$              call omp_set_lock(tlock(ig))
  1101: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1102: |||                       mormult(ig)=mormult(ig)+0.5d0
  1103: |||       c$              call omp_unset_lock(tlock(ig))
  1104: ||*----                 end do
  1105: ||                    end if
  1106: ||                  end if!nnje=1
  1107: |+-----           end do
  1108: P------         end do
  1109:           c$OMP END DO NOWAIT
  1110:           c$OMP END PARALLEL
  1111:           
  1112:                 return
  1113:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C_2
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (transfer.f:1003)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:1003)
    *** The number of VLOAD, VSTORE. :  0,  1. (transfer.f:1003)
  LOOP END

  LOOP BEGIN: (transfer.f:1008)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (transfer.f:1008)
    *** The number of VLOAD, VSTORE. :  0,  1. (transfer.f:1008)
  LOOP END

  LOOP BEGIN: (transfer.f:1014)
    <Unvectorized loop.>

    LOOP BEGIN: (transfer.f:1015)
      <Unvectorized loop.>

      LOOP BEGIN: (transfer.f:1047)
        <Unvectorized loop.>
      LOOP END
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:13 2021
FILE NAME: transfer.f

PROCEDURE NAME: TRANSFB_C_2
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 52 [s0-s2 s8-s12 s15-s16 s18-s44 s49-s63]
      Vector registers         :  2 [v62-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1376 bytes
      Register spill area      :  584 bytes
      Parameter area           :    8 bytes
      Register save area       :  176 bytes
      User data area           :  600 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (transfer.f:1003)
    *** Estimated execution cycle                       : 39
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END

  LOOP BEGIN: (transfer.f:1008)
    *** Estimated execution cycle                       : 39
  LOOP END

  LOOP BEGIN: (transfer.f:1014)
    *** Estimated execution cycle                       : 11
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 3
            Over basic blocks                           : 3
    *** The number of SCALAR REGISTER TRANSFER          : 1

    LOOP BEGIN: (transfer.f:1015)
      *** Estimated execution cycle                     : 218
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 53
              Across calls                              : 35
              Over basic blocks                         :  3
              Others                                    : 15
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 122
              Across calls                              :  38
              Over basic blocks                         :  72
              Others                                    :  12
      *** The number of SCALAR REGISTER TRANSFER        : 42

      LOOP BEGIN: (transfer.f:1047)
        *** Estimated execution cycle                   : 167
        *** The number of SCALAR REGISTER SPILL
              Total                                     : 8
                Across calls                            : 4
                Others                                  : 4
        *** The number of SCALAR REGISTER RESTORE
              Total                                     : 12
                Across calls                            :  4
                Over basic blocks                       :  5
                Others                                  :  3
        *** The number of SCALAR REGISTER TRANSFER      : 1
      LOOP END
    LOOP END
  LOOP END


