NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME : move.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    15: inl(1212): Source for routine not found.: NR_INIT_OMP
    18: par(1801): Parallel routine generated.: MOVE$1
    21: opt(3008): Reference within a conditional branch moved outside loop.
    21: vec( 110): Vectorization obstructive procedure reference.: COPY
    21: par(1803): Parallelized by "do".
    24: inl(1212): Source for routine not found.: COPY
    28: opt(3008): Reference within a conditional branch moved outside loop.
    28: vec( 101): Vectorized loop.
    49: opt(1418): Constant-length loop is expanded.
    50: opt(1418): Constant-length loop is expanded.
    73: vec( 110): Vectorization obstructive procedure reference.: COPY
    73: vec( 110): Vectorization obstructive procedure reference.: NCOPY
    73: par(1803): Parallelized by "do".
    73: vec( 103): Unvectorized loop.
    78: inl(1212): Source for routine not found.: NCOPY


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
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

     2:                 subroutine move
     3:           c---------------------------------------------------------------
     4:           c     move element to proper location in morton space filling curve
     5:           c---------------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:           
    10:                 integer i,iside,jface,iel,ntemp,ii1,ii2,n1,n2,cb
    11:           
    12:                 n2=2*6*nelt
    13:                 n1=n2*2
    14:           
    15:                 call nr_init_omp(sje_new,n1,0)
    16:                 call nr_init_omp(ijel_new,n2,0)
    17:           
    18:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iel,i,iside,jface,cb,ntemp,
    19:           c$OMP& ii1,ii2)
    20:           c$OMP DO
    21: P------>        do iel=1,nelt
    22: |                 i=mt_to_id(iel)
    23: |                 treenew(iel)=tree(i)
    24: |                 call copy(xc_new(1,iel),xc(1,i),8)
    25: |                 call copy(yc_new(1,iel),yc(1,i),8)
    26: |                 call copy(zc_new(1,iel),zc(1,i),8)
    27: |         
    28: |V----->          do iside=1,nsides
    29: ||                  jface = jjface(iside)
    30: ||                  cb=cbc(iside,i)
    31: ||                  xc_new(iside,iel)=xc(iside,i)
    32: ||                  yc_new(iside,iel)=yc(iside,i)
    33: ||                  zc_new(iside,iel)=zc(iside,i)
    34: ||                  cbc_new(iside,iel)=cb
    35: ||        
    36: ||                  if(cb.eq.2)then
    37: ||                    ntemp=sje(1,1,iside,i)
    38: ||                    ijel_new(1,iside,iel)=1
    39: ||                    ijel_new(2,iside,iel)=1
    40: ||      G             sje_new(1,1,iside,iel)=id_to_mt(ntemp)
    41: ||        
    42: ||                  else if(cb.eq.1) then
    43: ||                    ntemp=sje(1,1,iside,i)
    44: ||                    ijel_new(1,iside,iel)=ijel(1,iside,i)
    45: ||                    ijel_new(2,iside,iel)=ijel(2,iside,i)
    46: ||      G             sje_new(1,1,iside,iel)=id_to_mt(ntemp)
    47: ||         
    48: ||                  else if(cb.eq.3) then
    49: ||*---->              do ii2=1,2
    50: |||*--->                do ii1=1,2
    51: ||||                      ntemp=sje(ii1,ii2,iside,i)
    52: ||||                      ijel_new(1,iside,iel)=1
    53: ||||                      ijel_new(2,iside,iel)=1
    54: ||||    G                 sje_new(ii1,ii2,iside,iel)=id_to_mt(ntemp)
    55: |||*---                 end do
    56: ||*----               end do
    57: ||        
    58: ||                  else if(cb.eq.0)then
    59: ||                    sje_new(1,1,iside,iel)=0
    60: ||                    sje_new(1,2,iside,iel)=0
    61: ||                    sje_new(2,1,iside,iel)=0
    62: ||                    sje_new(2,2,iside,iel)=0
    63: ||                  end if
    64: ||        
    65: |V-----           end do
    66: |         
    67: |                 call copy(ta2(1,1,1,iel),ta1(1,1,1,i),nxyz)
    68: |         
    69: P------         end do
    70:           c$OMP ENDDO
    71:           
    72:           c$OMP DO
    73: P------>        do iel=1,nelt
    74: |                 call copy(xc(1,iel),xc_new(1,iel),8)
    75: |                 call copy(yc(1,iel),yc_new(1,iel),8)
    76: |                 call copy(zc(1,iel),zc_new(1,iel),8)
    77: |                 call copy(ta1(1,1,1,iel),ta2(1,1,1,iel),nxyz)
    78: |                 call ncopy(sje(1,1,1,iel),sje_new(1,1,1,iel),4*6)
    79: |                 call ncopy(ijel(1,1,iel),ijel_new(1,1,iel),2*6)
    80: |                 call ncopy(cbc(1,iel),cbc_new(1,iel),6)
    81: |                 mt_to_id(iel)=iel
    82: |                 id_to_mt(iel)=iel
    83: |                 tree(iel)=treenew(iel)
    84: P------         end do
    85:           c$OMP ENDDO
    86:           c$OMP END PARALLEL
    87:           
    88:                 return
    89:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
INLINE LIST

  ROOT: MOVE (move.f:2)
  -> NOINLINE: NR_INIT_OMP (move.f:15)
     *** Source for routine not found.
  -> NOINLINE: NR_INIT_OMP (move.f:16)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:24)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:25)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:26)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:67)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:74)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:75)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:76)
     *** Source for routine not found.
  -> NOINLINE: COPY (move.f:77)
     *** Source for routine not found.
  -> NOINLINE: NCOPY (move.f:78)
     *** Source for routine not found.
  -> NOINLINE: NCOPY (move.f:79)
     *** Source for routine not found.
  -> NOINLINE: NCOPY (move.f:80)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (move.f:21)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : COPY (move.f:21)

    LOOP BEGIN: (move.f:28)
      <Vectorized loop.>
      **  Runtime-tested short loop. (move.f:28)
      *** The number of VGT,   VSC.    :  6,  0. (move.f:28)
      *** The number of VLOAD, VSTORE. : 12, 20. (move.f:28)
      *** VGT generated (move.f:54)
      *** VGT generated (move.f:46)
      *** VGT generated (move.f:40)
    LOOP END
  LOOP END

  LOOP BEGIN: (move.f:73)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : COPY (move.f:73)
    *** Vectorization obstructive procedure reference. : NCOPY (move.f:73)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:37:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 47 [v17-v63]
      Vector mask registers    :  4 [vm12-vm15]
      VREG directive           :  0

  Routine stack
    Total size                 : 120160 bytes
      Register spill area      :  59944 bytes
      Parameter area           :     24 bytes
      Register save area       :    176 bytes
      User data area           :  60008 bytes
      Others                   :      8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (move.f:21)
    *** Estimated execution cycle                       : 2047
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 31
            Across calls                                : 10
            Not enough registers                        : 21
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 72
            Across calls                                : 15
            Over basic blocks                           : 57
    *** The number of SCALAR REGISTER TRANSFER          : 9
  LOOP END

  LOOP BEGIN: (move.f:73)
    *** Estimated execution cycle                       : 167
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 9
            Across calls                                : 6
            Over basic blocks                           : 3
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 30
            Across calls                                : 12
            Over basic blocks                           : 18
    *** The number of SCALAR REGISTER TRANSFER          : 8
  LOOP END


