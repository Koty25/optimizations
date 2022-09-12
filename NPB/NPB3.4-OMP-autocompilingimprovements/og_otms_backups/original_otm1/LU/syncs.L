NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME : syncs.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
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

    28:                 subroutine sync_init( jdim )
    29:           
    30:           c---------------------------------------------------------------------
    31:           c   Initialize sync-related variables
    32:           c---------------------------------------------------------------------
    33:           
    34:                 use syncs
    35:                 implicit none
    36:           
    37:                 integer jdim
    38:           
    39:           !$    integer, external :: omp_get_num_threads, omp_get_thread_num
    40:           
    41:                 mthreadnum = 0
    42:           !$    mthreadnum = omp_get_num_threads() - 1
    43:                 if (mthreadnum .gt. jdim) mthreadnum = jdim
    44:                 iam = 0
    45:           !$    iam = omp_get_thread_num()
    46:                 if (iam .le. mthreadnum) isync(1,iam) = 0
    47:           
    48:                 return
    49:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_INIT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 15 [s0 s8-s12 s15-s16 s23-s26 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    79: vec( 103): Unvectorized loop.
    79: vec( 109): Vectorization obstructive statement.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
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

    57:                 subroutine sync_left( ldmx, ldmy, ldmz, v )
    58:           
    59:           c---------------------------------------------------------------------
    60:           c   Thread synchronization for pipeline operation
    61:           c---------------------------------------------------------------------
    62:           
    63:                 use syncs
    64:                 implicit none
    65:           
    66:                 integer ldmx, ldmy, ldmz
    67:                 double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)
    68:           
    69:           c---------------------------------------------------------------------
    70:           c---------------------------------------------------------------------
    71:           
    72:                 integer neigh, iv
    73:           
    74:           
    75:                 if (iam .gt. 0 .and. iam .le. mthreadnum) then
    76:                    neigh = iam - 1
    77:           !$omp atomic read
    78:                    iv = isync(1,neigh)
    79: +------>           do while (iv .eq. 0)
    80: |         !$omp atomic read
    81: |                     iv = isync(1,neigh)
    82: +------            end do
    83:           !$omp atomic write
    84:                    isync(1,neigh) = 0
    85:                 endif
    86:           !$omp flush(isync,v)
    87:           
    88:           
    89:                 return
    90:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
VECTORIZATION LIST

  LOOP BEGIN: (syncs.f:79)
    <Unvectorized loop.>
    *** Vectorization obstructive statement. (syncs.f:79)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_LEFT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 20 [s0 s8-s12 s15-s16 s18-s19 s23-s29 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (syncs.f:79)
    *** Estimated execution cycle                       : 6
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   117: vec( 103): Unvectorized loop.
   117: vec( 109): Vectorization obstructive statement.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
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

    95:                 subroutine sync_right( ldmx, ldmy, ldmz, v )
    96:           
    97:           c---------------------------------------------------------------------
    98:           c   Thread synchronization for pipeline operation
    99:           c---------------------------------------------------------------------
   100:           
   101:                 use syncs
   102:                 implicit none
   103:           
   104:                 integer ldmx, ldmy, ldmz
   105:                 double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)
   106:           
   107:           c---------------------------------------------------------------------
   108:           c---------------------------------------------------------------------
   109:           
   110:                 integer iv
   111:           
   112:           
   113:           !$omp flush(isync,v)
   114:                 if (iam .lt. mthreadnum) then
   115:           !$omp atomic read
   116:                    iv = isync(1,iam)
   117: +------>           do while (iv .eq. 1)
   118: |         !$omp atomic read
   119: |                     iv = isync(1,iam)
   120: +------            end do
   121:           !$omp atomic write
   122:                    isync(1,iam) = 1
   123:                 endif
   124:           
   125:           
   126:                 return
   127:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
VECTORIZATION LIST

  LOOP BEGIN: (syncs.f:117)
    <Unvectorized loop.>
    *** Vectorization obstructive statement. (syncs.f:117)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:02 2021
FILE NAME: syncs.f

PROCEDURE NAME: SYNC_RIGHT
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 16 [s0 s8-s12 s15-s16 s23-s27 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   8 bytes
      Parameter area           :   0 bytes
      Register save area       : 176 bytes
      User data area           :   8 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (syncs.f:117)
    *** Estimated execution cycle                       : 14
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


