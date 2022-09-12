NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME : setiv.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setiv.f

PROCEDURE NAME: SETIV
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setiv.f

PROCEDURE NAME: SETIV
VECTORIZATION LIST

  LOOP BEGIN: (setiv.f:33)
    <Unvectorized loop.>

    LOOP BEGIN: (setiv.f:37)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : EXACT (setiv.f:37)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: setiv.f

PROCEDURE NAME: SETIV
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 52 [s0-s3 s8-s12 s15-s16 s18-s44 s50-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 688 bytes
      Register spill area      : 104 bytes
      Parameter area           :  32 bytes
      Register save area       : 176 bytes
      User data area           : 376 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (setiv.f:33)
    *** Estimated execution cycle                       : 78
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 1
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (setiv.f:37)
      *** Estimated execution cycle                     : 200
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 11
              Across calls                              : 11
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 36
              Across calls                              : 15
              Over basic blocks                         : 21
      *** The number of SCALAR REGISTER TRANSFER        : 33
    LOOP END
  LOOP END


