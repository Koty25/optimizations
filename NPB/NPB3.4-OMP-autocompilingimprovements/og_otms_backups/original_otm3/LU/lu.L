NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME : lu.f

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
FILE NAME: lu.f

PROCEDURE NAME: APPLU
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   171: vec( 103): Unvectorized loop.
   171: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   179: opt(3008): Reference within a conditional branch moved outside loop.
   179: vec( 103): Unvectorized loop.
   179: vec( 180): I/O statement obstructs vectorization.
   181: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: lu.f

PROCEDURE NAME: APPLU
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

    49:                 program applu
    50:           c---------------------------------------------------------------------
    51:           
    52:           c---------------------------------------------------------------------
    53:           c
    54:           c   driver for the performance evaluation of the solver for
    55:           c   five coupled parabolic/elliptic partial differential equations.
    56:           c
    57:           c---------------------------------------------------------------------
    58:           
    59:                 use lu_data
    60:                 implicit none
    61:           
    62:                 character class
    63:                 logical verified
    64:                 double precision mflops
    65:           
    66:                 double precision t, tmax, timer_read, trecs(t_last)
    67:                 external timer_read
    68:                 integer i
    69:                 character t_names(t_last)*8
    70:           
    71:           c---------------------------------------------------------------------
    72:           c     Setup info for timers
    73:           c---------------------------------------------------------------------
    74:           
    75:                 call check_timer_flag( timeron )
    76:                 if (timeron) then
    77:                    t_names(t_total) = 'total'
    78:                    t_names(t_rhsx) = 'rhsx'
    79:                    t_names(t_rhsy) = 'rhsy'
    80:                    t_names(t_rhsz) = 'rhsz'
    81:                    t_names(t_rhs) = 'rhs'
    82:                    t_names(t_jacld) = 'jacld'
    83:                    t_names(t_blts) = 'blts'
    84:                    t_names(t_jacu) = 'jacu'
    85:                    t_names(t_buts) = 'buts'
    86:                    t_names(t_add) = 'add'
    87:                    t_names(t_l2norm) = 'l2norm'
    88:                 endif
    89:           
    90:           c---------------------------------------------------------------------
    91:           c   read input data
    92:           c---------------------------------------------------------------------
    93:                 call read_input()
    94:           
    95:           c---------------------------------------------------------------------
    96:           c   set up domain sizes
    97:           c---------------------------------------------------------------------
    98:                 call domain()
    99:           
   100:                 call alloc_space
   101:           
   102:           c---------------------------------------------------------------------
   103:           c   set up coefficients
   104:           c---------------------------------------------------------------------
   105:                 call setcoeff()
   106:           
   107:           c---------------------------------------------------------------------
   108:           c   set the boundary values for dependent variables
   109:           c---------------------------------------------------------------------
   110:                 call setbv()
   111:           
   112:           c---------------------------------------------------------------------
   113:           c   set the initial values for dependent variables
   114:           c---------------------------------------------------------------------
   115:                 call setiv()
   116:           
   117:           c---------------------------------------------------------------------
   118:           c   compute the forcing term based on prescribed exact solution
   119:           c---------------------------------------------------------------------
   120:                 call erhs()
   121:           
   122:           c---------------------------------------------------------------------
   123:           c   perform one SSOR iteration to touch all data pages
   124:           c---------------------------------------------------------------------
   125:                 call ssor(1)
   126:           
   127:           c---------------------------------------------------------------------
   128:           c   reset the boundary and initial values
   129:           c---------------------------------------------------------------------
   130:                 call setbv()
   131:                 call setiv()
   132:           
   133:           c---------------------------------------------------------------------
   134:           c   perform the SSOR iterations
   135:           c---------------------------------------------------------------------
   136:                 call ssor(itmax)
   137:           
   138:           c---------------------------------------------------------------------
   139:           c   compute the solution error
   140:           c---------------------------------------------------------------------
   141:                 call error()
   142:           
   143:           c---------------------------------------------------------------------
   144:           c   compute the surface integral
   145:           c---------------------------------------------------------------------
   146:                 call pintgr()
   147:           
   148:           c---------------------------------------------------------------------
   149:           c   verification test
   150:           c---------------------------------------------------------------------
   151:                 call verify ( rsdnm, errnm, frc, class, verified )
   152:                 mflops = 1.0d-6*dble(itmax)*(1984.77*dble( nx0 )
   153:                >     *dble( ny0 )
   154:                >     *dble( nz0 )
   155:                >     -10923.3*(dble( nx0+ny0+nz0 )/3.)**2
   156:                >     +27770.9* dble( nx0+ny0+nz0 )/3.
   157:                >     -144010.)
   158:                >     / maxtime
   159:           
   160:                 call print_results('LU', class, nx0,
   161:                >  ny0, nz0, itmax,
   162:                >  maxtime, mflops, '          floating point', verified,
   163:                >  npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6,
   164:                >  '(none)')
   165:           
   166:           c---------------------------------------------------------------------
   167:           c      More timers
   168:           c---------------------------------------------------------------------
   169:                 if (.not.timeron) goto 999
   170:           
   171: +------>        do i=1, t_last
   172: |                  trecs(i) = timer_read(i)
   173: +------         end do
   174:                 tmax = maxtime
   175:                 if ( tmax .eq. 0. ) tmax = 1.0
   176:           
   177:                 write(*,800)
   178:            800  format('  SECTION     Time (secs)')
   179: +------>        do i=1, t_last
   180: |                  if (i.ne.t_jacld .and. i.ne.t_jacu) then
   181: |                     write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
   182: |                  endif
   183: |                  if (i.eq.t_rhs) then
   184: |                     t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
   185: |                     write(*,820) 'sub-rhs', t, t*100./tmax
   186: |                     t = trecs(i) - t
   187: |                     write(*,820) 'rest-rhs', t, t*100./tmax
   188: |                  endif
   189: |          810     format(2x,a8,':',f9.3,'  (',f6.2,'%)')
   190: |          820     format(5x,'--> ',a8,':',f9.3,'  (',f6.2,'%)')
   191: +------         end do
   192:           
   193:            999  continue
   194:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: lu.f

PROCEDURE NAME: APPLU
VECTORIZATION LIST

  LOOP BEGIN: (lu.f:171)
    <Unvectorized loop.>
    *** Vectorization obstructive procedure reference. : TIMER_READ (lu.f:171)
  LOOP END

  LOOP BEGIN: (lu.f:179)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (lu.f:179)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:28:59 2021
FILE NAME: lu.f

PROCEDURE NAME: APPLU
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 49 [s0-s12 s15-s16 s18-s35 s48-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1792 bytes
      Register spill area      :   64 bytes
      Parameter area           :  248 bytes
      Register save area       :  176 bytes
      User data area           : 1296 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (lu.f:171)
    *** Estimated execution cycle                       : 13
  LOOP END

  LOOP BEGIN: (lu.f:179)
    *** Estimated execution cycle                       : 147
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Across calls                                : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 8
            Across calls                                : 7
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 46
  LOOP END


