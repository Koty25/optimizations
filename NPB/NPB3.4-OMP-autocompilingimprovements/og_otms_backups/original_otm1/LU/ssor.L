NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME : ssor.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: ssor.f

PROCEDURE NAME: SSOR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    36: vec( 103): Unvectorized loop.
    36: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
    53: vec( 103): Unvectorized loop.
    53: vec( 110): Vectorization obstructive procedure reference.: TIMER_CLEAR
    61: opt(3008): Reference within a conditional branch moved outside loop.
    61: vec( 109): Vectorization obstructive statement.
    61: vec( 110): Vectorization obstructive procedure reference.: RHS
    61: vec( 110): Vectorization obstructive procedure reference.: L2NORM
    61: vec( 110): Vectorization obstructive procedure reference.: TIMER_START
    61: vec( 110): Vectorization obstructive procedure reference.: TIMER_STOP
    61: vec( 180): I/O statement obstructs vectorization.
    61: vec( 103): Unvectorized loop.
    73: par(1801): Parallel routine generated.: SSOR$1
    76: par(1805): Master-construct generated.
    81: par(1803): Parallelized by "do".
    81: opt(2495): Loop nest collapsed by "collapse" clause.
    81: vec( 101): Vectorized loop.
    83: vec( 101): Vectorized loop.
    84: opt(1418): Constant-length loop is expanded.
    91: par(1805): Master-construct generated.
    98: par(1809): Barrier synchronization.
   104: vec( 110): Vectorization obstructive procedure reference.: JACLD
   104: vec( 110): Vectorization obstructive procedure reference.: BLTS
   104: par(1803): Parallelized by "do".
   104: vec( 103): Unvectorized loop.
   126: par(1809): Barrier synchronization.
   127: par(1805): Master-construct generated.
   136: vec( 110): Vectorization obstructive procedure reference.: JACU
   136: vec( 110): Vectorization obstructive procedure reference.: BUTS
   136: par(1803): Parallelized by "do".
   136: vec( 103): Unvectorized loop.
   158: par(1809): Barrier synchronization.
   159: par(1805): Master-construct generated.
   167: par(1805): Master-construct generated.
   172: par(1803): Parallelized by "do".
   172: opt(2495): Loop nest collapsed by "collapse" clause.
   172: vec( 101): Vectorized loop.
   174: vec( 101): Vectorized loop.
   175: opt(1418): Constant-length loop is expanded.
   176: vec( 128): Fused multiply-add operation applied.
   183: par(1805): Master-construct generated.


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: ssor.f

PROCEDURE NAME: SSOR
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

     4:                 subroutine ssor(niter)
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:           c---------------------------------------------------------------------
    10:           c   to perform pseudo-time stepping SSOR iterations
    11:           c   for five nonlinear pde's.
    12:           c---------------------------------------------------------------------
    13:           
    14:                 use lu_data
    15:                 implicit none
    16:           
    17:                 integer niter
    18:           
    19:           c---------------------------------------------------------------------
    20:           c  local variables
    21:           c---------------------------------------------------------------------
    22:                 integer i, j, k, m, n
    23:                 integer istep
    24:                 double precision  tmp, tmp2
    25:                 double precision  delunm(5)
    26:           
    27:                 external timer_read
    28:                 double precision timer_read
    29:           
    30:            
    31:           c---------------------------------------------------------------------
    32:           c   begin pseudo-time stepping iterations
    33:           c---------------------------------------------------------------------
    34:                 tmp = 1.0d+00 / ( omega * ( 2.0d+00 - omega ) )
    35:           
    36: +------>        do i = 1, t_last
    37: |                  call timer_clear(i)
    38: +------         end do
    39:           
    40:           c---------------------------------------------------------------------
    41:           c   compute the steady-state residuals
    42:           c---------------------------------------------------------------------
    43:                 call rhs
    44:            
    45:           c---------------------------------------------------------------------
    46:           c   compute the L2 norms of newton iteration residuals
    47:           c---------------------------------------------------------------------
    48:                 call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
    49:                >             ist, iend, jst, jend,
    50:                >             rsd, rsdnm )
    51:           
    52:            
    53: +------>        do i = 1, t_last
    54: |                  call timer_clear(i)
    55: +------         end do
    56:                 call timer_start(1)
    57:            
    58:           c---------------------------------------------------------------------
    59:           c   the timestep loop
    60:           c---------------------------------------------------------------------
    61: +------>        do istep = 1, niter
    62: |         
    63: |                  if (mod ( istep, 20) .eq. 0 .or.
    64: |              >         istep .eq. itmax .or.
    65: |              >         istep .eq. 1) then
    66: |                     if (niter .gt. 1) write( *, 200) istep
    67: |          200        format(' Time step ', i4)
    68: |                  endif
    69: |          
    70: |         c---------------------------------------------------------------------
    71: |         c   perform SSOR iteration
    72: |         c---------------------------------------------------------------------
    73: |         !$omp parallel default(shared) private(i,j,k,m,tmp2)
    74: |         !$omp&  shared(ist,iend,jst,jend,nx,ny,nz,nx0,ny0,omega)
    75: |         
    76: |         !$omp master
    77: |                  if (timeron) call timer_start(t_rhs)
    78: |         !$omp end master
    79: |                  tmp2 = dt
    80: |         !$omp do schedule(static) collapse(2)
    81: |Y----->           do k = 2, nz - 1
    82: ||+---->              do j = jst, jend
    83: |||V--->                 do i = ist, iend
    84: ||||*-->                    do m = 1, 5
    85: |||||                          rsd(m,i,j,k) = tmp2 * rsd(m,i,j,k)
    86: ||||*--                     end do
    87: |||V---                  end do
    88: ||+----               end do
    89: |Y-----            end do
    90: |         !$omp end do nowait
    91: |         !$omp master
    92: |                  if (timeron) call timer_stop(t_rhs)
    93: |         
    94: |                  if (timeron) call timer_start(t_blts)
    95: |         !$omp end master
    96: |         
    97: |                  call sync_init( jend-jst )
    98: |         !$omp barrier
    99: |         
   100: |+----->           do k = 2, nz -1
   101: ||        
   102: ||                    call sync_left( isiz1, isiz2, isiz3, rsd )
   103: ||        !$omp do schedule(static)
   104: ||P---->              do j = jst, jend
   105: |||       
   106: |||       c---------------------------------------------------------------------
   107: |||       c   form the lower triangular part of the jacobian matrix
   108: |||       c---------------------------------------------------------------------
   109: |||                      call jacld(j, k)
   110: |||        
   111: |||       c---------------------------------------------------------------------
   112: |||       c   perform the lower triangular solution
   113: |||       c---------------------------------------------------------------------
   114: |||                      call blts( isiz1, isiz2, isiz3,
   115: |||            >                    nx, ny, nz,
   116: |||            >                    omega,
   117: |||            >                    rsd,
   118: |||            >                    a, b, c, d,
   119: |||            >                    ist, iend, j, k )
   120: |||       
   121: ||P----               end do
   122: ||        !$omp end do nowait
   123: ||                    call sync_right( isiz1, isiz2, isiz3, rsd )
   124: ||        
   125: |+-----            end do
   126: |         !$omp barrier
   127: |         !$omp master
   128: |                  if (timeron) call timer_stop(t_blts)
   129: |         
   130: |                  if (timeron) call timer_start(t_buts)
   131: |         !$omp end master
   132: |+----->           do k = nz - 1, 2, -1
   133: ||        
   134: ||                    call sync_left( isiz1, isiz2, isiz3, rsd )
   135: ||        !$omp do schedule(static)
   136: ||P---->              do j = jend, jst, -1
   137: |||       
   138: |||       c---------------------------------------------------------------------
   139: |||       c   form the strictly upper triangular part of the jacobian matrix
   140: |||       c---------------------------------------------------------------------
   141: |||                      call jacu(j, k)
   142: |||       
   143: |||       c---------------------------------------------------------------------
   144: |||       c   perform the upper triangular solution
   145: |||       c---------------------------------------------------------------------
   146: |||                      call buts( isiz1, isiz2, isiz3,
   147: |||            >                    nx, ny, nz,
   148: |||            >                    omega,
   149: |||            >                    rsd,
   150: |||            >                    d, a, b, c,
   151: |||            >                    ist, iend, j, k )
   152: |||       
   153: ||P----               end do
   154: ||        !$omp end do nowait
   155: ||                    call sync_right( isiz1, isiz2, isiz3, rsd )
   156: ||        
   157: |+-----            end do
   158: |         !$omp barrier
   159: |         !$omp master
   160: |                  if (timeron) call timer_stop(t_buts)
   161: |         !$omp end master
   162: |         
   163: |         c---------------------------------------------------------------------
   164: |         c   update the variables
   165: |         c---------------------------------------------------------------------
   166: |         
   167: |         !$omp master
   168: |                  if (timeron) call timer_start(t_add)
   169: |         !$omp end master
   170: |                  tmp2 = tmp
   171: |         !$omp do schedule(static) collapse(2)
   172: |Y----->           do k = 2, nz-1
   173: ||+---->              do j = jst, jend
   174: |||V--->                 do i = ist, iend
   175: ||||*-->                    do m = 1, 5
   176: |||||   F                      u( m, i, j, k ) = u( m, i, j, k )
   177: |||||          >                    + tmp2 * rsd( m, i, j, k )
   178: ||||*--                     end do
   179: |||V---                  end do
   180: ||+----               end do
   181: |Y-----            end do
   182: |         !$omp end do nowait
   183: |         !$omp master
   184: |                  if (timeron) call timer_stop(t_add)
   185: |         !$omp end master
   186: |         !$omp end parallel
   187: |          
   188: |         c---------------------------------------------------------------------
   189: |         c   compute the max-norms of newton iteration corrections
   190: |         c---------------------------------------------------------------------
   191: |                  if ( mod ( istep, inorm ) .eq. 0 ) then
   192: |                     if (timeron) call timer_start(t_l2norm)
   193: |                     call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
   194: |              >                   ist, iend, jst, jend,
   195: |              >                   rsd, delunm )
   196: |                     if (timeron) call timer_stop(t_l2norm)
   197: |         c            if ( ipr .eq. 1 ) then
   198: |         c                write (*,1006) ( delunm(m), m = 1, 5 )
   199: |         c            else if ( ipr .eq. 2 ) then
   200: |         c                write (*,'(i5,f15.6)') istep,delunm(5)
   201: |         c            end if
   202: |                  end if
   203: |          
   204: |         c---------------------------------------------------------------------
   205: |         c   compute the steady-state residuals
   206: |         c---------------------------------------------------------------------
   207: |                  call rhs
   208: |          
   209: |         c---------------------------------------------------------------------
   210: |         c   compute the max-norms of newton iteration residuals
   211: |         c---------------------------------------------------------------------
   212: |                  if ( ( mod ( istep, inorm ) .eq. 0 ) .or.
   213: |              >        ( istep .eq. itmax ) ) then
   214: |                     if (timeron) call timer_start(t_l2norm)
   215: |                     call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
   216: |              >                   ist, iend, jst, jend,
   217: |              >                   rsd, rsdnm )
   218: |                     if (timeron) call timer_stop(t_l2norm)
   219: |         c            if ( ipr .eq. 1 ) then
   220: |         c                write (*,1007) ( rsdnm(m), m = 1, 5 )
   221: |         c            end if
   222: |                  end if
   223: |         
   224: |         c---------------------------------------------------------------------
   225: |         c   check the newton-iteration residuals against the tolerance levels
   226: |         c---------------------------------------------------------------------
   227: |                  if ( ( rsdnm(1) .lt. tolrsd(1) ) .and.
   228: |              >        ( rsdnm(2) .lt. tolrsd(2) ) .and.
   229: |              >        ( rsdnm(3) .lt. tolrsd(3) ) .and.
   230: |              >        ( rsdnm(4) .lt. tolrsd(4) ) .and.
   231: |              >        ( rsdnm(5) .lt. tolrsd(5) ) ) then
   232: |         c            if (ipr .eq. 1 ) then
   233: |                        write (*,1004) istep
   234: |         c            end if
   235: |                     go to 900
   236: |                  end if
   237: |          
   238: +------         end do
   239:             900 continue
   240:            
   241:                 call timer_stop(1)
   242:                 maxtime= timer_read(1)
   243:            
   244:           
   245:           
   246:                 return
   247:            
   248:            1001 format (1x/5x,'pseudo-time SSOR iteration no.=',i4/)
   249:            1004 format (1x/1x,'convergence was achieved after ',i4,
   250:                >   ' pseudo-time steps' )
   251:            1006 format (1x/1x,'RMS-norm of SSOR-iteration correction ',
   252:                > 'for first pde  = ',1pe12.5/,
   253:                > 1x,'RMS-norm of SSOR-iteration correction ',
   254:                > 'for second pde = ',1pe12.5/,
   255:                > 1x,'RMS-norm of SSOR-iteration correction ',
   256:                > 'for third pde  = ',1pe12.5/,
   257:                > 1x,'RMS-norm of SSOR-iteration correction ',
   258:                > 'for fourth pde = ',1pe12.5/,
   259:                > 1x,'RMS-norm of SSOR-iteration correction ',
   260:                > 'for fifth pde  = ',1pe12.5)
   261:            1007 format (1x/1x,'RMS-norm of steady-state residual for ',
   262:                > 'first pde  = ',1pe12.5/,
   263:                > 1x,'RMS-norm of steady-state residual for ',
   264:                > 'second pde = ',1pe12.5/,
   265:                > 1x,'RMS-norm of steady-state residual for ',
   266:                > 'third pde  = ',1pe12.5/,
   267:                > 1x,'RMS-norm of steady-state residual for ',
   268:                > 'fourth pde = ',1pe12.5/,
   269:                > 1x,'RMS-norm of steady-state residual for ',
   270:                > 'fifth pde  = ',1pe12.5)
   271:            
   272:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: ssor.f

PROCEDURE NAME: SSOR
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (ssor.f:81)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ssor.f:81)
    *** The number of VLOAD, VSTORE. :  0,  2. (ssor.f:81)

    LOOP BEGIN: (ssor.f:83)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ssor.f:83)
      *** The number of VLOAD, VSTORE. :  5,  5. (ssor.f:83)
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:100)
    <Unvectorized loop.>

    LOOP BEGIN: (ssor.f:104)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : JACLD (ssor.f:104)
      *** Vectorization obstructive procedure reference. : BLTS (ssor.f:104)
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:132)
    <Unvectorized loop.>

    LOOP BEGIN: (ssor.f:136)
      <Unvectorized loop.>
      *** Vectorization obstructive procedure reference. : JACU (ssor.f:136)
      *** Vectorization obstructive procedure reference. : BUTS (ssor.f:136)
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:172)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (ssor.f:172)
    *** The number of VLOAD, VSTORE. :  0,  2. (ssor.f:172)

    LOOP BEGIN: (ssor.f:174)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  0,  0. (ssor.f:174)
      *** The number of VLOAD, VSTORE. : 10,  5. (ssor.f:174)
    LOOP END
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Thu Apr 15 12:29:00 2021
FILE NAME: ssor.f

PROCEDURE NAME: SSOR
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 61 [s0-s12 s15-s16 s18-s63]
      Vector registers         : 51 [v13-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1104 bytes
      Register spill area      :  336 bytes
      Parameter area           :  128 bytes
      Register save area       :  176 bytes
      User data area           :  464 bytes
      Others                   :    0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (ssor.f:81)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (ssor.f:81)
    *** Estimated execution cycle                       : 25
    *** The number of SCALAR REGISTER TRANSFER          : 5

    LOOP BEGIN: (ssor.f:83)
      *** Estimated execution cycle                     : 338
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:100)
    *** Estimated execution cycle                       : 85
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 15
            Across calls                                : 12
            Over basic blocks                           :  3
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 26
            Across calls                                :  6
            Over basic blocks                           : 20
    *** The number of SCALAR REGISTER TRANSFER          : 24

    LOOP BEGIN: (ssor.f:104)
      *** Estimated execution cycle                     : 81
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 15
              Across calls                              :  7
              Over basic blocks                         :  8
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:132)
    *** Estimated execution cycle                       : 90
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 15
            Across calls                                : 12
            Over basic blocks                           :  3
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 25
            Across calls                                :  6
            Over basic blocks                           : 19
    *** The number of SCALAR REGISTER TRANSFER          : 25

    LOOP BEGIN: (ssor.f:136)
      *** Estimated execution cycle                     : 83
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 16
              Across calls                              :  7
              Over basic blocks                         :  9
      *** The number of SCALAR REGISTER TRANSFER        : 2
    LOOP END
  LOOP END

  LOOP BEGIN: (ssor.f:172)
    *** Estimated execution cycle                       : 575
  LOOP END

  LOOP BEGIN: (ssor.f:172)
    *** Estimated execution cycle                       : 36
    *** The number of SCALAR REGISTER TRANSFER          : 6

    LOOP BEGIN: (ssor.f:174)
      *** Estimated execution cycle                     : 487
      *** The number of SCALAR REGISTER TRANSFER        : 4
    LOOP END
  LOOP END


