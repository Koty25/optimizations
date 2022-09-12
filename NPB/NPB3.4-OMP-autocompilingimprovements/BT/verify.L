NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME : verify.f

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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    35: inl(1212): Source for routine not found.: ERROR_NORM
    36: inl(1212): Source for routine not found.: COMPUTE_RHS
    38: inl(1212): Source for routine not found.: RHS_NORM
    40: opt(1418): Constant-length loop is expanded.
    48: opt(1418): Constant-length loop is expanded.
   304: opt(1418): Constant-length loop is expanded.
   341: opt(3008): Reference within a conditional branch moved outside loop.
   341: vec( 103): Unvectorized loop.
   341: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   341: vec( 180): I/O statement obstructs vectorization.
   343: opt(1118): This I/O statement inhibits optimization of loop.
   344: inl(1212): Source for routine not found.: IEEE_ARITHMETIC::IEEE_IS_NAN
   361: opt(3008): Reference within a conditional branch moved outside loop.
   361: vec( 103): Unvectorized loop.
   361: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   361: vec( 180): I/O statement obstructs vectorization.
   363: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
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

     5:                   subroutine verify(no_time_steps, class, verified)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c  verification routine
    12:           c---------------------------------------------------------------------
    13:           
    14:                   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    15:           
    16:                   use bt_data
    17:           
    18:                   implicit none
    19:           
    20:                   double precision xcrref(5),xceref(5),xcrdif(5),xcedif(5),
    21:                >                   epsilon, xce(5), xcr(5), dtref
    22:                   integer m, no_time_steps
    23:                   character class
    24:                   logical verified
    25:           
    26:           c---------------------------------------------------------------------
    27:           c   tolerance level
    28:           c---------------------------------------------------------------------
    29:                   epsilon = 1.0d-08
    30:           
    31:           
    32:           c---------------------------------------------------------------------
    33:           c   compute the error norm and the residual norm, and exit if not printing
    34:           c---------------------------------------------------------------------
    35:                   call error_norm(xce)
    36:                   call compute_rhs
    37:           
    38:                   call rhs_norm(xcr)
    39:           
    40: *------>          do m = 1, 5
    41: |                    xcr(m) = xcr(m) / dt
    42: *------           enddo
    43:           
    44:           
    45:                   class = 'U'
    46:                   verified = .true.
    47:           
    48: *------>          do m = 1,5
    49: |                    xcrref(m) = 1.0
    50: |                    xceref(m) = 1.0
    51: *------           end do
    52:           
    53:           c---------------------------------------------------------------------
    54:           c    reference data for 12X12X12 grids after 60 time steps, with DT = 1.0d-02
    55:           c---------------------------------------------------------------------
    56:                   if ( (grid_points(1)  .eq. 12     ) .and.
    57:                >       (grid_points(2)  .eq. 12     ) .and.
    58:                >       (grid_points(3)  .eq. 12     ) .and.
    59:                >       (no_time_steps   .eq. 60    ))  then
    60:           
    61:                      class = 'S'
    62:                      dtref = 1.0d-2
    63:           
    64:           c---------------------------------------------------------------------
    65:           c  Reference values of RMS-norms of residual.
    66:           c---------------------------------------------------------------------
    67:                    xcrref(1) = 1.7034283709541311d-01
    68:                    xcrref(2) = 1.2975252070034097d-02
    69:                    xcrref(3) = 3.2527926989486055d-02
    70:                    xcrref(4) = 2.6436421275166801d-02
    71:                    xcrref(5) = 1.9211784131744430d-01
    72:           
    73:           c---------------------------------------------------------------------
    74:           c  Reference values of RMS-norms of solution error.
    75:           c---------------------------------------------------------------------
    76:                    xceref(1) = 4.9976913345811579d-04
    77:                    xceref(2) = 4.5195666782961927d-05
    78:                    xceref(3) = 7.3973765172921357d-05
    79:                    xceref(4) = 7.3821238632439731d-05
    80:                    xceref(5) = 8.9269630987491446d-04
    81:           
    82:           c---------------------------------------------------------------------
    83:           c    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
    84:           c---------------------------------------------------------------------
    85:                   elseif ( (grid_points(1) .eq. 24) .and.
    86:                >           (grid_points(2) .eq. 24) .and.
    87:                >           (grid_points(3) .eq. 24) .and.
    88:                >           (no_time_steps . eq. 200) ) then
    89:           
    90:                      class = 'W'
    91:                      dtref = 0.8d-3
    92:           c---------------------------------------------------------------------
    93:           c  Reference values of RMS-norms of residual.
    94:           c---------------------------------------------------------------------
    95:                      xcrref(1) = 0.1125590409344d+03
    96:                      xcrref(2) = 0.1180007595731d+02
    97:                      xcrref(3) = 0.2710329767846d+02
    98:                      xcrref(4) = 0.2469174937669d+02
    99:                      xcrref(5) = 0.2638427874317d+03
   100:           
   101:           c---------------------------------------------------------------------
   102:           c  Reference values of RMS-norms of solution error.
   103:           c---------------------------------------------------------------------
   104:                      xceref(1) = 0.4419655736008d+01
   105:                      xceref(2) = 0.4638531260002d+00
   106:                      xceref(3) = 0.1011551749967d+01
   107:                      xceref(4) = 0.9235878729944d+00
   108:                      xceref(5) = 0.1018045837718d+02
   109:           
   110:           
   111:           c---------------------------------------------------------------------
   112:           c    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
   113:           c---------------------------------------------------------------------
   114:                   elseif ( (grid_points(1) .eq. 64) .and.
   115:                >           (grid_points(2) .eq. 64) .and.
   116:                >           (grid_points(3) .eq. 64) .and.
   117:                >           (no_time_steps . eq. 200) ) then
   118:           
   119:                      class = 'A'
   120:                      dtref = 0.8d-3
   121:           c---------------------------------------------------------------------
   122:           c  Reference values of RMS-norms of residual.
   123:           c---------------------------------------------------------------------
   124:                    xcrref(1) = 1.0806346714637264d+02
   125:                    xcrref(2) = 1.1319730901220813d+01
   126:                    xcrref(3) = 2.5974354511582465d+01
   127:                    xcrref(4) = 2.3665622544678910d+01
   128:                    xcrref(5) = 2.5278963211748344d+02
   129:           
   130:           c---------------------------------------------------------------------
   131:           c  Reference values of RMS-norms of solution error.
   132:           c---------------------------------------------------------------------
   133:                    xceref(1) = 4.2348416040525025d+00
   134:                    xceref(2) = 4.4390282496995698d-01
   135:                    xceref(3) = 9.6692480136345650d-01
   136:                    xceref(4) = 8.8302063039765474d-01
   137:                    xceref(5) = 9.7379901770829278d+00
   138:           
   139:           c---------------------------------------------------------------------
   140:           c    reference data for 102X102X102 grids after 200 time steps,
   141:           c    with DT = 3.0d-04
   142:           c---------------------------------------------------------------------
   143:                   elseif ( (grid_points(1) .eq. 102) .and.
   144:                >           (grid_points(2) .eq. 102) .and.
   145:                >           (grid_points(3) .eq. 102) .and.
   146:                >           (no_time_steps . eq. 200) ) then
   147:           
   148:                      class = 'B'
   149:                      dtref = 3.0d-4
   150:           
   151:           c---------------------------------------------------------------------
   152:           c  Reference values of RMS-norms of residual.
   153:           c---------------------------------------------------------------------
   154:                    xcrref(1) = 1.4233597229287254d+03
   155:                    xcrref(2) = 9.9330522590150238d+01
   156:                    xcrref(3) = 3.5646025644535285d+02
   157:                    xcrref(4) = 3.2485447959084092d+02
   158:                    xcrref(5) = 3.2707541254659363d+03
   159:           
   160:           c---------------------------------------------------------------------
   161:           c  Reference values of RMS-norms of solution error.
   162:           c---------------------------------------------------------------------
   163:                    xceref(1) = 5.2969847140936856d+01
   164:                    xceref(2) = 4.4632896115670668d+00
   165:                    xceref(3) = 1.3122573342210174d+01
   166:                    xceref(4) = 1.2006925323559144d+01
   167:                    xceref(5) = 1.2459576151035986d+02
   168:           
   169:           c---------------------------------------------------------------------
   170:           c    reference data for 162X162X162 grids after 200 time steps,
   171:           c    with DT = 1.0d-04
   172:           c---------------------------------------------------------------------
   173:                   elseif ( (grid_points(1) .eq. 162) .and.
   174:                >           (grid_points(2) .eq. 162) .and.
   175:                >           (grid_points(3) .eq. 162) .and.
   176:                >           (no_time_steps . eq. 200) ) then
   177:           
   178:                      class = 'C'
   179:                      dtref = 1.0d-4
   180:           
   181:           c---------------------------------------------------------------------
   182:           c  Reference values of RMS-norms of residual.
   183:           c---------------------------------------------------------------------
   184:                    xcrref(1) = 0.62398116551764615d+04
   185:                    xcrref(2) = 0.50793239190423964d+03
   186:                    xcrref(3) = 0.15423530093013596d+04
   187:                    xcrref(4) = 0.13302387929291190d+04
   188:                    xcrref(5) = 0.11604087428436455d+05
   189:           
   190:           c---------------------------------------------------------------------
   191:           c  Reference values of RMS-norms of solution error.
   192:           c---------------------------------------------------------------------
   193:                    xceref(1) = 0.16462008369091265d+03
   194:                    xceref(2) = 0.11497107903824313d+02
   195:                    xceref(3) = 0.41207446207461508d+02
   196:                    xceref(4) = 0.37087651059694167d+02
   197:                    xceref(5) = 0.36211053051841265d+03
   198:           
   199:           c---------------------------------------------------------------------
   200:           c    reference data for 408x408x408 grids after 250 time steps,
   201:           c    with DT = 0.2d-04
   202:           c---------------------------------------------------------------------
   203:                   elseif ( (grid_points(1) .eq. 408) .and.
   204:                >           (grid_points(2) .eq. 408) .and.
   205:                >           (grid_points(3) .eq. 408) .and.
   206:                >           (no_time_steps . eq. 250) ) then
   207:           
   208:                      class = 'D'
   209:                      dtref = 0.2d-4
   210:           
   211:           c---------------------------------------------------------------------
   212:           c  Reference values of RMS-norms of residual.
   213:           c---------------------------------------------------------------------
   214:                    xcrref(1) = 0.2533188551738d+05
   215:                    xcrref(2) = 0.2346393716980d+04
   216:                    xcrref(3) = 0.6294554366904d+04
   217:                    xcrref(4) = 0.5352565376030d+04
   218:                    xcrref(5) = 0.3905864038618d+05
   219:           
   220:           c---------------------------------------------------------------------
   221:           c  Reference values of RMS-norms of solution error.
   222:           c---------------------------------------------------------------------
   223:           
   224:                    xceref(1) = 0.3100009377557d+03
   225:                    xceref(2) = 0.2424086324913d+02
   226:                    xceref(3) = 0.7782212022645d+02
   227:                    xceref(4) = 0.6835623860116d+02
   228:                    xceref(5) = 0.6065737200368d+03
   229:           
   230:           c---------------------------------------------------------------------
   231:           c    reference data for 1020x1020x1020 grids after 250 time steps,
   232:           c    with DT = 0.4d-05
   233:           c---------------------------------------------------------------------
   234:                   elseif ( (grid_points(1) .eq. 1020) .and.
   235:                >           (grid_points(2) .eq. 1020) .and.
   236:                >           (grid_points(3) .eq. 1020) .and.
   237:                >           (no_time_steps . eq. 250) ) then
   238:           
   239:                      class = 'E'
   240:                      dtref = 0.4d-5
   241:           
   242:           c---------------------------------------------------------------------
   243:           c  Reference values of RMS-norms of residual.
   244:           c---------------------------------------------------------------------
   245:                    xcrref(1) = 0.9795372484517d+05
   246:                    xcrref(2) = 0.9739814511521d+04
   247:                    xcrref(3) = 0.2467606342965d+05
   248:                    xcrref(4) = 0.2092419572860d+05
   249:                    xcrref(5) = 0.1392138856939d+06
   250:           
   251:           c---------------------------------------------------------------------
   252:           c  Reference values of RMS-norms of solution error.
   253:           c---------------------------------------------------------------------
   254:           
   255:                    xceref(1) = 0.4327562208414d+03
   256:                    xceref(2) = 0.3699051964887d+02
   257:                    xceref(3) = 0.1089845040954d+03
   258:                    xceref(4) = 0.9462517622043d+02
   259:                    xceref(5) = 0.7765512765309d+03
   260:           
   261:           c---------------------------------------------------------------------
   262:           c    reference data for 2560x2560x2560 grids after 250 time steps,
   263:           c    with DT = 0.6d-06
   264:           c---------------------------------------------------------------------
   265:                   elseif ( (grid_points(1) .eq. 2560) .and.
   266:                >           (grid_points(2) .eq. 2560) .and.
   267:                >           (grid_points(3) .eq. 2560) .and.
   268:                >           (no_time_steps . eq. 250) ) then
   269:           
   270:                      class = 'F'
   271:                      dtref = 0.6d-6
   272:           
   273:           c---------------------------------------------------------------------
   274:           c  Reference values of RMS-norms of residual.
   275:           c---------------------------------------------------------------------
   276:                    xcrref(1) = 0.4240735175585d+06
   277:                    xcrref(2) = 0.4348701133212d+05
   278:                    xcrref(3) = 0.1078114688845d+06
   279:                    xcrref(4) = 0.9142160938556d+05
   280:                    xcrref(5) = 0.5879842143431d+06
   281:           
   282:           c---------------------------------------------------------------------
   283:           c  Reference values of RMS-norms of solution error.
   284:           c---------------------------------------------------------------------
   285:           
   286:                    xceref(1) = 0.5095577042351d+03
   287:                    xceref(2) = 0.4557065541652d+02
   288:                    xceref(3) = 0.1286632140581d+03
   289:                    xceref(4) = 0.1111419378722d+03
   290:                    xceref(5) = 0.8720011709356d+03
   291:           
   292:                   else
   293:                      verified = .false.
   294:                   endif
   295:           
   296:           c---------------------------------------------------------------------
   297:           c    verification test for residuals if gridsize is one of
   298:           c    the defined grid sizes above (class .ne. 'U')
   299:           c---------------------------------------------------------------------
   300:           
   301:           c---------------------------------------------------------------------
   302:           c    Compute the difference of solution values and the known reference values.
   303:           c---------------------------------------------------------------------
   304: *------>          do m = 1, 5
   305: |          
   306: |                    xcrdif(m) = dabs((xcr(m)-xcrref(m))/xcrref(m))
   307: |                    xcedif(m) = dabs((xce(m)-xceref(m))/xceref(m))
   308: |          
   309: *------           enddo
   310:           
   311:           c---------------------------------------------------------------------
   312:           c    Output the comparison of computed results to known cases.
   313:           c---------------------------------------------------------------------
   314:           
   315:                   if (class .ne. 'U') then
   316:                      write(*, 1990) class
   317:            1990      format(' Verification being performed for class ', a)
   318:                      write (*,2000) epsilon
   319:            2000      format(' accuracy setting for epsilon = ', E20.13)
   320:                      verified = (dabs(dt-dtref) .le. epsilon)
   321:                      if (.not.verified) then
   322:                         class = 'U'
   323:                         write (*,1000) dtref
   324:            1000         format(' DT does not match the reference value of ',
   325:                >                 E15.8)
   326:                      endif
   327:                   else
   328:                      write(*, 1995)
   329:            1995      format(' Unknown class')
   330:                   endif
   331:           
   332:           
   333:                   if (class .ne. 'U') then
   334:                      write (*, 2001)
   335:                   else
   336:                      write (*, 2005)
   337:                   endif
   338:           
   339:            2001   format(' Comparison of RMS-norms of residual')
   340:            2005   format(' RMS-norms of residual')
   341: +------>          do m = 1, 5
   342: |                    if (class .eq. 'U') then
   343: |                       write(*, 2015) m, xcr(m)
   344: |                    else if ((.not.ieee_is_nan(xcrdif(m))) .and.
   345: |              &              xcrdif(m) .le. epsilon) then
   346: |                       write (*,2011) m,xcr(m),xcrref(m),xcrdif(m)
   347: |                    else
   348: |                       verified = .false.
   349: |                       write (*,2010) m,xcr(m),xcrref(m),xcrdif(m)
   350: |                    endif
   351: +------           enddo
   352:           
   353:                   if (class .ne. 'U') then
   354:                      write (*,2002)
   355:                   else
   356:                      write (*,2006)
   357:                   endif
   358:            2002   format(' Comparison of RMS-norms of solution error')
   359:            2006   format(' RMS-norms of solution error')
   360:            
   361: +------>          do m = 1, 5
   362: |                    if (class .eq. 'U') then
   363: |                       write(*, 2015) m, xce(m)
   364: |                    else if ((.not.ieee_is_nan(xcedif(m))) .and.
   365: |              &              xcedif(m) .le. epsilon) then
   366: |                       write (*,2011) m,xce(m),xceref(m),xcedif(m)
   367: |                    else
   368: |                       verified = .false.
   369: |                       write (*,2010) m,xce(m),xceref(m),xcedif(m)
   370: |                    endif
   371: +------           enddo
   372:            
   373:            2010   format(' FAILURE: ', i2, E20.13, E20.13, E20.13)
   374:            2011   format('          ', i2, E20.13, E20.13, E20.13)
   375:            2015   format('          ', i2, E20.13)
   376:            
   377:                   if (class .eq. 'U') then
   378:                      write(*, 2022)
   379:                      write(*, 2023)
   380:            2022      format(' No reference values provided')
   381:            2023      format(' No verification performed')
   382:                   else if (verified) then
   383:                      write(*, 2020)
   384:            2020      format(' Verification Successful')
   385:                   else
   386:                      write(*, 2021)
   387:            2021      format(' Verification failed')
   388:                   endif
   389:           
   390:                   return
   391:           
   392:           
   393:                   end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
INLINE LIST

  ROOT: VERIFY (verify.f:5)
  -> NOINLINE: ERROR_NORM (verify.f:35)
     *** Source for routine not found.
  -> NOINLINE: COMPUTE_RHS (verify.f:36)
     *** Source for routine not found.
  -> NOINLINE: RHS_NORM (verify.f:38)
     *** Source for routine not found.
  -> NOINLINE: IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:344)
     *** Source for routine not found.
  -> NOINLINE: IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:364)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
VECTORIZATION LIST

  LOOP BEGIN: (verify.f:341)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (verify.f:341)
    *** Vectorization obstructive procedure reference. : IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:341)
  LOOP END

  LOOP BEGIN: (verify.f:361)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (verify.f:361)
    *** Vectorization obstructive procedure reference. : IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:361)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:47 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 54 [s0-s12 s15-s16 s18-s35 s43-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1728 bytes
      Register spill area      :  104 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1368 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (verify.f:341)
    *** Estimated execution cycle                       : 61
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 5
            Across calls                                : 4
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 53
  LOOP END

  LOOP BEGIN: (verify.f:361)
    *** Estimated execution cycle                       : 61
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 2
            Across calls                                : 2
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 5
            Across calls                                : 4
            Over basic blocks                           : 1
    *** The number of SCALAR REGISTER TRANSFER          : 53
  LOOP END


