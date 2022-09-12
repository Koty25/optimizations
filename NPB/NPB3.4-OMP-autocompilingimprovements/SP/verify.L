NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
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


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    35: inl(1212): Source for routine not found.: ERROR_NORM
    36: inl(1212): Source for routine not found.: COMPUTE_RHS
    38: inl(1212): Source for routine not found.: RHS_NORM
    40: opt(1418): Constant-length loop is expanded.
    47: opt(1418): Constant-length loop is expanded.
   303: opt(1418): Constant-length loop is expanded.
   340: opt(3008): Reference within a conditional branch moved outside loop.
   340: vec( 103): Unvectorized loop.
   340: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   340: vec( 180): I/O statement obstructs vectorization.
   342: opt(1118): This I/O statement inhibits optimization of loop.
   343: inl(1212): Source for routine not found.: IEEE_ARITHMETIC::IEEE_IS_NAN
   360: opt(3008): Reference within a conditional branch moved outside loop.
   360: vec( 103): Unvectorized loop.
   360: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   360: vec( 180): I/O statement obstructs vectorization.
   362: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
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
    16:                   use sp_data
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
    44:                   class = 'U'
    45:                   verified = .true.
    46:           
    47: *------>          do m = 1,5
    48: |                    xcrref(m) = 1.0
    49: |                    xceref(m) = 1.0
    50: *------           end do
    51:           
    52:           c---------------------------------------------------------------------
    53:           c    reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
    54:           c---------------------------------------------------------------------
    55:                   if ( (grid_points(1)  .eq. 12     ) .and.
    56:                >       (grid_points(2)  .eq. 12     ) .and.
    57:                >       (grid_points(3)  .eq. 12     ) .and.
    58:                >       (no_time_steps   .eq. 100    ))  then
    59:           
    60:                      class = 'S'
    61:                      dtref = 1.5d-2
    62:           
    63:           c---------------------------------------------------------------------
    64:           c    Reference values of RMS-norms of residual.
    65:           c---------------------------------------------------------------------
    66:                      xcrref(1) = 2.7470315451339479d-02
    67:                      xcrref(2) = 1.0360746705285417d-02
    68:                      xcrref(3) = 1.6235745065095532d-02
    69:                      xcrref(4) = 1.5840557224455615d-02
    70:                      xcrref(5) = 3.4849040609362460d-02
    71:           
    72:           c---------------------------------------------------------------------
    73:           c    Reference values of RMS-norms of solution error.
    74:           c---------------------------------------------------------------------
    75:                      xceref(1) = 2.7289258557377227d-05
    76:                      xceref(2) = 1.0364446640837285d-05
    77:                      xceref(3) = 1.6154798287166471d-05
    78:                      xceref(4) = 1.5750704994480102d-05
    79:                      xceref(5) = 3.4177666183390531d-05
    80:           
    81:           
    82:           c---------------------------------------------------------------------
    83:           c    reference data for 36X36X36 grids after 400 time steps, with DT = 1.5d-03
    84:           c---------------------------------------------------------------------
    85:                   elseif ( (grid_points(1) .eq. 36) .and.
    86:                >           (grid_points(2) .eq. 36) .and.
    87:                >           (grid_points(3) .eq. 36) .and.
    88:                >           (no_time_steps . eq. 400) ) then
    89:           
    90:                      class = 'W'
    91:                      dtref = 1.5d-3
    92:           
    93:           c---------------------------------------------------------------------
    94:           c    Reference values of RMS-norms of residual.
    95:           c---------------------------------------------------------------------
    96:                      xcrref(1) = 0.1893253733584d-02
    97:                      xcrref(2) = 0.1717075447775d-03
    98:                      xcrref(3) = 0.2778153350936d-03
    99:                      xcrref(4) = 0.2887475409984d-03
   100:                      xcrref(5) = 0.3143611161242d-02
   101:           
   102:           c---------------------------------------------------------------------
   103:           c    Reference values of RMS-norms of solution error.
   104:           c---------------------------------------------------------------------
   105:                      xceref(1) = 0.7542088599534d-04
   106:                      xceref(2) = 0.6512852253086d-05
   107:                      xceref(3) = 0.1049092285688d-04
   108:                      xceref(4) = 0.1128838671535d-04
   109:                      xceref(5) = 0.1212845639773d-03
   110:           
   111:           c---------------------------------------------------------------------
   112:           c    reference data for 64X64X64 grids after 400 time steps, with DT = 1.5d-03
   113:           c---------------------------------------------------------------------
   114:                   elseif ( (grid_points(1) .eq. 64) .and.
   115:                >           (grid_points(2) .eq. 64) .and.
   116:                >           (grid_points(3) .eq. 64) .and.
   117:                >           (no_time_steps . eq. 400) ) then
   118:           
   119:                      class = 'A'
   120:                      dtref = 1.5d-3
   121:           
   122:           c---------------------------------------------------------------------
   123:           c    Reference values of RMS-norms of residual.
   124:           c---------------------------------------------------------------------
   125:                      xcrref(1) = 2.4799822399300195d0
   126:                      xcrref(2) = 1.1276337964368832d0
   127:                      xcrref(3) = 1.5028977888770491d0
   128:                      xcrref(4) = 1.4217816211695179d0
   129:                      xcrref(5) = 2.1292113035138280d0
   130:           
   131:           c---------------------------------------------------------------------
   132:           c    Reference values of RMS-norms of solution error.
   133:           c---------------------------------------------------------------------
   134:                      xceref(1) = 1.0900140297820550d-04
   135:                      xceref(2) = 3.7343951769282091d-05
   136:                      xceref(3) = 5.0092785406541633d-05
   137:                      xceref(4) = 4.7671093939528255d-05
   138:                      xceref(5) = 1.3621613399213001d-04
   139:           
   140:           c---------------------------------------------------------------------
   141:           c    reference data for 102X102X102 grids after 400 time steps,
   142:           c    with DT = 1.0d-03
   143:           c---------------------------------------------------------------------
   144:                   elseif ( (grid_points(1) .eq. 102) .and.
   145:                >           (grid_points(2) .eq. 102) .and.
   146:                >           (grid_points(3) .eq. 102) .and.
   147:                >           (no_time_steps . eq. 400) ) then
   148:           
   149:                      class = 'B'
   150:                      dtref = 1.0d-3
   151:           
   152:           c---------------------------------------------------------------------
   153:           c    Reference values of RMS-norms of residual.
   154:           c---------------------------------------------------------------------
   155:                      xcrref(1) = 0.6903293579998d+02
   156:                      xcrref(2) = 0.3095134488084d+02
   157:                      xcrref(3) = 0.4103336647017d+02
   158:                      xcrref(4) = 0.3864769009604d+02
   159:                      xcrref(5) = 0.5643482272596d+02
   160:           
   161:           c---------------------------------------------------------------------
   162:           c    Reference values of RMS-norms of solution error.
   163:           c---------------------------------------------------------------------
   164:                      xceref(1) = 0.9810006190188d-02
   165:                      xceref(2) = 0.1022827905670d-02
   166:                      xceref(3) = 0.1720597911692d-02
   167:                      xceref(4) = 0.1694479428231d-02
   168:                      xceref(5) = 0.1847456263981d-01
   169:           
   170:           c---------------------------------------------------------------------
   171:           c    reference data for 162X162X162 grids after 400 time steps,
   172:           c    with DT = 0.67d-03
   173:           c---------------------------------------------------------------------
   174:                   elseif ( (grid_points(1) .eq. 162) .and.
   175:                >           (grid_points(2) .eq. 162) .and.
   176:                >           (grid_points(3) .eq. 162) .and.
   177:                >           (no_time_steps . eq. 400) ) then
   178:           
   179:                      class = 'C'
   180:                      dtref = 0.67d-3
   181:           
   182:           c---------------------------------------------------------------------
   183:           c    Reference values of RMS-norms of residual.
   184:           c---------------------------------------------------------------------
   185:                      xcrref(1) = 0.5881691581829d+03
   186:                      xcrref(2) = 0.2454417603569d+03
   187:                      xcrref(3) = 0.3293829191851d+03
   188:                      xcrref(4) = 0.3081924971891d+03
   189:                      xcrref(5) = 0.4597223799176d+03
   190:           
   191:           c---------------------------------------------------------------------
   192:           c    Reference values of RMS-norms of solution error.
   193:           c---------------------------------------------------------------------
   194:                      xceref(1) = 0.2598120500183d+00
   195:                      xceref(2) = 0.2590888922315d-01
   196:                      xceref(3) = 0.5132886416320d-01
   197:                      xceref(4) = 0.4806073419454d-01
   198:                      xceref(5) = 0.5483377491301d+00
   199:           
   200:           c---------------------------------------------------------------------
   201:           c    reference data for 408X408X408 grids after 500 time steps,
   202:           c    with DT = 0.3d-03
   203:           c---------------------------------------------------------------------
   204:                   elseif ( (grid_points(1) .eq. 408) .and.
   205:                >           (grid_points(2) .eq. 408) .and.
   206:                >           (grid_points(3) .eq. 408) .and.
   207:                >           (no_time_steps . eq. 500) ) then
   208:           
   209:                      class = 'D'
   210:                      dtref = 0.30d-3
   211:           
   212:           c---------------------------------------------------------------------
   213:           c    Reference values of RMS-norms of residual.
   214:           c---------------------------------------------------------------------
   215:                      xcrref(1) = 0.1044696216887d+05
   216:                      xcrref(2) = 0.3204427762578d+04
   217:                      xcrref(3) = 0.4648680733032d+04
   218:                      xcrref(4) = 0.4238923283697d+04
   219:                      xcrref(5) = 0.7588412036136d+04
   220:           
   221:           c---------------------------------------------------------------------
   222:           c    Reference values of RMS-norms of solution error.
   223:           c---------------------------------------------------------------------
   224:                      xceref(1) = 0.5089471423669d+01
   225:                      xceref(2) = 0.5323514855894d+00
   226:                      xceref(3) = 0.1187051008971d+01
   227:                      xceref(4) = 0.1083734951938d+01
   228:                      xceref(5) = 0.1164108338568d+02
   229:           
   230:           c---------------------------------------------------------------------
   231:           c    reference data for 1020X1020X1020 grids after 500 time steps,
   232:           c    with DT = 0.1d-03
   233:           c---------------------------------------------------------------------
   234:                   elseif ( (grid_points(1) .eq. 1020) .and.
   235:                >           (grid_points(2) .eq. 1020) .and.
   236:                >           (grid_points(3) .eq. 1020) .and.
   237:                >           (no_time_steps . eq. 500) ) then
   238:           
   239:                      class = 'E'
   240:                      dtref = 0.10d-3
   241:           
   242:           c---------------------------------------------------------------------
   243:           c    Reference values of RMS-norms of residual.
   244:           c---------------------------------------------------------------------
   245:                      xcrref(1) = 0.6255387422609d+05
   246:                      xcrref(2) = 0.1495317020012d+05
   247:                      xcrref(3) = 0.2347595750586d+05
   248:                      xcrref(4) = 0.2091099783534d+05
   249:                      xcrref(5) = 0.4770412841218d+05
   250:           
   251:           c---------------------------------------------------------------------
   252:           c    Reference values of RMS-norms of solution error.
   253:           c---------------------------------------------------------------------
   254:                      xceref(1) = 0.6742735164909d+02
   255:                      xceref(2) = 0.5390656036938d+01
   256:                      xceref(3) = 0.1680647196477d+02
   257:                      xceref(4) = 0.1536963126457d+02
   258:                      xceref(5) = 0.1575330146156d+03
   259:           
   260:           c---------------------------------------------------------------------
   261:           c    reference data for 2560X2560X2560 grids after 500 time steps,
   262:           c    with DT = 0.15d-04
   263:           c---------------------------------------------------------------------
   264:                   elseif ( (grid_points(1) .eq. 2560) .and.
   265:                >           (grid_points(2) .eq. 2560) .and.
   266:                >           (grid_points(3) .eq. 2560) .and.
   267:                >           (no_time_steps . eq. 500) ) then
   268:           
   269:                      class = 'F'
   270:                      dtref = 0.15d-4
   271:           
   272:           c---------------------------------------------------------------------
   273:           c    Reference values of RMS-norms of residual.
   274:           c---------------------------------------------------------------------
   275:                      xcrref(1) = 0.9281628449462d+05
   276:                      xcrref(2) = 0.2230152287675d+05
   277:                      xcrref(3) = 0.3493102358632d+05
   278:                      xcrref(4) = 0.3114096186689d+05
   279:                      xcrref(5) = 0.7424426448298d+05
   280:           
   281:           c---------------------------------------------------------------------
   282:           c    Reference values of RMS-norms of solution error.
   283:           c---------------------------------------------------------------------
   284:                      xceref(1) = 0.2683717702444d+03
   285:                      xceref(2) = 0.2030647554028d+02
   286:                      xceref(3) = 0.6734864248234d+02
   287:                      xceref(4) = 0.5947451301640d+02
   288:                      xceref(5) = 0.5417636652565d+03
   289:           
   290:           
   291:                   else
   292:                      verified = .false.
   293:                   endif
   294:           
   295:           c---------------------------------------------------------------------
   296:           c    verification test for residuals if gridsize is one of
   297:           c    the defined grid sizes above (class .ne. 'U')
   298:           c---------------------------------------------------------------------
   299:           
   300:           c---------------------------------------------------------------------
   301:           c    Compute the difference of solution values and the known reference values.
   302:           c---------------------------------------------------------------------
   303: *------>          do m = 1, 5
   304: |          
   305: |                    xcrdif(m) = dabs((xcr(m)-xcrref(m))/xcrref(m))
   306: |                    xcedif(m) = dabs((xce(m)-xceref(m))/xceref(m))
   307: |          
   308: *------           enddo
   309:           
   310:           c---------------------------------------------------------------------
   311:           c    Output the comparison of computed results to known cases.
   312:           c---------------------------------------------------------------------
   313:           
   314:                   if (class .ne. 'U') then
   315:                      write(*, 1990) class
   316:            1990      format(' Verification being performed for class ', a)
   317:                      write (*,2000) epsilon
   318:            2000      format(' accuracy setting for epsilon = ', E20.13)
   319:                      verified = (dabs(dt-dtref) .le. epsilon)
   320:                      if (.not.verified) then
   321:                         class = 'U'
   322:                         write (*,1000) dtref
   323:            1000         format(' DT does not match the reference value of ',
   324:                >                 E15.8)
   325:                      endif
   326:                   else
   327:                      write(*, 1995)
   328:            1995      format(' Unknown class')
   329:                   endif
   330:           
   331:           
   332:                   if (class .ne. 'U') then
   333:                      write (*, 2001)
   334:                   else
   335:                      write (*, 2005)
   336:                   endif
   337:           
   338:            2001   format(' Comparison of RMS-norms of residual')
   339:            2005   format(' RMS-norms of residual')
   340: +------>          do m = 1, 5
   341: |                    if (class .eq. 'U') then
   342: |                       write(*, 2015) m, xcr(m)
   343: |                    else if ((.not.ieee_is_nan(xcrdif(m))) .and.
   344: |              &              xcrdif(m) .le. epsilon) then
   345: |                       write (*,2011) m,xcr(m),xcrref(m),xcrdif(m)
   346: |                    else
   347: |                       verified = .false.
   348: |                       write (*,2010) m,xcr(m),xcrref(m),xcrdif(m)
   349: |                    endif
   350: +------           enddo
   351:           
   352:                   if (class .ne. 'U') then
   353:                      write (*,2002)
   354:                   else
   355:                      write (*,2006)
   356:                   endif
   357:            2002   format(' Comparison of RMS-norms of solution error')
   358:            2006   format(' RMS-norms of solution error')
   359:            
   360: +------>          do m = 1, 5
   361: |                    if (class .eq. 'U') then
   362: |                       write(*, 2015) m, xce(m)
   363: |                    else if ((.not.ieee_is_nan(xcedif(m))) .and.
   364: |              &              xcedif(m) .le. epsilon) then
   365: |                       write (*,2011) m,xce(m),xceref(m),xcedif(m)
   366: |                    else
   367: |                       verified = .false.
   368: |                       write (*,2010) m,xce(m),xceref(m),xcedif(m)
   369: |                    endif
   370: +------           enddo
   371:            
   372:            2010   format(' FAILURE: ', i2, E20.13, E20.13, E20.13)
   373:            2011   format('          ', i2, E20.13, E20.13, E20.13)
   374:            2015   format('          ', i2, E20.13)
   375:            
   376:                   if (class .eq. 'U') then
   377:                      write(*, 2022)
   378:                      write(*, 2023)
   379:            2022      format(' No reference values provided')
   380:            2023      format(' No verification performed')
   381:                   else if (verified) then
   382:                      write(*, 2020)
   383:            2020      format(' Verification Successful')
   384:                   else
   385:                      write(*, 2021)
   386:            2021      format(' Verification failed')
   387:                   endif
   388:           
   389:                   return
   390:           
   391:           
   392:                   end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
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
  -> NOINLINE: IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:343)
     *** Source for routine not found.
  -> NOINLINE: IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:363)
     *** Source for routine not found.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
VECTORIZATION LIST

  LOOP BEGIN: (verify.f:340)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (verify.f:340)
    *** Vectorization obstructive procedure reference. : IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:340)
  LOOP END

  LOOP BEGIN: (verify.f:360)
    <Unvectorized loop.>
    *** I/O statement obstructs vectorization. (verify.f:360)
    *** Vectorization obstructive procedure reference. : IEEE_ARITHMETIC::IEEE_IS_NAN (verify.f:360)
  LOOP END


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 55 [s0-s12 s15-s16 s18-s35 s42-s63]
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

  LOOP BEGIN: (verify.f:340)
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

  LOOP BEGIN: (verify.f:360)
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


