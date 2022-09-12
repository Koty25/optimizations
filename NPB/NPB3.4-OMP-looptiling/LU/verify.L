NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    36: opt(1418): Constant-length loop is expanded.
   332: opt(1418): Constant-length loop is expanded.
   371: vec( 103): Unvectorized loop.
   371: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   371: vec( 180): I/O statement obstructs vectorization.
   373: opt(1118): This I/O statement inhibits optimization of loop.
   374: inl(1212): Source for routine not found.: IEEE_ARITHMETIC::IEEE_IS_NAN
   391: vec( 103): Unvectorized loop.
   391: vec( 110): Vectorization obstructive procedure reference.: IEEE_ARITHMETIC::IEEE_IS_NAN
   391: vec( 180): I/O statement obstructs vectorization.
   393: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:23 2020
FILE NAME: verify.f

PROCEDURE NAME: VERIFY
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                   subroutine verify(xcr, xce, xci, class, verified)
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
    16:                   use lu_data
    17:           
    18:                   implicit none
    19:           
    20:                   double precision xcr(5), xce(5), xci
    21:                   double precision xcrref(5),xceref(5),xciref,
    22:                >                   xcrdif(5),xcedif(5),xcidif,
    23:                >                   epsilon, dtref
    24:                   integer m
    25:                   character class
    26:                   logical verified
    27:           
    28:           c---------------------------------------------------------------------
    29:           c   tolerance level
    30:           c---------------------------------------------------------------------
    31:                   epsilon = 1.0d-08
    32:           
    33:                   class = 'U'
    34:                   verified = .true.
    35:           
    36: *------>          do m = 1,5
    37: |                    xcrref(m) = 1.0
    38: |                    xceref(m) = 1.0
    39: *------           end do
    40:                   xciref = 1.0
    41:           
    42:                   if ( (nx0  .eq. 12     ) .and.
    43:                >       (ny0  .eq. 12     ) .and.
    44:                >       (nz0  .eq. 12     ) .and.
    45:                >       (itmax   .eq. 50    ))  then
    46:           
    47:                      class = 'S'
    48:                      dtref = 5.0d-1
    49:           c---------------------------------------------------------------------
    50:           c   Reference values of RMS-norms of residual, for the (12X12X12) grid,
    51:           c   after 50 time steps, with  DT = 5.0d-01
    52:           c---------------------------------------------------------------------
    53:                    xcrref(1) = 1.6196343210976702d-02
    54:                    xcrref(2) = 2.1976745164821318d-03
    55:                    xcrref(3) = 1.5179927653399185d-03
    56:                    xcrref(4) = 1.5029584435994323d-03
    57:                    xcrref(5) = 3.4264073155896461d-02
    58:           
    59:           c---------------------------------------------------------------------
    60:           c   Reference values of RMS-norms of solution error, for the (12X12X12) grid,
    61:           c   after 50 time steps, with  DT = 5.0d-01
    62:           c---------------------------------------------------------------------
    63:                    xceref(1) = 6.4223319957960924d-04
    64:                    xceref(2) = 8.4144342047347926d-05
    65:                    xceref(3) = 5.8588269616485186d-05
    66:                    xceref(4) = 5.8474222595157350d-05
    67:                    xceref(5) = 1.3103347914111294d-03
    68:           
    69:           c---------------------------------------------------------------------
    70:           c   Reference value of surface integral, for the (12X12X12) grid,
    71:           c   after 50 time steps, with DT = 5.0d-01
    72:           c---------------------------------------------------------------------
    73:                    xciref = 7.8418928865937083d+00
    74:           
    75:           
    76:                   elseif ( (nx0 .eq. 33) .and.
    77:                >           (ny0 .eq. 33) .and.
    78:                >           (nz0 .eq. 33) .and.
    79:                >           (itmax . eq. 300) ) then
    80:           
    81:                      class = 'W'   !SPEC95fp size
    82:                      dtref = 1.5d-3
    83:           c---------------------------------------------------------------------
    84:           c   Reference values of RMS-norms of residual, for the (33x33x33) grid,
    85:           c   after 300 time steps, with  DT = 1.5d-3
    86:           c---------------------------------------------------------------------
    87:                      xcrref(1) =   0.1236511638192d+02
    88:                      xcrref(2) =   0.1317228477799d+01
    89:                      xcrref(3) =   0.2550120713095d+01
    90:                      xcrref(4) =   0.2326187750252d+01
    91:                      xcrref(5) =   0.2826799444189d+02
    92:           
    93:           
    94:           c---------------------------------------------------------------------
    95:           c   Reference values of RMS-norms of solution error, for the (33X33X33) grid,
    96:           c---------------------------------------------------------------------
    97:                      xceref(1) =   0.4867877144216d+00
    98:                      xceref(2) =   0.5064652880982d-01
    99:                      xceref(3) =   0.9281818101960d-01
   100:                      xceref(4) =   0.8570126542733d-01
   101:                      xceref(5) =   0.1084277417792d+01
   102:           
   103:           
   104:           c---------------------------------------------------------------------
   105:           c   Reference value of surface integral, for the (33X33X33) grid,
   106:           c   after 300 time steps, with  DT = 1.5d-3
   107:           c---------------------------------------------------------------------
   108:                      xciref    =   0.1161399311023d+02
   109:           
   110:                   elseif ( (nx0 .eq. 64) .and.
   111:                >           (ny0 .eq. 64) .and.
   112:                >           (nz0 .eq. 64) .and.
   113:                >           (itmax . eq. 250) ) then
   114:           
   115:                      class = 'A'
   116:                      dtref = 2.0d+0
   117:           c---------------------------------------------------------------------
   118:           c   Reference values of RMS-norms of residual, for the (64X64X64) grid,
   119:           c   after 250 time steps, with  DT = 2.0d+00
   120:           c---------------------------------------------------------------------
   121:                    xcrref(1) = 7.7902107606689367d+02
   122:                    xcrref(2) = 6.3402765259692870d+01
   123:                    xcrref(3) = 1.9499249727292479d+02
   124:                    xcrref(4) = 1.7845301160418537d+02
   125:                    xcrref(5) = 1.8384760349464247d+03
   126:           
   127:           c---------------------------------------------------------------------
   128:           c   Reference values of RMS-norms of solution error, for the (64X64X64) grid,
   129:           c   after 250 time steps, with  DT = 2.0d+00
   130:           c---------------------------------------------------------------------
   131:                    xceref(1) = 2.9964085685471943d+01
   132:                    xceref(2) = 2.8194576365003349d+00
   133:                    xceref(3) = 7.3473412698774742d+00
   134:                    xceref(4) = 6.7139225687777051d+00
   135:                    xceref(5) = 7.0715315688392578d+01
   136:           
   137:           c---------------------------------------------------------------------
   138:           c   Reference value of surface integral, for the (64X64X64) grid,
   139:           c   after 250 time steps, with DT = 2.0d+00
   140:           c---------------------------------------------------------------------
   141:                    xciref = 2.6030925604886277d+01
   142:           
   143:           
   144:                   elseif ( (nx0 .eq. 102) .and.
   145:                >           (ny0 .eq. 102) .and.
   146:                >           (nz0 .eq. 102) .and.
   147:                >           (itmax . eq. 250) ) then
   148:           
   149:                      class = 'B'
   150:                      dtref = 2.0d+0
   151:           
   152:           c---------------------------------------------------------------------
   153:           c   Reference values of RMS-norms of residual, for the (102X102X102) grid,
   154:           c   after 250 time steps, with  DT = 2.0d+00
   155:           c---------------------------------------------------------------------
   156:                    xcrref(1) = 3.5532672969982736d+03
   157:                    xcrref(2) = 2.6214750795310692d+02
   158:                    xcrref(3) = 8.8333721850952190d+02
   159:                    xcrref(4) = 7.7812774739425265d+02
   160:                    xcrref(5) = 7.3087969592545314d+03
   161:           
   162:           c---------------------------------------------------------------------
   163:           c   Reference values of RMS-norms of solution error, for the (102X102X102)
   164:           c   grid, after 250 time steps, with  DT = 2.0d+00
   165:           c---------------------------------------------------------------------
   166:                    xceref(1) = 1.1401176380212709d+02
   167:                    xceref(2) = 8.1098963655421574d+00
   168:                    xceref(3) = 2.8480597317698308d+01
   169:                    xceref(4) = 2.5905394567832939d+01
   170:                    xceref(5) = 2.6054907504857413d+02
   171:           
   172:           c---------------------------------------------------------------------
   173:           c   Reference value of surface integral, for the (102X102X102) grid,
   174:           c   after 250 time steps, with DT = 2.0d+00
   175:           c---------------------------------------------------------------------
   176:                    xciref = 4.7887162703308227d+01
   177:           
   178:                   elseif ( (nx0 .eq. 162) .and.
   179:                >           (ny0 .eq. 162) .and.
   180:                >           (nz0 .eq. 162) .and.
   181:                >           (itmax . eq. 250) ) then
   182:           
   183:                      class = 'C'
   184:                      dtref = 2.0d+0
   185:           
   186:           c---------------------------------------------------------------------
   187:           c   Reference values of RMS-norms of residual, for the (162X162X162) grid,
   188:           c   after 250 time steps, with  DT = 2.0d+00
   189:           c---------------------------------------------------------------------
   190:                    xcrref(1) = 1.03766980323537846d+04
   191:                    xcrref(2) = 8.92212458801008552d+02
   192:                    xcrref(3) = 2.56238814582660871d+03
   193:                    xcrref(4) = 2.19194343857831427d+03
   194:                    xcrref(5) = 1.78078057261061185d+04
   195:           
   196:           c---------------------------------------------------------------------
   197:           c   Reference values of RMS-norms of solution error, for the (162X162X162)
   198:           c   grid, after 250 time steps, with  DT = 2.0d+00
   199:           c---------------------------------------------------------------------
   200:                    xceref(1) = 2.15986399716949279d+02
   201:                    xceref(2) = 1.55789559239863600d+01
   202:                    xceref(3) = 5.41318863077207766d+01
   203:                    xceref(4) = 4.82262643154045421d+01
   204:                    xceref(5) = 4.55902910043250358d+02
   205:           
   206:           c---------------------------------------------------------------------
   207:           c   Reference value of surface integral, for the (162X162X162) grid,
   208:           c   after 250 time steps, with DT = 2.0d+00
   209:           c---------------------------------------------------------------------
   210:                    xciref = 6.66404553572181300d+01
   211:           
   212:           c---------------------------------------------------------------------
   213:           c   Reference value of surface integral, for the (162X162X162) grid,
   214:           c   after 250 time steps, with DT = 2.0d+00
   215:           c---------------------------------------------------------------------
   216:                    xciref = 6.66404553572181300d+01
   217:           
   218:                   elseif ( (nx0 .eq. 408) .and.
   219:                >           (ny0 .eq. 408) .and.
   220:                >           (nz0 .eq. 408) .and.
   221:                >           (itmax . eq. 300) ) then
   222:           
   223:                      class = 'D'
   224:                      dtref = 1.0d+0
   225:           
   226:           c---------------------------------------------------------------------
   227:           c   Reference values of RMS-norms of residual, for the (408X408X408) grid,
   228:           c   after 300 time steps, with  DT = 1.0d+00
   229:           c---------------------------------------------------------------------
   230:                    xcrref(1) = 0.4868417937025d+05
   231:                    xcrref(2) = 0.4696371050071d+04
   232:                    xcrref(3) = 0.1218114549776d+05
   233:                    xcrref(4) = 0.1033801493461d+05
   234:                    xcrref(5) = 0.7142398413817d+05
   235:           
   236:           c---------------------------------------------------------------------
   237:           c   Reference values of RMS-norms of solution error, for the (408X408X408)
   238:           c   grid, after 300 time steps, with  DT = 1.0d+00
   239:           c---------------------------------------------------------------------
   240:                    xceref(1) = 0.3752393004482d+03
   241:                    xceref(2) = 0.3084128893659d+02
   242:                    xceref(3) = 0.9434276905469d+02
   243:                    xceref(4) = 0.8230686681928d+02
   244:                    xceref(5) = 0.7002620636210d+03
   245:           
   246:           c---------------------------------------------------------------------
   247:           c   Reference value of surface integral, for the (408X408X408) grid,
   248:           c   after 300 time steps, with DT = 1.0d+00
   249:           c---------------------------------------------------------------------
   250:                    xciref =    0.8334101392503d+02
   251:           
   252:                   elseif ( (nx0 .eq. 1020) .and.
   253:                >           (ny0 .eq. 1020) .and.
   254:                >           (nz0 .eq. 1020) .and.
   255:                >           (itmax . eq. 300) ) then
   256:           
   257:                      class = 'E'
   258:                      dtref = 0.5d+0
   259:           
   260:           c---------------------------------------------------------------------
   261:           c   Reference values of RMS-norms of residual, for the (1020X1020X1020) grid,
   262:           c   after 300 time steps, with  DT = 0.5d+00
   263:           c---------------------------------------------------------------------
   264:                    xcrref(1) = 0.2099641687874d+06
   265:                    xcrref(2) = 0.2130403143165d+05
   266:                    xcrref(3) = 0.5319228789371d+05
   267:                    xcrref(4) = 0.4509761639833d+05
   268:                    xcrref(5) = 0.2932360006590d+06
   269:           
   270:           c---------------------------------------------------------------------
   271:           c   Reference values of RMS-norms of solution error, for the (1020X1020X1020)
   272:           c   grid, after 300 time steps, with  DT = 0.5d+00
   273:           c---------------------------------------------------------------------
   274:                    xceref(1) = 0.4800572578333d+03
   275:                    xceref(2) = 0.4221993400184d+02
   276:                    xceref(3) = 0.1210851906824d+03
   277:                    xceref(4) = 0.1047888986770d+03
   278:                    xceref(5) = 0.8363028257389d+03
   279:           
   280:           c---------------------------------------------------------------------
   281:           c   Reference value of surface integral, for the (1020X1020X1020) grid,
   282:           c   after 300 time steps, with DT = 0.5d+00
   283:           c---------------------------------------------------------------------
   284:                    xciref =    0.9512163272273d+02
   285:           
   286:                   elseif ( (nx0 .eq. 2560) .and.
   287:                >           (ny0 .eq. 2560) .and.
   288:                >           (nz0 .eq. 2560) .and.
   289:                >           (itmax . eq. 300) ) then
   290:           
   291:                      class = 'F'
   292:                      dtref = 0.2d+0
   293:           
   294:           c---------------------------------------------------------------------
   295:           c   Reference values of RMS-norms of residual, for the (2560X2560X2560) grid,
   296:           c   after 300 time steps, with  DT = 0.2d+00
   297:           c---------------------------------------------------------------------
   298:                    xcrref(1) = 0.8505125358152d+06
   299:                    xcrref(2) = 0.8774655318044d+05
   300:                    xcrref(3) = 0.2167258198851d+06
   301:                    xcrref(4) = 0.1838245257371d+06
   302:                    xcrref(5) = 0.1175556512415d+07
   303:           
   304:           c---------------------------------------------------------------------
   305:           c   Reference values of RMS-norms of solution error, for the (2560X2560X2560)
   306:           c   grid, after 300 time steps, with  DT = 0.2d+00
   307:           c---------------------------------------------------------------------
   308:                    xceref(1) = 0.5293914132486d+03
   309:                    xceref(2) = 0.4784861621068d+02
   310:                    xceref(3) = 0.1337701281659d+03
   311:                    xceref(4) = 0.1154215049655d+03
   312:                    xceref(5) = 0.8956266851467d+03
   313:           
   314:           c---------------------------------------------------------------------
   315:           c   Reference value of surface integral, for the (2560X2560X2560) grid,
   316:           c   after 300 time steps, with DT = 0.2d+00
   317:           c---------------------------------------------------------------------
   318:                    xciref =    0.1002509436546d+03
   319:           
   320:                   else
   321:                      verified = .FALSE.
   322:                   endif
   323:           
   324:           c---------------------------------------------------------------------
   325:           c    verification test for residuals if gridsize is one of
   326:           c    the defined grid sizes above (class .ne. 'U')
   327:           c---------------------------------------------------------------------
   328:           
   329:           c---------------------------------------------------------------------
   330:           c    Compute the difference of solution values and the known reference values.
   331:           c---------------------------------------------------------------------
   332: *------>          do m = 1, 5
   333: |          
   334: |                    xcrdif(m) = dabs((xcr(m)-xcrref(m))/xcrref(m))
   335: |                    xcedif(m) = dabs((xce(m)-xceref(m))/xceref(m))
   336: |          
   337: *------           enddo
   338:                   xcidif = dabs((xci - xciref)/xciref)
   339:           
   340:           
   341:           c---------------------------------------------------------------------
   342:           c    Output the comparison of computed results to known cases.
   343:           c---------------------------------------------------------------------
   344:           
   345:                   if (class .ne. 'U') then
   346:                      write(*, 1990) class
   347:            1990      format(/, ' Verification being performed for class ', a)
   348:                      write (*,2000) epsilon
   349:            2000      format(' Accuracy setting for epsilon = ', E20.13)
   350:                      verified = (dabs(dt-dtref) .le. epsilon)
   351:                      if (.not.verified) then
   352:                         class = 'U'
   353:                         write (*,1000) dtref
   354:            1000         format(' DT does not match the reference value of ',
   355:                >                 E15.8)
   356:                      endif
   357:                   else
   358:                      write(*, 1995)
   359:            1995      format(' Unknown class')
   360:                   endif
   361:           
   362:           
   363:                   if (class .ne. 'U') then
   364:                      write (*, 2001)
   365:                   else
   366:                      write (*, 2005)
   367:                   endif
   368:           
   369:            2001   format(' Comparison of RMS-norms of residual')
   370:            2005   format(' RMS-norms of residual')
   371: +------>          do m = 1, 5
   372: |                    if (class .eq. 'U') then
   373: |                       write(*, 2015) m, xcr(m)
   374: |                    else if ((.not.ieee_is_nan(xcrdif(m))) .and.
   375: |              &              xcrdif(m) .le. epsilon) then
   376: |                       write (*,2011) m,xcr(m),xcrref(m),xcrdif(m)
   377: |                    else
   378: |                       verified = .false.
   379: |                       write (*,2010) m,xcr(m),xcrref(m),xcrdif(m)
   380: |                    endif
   381: +------           enddo
   382:           
   383:                   if (class .ne. 'U') then
   384:                      write (*,2002)
   385:                   else
   386:                      write (*,2006)
   387:                   endif
   388:            2002   format(' Comparison of RMS-norms of solution error')
   389:            2006   format(' RMS-norms of solution error')
   390:            
   391: +------>          do m = 1, 5
   392: |                    if (class .eq. 'U') then
   393: |                       write(*, 2015) m, xce(m)
   394: |                    else if ((.not.ieee_is_nan(xcedif(m))) .and.
   395: |              &              xcedif(m) .le. epsilon) then
   396: |                       write (*,2011) m,xce(m),xceref(m),xcedif(m)
   397: |                    else
   398: |                       verified = .false.
   399: |                       write (*,2010) m,xce(m),xceref(m),xcedif(m)
   400: |                    endif
   401: +------           enddo
   402:            
   403:            2010   format(' FAILURE: ', i2, 2x, E20.13, E20.13, E20.13)
   404:            2011   format('          ', i2, 2x, E20.13, E20.13, E20.13)
   405:            2015   format('          ', i2, 2x, E20.13)
   406:            
   407:                   if (class .ne. 'U') then
   408:                      write (*,2025)
   409:                   else
   410:                      write (*,2026)
   411:                   endif
   412:            2025   format(' Comparison of surface integral')
   413:            2026   format(' Surface integral')
   414:           
   415:           
   416:                   if (class .eq. 'U') then
   417:                      write(*, 2030) xci
   418:                   else if ((.not.ieee_is_nan(xcidif)) .and.
   419:                &           xcidif .le. epsilon) then
   420:                      write(*, 2032) xci, xciref, xcidif
   421:                   else
   422:                      verified = .false.
   423:                      write(*, 2031) xci, xciref, xcidif
   424:                   endif
   425:           
   426:            2030   format('          ', 4x, E20.13)
   427:            2031   format(' FAILURE: ', 4x, E20.13, E20.13, E20.13)
   428:            2032   format('          ', 4x, E20.13, E20.13, E20.13)
   429:           
   430:           
   431:           
   432:                   if (class .eq. 'U') then
   433:                      write(*, 2022)
   434:                      write(*, 2023)
   435:            2022      format(' No reference values provided')
   436:            2023      format(' No verification performed')
   437:                   else if (verified) then
   438:                      write(*, 2020)
   439:            2020      format(' Verification Successful')
   440:                   else
   441:                      write(*, 2021)
   442:            2021      format(' Verification failed')
   443:                   endif
   444:           
   445:                   return
   446:           
   447:           
   448:                   end


