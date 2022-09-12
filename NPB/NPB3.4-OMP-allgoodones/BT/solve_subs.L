NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: MATVEC_SUB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: MATVEC_SUB
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                 subroutine matvec_sub(ablock,avec,bvec)
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:           c---------------------------------------------------------------------
    11:           c     subtracts bvec=bvec - ablock*avec
    12:           c---------------------------------------------------------------------
    13:           
    14:                 implicit none
    15:           
    16:                 double precision ablock,avec,bvec
    17:                 dimension ablock(5,5),avec(5),bvec(5)
    18:           
    19:           c---------------------------------------------------------------------
    20:           c            rhs(i,ic,jc,kc) = rhs(i,ic,jc,kc)
    21:           c     $           - lhs(i,1,ablock,ia)*
    22:           c---------------------------------------------------------------------
    23:                    bvec(1) = bvec(1) - ablock(1,1)*avec(1)
    24:                >                     - ablock(1,2)*avec(2)
    25:                >                     - ablock(1,3)*avec(3)
    26:                >                     - ablock(1,4)*avec(4)
    27:                >                     - ablock(1,5)*avec(5)
    28:                    bvec(2) = bvec(2) - ablock(2,1)*avec(1)
    29:                >                     - ablock(2,2)*avec(2)
    30:                >                     - ablock(2,3)*avec(3)
    31:                >                     - ablock(2,4)*avec(4)
    32:                >                     - ablock(2,5)*avec(5)
    33:                    bvec(3) = bvec(3) - ablock(3,1)*avec(1)
    34:                >                     - ablock(3,2)*avec(2)
    35:                >                     - ablock(3,3)*avec(3)
    36:                >                     - ablock(3,4)*avec(4)
    37:                >                     - ablock(3,5)*avec(5)
    38:                    bvec(4) = bvec(4) - ablock(4,1)*avec(1)
    39:                >                     - ablock(4,2)*avec(2)
    40:                >                     - ablock(4,3)*avec(3)
    41:                >                     - ablock(4,4)*avec(4)
    42:                >                     - ablock(4,5)*avec(5)
    43:                    bvec(5) = bvec(5) - ablock(5,1)*avec(1)
    44:                >                     - ablock(5,2)*avec(2)
    45:                >                     - ablock(5,3)*avec(3)
    46:                >                     - ablock(5,4)*avec(4)
    47:                >                     - ablock(5,5)*avec(5)
    48:           
    49:           
    50:                 return
    51:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: MATMUL_SUB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: MATMUL_SUB
FORMAT LIST

 LINE   LOOP      STATEMENT

    56:                 subroutine matmul_sub(ablock, bblock, cblock)
    57:           
    58:           c---------------------------------------------------------------------
    59:           c---------------------------------------------------------------------
    60:           
    61:           c---------------------------------------------------------------------
    62:           c     subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
    63:           c---------------------------------------------------------------------
    64:           
    65:                 implicit none
    66:           
    67:                 double precision ablock, bblock, cblock
    68:                 dimension ablock(5,5), bblock(5,5), cblock(5,5)
    69:           
    70:           
    71:                    cblock(1,1) = cblock(1,1) - ablock(1,1)*bblock(1,1)
    72:                >                             - ablock(1,2)*bblock(2,1)
    73:                >                             - ablock(1,3)*bblock(3,1)
    74:                >                             - ablock(1,4)*bblock(4,1)
    75:                >                             - ablock(1,5)*bblock(5,1)
    76:                    cblock(2,1) = cblock(2,1) - ablock(2,1)*bblock(1,1)
    77:                >                             - ablock(2,2)*bblock(2,1)
    78:                >                             - ablock(2,3)*bblock(3,1)
    79:                >                             - ablock(2,4)*bblock(4,1)
    80:                >                             - ablock(2,5)*bblock(5,1)
    81:                    cblock(3,1) = cblock(3,1) - ablock(3,1)*bblock(1,1)
    82:                >                             - ablock(3,2)*bblock(2,1)
    83:                >                             - ablock(3,3)*bblock(3,1)
    84:                >                             - ablock(3,4)*bblock(4,1)
    85:                >                             - ablock(3,5)*bblock(5,1)
    86:                    cblock(4,1) = cblock(4,1) - ablock(4,1)*bblock(1,1)
    87:                >                             - ablock(4,2)*bblock(2,1)
    88:                >                             - ablock(4,3)*bblock(3,1)
    89:                >                             - ablock(4,4)*bblock(4,1)
    90:                >                             - ablock(4,5)*bblock(5,1)
    91:                    cblock(5,1) = cblock(5,1) - ablock(5,1)*bblock(1,1)
    92:                >                             - ablock(5,2)*bblock(2,1)
    93:                >                             - ablock(5,3)*bblock(3,1)
    94:                >                             - ablock(5,4)*bblock(4,1)
    95:                >                             - ablock(5,5)*bblock(5,1)
    96:                    cblock(1,2) = cblock(1,2) - ablock(1,1)*bblock(1,2)
    97:                >                             - ablock(1,2)*bblock(2,2)
    98:                >                             - ablock(1,3)*bblock(3,2)
    99:                >                             - ablock(1,4)*bblock(4,2)
   100:                >                             - ablock(1,5)*bblock(5,2)
   101:                    cblock(2,2) = cblock(2,2) - ablock(2,1)*bblock(1,2)
   102:                >                             - ablock(2,2)*bblock(2,2)
   103:                >                             - ablock(2,3)*bblock(3,2)
   104:                >                             - ablock(2,4)*bblock(4,2)
   105:                >                             - ablock(2,5)*bblock(5,2)
   106:                    cblock(3,2) = cblock(3,2) - ablock(3,1)*bblock(1,2)
   107:                >                             - ablock(3,2)*bblock(2,2)
   108:                >                             - ablock(3,3)*bblock(3,2)
   109:                >                             - ablock(3,4)*bblock(4,2)
   110:                >                             - ablock(3,5)*bblock(5,2)
   111:                    cblock(4,2) = cblock(4,2) - ablock(4,1)*bblock(1,2)
   112:                >                             - ablock(4,2)*bblock(2,2)
   113:                >                             - ablock(4,3)*bblock(3,2)
   114:                >                             - ablock(4,4)*bblock(4,2)
   115:                >                             - ablock(4,5)*bblock(5,2)
   116:                    cblock(5,2) = cblock(5,2) - ablock(5,1)*bblock(1,2)
   117:                >                             - ablock(5,2)*bblock(2,2)
   118:                >                             - ablock(5,3)*bblock(3,2)
   119:                >                             - ablock(5,4)*bblock(4,2)
   120:                >                             - ablock(5,5)*bblock(5,2)
   121:                    cblock(1,3) = cblock(1,3) - ablock(1,1)*bblock(1,3)
   122:                >                             - ablock(1,2)*bblock(2,3)
   123:                >                             - ablock(1,3)*bblock(3,3)
   124:                >                             - ablock(1,4)*bblock(4,3)
   125:                >                             - ablock(1,5)*bblock(5,3)
   126:                    cblock(2,3) = cblock(2,3) - ablock(2,1)*bblock(1,3)
   127:                >                             - ablock(2,2)*bblock(2,3)
   128:                >                             - ablock(2,3)*bblock(3,3)
   129:                >                             - ablock(2,4)*bblock(4,3)
   130:                >                             - ablock(2,5)*bblock(5,3)
   131:                    cblock(3,3) = cblock(3,3) - ablock(3,1)*bblock(1,3)
   132:                >                             - ablock(3,2)*bblock(2,3)
   133:                >                             - ablock(3,3)*bblock(3,3)
   134:                >                             - ablock(3,4)*bblock(4,3)
   135:                >                             - ablock(3,5)*bblock(5,3)
   136:                    cblock(4,3) = cblock(4,3) - ablock(4,1)*bblock(1,3)
   137:                >                             - ablock(4,2)*bblock(2,3)
   138:                >                             - ablock(4,3)*bblock(3,3)
   139:                >                             - ablock(4,4)*bblock(4,3)
   140:                >                             - ablock(4,5)*bblock(5,3)
   141:                    cblock(5,3) = cblock(5,3) - ablock(5,1)*bblock(1,3)
   142:                >                             - ablock(5,2)*bblock(2,3)
   143:                >                             - ablock(5,3)*bblock(3,3)
   144:                >                             - ablock(5,4)*bblock(4,3)
   145:                >                             - ablock(5,5)*bblock(5,3)
   146:                    cblock(1,4) = cblock(1,4) - ablock(1,1)*bblock(1,4)
   147:                >                             - ablock(1,2)*bblock(2,4)
   148:                >                             - ablock(1,3)*bblock(3,4)
   149:                >                             - ablock(1,4)*bblock(4,4)
   150:                >                             - ablock(1,5)*bblock(5,4)
   151:                    cblock(2,4) = cblock(2,4) - ablock(2,1)*bblock(1,4)
   152:                >                             - ablock(2,2)*bblock(2,4)
   153:                >                             - ablock(2,3)*bblock(3,4)
   154:                >                             - ablock(2,4)*bblock(4,4)
   155:                >                             - ablock(2,5)*bblock(5,4)
   156:                    cblock(3,4) = cblock(3,4) - ablock(3,1)*bblock(1,4)
   157:                >                             - ablock(3,2)*bblock(2,4)
   158:                >                             - ablock(3,3)*bblock(3,4)
   159:                >                             - ablock(3,4)*bblock(4,4)
   160:                >                             - ablock(3,5)*bblock(5,4)
   161:                    cblock(4,4) = cblock(4,4) - ablock(4,1)*bblock(1,4)
   162:                >                             - ablock(4,2)*bblock(2,4)
   163:                >                             - ablock(4,3)*bblock(3,4)
   164:                >                             - ablock(4,4)*bblock(4,4)
   165:                >                             - ablock(4,5)*bblock(5,4)
   166:                    cblock(5,4) = cblock(5,4) - ablock(5,1)*bblock(1,4)
   167:                >                             - ablock(5,2)*bblock(2,4)
   168:                >                             - ablock(5,3)*bblock(3,4)
   169:                >                             - ablock(5,4)*bblock(4,4)
   170:                >                             - ablock(5,5)*bblock(5,4)
   171:                    cblock(1,5) = cblock(1,5) - ablock(1,1)*bblock(1,5)
   172:                >                             - ablock(1,2)*bblock(2,5)
   173:                >                             - ablock(1,3)*bblock(3,5)
   174:                >                             - ablock(1,4)*bblock(4,5)
   175:                >                             - ablock(1,5)*bblock(5,5)
   176:                    cblock(2,5) = cblock(2,5) - ablock(2,1)*bblock(1,5)
   177:                >                             - ablock(2,2)*bblock(2,5)
   178:                >                             - ablock(2,3)*bblock(3,5)
   179:                >                             - ablock(2,4)*bblock(4,5)
   180:                >                             - ablock(2,5)*bblock(5,5)
   181:                    cblock(3,5) = cblock(3,5) - ablock(3,1)*bblock(1,5)
   182:                >                             - ablock(3,2)*bblock(2,5)
   183:                >                             - ablock(3,3)*bblock(3,5)
   184:                >                             - ablock(3,4)*bblock(4,5)
   185:                >                             - ablock(3,5)*bblock(5,5)
   186:                    cblock(4,5) = cblock(4,5) - ablock(4,1)*bblock(1,5)
   187:                >                             - ablock(4,2)*bblock(2,5)
   188:                >                             - ablock(4,3)*bblock(3,5)
   189:                >                             - ablock(4,4)*bblock(4,5)
   190:                >                             - ablock(4,5)*bblock(5,5)
   191:                    cblock(5,5) = cblock(5,5) - ablock(5,1)*bblock(1,5)
   192:                >                             - ablock(5,2)*bblock(2,5)
   193:                >                             - ablock(5,3)*bblock(3,5)
   194:                >                             - ablock(5,4)*bblock(4,5)
   195:                >                             - ablock(5,5)*bblock(5,5)
   196:           
   197:            
   198:                 return
   199:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: BINVCRHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: BINVCRHS
FORMAT LIST

 LINE   LOOP      STATEMENT

   206:                 subroutine binvcrhs( lhs,c,r )
   207:           
   208:           c---------------------------------------------------------------------
   209:           c---------------------------------------------------------------------
   210:           
   211:           c---------------------------------------------------------------------
   212:           c
   213:           c---------------------------------------------------------------------
   214:           
   215:                 implicit none
   216:           
   217:                 double precision pivot, coeff, lhs
   218:                 dimension lhs(5,5)
   219:                 double precision c(5,5), r(5)
   220:           
   221:           c---------------------------------------------------------------------
   222:           c
   223:           c---------------------------------------------------------------------
   224:           
   225:                 pivot = 1.00d0/lhs(1,1)
   226:                 lhs(1,2) = lhs(1,2)*pivot
   227:                 lhs(1,3) = lhs(1,3)*pivot
   228:                 lhs(1,4) = lhs(1,4)*pivot
   229:                 lhs(1,5) = lhs(1,5)*pivot
   230:                 c(1,1) = c(1,1)*pivot
   231:                 c(1,2) = c(1,2)*pivot
   232:                 c(1,3) = c(1,3)*pivot
   233:                 c(1,4) = c(1,4)*pivot
   234:                 c(1,5) = c(1,5)*pivot
   235:                 r(1)   = r(1)  *pivot
   236:           
   237:                 coeff = lhs(2,1)
   238:                 lhs(2,2)= lhs(2,2) - coeff*lhs(1,2)
   239:                 lhs(2,3)= lhs(2,3) - coeff*lhs(1,3)
   240:                 lhs(2,4)= lhs(2,4) - coeff*lhs(1,4)
   241:                 lhs(2,5)= lhs(2,5) - coeff*lhs(1,5)
   242:                 c(2,1) = c(2,1) - coeff*c(1,1)
   243:                 c(2,2) = c(2,2) - coeff*c(1,2)
   244:                 c(2,3) = c(2,3) - coeff*c(1,3)
   245:                 c(2,4) = c(2,4) - coeff*c(1,4)
   246:                 c(2,5) = c(2,5) - coeff*c(1,5)
   247:                 r(2)   = r(2)   - coeff*r(1)
   248:           
   249:                 coeff = lhs(3,1)
   250:                 lhs(3,2)= lhs(3,2) - coeff*lhs(1,2)
   251:                 lhs(3,3)= lhs(3,3) - coeff*lhs(1,3)
   252:                 lhs(3,4)= lhs(3,4) - coeff*lhs(1,4)
   253:                 lhs(3,5)= lhs(3,5) - coeff*lhs(1,5)
   254:                 c(3,1) = c(3,1) - coeff*c(1,1)
   255:                 c(3,2) = c(3,2) - coeff*c(1,2)
   256:                 c(3,3) = c(3,3) - coeff*c(1,3)
   257:                 c(3,4) = c(3,4) - coeff*c(1,4)
   258:                 c(3,5) = c(3,5) - coeff*c(1,5)
   259:                 r(3)   = r(3)   - coeff*r(1)
   260:           
   261:                 coeff = lhs(4,1)
   262:                 lhs(4,2)= lhs(4,2) - coeff*lhs(1,2)
   263:                 lhs(4,3)= lhs(4,3) - coeff*lhs(1,3)
   264:                 lhs(4,4)= lhs(4,4) - coeff*lhs(1,4)
   265:                 lhs(4,5)= lhs(4,5) - coeff*lhs(1,5)
   266:                 c(4,1) = c(4,1) - coeff*c(1,1)
   267:                 c(4,2) = c(4,2) - coeff*c(1,2)
   268:                 c(4,3) = c(4,3) - coeff*c(1,3)
   269:                 c(4,4) = c(4,4) - coeff*c(1,4)
   270:                 c(4,5) = c(4,5) - coeff*c(1,5)
   271:                 r(4)   = r(4)   - coeff*r(1)
   272:           
   273:                 coeff = lhs(5,1)
   274:                 lhs(5,2)= lhs(5,2) - coeff*lhs(1,2)
   275:                 lhs(5,3)= lhs(5,3) - coeff*lhs(1,3)
   276:                 lhs(5,4)= lhs(5,4) - coeff*lhs(1,4)
   277:                 lhs(5,5)= lhs(5,5) - coeff*lhs(1,5)
   278:                 c(5,1) = c(5,1) - coeff*c(1,1)
   279:                 c(5,2) = c(5,2) - coeff*c(1,2)
   280:                 c(5,3) = c(5,3) - coeff*c(1,3)
   281:                 c(5,4) = c(5,4) - coeff*c(1,4)
   282:                 c(5,5) = c(5,5) - coeff*c(1,5)
   283:                 r(5)   = r(5)   - coeff*r(1)
   284:           
   285:           
   286:                 pivot = 1.00d0/lhs(2,2)
   287:                 lhs(2,3) = lhs(2,3)*pivot
   288:                 lhs(2,4) = lhs(2,4)*pivot
   289:                 lhs(2,5) = lhs(2,5)*pivot
   290:                 c(2,1) = c(2,1)*pivot
   291:                 c(2,2) = c(2,2)*pivot
   292:                 c(2,3) = c(2,3)*pivot
   293:                 c(2,4) = c(2,4)*pivot
   294:                 c(2,5) = c(2,5)*pivot
   295:                 r(2)   = r(2)  *pivot
   296:           
   297:                 coeff = lhs(1,2)
   298:                 lhs(1,3)= lhs(1,3) - coeff*lhs(2,3)
   299:                 lhs(1,4)= lhs(1,4) - coeff*lhs(2,4)
   300:                 lhs(1,5)= lhs(1,5) - coeff*lhs(2,5)
   301:                 c(1,1) = c(1,1) - coeff*c(2,1)
   302:                 c(1,2) = c(1,2) - coeff*c(2,2)
   303:                 c(1,3) = c(1,3) - coeff*c(2,3)
   304:                 c(1,4) = c(1,4) - coeff*c(2,4)
   305:                 c(1,5) = c(1,5) - coeff*c(2,5)
   306:                 r(1)   = r(1)   - coeff*r(2)
   307:           
   308:                 coeff = lhs(3,2)
   309:                 lhs(3,3)= lhs(3,3) - coeff*lhs(2,3)
   310:                 lhs(3,4)= lhs(3,4) - coeff*lhs(2,4)
   311:                 lhs(3,5)= lhs(3,5) - coeff*lhs(2,5)
   312:                 c(3,1) = c(3,1) - coeff*c(2,1)
   313:                 c(3,2) = c(3,2) - coeff*c(2,2)
   314:                 c(3,3) = c(3,3) - coeff*c(2,3)
   315:                 c(3,4) = c(3,4) - coeff*c(2,4)
   316:                 c(3,5) = c(3,5) - coeff*c(2,5)
   317:                 r(3)   = r(3)   - coeff*r(2)
   318:           
   319:                 coeff = lhs(4,2)
   320:                 lhs(4,3)= lhs(4,3) - coeff*lhs(2,3)
   321:                 lhs(4,4)= lhs(4,4) - coeff*lhs(2,4)
   322:                 lhs(4,5)= lhs(4,5) - coeff*lhs(2,5)
   323:                 c(4,1) = c(4,1) - coeff*c(2,1)
   324:                 c(4,2) = c(4,2) - coeff*c(2,2)
   325:                 c(4,3) = c(4,3) - coeff*c(2,3)
   326:                 c(4,4) = c(4,4) - coeff*c(2,4)
   327:                 c(4,5) = c(4,5) - coeff*c(2,5)
   328:                 r(4)   = r(4)   - coeff*r(2)
   329:           
   330:                 coeff = lhs(5,2)
   331:                 lhs(5,3)= lhs(5,3) - coeff*lhs(2,3)
   332:                 lhs(5,4)= lhs(5,4) - coeff*lhs(2,4)
   333:                 lhs(5,5)= lhs(5,5) - coeff*lhs(2,5)
   334:                 c(5,1) = c(5,1) - coeff*c(2,1)
   335:                 c(5,2) = c(5,2) - coeff*c(2,2)
   336:                 c(5,3) = c(5,3) - coeff*c(2,3)
   337:                 c(5,4) = c(5,4) - coeff*c(2,4)
   338:                 c(5,5) = c(5,5) - coeff*c(2,5)
   339:                 r(5)   = r(5)   - coeff*r(2)
   340:           
   341:           
   342:                 pivot = 1.00d0/lhs(3,3)
   343:                 lhs(3,4) = lhs(3,4)*pivot
   344:                 lhs(3,5) = lhs(3,5)*pivot
   345:                 c(3,1) = c(3,1)*pivot
   346:                 c(3,2) = c(3,2)*pivot
   347:                 c(3,3) = c(3,3)*pivot
   348:                 c(3,4) = c(3,4)*pivot
   349:                 c(3,5) = c(3,5)*pivot
   350:                 r(3)   = r(3)  *pivot
   351:           
   352:                 coeff = lhs(1,3)
   353:                 lhs(1,4)= lhs(1,4) - coeff*lhs(3,4)
   354:                 lhs(1,5)= lhs(1,5) - coeff*lhs(3,5)
   355:                 c(1,1) = c(1,1) - coeff*c(3,1)
   356:                 c(1,2) = c(1,2) - coeff*c(3,2)
   357:                 c(1,3) = c(1,3) - coeff*c(3,3)
   358:                 c(1,4) = c(1,4) - coeff*c(3,4)
   359:                 c(1,5) = c(1,5) - coeff*c(3,5)
   360:                 r(1)   = r(1)   - coeff*r(3)
   361:           
   362:                 coeff = lhs(2,3)
   363:                 lhs(2,4)= lhs(2,4) - coeff*lhs(3,4)
   364:                 lhs(2,5)= lhs(2,5) - coeff*lhs(3,5)
   365:                 c(2,1) = c(2,1) - coeff*c(3,1)
   366:                 c(2,2) = c(2,2) - coeff*c(3,2)
   367:                 c(2,3) = c(2,3) - coeff*c(3,3)
   368:                 c(2,4) = c(2,4) - coeff*c(3,4)
   369:                 c(2,5) = c(2,5) - coeff*c(3,5)
   370:                 r(2)   = r(2)   - coeff*r(3)
   371:           
   372:                 coeff = lhs(4,3)
   373:                 lhs(4,4)= lhs(4,4) - coeff*lhs(3,4)
   374:                 lhs(4,5)= lhs(4,5) - coeff*lhs(3,5)
   375:                 c(4,1) = c(4,1) - coeff*c(3,1)
   376:                 c(4,2) = c(4,2) - coeff*c(3,2)
   377:                 c(4,3) = c(4,3) - coeff*c(3,3)
   378:                 c(4,4) = c(4,4) - coeff*c(3,4)
   379:                 c(4,5) = c(4,5) - coeff*c(3,5)
   380:                 r(4)   = r(4)   - coeff*r(3)
   381:           
   382:                 coeff = lhs(5,3)
   383:                 lhs(5,4)= lhs(5,4) - coeff*lhs(3,4)
   384:                 lhs(5,5)= lhs(5,5) - coeff*lhs(3,5)
   385:                 c(5,1) = c(5,1) - coeff*c(3,1)
   386:                 c(5,2) = c(5,2) - coeff*c(3,2)
   387:                 c(5,3) = c(5,3) - coeff*c(3,3)
   388:                 c(5,4) = c(5,4) - coeff*c(3,4)
   389:                 c(5,5) = c(5,5) - coeff*c(3,5)
   390:                 r(5)   = r(5)   - coeff*r(3)
   391:           
   392:           
   393:                 pivot = 1.00d0/lhs(4,4)
   394:                 lhs(4,5) = lhs(4,5)*pivot
   395:                 c(4,1) = c(4,1)*pivot
   396:                 c(4,2) = c(4,2)*pivot
   397:                 c(4,3) = c(4,3)*pivot
   398:                 c(4,4) = c(4,4)*pivot
   399:                 c(4,5) = c(4,5)*pivot
   400:                 r(4)   = r(4)  *pivot
   401:           
   402:                 coeff = lhs(1,4)
   403:                 lhs(1,5)= lhs(1,5) - coeff*lhs(4,5)
   404:                 c(1,1) = c(1,1) - coeff*c(4,1)
   405:                 c(1,2) = c(1,2) - coeff*c(4,2)
   406:                 c(1,3) = c(1,3) - coeff*c(4,3)
   407:                 c(1,4) = c(1,4) - coeff*c(4,4)
   408:                 c(1,5) = c(1,5) - coeff*c(4,5)
   409:                 r(1)   = r(1)   - coeff*r(4)
   410:           
   411:                 coeff = lhs(2,4)
   412:                 lhs(2,5)= lhs(2,5) - coeff*lhs(4,5)
   413:                 c(2,1) = c(2,1) - coeff*c(4,1)
   414:                 c(2,2) = c(2,2) - coeff*c(4,2)
   415:                 c(2,3) = c(2,3) - coeff*c(4,3)
   416:                 c(2,4) = c(2,4) - coeff*c(4,4)
   417:                 c(2,5) = c(2,5) - coeff*c(4,5)
   418:                 r(2)   = r(2)   - coeff*r(4)
   419:           
   420:                 coeff = lhs(3,4)
   421:                 lhs(3,5)= lhs(3,5) - coeff*lhs(4,5)
   422:                 c(3,1) = c(3,1) - coeff*c(4,1)
   423:                 c(3,2) = c(3,2) - coeff*c(4,2)
   424:                 c(3,3) = c(3,3) - coeff*c(4,3)
   425:                 c(3,4) = c(3,4) - coeff*c(4,4)
   426:                 c(3,5) = c(3,5) - coeff*c(4,5)
   427:                 r(3)   = r(3)   - coeff*r(4)
   428:           
   429:                 coeff = lhs(5,4)
   430:                 lhs(5,5)= lhs(5,5) - coeff*lhs(4,5)
   431:                 c(5,1) = c(5,1) - coeff*c(4,1)
   432:                 c(5,2) = c(5,2) - coeff*c(4,2)
   433:                 c(5,3) = c(5,3) - coeff*c(4,3)
   434:                 c(5,4) = c(5,4) - coeff*c(4,4)
   435:                 c(5,5) = c(5,5) - coeff*c(4,5)
   436:                 r(5)   = r(5)   - coeff*r(4)
   437:           
   438:           
   439:                 pivot = 1.00d0/lhs(5,5)
   440:                 c(5,1) = c(5,1)*pivot
   441:                 c(5,2) = c(5,2)*pivot
   442:                 c(5,3) = c(5,3)*pivot
   443:                 c(5,4) = c(5,4)*pivot
   444:                 c(5,5) = c(5,5)*pivot
   445:                 r(5)   = r(5)  *pivot
   446:           
   447:                 coeff = lhs(1,5)
   448:                 c(1,1) = c(1,1) - coeff*c(5,1)
   449:                 c(1,2) = c(1,2) - coeff*c(5,2)
   450:                 c(1,3) = c(1,3) - coeff*c(5,3)
   451:                 c(1,4) = c(1,4) - coeff*c(5,4)
   452:                 c(1,5) = c(1,5) - coeff*c(5,5)
   453:                 r(1)   = r(1)   - coeff*r(5)
   454:           
   455:                 coeff = lhs(2,5)
   456:                 c(2,1) = c(2,1) - coeff*c(5,1)
   457:                 c(2,2) = c(2,2) - coeff*c(5,2)
   458:                 c(2,3) = c(2,3) - coeff*c(5,3)
   459:                 c(2,4) = c(2,4) - coeff*c(5,4)
   460:                 c(2,5) = c(2,5) - coeff*c(5,5)
   461:                 r(2)   = r(2)   - coeff*r(5)
   462:           
   463:                 coeff = lhs(3,5)
   464:                 c(3,1) = c(3,1) - coeff*c(5,1)
   465:                 c(3,2) = c(3,2) - coeff*c(5,2)
   466:                 c(3,3) = c(3,3) - coeff*c(5,3)
   467:                 c(3,4) = c(3,4) - coeff*c(5,4)
   468:                 c(3,5) = c(3,5) - coeff*c(5,5)
   469:                 r(3)   = r(3)   - coeff*r(5)
   470:           
   471:                 coeff = lhs(4,5)
   472:                 c(4,1) = c(4,1) - coeff*c(5,1)
   473:                 c(4,2) = c(4,2) - coeff*c(5,2)
   474:                 c(4,3) = c(4,3) - coeff*c(5,3)
   475:                 c(4,4) = c(4,4) - coeff*c(5,4)
   476:                 c(4,5) = c(4,5) - coeff*c(5,5)
   477:                 r(4)   = r(4)   - coeff*r(5)
   478:           
   479:           
   480:                 return
   481:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: BINVRHS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:35 2020
FILE NAME: solve_subs.f

PROCEDURE NAME: BINVRHS
FORMAT LIST

 LINE   LOOP      STATEMENT

   488:                 subroutine binvrhs( lhs,r )
   489:           
   490:           c---------------------------------------------------------------------
   491:           c---------------------------------------------------------------------
   492:           
   493:           c---------------------------------------------------------------------
   494:           c
   495:           c---------------------------------------------------------------------
   496:           
   497:                 implicit none
   498:           
   499:                 double precision pivot, coeff, lhs
   500:                 dimension lhs(5,5)
   501:                 double precision r(5)
   502:           
   503:           c---------------------------------------------------------------------
   504:           c
   505:           c---------------------------------------------------------------------
   506:           
   507:           
   508:                 pivot = 1.00d0/lhs(1,1)
   509:                 lhs(1,2) = lhs(1,2)*pivot
   510:                 lhs(1,3) = lhs(1,3)*pivot
   511:                 lhs(1,4) = lhs(1,4)*pivot
   512:                 lhs(1,5) = lhs(1,5)*pivot
   513:                 r(1)   = r(1)  *pivot
   514:           
   515:                 coeff = lhs(2,1)
   516:                 lhs(2,2)= lhs(2,2) - coeff*lhs(1,2)
   517:                 lhs(2,3)= lhs(2,3) - coeff*lhs(1,3)
   518:                 lhs(2,4)= lhs(2,4) - coeff*lhs(1,4)
   519:                 lhs(2,5)= lhs(2,5) - coeff*lhs(1,5)
   520:                 r(2)   = r(2)   - coeff*r(1)
   521:           
   522:                 coeff = lhs(3,1)
   523:                 lhs(3,2)= lhs(3,2) - coeff*lhs(1,2)
   524:                 lhs(3,3)= lhs(3,3) - coeff*lhs(1,3)
   525:                 lhs(3,4)= lhs(3,4) - coeff*lhs(1,4)
   526:                 lhs(3,5)= lhs(3,5) - coeff*lhs(1,5)
   527:                 r(3)   = r(3)   - coeff*r(1)
   528:           
   529:                 coeff = lhs(4,1)
   530:                 lhs(4,2)= lhs(4,2) - coeff*lhs(1,2)
   531:                 lhs(4,3)= lhs(4,3) - coeff*lhs(1,3)
   532:                 lhs(4,4)= lhs(4,4) - coeff*lhs(1,4)
   533:                 lhs(4,5)= lhs(4,5) - coeff*lhs(1,5)
   534:                 r(4)   = r(4)   - coeff*r(1)
   535:           
   536:                 coeff = lhs(5,1)
   537:                 lhs(5,2)= lhs(5,2) - coeff*lhs(1,2)
   538:                 lhs(5,3)= lhs(5,3) - coeff*lhs(1,3)
   539:                 lhs(5,4)= lhs(5,4) - coeff*lhs(1,4)
   540:                 lhs(5,5)= lhs(5,5) - coeff*lhs(1,5)
   541:                 r(5)   = r(5)   - coeff*r(1)
   542:           
   543:           
   544:                 pivot = 1.00d0/lhs(2,2)
   545:                 lhs(2,3) = lhs(2,3)*pivot
   546:                 lhs(2,4) = lhs(2,4)*pivot
   547:                 lhs(2,5) = lhs(2,5)*pivot
   548:                 r(2)   = r(2)  *pivot
   549:           
   550:                 coeff = lhs(1,2)
   551:                 lhs(1,3)= lhs(1,3) - coeff*lhs(2,3)
   552:                 lhs(1,4)= lhs(1,4) - coeff*lhs(2,4)
   553:                 lhs(1,5)= lhs(1,5) - coeff*lhs(2,5)
   554:                 r(1)   = r(1)   - coeff*r(2)
   555:           
   556:                 coeff = lhs(3,2)
   557:                 lhs(3,3)= lhs(3,3) - coeff*lhs(2,3)
   558:                 lhs(3,4)= lhs(3,4) - coeff*lhs(2,4)
   559:                 lhs(3,5)= lhs(3,5) - coeff*lhs(2,5)
   560:                 r(3)   = r(3)   - coeff*r(2)
   561:           
   562:                 coeff = lhs(4,2)
   563:                 lhs(4,3)= lhs(4,3) - coeff*lhs(2,3)
   564:                 lhs(4,4)= lhs(4,4) - coeff*lhs(2,4)
   565:                 lhs(4,5)= lhs(4,5) - coeff*lhs(2,5)
   566:                 r(4)   = r(4)   - coeff*r(2)
   567:           
   568:                 coeff = lhs(5,2)
   569:                 lhs(5,3)= lhs(5,3) - coeff*lhs(2,3)
   570:                 lhs(5,4)= lhs(5,4) - coeff*lhs(2,4)
   571:                 lhs(5,5)= lhs(5,5) - coeff*lhs(2,5)
   572:                 r(5)   = r(5)   - coeff*r(2)
   573:           
   574:           
   575:                 pivot = 1.00d0/lhs(3,3)
   576:                 lhs(3,4) = lhs(3,4)*pivot
   577:                 lhs(3,5) = lhs(3,5)*pivot
   578:                 r(3)   = r(3)  *pivot
   579:           
   580:                 coeff = lhs(1,3)
   581:                 lhs(1,4)= lhs(1,4) - coeff*lhs(3,4)
   582:                 lhs(1,5)= lhs(1,5) - coeff*lhs(3,5)
   583:                 r(1)   = r(1)   - coeff*r(3)
   584:           
   585:                 coeff = lhs(2,3)
   586:                 lhs(2,4)= lhs(2,4) - coeff*lhs(3,4)
   587:                 lhs(2,5)= lhs(2,5) - coeff*lhs(3,5)
   588:                 r(2)   = r(2)   - coeff*r(3)
   589:           
   590:                 coeff = lhs(4,3)
   591:                 lhs(4,4)= lhs(4,4) - coeff*lhs(3,4)
   592:                 lhs(4,5)= lhs(4,5) - coeff*lhs(3,5)
   593:                 r(4)   = r(4)   - coeff*r(3)
   594:           
   595:                 coeff = lhs(5,3)
   596:                 lhs(5,4)= lhs(5,4) - coeff*lhs(3,4)
   597:                 lhs(5,5)= lhs(5,5) - coeff*lhs(3,5)
   598:                 r(5)   = r(5)   - coeff*r(3)
   599:           
   600:           
   601:                 pivot = 1.00d0/lhs(4,4)
   602:                 lhs(4,5) = lhs(4,5)*pivot
   603:                 r(4)   = r(4)  *pivot
   604:           
   605:                 coeff = lhs(1,4)
   606:                 lhs(1,5)= lhs(1,5) - coeff*lhs(4,5)
   607:                 r(1)   = r(1)   - coeff*r(4)
   608:           
   609:                 coeff = lhs(2,4)
   610:                 lhs(2,5)= lhs(2,5) - coeff*lhs(4,5)
   611:                 r(2)   = r(2)   - coeff*r(4)
   612:           
   613:                 coeff = lhs(3,4)
   614:                 lhs(3,5)= lhs(3,5) - coeff*lhs(4,5)
   615:                 r(3)   = r(3)   - coeff*r(4)
   616:           
   617:                 coeff = lhs(5,4)
   618:                 lhs(5,5)= lhs(5,5) - coeff*lhs(4,5)
   619:                 r(5)   = r(5)   - coeff*r(4)
   620:           
   621:           
   622:                 pivot = 1.00d0/lhs(5,5)
   623:                 r(5)   = r(5)  *pivot
   624:           
   625:                 coeff = lhs(1,5)
   626:                 r(1)   = r(1)   - coeff*r(5)
   627:           
   628:                 coeff = lhs(2,5)
   629:                 r(2)   = r(2)   - coeff*r(5)
   630:           
   631:                 coeff = lhs(3,5)
   632:                 r(3)   = r(3)   - coeff*r(5)
   633:           
   634:                 coeff = lhs(4,5)
   635:                 r(4)   = r(4)   - coeff*r(5)
   636:           
   637:           
   638:                 return
   639:                 end


