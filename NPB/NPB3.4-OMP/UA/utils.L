NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    12: par(1801): Parallel routine generated.: RECIPROCAL$1
    13: par(1803): Parallelized by "do".
    13: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: RECIPROCAL
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine reciprocal (a, n)
     3:           c------------------------------------------------------------------
     4:           c     initialize double precision array a with length of n
     5:           c------------------------------------------------------------------
     6:           
     7:                 implicit none
     8:           
     9:                 integer n, i
    10:                 double precision a(n)
    11:           
    12:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    13: Y------>        do i = 1, n
    14: |                 a(i) = 1.d0/a(i)
    15: Y------         end do
    16:           c$OMP END PARALLEL DO
    17:           
    18:                 return
    19:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    31: par(1801): Parallel routine generated.: R_INIT_OMP$1
    32: par(1803): Parallelized by "do".
    32: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT_OMP
FORMAT LIST

 LINE   LOOP      STATEMENT

    21:                 subroutine r_init_omp (a, n, const)
    22:           c------------------------------------------------------------------
    23:           c     initialize double precision array a with length of n
    24:           c------------------------------------------------------------------
    25:           
    26:                 implicit none
    27:           
    28:                 integer n, i
    29:                 double precision a(n), const
    30:           
    31:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    32: Y------>        do i = 1, n
    33: |                 a(i) = const
    34: Y------         end do
    35:           c$OMP END PARALLEL DO
    36:           
    37:                 return
    38:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    50: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: R_INIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    40:                 subroutine r_init (a, n, const)
    41:           c------------------------------------------------------------------
    42:           c     initialize double precision array a with length of n
    43:           c------------------------------------------------------------------
    44:           
    45:                 implicit none
    46:           
    47:                 integer n, i
    48:                 double precision a(n), const
    49:           
    50: V------>        do i = 1, n
    51: |                 a(i) = const
    52: V------         end do
    53:           
    54:                 return
    55:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    66: par(1801): Parallel routine generated.: NR_INIT_OMP$1
    67: par(1803): Parallelized by "do".
    67: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT_OMP
FORMAT LIST

 LINE   LOOP      STATEMENT

    57:                 subroutine nr_init_omp (a, n, const)
    58:           c------------------------------------------------------------------
    59:           c     initialize integer array a with length of n
    60:           c------------------------------------------------------------------
    61:           
    62:                 implicit none
    63:           
    64:                 integer n, i, a(n), const
    65:           
    66:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    67: Y------>        do i = 1, n
    68: |                 a(i) = const
    69: Y------         end do
    70:           c$OMP END PARALLEL DO
    71:           
    72:                 return
    73:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    85: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NR_INIT
FORMAT LIST

 LINE   LOOP      STATEMENT

    76:                 subroutine nr_init (a, n, const)
    77:           c------------------------------------------------------------------
    78:           c     initialize integer array a with length of n
    79:           c------------------------------------------------------------------
    80:           
    81:                 implicit none
    82:           
    83:                 integer n, i, a(n), const
    84:           
    85: V------>        do i = 1, n
    86: |                 a(i) = const
    87: V------         end do
    88:           
    89:                 return
    90:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   101: par(1801): Parallel routine generated.: L_INIT_OMP$1
   102: par(1803): Parallelized by "do".
   102: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: L_INIT_OMP
FORMAT LIST

 LINE   LOOP      STATEMENT

    92:                 subroutine l_init_omp (a, n, const)
    93:           c------------------------------------------------------------------
    94:           c     initialize integer array a with length of n
    95:           c------------------------------------------------------------------
    96:           
    97:                 implicit none
    98:                 integer n, i
    99:                 logical a(n), const
   100:           
   101:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
   102: Y------>        do i = 1, n
   103: |                 a(i) = const
   104: Y------         end do
   105:           c$OMP END PARALLEL DO
   106:           
   107:                 return
   108:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   121: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NCOPY
FORMAT LIST

 LINE   LOOP      STATEMENT

   111:                 subroutine ncopy (a,b,n)
   112:           c------------------------------------------------------------------
   113:           c     copy array of integers b to a, the length of array is n
   114:           c------------------------------------------------------------------
   115:           
   116:                 implicit none
   117:           
   118:                 integer n,i
   119:                 integer a(n), b(n)
   120:           
   121: V------>        do i = 1, n
   122: |                 a(i) = b(i)
   123: V------         end do
   124:           
   125:                 return
   126:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   139: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: COPY
FORMAT LIST

 LINE   LOOP      STATEMENT

   129:                 subroutine copy (a,b,n)
   130:           c------------------------------------------------------------------
   131:           c     copy double precision array b to a, the length of array is n
   132:           c------------------------------------------------------------------
   133:           
   134:                 implicit none
   135:           
   136:                 integer n,i
   137:                 double precision a(n), b(n)
   138:           
   139: V------>        do i = 1, n
   140: |                  a(i) = b(i)
   141: V------         end do
   142:           
   143:                 return
   144:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   155: par(1801): Parallel routine generated.: ADDS2M1$1
   156: par(1803): Parallelized by "do".
   156: vec( 101): Vectorized loop.
   157: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS2M1
FORMAT LIST

 LINE   LOOP      STATEMENT

   147:                 subroutine adds2m1(a,b,c1,n)
   148:           c-----------------------------------------------------------------
   149:           c     a=b*c1
   150:           c-----------------------------------------------------------------
   151:                 implicit none
   152:           
   153:                 integer n,i
   154:                 double precision a(n),b(n),c1
   155:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   156: Y------>        do i=1,n
   157: |       F         a(i)=a(i)+c1*b(i)
   158: Y------         end do
   159:           c$OMP END PARALLEL DO
   160:           
   161:                 return
   162:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   174: par(1801): Parallel routine generated.: ADDS1M1$1
   175: par(1803): Parallelized by "do".
   175: vec( 101): Vectorized loop.
   176: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADDS1M1
FORMAT LIST

 LINE   LOOP      STATEMENT

   165:                 subroutine adds1m1(a,b,c1,n )
   166:           c-----------------------------------------------------------------
   167:           c     a=c1*a+b
   168:           c-----------------------------------------------------------------
   169:           
   170:                 implicit none
   171:           
   172:                 integer n,i
   173:                 double precision a(n),b(n),c1
   174:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   175: Y------>        do i=1,n
   176: |       F         a(i)=c1*a(i)+b(i)
   177: Y------         end do
   178:           c$OMP END PARALLEL DO
   179:           
   180:                 return
   181:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   194: par(1801): Parallel routine generated.: COL2$1
   195: par(1803): Parallelized by "do".
   195: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: COL2
FORMAT LIST

 LINE   LOOP      STATEMENT

   184:                 subroutine col2(a,b,n)
   185:           c------------------------------------------------------------------
   186:           c     a=a*b
   187:           c------------------------------------------------------------------
   188:           
   189:                 implicit none
   190:           
   191:                 integer n,i
   192:                 double precision a(n),b(n)
   193:           
   194:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   195: Y------>        do i=1,n
   196: |                 a(i)=a(i)*b(i)
   197: Y------         end do
   198:           c$OMP END PARALLEL DO
   199:           
   200:                 return
   201:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   213: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: NRZERO
FORMAT LIST

 LINE   LOOP      STATEMENT

   204:                 subroutine nrzero (na,n)
   205:           c------------------------------------------------------------------
   206:           c     zero out array of integers
   207:           c------------------------------------------------------------------
   208:           
   209:                 implicit none
   210:           
   211:                 integer n,i,na(n)
   212:           
   213: V------>        do i = 1, n
   214: |                 na(i ) = 0
   215: V------         end do
   216:           
   217:                 return
   218:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   230: par(1801): Parallel routine generated.: ADD2$1
   231: par(1803): Parallelized by "do".
   231: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: ADD2
FORMAT LIST

 LINE   LOOP      STATEMENT

   221:                 subroutine add2(a,b,n)
   222:           c------------------------------------------------------------------
   223:           c     a=a+b
   224:           c------------------------------------------------------------------
   225:           
   226:                 implicit none
   227:           
   228:                 integer n,i
   229:                 double precision  a(n),b(n)
   230:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
   231: Y------>        do i=1,n
   232: |                 a(i)=a(i)+b(i)
   233: Y------         end do
   234:           c$OMP END PARALLEL DO
   235:           
   236:                 return
   237:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   252: par(1801): Parallel routine generated.: CALC_NORM$1
   255: par(1803): Parallelized by "do".
   255: vec( 103): Unvectorized loop.
   255: vec( 113): Overhead of loop division is too large.
   258: vec( 101): Vectorized loop.
   259: opt(1418): Constant-length loop is expanded.
   260: opt(1418): Constant-length loop is expanded.
   261: vec( 128): Fused multiply-add operation applied.
   261: vec( 126): Idiom detected.: SUM
   268: par(1807): Critical-construct generated.
   268: par(1809): Barrier synchronization.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: CALC_NORM
FORMAT LIST

 LINE   LOOP      STATEMENT

   240:                 double precision function calc_norm()
   241:           c------------------------------------------------------------------
   242:           c     calculate the integral of ta1 over the whole domain
   243:           c------------------------------------------------------------------
   244:           
   245:                 use ua_data
   246:                 implicit none
   247:           
   248:                 double precision total,ieltotal
   249:                 integer iel,k,j,i,isize
   250:           
   251:                 total=0.d0
   252:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,isize,ieltotal,iel)
   253:           c$OMP& REDUCTION(+:total)
   254:           
   255: P------>        do iel=1,nelt
   256: |                 ieltotal=0.d0
   257: |                 isize=size_e(iel)
   258: |V----->          do k=1,lx1
   259: ||*---->            do j=1,lx1
   260: |||*--->              do i=1,lx1
   261: ||||    F               ieltotal=ieltotal+ta1(i,j,k,iel)*w3m1(i,j,k)
   262: ||||           &                               *jacm1_s(i,j,k,isize)
   263: |||*---               end do
   264: ||*----             end do
   265: |V-----           end do
   266: |               total=total+ieltotal
   267: P------         end do
   268:           c$OMP END PARALLEL DO
   269:           
   270:                 calc_norm = total
   271:           
   272:                 return
   273:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   275: vec( 103): Unvectorized loop.
   275: vec( 108): Unvectorizable loop structure.
   293: vec( 109): Vectorization obstructive statement.
   293: vec( 103): Unvectorized loop.
   295: par(1801): Parallel routine generated.: PARALLEL_ADD$1
   296: par(1803): Parallelized by "do".
   298: vec( 101): Vectorized loop.
   308: par(1801): Parallel routine generated.: PARALLEL_ADD$2
   309: par(1803): Parallelized by "do".
   309: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: PARALLEL_ADD
FORMAT LIST

 LINE   LOOP      STATEMENT

   275:                 subroutine parallel_add(frontier)
   276:           c-----------------------------------------------------------------
   277:           c     input array frontier, perform (potentially) parallel add so that
   278:           c     the output frontier(i) has sum of frontier(1)+frontier(2)+...+frontier(i)
   279:           c-----------------------------------------------------------------
   280:           
   281:                 use ua_data
   282:                 implicit none
   283:           
   284:                 integer nellog,i,ahead,ii,ntemp,n1,ntemp1,frontier(lelt),iel
   285:           
   286:                 nellog=0
   287:                 iel=1
   288:              10 iel=iel*2
   289:                 nellog=nellog+1
   290:                 if (iel.lt.nelt) goto 10
   291:           
   292:                 ntemp=1
   293: +------>        do i=1,nellog
   294: |                 n1=ntemp*2
   295: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ahead,ii,iel)
   296: |P----->          do iel=n1, nelt,n1
   297: ||                  ahead=frontier(iel-ntemp)
   298: ||V---->            do ii=ntemp-1,0,-1
   299: |||                   frontier(iel-ii)=frontier(iel-ii)+ahead
   300: ||V----             end do
   301: |P-----           end do
   302: |         c$OMP END PARALLEL DO
   303: |         
   304: |                 iel=(nelt/n1+1)*n1
   305: |                 ntemp1=iel-nelt
   306: |                 if(ntemp1.lt.ntemp)then
   307: |                   ahead=frontier(iel-ntemp)
   308: |         c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)
   309: |Y----->            do ii=ntemp-1,ntemp1,-1
   310: ||                    frontier(iel-ii)=frontier(iel-ii)+ahead
   311: |Y-----             end do
   312: |         c$OMP END PARALLEL DO
   313: |                 end if
   314: |         
   315: |                 ntemp=n1
   316: +------         end do
   317:           
   318:                 return
   319:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: DSSUM
FORMAT LIST

 LINE   LOOP      STATEMENT

   322:                 subroutine dssum
   323:           
   324:           c------------------------------------------------------------------
   325:           c     Perform stiffness summation: element-mortar-element mapping
   326:           c------------------------------------------------------------------
   327:           
   328:                 use ua_data
   329:                 implicit none
   330:           
   331:                 call transfb(dpcmor,dpcelm)
   332:                 call transf (dpcmor,dpcelm)
   333:           
   334:                 return
   335:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   364: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: utils.f

PROCEDURE NAME: FACEV
FORMAT LIST

 LINE   LOOP      STATEMENT

   338:                 subroutine facev(a,iface,val)
   339:           c------------------------------------------------------------------
   340:           c     assign the value val to face(iface,iel) of array a.
   341:           c------------------------------------------------------------------
   342:           
   343:                 use ua_data
   344:                 implicit none
   345:           
   346:                 double precision a(lx1,lx1,lx1), val
   347:                 integer iface, kx1, kx2, ky1, ky2, kz1, kz2, ix, iy, iz
   348:           
   349:                 kx1=1
   350:                 ky1=1
   351:                 kz1=1
   352:                 kx2=lx1
   353:                 ky2=lx1
   354:                 kz2=lx1
   355:                 if (iface.eq.1) kx1=lx1
   356:                 if (iface.eq.2) kx2=1
   357:                 if (iface.eq.3) ky1=lx1
   358:                 if (iface.eq.4) ky2=1
   359:                 if (iface.eq.5) kz1=lx1
   360:                 if (iface.eq.6) kz2=1
   361:           
   362: +------>        do ix = kx1, kx2
   363: |+----->          do iy = ky1, ky2
   364: ||V---->            do iz = kz1, kz2
   365: |||                   a(ix,iy,iz)=val
   366: ||V----             end do
   367: |+-----           end do
   368: +------         end do
   369:           
   370:                 return
   371:                 end


