NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:28 2020
FILE NAME: adapt.f

PROCEDURE NAME: ADAPTATION
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    13: inl(1212): Source for routine not found.: TIMER_START
    32: inl(1388): Inlining inhibited: OpenMP or parallel directive.: FIND_REFINE
    35: vec( 103): Unvectorized loop.
    35: vec( 110): Vectorization obstructive procedure reference.: CHECK_REFINE
    39: inl(1388): Inlining inhibited: OpenMP or parallel directive.: CHECK_REFINE
    43: inl(1388): Inlining inhibited: OpenMP or parallel directive.: DO_REFINE
    56: par(1801): Parallel routine generated.: ADAPTATION$1
    57: par(1803): Parallelized by "do".
    57: vec( 101): Vectorized loop.
    62: vec( 103): Unvectorized loop.
    62: vec( 110): Vectorization obstructive procedure reference.: FIND_COARSEN
    62: vec( 110): Vectorization obstructive procedure reference.: DO_COARSEN
    68: inl(1388): Inlining inhibited: OpenMP or parallel directive.: FIND_COARSEN
    74: inl(1388): Inlining inhibited: OpenMP or parallel directive.: DO_COARSEN
    91: par(1801): Parallel routine generated.: ADAPTATION$2
    92: par(1803): Parallelized by "do".
    92: vec( 101): Vectorized loop.
   100: inl(1212): Source for routine not found.: MOVE
   105: inl(1212): Source for routine not found.: MORTAR
   106: inl(1212): Source for routine not found.: PREPWORK
   108: inl(1212): Source for routine not found.: TIMER_STOP


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:28 2020
FILE NAME: adapt.f

PROCEDURE NAME: ADAPTATION
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine adaptation (ifmortar,step)
     3:           c-----------------------------------------------------------
     4:           c     For 3-D mesh adaptation (refinement+ coarsening)
     5:           c-----------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:            
    10:                 logical if_coarsen,if_refine,ifmortar,ifrepeat
    11:                 integer iel,miel,irefine,icoarsen,neltold,step
    12:           
    13:                 if (timeron) call timer_start(t_adaptation)
    14:                 ifmortar=.false.
    15:           c.....compute heat source center(x0,y0,z0)
    16:                 x0=x00+velx*time
    17:                 y0=y00+vely*time
    18:                 z0=z00+velz*time
    19:           
    20:           c.....Search elements to be refined. Check with restrictions. Perform
    21:           c     refinement repeatedly until all desired refinements are done.
    22:           
    23:           c.....ich(iel)=0 no grid change on element iel
    24:           c.....ich(iel)=2 iel is marked to be coarsened
    25:           c.....ich(iel)=4 iel is marked to be refined
    26:           
    27:           c.....irefine records how many elements got refined
    28:                 irefine=0
    29:           
    30:           c.....check whether elements need to be refined because they have overlap
    31:           c     with the  heat source
    32:           4     call find_refine(if_refine)
    33:           
    34:                 if(if_refine) then
    35:                   ifrepeat=.true.
    36:           2       if(ifrepeat) then
    37:           c.........Check with restriction, unmark elements that cannot be refined.
    38:           c         Elements preventing desired refinement will be marked to be refined.
    39:                     call check_refine(ifrepeat)
    40:                     go to 2
    41:                   end if
    42:           c.......perform refinement
    43:                   call do_refine(ifmortar,irefine)
    44:                   goto 4
    45:                 endif
    46:           
    47:           c.....Search for elements to be coarsened. Check with restrictions,
    48:           c     Perform coarsening repeatedly until all possible coarsening
    49:           c     is done.
    50:           
    51:           c.....icoarsen records how many elements got coarsened
    52:                 icoarsen=0
    53:           
    54:           c.....skip(iel)=.true. indicates an element no longer exists (because it
    55:           c     got merged)
    56:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel)
    57: Y------>         do iel=1,nelt
    58: |                 skip(iel)=.false.
    59: Y------         end do
    60:           c$OMP END PARALLEL DO
    61:           
    62:                 neltold=nelt
    63:           
    64:           c.....Check whether elements need to be coarsened because they don't have
    65:           c     overlap with the heat source. Only elements that don't have a larger
    66:           c     size neighbor can be marked to be coarsened
    67:           
    68:           5     call find_coarsen(if_coarsen,neltold)
    69:           
    70:                 if(if_coarsen) then
    71:           c.......Perform coarsening, however subject to restriction. Only possible
    72:           c       coarsening will be performed. if_coarsen=.true. indicates that
    73:           c       actual coarsening happened
    74:                   call do_coarsen(if_coarsen,icoarsen,neltold)
    75:                   if(if_coarsen) then
    76:           c.........ifmortar=.true. indicates the grid changed, i.e. the mortar points
    77:           c         indices need to be regenerated on the new grid.
    78:                     ifmortar=.true.
    79:                     go to 5
    80:                   end if
    81:                 end if
    82:           
    83:                 write(*,1000) step, irefine, icoarsen, nelt
    84:            1000 format('Step ',i4, ': elements refined, merged, total: ',
    85:                &       i7, 1X , i7, 1X, i7)
    86:           
    87:           c.....mt_to_id(miel) takes as argument the morton index  and returns the actual
    88:           c                    element index
    89:           c.....id_to_mt(iel)  takes as argument the actual element index and returns the
    90:           c                    morton index
    91:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel,iel)
    92: Y------>        do miel=1,nelt
    93: |                 iel=mt_to_id(miel)
    94: |       C         id_to_mt(iel)=miel
    95: Y------         end do
    96:           c$OMP END PARALLEL DO
    97:           
    98:           c.....Reorder the elements in the order of the morton curve. After the move
    99:           c     subroutine the element indices are  the same as the morton indices
   100:                 call move
   101:           
   102:           c.....if the grid changed, regenerate mortar indices and update variables
   103:           c     associated to grid.
   104:                 if (ifmortar) then
   105:                   call mortar
   106:                   call prepwork
   107:                 endif
   108:                 if (timeron) call timer_stop(t_adaptation)
   109:           
   110:                 return
   111:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:28 2020
FILE NAME: adapt.f

PROCEDURE NAME: DO_COARSEN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   138: par(1801): Parallel routine generated.: DO_COARSEN$1
   140: par(1803): Parallelized by "do".
   140: vec( 101): Vectorized loop.
   146: par(1803): Parallelized by "do".
   146: vec( 101): Vectorized loop.
   155: par(1801): Parallel routine generated.: DO_COARSEN$2
   158: par(1803): Parallelized by "do".
   158: vec( 108): Unvectorizable loop structure.
   186: vec( 101): Vectorized loop.
   193: vec( 101): Vectorized loop.
   200: vec( 101): Vectorized loop.
   201: inl(1222): Inlined: ICHECK
   217: vec( 101): Vectorized loop.
   221: vec( 101): Vectorized loop.
   233: inl(1212): Source for routine not found.: PARALLEL_ADD
   242: par(1801): Parallel routine generated.: DO_COARSEN$3
   243: par(1803): Parallelized by "do".
   243: vec( 101): Vectorized loop.
   258: par(1801): Parallel routine generated.: DO_COARSEN$4
   259: vec( 110): Vectorization obstructive procedure reference.: MERGING
   259: par(1803): Parallelized by "do".
   259: vec( 103): Unvectorized loop.
   272: inl(1214): Expansion routine is too big for automatic expansion.: MERGING


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:28 2020
FILE NAME: adapt.f

PROCEDURE NAME: DO_COARSEN
FORMAT LIST

 LINE   LOOP      STATEMENT

   115:                 subroutine do_coarsen(if_coarsen,icoarsen,neltold)
   116:           c---------------------------------------------------------------
   117:           c     Coarsening procedure:
   118:           c     1) check with restrictions
   119:           c     2) perform coarsening
   120:           c---------------------------------------------------------------
   121:           
   122:                 use ua_data
   123:                 implicit none
   124:           
   125:                 logical if_coarsen, icheck,test,test1,test2,test3
   126:                 integer iel, ntp(8), ntempmin, ic, parent, mielnew, miel,
   127:                &        icoarsen, i, index, num_coarsen, ntemp, ii, ntemp1,
   128:                &        neltold
   129:            
   130:                 if_coarsen=.false.
   131:           
   132:           c.....If an element has been merged, it will be skipped afterwards
   133:           c     skip(iel)=.true. for elements that will be skipped.
   134:           c     ifcoa_id(iel)=.true. indicates that element iel will be coarsened
   135:           c     ifcoa(miel)=.true. refers to element miel(mortar index) will be
   136:           c                        coarsened
   137:           
   138:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iel)
   139:           c$OMP DO
   140: Y------>        do iel=1,nelt
   141: |                 mt_to_id_old(iel)=mt_to_id(iel)
   142: |                 mt_to_id(iel)=0
   143: Y------         end do
   144:           c$OMP END DO nowait
   145:           c$OMP DO
   146: Y------>        do iel=1,neltold
   147: |                 ifcoa_id(iel)=.false.
   148: Y------         end do
   149:           c$OMP END DO nowait
   150:           c$OMP END PARALLEL
   151:           
   152:           c.....Check whether the potential coarsening will make neighbor,
   153:           c     and neighbor's neighbor....break grid restriction
   154:           
   155:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel,iel,ic,
   156:           c$OMP& ntp,parent,test,test1,i,test2,test3)
   157:           c$OMP& SHARED(if_coarsen)
   158: P------>        do miel=1,nelt
   159: |                 ifcoa(miel)=.false.
   160: |                 front(miel)=0
   161: |                 iel=mt_to_id_old(miel)
   162: |         c.......if an element is marked to be coarsened
   163: |                 if(ich(iel).eq.2) then
   164: |         
   165: |         c.........If the current  element is the "first" child (front-left-
   166: |         c         bottom) of its parent (tree(iel) mod 8 equals 0), then
   167: |         c         find all its neighbors. Check whether they are from the same
   168: |         c         parent.
   169: |         
   170: |                   ic=tree(iel)
   171: |                   if(.not.btest(ic,0).and..not.btest(ic,1).and.
   172: |              &       .not.btest(ic,2)) then
   173: |                     ntp(1)=iel
   174: |                     ntp(2)=sje(1,1,1,iel)
   175: |                     ntp(3)=sje(1,1,3,iel)
   176: |                     ntp(4)=sje(1,1,1,ntp(3))
   177: |                     ntp(5)=sje(1,1,5,iel)
   178: |                     ntp(6)=sje(1,1,1,ntp(5))
   179: |                     ntp(7)=sje(1,1,3,ntp(5))
   180: |                     ntp(8)=sje(1,1,1,ntp(7))
   181: |          
   182: |                     parent=ishft(tree(iel),-3)
   183: |                     test=.false.
   184: |         
   185: |                     test1=.true.
   186: |V----->              do i=1,8
   187: ||      G               if(ishft(tree(ntp(i)),-3).ne.parent)test1=.false.
   188: |V-----               end do
   189: |         
   190: |         c...........check whether all child elements are marked to be coarsened
   191: |                     if(test1)then
   192: |                       test2=.true.
   193: |V----->                do i=1,8
   194: ||      G                 if(ich(ntp(i)).ne.2)test2=.false.
   195: |V-----                 end do
   196: |         
   197: |         c.............check whether all child elements can be coarsened or not.
   198: |                       if(test2)then
   199: |                         test3=.true.
   200: |V----->                  do i=1,8
   201: ||      I                   if(.not.icheck(ntp(i),i))test3=.false.
   202: |V-----                   end do
   203: |                         if(test3)test=.true.
   204: |                       end if
   205: |                     end if
   206: |         c...........if the eight child elements are eligible to be coarsened
   207: |         c           mark the first children ifcoa(miel)=.true.
   208: |         c           mark them all ifcoa_id()=.true.
   209: |         c           front(miel) will be used to calculate (potentially in parallel)
   210: |         c                       how many elements with seuqnece numbers less than
   211: |         c                       miel will be coarsened.
   212: |         c           skip()      marks that an element will no longer exist after merge.
   213: |         
   214: |                     if(test)then
   215: |         
   216: |                       ifcoa(miel)=.true.
   217: |V----->                do i=1,8
   218: ||      C                 ifcoa_id(ntp(i))=.true.
   219: |V-----                 end do
   220: |                       front(miel)=1
   221: |V----->                do i=1,7
   222: ||      C                  skip(ntp(i+1))=.true.
   223: |V-----                 end do
   224: |                       if(.not.if_coarsen) if_coarsen=.true.
   225: |                     end if
   226: |                   end if
   227: |                 end if
   228: P------         end do
   229:           c$OMP END PARALLEL DO
   230:           
   231:           c.....compute front(iel), how many elements will be coarsened before iel
   232:           c     (including iel)
   233:                 call parallel_add(front)
   234:           
   235:           c.....num_coarsen is the total number of elements that will be coarsened
   236:                 num_coarsen=front(nelt)
   237:           
   238:           c.....action(i) records the morton index of the i'th element (if it is an
   239:           c     element's front-left-bottom-child) to be coarsened.
   240:           
   241:           c.....create array mt_to_id to convert actual element index to morton index
   242:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel,iel,mielnew)
   243: Y------>        do miel=1,nelt
   244: |                 iel=mt_to_id_old(miel)
   245: |       G         if(.not.skip(iel))then
   246: |                   if(ifcoa(miel))then
   247: |       C             action(front(miel))=miel
   248: |                     mielnew=miel-(front(miel)-1)*7
   249: |                   else
   250: |                     mielnew=miel-front(miel)*7
   251: |                   end if
   252: |       C           mt_to_id(mielnew)=iel
   253: |                 end if
   254: Y------         end do
   255:           c$OMP END PARALLEL DO
   256:           
   257:           c.....perform the coarsening procedure (potentially in parallel)
   258:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(index,miel,iel,ntp)
   259: P------>        do index=1,num_coarsen
   260: |                 miel=action(index)
   261: |                 iel=mt_to_id_old(miel)
   262: |         c.......find eight child elements to be coarsened
   263: |                 ntp(1)=iel
   264: |                 ntp(2)=sje(1,1,1,iel)
   265: |                 ntp(3)=sje(1,1,3,iel)
   266: |                 ntp(4)=sje(1,1,1,ntp(3))
   267: |                 ntp(5)=sje(1,1,5,iel)
   268: |                 ntp(6)=sje(1,1,1,ntp(5))
   269: |                 ntp(7)=sje(1,1,3,ntp(5))
   270: |                 ntp(8)=sje(1,1,1,ntp(7))
   271: |         c.......merge them to be the parent
   272: |                 call merging(ntp)
   273: P------         end do
   274:           c$OMP END PARALLEL DO
   275:                 nelt=nelt-num_coarsen*7
   276:                 icoarsen=icoarsen+num_coarsen*8
   277:           
   278:                 return
   279:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: DO_REFINE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   302: par(1801): Parallel routine generated.: DO_REFINE$1
   303: par(1803): Parallelized by "do".
   303: vec( 101): Vectorized loop.
   317: inl(1212): Source for routine not found.: PARALLEL_ADD
   323: par(1801): Parallel routine generated.: DO_REFINE$2
   324: par(1803): Parallelized by "do".
   324: vec( 101): Vectorized loop.
   335: par(1801): Parallel routine generated.: DO_REFINE$3
   336: par(1803): Parallelized by "do".
   336: vec( 101): Vectorized loop.
   362: par(1801): Parallel routine generated.: DO_REFINE$4
   366: par(1803): Parallelized by "do".
   374: vec( 101): Vectorized loop.
   379: vec( 101): Vectorized loop.
   381: opt(1418): Constant-length loop is expanded.
   382: opt(1418): Constant-length loop is expanded.
   387: inl(1212): Source for routine not found.: COPY
   392: inl(1212): Source for routine not found.: NR_INIT
   395: inl(1212): Source for routine not found.: R_INIT
   399: vec( 103): Unvectorized loop.
   399: vec( 110): Vectorization obstructive procedure reference.: R_INIT
   399: vec( 110): Vectorization obstructive procedure reference.: NR_INIT
   402: opt(1025): Reference to this procedure inhibits optimization.: NR_INIT
   411: vec( 101): Vectorized loop.
   425: opt(1418): Constant-length loop is expanded.
   426: opt(1418): Constant-length loop is expanded.
   432: opt(1418): Constant-length loop is expanded.
   433: opt(1418): Constant-length loop is expanded.
   439: opt(1418): Constant-length loop is expanded.
   444: opt(1418): Constant-length loop is expanded.
   452: opt(1418): Constant-length loop is expanded.
   460: opt(1418): Constant-length loop is expanded.
   461: opt(1418): Constant-length loop is expanded.
   467: opt(1418): Constant-length loop is expanded.
   475: opt(1418): Constant-length loop is expanded.
   479: opt(1418): Constant-length loop is expanded.
   484: opt(1418): Constant-length loop is expanded.
   485: opt(1418): Constant-length loop is expanded.
   490: opt(1418): Constant-length loop is expanded.
   491: opt(1418): Constant-length loop is expanded.
   496: opt(1418): Constant-length loop is expanded.
   509: opt(1033): Potential multiple store conflict -- use directive if OK.: SJE
   509: opt(1033): Potential multiple store conflict -- use directive if OK.: CBC
   509: vec( 103): Unvectorized loop.
   509: vec( 107): Iteration count is too small.
   516: opt(1418): Constant-length loop is expanded.
   527: opt(1418): Constant-length loop is expanded.
   529: opt(1033): Potential multiple store conflict -- use directive if OK.: SJE
   536: opt(1019): Feedback of scalar value from one loop pass to another.: NTEMP
   544: opt(1418): Constant-length loop is expanded.
   570: opt(1418): Constant-length loop is expanded.
   581: opt(1418): Constant-length loop is expanded.
   603: opt(1418): Constant-length loop is expanded.
   609: opt(1062): Too many data dependency problems.
   612: inl(1214): Expansion routine is too big for automatic expansion.: REMAP


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: DO_REFINE
FORMAT LIST

 LINE   LOOP      STATEMENT

   282:                 subroutine do_refine(ifmortar,irefine)
   283:           c-------------------------------------------------------
   284:           c     Refinement procedure
   285:           c--------------------------------------------------------
   286:           
   287:                 use ua_data
   288:                 implicit none
   289:           
   290:                 logical ifmortar
   291:                 double precision xctemp(8), yctemp(8), zctemp(8), xleft, xright,
   292:                &       yleft, yright, zleft, zright, ta1temp(lx1,lx1,lx1),
   293:                &       xhalf, yhalf, zhalf
   294:                 integer iel, i, ii, jj, j, jface,
   295:                &        ntemp, ndir, facedir, k, le(4), ne(4), mielnew,
   296:                &        miel, irefine,ntemp1, num_refine, index, treetemp,
   297:                &        sjetemp(2,2,6), n1, n2, nelttemp,
   298:                &        cb, cbctemp(6)
   299:           
   300:           c.....initialize
   301:           
   302:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel)
   303: Y------>        do miel=1,nelt
   304: |                 mt_to_id_old(miel)=mt_to_id(miel)
   305: |                 mt_to_id(miel)=0
   306: |                 action(miel)=0
   307: |       G         if(ich(mt_to_id_old(miel)).ne.4)then
   308: |                   front(miel)=0
   309: |                 else
   310: |                   front(miel)=1
   311: |                 end if
   312: Y------         end do
   313:           c$OMP END PARALLEL DO
   314:           
   315:           c.....front(iel) records how many elements with sequence numbers less than
   316:           c     or equal to iel will be refined
   317:                 call parallel_add(front)
   318:           
   319:           c.....num_refine is the total number of elements that will be refined
   320:                 num_refine=front(nelt)
   321:           
   322:           c.....action(i) records the morton index of the  i'th element to be refined
   323:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel,iel)
   324: Y------>        do miel=1,nelt
   325: |                 iel=mt_to_id_old(miel)
   326: |       G         if(ich(iel).eq.4)then
   327: |       C           action(front(miel))=miel
   328: |                 end if
   329: Y------         end do
   330:           c$OMP END PARALLEL DO
   331:           
   332:           c.....Compute array mt_to_id to convert the element index to morton index.
   333:           c     ref_front_id(iel) records how many elements with index less than
   334:           c     iel (actual element index, not morton index), will be refined.
   335:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(miel,iel,ntemp,mielnew)
   336: Y------>        do miel=1,nelt
   337: |                 iel=mt_to_id_old(miel)
   338: |       G         if(ich(iel).eq.4)then
   339: |                   ntemp=(front(miel)-1)*7
   340: |                   mielnew=miel+ntemp
   341: |                 else
   342: |                   ntemp=front(miel)*7
   343: |                   mielnew=miel+ntemp
   344: |                 end if
   345: |         
   346: |       C         mt_to_id(mielnew)=iel
   347: |       C         ref_front_id(iel)=nelt+ntemp
   348: Y------         end do
   349:           c$OMP END PARALLEL DO
   350:           
   351:           
   352:           c.....Perform refinement (potentially in parallel):
   353:           c       - Cut an element into eight children.
   354:           c       - Assign them element index  as iel, nelt+1,...., nelt+7.
   355:           c       - Update neighboring information.
   356:           
   357:                 nelttemp=nelt
   358:           
   359:                 if (num_refine .gt. 0) then
   360:                   ifmortar=.true.
   361:                 endif
   362:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(index,miel,mielnew,iel,nelt,
   363:           c$OMP& treetemp,xctemp,yctemp,zctemp,cbctemp,sjetemp,ta1temp,
   364:           c$OMP& ii,jj,ntemp,xleft,xright,xhalf,yleft,yright,yhalf,zleft,zright,
   365:           c$OMP& zhalf,ndir,facedir,jface,cb,le,ne,n1,n2,i,j,k)
   366: P------>        do index=1, num_refine
   367: |         c.......miel is old morton index and mielnew is new morton index after refinement.
   368: |                 miel=action(index)
   369: |                 mielnew=miel+(front(miel)-1)*7
   370: |                 iel=mt_to_id_old(miel)
   371: |                 nelt=nelttemp+(front(miel)-1)*7
   372: |         c.......save iel's information in a temporary array
   373: |                 treetemp=tree(iel)
   374: |V----->          do i=1,8
   375: ||                  xctemp(i)=xc(i,iel)
   376: ||                  yctemp(i)=yc(i,iel)
   377: ||                  zctemp(i)=zc(i,iel)
   378: |V-----           end do
   379: |V----->          do i=1,6
   380: ||                  cbctemp(i)=cbc(i,iel)
   381: ||*---->            do jj=1,2
   382: |||*--->              do ii=1,2
   383: ||||                    sjetemp(ii,jj,i)=sje(ii,jj,i,iel)
   384: |||*---               end do
   385: ||*----             end do
   386: |V-----           end do
   387: |                 call copy(ta1temp,ta1(1,1,1,iel),nxyz)
   388: |         
   389: |         
   390: |         c.......zero out iel here
   391: |                 tree(iel)=0
   392: |                 call nr_init(cbc(1,iel),6,0)
   393: |                 call nr_init(sje(1,1,1,iel),24,0)
   394: |                 call nr_init(ijel(1,1,iel),12,0)
   395: |                 call r_init(ta1(1,1,1,iel),nxyz,0.d0)
   396: |         
   397: |         
   398: |         c.......initialize new child elements:iel and nelt+1~nelt+7
   399: |+----->          do j=1,7
   400: ||                  mt_to_id(mielnew+j)=nelt+j
   401: ||                  tree(nelt+j)=0
   402: ||                  call nr_init(cbc(1,nelt+j),6,0)
   403: ||                  call nr_init(sje(1,1,1,nelt+j),24,0)
   404: ||                  call nr_init(ijel(1,1,nelt+j),12,0)
   405: ||                  call r_init(ta1(1,1,1,nelt+j),nxyz,0.d0)
   406: |+-----           end do
   407: |          
   408: |         c.......update the tree()
   409: |                 ntemp=ishft(treetemp,3)
   410: |                 tree(iel)=ntemp
   411: |V----->          do i=1,7
   412: ||                  tree(nelt+i)=ntemp+mod(i,8)
   413: |V-----           end do
   414: |         c.......update the children's vertices' coordinates
   415: |                 xhalf=xctemp(1)+(xctemp(2)-xctemp(1))/2.d0
   416: |                 xleft=xctemp(1)
   417: |                 xright=xctemp(2)
   418: |                 yhalf=yctemp(1)+(yctemp(3)-yctemp(1))/2.d0
   419: |                 yleft=yctemp(1)
   420: |                 yright=yctemp(3)
   421: |                 zhalf=zctemp(1)+(zctemp(5)-zctemp(1))/2.d0
   422: |                 zleft=zctemp(1)
   423: |                 zright=zctemp(5)
   424: |          
   425: |*----->          do j=1,7,2
   426: ||*---->            do i=1,7,2
   427: |||                   xc(i,nelt+j)     = xhalf
   428: |||                   xc(i+1,nelt+j)   = xright
   429: ||*----             end do
   430: |*-----           end do
   431: |         
   432: |*----->          do j=2,6,2
   433: ||*---->            do i=1,7,2
   434: |||                   xc(i,nelt+j)   = xleft
   435: |||                   xc(i+1,nelt+j) = xhalf
   436: ||*----             end do
   437: |*-----           end do
   438: |          
   439: |*----->          do i=1,7,2
   440: ||                  xc(i,iel)=xleft
   441: ||                  xc(i+1,iel)=xhalf
   442: |*-----           end do
   443: |         
   444: |*----->          do i=1,2
   445: ||                  yc(i,nelt+1)=yleft
   446: ||                  yc(i,nelt+4)=yleft
   447: ||                  yc(i,nelt+5)=yleft
   448: ||                  yc(i+4,nelt+1)=yleft
   449: ||                  yc(i+4,nelt+4)=yleft
   450: ||                  yc(i+4,nelt+5)=yleft
   451: |*-----           enddo
   452: |*----->          do i=3,4
   453: ||                  yc(i,nelt+1)=yhalf
   454: ||                  yc(i,nelt+4)=yhalf
   455: ||                  yc(i,nelt+5)=yhalf
   456: ||                  yc(i+4,nelt+1)=yhalf
   457: ||                  yc(i+4,nelt+4)=yhalf
   458: ||                  yc(i+4,nelt+5)=yhalf
   459: |*-----           end do
   460: |*----->          do j=2,3
   461: ||*---->            do i=1,2
   462: |||                   yc(i,nelt+j)=yhalf
   463: |||                   yc(i,nelt+j+4)=yhalf
   464: |||                   yc(i+4,nelt+j)=yhalf
   465: |||                   yc(i+4,nelt+j+4)=yhalf
   466: ||*----             end do
   467: ||*---->            do i=3,4
   468: |||                   yc(i,nelt+j)=yright
   469: |||                   yc(i,nelt+j+4)=yright
   470: |||                   yc(i+4,nelt+j)=yright
   471: |||                   yc(i+4,nelt+j+4)=yright
   472: ||*----             end do
   473: |*-----           end do
   474: |          
   475: |*----->          do i=1,2
   476: ||                  yc(i,iel)=yleft
   477: ||                  yc(i+4,iel)=yleft
   478: |*-----           end do
   479: |*----->          do i=3,4
   480: ||                  yc(i,iel)=yhalf
   481: ||                  yc(i+4,iel)=yhalf
   482: |*-----           end do
   483: |         
   484: |*----->          do j=1,3
   485: ||*---->            do i=1,4
   486: |||                   zc(i,nelt+j)=zleft
   487: |||                   zc(i+4,nelt+j)=zhalf
   488: ||*----             end do
   489: |*-----           end do
   490: |*----->          do j=4,7
   491: ||*---->            do i=1,4
   492: |||                   zc(i,nelt+j)=zhalf
   493: |||                   zc(i+4,nelt+j)=zright
   494: ||*----             end do
   495: |*-----           end do
   496: |*----->          do i=1,4
   497: ||                  zc(i,iel)=zleft
   498: ||                  zc(i+4,iel)=zhalf
   499: |*-----           end do
   500: |         
   501: |         c.......update the children's neighbor information
   502: |         
   503: |         c.......ndir refers to the x,y,z directions, respectively.
   504: |         c       facedir refers to the orientation of the face in each direction,
   505: |         c       e.g. ndir=1, facedir=0 refers to face 1,
   506: |         c       and ndir =1, facedir=1 refers to face 2.
   507: |         
   508: |+----->          do ndir = 1, 3
   509: ||+---->            do facedir = 0, 1
   510: |||                   i=2*ndir-1+facedir
   511: |||                   jface=jjface(i)
   512: |||                   cb=cbctemp(i)
   513: |||       
   514: |||       c...........find the new element indices of the four children on each
   515: |||       c           face of the parent element
   516: |||*--->              do k = 1, 4
   517: ||||                    le(k) = le_arr(k,facedir,ndir)+nelt
   518: ||||                    ne(k) = le_arr(k,1-facedir,ndir)+nelt
   519: |||*---               end do
   520: |||                   if(facedir.eq.0)then
   521: |||                     le(1)=iel
   522: |||                   else
   523: |||                     ne(1)=iel
   524: |||                   end if
   525: |||       c...........update neighbor information of the four child elements on each
   526: |||       c           face of the parent element
   527: |||*--->              do k=1,4
   528: ||||                    cbc(i,le(k))=2
   529: ||||                    sje(1,1,i,le(k))=ne(k)
   530: ||||                    ijel(1,i,le(k))=1
   531: ||||                    ijel(2,i,le(k))=1
   532: |||*---               end do
   533: |||       
   534: |||       c...........if the face type of the parent element is type 2
   535: |||                   if(cb.eq.2) then
   536: |||                     ntemp=sjetemp(1,1,i)
   537: |||       
   538: |||       c.............if the neighbor ntemp is not marked to be refined
   539: |||                     if(ich(ntemp).ne.4)then
   540: |||                       cbc(jface,ntemp)=3
   541: |||                       ijel(1,jface,ntemp)=1
   542: |||                       ijel(2,jface,ntemp)=1
   543: |||        
   544: |||*--->                  do k=1,4
   545: ||||                        cbc(i,ne(k))=1
   546: ||||                        sje(1,1,i,ne(k))=ntemp
   547: ||||                        if(k.eq.1) then
   548: ||||                          ijel(1,i,ne(k))=1
   549: ||||                          ijel(2,i,ne(k))=1
   550: ||||                          sje(1,1,jface,ntemp)=ne(k)
   551: ||||                        elseif(k.eq.2) then
   552: ||||                          ijel(1,i,ne(k))=1
   553: ||||                          ijel(2,i,ne(k))=2
   554: ||||                          sje(1,2,jface,ntemp)=ne(k)
   555: ||||                        elseif(k.eq.3) then
   556: ||||                          ijel(1,i,ne(k))=2
   557: ||||                          ijel(2,i,ne(k))=1
   558: ||||                          sje(2,1,jface,ntemp)=ne(k)
   559: ||||                        elseif(k.eq.4) then
   560: ||||                          ijel(1,i,ne(k))=2
   561: ||||                          ijel(2,i,ne(k))=2
   562: ||||                          sje(2,2,jface,ntemp)=ne(k)
   563: ||||                        end if
   564: |||*---                   end do
   565: |||       
   566: |||       c.............if the neighbor ntemp is also marked to be refined
   567: |||                     else
   568: |||                       n1=ref_front_id(ntemp)
   569: |||        
   570: |||*--->                  do k=1,4
   571: ||||                        cbc(i,ne(k))=2
   572: ||||                        n2=n1+le_arr(k,facedir,ndir)
   573: ||||                        if(n2.eq.n1+8)n2=ntemp
   574: ||||                        sje(1,1,i,ne(k))=n2
   575: ||||                        ijel(1,i,ne(k))=1
   576: |||*---                   end do
   577: |||       
   578: |||                     endif
   579: |||       c...........if the face type of the parent element is type 3
   580: |||                   elseif(cb.eq.3) then
   581: |||*--->                do k=1,4
   582: ||||                      cbc(i,ne(k))=2
   583: ||||                      if(k.eq.1) then
   584: ||||                        ntemp=sjetemp(1,1,i)
   585: ||||                      elseif (k.eq.2) then
   586: ||||                        ntemp=sjetemp(1,2,i)
   587: ||||                      elseif(k.eq.3) then
   588: ||||                        ntemp=sjetemp(2,1,i)
   589: ||||                      elseif(k.eq.4) then
   590: ||||                        ntemp=sjetemp(2,2,i)
   591: ||||                      end if
   592: ||||                      ijel(1,i,ne(k))=1
   593: ||||                      ijel(2,i,ne(k))=1
   594: ||||                      sje(1,1,i,ne(k))=ntemp
   595: ||||                      cbc(jface,ntemp)=2
   596: ||||                      sje(1,1,jface,ntemp)=ne(k)
   597: ||||                      ijel(1,jface,ntemp)=1
   598: ||||                      ijel(2,jface,ntemp)=1
   599: |||*---                 end do
   600: |||       
   601: |||       c...........if the face type of the parent element is type 0
   602: |||                   elseif(cb.eq.0) then
   603: |||*--->                do k=1,4
   604: ||||                      cbc(i,ne(k))=cb
   605: |||*---                 end do
   606: |||                   end if
   607: |||       
   608: ||+----             end do
   609: |+-----           end do
   610: |         
   611: |         c.......map solution from parent element to children
   612: |                 call remap(ta1(1,1,1,iel),ta1(1,1,1,ref_front_id(iel)+1),
   613: |              &             ta1temp(1,1,1))
   614: P------         end do
   615:           c$OMP ENDPARALLEL DO
   616:           
   617:                 nelt=nelttemp+num_refine*7
   618:                 irefine=irefine+num_refine
   619:                 ntot=nelt*lx1*lx1*lx1
   620:                 return
   621:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: IFCOR
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   639: inl(1212): Source for routine not found.: IFSAME
   639: opt(1025): Reference to this procedure inhibits optimization.: IFSAME


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: IFCOR
FORMAT LIST

 LINE   LOOP      STATEMENT

   624:                  logical function ifcor(n1,n2,i,iface)
   625:           c-----------------------------------------------------------
   626:           c      returns whether element n1's face i and element n2's
   627:           c      jjface(iface) have intersections, i.e. whether n1 and
   628:           c      n2 are neighbored by an edge.
   629:           c-----------------------------------------------------------
   630:           
   631:                  use ua_data
   632:                  implicit none
   633:           
   634:                  integer n1,n2,i,iface
   635:                  logical ifsame
   636:           
   637:                  ifcor=.false.
   638:           
   639:                  if(ifsame(n1,e1v1(iface,i),n2,e2v1(iface,i)).or.
   640:                &    ifsame(n1,e1v2(iface,i),n2,e2v2(iface,i))) then
   641:                     ifcor=.true.
   642:                  end if
   643:           
   644:                  return
   645:                  end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: ICHECK
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: ICHECK
FORMAT LIST

 LINE   LOOP      STATEMENT

   648:                 logical function icheck(ie,n)
   649:           c-----------------------------------------------------------
   650:           c     Check whether element ie's three faces (sharing vertex n)
   651:           c     are nonconforming. This will prevent it from being coarsened.
   652:           c     Also check ie's neighbors on those three faces, whether ie's
   653:           c     neighbors by only an edge have a size smaller than ie's,
   654:           c     which also prevents ie from being coarsened.
   655:           c-----------------------------------------------------------
   656:           
   657:                 use ua_data
   658:                 implicit none
   659:           
   660:                 integer ie, n, iside, ntemp1, ntemp2, ntemp3, n1, n2, n3,
   661:                &cb2_1,cb3_1,cb1_2,cb3_2,cb1_3,cb2_3
   662:           
   663:                 icheck=.true.
   664:                 cb2_1=0
   665:                 cb3_1=0
   666:                 cb1_2=0
   667:                 cb3_2=0
   668:                 cb1_3=0
   669:                 cb2_3=0
   670:           
   671:                 n1=f_c(1,n)
   672:                 n2=f_c(2,n)
   673:                 n3=f_c(3,n)
   674:         G       if((cbc(n1,ie).eq.3) .or. (cbc(n2,ie).eq.3) .or.
   675:                &   (cbc(n3,ie).eq.3)) then
   676:                    icheck=.false.
   677:                 else
   678:         G         ntemp1=sje(1,1,n1,ie)
   679:         G         ntemp2=sje(1,1,n2,ie)
   680:         G         ntemp3=sje(1,1,n3,ie)
   681:                   if(ntemp1.ne.0)then
   682:         G            cb2_1=cbc(n2,ntemp1)
   683:         G            cb3_1=cbc(n3,ntemp1)
   684:                   end if
   685:                   if(ntemp2.ne.0)then
   686:         G            cb3_2=cbc(n3,ntemp2)
   687:         G            cb1_2=cbc(n1,ntemp2)
   688:                   end if
   689:                   if(ntemp3.ne.0)then
   690:         G            cb1_3=cbc(n1,ntemp3)
   691:         G            cb2_3=cbc(n2,ntemp3)
   692:                   end if
   693:         G         if((cbc(n1,ie).eq.2.and.(cb2_1.eq.3.or.
   694:                &                               cb3_1.eq.3)).or.
   695:                &     (cbc(n2,ie).eq.2.and.(cb3_2.eq.3.or.
   696:                &                               cb1_2.eq.3)).or.
   697:                &     (cbc(n3,ie).eq.2.and.(cb1_3.eq.3.or.
   698:                &                              cb2_3.eq.3)))then
   699:                     icheck=.false.
   700:                   end if
   701:                 end if
   702:           
   703:                 return
   704:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: FIND_COARSEN
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   722: par(1801): Parallel routine generated.: FIND_COARSEN$1
   724: par(1803): Parallelized by "do".
   724: vec( 108): Unvectorizable loop structure.
   727: inl(1222): Inlined: IFTOUCH
   729: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: FIND_COARSEN
FORMAT LIST

 LINE   LOOP      STATEMENT

   707:                 subroutine find_coarsen(if_coarsen,neltold)
   708:           c-----------------------------------------------------------
   709:           c     Search elements to be coarsened. Check with restrictions.
   710:           c     This subroutine only checks the element itself, not its
   711:           c     neighbors.
   712:           c-----------------------------------------------------------
   713:           
   714:                 use ua_data
   715:                 implicit none
   716:           
   717:                 logical if_coarsen, iftemp, iftouch
   718:                 integer iel,i,neltold
   719:           
   720:                 if_coarsen=.false.
   721:           
   722:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,i,iftemp)
   723:           c$OMP& SHARED(if_coarsen)
   724: P------>        do iel=1,neltold
   725: |                 if(.not.skip(iel))then
   726: |                   ich(iel)=0
   727: |       I           if(.not.iftouch(iel)) then
   728: |                     iftemp=.false.
   729: |V----->              do i=1,nsides
   730: ||        c.............if iel has a larger size than its face neighbors, it
   731: ||        c             can not be coarsened
   732: ||                      if(cbc(i,iel).eq.3) then
   733: ||                        iftemp=.true.
   734: ||                      endif
   735: |V-----               enddo
   736: |                     if(.not.iftemp) then
   737: |                       if(.not.if_coarsen) if_coarsen=.true.
   738: |                       ich(iel)=2
   739: |                     end if
   740: |                   end if
   741: |                 endif
   742: P------         enddo
   743:           c$OMP END PARALLEL DO
   744:           
   745:                 return
   746:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: FIND_REFINE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   763: par(1801): Parallel routine generated.: FIND_REFINE$1
   765: par(1803): Parallelized by "do".
   765: vec( 103): Unvectorized loop.
   765: vec( 108): Unvectorizable loop structure.
   767: inl(1222): Inlined: IFTOUCH
   769: vec( 120): Unvectorizable dependency.: IF_REFINE$374


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: FIND_REFINE
FORMAT LIST

 LINE   LOOP      STATEMENT

   749:                 subroutine find_refine(if_refine)
   750:           c-----------------------------------------------------------
   751:           c     search elements to be refined based on whether they
   752:           c     have overlap with the heat source
   753:           c-----------------------------------------------------------
   754:           
   755:                 use ua_data
   756:                 implicit none
   757:           
   758:                 logical if_refine, iftouch
   759:                 integer iel
   760:           
   761:                 if_refine=.false.
   762:           
   763:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel)
   764:           c$OMP& SHARED(if_refine)
   765: P------>        do iel=1,nelt
   766: |                 ich(iel)=0
   767: |       I         if(iftouch(iel)) then
   768: |                   if((xc(2,iel)-xc(1,iel)).gt.dlmin) then
   769: |                     if(.not.if_refine) if_refine=.true.
   770: |                     ich(iel)=4
   771: |                   end if
   772: |                 end if
   773: P------         enddo
   774:           c$OMP END PARALLEL DO
   775:           
   776:                 return
   777:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: CHECK_REFINE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   796: par(1801): Parallel routine generated.: CHECK_REFINE$1
   798: par(1803): Parallelized by "do".
   816: vec( 103): Unvectorized loop.
   816: vec( 110): Vectorization obstructive procedure reference.: IFSAME
   822: inl(1222): Inlined: IFCOR
   831: vec( 103): Unvectorized loop.
   831: vec( 113): Overhead of loop division is too large.
   835: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: ICH
   837: vec( 120): Unvectorizable dependency.: IFREPEAT$420


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: CHECK_REFINE
FORMAT LIST

 LINE   LOOP      STATEMENT

   780:                 subroutine check_refine(ifrepeat)
   781:           c-----------------------------------------------------------------
   782:           c     Check whether the potential refinement will violate the
   783:           c     restriction. If so, mark the neighbor and unmark the
   784:           c     original element, and set ifrepeat true. i.e. this procedure
   785:           c     needs to be repeated until no further check is needed
   786:           c-----------------------------------------------------------------
   787:           
   788:                 use ua_data
   789:                 implicit none
   790:            
   791:                 logical ifrepeat,ifcor
   792:                 integer iel,iface,ntemp,nntemp,i,jface
   793:           
   794:                 ifrepeat=.false.
   795:           
   796:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iel,i,jface,ntemp,
   797:           c$OMP& iface,nntemp) SHARED(ifrepeat)
   798: P------>        do iel=1,nelt
   799: |         c.......if iel is marked to be refined
   800: |                 if(ich(iel).eq.4) then
   801: |         c.........check its six faces
   802: |+----->            do i=1,nsides
   803: ||                    jface=jjface(i)
   804: ||                    ntemp=sje(1,1,i,iel)
   805: ||        c...........if one face neighbor is larger in size than iel
   806: ||                    if(cbc(i,iel).eq.1) then
   807: ||        c.............unmark iel
   808: ||                      ich(iel)=0
   809: ||        c.............the large size neighbor ntemp is marked to be refined
   810: ||                      if(ich(ntemp).ne.4) then
   811: ||                        if(.not.ifrepeat) ifrepeat=.true.
   812: ||                        ich(ntemp)=4
   813: ||                      end if
   814: ||        c.............check iel's neighbor, neighbored by an edge on face i, which
   815: ||        c             must be a face neighbor of ntemp
   816: ||+---->                do iface=1,nsides
   817: |||                       if(iface.ne.i.and.iface.ne.jface) then
   818: |||       c................if edge neighbors are larger than iel, mark them to be refined
   819: |||                         if(cbc(iface,ntemp).eq.2) then
   820: |||                           nntemp=sje(1,1,iface,ntemp)
   821: |||       c..................ifcor is to make sure the edge neighbor exist
   822: |||     I                     if(ich(nntemp).ne.4.and.
   823: |||            &                 ifcor(iel,nntemp,i,iface))then
   824: |||                             ich(nntemp)=4
   825: |||                           end if
   826: |||                         end if
   827: |||                       end if
   828: ||+----                 end do
   829: ||        c...........if face neighbor are of the same size of iel, check edge neighbors
   830: ||                    elseif(cbc(i,iel).eq.2)then
   831: ||+---->                do iface=1,nsides
   832: |||                       if(iface.ne.i.and.iface.ne.jface) then
   833: |||                         if(cbc(iface,ntemp).eq.1)then
   834: |||                           nntemp=sje(1,1,iface,ntemp)
   835: |||                           ich(nntemp)=4
   836: |||                           ich(iel)=0
   837: |||                           if(.not.ifrepeat) ifrepeat=.true.
   838: |||                         end if
   839: |||                       end if
   840: ||+----                 end do
   841: ||                    end if
   842: |+-----             enddo
   843: |                 end if
   844: P------         end do
   845:           c$OMP END PARALLEL DO
   846:           
   847:                 return
   848:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: IFTOUCH
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: IFTOUCH
FORMAT LIST

 LINE   LOOP      STATEMENT

   851:                 logical function iftouch(iel)
   852:           c-----------------------------------------------------------------
   853:           c     check whether element iel has overlap with the heat source
   854:           c-----------------------------------------------------------------
   855:           
   856:                 use ua_data
   857:                 implicit none
   858:           
   859:                 double precision dis, dis1, dis2, dis3, alpha2
   860:                 integer iel
   861:           
   862:                 alpha2 = alpha*alpha
   863:           
   864:                 if     (x0 .lt. xc(1,iel)) then
   865:                   dis1 = xc(1,iel) - x0
   866:                 elseif (x0 .gt. xc(2,iel)) then
   867:                   dis1 = x0 - xc(2,iel)
   868:                 else
   869:                   dis1 = 0.d0
   870:                 endif
   871:           
   872:                 if     (y0 .lt. yc(1,iel)) then
   873:                   dis2 = yc(1,iel) - y0
   874:                 elseif (y0 .gt. yc(3,iel)) then
   875:                   dis2 = y0 - yc(3,iel)
   876:                 else
   877:                   dis2 = 0.d0
   878:                 endif
   879:           
   880:                 if     (z0 .lt. zc(1,iel)) then
   881:                   dis3 = zc(1,iel) - z0
   882:                 elseif (z0 .gt. zc(5,iel)) then
   883:                   dis3 = z0 - zc(5,iel)
   884:                 else
   885:                  dis3 = 0.d0
   886:                 endif
   887:           
   888:                 dis = dis1**2+dis2**2+dis3**2
   889:           
   890:                 if (dis .lt. alpha2) then
   891:                  iftouch=.true.
   892:                 else
   893:                  iftouch=.false.
   894:                 end if
   895:           
   896:                 return
   897:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAP
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   916: inl(1212): Source for routine not found.: R_INIT
   923: vec( 101): Vectorized loop.
   924: opt(1418): Constant-length loop is expanded.
   925: vec( 128): Fused multiply-add operation applied.
   926: vec( 128): Fused multiply-add operation applied.
   932: vec( 101): Vectorized loop.
   933: opt(1418): Constant-length loop is expanded.
   934: vec( 128): Fused multiply-add operation applied.
   936: vec( 128): Fused multiply-add operation applied.
   938: vec( 128): Fused multiply-add operation applied.
   940: vec( 128): Fused multiply-add operation applied.
   950: vec( 101): Vectorized loop.
   951: opt(1800): Idiom detected (matrix multiply).
   953: vec( 128): Fused multiply-add operation applied.
   955: vec( 128): Fused multiply-add operation applied.
   957: vec( 128): Fused multiply-add operation applied.
   959: vec( 128): Fused multiply-add operation applied.
   961: vec( 128): Fused multiply-add operation applied.
   963: vec( 128): Fused multiply-add operation applied.
   965: opt(1800): Idiom detected (matrix multiply).


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAP
FORMAT LIST

 LINE   LOOP      STATEMENT

   901:                 subroutine remap (y,y1,x)
   902:           c-----------------------------------------------------------------
   903:           c     After a refinement, map the solution  from the parent (x) to
   904:           c     the eight children. y is the solution on the first child
   905:           c     (front-bottom-left) and y1 is the solution on the next 7
   906:           c     children.
   907:           c-----------------------------------------------------------------
   908:           
   909:                 use ua_data
   910:                 implicit none
   911:           
   912:                 double precision x(lx1,lx1,lx1),y(lx1,lx1,lx1),y1(lx1,lx1,lx1,7),
   913:                &       yone(lx1,lx1,lx1,2), ytwo(lx1,lx1,lx1,4)
   914:                 integer i, iz, ii, jj, kk
   915:           
   916:                 call r_init(y,lx1*lx1*lx1,0.d0)
   917:                 call r_init(y1,lx1*lx1*lx1*7,0.d0)
   918:                 call r_init(yone,lx1*lx1*lx1*2,0.d0)
   919:                 call r_init(ytwo,lx1*lx1*lx1*4,0.d0)
   920:           
   921: +------>        do  i=1,lx1
   922: |+----->          do kk = 1, lx1
   923: ||V---->            do jj = 1, lx1
   924: |||*--->              do ii = 1, lx1
   925: ||||    F               yone(ii,jj,i,1) = yone(ii,jj,i,1) +ixmc1(ii,kk)*x(kk,jj,i)
   926: ||||    F               yone(ii,jj,i,2) = yone(ii,jj,i,2) +ixmc2(ii,kk)*x(kk,jj,i)
   927: |||*---               end do
   928: ||V----             end do
   929: |+-----           end do
   930: |         
   931: |+----->          do kk = 1, lx1
   932: ||V---->            do jj = 1, lx1
   933: |||*--->              do ii = 1, lx1
   934: ||||    F               ytwo(ii,i,jj,1) = ytwo(ii,i,jj,1) +
   935: ||||           &                          yone(ii,kk,i,1)*ixtmc1(kk,jj)
   936: ||||    F               ytwo(ii,i,jj,2) = ytwo(ii,i,jj,2) +
   937: ||||           &                          yone(ii,kk,i,1)*ixtmc2(kk,jj)
   938: ||||    F               ytwo(ii,i,jj,3) = ytwo(ii,i,jj,3) +
   939: ||||           &                          yone(ii,kk,i,2)*ixtmc1(kk,jj)
   940: ||||    F               ytwo(ii,i,jj,4) = ytwo(ii,i,jj,4) +
   941: ||||           &                          yone(ii,kk,i,2)*ixtmc2(kk,jj)
   942: |||*---               end do
   943: ||V----             end do
   944: |+-----           end do
   945: +------         end do
   946:           
   947: +------>        do  iz=1,lx1
   948: |+----->          do kk = 1, lx1
   949: ||+---->            do jj = 1, lx1
   950: |||V--->              do ii = 1, lx1
   951: ||||    M               y(ii,iz,jj) = y(ii,iz,jj) +
   952: ||||           &                        ytwo(ii,kk,iz,1)*ixtmc1(kk,jj)
   953: ||||    F               y1(ii,iz,jj,1) = y1(ii,iz,jj,1) +
   954: ||||           &                        ytwo(ii,kk,iz,3)*ixtmc1(kk,jj)
   955: ||||    F               y1(ii,iz,jj,2) = y1(ii,iz,jj,2) +
   956: ||||           &                        ytwo(ii,kk,iz,2)*ixtmc1(kk,jj)
   957: ||||    F               y1(ii,iz,jj,3) = y1(ii,iz,jj,3) +
   958: ||||           &                        ytwo(ii,kk,iz,4)*ixtmc1(kk,jj)
   959: ||||    F               y1(ii,iz,jj,4) = y1(ii,iz,jj,4) +
   960: ||||           &                        ytwo(ii,kk,iz,1)*ixtmc2(kk,jj)
   961: ||||    F               y1(ii,iz,jj,5) = y1(ii,iz,jj,5) +
   962: ||||           &                        ytwo(ii,kk,iz,3)*ixtmc2(kk,jj)
   963: ||||    F               y1(ii,iz,jj,6) = y1(ii,iz,jj,6) +
   964: ||||           &                        ytwo(ii,kk,iz,2)*ixtmc2(kk,jj)
   965: ||||    M               y1(ii,iz,jj,7) = y1(ii,iz,jj,7) +
   966: ||||           &                        ytwo(ii,kk,iz,4)*ixtmc2(kk,jj)
   967: |||V---               end do
   968: ||+----             end do
   969: |+-----           end do
   970: +------         end do
   971:           
   972:                 return
   973:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: MERGING
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1003: opt(1418): Constant-length loop is expanded.
  1006: opt(1418): Constant-length loop is expanded.
  1009: opt(1418): Constant-length loop is expanded.
  1013: opt(1418): Constant-length loop is expanded.
  1017: opt(1418): Constant-length loop is expanded.
  1020: opt(1418): Constant-length loop is expanded.
  1025: vec( 102): Partially vectorized loop.
  1028: opt(1418): Constant-length loop is expanded.
  1032: opt(1036): Potential feedback - use directive or compiler option if OK.: CBC
  1036: opt(1036): Potential feedback - use directive or compiler option if OK.: SJE
  1042: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1045: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1049: opt(1418): Constant-length loop is expanded.
  1051: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1052: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1053: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1054: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1055: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1056: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1057: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1059: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1065: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1066: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1067: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1068: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1069: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1070: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1071: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1073: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1074: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1075: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: IJEL
  1076: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1079: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: CBC
  1080: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1081: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1082: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1083: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: SJE
  1086: opt(1062): Too many data dependency problems.
  1089: inl(1222): Inlined: REMAP2


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: MERGING
FORMAT LIST

 LINE   LOOP      STATEMENT

   977:                 subroutine merging(iela)
   978:           c-----------------------------------------------------------------------
   979:           c     This subroutine is to merge the eight child elements and map
   980:           c     the solution from eight children to the  merged element.
   981:           c     iela array records the eight elements to be merged.
   982:           c-----------------------------------------------------------------------
   983:           
   984:                 use ua_data
   985:                 implicit none
   986:           
   987:                 double precision x1,x2,y1,y2,z1,z2
   988:                 integer ielnew,i,ntemp,jface,ii,cb,ntempa(4),iela(8),ielold,
   989:                &        ntema(4)
   990:           
   991:                 ielnew=iela(1)
   992:           
   993:                 tree(ielnew)=ishft(tree(ielnew),-3)
   994:           
   995:           c.....element vertices
   996:                 x1=xc(1,iela(1))
   997:                 x2=xc(2,iela(2))
   998:                 y1=yc(1,iela(1))
   999:                 y2=yc(3,iela(3))
  1000:                 z1=zc(1,iela(1))
  1001:                 z2=zc(5,iela(5))
  1002:           
  1003: *------>        do i=1,7,2
  1004: |                 xc(i,ielnew)=x1
  1005: *------         end do
  1006: *------>        do i=2,8,2
  1007: |                 xc(i,ielnew)=x2
  1008: *------         end do
  1009: *------>        do i=1,2
  1010: |                 yc(i,ielnew)=y1
  1011: |                 yc(i+4,ielnew)=y1
  1012: *------         end do
  1013: *------>        do i=3,4
  1014: |                 yc(i,ielnew)=y2
  1015: |                 yc(i+4,ielnew)=y2
  1016: *------         end do
  1017: *------>        do i=1,4
  1018: |                 zc(i,ielnew)=z1
  1019: *------         end do
  1020: *------>        do i=5,8
  1021: |                 zc(i,ielnew)=z2
  1022: *------         end do
  1023:           
  1024:           c.....update neighboring information
  1025: S------>        do i=1,nsides
  1026: |                 jface=jjface(i)
  1027: |       G         ielold=iela(children(1,i))
  1028: |*----->          do ii=1,4
  1029: ||      G           ntempa(ii)=iela(children(ii,i))
  1030: |*-----           end do
  1031: |         
  1032: |                 cb=cbc(i,ielold)
  1033: |          
  1034: |                 if (cb.eq.2) then
  1035: |         c.........if the neighbor elements also will be coarsened
  1036: |                   if(ifcoa_id(sje(1,1,i,ielold)))then
  1037: |                     if (i.eq.2 .or. i.eq. 4 .or. i.eq.6) then
  1038: |                       ntemp=sje(1,1,i,sje(1,1,i,ntempa(1)))
  1039: |                     else
  1040: |                       ntemp=sje(1,1,i,ntempa(1))
  1041: |                     end if
  1042: |                     sje(1,1,i,ielnew)=ntemp
  1043: |                     ijel(1,i,ielnew)=1
  1044: |                     ijel(2,i,ielnew)=1
  1045: |                     cbc(i,ielnew)=2
  1046: |         
  1047: |         c.........if the neighbor elements will not be coarsened
  1048: |                   else
  1049: |*----->              do ii=1,4
  1050: ||                      ntema(ii)=sje(1,1,i,ntempa(ii))
  1051: ||                      cbc(jface,ntema(ii))=1
  1052: ||                      sje(1,1,jface,ntema(ii))=ielnew
  1053: ||                      ijel(1,jface,ntema(ii))=iijj(1,ii)
  1054: ||                      ijel(2,jface,ntema(ii))=iijj(2,ii)
  1055: ||                      sje(iijj(1,ii),iijj(2,ii),i,ielnew)=ntema(ii)
  1056: ||                      ijel(1,i,ielnew)=1
  1057: ||                      ijel(2,i,ielnew)=1
  1058: |*-----               end do
  1059: |                     cbc(i,ielnew)=3
  1060: |                   end if
  1061: |         
  1062: |                 else if(cb.eq.1)then
  1063: |         
  1064: |                   ntemp=sje(1,1,i,ielold)
  1065: |                   cbc(jface,ntemp)=2
  1066: |                   ijel(1,jface,ntemp)=1
  1067: |                   ijel(2,jface,ntemp)=1
  1068: |                   sje(1,1,jface,ntemp)=ielnew
  1069: |                   sje(1,2,jface,ntemp)=0
  1070: |                   sje(2,1,jface,ntemp)=0
  1071: |                   sje(2,2,jface,ntemp)=0
  1072: |          
  1073: |                   cbc(i,ielnew)=2
  1074: |                   ijel(1,i,ielnew)=1
  1075: |                   ijel(2,i,ielnew)=1
  1076: |                   sje(1,1,i,ielnew)=ntemp
  1077: |          
  1078: |                 else if(cb.eq.0)then
  1079: |                   cbc(i,ielnew)=0
  1080: |                   sje(1,1,i,ielnew)=0
  1081: |                   sje(1,2,i,ielnew)=0
  1082: |                   sje(2,1,i,ielnew)=0
  1083: |                   sje(2,2,i,ielnew)=0
  1084: |                 endif
  1085: |         
  1086: S------         end do
  1087:           
  1088:           c.....map solution from children to the merged element
  1089:         I       call remap2(iela, ielnew)
  1090:            
  1091:                 return
  1092:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAP2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1111: inl(1222): Inlined: REMAPX
  1112: inl(1222): Inlined: REMAPX
  1113: inl(1222): Inlined: REMAPX
  1114: inl(1235): Inline halted: code size exceeds limit.: REMAP2
  1114: inl(1222): Inlined: REMAPX
  1115: inl(1222): Inlined: REMAPY
  1116: inl(1222): Inlined: REMAPY
  1117: inl(1222): Inlined: REMAPZ


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAP2
FORMAT LIST

 LINE   LOOP      STATEMENT

  1095:                 subroutine remap2(iela, ielnew)
  1096:           c-----------------------------------------------------------------
  1097:           c     Map the solution from the children to the parent.
  1098:           c     iela array records the eight elements to be merged.
  1099:           c     ielnew is the element index of the merged element.
  1100:           c-----------------------------------------------------------------
  1101:           
  1102:                 use ua_data
  1103:                 implicit none
  1104:           
  1105:                 integer iela(8), ielnew
  1106:           
  1107:                 double precision temp1(lx1,lx1,lx1),
  1108:                &       temp2(lx1,lx1,lx1),temp3(lx1,lx1,lx1),temp4(lx1,lx1,lx1),
  1109:                &       temp5(lx1,lx1,lx1),temp6(lx1,lx1,lx1)
  1110:           
  1111:         I       call remapx(ta1(1,1,1,iela(1)),ta1(1,1,1,iela(2)),temp1)
  1112:         I       call remapx(ta1(1,1,1,iela(3)),ta1(1,1,1,iela(4)),temp2)
  1113:         I       call remapx(ta1(1,1,1,iela(5)),ta1(1,1,1,iela(6)),temp3)
  1114:         I       call remapx(ta1(1,1,1,iela(7)),ta1(1,1,1,iela(8)),temp4)
  1115:         I       call remapy(temp1,temp2,temp5)
  1116:         I       call remapy(temp3,temp4,temp6)
  1117:         I       call remapz(temp5,temp6,ta1(1,1,1,ielnew))
  1118:           
  1119:                 return
  1120:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPZ
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1135: opt(1593): Loop nest collapsed into one loop.: IY
  1135: vec( 101): Vectorized loop.
  1140: opt(1418): Constant-length loop is expanded.
  1141: vec( 128): Fused multiply-add operation applied.
  1147: opt(1418): Constant-length loop is expanded.
  1148: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPZ
FORMAT LIST

 LINE   LOOP      STATEMENT

  1123:                 subroutine remapz(x1,x2,y)
  1124:           c-----------------------------------------------------------------
  1125:           c     z direction mapping after the merge.
  1126:           c     Map solution from x1 & x2 to y.
  1127:           c-----------------------------------------------------------------
  1128:           
  1129:                 use ua_data
  1130:                 implicit none
  1131:           
  1132:                 double precision x1(lx1,lx1,lx1),x2(lx1,lx1,lx1),y(lx1,lx1,lx1)
  1133:                 integer ix, iy, ip
  1134:           
  1135: W------>        do iy=1,lx1
  1136: |*----->          do ix=1,lx1
  1137: ||                  y(ix,iy,1)=x1(ix,iy,1)
  1138: ||        
  1139: ||                  y(ix,iy,2)=0.d0
  1140: ||*---->            do ip=1,lx1
  1141: |||     F             y(ix,iy,2)=y(ix,iy,2)+map2(ip)*x1(ix,iy,ip)
  1142: ||*----             end do
  1143: ||        
  1144: ||                  y(ix,iy,3)=x1(ix,iy,lx1)
  1145: ||        
  1146: ||                  y(ix,iy,4)=0.d0
  1147: ||*---->            do ip=1,lx1
  1148: |||     F             y(ix,iy,4)=y(ix,iy,4)+map4(ip)*x2(ix,iy,ip)
  1149: ||*----             end do
  1150: ||        
  1151: ||                  y(ix,iy,lx1)=x2(ix,iy,lx1)
  1152: |*-----           end do
  1153: W------         end do
  1154:           
  1155:                 return
  1156:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPY
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1172: vec( 101): Vectorized loop.
  1176: opt(1418): Constant-length loop is expanded.
  1177: vec( 128): Fused multiply-add operation applied.
  1183: opt(1418): Constant-length loop is expanded.
  1184: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPY
FORMAT LIST

 LINE   LOOP      STATEMENT

  1159:                 subroutine remapy(x1,x2,y)
  1160:           c-----------------------------------------------------------------
  1161:           c     y direction mapping after the merge.
  1162:           c     Map solution from x1 & x2 to y.
  1163:           c-----------------------------------------------------------------
  1164:           
  1165:                 use ua_data
  1166:                 implicit none
  1167:           
  1168:                 double precision x1(lx1,lx1,lx1),x2(lx1,lx1,lx1),y(lx1,lx1,lx1)
  1169:                 integer ix, iz, ip
  1170:           
  1171: +------>        do iz=1,lx1
  1172: |V----->          do ix=1,lx1
  1173: ||                  y(ix,1,iz)=x1(ix,1,iz)
  1174: ||        
  1175: ||                  y(ix,2,iz)=0.d0
  1176: ||*---->            do ip=1,lx1
  1177: |||     F             y(ix,2,iz)=y(ix,2,iz)+map2(ip)*x1(ix,ip,iz)
  1178: ||*----             end do
  1179: ||        
  1180: ||                  y(ix,3,iz)=x1(ix,lx1,iz)
  1181: ||        
  1182: ||                  y(ix,4,iz)=0.d0
  1183: ||*---->            do ip=1,lx1
  1184: |||     F             y(ix,4,iz)=y(ix,4,iz)+map4(ip)*x2(ix,ip,iz)
  1185: ||*----             end do
  1186: ||        
  1187: ||                  y(ix,lx1,iz)=x2(ix,lx1,iz)
  1188: |V-----           end do
  1189: +------         end do
  1190:           
  1191:                 return
  1192:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPX
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  1207: opt(1593): Loop nest collapsed into one loop.: IZ
  1207: opt(1593): Loop nest collapsed into one loop.: IZ
  1207: vec( 101): Vectorized loop.
  1212: opt(1418): Constant-length loop is expanded.
  1213: vec( 128): Fused multiply-add operation applied.
  1219: opt(1418): Constant-length loop is expanded.
  1220: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:29 2020
FILE NAME: adapt.f

PROCEDURE NAME: REMAPX
FORMAT LIST

 LINE   LOOP      STATEMENT

  1195:                 subroutine remapx(x1,x2,y)
  1196:           c-----------------------------------------------------------------
  1197:           c     x direction mapping after the merge.
  1198:           c     Map solution from x1 & x2 to y.
  1199:           c-----------------------------------------------------------------
  1200:           
  1201:                 use ua_data
  1202:                 implicit none
  1203:           
  1204:                 double precision x1(lx1,lx1,lx1),x2(lx1,lx1,lx1),y(lx1,lx1,lx1)
  1205:                 integer iy, iz, ip
  1206:           
  1207: W------>        do iz=1,lx1
  1208: |*----->          do iy=1,lx1
  1209: ||                  y(1,iy,iz)=x1(1,iy,iz)
  1210: ||        
  1211: ||                  y(2,iy,iz)=0.d0
  1212: ||*---->            do ip=1,lx1
  1213: |||     F             y(2,iy,iz)=y(2,iy,iz)+map2(ip)*x1(ip,iy,iz)
  1214: ||*----             end do
  1215: ||        
  1216: ||                  y(3,iy,iz)=x1(lx1,iy,iz)
  1217: ||        
  1218: ||                  y(4,iy,iz)=0.d0
  1219: ||*---->            do ip=1,lx1
  1220: |||     F             y(4,iy,iz)=y(4,iy,iz)+map4(ip)*x2(ip,iy,iz)
  1221: ||*----             end do
  1222: ||        
  1223: ||                  y(lx1,iy,iz)=x2(lx1,iy,iz)
  1224: |*-----           end do
  1225: W------         end do
  1226:           
  1227:                 return
  1228:                 end


