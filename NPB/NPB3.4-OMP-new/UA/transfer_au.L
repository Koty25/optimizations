NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: INIT_LOCKS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    14: par(1801): Parallel routine generated.: INIT_LOCKS$1
    15: par(1803): Parallelized by "do".
    15: vec( 103): Unvectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: INIT_LOCKS
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine init_locks
     3:           c------------------------------------------------------------------
     4:           c     This version uses ATOMIC for atomic updates,
     5:           c     but locks are still used in get_emo (mason.f).
     6:           c------------------------------------------------------------------
     7:           
     8:                 use ua_data
     9:                 implicit none
    10:           
    11:                 integer i
    12:           
    13:           c.....initialize locks in parallel
    14:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    15: P------>  c$    do i=1,8*lelt
    16: |         c$      call omp_init_lock(tlock(i))
    17: P------   c$    end do
    18:           
    19:                 return
    20:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSF
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    40: par(1801): Parallel routine generated.: TRANSF$1
    42: par(1803): Parallelized by "do".
    43: vec( 110): Vectorization obstructive procedure reference.: R_INIT
    76: opt(1025): Reference to this procedure inhibits optimization.: R_INIT
    79: vec( 102): Partially vectorized loop.
    80: opt(1036): Potential feedback - use directive or compiler option if OK.: TMP
    80: vec( 102): Partially vectorized loop.
    86: opt(1033): Potential multiple store conflict -- use directive if OK.: TMP
    86: opt(1036): Potential feedback - use directive or compiler option if OK.: TMP
    86: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: TMP$57
    90: opt(1418): Constant-length loop is expanded.
    92: opt(1418): Constant-length loop is expanded.
    94: opt(1036): Potential feedback - use directive or compiler option if OK.: TMP
    94: opt(1019): Feedback of scalar value from one loop pass to another.: TMP
   106: vec( 108): Unvectorizable loop structure.
   114: opt(1418): Constant-length loop is expanded.
   116: opt(1418): Constant-length loop is expanded.
   117: opt(1036): Potential feedback - use directive or compiler option if OK.
   117: opt(1019): Feedback of scalar value from one loop pass to another.
   123: vec( 103): Unvectorized loop.
   123: vec( 107): Iteration count is too small.
   129: opt(1036): Potential feedback - use directive or compiler option if OK.
   133: opt(1418): Constant-length loop is expanded.
   135: opt(1418): Constant-length loop is expanded.
   143: opt(1418): Constant-length loop is expanded.
   145: opt(1418): Constant-length loop is expanded.
   146: opt(1019): Feedback of scalar value from one loop pass to another.
   156: opt(1418): Constant-length loop is expanded.
   157: opt(1418): Constant-length loop is expanded.
   169: vec( 103): Unvectorized loop.
   169: vec( 107): Iteration count is too small.
   171: opt(1418): Constant-length loop is expanded.
   172: opt(1418): Constant-length loop is expanded.
   174: opt(1036): Potential feedback - use directive or compiler option if OK.
   181: opt(1418): Constant-length loop is expanded.
   190: vec( 103): Unvectorized loop.
   190: vec( 107): Iteration count is too small.
   192: opt(1418): Constant-length loop is expanded.
   193: opt(1418): Constant-length loop is expanded.
   195: opt(1036): Potential feedback - use directive or compiler option if OK.
   202: opt(1418): Constant-length loop is expanded.
   211: vec( 103): Unvectorized loop.
   211: vec( 107): Iteration count is too small.
   213: opt(1418): Constant-length loop is expanded.
   214: opt(1418): Constant-length loop is expanded.
   216: opt(1036): Potential feedback - use directive or compiler option if OK.
   223: opt(1418): Constant-length loop is expanded.
   232: vec( 103): Unvectorized loop.
   232: vec( 107): Iteration count is too small.
   234: opt(1418): Constant-length loop is expanded.
   235: opt(1418): Constant-length loop is expanded.
   237: opt(1036): Potential feedback - use directive or compiler option if OK.
   243: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSF
FORMAT LIST

 LINE   LOOP      STATEMENT

    24:                 subroutine transf(tmor,tx)
    25:           c------------------------------------------------------------------
    26:           c     Map values from mortar(tmor) to element(tx)
    27:           c------------------------------------------------------------------
    28:           
    29:                 use ua_data
    30:                 implicit none
    31:           
    32:                 double precision tmor(*),tx(*), tmp(lx1,lx1,2)
    33:                 integer ig1,ig2,ig3,ig4,ie,iface,il1,il2,il3,il4,
    34:                &        nnje,ije1,ije2,col,i,j,ig,il
    35:           
    36:           
    37:           c.....zero out tx on element boundaries
    38:                 call col2(tx,tmult,ntot)
    39:           
    40:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(il,j,ig,i,col,ije2,ije1,ig4,
    41:           c$OMP& ig3,ig2,ig1,nnje,il4,il3,il2,il1,iface,ie,tmp)
    42: P------>        do ie=1,nelt
    43: |+----->          do iface=1,nsides
    44: ||        
    45: ||        c.........get the collocation point index of the four local corners on the
    46: ||        c         face iface of element ie
    47: ||                  il1=idel(1,1,iface,ie)
    48: ||                  il2=idel(lx1,1,iface,ie)
    49: ||                  il3=idel(1,lx1,iface,ie)
    50: ||                  il4=idel(lx1,lx1,iface,ie)
    51: ||        
    52: ||        c.........get the mortar indices of the four local corners
    53: ||                  ig1= idmo(1,  1  ,1,1,iface,ie)
    54: ||                  ig2= idmo(lx1,1  ,1,2,iface,ie)
    55: ||                  ig3= idmo(1,  lx1,2,1,iface,ie)
    56: ||                  ig4= idmo(lx1,lx1,2,2,iface,ie)
    57: ||         
    58: ||        c.........copy the value from tmor to tx for these four local corners
    59: ||                  tx(il1) = tmor(ig1)
    60: ||                  tx(il2) = tmor(ig2)
    61: ||                  tx(il3) = tmor(ig3)
    62: ||                  tx(il4) = tmor(ig4)
    63: ||         
    64: ||        c.........nnje=1 for conforming faces, nnje=2 for nonconforming faces
    65: ||                  if(cbc(iface,ie).eq.3) then
    66: ||                    nnje=2
    67: ||                  else
    68: ||                    nnje=1
    69: ||                  end if
    70: ||        
    71: ||        c.........for nonconforming faces
    72: ||                  if(nnje.eq.2) then
    73: ||        
    74: ||        c...........nonconforming faces have four pieces of mortar, first map them to
    75: ||        c           two intermediate mortars, stored in tmp
    76: ||                    call r_init(tmp,lx1*lx1*2,0.d0)
    77: ||         
    78: ||+---->              do ije1=1,nnje
    79: |||S--->                do ije2=1,nnje
    80: ||||S-->                  do col=1,lx1
    81: |||||     
    82: |||||     c.................in each row col, when coloumn i=1 or lx1, the value
    83: |||||     c                 in tmor is copied to tmp
    84: |||||                       i = v_end(ije2)
    85: |||||                       ig=idmo(i,col,ije1,ije2,iface,ie)
    86: |||||   G                   tmp(i,col,ije1)=tmor(ig)
    87: |||||     
    88: |||||     c.................in each row col, value in the interior three collocation
    89: |||||     c                 points is computed by apply mapping matrix qbnew to tmor
    90: |||||*->                    do i=2,lx1-1
    91: ||||||                        il= idel(i,col,iface,ie)
    92: ||||||*>                      do j=1,lx1
    93: |||||||                         ig=idmo(j,col,ije1,ije2,iface,ie)
    94: ||||||| G                       tmp(i,col,ije1) = tmp(i,col,ije1) +
    95: |||||||        &                qbnew(i-1,j,ije2)*tmor(ig)
    96: ||||||*                       end do
    97: |||||*-                     end do
    98: |||||     
    99: ||||S--                   end do
   100: |||S---                 end do
   101: ||+----               end do
   102: ||         
   103: ||        c...........mapping from two pieces of intermediate mortar tmp to element
   104: ||        c           face tx
   105: ||        
   106: ||+---->              do ije1=1, nnje
   107: |||       
   108: |||       c.............the first column, col=1, is an edge of face iface.
   109: |||       c             the value on the three interior collocation points, tx, is
   110: |||       c             computed by applying mapping matrices qbnew to tmp.
   111: |||       c             the mapping result is divided by 2, because there will be
   112: |||       c             duplicated contribution from another face sharing this edge.
   113: |||                     col=1
   114: |||*--->                do i=2,lx1-1
   115: ||||                      il= idel(col,i,iface,ie)
   116: ||||*-->                  do j=1,lx1
   117: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*
   118: |||||          &                       tmp(col,j,ije1)*0.5d0
   119: ||||*--                   end do
   120: |||*---                 end do
   121: |||       
   122: |||       c.............for column 2 ~ lx-1
   123: |||+--->                do col=2,lx1-1
   124: ||||      
   125: ||||      c...............when i=1 or lx1, the collocation points are also on an edge of
   126: ||||      c               the face, so the mapping result also needs to be divided by 2
   127: ||||                      i = v_end(ije1)
   128: ||||                      il= idel(col,i,iface,ie)
   129: ||||                      tx(il)=tx(il)+tmp(col,i,ije1)*0.5d0
   130: ||||      
   131: ||||      c...............compute the value at interior collocation points in
   132: ||||      c               columns 2 ~ lx1
   133: ||||*-->                  do i=2,lx1-1
   134: |||||                       il= idel(col,i,iface,ie)
   135: |||||*->                    do j=1,lx1
   136: ||||||                        tx(il) = tx(il) + qbnew(i-1,j,ije1)* tmp(col,j,ije1)
   137: |||||*-                     end do
   138: ||||*--                   end do
   139: |||+---                 end do
   140: |||       
   141: |||       c.............same as col=1
   142: |||                     col=lx1
   143: |||*--->                do  i=2,lx1-1
   144: ||||                      il= idel(col,i,iface,ie)
   145: ||||*-->                  do j=1,lx1
   146: |||||                       tx(il) = tx(il) + qbnew(i-1,j,ije1)*
   147: |||||          &                     tmp(col,j,ije1)*0.5d0
   148: ||||*--                   end do
   149: |||*---                 end do
   150: ||+----               end do
   151: ||        
   152: ||        c.........for conforming faces
   153: ||                  else
   154: ||        
   155: ||        c.........face interior
   156: ||*---->              do col=2,lx1-1
   157: |||*--->                do i=2,lx1-1
   158: ||||                      il= idel(i,col,iface,ie)
   159: ||||                      ig= idmo(i,col,1,1,iface,ie)
   160: ||||                      tx(il)=tmor(ig)
   161: |||*---                 end do
   162: ||*----               end do
   163: ||        
   164: ||         
   165: ||        c...........edges of conforming faces
   166: ||        
   167: ||        c...........if local edge 1 is a nonconforming edge
   168: ||                    if(idmo(lx1,1,1,1,iface,ie).ne.0)then
   169: ||+---->                do i=2,lx1-1
   170: |||                       il= idel(i,1,iface,ie)
   171: |||*--->                  do ije1=1,2
   172: ||||*-->                    do j=1,lx1
   173: |||||                         ig=idmo(j,1,1,ije1,iface,ie)
   174: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   175: ||||*--                     end do
   176: |||*---                   end do
   177: ||+----                 end do
   178: ||        
   179: ||        c...........if local edge 1 is a conforming edge
   180: ||                    else
   181: ||*---->                do i=2,lx1-1
   182: |||                       il= idel(i,1,iface,ie)
   183: |||                       ig= idmo(i,1,1,1,iface,ie)
   184: |||                       tx(il)=tmor(ig)
   185: ||*----                 end do
   186: ||                    end if
   187: ||        
   188: ||        c...........if local edge 2 is a nonconforming edge
   189: ||                    if(idmo(lx1,2,1,2,iface,ie).ne.0)then
   190: ||+---->                do i=2,lx1-1
   191: |||                       il= idel(lx1,i,iface,ie)
   192: |||*--->                  do ije1=1,2
   193: ||||*-->                    do j=1,lx1
   194: |||||                         ig=idmo(lx1,j,ije1,2,iface,ie)
   195: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   196: ||||*--                     end do
   197: |||*---                   end do
   198: ||+----                 end do
   199: ||        
   200: ||        c...........if local edge 2 is a conforming edge
   201: ||                    else
   202: ||*---->                do i=2,lx1-1
   203: |||                       il= idel(lx1,i,iface,ie)
   204: |||                       ig= idmo(lx1,i,1,1,iface,ie)
   205: |||                       tx(il)=tmor(ig)
   206: ||*----                 end do
   207: ||                    end if
   208: ||        
   209: ||        c...........if local edge 3 is a nonconforming edge
   210: ||                    if(idmo(2,lx1,2,1,iface,ie).ne.0)then
   211: ||+---->                do  i=2,lx1-1
   212: |||                       il= idel(i,lx1,iface,ie)
   213: |||*--->                  do ije1=1,2
   214: ||||*-->                    do j=1,lx1
   215: |||||                         ig=idmo(j,lx1,2,ije1,iface,ie)
   216: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   217: ||||*--                     end do
   218: |||*---                   end do
   219: ||+----                 end do
   220: ||        
   221: ||        c...........if local edge 3 is a conforming edge
   222: ||                    else
   223: ||*---->                do i=2,lx1-1
   224: |||                       il= idel(i,lx1,iface,ie)
   225: |||                       ig= idmo(i,lx1,1,1,iface,ie)
   226: |||                       tx(il)=tmor(ig)
   227: ||*----                 end do
   228: ||                    end if
   229: ||        
   230: ||        c...........if local edge 4 is a nonconforming edge
   231: ||                    if(idmo(1,lx1,1,1,iface,ie).ne.0)then
   232: ||+---->                do i=2,lx1-1
   233: |||                       il= idel(1,i,iface,ie)
   234: |||*--->                  do ije1=1,2
   235: ||||*-->                    do j=1,lx1
   236: |||||                         ig=idmo(1,j,ije1,1,iface,ie)
   237: |||||                         tx(il) = tx(il) + qbnew(i-1,j,ije1)*tmor(ig)*0.5d0
   238: ||||*--                     end do
   239: |||*---                   end do
   240: ||+----                 end do
   241: ||        c...........if local edge 4 is a conforming edge
   242: ||                    else
   243: ||*---->                do i=2,lx1-1
   244: |||                       il= idel(1,i,iface,ie)
   245: |||                       ig= idmo(1,i,1,1,iface,ie)
   246: |||                       tx(il)=tmor(ig)
   247: ||*----                 end do
   248: ||                    end if
   249: ||                  end if
   250: ||         
   251: |+-----           end do
   252: P------         end do
   253:           c$OMP END PARALLEL DO
   254:           
   255:                 return
   256:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   278: par(1801): Parallel routine generated.: TRANSFB$1
   283: par(1803): Parallelized by "do".
   283: vec( 101): Vectorized loop.
   289: par(1803): Parallelized by "do".
   313: par(3007): Atomic-construct generated.
   315: par(3007): Atomic-construct generated.
   317: par(3007): Atomic-construct generated.
   319: par(3007): Atomic-construct generated.
   329: vec( 102): Partially vectorized loop.
   329: vec( 113): Overhead of loop division is too large.
   331: vec( 102): Partially vectorized loop.
   335: opt(1033): Potential multiple store conflict -- use directive if OK.: TEMP
   335: opt(1036): Potential feedback - use directive or compiler option if OK.: TEMP
   335: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: TEMP$105
   341: opt(1418): Constant-length loop is expanded.
   343: vec( 128): Fused multiply-add operation applied.
   346: opt(1033): Potential multiple store conflict -- use directive if OK.: TOP
   350: opt(1418): Constant-length loop is expanded.
   352: opt(1418): Constant-length loop is expanded.
   354: vec( 128): Fused multiply-add operation applied.
   356: opt(1036): Potential feedback - use directive or compiler option if OK.: TEMP
   375: par(3007): Atomic-construct generated.
   382: vec( 103): Unvectorized loop.
   382: vec( 109): Vectorization obstructive statement.
   384: opt(1418): Constant-length loop is expanded.
   388: par(3007): Atomic-construct generated.
   400: par(3007): Atomic-construct generated.
   402: vec( 103): Unvectorized loop.
   402: vec( 109): Vectorization obstructive statement.
   405: opt(1418): Constant-length loop is expanded.
   410: par(3007): Atomic-construct generated.
   420: vec( 103): Unvectorized loop.
   420: vec( 109): Vectorization obstructive statement.
   421: opt(1418): Constant-length loop is expanded.
   424: par(3007): Atomic-construct generated.
   433: opt(1418): Constant-length loop is expanded.
   434: vec( 103): Unvectorized loop.
   434: vec( 109): Vectorization obstructive statement.
   436: opt(1418): Constant-length loop is expanded.
   441: par(3007): Atomic-construct generated.
   448: opt(1418): Constant-length loop is expanded.
   451: par(3007): Atomic-construct generated.
   458: opt(1418): Constant-length loop is expanded.
   459: vec( 103): Unvectorized loop.
   459: vec( 109): Vectorization obstructive statement.
   461: opt(1418): Constant-length loop is expanded.
   466: par(3007): Atomic-construct generated.
   473: opt(1418): Constant-length loop is expanded.
   476: par(3007): Atomic-construct generated.
   483: opt(1418): Constant-length loop is expanded.
   484: vec( 103): Unvectorized loop.
   484: vec( 109): Vectorization obstructive statement.
   486: opt(1418): Constant-length loop is expanded.
   491: par(3007): Atomic-construct generated.
   498: opt(1418): Constant-length loop is expanded.
   501: par(3007): Atomic-construct generated.
   508: opt(1418): Constant-length loop is expanded.
   509: vec( 103): Unvectorized loop.
   509: vec( 109): Vectorization obstructive statement.
   511: opt(1418): Constant-length loop is expanded.
   516: par(3007): Atomic-construct generated.
   523: opt(1418): Constant-length loop is expanded.
   526: par(3007): Atomic-construct generated.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB
FORMAT LIST

 LINE   LOOP      STATEMENT

   260:                 subroutine transfb(tmor,tx)
   261:           c------------------------------------------------------------------
   262:           c     Map from element(tx) to mortar(tmor).
   263:           c     tmor sums contributions from all elements.
   264:           c------------------------------------------------------------------
   265:           
   266:                 use ua_data
   267:                 implicit none
   268:           
   269:                 double precision third
   270:                 parameter (third=1.d0/3.d0)
   271:                 integer shift
   272:           
   273:                 double precision tmp,tmp1,tx(*),tmor(*),temp(lx1,lx1,2),
   274:                &                 top(lx1,2)
   275:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,nnje,
   276:                &        ije1,ije2,col,i,j,ije,ig,il
   277:           
   278:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(il,j,ig,i,col,ije2,ije1,ig4,
   279:           c$OMP& ig3,ig2,ig1,nnje,il4,il3,il2,il1,iface,ie,ije,
   280:           c$OMP& tmp,shift,temp,top,tmp1)
   281:           
   282:           c$OMP DO
   283: Y------>        do ie=1,nmor
   284: |                 tmor(ie)=0.d0
   285: Y------         end do
   286:           c$OMP END DO
   287:           
   288:           c$OMP DO
   289: P------>        do ie=1,nelt
   290: |+----->          do iface=1,nsides
   291: ||        c.........nnje=1 for conforming faces, nnje=2 for nonconforming faces
   292: ||                  if(cbc(iface,ie).eq.3) then
   293: ||                    nnje=2
   294: ||                  else
   295: ||                    nnje=1
   296: ||                  end if
   297: ||        
   298: ||        c.........get collocation point index of four local corners on the face
   299: ||                  il1 = idel(1,  1,  iface,ie)
   300: ||                  il2 = idel(lx1,1,  iface,ie)
   301: ||                  il3 = idel(1,  lx1,iface,ie)
   302: ||                  il4 = idel(lx1,lx1,iface,ie)
   303: ||        
   304: ||        c.........get the mortar indices of the four local corners
   305: ||                  ig1 = idmo(1,  1,  1,1,iface,ie)
   306: ||                  ig2 = idmo(lx1,1,  1,2,iface,ie)
   307: ||                  ig3 = idmo(1,  lx1,2,1,iface,ie )
   308: ||                  ig4 = idmo(lx1,lx1,2,2,iface,ie)
   309: ||        
   310: ||        c.........sum the values from tx to tmor for these four local corners
   311: ||        c         only 1/3 of the value is summed, since there will be two duplicated
   312: ||        c         contributions from the other two faces sharing this vertex
   313: ||        c$OMP ATOMIC
   314: ||                  tmor(ig1) = tmor(ig1)+tx(il1)*third
   315: ||        c$OMP ATOMIC
   316: ||                  tmor(ig2) = tmor(ig2)+tx(il2)*third
   317: ||        c$OMP ATOMIC
   318: ||                  tmor(ig3) = tmor(ig3)+tx(il3)*third
   319: ||        c$OMP ATOMIC
   320: ||                  tmor(ig4) = tmor(ig4)+tx(il4)*third
   321: ||        
   322: ||        c.........for nonconforming faces
   323: ||                  if(nnje.eq.2) then
   324: ||                    call r_init(temp,lx1*lx1*2,0.d0)
   325: ||        
   326: ||        c...........nonconforming faces have four pieces of mortar, first map tx to
   327: ||        c           two intermediate mortars stored in temp
   328: ||        
   329: ||S---->              do ije2 = 1, nnje
   330: |||                     shift = ije2-1
   331: |||S--->                do col=1,lx1
   332: ||||      c...............For mortar points on face edge (top and bottom), copy the
   333: ||||      c               value from tx to temp
   334: ||||                      il=idel(col,v_end(ije2),iface,ie)
   335: ||||    G                 temp(col,v_end(ije2),ije2)=tx(il)
   336: ||||      
   337: ||||      c...............For mortar points on face edge (top and bottom), calculate
   338: ||||      c               the interior points' contribution to them, i.e. top()
   339: ||||                      j = v_end(ije2)
   340: ||||                      tmp=0.d0
   341: ||||*-->                  do i=2,lx1-1
   342: |||||                       il=idel(col,i,iface,ie)
   343: |||||   G                   tmp = tmp + qbnew(i-1,j,ije2)*tx(il)
   344: ||||*--                   end do
   345: ||||      
   346: ||||                      top(col,ije2)=tmp
   347: ||||      
   348: ||||      c...............Use mapping matrices qbnew to map the value from tx to temp
   349: ||||      c               for mortar points not on the top bottom face edge.
   350: ||||*-->                  do j=2-shift,lx1-shift
   351: |||||                       tmp=0.d0
   352: |||||*->                    do i=2,lx1-1
   353: ||||||                        il=idel(col,i,iface,ie)
   354: ||||||  F                     tmp = tmp + qbnew(i-1,j,ije2)*tx(il)
   355: |||||*-                     end do
   356: |||||                       temp(col,j,ije2) = tmp + temp(col,j,ije2)
   357: ||||*--                   end do
   358: |||S---                 end do
   359: ||S----               end do
   360: ||        
   361: ||        c...........mapping from temp to tmor
   362: ||        
   363: ||+---->              do ije1=1, nnje
   364: |||                     shift = ije1-1
   365: |||+--->                do ije2=1,nnje
   366: ||||      
   367: ||||      c...............for each column of collocation points on a piece of mortar
   368: ||||+-->                  do col=2-shift,lx1-shift
   369: |||||     
   370: |||||     c.................For the end point, which is on an edge (local edge 2,4),
   371: |||||     c                 the contribution is halved since there will be duplicated
   372: |||||     c                 contribution from another face sharing this edge.
   373: |||||     
   374: |||||                       ig=idmo(v_end(ije2),col,ije1,ije2,iface,ie)
   375: |||||     c$OMP ATOMIC
   376: |||||                       tmor(ig)=tmor(ig)+temp(v_end(ije2),col,ije1)*0.5d0
   377: |||||     
   378: |||||     c.................In each row of collocation points on a piece of mortar,
   379: |||||     c                 sum the contributions from interior collocation points
   380: |||||     c                 (i=2,lx1-1)
   381: |||||     
   382: |||||+->                    do  j=1,lx1
   383: ||||||                        tmp=0.d0
   384: ||||||*>                      do i=2,lx1-1
   385: |||||||                         tmp = tmp + qbnew(i-1,j,ije2) * temp(i,col,ije1)
   386: ||||||*                       end do
   387: ||||||                        ig=idmo(j,col,ije1,ije2,iface,ie)
   388: ||||||    c$OMP ATOMIC
   389: ||||||                        tmor(ig)=tmor(ig)+tmp
   390: |||||+-                     end do
   391: ||||+--                   end do
   392: ||||      
   393: ||||      c...............For tmor on local edge 1 and 3, tmp is the contribution from
   394: ||||      c               an edge, so it is halved because of duplicated contribution
   395: ||||      c               from another face sharing this edge. tmp1 is contribution
   396: ||||      c               from face interior.
   397: ||||      
   398: ||||                      col = v_end(ije1)
   399: ||||                      ig=idmo(v_end(ije2),col,ije1,ije2,iface,ie)
   400: ||||      c$OMP ATOMIC
   401: ||||                      tmor(ig)=tmor(ig)+top(v_end(ije2),ije1)*0.5d0
   402: ||||+-->                  do  j=1,lx1
   403: |||||                       tmp=0.d0
   404: |||||                       tmp1=0.d0
   405: |||||*->                    do i=2,lx1-1
   406: ||||||                        tmp  = tmp  + qbnew(i-1,j,ije2) * temp(i,col,ije1)
   407: ||||||                        tmp1 = tmp1 + qbnew(i-1,j,ije2) * top(i,ije1)
   408: |||||*-                     end do
   409: |||||                       ig=idmo(j,col,ije1,ije2,iface,ie)
   410: |||||     c$OMP ATOMIC
   411: |||||                       tmor(ig)=tmor(ig)+tmp*0.5d0+tmp1
   412: ||||+--                   end do
   413: |||+---                 end do
   414: ||+----               end do
   415: ||        
   416: ||        c.........for conforming faces
   417: ||                  else
   418: ||        
   419: ||        c.........face interior
   420: ||+---->              do col=2,lx1-1
   421: |||*--->                do j=2,lx1-1
   422: ||||                      il=idel(j,col,iface,ie)
   423: ||||                      ig=idmo(j,col,1,1,iface,ie)
   424: ||||      c$OMP ATOMIC
   425: ||||                      tmor(ig)=tmor(ig)+tx(il)
   426: |||*---                 end do
   427: ||+----               end do
   428: ||        
   429: ||        c...........edges of conforming faces
   430: ||        
   431: ||        c...........if local edge 1 is a nonconforming edge
   432: ||                    if(idmo(lx1,1,1,1,iface,ie).ne.0)then
   433: ||*---->                do ije=1,2
   434: |||+--->                  do j=1,lx1
   435: ||||                        tmp=0.d0
   436: ||||*-->                    do i=2,lx1-1
   437: |||||                         il=idel(i,1,iface,ie)
   438: |||||                         tmp= tmp + qbnew(i-1,j,ije)*tx(il)
   439: ||||*--                     end do
   440: ||||                        ig=idmo(j,1,1,ije,iface,ie)
   441: ||||      c$OMP ATOMIC
   442: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   443: |||+---                   end do
   444: ||*----                 end do
   445: ||        
   446: ||        c...........if local edge 1 is a conforming edge
   447: ||                    else
   448: ||*---->                do j=2,lx1-1
   449: |||                       il=idel(j,1,iface,ie)
   450: |||                       ig=idmo(j,1,1,1,iface,ie)
   451: |||       c$OMP ATOMIC
   452: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   453: ||*----                 end do
   454: ||                    end if
   455: ||        
   456: ||        c...........if local edge 2 is a nonconforming edge
   457: ||                    if(idmo(lx1,2,1,2,iface,ie).ne.0)then
   458: ||*---->                do ije=1,2
   459: |||+--->                  do j=1,lx1
   460: ||||                        tmp=0.d0
   461: ||||*-->                    do i=2,lx1-1
   462: |||||                         il=idel(lx1,i,iface,ie)
   463: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   464: ||||*--                     end do
   465: ||||                        ig=idmo(lx1,j,ije,2,iface,ie)
   466: ||||      c$OMP ATOMIC
   467: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   468: |||+---                   end do
   469: ||*----                 end do
   470: ||        
   471: ||        c...........if local edge 2 is a conforming edge
   472: ||                    else
   473: ||*---->                do j=2,lx1-1
   474: |||                       il=idel(lx1,j,iface,ie)
   475: |||                       ig=idmo(lx1,j,1,1,iface,ie)
   476: |||       c$OMP ATOMIC
   477: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   478: ||*----                 end do
   479: ||                    end if
   480: ||        
   481: ||        c...........if local edge 3 is a nonconforming edge
   482: ||                    if(idmo(2,lx1,2,1,iface,ie).ne.0)then
   483: ||*---->                do ije=1,2
   484: |||+--->                  do j=1,lx1
   485: ||||                        tmp=0.d0
   486: ||||*-->                    do i=2,lx1-1
   487: |||||                         il=idel(i,lx1,iface,ie)
   488: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   489: ||||*--                     end do
   490: ||||                        ig=idmo(j,lx1,2,ije,iface,ie)
   491: ||||      c$OMP ATOMIC
   492: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   493: |||+---                   end do
   494: ||*----                 end do
   495: ||        
   496: ||        c...........if local edge 3 is a conforming edge
   497: ||                    else
   498: ||*---->                do j=2,lx1-1
   499: |||                       il=idel(j,lx1,iface,ie)
   500: |||                       ig=idmo(j,lx1,1,1,iface,ie)
   501: |||       c$OMP ATOMIC
   502: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   503: ||*----                 end do
   504: ||                    end if
   505: ||        
   506: ||        c...........if local edge 4 is a nonconforming edge
   507: ||                    if(idmo(1,lx1,1,1,iface,ie).ne.0)then
   508: ||*---->                do ije=1,2
   509: |||+--->                  do j=1,lx1
   510: ||||                        tmp=0.d0
   511: ||||*-->                    do i=2,lx1-1
   512: |||||                         il=idel(1,i,iface,ie)
   513: |||||                         tmp = tmp + qbnew(i-1,j,ije)*tx(il)
   514: ||||*--                     end do
   515: ||||                        ig=idmo(1,j,ije,1,iface,ie)
   516: ||||      c$OMP ATOMIC
   517: ||||                        tmor(ig)=tmor(ig)+tmp*0.5d0
   518: |||+---                   end do
   519: ||*----                 end do
   520: ||        
   521: ||        c...........if local edge 4 is a conforming edge
   522: ||                    else
   523: ||*---->                do j=2,lx1-1
   524: |||                       il=idel(1,j,iface,ie)
   525: |||                       ig=idmo(1,j,1,1,iface,ie)
   526: |||       c$OMP ATOMIC
   527: |||                       tmor(ig)=tmor(ig)+tx(il)*0.5d0
   528: ||*----                 end do
   529: ||                    end if
   530: ||                  end if
   531: |+-----           end do
   532: P------         end do
   533:           c$OMP END DO NOWAIT
   534:           c$OMP END PARALLEL
   535:           
   536:                 return
   537:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_COR_E
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   560: opt(1418): Constant-length loop is expanded.
   565: opt(1418): Constant-length loop is expanded.
   571: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_COR_E
FORMAT LIST

 LINE   LOOP      STATEMENT

   541:                 subroutine transfb_cor_e(n,tmor,tx)
   542:           c--------------------------------------------------------------
   543:           c     This subroutine performs the edge to mortar mapping and
   544:           c     calculates the mapping result on the mortar point at a vertex
   545:           c     under situation 1,2, or 3.
   546:           c     n refers to the configuration of three edges sharing a vertex,
   547:           c     n = 1: only one edge is nonconforming
   548:           c     n = 2: two edges are nonconforming
   549:           c     n = 3: three edges are nonconforming
   550:           c-------------------------------------------------------------------
   551:           
   552:                 use ua_data
   553:                 implicit none
   554:           
   555:                 double precision tmor,tx(lx1,lx1,lx1),tmp
   556:                 integer i,n
   557:           
   558:                 tmor=tx(1,1,1)
   559:           
   560: *------>        do i=2,lx1-1
   561: |                 tmor= tmor + qbnew(i-1,1,1)*tx(i,1,1)
   562: *------         end do
   563:           
   564:                 if(n.gt.1)then
   565: *------>          do i=2,lx1-1
   566: |                   tmor= tmor + qbnew(i-1,1,1)*tx(1,i,1)
   567: *------           end do
   568:                 end if
   569:           
   570:                 if(n.eq.3)then
   571: *------>          do i=2,lx1-1
   572: |                   tmor= tmor + qbnew(i-1,1,1)*tx(1,1,i)
   573: *------           end do
   574:                 end if
   575:           
   576:                 return
   577:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_COR_F
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   599: opt(1418): Constant-length loop is expanded.
   601: opt(1418): Constant-length loop is expanded.
   607: opt(1418): Constant-length loop is expanded.
   612: opt(1418): Constant-length loop is expanded.
   619: opt(1418): Constant-length loop is expanded.
   620: opt(1418): Constant-length loop is expanded.
   625: opt(1418): Constant-length loop is expanded.
   632: opt(1418): Constant-length loop is expanded.
   633: opt(1418): Constant-length loop is expanded.
   637: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_COR_F
FORMAT LIST

 LINE   LOOP      STATEMENT

   580:                 subroutine transfb_cor_f(n,tmor,tx)
   581:           c--------------------------------------------------------------
   582:           c     This subroutine performs the mapping from face to mortar.
   583:           c     Output tmor is the mapping result on a mortar vertex
   584:           c     of situations of three edges and three faces sharing a vertex:
   585:           c     n=4: only one face is nonconforming
   586:           c     n=5: one face and one edge are nonconforming
   587:           c     n=6: two faces are nonconforming
   588:           c     n=7: three faces are nonconforming
   589:           c--------------------------------------------------------------
   590:           
   591:                 use ua_data
   592:                 implicit none
   593:           
   594:                 double precision tx(lx1,lx1,lx1),tmor,temp(lx1)
   595:                 integer col,i,n
   596:           
   597:                 call r_init(temp,lx1,0.d0)
   598:           
   599: *------>        do col=1,lx1
   600: |                 temp(col)=tx(col,1,1)
   601: |*----->          do i=2,lx1-1
   602: ||                  temp(col) = temp(col) + qbnew(i-1,1,1)*tx(col,i,1)
   603: |*-----           end do
   604: *------         end do
   605:                 tmor=temp(1)
   606:           
   607: *------>        do i=2,lx1-1
   608: |                 tmor = tmor + qbnew(i-1,1,1) *temp(i)
   609: *------         end do
   610:           
   611:                 if(n.eq.5)then
   612: *------>          do i=2,lx1-1
   613: |                   tmor = tmor + qbnew(i-1,1,1) *tx(1,1,i)
   614: *------           end do
   615:                 end if
   616:            
   617:                 if(n.ge.6)then
   618:                   call r_init(temp,lx1,0.d0)
   619: *------>          do col=1,lx1
   620: |*----->            do i=2,lx1-1
   621: ||                    temp(col) = temp(col) + qbnew(i-1,1,1)*tx(col,1,i)
   622: |*-----             end do
   623: *------           end do
   624:                   tmor=tmor+temp(1)
   625: *------>          do i=2,lx1-1
   626: |                   tmor = tmor +qbnew(i-1,1,1) *temp(i)
   627: *------           end do
   628:                 end if
   629:            
   630:                 if(n.eq.7)then
   631:                   call r_init(temp,lx1,0.d0)
   632: *------>          do col=2,lx1-1
   633: |*----->            do i=2,lx1-1
   634: ||                    temp(col) = temp(col) + qbnew(i-1,1,1)*tx(1,col,i)
   635: |*-----             end do
   636: *------           end do
   637: *------>          do i=2,lx1-1
   638: |                   tmor = tmor + qbnew(i-1,1,1) *temp(i)
   639: *------           end do
   640:                 end if
   641:           
   642:                 return
   643:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSF_NC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   663: vec( 101): Vectorized loop.
   666: opt(1418): Constant-length loop is expanded.
   667: opt(1418): Constant-length loop is expanded.
   668: vec( 128): Fused multiply-add operation applied.
   673: vec( 101): Vectorized loop.
   676: opt(1418): Constant-length loop is expanded.
   677: opt(1418): Constant-length loop is expanded.
   678: vec( 128): Fused multiply-add operation applied.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSF_NC
FORMAT LIST

 LINE   LOOP      STATEMENT

   647:                 subroutine transf_nc(tmor,tx)
   648:           c------------------------------------------------------------------------
   649:           c     Perform mortar to element mapping on a nonconforming face.
   650:           c     This subroutin is used when all entries in tmor are zero except
   651:           c     one tmor(i,j)=1. So this routine is simplified. Only one piece of
   652:           c     mortar  (tmor only has two indices) and one piece of intermediate
   653:           c     mortar (tmp) are involved.
   654:           c------------------------------------------------------------------------
   655:           
   656:                 use ua_data
   657:                 implicit none
   658:           
   659:                 double precision tmor(lx1,lx1), tx(lx1,lx1), tmp(lx1,lx1)
   660:                 integer col,i,j
   661:           
   662:                 call r_init(tmp,lx1*lx1,0.d0)
   663: V------>        do col=1,lx1
   664: |                 i = 1
   665: |                 tmp(i,col)=tmor(i,col)
   666: |*----->          do i=2,lx1-1
   667: ||*---->            do j=1,lx1
   668: |||     F             tmp(i,col) = tmp(i,col) + qbnew(i-1,j,1)*tmor(j,col)
   669: ||*----             end do
   670: |*-----           end do
   671: V------         end do
   672:           
   673: V------>        do col=1,lx1
   674: |                 i = 1
   675: |                 tx(col,i)   = tx(col,i)   + tmp(col,i)
   676: |*----->          do i=2,lx1-1
   677: ||*---->            do j=1,lx1
   678: |||     F             tx(col,i) = tx(col,i) + qbnew(i-1,j,1)*tmp(col,j)
   679: ||*----             end do
   680: |*-----           end do
   681: V------         end do
   682:           
   683:                 return
   684:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC0
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   700: opt(1418): Constant-length loop is expanded.
   701: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC0
FORMAT LIST

 LINE   LOOP      STATEMENT

   687:                 subroutine transfb_nc0(tmor,tx)
   688:           c------------------------------------------------------------------------
   689:           c     Performs mapping from element to mortar when the nonconforming
   690:           c     edges are shared by two conforming faces of an element.
   691:           c------------------------------------------------------------------------
   692:           
   693:                 use ua_data
   694:                 implicit none
   695:           
   696:                 double precision tmor(lx1,lx1),tx(lx1,lx1,lx1)
   697:                 integer i,j
   698:           
   699:                 call r_init(tmor,lx1*lx1,0.d0)
   700: *------>        do j=1,lx1
   701: |*----->          do i=2,lx1-1
   702: ||                  tmor(j,1)= tmor(j,1) + qbnew(i-1,j  ,1)*tx(i,1,1)
   703: |*-----           end do
   704: *------         end do
   705:           
   706:                 return
   707:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   731: vec( 101): Vectorized loop.
   735: opt(1418): Constant-length loop is expanded.
   736: vec( 128): Fused multiply-add operation applied.
   739: opt(1418): Constant-length loop is expanded.
   740: opt(1418): Constant-length loop is expanded.
   741: vec( 128): Fused multiply-add operation applied.
   751: opt(1418): Constant-length loop is expanded.
   752: opt(1418): Constant-length loop is expanded.
   758: vec( 103): Unvectorized loop.
   758: vec( 107): Iteration count is too small.
   760: opt(1418): Constant-length loop is expanded.
   761: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC2
FORMAT LIST

 LINE   LOOP      STATEMENT

   710:                 subroutine transfb_nc2(tmor,tx)
   711:           c------------------------------------------------------------------------
   712:           c     Maps values from element to mortar when the nonconforming edges are
   713:           c     shared by two nonconforming faces of an element.
   714:           c     Although each face shall have four pieces of mortar, only value in
   715:           c     one piece (location (1,1)) is used in the calling routine so only
   716:           c     the value in the first mortar is calculated in this subroutine.
   717:           c------------------------------------------------------------------------
   718:           
   719:                 use ua_data
   720:                 implicit none
   721:           
   722:                 double precision tx(lx1,lx1),tmor(lx1,lx1),bottom(lx1),
   723:                &                 temp(lx1,lx1)
   724:                 integer col,j,i
   725:           
   726:                 call r_init(tmor,lx1*lx1,0.d0)
   727:                 call r_init(temp,lx1*lx1,0.d0)
   728:                 tmor(1,1)=tx(1,1)
   729:           
   730:           c.....mapping from tx to intermediate mortar temp + bottom
   731: V------>        do col=1,lx1
   732: |                 temp(col,1)=tx(col,1)
   733: |                 j=1
   734: |                 bottom(col)= 0.d0
   735: |*----->          do i=2,lx1-1
   736: ||      F           bottom(col) = bottom(col) + qbnew(i-1,j,1)*tx(col,i)
   737: |*-----           end do
   738: |         
   739: |*----->          do j=2,lx1
   740: ||*---->            do i=2,lx1-1
   741: |||     F             temp(col,j) = temp(col,j) + qbnew(i-1,j,1)*tx(col,i)
   742: ||*----             end do
   743: |*-----           end do
   744: V------         end do
   745:           
   746:           c.....from intermediate mortar to mortar
   747:           
   748:           c.....On the nonconforming edge, temp is divided by 2 as there will be
   749:           c     a duplicate contribution from another face sharing this edge
   750:                 col=1
   751: *------>        do j=1,lx1
   752: |*----->          do i=2,lx1-1
   753: ||                  tmor(j,col)=tmor(j,col)+ qbnew(i-1,j,1) * bottom(i) +
   754: ||             &                             qbnew(i-1,j,1) * temp(i,col) * 0.5d0
   755: |*-----           end do
   756: *------         end do
   757:           
   758: +------>        do col=2,lx1
   759: |                 tmor(1,col)=tmor(1,col)+temp(1,col)
   760: |*----->          do j=1,lx1
   761: ||*---->            do i=2,lx1-1
   762: |||                   tmor(j,col) = tmor(j,col) + qbnew(i-1,j,1) *temp(i,col)
   763: ||*----             end do
   764: |*-----           end do
   765: +------         end do
   766:           
   767:                 return
   768:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC1
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   793: vec( 101): Vectorized loop.
   797: opt(1418): Constant-length loop is expanded.
   798: vec( 128): Fused multiply-add operation applied.
   801: opt(1418): Constant-length loop is expanded.
   802: opt(1418): Constant-length loop is expanded.
   803: vec( 128): Fused multiply-add operation applied.
   811: opt(1418): Constant-length loop is expanded.
   812: opt(1418): Constant-length loop is expanded.
   822: vec( 103): Unvectorized loop.
   822: vec( 107): Iteration count is too small.
   824: opt(1418): Constant-length loop is expanded.
   825: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_NC1
FORMAT LIST

 LINE   LOOP      STATEMENT

   772:                 subroutine transfb_nc1(tmor,tx)
   773:           c------------------------------------------------------------------------
   774:           c     Maps values from element to mortar when the nonconforming edges are
   775:           c     shared by a nonconforming face and a conforming face of an element
   776:           c------------------------------------------------------------------------
   777:           
   778:                 use ua_data
   779:                 implicit none
   780:           
   781:                 double precision tx(lx1,lx1),tmor(lx1,lx1),bottom(lx1),
   782:                &                 temp(lx1,lx1)
   783:                 integer col,j,i
   784:           
   785:                 call r_init(tmor,lx1*lx1,0.d0)
   786:                 call r_init(temp,lx1*lx1,0.d0)
   787:           
   788:                 tmor(1,1)=tx(1,1)
   789:           c.....Contribution from the nonconforming faces
   790:           c     Since the calling subroutine is only interested in the value on the
   791:           c     mortar (location (1,1)), only this piece of mortar is calculated.
   792:           
   793: V------>        do col=1,lx1
   794: |                 temp(col,1)=tx(col,1)
   795: |                 j = 1
   796: |                 bottom(col)= 0.d0
   797: |*----->          do i=2,lx1-1
   798: ||      F           bottom(col)=bottom(col) + qbnew(i-1,j,1)*tx(col,i)
   799: |*-----           end do
   800: |         
   801: |*----->          do j=2,lx1
   802: ||*---->            do i=2,lx1-1
   803: |||     F             temp(col,j) = temp(col,j) + qbnew(i-1,j,1)*tx(col,i)
   804: ||*----             end do
   805: ||        
   806: |*-----           end do
   807: V------         end do
   808:           
   809:                 col=1
   810:                 tmor(1,col)=tmor(1,col)+bottom(1)
   811: *------>        do j=1,lx1
   812: |*----->          do i=2,lx1-1
   813: ||        
   814: ||        c.........temp is not divided by 2 here. It includes the contribution
   815: ||        c         from the other conforming face.
   816: ||        
   817: ||                  tmor(j,col)=tmor(j,col) + qbnew(i-1,j,1) *bottom(i) +
   818: ||             &                              qbnew(i-1,j,1) *temp(i,col)
   819: |*-----           end do
   820: *------         end do
   821:           
   822: +------>        do col=2,lx1
   823: |                 tmor(1,col)=tmor(1,col)+temp(1,col)
   824: |*----->          do j=1,lx1
   825: ||*---->            do i=2,lx1-1
   826: |||                   tmor(j,col) = tmor(j,col) + qbnew(i-1,j,1) *temp(i,col)
   827: ||*----             end do
   828: |*-----           end do
   829: +------         end do
   830:           
   831:                 return
   832:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_C
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   850: par(1801): Parallel routine generated.: TRANSFB_C$1
   854: par(1803): Parallelized by "do".
   854: vec( 101): Vectorized loop.
   860: par(1803): Parallelized by "do".
   872: par(3007): Atomic-construct generated.
   874: par(3007): Atomic-construct generated.
   876: par(3007): Atomic-construct generated.
   878: par(3007): Atomic-construct generated.
   881: vec( 103): Unvectorized loop.
   881: vec( 109): Vectorization obstructive statement.
   882: opt(1418): Constant-length loop is expanded.
   885: par(3007): Atomic-construct generated.
   891: opt(1418): Constant-length loop is expanded.
   894: par(3007): Atomic-construct generated.
   900: opt(1418): Constant-length loop is expanded.
   903: par(3007): Atomic-construct generated.
   909: opt(1418): Constant-length loop is expanded.
   912: par(3007): Atomic-construct generated.
   918: opt(1418): Constant-length loop is expanded.
   921: par(3007): Atomic-construct generated.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_C
FORMAT LIST

 LINE   LOOP      STATEMENT

   836:                 subroutine transfb_c(tx)
   837:           c-------------------------------------------------------------------
   838:           c     Prepare initial guess for cg. All values from conforming
   839:           c     boundary are copied and summed on tmor.
   840:           c-------------------------------------------------------------------
   841:           
   842:                 use ua_data
   843:                 implicit none
   844:           
   845:                 double precision third
   846:                 parameter (third = 1.d0/3.d0)
   847:                 double precision tx(*)
   848:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,col,j,ig,il
   849:           
   850:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IE,IFACE,IL1,IL2,
   851:           c$OMP& IL3,IL4,IG1,IG2,IG3,IG4,COL,J,IG,IL)
   852:           
   853:           c$OMP DO
   854: Y------>        do j=1,nmor
   855: |                 tmort(j)=0.d0
   856: Y------         end do
   857:           c$OMP END DO
   858:           
   859:           c$OMP DO
   860: P------>        do ie=1,nelt
   861: |+----->          do iface=1,nsides
   862: ||                  if(cbc(iface,ie).ne.3)then
   863: ||                    il1 = idel(1,1,iface,ie)
   864: ||                    il2 = idel(lx1,1,iface,ie)
   865: ||                    il3 = idel(1,lx1,iface,ie)
   866: ||                    il4 = idel(lx1,lx1,iface,ie)
   867: ||                    ig1 = idmo(1,  1,  1,1,iface,ie)
   868: ||                    ig2 = idmo(lx1,1,  1,2,iface,ie)
   869: ||                    ig3 = idmo(1,  lx1,2,1,iface,ie)
   870: ||                    ig4 = idmo(lx1,lx1,2,2,iface,ie)
   871: ||        
   872: ||        c$OMP ATOMIC
   873: ||                    tmort(ig1) = tmort(ig1)+tx(il1)*third
   874: ||        c$OMP ATOMIC
   875: ||                    tmort(ig2) = tmort(ig2)+tx(il2)*third
   876: ||        c$OMP ATOMIC
   877: ||                    tmort(ig3) = tmort(ig3)+tx(il3)*third
   878: ||        c$OMP ATOMIC
   879: ||                    tmort(ig4) = tmort(ig4)+tx(il4)*third
   880: ||        
   881: ||+---->              do  col=2,lx1-1
   882: |||*--->                do j=2,lx1-1
   883: ||||                      il=idel(j,col,iface,ie)
   884: ||||                      ig=idmo(j,col,1,1,iface,ie)
   885: ||||      c$OMP ATOMIC
   886: ||||                      tmort(ig)=tmort(ig)+tx(il)
   887: |||*---                 end do
   888: ||+----               end do
   889: ||        
   890: ||                    if(idmo(lx1,1,1,1,iface,ie).eq.0)then
   891: ||*---->                do j=2,lx1-1
   892: |||                       il=idel(j,1,iface,ie)
   893: |||                       ig=idmo(j,1,1,1,iface,ie)
   894: |||       c$OMP ATOMIC
   895: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   896: ||*----                 end do
   897: ||                    end if
   898: ||        
   899: ||                    if(idmo(lx1,2,1,2,iface,ie).eq.0)then
   900: ||*---->                do j=2,lx1-1
   901: |||                       il=idel(lx1,j,iface,ie)
   902: |||                       ig=idmo(lx1,j,1,1,iface,ie)
   903: |||       c$OMP ATOMIC
   904: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   905: ||*----                 end do
   906: ||                    end if
   907: ||        
   908: ||                    if(idmo(2,lx1,2,1,iface,ie).eq.0)then
   909: ||*---->                do j=2,lx1-1
   910: |||                       il=idel(j,lx1,iface,ie)
   911: |||                       ig=idmo(j,lx1,1,1,iface,ie)
   912: |||       c$OMP ATOMIC
   913: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   914: ||*----                 end do
   915: ||                    end if
   916: ||        
   917: ||                    if(idmo(1,lx1,1,1,iface,ie).eq.0)then
   918: ||*---->                do j=2,lx1-1
   919: |||                       il=idel(1,j,iface,ie)
   920: |||                       ig=idmo(1,j,1,1,iface,ie)
   921: |||       c$OMP ATOMIC
   922: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
   923: ||*----                 end do
   924: ||                    end if
   925: ||                  end if!
   926: |+-----           end do
   927: P------         end do
   928:           c$OMP END DO NOWAIT
   929:           c$OMP END PARALLEL
   930:                 return
   931:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_C_2
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   949: par(1801): Parallel routine generated.: TRANSFB_C_2$1
   953: par(1803): Parallelized by "do".
   953: vec( 101): Vectorized loop.
   958: par(1803): Parallelized by "do".
   958: vec( 101): Vectorized loop.
   964: par(1803): Parallelized by "do".
   977: par(3007): Atomic-construct generated.
   979: par(3007): Atomic-construct generated.
   981: par(3007): Atomic-construct generated.
   983: par(3007): Atomic-construct generated.
   985: par(3007): Atomic-construct generated.
   987: par(3007): Atomic-construct generated.
   989: par(3007): Atomic-construct generated.
   991: par(3007): Atomic-construct generated.
   994: vec( 103): Unvectorized loop.
   994: vec( 109): Vectorization obstructive statement.
   995: opt(1418): Constant-length loop is expanded.
   998: par(3007): Atomic-construct generated.
  1000: par(3007): Atomic-construct generated.
  1006: opt(1418): Constant-length loop is expanded.
  1009: par(3007): Atomic-construct generated.
  1011: par(3007): Atomic-construct generated.
  1017: opt(1418): Constant-length loop is expanded.
  1020: par(3007): Atomic-construct generated.
  1022: par(3007): Atomic-construct generated.
  1028: opt(1418): Constant-length loop is expanded.
  1031: par(3007): Atomic-construct generated.
  1033: par(3007): Atomic-construct generated.
  1039: opt(1418): Constant-length loop is expanded.
  1042: par(3007): Atomic-construct generated.
  1044: par(3007): Atomic-construct generated.


NEC Fortran Compiler (3.0.4) for Vector Engine   Mon Jul 20 02:51:54 2020
FILE NAME: transfer_au.f

PROCEDURE NAME: TRANSFB_C_2
FORMAT LIST

 LINE   LOOP      STATEMENT

   934:                 subroutine transfb_c_2(tx)
   935:           c-------------------------------------------------------------------
   936:           c     Prepare initial guess for CG. All values from conforming
   937:           c     boundary are copied and summed in tmort.
   938:           c     mormult is multiplicity, which is used to average tmort.
   939:           c-------------------------------------------------------------------
   940:           
   941:                 use ua_data
   942:                 implicit none
   943:           
   944:                 double precision third
   945:                 parameter (third = 1.d0/3.d0)
   946:                 double precision tx(*)
   947:                 integer il1,il2,il3,il4,ig1,ig2,ig3,ig4,ie,iface,col,j,ig,il
   948:           
   949:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IE,IFACE,IL1,IL2,
   950:           c$OMP& IL3,IL4,IG1,IG2,IG3,IG4,COL,J,IG,IL)
   951:           
   952:           c$OMP DO
   953: Y------>        do j=1,nmor
   954: |                 tmort(j)=0.d0
   955: Y------         end do
   956:           c$OMP END DO nowait
   957:           c$OMP DO
   958: Y------>        do j=1,nmor
   959: |                 mormult(j)=0.d0
   960: Y------         end do
   961:           c$OMP END DO
   962:           
   963:           c$OMP DO
   964: P------>        do ie=1,nelt
   965: |+----->          do iface=1,nsides
   966: ||         
   967: ||                  if(cbc(iface,ie).ne.3)then
   968: ||                    il1 = idel(1,  1,  iface,ie)
   969: ||                    il2 = idel(lx1,1,  iface,ie)
   970: ||                    il3 = idel(1,  lx1,iface,ie)
   971: ||                    il4 = idel(lx1,lx1,iface,ie)
   972: ||                    ig1 = idmo(1,  1,  1,1,iface,ie)
   973: ||                    ig2 = idmo(lx1,1,  1,2,iface,ie)
   974: ||                    ig3 = idmo(1,  lx1,2,1,iface,ie)
   975: ||                    ig4 = idmo(lx1,lx1,2,2,iface,ie)
   976: ||        
   977: ||        c$OMP ATOMIC
   978: ||                    tmort(ig1) = tmort(ig1)+tx(il1)*third
   979: ||        c$OMP ATOMIC
   980: ||                    tmort(ig2) = tmort(ig2)+tx(il2)*third
   981: ||        c$OMP ATOMIC
   982: ||                    tmort(ig3) = tmort(ig3)+tx(il3)*third
   983: ||        c$OMP ATOMIC
   984: ||                    tmort(ig4) = tmort(ig4)+tx(il4)*third
   985: ||        c$OMP ATOMIC
   986: ||                    mormult(ig1) = mormult(ig1)+third
   987: ||        c$OMP ATOMIC
   988: ||                    mormult(ig2) = mormult(ig2)+third
   989: ||        c$OMP ATOMIC
   990: ||                    mormult(ig3) = mormult(ig3)+third
   991: ||        c$OMP ATOMIC
   992: ||                    mormult(ig4) = mormult(ig4)+third
   993: ||        
   994: ||+---->              do  col=2,lx1-1
   995: |||*--->                do j=2,lx1-1
   996: ||||                      il=idel(j,col,iface,ie)
   997: ||||                      ig=idmo(j,col,1,1,iface,ie)
   998: ||||      c$OMP ATOMIC
   999: ||||                      tmort(ig)=tmort(ig)+tx(il)
  1000: ||||      c$OMP ATOMIC
  1001: ||||                      mormult(ig)=mormult(ig)+1.d0
  1002: |||*---                 end do
  1003: ||+----               end do
  1004: ||        
  1005: ||                    if(idmo(lx1,1,1,1,iface,ie).eq.0)then
  1006: ||*---->                do j=2,lx1-1
  1007: |||                       il=idel(j,1,iface,ie)
  1008: |||                       ig=idmo(j,1,1,1,iface,ie)
  1009: |||       c$OMP ATOMIC
  1010: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1011: |||       c$OMP ATOMIC
  1012: |||                       mormult(ig)=mormult(ig)+0.5d0
  1013: ||*----                 end do
  1014: ||                    end if
  1015: ||        
  1016: ||                    if(idmo(lx1,2,1,2,iface,ie).eq.0)then
  1017: ||*---->                do j=2,lx1-1
  1018: |||                       il=idel(lx1,j,iface,ie)
  1019: |||                       ig=idmo(lx1,j,1,1,iface,ie)
  1020: |||       c$OMP ATOMIC
  1021: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1022: |||       c$OMP ATOMIC
  1023: |||                       mormult(ig)=mormult(ig)+0.5d0
  1024: ||*----                 end do
  1025: ||                    end if
  1026: ||        
  1027: ||                    if(idmo(2,lx1,2,1,iface,ie).eq.0)then
  1028: ||*---->                do j=2,lx1-1
  1029: |||                       il=idel(j,lx1,iface,ie)
  1030: |||                       ig=idmo(j,lx1,1,1,iface,ie)
  1031: |||       c$OMP ATOMIC
  1032: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1033: |||       c$OMP ATOMIC
  1034: |||                       mormult(ig)=mormult(ig)+0.5d0
  1035: ||*----                  end do
  1036: ||                    end if
  1037: ||        
  1038: ||                    if(idmo(1,lx1,1,1,iface,ie).eq.0)then
  1039: ||*---->                do j=2,lx1-1
  1040: |||                       il=idel(1,j,iface,ie)
  1041: |||                       ig=idmo(1,j,1,1,iface,ie)
  1042: |||       c$OMP ATOMIC
  1043: |||                       tmort(ig)=tmort(ig)+tx(il)*0.5d0
  1044: |||       c$OMP ATOMIC
  1045: |||                       mormult(ig)=mormult(ig)+0.5d0
  1046: ||*----                 end do
  1047: ||                    end if
  1048: ||                  end if!nnje=1
  1049: |+-----           end do
  1050: P------         end do
  1051:           c$OMP END DO NOWAIT
  1052:           c$OMP END PARALLEL
  1053:           
  1054:                 return
  1055:                 end


