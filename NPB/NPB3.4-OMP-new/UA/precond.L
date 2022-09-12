NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETUPPC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    16: opt(1418): Constant-length loop is expanded.
    17: opt(1418): Constant-length loop is expanded.
    24: par(1801): Parallel routine generated.: SETUPPC$1
    25: vec( 110): Vectorization obstructive procedure reference.: R_INIT
    25: par(1803): Parallelized by "do".
    29: vec( 101): Vectorized loop.
    30: opt(1418): Constant-length loop is expanded.
    31: opt(1418): Constant-length loop is expanded.
    32: vec( 128): Fused multiply-add operation applied.
    37: vec( 128): Fused multiply-add operation applied.
    54: par(1801): Parallel routine generated.: SETUPPC$2
    55: par(1803): Parallelized by "do".
    55: vec( 101): Vectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETUPPC
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine setuppc
     3:           c------------------------------------------------------------------
     4:           c     Generate diagonal preconditioner for CG.
     5:           c     Preconditioner computed in this subroutine is correct only
     6:           c     for collocation point in element interior, on conforming face
     7:           c     interior and conforming edge.
     8:           c------------------------------------------------------------------
     9:           
    10:                 use ua_data
    11:                 implicit none
    12:           
    13:                 double precision dxtm1_2(lx1,lx1), rdtime
    14:                 integer ie,k,i,j,q,isize
    15:           
    16: *------>        do j=1,lx1
    17: |*----->          do i=1,lx1
    18: ||                  dxtm1_2(i,j)=dxtm1(i,j)**2
    19: |*-----           end do
    20: *------         end do
    21:           
    22:                 rdtime=1.d0/dtime
    23:           
    24:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ie,isize,i,j,k,q)
    25: P------>        do ie = 1, nelt
    26: |                 call r_init(dpcelm(1,1,1,ie),nxyz,0.d0)
    27: |                 isize=size_e(ie)
    28: |+----->          do k = 1, lx1
    29: ||V---->            do j = 1, lx1
    30: |||*--->              do i = 1, lx1
    31: ||||*-->                do q = 1, lx1
    32: |||||   F                 dpcelm(i,j,k,ie) = dpcelm(i,j,k,ie) +
    33: |||||          &                        g1m1_s(q,j,k,isize) * dxtm1_2(i,q) +
    34: |||||          &                        g1m1_s(i,q,k,isize) * dxtm1_2(j,q) +
    35: |||||          &                        g1m1_s(i,j,q,isize) * dxtm1_2(k,q)
    36: ||||*--                 end do
    37: ||||    F               dpcelm(i,j,k,ie)=visc*dpcelm(i,j,k,ie)+
    38: ||||           &                      rdtime*bm1_s(i,j,k,isize)
    39: |||*---               end do
    40: ||V----             end do
    41: |+-----           end do
    42: P------         end do
    43:           c$OMP END PARALLEL DO
    44:           
    45:           c.....do the stiffness summation
    46:                 call dssum
    47:           
    48:           c.....take inverse.
    49:           
    50:                 call reciprocal(dpcelm,ntot)
    51:           
    52:           c.....compute preconditioner on mortar points. NOTE:  dpcmor for
    53:           c     nonconforming cases will be corrected in subroutine setpcmo
    54:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    55: Y------>        do i=1,nmor
    56: |                 dpcmor(i)=1.d0/dpcmor(i)
    57: Y------         end do
    58:           c$OMP END PARALLEL DO
    59:           
    60:                 return
    61:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETPCMO_PRE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    90: opt(1418): Constant-length loop is expanded.
    91: opt(1418): Constant-length loop is expanded.
    96: vec( 101): Vectorized loop.
    99: opt(1418): Constant-length loop is expanded.
   100: opt(1418): Constant-length loop is expanded.
   101: vec( 128): Fused multiply-add operation applied.
   107: par(1801): Parallel routine generated.: SETPCMO_PRE$1
   109: par(1803): Parallelized by "do".
   119: vec( 103): Unvectorized loop.
   119: vec( 110): Vectorization obstructive procedure reference.: LAPLACIAN
   119: vec( 110): Vectorization obstructive procedure reference.: R_INIT
   120: opt(1025): Reference to this procedure inhibits optimization.: R_INIT
   143: vec( 103): Unvectorized loop.
   143: vec( 110): Vectorization obstructive procedure reference.: LAPLACIAN
   143: vec( 110): Vectorization obstructive procedure reference.: TRANSFB_NC1
   143: vec( 110): Vectorization obstructive procedure reference.: TRANSF_NC
   143: vec( 110): Vectorization obstructive procedure reference.: R_INIT
   144: opt(1025): Reference to this procedure inhibits optimization.: R_INIT
   173: vec( 103): Unvectorized loop.
   173: vec( 110): Vectorization obstructive procedure reference.: TRANSFB_NC1
   173: vec( 110): Vectorization obstructive procedure reference.: TRANSFB_NC2
   173: vec( 110): Vectorization obstructive procedure reference.: LAPLACIAN
   173: vec( 110): Vectorization obstructive procedure reference.: TRANSFB_NC0
   173: vec( 110): Vectorization obstructive procedure reference.: TRANSF_NC
   173: vec( 110): Vectorization obstructive procedure reference.: R_INIT
   174: opt(1025): Reference to this procedure inhibits optimization.: R_INIT
   184: opt(1418): Constant-length loop is expanded.
   188: opt(1418): Constant-length loop is expanded.
   230: opt(1418): Constant-length loop is expanded.
   239: opt(1418): Constant-length loop is expanded.
   240: opt(1418): Constant-length loop is expanded.
   261: opt(1418): Constant-length loop is expanded.
   269: opt(1418): Constant-length loop is expanded.
   276: opt(1418): Constant-length loop is expanded.
   294: opt(1418): Constant-length loop is expanded.
   295: opt(1418): Constant-length loop is expanded.
   308: opt(1418): Constant-length loop is expanded.
   336: opt(1418): Constant-length loop is expanded.
   344: opt(1418): Constant-length loop is expanded.
   353: opt(1418): Constant-length loop is expanded.
   363: opt(1418): Constant-length loop is expanded.
   364: opt(1418): Constant-length loop is expanded.
   373: opt(1418): Constant-length loop is expanded.
   374: opt(1418): Constant-length loop is expanded.
   378: opt(1418): Constant-length loop is expanded.
   386: vec( 103): Unvectorized loop.
   386: vec( 113): Overhead of loop division is too large.
   387: opt(1418): Constant-length loop is expanded.
   388: opt(1034): Multiple store conflict.: P0
   388: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: P0$77
   396: vec( 103): Unvectorized loop.
   396: vec( 113): Overhead of loop division is too large.
   397: opt(1418): Constant-length loop is expanded.
   398: opt(1034): Multiple store conflict.: P0
   398: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: P0$77


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETPCMO_PRE
FORMAT LIST

 LINE   LOOP      STATEMENT

    65:                 subroutine setpcmo_pre
    66:           c--------------------------------------------------------------
    67:           c     pre-compute elemental contribution to preconditioner
    68:           c     for all situations
    69:           c--------------------------------------------------------------
    70:            
    71:                 use ua_data
    72:                 implicit none
    73:           
    74:                 integer element_size, i, j, ii, jj, col
    75:                 double precision
    76:                &       p(lx1,lx1,lx1), p0(lx1,lx1,lx1), mtemp(lx1,lx1),
    77:                &       temp(lx1,lx1,lx1), temp1(lx1,lx1), tmp(lx1,lx1),tig(lx1)
    78:           
    79:           c.....corners on face of type 3
    80:           
    81:                 call r_init(tcpre,lx1*lx1,0.d0)
    82:                 call r_init(tmp,lx1*lx1,0.d0)
    83:                 call r_init(tig,5,0.d0)
    84:                 tig(1)   =1.d0
    85:                 tmp(1,1) =1.d0
    86:           
    87:           c.....tcpre results from mapping a unit spike field (unity at
    88:           c     collocation point (1,1), zero elsewhere) on an entire element
    89:           c     face to the (1,1) segment of a nonconforming face
    90: *------>        do i=2,lx1-1
    91: |*----->          do j=1,lx1
    92: ||                  tmp(i,1) = tmp(i,1)+ qbnew(i-1,j,1)*tig(j)
    93: |*-----           end do
    94: *------         end do
    95:            
    96: V------>        do col=1,lx1
    97: |                 tcpre(col,1)=tmp(col,1)
    98: |         
    99: |*----->          do j=2,lx1-1
   100: ||*---->            do i=1,lx1
   101: |||     F             tcpre(col,j) = tcpre(col,j) + qbnew(j-1,i,1)*
   102: |||            &                                     tmp(col,i)
   103: ||*----             end do
   104: |*-----           end do
   105: V------         end do
   106:           
   107:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(element_size,i,j,p,temp,
   108:           c$OMP& mtemp,temp1,p0,ii,jj)
   109: P------>        do element_size=1,refine_max
   110: |         
   111: |         c.......for conforming cases
   112: |         
   113: |         c.......pcmor_c (i,j,element_size) records the intermediate value
   114: |         c       (preconditioner=1/pcmor_c) of the preconditor on collocation
   115: |         c       point (i,j) on a conforming face of an element of size
   116: |         c       element_size.
   117: |         
   118: |+----->          do j=1,lx1/2+1
   119: ||+---->            do i=j,lx1/2+1
   120: |||                   call r_init(p,nxyz,0.d0)
   121: |||                   p(i,j,1)=1.d0
   122: |||                   call laplacian(temp,p,element_size)
   123: |||                   pcmor_c(i,j,element_size)=temp(i,j,1)
   124: |||                   pcmor_c(lx1+1-i,j,element_size)=temp(i,j,1)
   125: |||                   pcmor_c(j,i,element_size)=temp(i,j,1)
   126: |||                   pcmor_c(lx1+1-j,i,element_size)=temp(i,j,1)
   127: |||                   pcmor_c(j,lx1+1-i,element_size)=temp(i,j,1)
   128: |||                   pcmor_c(lx1+1-j,lx1+1-i,element_size)=temp(i,j,1)
   129: |||                   pcmor_c(i,lx1+1-j,element_size)=temp(i,j,1)
   130: |||                   pcmor_c(lx1+1-i,lx1+1-j,element_size)=temp(i,j,1)
   131: ||+----             end do
   132: |+-----           end do
   133: |         
   134: |         c.......for nonconforming cases
   135: |         
   136: |         c.......nonconforming face interior
   137: |         
   138: |         c.......pcmor_nc1(i,j,ii,jj,element_size) records the intermediate
   139: |         c       preconditioner value on collocation point (i,j) on mortar
   140: |         c       (ii,jj)  on a nonconforming face of an element of size element_
   141: |         c       size
   142: |+----->          do j=2,lx1
   143: ||+---->            do i=j,lx1
   144: |||                   call r_init(mtemp,lx1*lx1,0.d0)
   145: |||                   call r_init(p,nxyz,0.d0)
   146: |||                   mtemp(i,j)=1.d0
   147: |||       c...........when i, j=lx1, mortar points are duplicated, so mtemp needs
   148: |||       c           to be doubled.
   149: |||                   if(i.eq.lx1)mtemp(i,j)=mtemp(i,j)*2.d0
   150: |||                   if(j.eq.lx1)mtemp(i,j)=mtemp(i,j)*2.d0
   151: |||                   call transf_nc(mtemp,p)
   152: |||                   call laplacian(temp,p,element_size)
   153: |||                   call transfb_nc1(temp1,temp)
   154: |||       
   155: |||       c...........values at points (i,j) and (j,i) are the same
   156: |||                   pcmor_nc1(i,j,1,1,element_size)=temp1(i,j)
   157: |||                   pcmor_nc1(j,i,1,1,element_size)=temp1(i,j)
   158: ||+----             end do
   159: ||        
   160: ||        c.........when i, j=lx1, mortar points are duplicated. so pcmor_nc1 needs
   161: ||        c         to be doubled on those points
   162: ||                  pcmor_nc1(lx1,j,1,1,element_size)=
   163: ||             &          pcmor_nc1(lx1,j,1,1,element_size)*2.d0
   164: ||                  pcmor_nc1(j,lx1,1,1,element_size)=
   165: ||             &          pcmor_nc1(lx1,j,1,1,element_size)
   166: ||        
   167: |+-----           end do
   168: |                 pcmor_nc1(lx1,lx1,1,1,element_size)=
   169: |              &      pcmor_nc1(lx1,lx1,1,1,element_size)*2.d0
   170: |         
   171: |         c.......nonconforming edges
   172: |                 j=1
   173: |+----->          do i=2,lx1
   174: ||                  call r_init(mtemp,lx1*lx1,0.d0)
   175: ||                  call r_init(p,nxyz,0.d0)
   176: ||                  call r_init(p0,nxyz,0.d0)
   177: ||                  mtemp(i,j)=1.d0
   178: ||                  if(i.eq.lx1)mtemp(i,j)=2.d0
   179: ||                  call transf_nc(mtemp,p)
   180: ||                  call laplacian(temp,p,element_size)
   181: ||                  call transfb_nc1(temp1,temp)
   182: ||                  pcmor_nc1(i,j,1,1,element_size)=temp1(i,j)
   183: ||                  pcmor_nc1(j,i,1,1,element_size)=temp1(i,j)
   184: ||*---->            do ii=1,lx1
   185: |||       c...........p0 is for the case that a nonconforming edge is shared by
   186: |||       c           two conforming faces
   187: |||                   p0(ii,1,1)=p(ii,1,1)
   188: |||*--->              do jj=1,lx1
   189: ||||      c.............now p is for the case that a nonconforming edge is shared
   190: ||||      c             by nonconforming faces
   191: ||||                    p(ii,1,jj)=p(ii,jj,1)
   192: |||*---               end do
   193: ||*----             end do
   194: ||        
   195: ||                  call laplacian(temp,p,element_size)
   196: ||                  call transfb_nc2(temp1,temp)
   197: ||        
   198: ||        c.........pcmor_nc2(i,j,ii,jj,element_size) gives the intermediate
   199: ||        c         preconditioner value on collocation point (i,j) on a
   200: ||        c         nonconforming face of an element with size size_element
   201: ||        
   202: ||                  pcmor_nc2(i,j,1,1,element_size)=temp1(i,j)*2.d0
   203: ||                  pcmor_nc2(j,i,1,1,element_size)=
   204: ||             &          pcmor_nc2(i,j,1,1,element_size)
   205: ||        
   206: ||                  call laplacian(temp,p0,element_size)
   207: ||                  call transfb_nc0(temp1,temp)
   208: ||        
   209: ||        c.........pcmor_nc0(i,j,ii,jj,element_size) gives the intermediate
   210: ||        c         preconditioner value on collocation point (i,j) on a
   211: ||        c         conforming face of an element, which shares a nonconforming
   212: ||        c         edge with another conforming face
   213: ||                  pcmor_nc0(i,j,1,1,element_size)=temp1(i,j)
   214: ||                  pcmor_nc0(j,i,1,1,element_size)=temp1(i,j)
   215: |+-----           end do
   216: |                 pcmor_nc1(lx1,j,1,1,element_size)=
   217: |              &        pcmor_nc1(lx1,j,1,1,element_size)*2.d0
   218: |                 pcmor_nc1(j,lx1,1,1,element_size)=
   219: |              &        pcmor_nc1(lx1,j,1,1,element_size)
   220: |                 pcmor_nc2(lx1,j,1,1,element_size)=
   221: |              &        pcmor_nc2(lx1,j,1,1,element_size)*2.d0
   222: |                 pcmor_nc2(j,lx1,1,1,element_size)=
   223: |              &        pcmor_nc2(lx1,j,1,1,element_size)
   224: |                 pcmor_nc0(lx1,j,1,1,element_size)=
   225: |              &        pcmor_nc0(lx1,j,1,1,element_size)*2.d0
   226: |                 pcmor_nc0(j,lx1,1,1,element_size)=
   227: |              &        pcmor_nc0(lx1,j,1,1,element_size)
   228: |         
   229: |         c.......symmetrical copy
   230: |*----->          do i=1,lx1-1
   231: ||                  pcmor_nc1(i,j,1,2,element_size)=
   232: ||             &          pcmor_nc1(lx1+1-i,j,1,1,element_size)
   233: ||                  pcmor_nc0(i,j,1,2,element_size)=
   234: ||             &          pcmor_nc0(lx1+1-i,j,1,1,element_size)
   235: ||                  pcmor_nc2(i,j,1,2,element_size)=
   236: ||             &          pcmor_nc2(lx1+1-i,j,1,1,element_size)
   237: |*-----           end do
   238: |         
   239: |*----->          do j=2,lx1
   240: ||*---->            do i=1,lx1-1
   241: |||                   pcmor_nc1(i,j,1,2,element_size)=
   242: |||            &            pcmor_nc1(lx1+1-i,j,1,1,element_size)
   243: ||*----             end do
   244: ||                  i=lx1
   245: ||                  pcmor_nc1(i,j,1,2,element_size)=
   246: ||             &          pcmor_nc1(lx1+1-i,j,1,1,element_size)
   247: ||                  pcmor_nc0(i,j,1,2,element_size)=
   248: ||             &          pcmor_nc0(lx1+1-i,j,1,1,element_size)
   249: ||                  pcmor_nc2(i,j,1,2,element_size)=
   250: ||             &          pcmor_nc2(lx1+1-i,j,1,1,element_size)
   251: |*-----           end do
   252: |         
   253: |                 j=1
   254: |                 i=1
   255: |                 pcmor_nc1(i,j,2,1,element_size)=
   256: |              &        pcmor_nc1(i,lx1+1-j,1,1,element_size)
   257: |                 pcmor_nc0(i,j,2,1,element_size)=
   258: |              &        pcmor_nc0(i,lx1+1-j,1,1,element_size)
   259: |                 pcmor_nc2(i,j,2,1,element_size)=
   260: |              &        pcmor_nc2(i,lx1+1-j,1,1,element_size)
   261: |*----->          do j=2,lx1-1
   262: ||                  i=1
   263: ||                  pcmor_nc1(i,j,2,1,element_size)=
   264: ||             &          pcmor_nc1(i,lx1+1-j,1,1,element_size)
   265: ||                  pcmor_nc0(i,j,2,1,element_size)=
   266: ||             &          pcmor_nc0(i,lx1+1-j,1,1,element_size)
   267: ||                  pcmor_nc2(i,j,2,1,element_size)=
   268: ||             &          pcmor_nc2(i,lx1+1-j,1,1,element_size)
   269: ||*---->            do i=2,lx1
   270: |||                   pcmor_nc1(i,j,2,1,element_size)=
   271: |||            &            pcmor_nc1(i,lx1+1-j,1,1,element_size)
   272: ||*----             end do
   273: |*-----           end do
   274: |         
   275: |                 j=lx1
   276: |*----->          do i=2,lx1
   277: ||                  pcmor_nc1(i,j,2,1,element_size)=
   278: ||             &          pcmor_nc1(i,lx1+1-j,1,1,element_size)
   279: ||                  pcmor_nc0(i,j,2,1,element_size)=
   280: ||             &          pcmor_nc0(i,lx1+1-j,1,1,element_size)
   281: ||                  pcmor_nc2(i,j,2,1,element_size)=
   282: ||             &          pcmor_nc2(i,lx1+1-j,1,1,element_size)
   283: |*-----           end do
   284: |         
   285: |                 j=1
   286: |                 i=lx1
   287: |                 pcmor_nc1(i,j,2,2,element_size)=
   288: |              &        pcmor_nc1(lx1+1-i,lx1+1-j,1,1,element_size)
   289: |                 pcmor_nc0(i,j,2,2,element_size)=
   290: |              &        pcmor_nc0(lx1+1-i,lx1+1-j,1,1,element_size)
   291: |                 pcmor_nc2(i,j,2,2,element_size)=
   292: |              &        pcmor_nc2(lx1+1-i,lx1+1-j,1,1,element_size)
   293: |          
   294: |*----->          do j=2,lx1-1
   295: ||*---->            do i=2,lx1-1
   296: |||                   pcmor_nc1(i,j,2,2,element_size)=
   297: |||            &            pcmor_nc1(lx1+1-i,lx1+1-j,1,1,element_size)
   298: ||*----             end do
   299: ||                  i=lx1
   300: ||                  pcmor_nc1(i,j,2,2,element_size)=
   301: ||             &          pcmor_nc1(lx1+1-i,lx1+1-j,1,1,element_size)
   302: ||                  pcmor_nc0(i,j,2,2,element_size)=
   303: ||             &          pcmor_nc0(lx1+1-i,lx1+1-j,1,1,element_size)
   304: ||                  pcmor_nc2(i,j,2,2,element_size)=
   305: ||             &          pcmor_nc2(lx1+1-i,lx1+1-j,1,1,element_size)
   306: |*-----           end do
   307: |                 j=lx1
   308: |*----->          do i=2,lx1-1
   309: ||                  pcmor_nc1(i,j,2,2,element_size)=
   310: ||             &          pcmor_nc1(lx1+1-i,lx1+1-j,1,1,element_size)
   311: ||                  pcmor_nc0(i,j,2,2,element_size)=
   312: ||             &          pcmor_nc0(lx1+1-i,lx1+1-j,1,1,element_size)
   313: ||                  pcmor_nc2(i,j,2,2,element_size)=
   314: ||             &          pcmor_nc2(lx1+1-i,lx1+1-j,1,1,element_size)
   315: |*-----           end do
   316: |         
   317: |         
   318: |         c.......vertices shared by at least one nonconforming face or edge
   319: |         
   320: |         c.......Among three edges and three faces sharing a vertex on an element
   321: |         c       situation 1: only one edge is nonconforming
   322: |         c       situation 2: two edges are nonconforming
   323: |         c       situation 3: three edges are nonconforming
   324: |         c       situation 4: one face is nonconforming
   325: |         c       situation 5: one face and one edge are nonconforming
   326: |         c       situation 6: two faces are nonconforming
   327: |         c       situation 7: three faces are nonconforming
   328: |         
   329: |                 call r_init(p0,nxyz,0.d0)
   330: |                 p0(1,1,1)=1.d0
   331: |                 call laplacian(temp,p0,element_size)
   332: |                 pcmor_cor(8,element_size)=temp(1,1,1)
   333: |         
   334: |         c.......situation 1
   335: |                 call r_init(p0,nxyz,0.d0)
   336: |*----->          do i=1,lx1
   337: ||                   p0(i,1,1)=tcpre(i,1)
   338: |*-----           end do
   339: |                 call laplacian(temp,p0,element_size)
   340: |                 call transfb_cor_e(1,pcmor_cor(1,element_size),temp)
   341: |         
   342: |         c.......situation 2
   343: |                 call r_init(p0,nxyz,0.d0)
   344: |*----->          do i=1,lx1
   345: ||                   p0(i,1,1)=tcpre(i,1)
   346: ||                   p0(1,i,1)=tcpre(i,1)
   347: |*-----           end do
   348: |                 call laplacian(temp,p0,element_size)
   349: |                 call transfb_cor_e(2,pcmor_cor(2,element_size),temp)
   350: |         
   351: |         c.......situation 3
   352: |                 call r_init(p0,nxyz,0.d0)
   353: |*----->          do i=1,lx1
   354: ||                   p0(i,1,1)=tcpre(i,1)
   355: ||                   p0(1,i,1)=tcpre(i,1)
   356: ||                   p0(1,1,i)=tcpre(i,1)
   357: |*-----           end do
   358: |                 call laplacian(temp,p0,element_size)
   359: |                 call transfb_cor_e(3,pcmor_cor(3,element_size),temp)
   360: |         
   361: |         c.......situation 4
   362: |                 call r_init(p0,nxyz,0.d0)
   363: |*----->          do j=1,lx1
   364: ||*---->            do i=1,lx1
   365: |||                   p0(i,j,1)=tcpre(i,j)
   366: ||*----             end do
   367: |*-----           end do
   368: |                 call laplacian(temp,p0,element_size)
   369: |                 call transfb_cor_f(4,pcmor_cor(4,element_size),temp)
   370: |         
   371: |         c.......situation 5
   372: |                 call r_init(p0,nxyz,0.d0)
   373: |*----->          do j=1,lx1
   374: ||*---->            do i=1,lx1
   375: |||                   p0(i,j,1)=tcpre(i,j)
   376: ||*----             end do
   377: |*-----           end do
   378: |*----->          do i=1,lx1
   379: ||                   p0(1,1,i)=tcpre(i,1)
   380: |*-----           end do
   381: |                 call laplacian(temp,p0,element_size)
   382: |                 call transfb_cor_f(5,pcmor_cor(5,element_size),temp)
   383: |          
   384: |         c.......situation 6
   385: |                 call r_init(p0,nxyz,0.d0)
   386: |+----->          do j=1,lx1
   387: ||*---->            do i=1,lx1
   388: |||                   p0(i,j,1)=tcpre(i,j)
   389: |||                   p0(i,1,j)=tcpre(i,j)
   390: ||*----             end do
   391: |+-----           end do
   392: |                 call laplacian(temp,p0,element_size)
   393: |                 call transfb_cor_f(6,pcmor_cor(6,element_size),temp)
   394: |         
   395: |         c.......situation 7
   396: |+----->          do j=1,lx1
   397: ||*---->            do i=1,lx1
   398: |||                   p0(i,j,1)=tcpre(i,j)
   399: |||                   p0(i,1,j)=tcpre(i,j)
   400: |||                   p0(1,i,j)=tcpre(i,j)
   401: ||*----             end do
   402: |+-----           end do
   403: |                 call laplacian(temp,p0,element_size)
   404: |                 call transfb_cor_f(7,pcmor_cor(7,element_size),temp)
   405: |         
   406: P------         end do
   407:           c$OMP END PARALLEL DO
   408:                 return
   409:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETPCMO
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   425: par(1801): Parallel routine generated.: SETPCMO$1
   427: par(1803): Parallelized by "do".
   427: vec( 101): Vectorized loop.
   433: par(1803): Parallelized by "do".
   434: vec( 101): Vectorized loop.
   435: opt(1418): Constant-length loop is expanded.
   443: par(1801): Parallel routine generated.: SETPCMO$2
   446: par(1803): Parallelized by "do".
   482: vec( 103): Unvectorized loop.
   482: vec( 110): Vectorization obstructive procedure reference.: COM_DPC
   531: vec( 103): Unvectorized loop.
   531: vec( 107): Iteration count is too small.
   532: opt(1418): Constant-length loop is expanded.
   533: opt(1418): Constant-length loop is expanded.
   534: opt(1418): Constant-length loop is expanded.
   536: opt(1033): Potential multiple store conflict -- use directive if OK.: DPCMOR
   546: opt(1418): Constant-length loop is expanded.
   559: opt(1418): Constant-length loop is expanded.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: SETPCMO
FORMAT LIST

 LINE   LOOP      STATEMENT

   413:                 subroutine setpcmo
   414:           c------------------------------------------------------------------------
   415:           c     compute the preconditioner by identifying its geometry configuration
   416:           c     and sum the values from the precomputed elemental contributions
   417:           c------------------------------------------------------------------------
   418:            
   419:                 use ua_data
   420:                 implicit none
   421:           
   422:                 integer face2, nb1, nb2, sizei, imor, enum, i,j,
   423:                &        iel, iside, nn1, nn2
   424:           
   425:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IMOR,IEL,ISIDE,I)
   426:           c$OMP DO
   427: Y------>        do imor=1,nvertex
   428: |                ifpcmor(imor)=.false.
   429: Y------         end do
   430:           c$OMP END DO nowait
   431:            
   432:           c$OMP DO
   433: P------>        do iel=1,nelt
   434: |V----->          do iside=1,nsides
   435: ||*---->            do i=1,4
   436: |||                   edgevis(i,iside,iel)=.false.
   437: ||*----             end do
   438: |V-----           end do
   439: P------         end do
   440:           c$OMP END DO
   441:           c$OMP END PARALLEL
   442:           
   443:           c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IEL,iside,sizei,
   444:           c$OMP& imor,enum,face2,nb1,nb2,i,j,nn1,nn2)
   445:           
   446: P------>        do iel=1,nelt
   447: |+----->          do iside=1,nsides
   448: ||        c.........for nonconforming faces
   449: ||                  if(cbc(iside,iel).eq.3)then
   450: ||                    sizei=size_e(iel)
   451: ||        
   452: ||        c...........vertices
   453: ||        
   454: ||        c...........ifpcmor(imor)=.true. indicates that mortar point imor has
   455: ||        c           been visited
   456: ||                    imor=idmo(1,1,1,1,iside,iel)
   457: ||                    if(.not.ifpcmor(imor))then
   458: ||        c.............compute the preconditioner on mortar point imor
   459: ||                      call pc_corner(imor)
   460: ||                      ifpcmor(imor)=.true.
   461: ||                    end if
   462: ||        
   463: ||                    imor=idmo(lx1,1,1,2,iside,iel)
   464: ||                    if(.not.ifpcmor(imor))then
   465: ||                      call pc_corner(imor)
   466: ||                      ifpcmor(imor)=.true.
   467: ||                    end if
   468: ||        
   469: ||                    imor=idmo(1,lx1,2,1,iside,iel)
   470: ||                    if(.not.ifpcmor(imor))then
   471: ||                      call pc_corner(imor)
   472: ||                      ifpcmor(imor)=.true.
   473: ||                    end if
   474: ||        
   475: ||                    imor=idmo(lx1,lx1,2,2,iside,iel)
   476: ||                    if(.not.ifpcmor(imor))then
   477: ||                      call pc_corner(imor)
   478: ||                      ifpcmor(imor)=.true.
   479: ||                    end if
   480: ||        
   481: ||        c...........edges on nonconforming faces, enum is local edge number
   482: ||+---->              do enum=1,4
   483: |||       
   484: |||       c.............edgevis(enum,iside,iel)=.true. indicates that local edge
   485: |||       c             enum of face iside of iel has been visited
   486: |||                     if(.not.edgevis(enum,iside,iel))then
   487: |||                       edgevis(enum,iside,iel)=.true.
   488: |||       
   489: |||       c...............Examing neighbor element information,
   490: |||       c               calculateing the preconditioner value.
   491: |||                       face2= f_e_ef(enum,iside)
   492: |||                       if(cbc(face2,iel).eq.2)then
   493: |||                         nb1=sje(1,1,face2,iel)
   494: |||                         if(cbc(iside,nb1).eq.2)then
   495: |||       
   496: |||       c...................Compute the preconditioner on local edge enum on face
   497: |||       c                   iside of element iel, 1 is neighborhood information got
   498: |||       c                   by examing neighbors(nb1). For detailed meaning of 1,
   499: |||       c                   see subroutine com_dpc.
   500: |||       
   501: |||                           call com_dpc(iside,iel,enum,1,sizei)
   502: |||                           nb2=sje(1,1,iside,nb1)
   503: |||                           edgevis(op(e_face2(enum,iside)),
   504: |||            &                      jjface(face2),nb2)=.true.
   505: |||       
   506: |||                         elseif(cbc(iside,nb1).eq.3)then
   507: |||                           call com_dpc(iside,iel,enum,2,sizei)
   508: |||                           edgevis(op(enum),iside,nb1)=.true.
   509: |||                         end if
   510: |||       
   511: |||                       elseif(cbc(face2,iel).eq.3)then
   512: |||                         edgevis(e_face2(enum,iside),face2,iel)=.true.
   513: |||                         nb1=sje(1,2,face2,iel)
   514: |||                         if(cbc(iside,nb1).eq.1)then
   515: |||                           call com_dpc(iside,iel,enum,3,sizei)
   516: |||                           nb2=sje(1,1,iside,nb1)
   517: |||                           edgevis(op(enum),jjface(iside),nb2)=.true.
   518: |||                           edgevis(op(e_face2(enum,iside)),
   519: |||            &                      jjface(face2),nb2)=.true.
   520: |||                         elseif(cbc(iside,nb1).eq.2)then
   521: |||                           call com_dpc(iside,iel,enum,4,sizei)
   522: |||                         end if
   523: |||                       else if (cbc(face2,iel).eq.0)then
   524: |||                         call com_dpc(iside,iel,enum,0,sizei)
   525: |||                       end if
   526: |||                     end if
   527: ||+----               end do
   528: ||        
   529: ||        c...........mortar element interior (not edge of mortar)
   530: ||        
   531: ||+---->              do nn1=1,2
   532: |||*--->                do nn2=1,2
   533: ||||*-->                  do j=2,lx1-1
   534: |||||*->                    do i=2,lx1-1
   535: ||||||                        imor=idmo(i,j,nn1,nn2,iside,iel)
   536: ||||||                        dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,nn1,nn2,sizei)+
   537: ||||||         &                                pcmor_c(i,j,sizei+1))
   538: |||||*-                     end do
   539: ||||*--                   end do
   540: |||*---                 end do
   541: ||+----               end do
   542: ||        
   543: ||        c...........for i,j=lx1 there are duplicated mortar points, so
   544: ||        c           pcmor_c needs to be doubled or quadrupled
   545: ||                    i=lx1
   546: ||*---->              do j=2,lx1-1
   547: |||                     imor=idmo(i,j,1,1,iside,iel)
   548: |||                     dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,1,1,sizei)+
   549: |||            &                          pcmor_c(i,j,sizei+1)*2.d0)
   550: |||                     imor=idmo(i,j,2,1,iside,iel)
   551: |||                     dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,2,1,sizei)+
   552: |||            &                          pcmor_c(i,j,sizei+1)*2.d0)
   553: ||*----               end do
   554: ||        
   555: ||                    j=lx1
   556: ||                    imor=idmo(i,j,1,1,iside,iel)
   557: ||                    dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,1,1,sizei)+
   558: ||             &                        pcmor_c(i,j,sizei+1)*4.d0)
   559: ||*---->              do i=2,lx1-1
   560: |||                     imor=idmo(i,j,1,1,iside,iel)
   561: |||                     dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,1,1,sizei)+
   562: |||            &                          pcmor_c(i,j,sizei+1)*2.d0)
   563: |||                     imor=idmo(i,j,1,2,iside,iel)
   564: |||                     dpcmor(imor) = 1.d0/(pcmor_nc1(i,j,1,2,sizei)+
   565: |||            &                          pcmor_c(i,j,sizei+1)*2.d0)
   566: ||*----               end do
   567: ||        
   568: ||                  end if
   569: |+-----           end do
   570: P------         end do
   571:           c$OMP END PARALLEL DO
   572:           
   573:                 return
   574:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: PC_CORNER
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   591: vec( 102): Partially vectorized loop.
   597: opt(1418): Constant-length loop is expanded.
   607: opt(1418): Constant-length loop is expanded.
   617: opt(1019): Feedback of scalar value from one loop pass to another.: N
   617: vec( 121): Unvectorizable dependency.
   619: vec( 121): Unvectorizable dependency.
   621: vec( 121): Unvectorizable dependency.
   623: vec( 121): Unvectorizable dependency.
   627: vec( 121): Unvectorizable dependency.
   629: vec( 121): Unvectorizable dependency.
   632: vec( 121): Unvectorizable dependency.
   634: vec( 121): Unvectorizable dependency.
   639: vec( 126): Idiom detected.: SUM


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: PC_CORNER
FORMAT LIST

 LINE   LOOP      STATEMENT

   577:                 subroutine pc_corner(imor)
   578:           c------------------------------------------------------------------------
   579:           c     calculate preconditioner value for vertex with mortar index imor
   580:           c------------------------------------------------------------------------
   581:           
   582:                 use ua_data
   583:                 implicit none
   584:           
   585:                 double precision tmortemp
   586:                 integer imor, inemo,ie, sizei,cornernumber,
   587:                &        sface,sedge,iiface,iface,iiedge,iedge,n
   588:           
   589:                 tmortemp=0.d0
   590:           c.....loop over all elements sharing this vertex
   591: S------>        do inemo=1,nemo(imor)
   592: |                 ie=emo(1,inemo,imor)
   593: |       G         sizei=size_e(ie)
   594: |                 cornernumber=emo(2,inemo,imor)
   595: |                 sface=0
   596: |                 sedge=0
   597: |*----->          do iiface=1,3
   598: ||      G           iface=f_c(iiface,cornernumber)
   599: ||        c.........sface sums the number of nonconforming faces sharing this vertex on
   600: ||        c         one element
   601: ||      G           if(cbc(iface,ie).eq.3)then
   602: ||                    sface=sface+1
   603: ||                  end if
   604: |*-----           end do
   605: |         c.......sedge sums the number of nonconforming edges sharing this vertex on
   606: |         c       one element
   607: |*----->          do iiedge=1,3
   608: ||      G           iedge=e_c(iiedge,cornernumber)
   609: ||      G           if(ncon_edge(iedge,ie))sedge=sedge+1
   610: |*-----           end do
   611: |         
   612: |         c.......each n indicates how many nonconforming faces and nonconforming
   613: |         c       edges share this vertex on an element,
   614: |         
   615: |                 if(sface.eq.0)then
   616: |                   if(sedge.eq.0)then
   617: |                      n=8
   618: |                   elseif(sedge.eq.1)then
   619: |                      n=1
   620: |                   elseif(sedge.eq.2)then
   621: |                      n=2
   622: |                   elseif(sedge.eq.3)then
   623: |                      n=3
   624: |                   end if
   625: |                 elseif (sface.eq.1)then
   626: |                   if (sedge.eq.1)then
   627: |                    n=5
   628: |                   else
   629: |                    n=4
   630: |                   end if
   631: |                 else if (sface.eq.2)then
   632: |                    n=6
   633: |                 else if(sface.eq.3)then
   634: |                    n=7
   635: |                 end if
   636: |          
   637: |         c.......sum the intermediate pre-computed preconditioner values for
   638: |         c       all elements
   639: |       G         tmortemp=tmortemp+pcmor_cor(n,sizei)
   640: |         
   641: S------         end do
   642:           
   643:           c.....dpcmor(imor) is the value of the preconditioner on mortar point imor
   644:                 dpcmor(imor)=1.d0/tmortemp
   645:           
   646:                 return
   647:                 end


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: COM_DPC
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   751: vec( 101): Vectorized loop.
   753: vec( 128): Fused multiply-add operation applied.
   757: opt(1033): Potential multiple store conflict -- use directive if OK.: DPCMOR


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:02 2021
FILE NAME: precond.f

PROCEDURE NAME: COM_DPC
FORMAT LIST

 LINE   LOOP      STATEMENT

   650:                 subroutine com_dpc(iside,iel,enumber,n,isize)
   651:           c------------------------------------------------------------------------
   652:           c     Compute preconditioner for local edge enumber of face iside
   653:           c     on element iel.
   654:           c     isize is element size,
   655:           c     n is one of five different configurations
   656:           c     anc1, ac, anc2, anc0 are coefficients for different edges.
   657:           c     nc0 refers to nonconforming edge shared by two conforming faces
   658:           c     nc1 refers to nonconforming edge shared by one nonconforming face
   659:           c     nc2 refers to nonconforming edges shared by two nonconforming faces
   660:           c     c refers to conforming edge
   661:           c------------------------------------------------------------------------
   662:           
   663:                 use ua_data
   664:                 implicit none
   665:           
   666:                 integer n, isize,iside,iel, enumber, nn1start, nn1end, nn2start,
   667:                &        nn2end, jstart, jend, istart, iend, i, j, nn1, nn2, imor
   668:                 double precision anc1,ac,anc2,anc0,temp
   669:           
   670:           c.....different local edges have different loop ranges
   671:                 if(enumber.eq.1)then
   672:                   nn1start=1
   673:                   nn1end=1
   674:                   nn2start=1
   675:                   nn2end=2
   676:                   jstart=1
   677:                   jend=1
   678:                   istart=2
   679:                   iend=lx1-1
   680:                 elseif (enumber.eq.2) then
   681:                   nn1start=1
   682:                   nn1end=2
   683:                   nn2start=2
   684:                   nn2end=2
   685:                   jstart=2
   686:                   jend=lx1-1
   687:                   istart=lx1
   688:                   iend=lx1
   689:                 elseif (enumber.eq.3) then
   690:                   nn1start=2
   691:                   nn1end=2
   692:                   nn2start=1
   693:                   nn2end=2
   694:                   jstart=lx1
   695:                   jend=lx1
   696:                   istart=2
   697:                   iend=lx1-1
   698:                 elseif (enumber.eq.4) then
   699:                   nn1start=1
   700:                   nn1end=2
   701:                   nn2start=1
   702:                   nn2end=1
   703:                   jstart=2
   704:                   jend=lx1-1
   705:                   istart=1
   706:                   iend=1
   707:                 end if
   708:           
   709:           c.....among the four elements sharing this edge
   710:           
   711:           c.....one has a smaller size
   712:                 if(n.eq.1)then
   713:                   anc1=2.d0
   714:                   ac=1.d0
   715:                   anc0=1.d0
   716:                   anc2=0.d0
   717:           
   718:           c.....two (neighbored by a face) are of  smaller size
   719:                 else if (n.eq.2)then
   720:                   anc1=2.d0
   721:                   ac=2.d0
   722:                   anc0=0.d0
   723:                   anc2=0.d0
   724:           
   725:           c.....two (neighbored by an edge) are of smaller size
   726:                 else if (n.eq.3)then
   727:                   anc2=2.d0
   728:                   ac=2.d0
   729:                   anc1=0.d0
   730:                   anc0=0.d0
   731:           
   732:           c.....three are of smaller size
   733:                 else if (n.eq.4)then
   734:                   anc1=0.d0
   735:                   ac=3.d0
   736:                   anc2=1.d0
   737:                   anc0=0.d0
   738:           
   739:           c.....on the boundary
   740:                 else if (n.eq.0)then
   741:                   anc1=1.d0
   742:                   ac=1.d0
   743:                   anc2=0.d0
   744:                   anc0=0.d0
   745:                 end if
   746:           
   747:           c.....edge interior
   748: +------>        do nn2=nn2start,nn2end
   749: |+----->          do nn1=nn1start,nn1end
   750: ||+---->            do j=jstart,jend
   751: |||V--->              do i=istart,iend
   752: ||||                    imor=idmo(i,j,nn1,nn2,iside,iel)
   753: ||||    F               temp=anc1* pcmor_nc1(i,j,nn1,nn2,isize) +
   754: ||||           &             ac*  pcmor_c(i,j,isize+1)+
   755: ||||           &             anc0*  pcmor_nc0(i,j,nn1,nn2,isize)+
   756: ||||           &             anc2*pcmor_nc2(i,j,nn1,nn2,isize)
   757: ||||    C                 dpcmor(imor)=1.d0/temp
   758: |||V---                 end do
   759: ||+----               end do
   760: |+-----             end do
   761: +------           end do
   762:           
   763:           c.......local edge 1
   764:                   if (enumber.eq.1) then
   765:                     imor=idmo(lx1,1,1,1,iside,iel)
   766:                     temp=anc1* pcmor_nc1(lx1,1,1,1,isize) +
   767:                &         ac*  pcmor_c(lx1,1,isize+1)*2.d0+
   768:                &         anc0*  pcmor_nc0(lx1,1,1,1,isize)+
   769:                &         anc2*pcmor_nc2(lx1,1,1,1,isize)
   770:           c.......local edge 2
   771:                   elseif (enumber.eq.2) then
   772:                     imor=idmo(lx1,lx1,1,2,iside,iel)
   773:                     temp=anc1* pcmor_nc1(lx1,lx1,1,2,isize) +
   774:                &         ac*  pcmor_c(lx1,lx1,isize+1)*2.d0+
   775:                &         anc0*  pcmor_nc0(lx1,lx1,1,2,isize)+
   776:                &         anc2*pcmor_nc2(lx1,lx1,1,2,isize)
   777:           c.......local edge 3
   778:                   elseif (enumber.eq.3) then
   779:                     imor=idmo(lx1,lx1,2,1,iside,iel)
   780:                     temp=anc1* pcmor_nc1(lx1,lx1,2,1,isize) +
   781:                &         ac*  pcmor_c(lx1,lx1,isize+1)*2.d0+
   782:                &         anc0*  pcmor_nc0(lx1,lx1,2,1,isize)+
   783:                &         anc2*pcmor_nc2(lx1,lx1,2,1,isize)
   784:           c.......local edge 4
   785:                   elseif (enumber.eq.4) then
   786:                     imor=idmo(1,lx1,1,1,iside,iel)
   787:                     temp=anc1* pcmor_nc1(1,lx1,1,1,isize) +
   788:                &         ac*  pcmor_c(1,lx1,isize+1)*2.d0+
   789:                &         anc0*  pcmor_nc0(1,lx1,1,1,isize)+
   790:                &         anc2*pcmor_nc2(1,lx1,1,1,isize)
   791:                   end if
   792:           
   793:                   dpcmor(imor)=1.d0/temp
   794:           
   795:                 return
   796:                 end


