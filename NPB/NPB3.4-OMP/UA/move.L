NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    18: par(1801): Parallel routine generated.: MOVE$1
    21: vec( 110): Vectorization obstructive procedure reference.: COPY
    21: par(1803): Parallelized by "do".
    28: vec( 101): Vectorized loop.
    49: opt(1418): Constant-length loop is expanded.
    50: opt(1418): Constant-length loop is expanded.
    73: vec( 110): Vectorization obstructive procedure reference.: NCOPY
    73: vec( 110): Vectorization obstructive procedure reference.: COPY
    73: par(1803): Parallelized by "do".
    73: vec( 103): Unvectorized loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:51:00 2021
FILE NAME: move.f

PROCEDURE NAME: MOVE
FORMAT LIST

 LINE   LOOP      STATEMENT

     2:                 subroutine move
     3:           c---------------------------------------------------------------
     4:           c     move element to proper location in morton space filling curve
     5:           c---------------------------------------------------------------
     6:           
     7:                 use ua_data
     8:                 implicit none
     9:           
    10:                 integer i,iside,jface,iel,ntemp,ii1,ii2,n1,n2,cb
    11:           
    12:                 n2=2*6*nelt
    13:                 n1=n2*2
    14:           
    15:                 call nr_init_omp(sje_new,n1,0)
    16:                 call nr_init_omp(ijel_new,n2,0)
    17:           
    18:           c$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iel,i,iside,jface,cb,ntemp,
    19:           c$OMP& ii1,ii2)
    20:           c$OMP DO
    21: P------>        do iel=1,nelt
    22: |                 i=mt_to_id(iel)
    23: |                 treenew(iel)=tree(i)
    24: |                 call copy(xc_new(1,iel),xc(1,i),8)
    25: |                 call copy(yc_new(1,iel),yc(1,i),8)
    26: |                 call copy(zc_new(1,iel),zc(1,i),8)
    27: |         
    28: |V----->          do iside=1,nsides
    29: ||                  jface = jjface(iside)
    30: ||                  cb=cbc(iside,i)
    31: ||                  xc_new(iside,iel)=xc(iside,i)
    32: ||                  yc_new(iside,iel)=yc(iside,i)
    33: ||                  zc_new(iside,iel)=zc(iside,i)
    34: ||                  cbc_new(iside,iel)=cb
    35: ||        
    36: ||                  if(cb.eq.2)then
    37: ||                    ntemp=sje(1,1,iside,i)
    38: ||                    ijel_new(1,iside,iel)=1
    39: ||                    ijel_new(2,iside,iel)=1
    40: ||      G             sje_new(1,1,iside,iel)=id_to_mt(ntemp)
    41: ||        
    42: ||                  else if(cb.eq.1) then
    43: ||                    ntemp=sje(1,1,iside,i)
    44: ||                    ijel_new(1,iside,iel)=ijel(1,iside,i)
    45: ||                    ijel_new(2,iside,iel)=ijel(2,iside,i)
    46: ||      G             sje_new(1,1,iside,iel)=id_to_mt(ntemp)
    47: ||         
    48: ||                  else if(cb.eq.3) then
    49: ||*---->              do ii2=1,2
    50: |||*--->                do ii1=1,2
    51: ||||                      ntemp=sje(ii1,ii2,iside,i)
    52: ||||                      ijel_new(1,iside,iel)=1
    53: ||||                      ijel_new(2,iside,iel)=1
    54: ||||    G                 sje_new(ii1,ii2,iside,iel)=id_to_mt(ntemp)
    55: |||*---                 end do
    56: ||*----               end do
    57: ||        
    58: ||                  else if(cb.eq.0)then
    59: ||                    sje_new(1,1,iside,iel)=0
    60: ||                    sje_new(1,2,iside,iel)=0
    61: ||                    sje_new(2,1,iside,iel)=0
    62: ||                    sje_new(2,2,iside,iel)=0
    63: ||                  end if
    64: ||        
    65: |V-----           end do
    66: |         
    67: |                 call copy(ta2(1,1,1,iel),ta1(1,1,1,i),nxyz)
    68: |         
    69: P------         end do
    70:           c$OMP ENDDO
    71:           
    72:           c$OMP DO
    73: P------>        do iel=1,nelt
    74: |                 call copy(xc(1,iel),xc_new(1,iel),8)
    75: |                 call copy(yc(1,iel),yc_new(1,iel),8)
    76: |                 call copy(zc(1,iel),zc_new(1,iel),8)
    77: |                 call copy(ta1(1,1,1,iel),ta2(1,1,1,iel),nxyz)
    78: |                 call ncopy(sje(1,1,1,iel),sje_new(1,1,1,iel),4*6)
    79: |                 call ncopy(ijel(1,1,iel),ijel_new(1,1,iel),2*6)
    80: |                 call ncopy(cbc(1,iel),cbc_new(1,iel),6)
    81: |                 mt_to_id(iel)=iel
    82: |                 id_to_mt(iel)=iel
    83: |                 tree(iel)=treenew(iel)
    84: P------         end do
    85:           c$OMP ENDDO
    86:           c$OMP END PARALLEL
    87:           
    88:                 return
    89:                 end


