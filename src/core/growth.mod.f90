!    EmbryoMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.



!TO DO: 

!-reallocate the matrix node and all the variables that have a nd or nda dimension

module growth  !MIQUEL MADE MODULE
use general
use neighboring
use io
!use mitosi	!>>>>>>>>Miquel 21-3-13
!use del         !>>>>>>>>>>Roland
use aleas
use shell ! miguel4-11-13
use energy ! >>> Is 5-6-14

real*8,parameter::third=1d0/3d0
integer aax,aaw,aaz

contains

!******************************************************************

subroutine should_I_grow !>>>>>>>> Miquel21-10-13   only one node per cell grows at a time. The amount of %eqd gained by this node depends on the expression genes among the entire cell
integer::i,j,k,ii,jj,kk,iii,jjj,kkk,ci
real*8::grate,grac
integer nnod
real*8 a,b,c,d,aa,bb,ra,rb,rc,ja,jb,jc,jd,jaa
real*8 maxreq,differ

el: do ci=1,ncels

      ! we first check if the nodes of the cell are too compressed
      a=0.0d0
      do i=1,cels(ci)%nunodes
        ii=cels(ci)%node(i)
        a=a+fmeanl(ii)
      end do
      b=cels(ci)%nunodes
      a=a/b
!print*,ci,ii,"a",a,"mincomp",min_comp
      if (a<-epsilod.and.abs(a)>abs(min_comp)) cycle ! if the compresion per node in the cell is too large (in negative) then we do not add nodes

      if(node(cels(ci)%node(1))%fix>0) cycle !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
      grac=0d0
      if(cels(ci)%ctipus<3)then  !>>Miquel23-1-14
        maxreq=2*mmae !
        nnod=cels(ci)%nunodes      !
        iii=0 !the node that doesn't have max size
        grate=0d0 ; grac=0d0
        do ii=1,nnod
          j=cels(ci)%node(ii)
          differ=1-node(j)%dif
          if(node(j)%tipus==2) cycle   !>>Miquel23-1-14
          jjj=node(j)%altre             !>>Miquel23-1-14
          do jj=1,npag(nparam_per_node+1)    !number of genes affecting growth
            k=whonpag(nparam_per_node+1,jj)  !which are those genes
            grate=grate+(gex(j,k)+gex(jjj,k))*gen(k)%e(nparam_per_node+1)*differ !this is the differential of growth for the node !>>Miquel23-1-13
          end do                                                                  ! wa in space-req units
          grate=grate*delta
          a=node(j)%grd+node(jjj)%grd !;print*,"req",node(j)%eqd !>>Miquel28-1-14
          if(a<maxreq) iii=j  !this is the node that has to grow
        end do
      else                       !
        maxreq=mmae   !
        nnod=cels(ci)%nunodes      !
        iii=0 !the node that doesn't have max size
        grate=0d0 ; grac=0d0
        do ii=1,nnod
          j=cels(ci)%node(ii)
          differ=1-node(j)%dif
          do jj=1,npag(nparam_per_node+1)    !number of genes affecting growth
            k=whonpag(nparam_per_node+1,jj)  !which are those genes
            grate=grate+gex(j,k)*gen(k)%e(nparam_per_node+1)*differ !this is the differential of growth for the node !>>Miquel23-1-13
          end do
          grate=grate*delta
          a=node(j)%grd !;print*,"req",node(j)%eqd !>>Miquel23-1-14
          if(a<maxreq) iii=j  !this is the node that has to grow
        end do
      end if                     


      if(grate==0) cycle !there is no growth

      if(iii==0) then;grac=grate;goto 123;end if !all nodes are full, goes directly to add a new node

      if(node(iii)%tipus<3) grate=0.5d0*grate  !since there are two nodes growing, the rate per node is half

      a=node(iii)%grd ; !b=node(iii)%codel-a !the node increases in size
      !c=node(iii)%add-a
      aa=a+grate          !this is a growth in node's "radius"
      node(iii)%grd=aa
      !node(iii)%add=aa+c  !%add is also growing in the same amount (but not in the same proportion) than req

      
      nodeo(iii)%grd=nodeo(iii)%grd+grate   !?????????
!when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
      !nodeo(iii)%add=nodeo(iii)%add+grate    !values have to be updated, otherwise it will create a conflict with nexus

      if (node(iii)%tipus<3) then  !epithelial: the node beneath also grows in the same amount
        jj=node(iii)%altre
        ja=node(jj)%grd ; !jb=node(jj)%codel-a
        !jc=node(jj)%add-a
        jaa=ja+grate
        node(jj)%grd=jaa
        !node(jj)%add=jaa+jc

        nodeo(jj)%grd=nodeo(jj)%grd+grate   
!???????????????????????????????????????
!when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
        !nodeo(jj)%add=nodeo(jj)%add+grate    !values have to be updated, otherwise it will create a conflict with nexus
      end if

      if(aa>maxreq)then   !if the node has reax max size, we add a new node
        d=aa-maxreq
        node(iii)%grd=maxreq
        !node(iii)%add=maxreq+c
        grac=grac+d   !if the node is "full" we save the growth to apply it to another node

        nodeo(iii)%grd=maxreq    !???????????????????????????????????????
!when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
        !nodeo(iii)%add=maxreq+c    !values have to be updated, otherwise it will create a conflict with nexus

        if (node(iii)%tipus<3) then  !epithelial: the node below
          jd=jaa-maxreq
          node(jj)%grd=maxreq
          !node(jj)%add=maxreq+jc
          grac=grac+jd   !if the node is "full" we save the growth to apply it to another node

          nodeo(jj)%eqd=maxreq   !?????????????????????????????????????
!when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
          nodeo(jj)%add=maxreq+jc  !?????????????????????????????????????????????????????????  
!values have to be updated, otherwise it will create a conflict with nexus
        end if


        !Adding a new node
123     if (ffu(1)==1) then !ACHTUNG
        if(cels(ci)%ctipus<3)then             ! We add at most a node per cell per iteration
          !epithelial
          a=0 ; ii=0
          do i=1,cels(ci)%nunodes    !>> Miquel6-8-14
            j=cels(ci)%node(i)
            do jj=1,npag(nparam_per_node+9)
              k=whonpag(nparam_per_node+9,jj)   !IT'S IMPORTANT THAT THIS GENE DIFFUSES WITHIN THE CELL
              if (gex(j,k)>0) then      !POLAR GROWTH FORBIDS NON_POLAR GROWTH (makes the algorithm not to work) 
               a=a+gen(k)%e(nparam_per_node+9)*gex(j,k) ; ii=ii+1 !>> Miquel6-8-14
              end if
            end do
          end do
          a=a/real(ii)!>> Miquel6-8-14
          call random_number(ga)                      ! wa in units of probability
          if(ga<a)then !the intensity of the parameter determines the likelyhood of polar growth at a given time  !>>>>Miquel14-12-13
            call creixement_polar(ci)             ! we add a node (it does not matter which since the subroutine choses)
            goto 13                     
          else
            goto 12
          end if
12        call creixement(ci)   ! we add a node (it does not matter which since the subroutine choses)
13        if(grac<2*reqmin) then
            a=reqmin
            grac=0.0d0
          else
            if (grac>2*maxreq) then ; a=maxreq ; else ; a=grac*0.5d0 ; end if
          end if
          !c=node(cels(ci)%node(1))%add-node(cels(ci)%node(1))%eqd
          node(nd)%grd=a ; nodeo(nd)%grd=a
          !node(nd)%add=a+c   
          k=node(nd)%altre
          node(k)%grd=a ; nodeo(k)%grd=a
          !node(k)%add=a+c   
          grac=grac-a
        else
          !mesenchymal
          call random_number(ga)
          jjj=int(ga*nnod)+1
          a=0 ; ii=0
          do i=1,cels(ci)%nunodes    !>> Miquel6-8-14
            j=cels(ci)%node(i)
            do jj=1,npag(nparam_per_node+9)
              k=whonpag(nparam_per_node+9,jj)   !IT'S IMPORTANT THAT THIS GENE DIFFUSES WITHIN THE CELL
              if (gex(j,k)>0) then      !POLAR GROWTH FORBIDS NON_POLAR GROWTH (makes the algorithm not to work) 
               a=a+gen(k)%e(nparam_per_node+9); ii=ii+1 !>> Miquel6-8-14
              end if
            end do
          end do
          a=a/real(ii)!>> Miquel6-8-14
          call random_number(ga)
          if(ga<a)then !the intensity of the parameter determines the likelyhood of polar growth at a given time !>>>>Miquel14-12-13
!            print*,"CREIXEMENT POLAR"
            call creixement_polar(ci)   ! we add a node (it does not matter which since the subroutine choses)
            goto 14
          else
!            print*,"CREIXEMENT normal"

            goto 15
          end if
15        call creixement(ci)   ! we add a node (it does not matter which since the subroutine choses)
!14        c=node(cels(ci)%node(1))%add-node(cels(ci)%node(1))%eqd
14        if(grac<reqmin) then !;print*,"grac pre",grac
            a=reqmin
            grac=0.0d0
          else
            if (grac>maxreq) then ; a=maxreq ; else ; a=grac ; end if
          end if
          node(nd)%grd=a ; nodeo(nd)%grd=a !??????????????????????????????!;print*,"grac",grac,"a",a
          !node(nd)%add=a+c
          grac=grac-a
        end if
        if (ffu(13)==1.and.grate>0.and.grate>reqmin) then

          ! we first check if the nodes of the cell are too compressed
          a=0.0d0
          do i=1,cels(ci)%nunodes
            ii=cels(ci)%node(i)
            a=a+fmeanl(ii)
          end do
          b=cels(ci)%nunodes
          a=a/b
          if (a<min_comp) cycle ! if the compresion per node in the cell is too large (in negative) then we do not add nodes

          print *,"WARNING: we are adding more than one node per cell per iteration"
          goto 123  ! >>> Is 1-7-14 this makes that one can add more than one per cell and iteration
        end if      ! THIS IS NECESSARY OTHERWISE THE BEHAVIOUR OF THE MODEL DEPENDS ON WHETHER WE USE DYNAMIC DELTA OR NOT
end if
      end if
    end do el 
end subroutine


!************************************************************************************************************
! Is >>> 17-6-14
subroutine creixement(ce)
  integer nods,ce

  if (cels(ce)%ctipus<3) then
    ndepi=ndepi+2
    nd=nd+2
    nods=cels(ce)%node(1)
if (ffu(1)==1)  call addanode(nods) !ACHTUNG

    if (cels(ce)%nunodes==4) then
      call findepixyz_single(ce)
    else
      if (cels(ce)%nunodes==6) then
        call findepixyz_double(ce)
      else
        call findepixyz_nodel(ce)   !!>> HC 10-8-2021 redundant code removed here
      end if
    end if
  else
    ndmes=ndmes+1
    nd=nd+1
    nods=cels(ce)%node(1)
if (ffu(1)==1)  call addanode(nods) !ACHTUNG
!    call addanode(nods)
    if (cels(ce)%nunodes==2) then
      call findmesxyz_single(ce)
    else
      call findepixyz_nodel(ce)     !!>> HC 10-8-2021 redundant code removed here
    end if
  end if
end subroutine

!************************************************************************************************************
! Is >>> 25-6-14
subroutine creixement_polar(ce)
  integer nods,ce

  if (cels(ce)%ctipus<3) then
    ndepi=ndepi+2
    nd=nd+2

    if (cels(ce)%nunodes==4) then
      call findepixyz_single_polar(ce)

    else
      if (cels(ce)%nunodes==6) then
        call findepixyz_double_polar(ce)
      else
        call findepixyz_polar(ce)
      end if
    end if
  else
    ndmes=ndmes+1
    nd=nd+1
!    nods=cels(ce)%node(1)
!    call addanode(nods)
    if (cels(ce)%nunodes==2) then
      call findmesxyz_single_polar(ce)
    else
      call findepixyz_polar(ce)
    end if
  end if
  nods=cels(ce)%node(1)
if (ffu(1)==1)  call addanode(nods) !ACHTUNG
!  call addanode(nods)
end subroutine


!************************************************************************************************************

subroutine creixement_polar_old(ce)  !AQUI CAL MIRAR QUE PASA QUAN NUNODES==1 !!!!!!!!!!!!!!!!!!!!
  integer nods,ce

  if (cels(ce)%ctipus<3) then
    ndepi=ndepi+2
    nd=nd+2
    call findepixyz_polar(ce)
  else
    ndmes=ndmes+1
    nd=nd+1
    call findepixyz_polar(ce)
  end if
  nods=cels(ce)%node(1)
  call addanode(nods)
end subroutine

!*******************************************

subroutine findmesxyz_single(ce) ! >>> Is 17-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,ii

  ii=cels(ce)%node(1)

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  node(nd)%x=node(ii)%x+particions_esfera(k,1)*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd)%y=node(ii)%y+particions_esfera(k,2)*1d-4
  node(nd)%z=node(ii)%z+particions_esfera(k,3)*1d-4 
end subroutine

!*******************************************

subroutine findmesxyz_single_polar(ce) ! >>> Is 25-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,ii

  ii=cels(ce)%node(1) ;!print*,"cel",ce,"ii",ii

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  node(nd)%x=node(ii)%x+cels(ce)%polx*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd)%y=node(ii)%y+cels(ce)%poly*1d-4
  node(nd)%z=node(ii)%z+cels(ce)%polz*1d-4 
end subroutine


!*******************************************

subroutine findepixyz_single(ce)  ! >>> Is 17-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,i,ii
 
  ii=cels(ce)%node(1)
  if (node(ii)%tipus==1) then
    i=cels(ce)%node(2)
  else
    ii=cels(ce)%node(2)
    i=cels(ce)%node(1)
  end if

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  node(nd-1)%x=node(ii)%x+particions_esfera(k,1)*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd-1)%y=node(ii)%y+particions_esfera(k,2)*1d-4
  node(nd-1)%z=node(ii)%z+particions_esfera(k,3)*1d-4 
  node(nd)%x=node(i)%x+particions_esfera(k,1)*1d-4  
  node(nd)%y=node(i)%y+particions_esfera(k,2)*1d-4
  node(nd)%z=node(i)%z+particions_esfera(k,3)*1d-4 

end subroutine

!*******************************************

subroutine findepixyz_single_polar(ce)  ! >>> Is 25-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,i,ii
 
  ii=cels(ce)%node(1)
  if (node(ii)%tipus==1) then
    i=cels(ce)%node(2)
  else
    ii=cels(ce)%node(2)
    i=cels(ce)%node(1)
  end if

  ! we simply add the node in the direction of the polarity

  node(nd-1)%x=node(ii)%x+cels(ce)%polx*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd-1)%y=node(ii)%y+cels(ce)%poly*1d-4
  node(nd-1)%z=node(ii)%z+cels(ce)%polz*1d-4 
  node(nd)%x=node(i)%x+cels(ce)%polx*1d-4  
  node(nd)%y=node(i)%y+cels(ce)%poly*1d-4
  node(nd)%z=node(i)%z+cels(ce)%polz*1d-4 
end subroutine


!*******************************************

subroutine findepixyz_double(ce) ! >>> Is 17-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,kk,i,ii

11 call random_number(a) ! we chose one of the ellipses at random >>> Is 24-6-14
  kk=int(a*4)+1

  ii=cels(ce)%node(kk)
  if (node(ii)%tipus==1) then
    i=node(cels(ce)%node(kk))%altre
  else
    ii=node(cels(ce)%node(kk))%altre
    i=cels(ce)%node(kk)
  end if

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  node(nd-1)%x=node(ii)%x+particions_esfera(k,1)*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd-1)%y=node(ii)%y+particions_esfera(k,2)*1d-4
  node(nd-1)%z=node(ii)%z+particions_esfera(k,3)*1d-4 
  node(nd)%x=node(i)%x+particions_esfera(k,1)*1d-4  
  node(nd)%y=node(i)%y+particions_esfera(k,2)*1d-4
  node(nd)%z=node(i)%z+particions_esfera(k,3)*1d-4 
end subroutine

!*******************************************

subroutine findepixyz_double_polar(ce) ! >>> Is 25-6-14
  integer ce
  real*8 a,aaa,bbb,ccc
  integer k,kk,i,ii

11 call random_number(a) ! we chose one of the ellipses at random >>> Is 24-6-14
  kk=int(a*4)+1

  ii=cels(ce)%node(kk)
  if (node(ii)%tipus==1) then
    i=node(cels(ce)%node(kk))%altre
  else
    ii=node(cels(ce)%node(kk))%altre
    i=cels(ce)%node(kk)
  end if

  node(nd-1)%x=node(ii)%x+cels(ce)%polx*1d-4  ! WARNING ACHTUNG this is a very small arbitrary value but it may still give problems if cells are very dense in nodes
  node(nd-1)%y=node(ii)%y+cels(ce)%poly*1d-4
  node(nd-1)%z=node(ii)%z+cels(ce)%polz*1d-4 
  node(nd)%x=node(i)%x+cels(ce)%polx*1d-4  
  node(nd)%y=node(i)%y+cels(ce)%poly*1d-4
  node(nd)%z=node(i)%z+cels(ce)%polz*1d-4 
end subroutine


!************************************************************************************************************ 

!subroutine findepixyz(ce)  !it uses delaunay
!  integer i,ce
!  integer tipi,iik
!  real*8 aream
!  integer::j,ic,icc,ii,jj,kk,kkk,kkkk,kkkkk,um
!  real*8 ::i1x,i1y,i2x,i2y,i2z,i3x,i3y,i3z,i4x,i4y,i4z,j1x,j1y,area
!  
!  um=6
!  if (cels(ce)%ctipus==3) um=3
!  if (cels(ce)%nunodes<um.or.cels(ce)%ctipus==3) then
!    call findepixyz_nodel(ce)
!    return
!  end if
!
!!  iik=node(i)%icel
!  tipi=cels(ce)%ctipus
!  call delaunay(ce)  !triangulation algorithm
!
!  !so, we have to calculate the areas of all triangles
!  aream=0
!  aax=1; aaw=1; aaz=1
!  !iii=numnod  ! >>> Is 7-6-14
!  numnod=cels(ce)%nunodes  !>>> Is 4-5-14
!  do iii=1,numnod !trmax              !THIS IS KIND OF CRAP CODE
!    do jjj=1,numnod !trmax                                          ! miguel 14-5-13
!      if(triang(iii,1)==nlist(jjj,4))then                           ! miguel 14-5-13
!        i2x=nlist(jjj,1); i2y=nlist(jjj,2); !i2z=nlist(jjj,3)       ! miguel 14-5-13 
!      end if                                                        ! miguel 14-5-13
!      if(triang(iii,2)==nlist(jjj,4))then                           ! miguel 14-5-13
!        i3x=nlist(jjj,1); i3y=nlist(jjj,2); !i3z=nlist(jjj,3)       ! miguel 14-5-13 
!      end if                                                        ! miguel 14-5-13
!      if(triang(iii,3)==nlist(jjj,4))then                           ! miguel 14-5-13
!        i4x=nlist(jjj,1); i4y=nlist(jjj,2); !i4z=nlist(jjj,3)       ! miguel 14-5-13
!      end if                                                        ! miguel 14-5-13
!    end do                                                          ! miguel 14-5-13
!    i1x=i2x-i3x; i1y=i2y-i3y; j1x=i2x-i4x; j1y=i2y-i4y
!    area=0.5d0*abs(((i1x*j1y)-(i1y*j1x)))
!    if(area.gt.aream) then
!      aream=area
!      aax=iii
!      aax=triang(iii,1) ; aaw=triang(iii,2) ; aaz=triang(iii,3) ! miguel 14-5-13 substitution of tri(:,:) by triang(:,:)
!    end if
!  end do
!
!  !THE NEW NODES are now put un the right triangle
!    node(nd)%x=(node(aax)%x+third*((node(aaz)%x-node(aax)%x)+(node(aaw)%x-node(aax)%x)))
!    node(nd)%y=(node(aax)%y+third*((node(aaz)%y-node(aax)%y)+(node(aaw)%y-node(aax)%y)))
!    node(nd)%z=(node(aax)%z+third*((node(aaz)%z-node(aax)%z)+(node(aaw)%z-node(aax)%z)))
!    
!    iii=node(aax)%altre;jjj=node(aaz)%altre;kkk=node(aaw)%altre
!    node(nd-1)%x=(node(iii)%x+third*((node(jjj)%x-node(iii)%x)+(node(kkk)%x-node(iii)%x)))
!    node(nd-1)%y=(node(iii)%y+third*((node(jjj)%y-node(iii)%y)+(node(kkk)%y-node(iii)%y)))
!    node(nd-1)%z=(node(iii)%z+third*((node(jjj)%z-node(iii)%z)+(node(kkk)%z-node(iii)%z)))
!    
!
!!  if(tipi==2)then                          !X>>>>>>> Miquel21-10-13
!!    node(nd)%x=(node(aax)%x+third*((node(aaz)%x-node(aax)%x)+(node(aaw)%x-node(aax)%x)))
!!    node(nd)%y=(node(aax)%y+third*((node(aaz)%y-node(aax)%y)+(node(aaw)%y-node(aax)%y)))
!!    node(nd)%z=(node(aax)%z+third*((node(aaz)%z-node(aax)%z)+(node(aaw)%z-node(aax)%z)))
!!  else
!!    iii=node(aax)%altre;jjj=node(aaz)%altre;kkk=node(aaw)%altre
!!    node(nd)%x=(node(iii)%x+third*((node(jjj)%x-node(iii)%x)+(node(kkk)%x-node(iii)%x)))
!!    node(nd)%y=(node(iii)%y+third*((node(jjj)%y-node(iii)%y)+(node(kkk)%y-node(iii)%y)))
!!    node(nd)%z=(node(iii)%z+third*((node(jjj)%z-node(iii)%z)+(node(kkk)%z-node(iii)%z)))
!!  end if
!
!!  if(tipi==1)then
!!    node(nd-1)%x=(node(aax)%x+third*((node(aaz)%x-node(aax)%x)+(node(aaw)%x-node(aax)%x)))
!!    node(nd-1)%y=(node(aax)%y+third*((node(aaz)%y-node(aax)%y)+(node(aaw)%y-node(aax)%y)))
!!    node(nd-1)%z=(node(aax)%z+third*((node(aaz)%z-node(aax)%z)+(node(aaw)%z-node(aax)%z)))
!!  else
!!    iii=node(aax)%altre;jjj=node(aaz)%altre;kkk=node(aaw)%altre
!!    node(nd-1)%x=(node(iii)%x+third*((node(jjj)%x-node(iii)%x)+(node(kkk)%x-node(iii)%x)))
!!    node(nd-1)%y=(node(iii)%y+third*((node(jjj)%y-node(iii)%y)+(node(kkk)%y-node(iii)%y)))
!!    node(nd-1)%z=(node(iii)%z+third*((node(jjj)%z-node(iii)%z)+(node(kkk)%z-node(iii)%z)))
!end subroutine

!***********************************************************

subroutine findepixyz_polar(iik)  !AQUI CAL MIRAR QUE PASA QUAN NUNODES==1 !!!!!!!!!!!!!!!!!!!!
  integer no
  integer iik,tipi,nvaloq,k
  real*8 a,b,c,aa,bb,cc,aaa,bbb,ccc,ded

  tipi=cels(iik)%ctipus

  ! we take a position for the new node in the direction of polarization
  ! this position is at a distance from the centroid similar to the more distant node existing in the direction
  ! of polarization

  ! so we make all the dot products between the polarization and each vector between the centroid and each node

  aax=1 ; aaw=1 ; aaz=1

  aa=cels(iik)%cex   ; bb=cels(iik)%cey   ; cc=cels(iik)%cez  
  aaa=cels(iik)%polx ; bbb=cels(iik)%poly ; ccc=cels(iik)%polz

  ded=sqrt(aaa**2+bbb**2+ccc**2)      ! >>> Is 28-6-14
!  if (ded==0.0d0) then
!    call random_number(a)   ! if the cell is unpolarized then the direction of mitosis is at random
!    k=int(a*nvaloq)+1
!    aaa=particions_esfera(k,1)
!    bbb=particions_esfera(k,2)
!    ccc=particions_esfera(k,3)
!    ded=sqrt(aaa**2+bbb**2+ccc**2)      ! >>> Is 28-6-14
!  end if

  ded=1d0/ded
  aaa=aaa*ded ; bbb=bbb*ded ; ccc=ccc*ded ! >>> Is 28-6-14

  dd=-1d10
  if (node(cels(iik)%node(1))%tipus==3) then
    do i=1,cels(iik)%nunodes-2
      ii=cels(iik)%node(i)
      a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
      d=sqrt(a**2+b**2+c**2)    ! >>> Is 1-8-14
      if (d<1d-2) cycle  ! >>> Is 1-8-14 ACHTUNG arbitrary value 
      d=1d0/d
      a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel4-7-14
      d=a*aaa+b*bbb+c*ccc
      if (d>dd) then ; dd=d ; jj=ii ; end if
    end do 
  else
    do i=1,cels(iik)%nunodes-1
      ii=cels(iik)%node(i)
      if (node(ii)%tipus==1) then
        a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
        d=sqrt(a**2+b**2+c**2)    ! >>> Is 1-8-14
        if (d<1d-2) cycle  ! >>> Is 1-8-14 ACHTUNG arbitrary value 
        d=1d0/d
        a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel4-7-14
        d=a*aaa+b*bbb+c*ccc
        if (d>dd) then ; dd=d ; jj=ii ; end if
      end if
    end do 
  end if

  if (dd/=0.0) then 
!    dd=dd*0.8d0  !THIS IS TO AVOID THAT CELLS CROSS EACH OTHER
    dd=0.8d0
    aa=aa+dd*(node(jj)%x-aa)
    bb=bb+dd*(node(jj)%y-bb)
    cc=cc+dd*(node(jj)%z-cc)
  else             ! the cell was apolar so we take the direction from the center to a random node in the cell
    if (cels(iik)%nunodes==1) then
      call random_number(a)
      a=a*node(cels(iik)%node(1))%add
      aa=cels(iik)%cex+a    !it is not spherically aleatory but this should happen very rarely
      bb=cels(iik)%cey+a
      cc=cels(iik)%cez+a
    else
      call random_number(a)   !that is the part that is different from the findepixyz_polar
      k=int(a*nvaloq)+1
      aaa=particions_esfera(k,1)
      bbb=particions_esfera(k,2)
      ccc=particions_esfera(k,3)
      dd=-1d10
      if (node(cels(iik)%node(1))%tipus==3) then
        do i=1,cels(iik)%nunodes
          ii=cels(iik)%node(i)
          a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
          d=1d0/sqrt(a**2+b**2+c**2)
          a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel9-7-14
          d=a*aaa+b*bbb+c*ccc
          if (d>dd) then ; dd=d ; jj=ii ; end if
        end do 
      else
        do i=1,cels(iik)%nunodes
          ii=cels(iik)%node(i)
          if (node(ii)%tipus==1) then
            a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
            d=1d0/sqrt(a**2+b**2+c**2)
            a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel9-7-14
            d=a*aaa+b*bbb+c*ccc
            if (d>dd) then ; dd=d ; jj=ii ; end if
          end if
        end do 
      end if
      aa=aa+dd*(node(jj)%x-aa)
      bb=bb+dd*(node(jj)%y-bb)
      cc=cc+dd*(node(jj)%z-cc)
    end if
  end if

  if (tipi>2) then
    if (tipi==3) then
      node(nd)%x=aa                
      node(nd)%y=bb
      node(nd)%z=cc
    end if
  else

    node(nd-1)%x=aa                
    node(nd-1)%y=bb
    node(nd-1)%z=cc

    ! now to take its altre we take the three nodes that are closer to the upper one
    a=extre ; b=extre ; c=extre
    do ii=1,cels(iik)%nunodes
      i=cels(iik)%node(ii)
      if (node(i)%tipus==2) then
        d=sqrt((aa-node(i)%x)**2+(bb-node(i)%y)**2+(cc-node(i)%z)**2)
        if (d<a) then
          if (a<b) then
            if (b<c) then ; c=b ; aaz=aaw ; end if      
            b=a ; aaw=aax
          end if       
          a=d ; aax=i
        else
          if (d<b) then
            if (b<c) then ; c=b ; aaz=aaw ; end if      
            b=d ; aaw=i
          else
            if (d<c) then ; b=d ; aaz=i ; end if
          end if
        end if
      end if
    end do
    !now we take the average of the vectors going from each of those three to their altre

    a=node(aax)%x-node(node(aax)%altre)%x
    b=node(aax)%y-node(node(aax)%altre)%y
    c=node(aax)%z-node(node(aax)%altre)%z
    a=a+node(aaw)%x-node(node(aaw)%altre)%x
    b=b+node(aaw)%y-node(node(aaw)%altre)%y
    c=c+node(aaw)%z-node(node(aaw)%altre)%z
    a=a+node(aaz)%x-node(node(aaz)%altre)%x
    b=b+node(aaz)%y-node(node(aaz)%altre)%y
    c=c+node(aaz)%z-node(node(aaz)%altre)%z
    node(nd)%x=node(nd-1)%x+third*a
    node(nd)%y=node(nd-1)%y+third*b
    node(nd)%z=node(nd-1)%z+third*c
  end if
end subroutine


!***********************************************************

subroutine findepixyz_nodel(iik)   ! it is like findepixyz but putting the new node in the direction of a randomly chosen existing node
  integer no
  integer iik,tipi,kul
  real*8 a,b,c,aa,bb,cc,aaa,bbb,ccc,dx,dy,dz,dd,ddd,d

  tipi=cels(iik)%ctipus

  ! we take a position for the new node in the direction of polarization
  ! this position is at a distance from the centroid similar to the more distant node existing in the direction
  ! of polarization

  ! so we make all the dot products between the polarization and each vector between the centroid and each node

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  aaa=particions_esfera(k,1)
  bbb=particions_esfera(k,2)
  ccc=particions_esfera(k,3) 

  aax=1 ; aaw=1 ; aaw=1

  aa=cels(iik)%cex   ; bb=cels(iik)%cey   ; cc=cels(iik)%cez   

  dd=-1d10
  if (node(cels(iik)%node(1))%tipus==3) then
    do i=1,cels(iik)%nunodes
      ii=cels(iik)%node(i)
      if (ii==nd) cycle  ! >>> Is 7-6-14
      a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
      d=sqrt(a**2+b**2+c**2)    ! >>> Is 1-8-14
      if (d<1d-2) cycle  ! >>> Is 1-8-14 ACHTUNG arbitrary value 
      a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel4-7-14
      d=a*aaa+b*bbb+c*ccc
      if (d>dd) then ; dd=d ; jj=ii ; end if
    end do  
  else
    do i=1,cels(iik)%nunodes
      ii=cels(iik)%node(i)
      if (ii==nd-1) cycle  ! >>> Is 7-6-14
      if (ii==nd) cycle    ! >>> Is 7-6-14
      if (node(ii)%tipus==1) then
        a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
        d=sqrt(a**2+b**2+c**2)    ! >>> Is 1-8-14
        if (d<1d-2) cycle  ! >>> Is 1-8-14 ACHTUNG arbitrary value 
        d=1d0/d
        a=a*d;b=b*d;c=c*d  !we need to make the unit vector !>>Miquel4-7-14
        d=a*aaa+b*bbb+c*ccc
        if (d>dd) then ; dd=d ; jj=ii ; end if
      end if
    end do  
  end if

  if (dd/=0.0) then  

    call random_number(a)
    dd=dd*a

!    dd=dd*0.9d0  !THIS IS TO AVOID THAT CELLS CROSS EACH OTHER

    aa=aa+dd*(node(jj)%x-aa)
    bb=bb+dd*(node(jj)%y-bb)
    cc=cc+dd*(node(jj)%z-cc)

  else           
    print *,"the position of nodes in cel",iik," is fucked up, they all in the same place or something"
  end if

  if (tipi>2) then
    if (tipi==3) then
      node(nd)%x=aa                 
      node(nd)%y=bb 
      node(nd)%z=cc 
      return  ! >>> Is 27-6-14
    end if
  else

    node(nd-1)%x=aa                 
    node(nd-1)%y=bb 
    node(nd-1)%z=cc
    a=extre ; b=extre ; c=extre

    ! now to take its altre we take the three nodes that are closer to the upper one
    do ii=1,cels(iik)%nunodes
      i=cels(iik)%node(ii)
      if (i==nd) cycle      ! >>> 15-6-14 this is because we are addin the node before 
      if (i==nd-1) cycle    ! >>> 15-6-14 
      if (node(i)%tipus==1) then ! >>> 15-6-14
        d=sqrt((aa-node(i)%x)**2+(bb-node(i)%y)**2+(cc-node(i)%z)**2)
        if (d<a) then
          if (a<b) then
            if (b<c) then ; c=b ; aaz=aaw ; end if       
            b=a ; aaw=aax 
          end if        
          a=d ; aax=i
        else
          if (d<b) then
            if (b<c) then ; c=b ; aaz=aaw ; end if       
            b=d ; aaw=i
          else
            if (d<c) then ; b=d ; aaz=i ; end if
          end if
        end if
      end if
    end do

    aax=node(aax)%altre  ! >>> Is 15-6-14
    aaw=node(aaw)%altre  ! >>> Is 15-6-14
    aaz=node(aaz)%altre  ! >>> Is 15-6-14

  end if
    ! now we take the average of the vectors going from each of those three to their altre 
    ! we interpolate so it means that we take the average weighted by the distance between the nd position and the closest and
    ! and then interpolate that in the other side node

  !d=(sqrt((aa-node(aax)%x)**2+(bb-node(aax)%y)**2+(cc-node(aax)%z)**2))
  !dd=(sqrt((aa-node(aaw)%x)**2+(bb-node(aaw)%y)**2+(cc-node(aaw)%z)**2))
  !ddd=(sqrt((aa-node(aax)%z)**2+(bb-node(aaz)%y)**2+(cc-node(aaz)%z)**2))
  !d=1
  !dd=1
  !ddd=1
  !   d=d/(d+dd+ddd)
  !   dd=dd/(d+dd+ddd)
  !   ddd=ddd/(d+dd+ddd)
    a=(node(aax)%x-node(node(aax)%altre)%x) !*d
    b=(node(aax)%y-node(node(aax)%altre)%y) !*d
    c=(node(aax)%z-node(node(aax)%altre)%z) !*d
    a=a+(node(aaw)%x-node(node(aaw)%altre)%x) !*dd
    b=b+(node(aaw)%y-node(node(aaw)%altre)%y) !*dd
    c=c+(node(aaw)%z-node(node(aaw)%altre)%z) !*dd
    a=a+(node(aaz)%x-node(node(aaz)%altre)%x) !*ddd
    b=b+(node(aaz)%y-node(node(aaz)%altre)%y) !*ddd
    c=c+(node(aaz)%z-node(node(aaz)%altre)%z) !*ddd
    node(nd)%x=node(nd-1)%x+a*third
    node(nd)%y=node(nd-1)%y+b*third
    node(nd)%z=node(nd-1)%z+c*third
end subroutine

!********************************************************************************************

subroutine findepixyz_out(iik)   ! it is like findepixyz but putting the new node in the direction of a randomly chosen existing node
  integer no,i,j,k,ii,jj,kk,iii,jjj,kkk
  integer iik,tipi
  real*8 a,b,c,aa,bb,cc,aaa,bbb,ccc,d,dd,cex,cey,cez
  integer sol,ntrials
  real*8 enepos(cels(iik)%nunodes),pos(cels(iik)%nunodes,2,3)

  ntrials=cels(iik)%nunodes
  enepos=1d10
  pos=0
  tipi=node(cels(iik)%node(1))%tipus

  a=0 ; b=0 ; c=0
  do j=1,cels(iik)%nunodes
    k=cels(iik)%node(j)
    if((node(k)%tipus==2.or.node(k)%tipus==3).and.k<nd-1)then
      a=a+node(k)%x ; b=b+node(k)%y ; c=c+node(k)%z
    end if
  end do
  d=1d0/real(cels(iik)%nunodes)
  if(node(cels(iik)%node(1))%tipus<3) d=2d0*d !if it's epithelial
  cex=a*d ; cey=b*d ; cez=c*d

  do sol=1,ntrials !the arbitrary is the number of different locations of the new nodes we try to decide which is the best energetically

    jj=cels(iik)%node(sol)
    enepos(sol)=1d10
    if (node(jj)%tipus==1)  cycle  !>>> Is 15-6-14
    if (jj==nd) cycle
    if (node(jj)%tipus<3.and.jj==nd-1) cycle ! >>> Is 28-6-14

    dd=0.01d0  !THIS IS THE PART THAT IS DIFFRENT ARBITRARY VALUE!!!!!!!!!!!!!!!!!!!! Is 

    aa=node(jj)%x+dd*(node(jj)%x-cex)
    bb=node(jj)%y+dd*(node(jj)%y-cey)
    cc=node(jj)%z+dd*(node(jj)%z-cez)
    if (tipi>2) then
      if (tipi==3) then
        node(nd)%x=aa                 
        node(nd)%y=bb 
        node(nd)%z=cc 
        node(nd)%eqd=node(jj)%eqd
        node(nd)%add=node(jj)%add
        node(nd)%grd=node(jj)%grd
        if(nd>1) call neighbor_build

        call energia(nd)
        enepos(sol)=node(nd)%e
        pos(sol,1,1)=aa ; pos(sol,1,2)=bb ; pos(sol,1,3)=cc
      end if
    else
  
      node(nd)%x=aa                 
      node(nd)%y=bb 
      node(nd)%z=cc

      !node(nd-1)%x=aa+1d-7 ! >>> Is 15-6-14 this value is totally arbitrary and it is simply so that iniboxesll does no crush                 
      !node(nd-1)%y=bb+1d-7 ! >>> Is 15-6-14 soon we give the proper values to that                
      !node(nd-1)%z=cc+1d-7 ! >>> Is 15-6-14                 

      node(nd)%eqd=node(jj)%eqd
      node(nd)%add=node(jj)%add
      node(nd)%grd=node(jj)%grd
      pos(sol,1,1)=aa ; pos(sol,1,2)=bb ; pos(sol,1,3)=cc

      ! now we have to find the position of the altre node
      ! we do it by taking the vector from the upper to the lower for the already existing nodes in the cell

      ii=jj
      a=(node(node(ii)%altre)%x-node(ii)%x)
      b=(node(node(ii)%altre)%y-node(ii)%y)
      c=(node(node(ii)%altre)%z-node(ii)%z)

      node(nd-1)%x=aa+a                 
      node(nd-1)%y=bb+b 
      node(nd-1)%z=cc+c

      node(nd-1)%eqd=node(node(jj)%altre)%eqd
      node(nd-1)%add=node(node(jj)%altre)%add
      node(nd-1)%grd=node(node(jj)%altre)%eqd

      if(nd>1) call neighbor_build

      call energia(nd)
      enepos(sol)=node(nd)%e
      call energia(nd-1)
      enepos(sol)=enepos(sol)+node(nd-1)%e

      pos(sol,2,1)=aa+a ; pos(sol,2,2)=bb+b ; pos(sol,2,3)=cc+c

    end if
  end do

  ! now we chose the position with the lower energy
  if (node(nd)%tipus==3) then
    a=huge(a)
    do i=1,ntrials
      if (a>enepos(i)) then ; j=i ; a=enepos(i) ; end if
    end do
    node(nd)%x=pos(j,1,1)                 
    node(nd)%y=pos(j,1,2)
    node(nd)%z=pos(j,1,3)
  else
    a=1d10
    do i=1,ntrials
      if (a>enepos(i)) then ; j=i ; a=enepos(i) ; end if
    end do

    node(nd)%x=pos(j,1,1)                 
    node(nd)%y=pos(j,1,2)
    node(nd)%z=pos(j,1,3)
    node(nd-1)%x=pos(j,2,1)                 
    node(nd-1)%y=pos(j,2,2)
    node(nd-1)%z=pos(j,2,3)
  end if

  if(nd>1) call neighbor_build

end subroutine

!*****************************************************************************
! >>> Is 5-6-14
subroutine findepixyz_out_no_energy(iik)   ! it is like findepixyz but putting the new node in the direction of a randomly chosen existing node
  integer no
  integer iik,tipi
  real*8 a,b,c,aa,bb,cc,aaa,bbb,ccc

  tipi=cels(iik)%ctipus

  ! we take a position for the new node in the direction of polarization
  ! this position is at a distance from the centroid similar to the more distant node existing in the direction
  ! of polarization

  ! so we make all the dot products between the polarization and each vector between the centroid and each node

  call random_number(a)   !that is the part that is different from the findepixyz_polar
  k=int(a*nvaloq)+1
  aaa=particions_esfera(k,1)
  bbb=particions_esfera(k,2)
  ccc=particions_esfera(k,3) 

  aax=1 ; aaw=1 ; aaw=1

  aa=cels(iik)%cex   ; bb=cels(iik)%cey   ; cc=cels(iik)%cez   

  dd=-1d10
  if (node(cels(iik)%node(1))%tipus==3) then
    do i=1,cels(iik)%nunodes
      ii=cels(iik)%node(i)
      a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
      d=a*aaa+b*bbb+c*ccc
      if (d>dd) then ; dd=d ; jj=ii ; end if
    end do  
  else
    do i=1,cels(iik)%nunodes
      ii=cels(iik)%node(i)
      if (node(ii)%tipus==1) then
        a=node(ii)%x-aa ; b=node(ii)%y-bb ; c=node(ii)%z-cc
        d=a*aaa+b*bbb+c*ccc
        if (d>dd) then ; dd=d ; jj=ii ; end if
      end if
    end do  
  end if

  if (dd/=0.0) then  

    call random_number(a)
!    dd=dd*a

    dd=dd*0.01d0  !THIS IS THE PART THAT IS DIFFRENT

    if (node(jj)%tipus==1) then
      jj=node(jj)%altre
    end if

    aa=node(jj)%x+dd*(node(jj)%x-aa)
    bb=node(jj)%y+dd*(node(jj)%y-bb)
    cc=node(jj)%z+dd*(node(jj)%z-cc)

  else             ! the cell was apolar so we take the direction from the center to a random node in the cell 

    if (cels(iik)%nunodes==1) then
      call random_number(a)
      a=a*node(cels(iik)%node(1))%add
      aa=cels(iik)%cex+a    !it is not spherically aleatory but this should happen very rarely
      bb=cels(iik)%cey+a
      cc=cels(iik)%cez+a
    else
13    call random_number(a)
      k=int(a*cels(iik)%nunodes)+1
      aa=(cels(iik)%cex+node(k)%x)*0.5
      bb=(cels(iik)%cey+node(k)%y)*0.5
      cc=(cels(iik)%cez+node(k)%z)*0.5
      if (sqrt((node(ii)%x-aa)**2+(node(ii)%y-bb)**2+(node(ii)%z-cc)**2)==0.0d0) goto 13 !this could happen sometimes
    end if
  end if

  if (tipi>2) then
    if (tipi==3) then
      node(nd)%x=aa                 
      node(nd)%y=bb 
      node(nd)%z=cc 
    end if
  else

    node(nd)%x=aa                 
    node(nd)%y=bb 
    node(nd)%z=cc

    ! now we have to find the position of the altre node
    ! we do it by taking the vector from the upper to the lower for the already existing nodes in the cell

    ii=jj
    a=(node(node(ii)%altre)%x-node(ii)%x)
    b=(node(node(ii)%altre)%y-node(ii)%y)
    c=(node(node(ii)%altre)%z-node(ii)%z)

    node(nd-1)%x=aa+a                 
    node(nd-1)%y=bb+b 
    node(nd-1)%z=cc+c

  end if
end subroutine

!***********************************************************

subroutine addanode(i)
integer::i,j,ic,icc,ii,jj,kk,kkk,kkkk,kkkkk,tipi,celi,cont
integer,allocatable::clist(:)
real*8,allocatable::cgex(:,:)
real*8::ix,iy,iz,cx,cy,cz,ccx,ccy,ccz,d,dd,pesc,maxd,ee,aa,bb,cc,aaa,bbb,ccc,b
real*8 ox,oy,oz !>>> Is 2-9-13
integer:: imaxd,iimaxd
integer:: cnodecel(cels(node(i)%icel)%nodela)
type(nod),allocatable :: cnode(:),cpnode(:,:),cnodeo(:)
real*8,allocatable::cpx(:),cpy(:),cpz(:),cdex(:)               !>>>>>Miquel 17-6-13
real*8,allocatable::cvcil(:),cvtor(:),cvstor(:),cvspr(:)              !>>>>>Miquel 20-6-13
integer,allocatable::cneigh(:,:),cnneigh(:)                    !<<< Is 30-4-15
real*8,allocatable::cdneigh(:,:)                               !>>> Is 30-4-15
real*8,allocatable::cfmeanv(:),cfmeanl(:),cdidpol(:,:)

  iik=node(i)%icel
  tipi=node(i)%tipus

  !if(ffu(8)==0)then
  !  ! now we have to put the new nodes in the boxes
  !  ic=boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))!position in boxes
  !  if(ic==0)then
  !    !the cube was empty
  !    boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))=nd
  !    list(nd)=0
  !  else		
  !    !the cube was not empty so we search for the end																
  !    do
  !      icc=ic; ic=list(icc) ; if(ic==0)then ; list(icc)=nd ; list(nd)=0 ; exit ; end if
  !    end do
  !  end if
  !end if

  !copy the data in the new node from its parent
  ox=node(nd)%x ; oy=node(nd)%y ; oz=node(nd)%z  !>>> Is 2-9-13
  node(nd)=node(i)   ! we copy everything        !>>> Is 2-9-13
  node(nd)%tipus=node(i)%tipus
  node(nd)%x=ox ; node(nd)%y=oy ; node(nd)%z=oz  !except the position !>>> Is 2-9-13 
  node(nd)%marge=1   !so the nucleus doesn't duplicate >>>Miquel4-10-13

  !now nodeo
  ox=nodeo(nd)%x ; oy=nodeo(nd)%y ; oz=nodeo(nd)%z  !>>> Is 2-9-13
  nodeo(nd)=nodeo(i)   ! we copy everything        !>>> Is 2-9-13
  nodeo(nd)%tipus=node(i)%tipus
  nodeo(nd)%x=ox ; nodeo(nd)%y=oy ; nodeo(nd)%z=oz  !except the position !>>> Is 2-9-13 
  nodeo(nd)%marge=1   !so the nucleus doesn't duplicate >>>Miquel4-10-13

  !the new node is negligibly small and and empty of any gene product

  !if(tipi==2)then
    !gex(nd,:)=(gex(aax,:)+gex(aaz,:)+gex(aaw,:))*third                            !>>>>>>>>>>> Miquel 3-6-13
    !gex(nd,:)=0d0                                                       !>>>>> Miquel 9-10-13
  !else if(tipi==1)then
    !iii=node(aax)%altre;jjj=node(aaz)%altre;kkk=node(aaw)%altre
    !gex(nd,:)=(gex(iii,:)+gex(jjj,:)+gex(kkk,:))*third
    !gex(nd,:)=0d0                                                       !>>>>> Miquel 9-10-13
  !else
    !gex(nd,:)=gex(i,:)
    !gex(nd,:)=0d0                                                       !>>>>> Miquel 9-10-13
  !end if

  if(ffu(1)==1)then
    !agex(nd,:)=0d0
    !gex(nd,:)=0d0  !>>> Is 13-5-14 THIS IS TRICKY since this gex values is used for the rest of the nexus module 
    agex(nd,:)=agex(i,:)
    gex(nd,:)=gex(i,:)
  else     
        !NO in 'single node per cell' mode , growth equals division, so the gene products split between the two nodes
!    agex(i,:)=0.5*agex(i,:)  
!    agex(i,:)=agex(i,:)
    agex(nd,:)=agex(i,:)
!    gex(i,:)=0.5*gex(i,:)
!    gex(i,:)=gex(i,:)
    gex(nd,:)=gex(i,:)    
  end if
  
  if(tipi<3)then													!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
    node(nd)%altre=nd-1
  else															!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
    node(nd)%altre=0											!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
  end if															!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13

  nneigh(nd)=0 !initialize the neighbors to not get errors in other parts
  neigh(nd,:)=0 !>>>Miquel14-3-14
  fmeanl(nd)=0d0 !>>>Miquel20-4-15
  
  if(tipi<3)then													!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13

    !if(ffu(8)==0)then
    !  ic=boxes(nint(node(nd-1)%x*urv),nint(node(nd-1)%y*urv),nint(node(nd-1)%z*urv))!position in boxes
    !  if(ic==0)then
    !    !the cube was empty
    !    boxes(nint(node(nd-1)%x*urv),nint(node(nd-1)%y*urv),nint(node(nd-1)%z*urv))=nd
    !    list(nd-1)=0
    !  else		
    !    !the cube was not empty so we search for the end																
    !    do
    !      icc=ic ; ic=list(icc) ; if(ic==0)then ; list(icc)=nd-1 ; list(nd)=0 ; exit ; end if
    !    end do
    !  end if
    !end if

    !copy the data in the new node from its parent
    node(nd-1)%icel=iik  
    j=node(i)%altre
    ox=node(nd-1)%x ; oy=node(nd-1)%y ; oz=node(nd-1)%z  !>>> Is 2-9-13
    node(nd-1)=node(j)   ! we copy everything        !>>> Is 2-9-13
    node(nd-1)%x=ox ; node(nd-1)%y=oy ; node(nd-1)%z=oz  !except the position !>>> Is 2-9-13 
    node(nd-1)%marge=1   !so the nucleus doesn't duplicate >>>Miquel4-10-13

    node(nd-1)%altre=nd

    !copy the data in the new node from its parent ALSO FOR nodeo so that it can be modified by nexus from that value
    !nodeo(nd-1)%icel=iik  
    j=nodeo(i)%altre
    ox=nodeo(nd-1)%x ; oy=nodeo(nd-1)%y ; oz=nodeo(nd-1)%z  !>>> Is 2-9-13

    nodeo(nd-1)=node(j)   ! we copy everything        !>>> Is 2-9-13
    nodeo(nd-1)=nodeo(j)  ! we copy everything        !>>> Is 2-9-13
    nodeo(nd-1)%x=ox ; nodeo(nd-1)%y=oy ; nodeo(nd-1)%z=oz  !except the position !>>> Is 2-9-13 
    nodeo(nd-1)%marge=1   !so the nucleus doesn't duplicate >>>Miquel4-10-13

    node(nd-1)%altre=nd

  
    nneigh(nd-1)=0 !initialize the neighbors to not get errors in other parts
    neigh(nd-1,:)=0 !>>>Miquel14-3-14
    fmeanl(nd-1)=0d0 !>>>Miquel20-4-15

    !if(tipi==1)then
      !gex(nd-1,:)=(gex(aax,:)+gex(aaz,:)+gex(aaw,:))*third                            !>>>>>>>>>>> Miquel 3-6-13
      !gex(nd-1,:)=0d0                                                       !>>>>> Miquel 9-10-13
    !else
      !iii=node(aax)%altre;jjj=node(aaz)%altre;kkk=node(aaw)%altre
      !gex(nd-1,:)=(gex(iii,:)+gex(jjj,:)+gex(kkk,:))*third
      !gex(nd-1,:)=0d0                                                       !>>>>> Miquel 9-10-13
    !end if
    j=node(i)%altre    !!>>HC 4-10-2022 We are going to copy the gene concentration from the mother node j to the daughter nd-1
    if(ffu(1)==1)then
      !agex(nd-1,:)=0d0
      !gex(nd-1,:)=0d0   !>>> Is 13-5-14 THIS IS TRICKY since this gex values is used for the rest of the nexus module      
      agex(nd-1,:)=agex(j,:) !>>> Is 10-10-14
      gex(nd-1,:)=gex(j,:)   !>>> Is 10-10-14
    else     !in 'single node per cell' mode , growth equals division, so the gene products split between the two nodes
      !agex(j,:)=0.5*agex(j,:)
      agex(nd-1,:)=agex(j,:)     !>>> Is 10-10-14
      !gex(j,:)=0.5*gex(j,:)
      gex(nd-1,:)=gex(j,:)   !>>> Is 10-10-14
    end if

  end if										!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13

  !reallocatation of the cels()%node matrix
  if(tipi<3)then
    kk=cels(iik)%nunodes
    cels(iik)%nunodes=kk+2
    if (kk+4>=cels(iik)%nodela) cels(iik)%nodela=kk+20   
    cnodecel(1:kk)=cels(iik)%node(1:kk)
    deallocate(cels(iik)%node)
    allocate(cels(iik)%node(cels(iik)%nodela))
    cels(iik)%node(1:kk)=cnodecel(1:kk)
    cels(iik)%node(kk+1)=nd-1
    cels(iik)%node(kk+2)=nd
  else
    kk=cels(iik)%nunodes
    cels(iik)%nunodes=kk+1
    if (kk+2>=cels(iik)%nodela) cels(iik)%nodela=kk+10   
    cnodecel(1:kk)=cels(iik)%node(1:kk)
    deallocate(cels(iik)%node)
    allocate(cels(iik)%node(cels(iik)%nodela))
    cels(iik)%node(1:kk)=cnodecel(1:kk)
    cels(iik)%node(kk+1)=nd
  end if

  if(nd+2>=nda)then	!let's enhance the node matrix and list matrix
    nda=nda+10
    allocate(cnode(nd))
    cnode(1:nd)=node(1:nd)
    deallocate(node)
    allocate(node(nda))
    node(1:nd)=cnode(1:nd)  
    deallocate(cnode)

    allocate(cpnode(mamax,nd))	!!>>Miquel 16-10-12
    cpnode(:,nd)=pnode(:,nd)
    deallocate(pnode)
    allocate(pnode(mamax,nda))
    pnode(:,1:nd)=cpnode(:,1:nd)
    deallocate(cpnode)
    allocate(cgex(nda,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
    cgex(:nd,:ng)=gex(:nd,:ng)
    deallocate(gex);allocate(gex(nda,ng))
    gex=0
    gex(:nd,:ng)=cgex(:nd,:ng)
    
    cgex(:nd,:ng)=agex(:nd,:ng)
    deallocate(agex) ; allocate(agex(nda,ng)) ;agex=0
    agex(:nd,:ng)=cgex(:nd,:ng)
    
    deallocate(cgex)
    allocate(cpx(nd),cpy(nd),cpz(nd))
    cpx(1:nd)=px(1:nd)
    cpy(1:nd)=py(1:nd)
    cpz(1:nd)=pz(1:nd)
    deallocate(px,py,pz);allocate(px(nda),py(nda),pz(nda))
    px=0;py=0;pz=0
    px(1:nd)=cpx(1:nd)
    py(1:nd)=cpy(1:nd)
    pz(1:nd)=cpz(1:nd)
    allocate(cdex(nd))
    cdex(1:nd)=dex(1:nd)
    deallocate(dex);allocate(dex(nda))
    dex=0
    dex(1:nd)=cdex(1:nd)

    allocate(cvcil(nd),cvtor(nd),cvstor(nd),cvspr(nd))  !force components storing arrays !>>>>>>>>>>>>>> Miquel 20-6-13
    allocate(cneigh(nda,omnn))
    cneigh(1:nd,:)=neigh(1:nd,:)
    deallocate(neigh)
    allocate(neigh(nda,omnn))
    neigh(:nd,:)=cneigh(:nd,:)
    deallocate(cneigh)
    allocate(cdneigh(nda,omnn))
    cdneigh(1:nd,:)=dneigh(1:nd,:)
    deallocate(dneigh)
    allocate(dneigh(nda,omnn))
    dneigh(:nd,:)=cdneigh(:nd,:)
    deallocate(cdneigh)
    allocate(cnneigh(nd))
    cnneigh(1:nd)=nneigh(1:nd)
    deallocate(nneigh)
    allocate(nneigh(nda))
    nneigh(1:nd)=cnneigh(1:nd)
    deallocate(cnneigh)

    if(ffu(8)==0)then
      !allocate(cnneigh(nd))
      !cnneigh(1:nd)=dif_nneigh(1:nd)
      !deallocate(dif_nneigh)
      !allocate(dif_nneigh(nda))
      !dif_nneigh(1:nd)=cnneigh(1:nd)
      !deallocate(cnneigh)
      
      allocate(clist(nd))
      clist(1:nd)=list(1:nd)
      deallocate(list)
      allocate(list(nda))
      list=0
      list(1:nd)=clist(1:nd)
      deallocate(clist)
    end if

    cvcil(1:nd)=vcilx(1:nd)
    cvtor(1:nd)=vtorx(1:nd)
    cvstor(1:nd)=vstorx(1:nd)
    cvspr(1:nd)=vsprx(1:nd)
    deallocate(vcilx,vtorx,vstorx,vsprx)
    allocate(vcilx(nda),vtorx(nda),vstorx(nda),vsprx(nda))
    vcilx(1:nd)=cvcil(1:nd)
    vtorx(1:nd)=cvtor(1:nd)
    vstorx(1:nd)=cvstor(1:nd)
    vsprx(1:nd)=cvspr(1:nd)
    cvcil(1:nd)=vcily(1:nd)
    cvtor(1:nd)=vtory(1:nd)
    cvstor(1:nd)=vstory(1:nd)
    deallocate(vcily,vtory,vstory,vspry)
    allocate(vcily(nda),vtory(nda),vstory(nda),vspry(nda))
    vcily(1:nd)=cvcil(1:nd)
    vtory(1:nd)=cvtor(1:nd)
    vstory(1:nd)=cvstor(1:nd)
    vspry(1:nd)=cvspr(1:nd)
    cvcil(1:nd)=vcilz(1:nd)
    cvtor(1:nd)=vtorz(1:nd)
    cvstor(1:nd)=vstorz(1:nd)
    deallocate(vcilz,vtorz,vstorz,vsprz)
    allocate(vcilz(nda),vtorz(nda),vstorz(nda),vsprz(nda))
    vcilz(1:nd)=cvcil(1:nd)
    vtorz(1:nd)=cvtor(1:nd)
    vstorz(1:nd)=cvstor(1:nd)
    vsprz(1:nd)=cvspr(1:nd)
    deallocate(cvcil,cvtor,cvstor,cvspr)
    if (allocated(erep)) deallocate(erep)
    if (allocated(erepcel)) deallocate(erepcel)
    if (allocated(eadh)) deallocate(eadh)
    if (allocated(eyou)) deallocate(eyou)
    if (allocated(espring)) deallocate(espring)
    if (allocated(eadh)) deallocate(eadh)
    if (allocated(etor)) deallocate(etor)
    allocate(erep(nda),erepcel(nda),eadh(nda),eyou(nda),espring(nda),etor(nda))

    allocate(cfmeanv(nd),cfmeanl(nd))
    cfmeanv(1:nd)=fmeanv(1:nd) ; cfmeanl(1:nd)=fmeanl(1:nd)
    deallocate(fmeanv,fmeanl)
    allocate(fmeanv(nda),fmeanl(nda))
    fmeanv(1:nd)=cfmeanv(1:nd) ; fmeanl(1:nd)=cfmeanl(1:nd)
    deallocate(cfmeanv,cfmeanl)

    ! we have to re-size nodeo too so that nexus can modify the new nodes on the bases of their values when they first
    ! arose
    allocate(cnodeo(nd))
    cnodeo(:nd)=nodeo(:nd)
    deallocate(nodeo)
    allocate(nodeo(nda))
    nodeo(:nd)=cnodeo(:nd)
    nodeo(:nd+1:nda)%ecm=0.0d0
    node(:nd+1:nda)%ecm=0.0d0
    nodeo(:nd+1:nda)%fix=0.0d0
    node(:nd+1:nda)%fix=0.0d0
  end if		

  !for nodeo
  if (tipi<3) then
    if (ffu(1)==0) node(nd)%marge=0 ! >>> Is 6-10-14 it is necessary if each mesenchymal cell is a single node
  end if

end subroutine


end module growth
