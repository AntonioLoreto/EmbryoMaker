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

module editor
use general
use ic
!use opengl
use opengl_gl
use opengl_glu
!use opengl_glut ! not used RZ 17-11-14
use opengl_kinds
use growth
!use del
use nexus
use io
use neighboring
use mitosis

integer,public :: nodeindex, cellid,nodetarget,target
public :: addnode,passvalues,selectnode,pastenode,getcoord,update,addcell,addepi,&
          selectcell,pastecell,deletenode,deletecell

integer, private :: and
type(nod), public :: tempnod !Tommi 23.9.2013
type(cel), public :: tempcel !Tommi 23.9.2013


contains

!***************************************************************************************************

subroutine addecm(x,y,z) !adds a basic mesenchymal node with similar properties than the previuos node in the node matrix
  real*8 x,y,z

  call ecm_secretion(1)

  node(nd)%tipus=4
  node(nd)%x=x
  node(nd)%y=y
  node(nd)%z=z
end subroutine

!********************************************************************************************************

subroutine addmescell(x,y,z)
  real*8 x,y,z
  integer i

  call addcell(x,y,z)
  cels(:)%temt=0.0d0
  cels(ncels)%temt=2
  call emt
  cels(ncels)%ctipus=3
!  do i=1,cels(ncels)%nunodes
!    node(i)%tipus=3
!  end do
end subroutine

!***************************************************************************************************

subroutine addnode(x,y,z) !adds a basic node 
real*8 x,y,z,ax,ay,az
integer dal,sot ! RZ 17-11-14 added

if (cellid<1.or.cellid>ncels) then
  print *,"you should select a cell first"
  return
end if

if (cels(cellid)%ctipus==3) then
  call adding(x,y,z)
  node(nd)%tipus=3
else
  call adding(x,y,z)
  node(nd)%tipus=1
print *,x,y,z
  !now we estimate the position of the other side of the epithelium based on the vector going from one side to the other for node 1 in that cell
  dal=int(cels(cellid)%node(1)) ! RZ 17-11-14 INT
  sot=node(dal)%altre
  ax=node(sot)%x-node(dal)%x
  ay=node(sot)%y-node(dal)%y
  az=node(sot)%z-node(dal)%z
print *,x+ax,y+ay,z+az
  call adding(x+ax,y+ay,z+az)
  node(nd)%tipus=2
  node(nd)%altre=nd-1
  node(nd-1)%altre=nd
end if

end subroutine

!*******************************************************************************************************************

subroutine adding(x,y,z) !adds a node 
real*8, intent(in) :: x,y,z
real*8::you,adh, rep, repcel,req,reqcel,reqs,ke,tor,stor,x1,y1,z1
integer::fase,tipus,icel,altre,marge,ic,icc,cell,and
type(nod),allocatable :: cnode(:),cpnode(:,:),cnodeo(:)
integer,allocatable::clist(:),ccelnode(:)
real*8,allocatable::cgex(:,:)


and=nd


  if(nd+1>=nda)then	!let's enhance the node matrix and list matrix
		nda=nda+10
		allocate(cnode(nd))
		cnode=node
		deallocate(node)
		allocate(node(nda))
		node(1:nd)=cnode(1:nd)
		deallocate(cnode)


		allocate(clist(nd))
		clist=list
		deallocate(list)
		allocate(list(nda))
		list=0
		list(1:nd)=clist(1:nd)
		deallocate(clist)

		allocate(cpnode(mamax,nd))	!!>>Miquel 16-10-12
		cpnode=pnode
		deallocate(pnode)
		allocate(pnode(mamax,nda))
		pnode(:,1:nd)=cpnode(:,1:nd)
		deallocate(cpnode)


        allocate(cgex(nd,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
        cgex(1:nd,1:ng)=gex(1:nd,1:ng)
        deallocate(gex);allocate(gex(nda,ng))
        gex=0
        gex(1:nd,1:ng)=cgex(1:nd,1:ng)
        deallocate(cgex)
  end if

  nd=nd+1
  ndmes=ndmes+1
  nda=nda+1
print *,nd,nda,"nd nda",nd-1,ubound(node)
  node(nd)=node(nd-1)!mesenchymal cell gets epithelial properties, since nd-1 is epithelial.
  node(nd)%x=x
  node(nd)%y=y
  node(nd)%z=z
  node(nd)%altre=0
  node(nd)%icel=cellid


  !call iniboxes
  !ic=boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))!position in boxes
  !if(ic==0)then
  !  !the cube was empty
  !  boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))=nd
  !  list(nd)=0
  !else		
  !  !the cube was not empty so we search for the end																
  !  do
  !    icc=ic; ic=list(icc) ; if(ic==0)then ; list(icc)=nd ; list(nd)=0 ; exit ; end if
  !  end do
  !end if

  gex(nd,:)=gex(nd-1,:) 



  cell=cellid
  j=cels(cell)%nunodes !<<<<< 11.9.2013
  i=cels(cell)%nodela
  allocate(ccelnode(i))
  ccelnode(1:i)=cels(cell)%node(1:i)
  deallocate(cels(cell)%node)
  cels(cell)%nunodes=j+1
  cels(cell)%nodela=i+1
  allocate(cels(cell)%node(i+1))
  cels(cell)%node(1:i)=ccelnode(1:i)
  cels(cell)%node(j+1)=nd
  cels(cell)%node(i+1)=0 ! >>>>> 11.9.2013

  call remade
  if(nd>1) call neighbor_build

end subroutine

!***********************************************************************************************************************

subroutine remade
type(nod),allocatable :: cnodeo(:)
real*8, allocatable   :: cgex(:,:)
deallocate(nneigh)
deallocate(neigh)
deallocate(dneigh)
allocate(nneigh(nda))
allocate(dneigh(nda,omnn))
allocate(neigh(nda,omnn))

deallocate(fmeanl)
allocate(fmeanl(nda))
deallocate(fmeanv)
allocate(fmeanv(nda))

deallocate(vsprx)
allocate(vsprx(nda))
deallocate(vspry)
allocate(vspry(nda))
deallocate(vsprz)
allocate(vsprz(nda))

deallocate(vcilx)
allocate(vcilx(nda))
deallocate(vcily)
allocate(vcily(nda))
deallocate(vcilz)
allocate(vcilz(nda))

deallocate(vtorx)
allocate(vtorx(nda))
deallocate(vtory)
allocate(vtory(nda))
deallocate(vtorz)
allocate(vtorz(nda))

deallocate(vstorx)
allocate(vstorx(nda))
deallocate(vstory)
allocate(vstory(nda))
deallocate(vstorz)
allocate(vstorz(nda))

deallocate(dex)
allocate(dex(nda))

deallocate(px)
allocate(px(nda))
deallocate(py)
allocate(py(nda))
deallocate(pz)
allocate(pz(nda))

allocate(cnodeo(nda))
cnodeo(:and)=nodeo(:and)
deallocate(nodeo)
allocate(nodeo(nda))
nodeo(:and)=cnodeo(:and)
print *,and,nd,nda,"and nd nda"
nodeo(and+1:nd)=node(and+1:nd)

allocate(cgex(nda,ng))
cgex(:and,:)=gex(:and,:)
deallocate(gex)
allocate(gex(nda,ng))
gex(:and,:)=cgex(:and,:)

!call neighbor_build

end subroutine

!***********************************************************************************************************************
subroutine selectnode(x,y,z) !selects a node nearest to the given coordinates
real*8, intent(in) :: x,y,z
real*8 :: dummyx,dummyy,dummyz, minimumvalue
integer :: i, sizelist
!logical :: inorder
type :: nodelist
     real*8 :: distance
     integer :: nodeid
end type
type(nodelist), allocatable :: nodescalc(:) !comparelist
cellid=0
allocate(nodescalc(nd))
do i=1,nd
  dummyx=node(i)%x
  dummyy=node(i)%y
  dummyz=node(i)%z
  nodescalc(i)%distance=sqrt(((x-dummyx)**2)+((y-dummyy)**2)+((z-dummyz)**2)) !d = sqrt( (x[1]-x[2])² + (y[1]-y[2])² + (z[1]-z[2])², euclidian distance
  nodescalc(i)%nodeid=i
end do   

minimumvalue=minval(nodescalc%distance)
do i=1,nd
  if(minimumvalue==nodescalc(i)%distance) then
    if (target==0) then
      nodeindex = nodescalc(i)%nodeid
      tempnod=node(nodeindex)
    else
      nodetarget = nodescalc(i)%nodeid
    end if
  endif
end do
!return nodeindex
if (cellid == 0) then
print*,"The selected node is node number ", nodeindex!<<<<<<<Tommi 17.9.2013
print*,"Position of the node",node(nodeindex)%x,node(nodeindex)%y,node(nodeindex)%z
print*,"Energy of the selected node",node(nodeindex)%e
print*,"Tipus of the selected node",node(nodeindex)%tipus
print*,"Cell of the selected node",node(nodeindex)%icel
print *,"Altre of the selected node",node(nodeindex)%altre!>>>>>>>Tommi

!we copy its data

endif
end subroutine selectnode

!*************************************************************************************************

subroutine pastenode(x,y,z) !pastes a selected node. The pasted node will have the same properties as the selected node, except for position and altre.
  real*8 x,y,z,ax,ay,az,a,b,c
  integer cella,ol,dal,sot ! RZ 17-11-14 added dal,sot

  if (nodeindex<=0) then
    print*,"No node selected!"
  else
    if (node(nodeindex)%tipus<3) then
      ndepi=ndepi+1
      ol=nodeindex
      call pasti(x,y,z)
      ndepi=ndepi+1
      nodeindex=node(ol)%altre
      !now we estimate the position of the other side of the epithelium based on the vector going from one side to the other for node 1 in that cell
      cella=node(nodeindex)%icel
      dal=cels(cella)%node(1)
      sot=node(dal)%altre
      ax=node(sot)%x-node(dal)%x
      ay=node(sot)%y-node(dal)%y
      az=node(sot)%z-node(dal)%z
      call pasti(x+ax,y+ay,z+az)
!      node(nd)%tipus=2
      node(nd)%altre=nd-1
      node(nd-1)%altre=nd
    else 
      ndmes=ndmes+1
      call pasti(x,y,z)
    end if
  endif

end subroutine

!*************************************************************************************************

subroutine pasti(x,y,z) !pastes a selected node. The pasted node will have the same properties as the selected node, except for position and altre.
real*8, intent(in) :: x,y,z
integer::fase,tipus,icel,altre,marge,ic,icc,cell ! RZ 17-11-14 added cell
type(nod),allocatable :: cnode(:),cpnode(:,:)
integer,allocatable::clist(:),ccelnode(:)
real*8,allocatable::cgex(:,:)
type(nod) :: temp 

  nd=nd+1
  node(nd)=node(nodeindex)
  node(nd)%x=x
  node(nd)%y=y
  node(nd)%z=z
  cell=node(nd)%icel

  gex(nd,:)=gex(nodeindex,:) 

  if(nd+1>=nda)then	!let's enhance the node matrix and list matrix
		nda=nda+10
		allocate(cnode(nd))
		cnode=node
		deallocate(node)
		allocate(node(nda))
		node(1:nd)=cnode(1:nd)
		deallocate(cnode)


		allocate(clist(nd))
		clist=list
		deallocate(list)
		allocate(list(nda))
		list=0
		list(1:nd)=clist(1:nd)
		deallocate(clist)

		allocate(cpnode(mamax,nd))	!!>>Miquel 16-10-12
		cpnode=pnode
		deallocate(pnode)
		allocate(pnode(mamax,nda))
		pnode(:,1:nd)=cpnode(:,1:nd)
		deallocate(cpnode)


        allocate(cgex(nd,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
        cgex(1:nd,1:ng)=gex(1:nd,1:ng)
        deallocate(gex);allocate(gex(nda,ng))
        gex=0
        gex(1:nd,1:ng)=cgex(1:nd,1:ng)
        deallocate(cgex)
  end if

  j=cels(cell)%nunodes !<<<<<<20.9.2013
  i=cels(cell)%nodela
  allocate(ccelnode(i))
  ccelnode(1:i)=cels(cell)%node(1:i)
  cels(cell)%nunodes=cels(cell)%nunodes+1
  cels(cell)%nodela=cels(cell)%nodela+1
  deallocate(cels(cell)%node)
  allocate(cels(cell)%node(cels(cell)%nodela))
  cels(cell)%node(1:i)=ccelnode(1:i)
  cels(cell)%node(j+1)=nd
  cels(cell)%node(i+1)=0 !>>>>20.9.2013
  nodeindex=0

call remade
if(nd>1) call neighbor_build

end subroutine pasti

!**********************************************************************************************************************
  
subroutine addcell(x,ycoord,z) !adds a basic epithelial cell to the desired plane 11.9.2013
real*8, intent(in) :: x,ycoord,z
integer plane
type(nod),allocatable :: cnodeo(:)
integer::i,radi,num,j,nodes,and
type(cel),allocatable :: ccels(:)
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi,xepi,yepi,zepi2,xstart,ystart,zstart,dy,d
integer kj

and=nd

xstart=x!+0.5d0
y=ycoord
de=0.5d0
di=0.5d0
vector=0.0d0
beta=2d0*pi/6d0
zepi=z
zepi2=z+di
alt=zepi
d=x!+de
dy=y

ncels=ncels+1
ncelsepi=ncelsepi+1
cels(ncels)=cels(ncels-1)
cels(ncels)%ctipus=1
cels(ncels)%nunodes=0
nodes=cels(ncels-1)%nodela
	!deallocate(cels(ncels)%node)
	!cels(ncels)%nodela=nodes
	!allocate(cels(ncels)%node(nodes))
	!cels(ncels)%node(:)=0

call addepi(xstart,dy,zepi)!central node down
node(nd)%altre=(nd+7)
node(nd)%icel=ncels
node(nd)%tipus=2
cels(ncels)%node(1)=nd 

call addepi(d+de,dy,zepi)!node 1
node(nd)%altre=(nd+7)
node(nd)%icel=ncels
node(nd)%tipus=2
!node1=node(nd)
cels(ncels)%node(2)=nd

angle=0
p1(1)=d;p1(2)=dy;p1(3)=alt
do kj=3,7 !nodes 3-7
angle=angle+beta
!angle=beta
p2(1)=d+de*dcos(angle);p2(2)=dy+de*dsin(angle);p2(3)=alt
p1=p2
call addepi(p1(1),p1(2),p1(3))
node(nd)%altre=(nd+7)
node(nd)%icel=ncels
node(nd)%tipus=2
cels(ncels)%node(kj)=nd
p2(1)=d+de*dcos(angle);p2(2)=dy+de*dsin(angle);p2(3)=alt
end do

alt=zepi2

call addepi(xstart,dy,zepi2)!central node up
node(nd)%altre=(nd-7)
node(nd)%icel=ncels
node(nd)%tipus=1
cels(ncels)%node(8)=nd

angle=0

call addepi(d+de,dy,zepi2)!node 8
node(nd)%altre=(nd-7)
node(nd)%icel=ncels
node(nd)%tipus=1
cels(ncels)%node(9)=nd

p1(1)=d;p1(2)=dy;p1(3)=alt
do kj=10,14 !nodes 10-14
angle=angle+beta
p2(1)=d+de*dcos(angle);p2(2)=dy+de*dsin(angle);p2(3)=alt
p1=p2
call addepi(p1(1),p1(2),p1(3))
node(nd)%altre=(nd-7)
node(nd)%icel=ncels
node(nd)%tipus=1
cels(ncels)%node(kj)=nd
p2(1)=d+de*dcos(angle);p2(2)=dy+de*dsin(angle);p2(3)=alt
end do !>>>>> 11.9.2013


do i=1,cels(ncels)%nunodes
  j=cels(ncels)%node(i)
  print*,j,"altre",node(j)%altre
end do

cels(ncels)%fase=0
 
!define the cell's centroid
jj=cels(ncels)%nunodes !to deal with if the previous cell has more than 14 nodes, since this basic type only creates 14 nodes in a cell.
do j=1,cels(ncels)%nunodes
   if (cels(ncels)%node(j) == 0) then
      jj=jj-1
   endif
end do
cels(ncels)%nunodes=jj   
 
		cels(ncels)%ctipus=1
		cels(ncels)%cex=0;cels(ncels)%cey=0;cels(ncels)%cez=0
		do j=1,cels(ncels)%nunodes
			k=cels(ncels)%node(j)
           		 if(node(k)%tipus==1)then
	 			cels(ncels)%cex=cels(ncels)%cex+node(ncels)%x
 				cels(ncels)%cey=cels(ncels)%cey+node(ncels)%y
				cels(ncels)%cez=cels(ncels)%cez+node(ncels)%z
			end if
		end do
		cels(ncels)%cex=2*cels(ncels)%cex/real(cels(ncels)%nunodes)
		cels(ncels)%cey=2*cels(ncels)%cey/real(cels(ncels)%nunodes)
		cels(ncels)%cez=2*cels(ncels)%cez/real(cels(ncels)%nunodes)



if(ncels+1>=ncals) then
	allocate(ccels(ncels))
	do i=1,ncels
		allocate(ccels(i)%node(cels(i)%nodela))
		ccels(i)%nunodes=cels(i)%nunodes
		ccels(i)%node(1:ccels(i)%nunodes)=cels(i)%node(1:cels(i)%nunodes)
		ccels(i)%nodela=cels(i)%nodela
		ccels(i)%cex=cels(i)%cex ; ccels(i)%cey=cels(i)%cey ; ccels(i)%cez=cels(i)%cez
		ccels(i)%polx=cels(i)%polx ; ccels(i)%poly=cels(i)%poly ; ccels(i)%polz=cels(i)%polz
!		print*,"i",cels(i)%nunodes,"nodes",cels(i)%node
!		print*,"i",ccels(i)%nunodes,"nodes",ccels(i)%node
	end do
	deallocate(cels)

	ncals=ncals+20


	allocate(cels(ncals))
	do i=1,ncels
		allocate(cels(i)%node(ccels(i)%nodela))
		cels(i)%node=0
		cels(i)%nunodes=0
		cels(i)%node(1:ccels(i)%nunodes)=ccels(i)%node(1:ccels(i)%nunodes)
		cels(i)%nunodes=ccels(i)%nunodes
		cels(i)%nodela=ccels(i)%nodela
		cels(i)%cex=ccels(i)%cex ; cels(i)%cey=ccels(i)%cey ; cels(i)%cez=ccels(i)%cez
		cels(i)%polx=ccels(i)%polx ; cels(i)%poly=ccels(i)%poly ; cels(i)%polz=ccels(i)%polz
	end do
	deallocate(ccels)
endif






call remade
if(nd>1) call neighbor_build

end subroutine addcell
      
!*************************************************************************************

subroutine update
integer:: ic,icc
type(nod),allocatable :: cnode(:),cpnode(:,:)
integer,allocatable::clist(:)
real*8,allocatable::cgex(:,:)

!call iniboxes
!ic=boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))!position in boxes
!  if(ic==0)then
!    !the cube was empty
!    boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))=nd
!    list(nd)=0
!  else		
!    !the cube was not empty so we search for the end																
!   do
!      icc=ic; ic=list(icc) ; if(ic==0)then ; list(icc)=nd ; list(nd)=0 ; exit ; end if
!    end do
!  end if

gex(nd,:)=gex(nd-1,:) 


if(nd+1>=nda)then	!let's enhance the node matrix and list matrix
		nda=nda+10
		allocate(cnode(nd))
		cnode=node
		deallocate(node)
		allocate(node(nda))
		node(1:nd)=cnode(1:nd)
		deallocate(cnode)


		allocate(clist(nd))
		clist=list
		deallocate(list)
		allocate(list(nda))
		list=0
		list(1:nd)=clist(1:nd)
		deallocate(clist)

		allocate(cpnode(mamax,nd))	!!>>Miquel 16-10-12
		cpnode=pnode
		deallocate(pnode)
		allocate(pnode(mamax,nda))
		pnode(:,1:nd)=cpnode(:,1:nd)
		deallocate(cpnode)


        allocate(cgex(nd,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
        cgex(1:nd,1:ng)=gex(1:nd,1:ng)
        deallocate(gex);allocate(gex(nda,ng))
        gex=0
        gex(1:nd,1:ng)=cgex(1:nd,1:ng)
        deallocate(cgex)
end if

call remade
if(nd>1) call neighbor_build

end subroutine update

!*************************************************************************************

subroutine addepi(x,y,z) !Adds an epithelial node
!logical,intent(in) :: altreswitch
real*8, intent(in) :: x,y,z
real*8::you,adh, rep, repcel,req,reqcel,reqs,ke,tor,stor,x1,y1,z1
integer::fase,tipus,icel,altre,marge,ic,icc
type(nod),allocatable :: cnode(:),cpnode(:,:)
integer,allocatable::clist(:),ccelnode(:)
real*8,allocatable::cgex(:,:),cagex(:,:)

nd=nd+1
ndepi=ndepi+1
node(nd)=node(nd-1)
node(nd)%x=x
node(nd)%y=y
node(nd)%z=z
node(nd)%tipus=1
node(nd)%altre=0
cell=ncels

!call iniboxes
!ic=boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))!position in boxes
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

gex(nd,:)=gex(nd-1,:) 


if(nd+1>=nda)then	!let's enhance the node matrix and list matrix
		nda=nda+10
		allocate(cnode(nd))
		cnode=node
		deallocate(node)
		allocate(node(nda))
		node(1:nd)=cnode(1:nd)
		deallocate(cnode)


		allocate(clist(nd))
		clist=list
		deallocate(list)
		allocate(list(nda))
		list=0
		list(1:nd)=clist(1:nd)
		deallocate(clist)

		allocate(cpnode(mamax,nd))	!!>>Miquel 16-10-12
		cpnode=pnode
		deallocate(pnode)
		allocate(pnode(mamax,nda))
		pnode(:,1:nd)=cpnode(:,1:nd)
		deallocate(cpnode)


        allocate(cgex(nd,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
        cgex(1:nd,1:ng)=gex(1:nd,1:ng)
        deallocate(gex);allocate(gex(nda,ng))
        gex=0
        gex(1:nd,1:ng)=cgex(1:nd,1:ng)
        deallocate(cgex)

        allocate(cagex(nd,ng))              !>>>>>>>>>>>>>>>>>>Miquel 3-6-13
        cagex(1:nd,1:ng)=agex(1:nd,1:ng)
        deallocate(agex);allocate(agex(nda,ng))
        agex=0
        agex(1:nd,1:ng)=cagex(1:nd,1:ng)
        deallocate(cagex)
end if
j=cels(cell)%nunodes !<<<<< 11.9.2013
i=cels(cell)%nodela
allocate(ccelnode(i))
ccelnode(1:i)=cels(cell)%node(1:i)
deallocate(cels(cell)%node)
cels(cell)%nunodes=j+1
cels(cell)%nodela=i+1
allocate(cels(cell)%node(i+1))
cels(cell)%node(1:i)=ccelnode(1:i)
cels(cell)%node(j+1)=nd
cels(cell)%node(i+1)=0 ! >>>>> 11.9.2013

call remade

end subroutine addepi

!*********************************************************************************************

subroutine selectcell(x,y,z) !Selects a cell nearest to the given coordinates
real*8, intent(in) :: x,y,z
call selectnode(x,y,z)
cellid = node(nodeindex)%icel
print*,"The selected cell is cell number ", cellid !<<<<<<<Tommi 17.9.2013
print*,"Tipus of the selected cell", cels(cellid)%ctipus
print*,"Fase of the selected cell",cels(cellid)%fase
k=cels(cellid)%nunodes
print*,"The cell has ",k," nodes"
print*,"which are:"
do i = 1,k
  print*,"Node ",i, "is", cels(cellid)%node(i), "and its altre", node(cels(cellid)%node(i))%altre
end do!>>>>>>>Tommi
end subroutine selectcell

!***********************************************************************************************
subroutine pastecell(x,y,z) !<<<<< 11.9.2013 ! Pastes a cell from selection to the chosen plane
real*8, intent(in) :: x,y,z
integer :: plane
integer :: pasteselection
integer::i,j,numbernodes,nodestart,tempnod,nodenext, testvar, nodetemp
logical::down
type(cel),allocatable :: ccels(:)
real*8             ::planedistvect, minimumvaluex,maximumvaluez,minimumvaluey,& 
                     minx,distance,posx,newx,newy,newz,vectorx,vectory,vectorz,posy,posz !The pasting of the cell starts from the node in upper left corner of the selected cell into the the selected coordinates
type(cel) :: copycell
!type(nod) :: tempnod
type :: nodelist
     real*8 :: xvalue
     real*8 :: yvalue
     real*8 :: zvalue
     integer :: nodeid
end type
type(nodelist), allocatable :: nodecopy(:),nodex(:), nodey(:), nodey1(:)&
                              ,nodey2(:), nodez(:),nodeslayer(:), cnode(:),ccnode(:), nodeorg(:) !comparelist
type :: vectors
     real*8 :: vectorx
     real*8 :: vectory
     real*8  :: vectorz
     integer :: nodeid
end type
type(vectors), allocatable :: positions(:)

if (cellid<=0) then
print*,"No cell selected!"
else
if (cels(cellid)%ctipus<2) then
ncelsepi=ncelsepi+1
else
ncelsmes=ncelsmes+1
endif

ncels=ncels+1
cels(ncels)=cels(cellid)

allocate(nodecopy(cels(cellid)%nunodes))
do i=1,cels(cellid)%nunodes
  tempnod=cels(cellid)%node(i)
  nodecopy(i)%xvalue=node(tempnod)%x
  nodecopy(i)%zvalue=node(tempnod)%z
  nodecopy(i)%yvalue=node(tempnod)%y
  nodecopy(i)%nodeid=tempnod
end do   

minimumvaluex=maxval(nodecopy%zvalue)
tempnod=cels(cellid)%node(1)
planedistvect=minimumvaluex-node(tempnod)%z
posy=y
posx=x
posz=z!-(planedistvect)
nodestart=cels(cellid)%node(1)

allocate(positions(cels(cellid)%nunodes-1))

do i=2,cels(cellid)%nunodes
nodenext=cels(cellid)%node(i)
vectorx=node(nodenext)%x-node(nodestart)%x
vectory=node(nodenext)%y-node(nodestart)%y
vectorz=node(nodenext)%z-node(nodestart)%z
positions(i-1)%vectorx=vectorx
positions(i-1)%vectory=vectory
positions(i-1)%vectorz=vectorz
positions(i-1)%nodeid=nodenext
end do

cels(ncels)%nunodes=0
 
if (node(nodestart)%tipus<3) then
call addepi(posx,posy,posz)
else
call addnode(posx,posy,posz)
endif
node(nd)=node(nodestart)
node(nd)%x=posx;node(nd)%y=posy;node(nd)%z=posz;
node(nd)%icel=ncels
cels(ncels)%node(1)=nd
kk=nodestart-node(nodestart)%altre
if (node(nodestart)%altre==0) then 
   node(nd)%altre=0
else
node(nd)%altre=nd-kk
endif

do i=1,size(positions)
newx=posx+positions(i)%vectorx
newy=posy+positions(i)%vectory
newz=posz+positions(i)%vectorz
nodenext=positions(i)%nodeid
if (node(nodenext)%tipus<3) then
call addepi(newx,newy,newz)
else
call addnode(newx,newy,newz)
endif
node(nd)=node(nodenext)
node(nd)%x=newx;node(nd)%y=newy;node(nd)%z=newz;
node(nd)%icel=ncels
cels(ncels)%node(i+1)=nd
kk=nodenext-node(nodenext)%altre
if (node(nodenext)%altre==0) then 
   node(nd)%altre=0
else
node(nd)%altre=nd-kk
endif
enddo

cels(ncels)%cex=0;cels(ncels)%cey=0;cels(ncels)%cez=0
		do j=1,cels(ncels)%nunodes
			k=cels(ncels)%node(j)
           		if(node(k)%tipus==1)then
	 			cels(ncels)%cex=cels(ncels)%cex+node(ncels)%x
 				cels(ncels)%cey=cels(ncels)%cey+node(ncels)%y
				cels(ncels)%cez=cels(ncels)%cez+node(ncels)%z
			end if
		end do
		cels(ncels)%cex=2*cels(ncels)%cex/real(cels(ncels)%nunodes)
		cels(ncels)%cey=2*cels(ncels)%cey/real(cels(ncels)%nunodes)
		cels(ncels)%cez=2*cels(ncels)%cez/real(cels(ncels)%nunodes)
if(ncels+1>=ncals) then
	allocate(ccels(ncels))
	do i=1,ncels
		allocate(ccels(i)%node(cels(i)%nodela))
		ccels(i)%nunodes=cels(i)%nunodes
		ccels(i)%node(1:ccels(i)%nunodes)=cels(i)%node(1:cels(i)%nunodes)
		ccels(i)%nodela=cels(i)%nodela
		ccels(i)%cex=cels(i)%cex ; ccels(i)%cey=cels(i)%cey ; ccels(i)%cez=cels(i)%cez
		ccels(i)%polx=cels(i)%polx ; ccels(i)%poly=cels(i)%poly ; ccels(i)%polz=cels(i)%polz
	end do
	deallocate(cels)

	ncals=ncals+20


	allocate(cels(ncals))
	do i=1,ncels
		allocate(cels(i)%node(ccels(i)%nodela))
		cels(i)%node=0
		cels(i)%nunodes=0
		cels(i)%node(1:ccels(i)%nunodes)=ccels(i)%node(1:ccels(i)%nunodes)
		cels(i)%nunodes=ccels(i)%nunodes
		cels(i)%nodela=ccels(i)%nodela
		cels(i)%cex=ccels(i)%cex ; cels(i)%cey=ccels(i)%cey ; cels(i)%cez=ccels(i)%cez
		cels(i)%polx=ccels(i)%polx ; cels(i)%poly=ccels(i)%poly ; cels(i)%polz=ccels(i)%polz
	end do
	deallocate(ccels)

end if
cellid=0
end if

call remade
if(nd>1) call neighbor_build

end subroutine pastecell !>>>>11.9.2013

!************************************************************************************************

subroutine deletenode !<<<<<20.9.2013
integer laltre,icel
!integer tipi, altre !>>> TT 15-6-2020

icel=node(nodeindex)%icel



call apoptosis(nodeindex)



if(cels(icel)%nunodes==0)then
  !print*,"ncels",ncels
  call eliminate_cell(icel)
  !print*,"ncels",ncels
end if

end subroutine

!*****************************************************************************************

subroutine deletecell(cell) !<<<< 11.9.2013 Deletes the selected cell, doesn't work yet
integer :: cell, nodes,ii,jj

 
 print *,cell,"cell"
 nodes=cels(cell)%nunodes
 if(cels(cell)%ctipus<3) nodes=nodes/2

 do ii=1,nodes
   jj=cels(cell)%node(1)
   call apoptosis(jj)
 end do 

 print *,nd,ncels
 !read(*,*)

 call eliminate_cell(cell)

if(nd>1) call neighbor_build

end subroutine


!subroutine updatecellnodes(cel)!subroutine that updates every cells node matrix, so that if a node is deleted, also cells node matrices are uptodate. Receives the index of the cell of deleted node as input
!integer, intent(in) :: cel
!integer :: cellindex
 

end module editor

