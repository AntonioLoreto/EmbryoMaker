!    GNOMO software (General Node Model)
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
module fitmo
  use general
  use neighboring
  use io
  !use conservative6 !pfh procrustes, not currently used
  !use complexity6 no functions currently used!
  use conservative_R ! module created by Renske, contains simple EMD and PF H measures of similarity to ancestor
  
  implicit none
  
  real*8 scale
  integer nbo
  integer, allocatable :: bo(:,:,:)  
  integer, allocatable :: bom(:,:,:,:)
  integer, allocatable :: outside (:)
  real*8,  allocatable :: boave(:,:,:,:),altre_boave(:,:,:,:)
  integer, allocatable :: nboave(:,:,:)

contains

!*********************************************************************************************************************
subroutine fit(ind,refit,val,pernofi2, tarpernofi) !!>> HC 30-6-2023 Save the EMD value
  integer ind
  real*8 refit,com,val,distfit !!>> HC 30-6-2023 Save the EMD value
  character*140 :: pernofi2, tarpernofi

  !print *, "going to calculate fitness"
 ! call readsnap(pernofi2)  !!>> HC 9-5-2023 Not needed, already done in fit.f90

  !print *, "readsnap in fit has run"
  if (allocated(outside)) deallocate(outside)
  allocate(outside(nd))


  call neighbor_build
  com=0.0
  
  call simple_EMD(1, com, pernofi2, tarpernofi)
  if (distfitscale==0.0d0 .or. distfitmag==0.0d0) then
    print *, "error in fitness: parameters distfit and distfitmag undefined"
    stop
  end if
  val=com   !!>> HC 30-6-2023 Save the EMD value
  distfit=distfitmag/(1.+(com**3)/(distfitscale**3)) !first order declining hill function
  !have to rescale the com from simple_EMD: this is the Euclidean distance, so should get the highest fitness when distance is 0. 
  !print *, "fitness is ", distfit
  refit=distfit

end subroutine


subroutine fit_old(ind,refit,pernofi2)
  integer :: ind,i,j,k,kk,kkk,jj
  real*8 :: refit,com,sumcom,d,dbas,dbasmin,dapi,dapimin,dd,facom
  character*140 :: pernofi2
!  real*8, dimension(nd,3) :: finalcloud4
  integer :: k2, kk2,  kkk2
  
 call readsnap(pernofi2)
 
 if (allocated(outside)) deallocate(outside)
 allocate(outside(nd))

  outside=0
  !no_epi_break: it checks if the epithelium is broken, it just checks if some
  ! epitelial node has so few neighbors that is means that the epithelium is broken
  ! it just check if more than 10% of the nodes are in an open side
  call neighbor_build

  kk=0 ; kk2=0
  kkk=0 ; kkk2=0
  do i=1,nd
    if (node(i)%tipus==1) then
      k=0
      kkk=kkk+1
      do j=1,nneigh(i)
        ii=neigh(i,j)
        if (node(ii)%tipus==1) then
          d=sqrt((node(i)%x-node(ii)%x)**2+(node(i)%y-node(ii)%y)**2+(node(i)%z-node(ii)%z)**2)
          if (d<(node(i)%add+node(ii)%add)*1.2d0) then !!>> HC 30-6-2020
            k=k+1
          end if
        end if
      end do
      if (k<5) then 
        kk=kk+1 !; print*,"nd",i,"number neigh 1",k
      end if
    end if
    if (node(i)%tipus==2) then
      k2=0
      kkk2=kkk2+1
      do j=1,nneigh(i)
        ii=neigh(i,j)
        if (node(ii)%tipus==2) then
          d=sqrt((node(i)%x-node(ii)%x)**2+(node(i)%y-node(ii)%y)**2+(node(i)%z-node(ii)%z)**2)
          if (d<(node(i)%add+node(ii)%add)*1.2d0) then !!>> HC 30-6-2020
            k2=k2+1
          end if
        end if
      end do
      if (k2<5) then 
        kk2=kk2+1 !; print*,"nd",i,"number neigh 2",k2
      end if
    end if
  end do
   
  print *,"open nodes type 1",kk,"nd*0.2",kkk*0.2,"nd",kkk
  if (kk>kkk*0.2.or.kkk==0) then ; print *,"open epithelium 1" ; refit=0.0 ; return ; end if
  
  print *,"open nodes type 2",kk2,"nd*0.2",kkk2*0.2,"nd",kkk2
  if (kk2>kkk2*0.2.or.kkk==0) then ; print *,"open epithelium 2" ; refit=0.0 ; return ; end if
  !***** pfh
 k=0
 kk=0
 do i=1,nd   
   if (node(i)%tipus==1) then   
    ii=node(i)%altre  
    dd=sqrt((node(i)%x-node(ii)%x)**2+(node(i)%y-node(ii)%y)**2+(node(i)%z-node(ii)%z)**2) !dd:distance between apical and basal nodes of 1 epith cell
      do j=1,nneigh(ii)	
          iii=neigh(ii,j)  
          if(i==iii)cycle
     !     if(node(iii)%tipus.ne.2)then;k=k+1;cycle;endif
          kk=kk+1
          d=sqrt((node(i)%x-node(iii)%x)**2+(node(i)%y-node(iii)%y)**2+(node(i)%z-node(iii)%z)**2) !distance between apical node and neigh basal cells
          if (d.le.dd*0.75) then
             k=k+1
          end if
      end do
   end if 
 end do  
a=real(k)/nd
print*,"FU?",a !;stop
if (a>0.05) then ; print *,"fucked up epithelium" ; refit=0.0 ; return ; end if    
!average distance to neigh from apical and basal side more or less the same
kkk=0
do i=1,nd
   k=0;kk=0;d=0;dd=0
   do j=1,nneigh(i)
     iii=neigh(i,j)  
     k=k+1	
     if(i==iii)cycle
     dd=dd+sqrt((node(i)%x-node(iii)%x)**2+(node(i)%y-node(iii)%y)**2+(node(i)%z-node(iii)%z)**2)
   enddo      
   ii=node(i)%altre 
   do j=1,nneigh(ii)	
     iii=neigh(ii,j) 
     kk=kk+1 
     if(iii==i)cycle 
     d=d+sqrt((node(i)%x-node(iii)%x)**2+(node(i)%y-node(iii)%y)**2+(node(i)%z-node(iii)%z)**2) 
   enddo 
   dd=dd/k;d=d/kk
   if(dd/d>1.1.or.dd/d<0.1)kkk=kkk+1   
enddo
print*,"kkk",kkk,"nd",nd,"limit",nd/5
if(kkk>nd/5)then;refit=0d0;print*,"CAAOOOS";refit=0.0; return ;endif 


  com=0.0
  
!  call total(ind,com,10)
!  call node_to_node(ind,com)
call functional_complexity2(ind,com,pernofi2) 
!**call functional_complexity(ind,com,pernofi2)  
!call volumenfit(ind,com,pernofi2)


!*! print *,"complextiat",com
!**call surfacearea(ind,com)
sumcom=com
!  com=maxval(node(:nd)%x)-minval(node(:nd)%x)
!print *,com,"com"
!  com=0.0
!  call ord4_pco_per_logpco_forcells(ind,com)
!print *,refit,com
!  sumcom=sumcom+com
  com=1.0d0
!print *,sumcom,"pre other control"
  facom=1
  call other_control(ind,facom)     ! this is to give fitness 0 to the fusion of cells in a syncitial ball
!print *,com,"coms"
!  call check_radial_sim(ind,com)
!**print *,com,"fau",facom,"facom"
!print*,"simetry",com;stop
  com=sumcom*com*facom  !yes here it has to be multiplied


  refit=com
 print *,com,sumcom,"final similarity"

end subroutine

!*********************************************************************************************************************
subroutine volumenfit(ind,com,pernofi3) !change in conservative6 if shared volumen or emd is wanted
integer ind
real*8, dimension(nd,3) :: finalcloud4
real*8 com
 character*140 :: pernofi3
!print*,"1",node(1)%x
!  do i=1,nd
!   finalcloud4(i,1)=node(i)%x
!   finalcloud4(i,2)=node(i)%y
!   finalcloud4(i,3)=node(i)%z
!  end do 
!print*,"2",finalcloud4(1,1)
!open(45,file="emd.dat")
  print*,pernofi3,"*/*/*/*/*/*/*"
  finalcloud4=0
 ! call conservative(com,finalcloud4,pernofi3)
 ! call conservative(com,pernofi3)
  com=com*1d5
  com=1/(1+com) 
  print*,"com",com
  if(com.ge.0.99)print*,"OPTIMUM ACHIEVED!!!, ind",ind!;stop;endif

! close(45)
!  e=2.718281828459045
!  com=com*com*com*com
!  com=-1*com/2
!  com=e**com
!  print*,"com",com
!  if(com.ge.0.99)then;print*,"OPTIMUM ACHIEVED!!!, ind",ind;stop;endif
!stop

end subroutine volumenfit
!*********************************************************************************************************************

subroutine surfacearea(ind,com)
real*8 aa,surfacesphere,surfaceembryo,com,extrea
integer ind
!integer surfaceembryo
!call extrem!no, solo epi, no ecm
! aa=0d0
surfaceembryo=0
surfacesphere=0
 do i=1,nd
   if(node(i)%tipus==2) then
      aa=sqrt((node(i)%x)**2+(node(i)%y)**2+(node(i)%z)**2)
      if (aa>extrea) then ; extrea=aa ; end if 
      surfaceembryo=node(i)%eqd*node(i)%eqd*pi+surfaceembryo !!>> HC 30-6-2020
   endif
  enddo
surfacesphere=4*pi*extrea*extrea
print*,"surfacesphere",surfacesphere
print*,"surfaceembryo",surfaceembryo
!com=surfaceembryo*surfaceembryo-surfacesphere
com=surfaceembryo-surfacesphere
com=com*com
print*,"complejidad",com

end subroutine surfacearea
!*********************************************************************************************************************
subroutine check_radial_sim(ind,com)

integer ind
integer imaxd,imind
integer qua(nd)
real*8 noux(nd),nouy(nd),nouz(nd)
real*8 cnoux(nd),cnouy(nd),cnouz(nd)
real*8 com
real*8 a,b,c,mix,miy,miz
real*8 covmat(3,3),eigenvec(3,3),eigenva(3)
real*8 prodesc
real*8 mind,mindd,angle,x,y,z,sumind
integer nparti

nparti=100 !number of different angles checked

! find centroid
mix=0.0d0 ; miy=0.0d0 ; miz=0.0d0

do i=1,nd
  mix=mix+node(i)%x  
  miy=miy+node(i)%y  
  miz=miz+node(i)%z  
end do

a=1.0d0/nd

mix=mix*a ; miy=miy*a ; miz=miz*a

node(:nd)%x=node(:nd)%x-mix
node(:nd)%y=node(:nd)%y-miy
node(:nd)%z=node(:nd)%z-miz

noux=node(:nd)%x
nouy=node(:nd)%y
nouz=node(:nd)%z

 cnoux=noux ; cnouy=nouy ; cnouz=nouz

 c=0
mind=1.0d10
sumind=0.0d0
do ii=1,nparti !we change the angles: we rotate: how many times we rotate (total rotation is always 2*pi)

  b=2*3.1416/100.0d0
  angle=b*ii

  do i=1,nd
    noux(i)=cnoux(i)*cos(angle)-cnouy(i)*sin(angle)
    nouy(i)=cnouy(i)*cos(angle)-cnoux(i)*sin(angle)
  end do

  ! EMD distance between between nou and cnou for that angle
  a=0.0d0
  do i=1,nd
    mindd=1.0d10
    do j=1,nd
      if (i==j) cycle
      d=sqrt((noux(i)-cnoux(j))**2+(nouy(i)-cnouy(j))**2+(nouz(i)-cnouz(j))**2)
      if (d<mindd) then ; mindd=d ;kk=j;end if
    end do
    a=a+mindd
  end do

  do i=1,nd
    mindd=1.0d10
    do j=1,nd
      if (i==j) cycle
      d=sqrt((cnoux(i)-noux(j))**2+(cnouy(i)-nouy(j))**2+(cnouz(i)-nouz(j))**2)
      if (d<mindd) then ; mindd=d ;kk=j;end if
    end do
    a=a+mindd
  end do
  a=a/real(2*nd**2)
  sumind=sumind+a
  if (mind>a) then
    mind=a
  end if
end do


print *,mind,"mind"
if (mind>8.0d-4) then ;  !this is a reasonable value
  com=0   !not radial enough
else
  com=1
end if

end subroutine check_radial_sim

!*********************************************************************************************************************
!*********************************************************************************************************************
subroutine fit_radial_sim(ind,com)

integer ind,axis
integer imaxd,imind
integer qua(nd)
real*8 noux(nd),nouy(nd),nouz(nd)
real*8 cnoux(nd),cnouy(nd),cnouz(nd)
real*8 com
real*8 a,b,c,mix,miy,miz
real*8 covmat(3,3),eigenvec(3,3),eigenva(3)
real*8 prodesc
real*8 mind,mindd,angle,x,y,z,sumind
integer nparti

nparti=100 !number of different angles checked
axis=nparti
! find centroid
mix=0.0d0 ; miy=0.0d0 ; miz=0.0d0

do i=1,nd
  mix=mix+node(i)%x  
  miy=miy+node(i)%y  
  miz=miz+node(i)%z  
end do

a=1.0d0/nd

mix=mix*a ; miy=miy*a ; miz=miz*a

node(:nd)%x=node(:nd)%x-mix
node(:nd)%y=node(:nd)%y-miy
node(:nd)%z=node(:nd)%z-miz



noux=node(:nd)%x
nouy=node(:nd)%y
nouz=node(:nd)%z


 cnoux=noux ; cnouy=nouy ; cnouz=nouz

 c=0
mind=1.0d10
sumind=0.0d0
do ii=1,nparti !we change the angles: we rotate: how many times we rotate (total rotation is always 2*pi)

  b=2*3.1416/100.0d0
  angle=b*ii

  do i=1,nd
    noux(i)=cnoux(i)*cos(angle)-cnouy(i)*sin(angle)
    nouy(i)=cnouy(i)*cos(angle)-cnoux(i)*sin(angle)
  end do

  ! EMD distance between between nou and cnou for that angle
  a=0.0d0
  do i=1,nd
    mindd=1.0d10
    do j=1,nd
      if (i==j) cycle
      d=sqrt((noux(i)-cnoux(j))**2+(nouy(i)-cnouy(j))**2+(nouz(i)-cnouz(j))**2)
      if (d<mindd) then ; mindd=d ;kk=j;end if
    end do
!print *,i,kk,mindd
    a=a+mindd
  end do

  do i=1,nd
    mindd=1.0d10
    do j=1,nd
      if (i==j) cycle
      d=sqrt((cnoux(i)-noux(j))**2+(cnouy(i)-nouy(j))**2+(cnouz(i)-nouz(j))**2)
      if (d<mindd) then ; mindd=d ;kk=j;end if
    end do
    a=a+mindd

  end do
  a=a/real(2*nd**2)
!axis of simmetry
if (a>2.0d-3) axis=axis-1
  sumind=sumind+a 
  if (mind>a) then
    mind=a 
  end if
end do


print *,mind,"mind"
if (mind>8.0d-4 .or. axis==0) then ;  !this is a reasonable value
  com=0   !not radial enough
else
  com=real(1/axis)
end if

end subroutine fit_radial_sim

!*********************************************************************************************************************
subroutine check_radial_simd(ind,com)

integer ind
integer imaxd,imind
integer qua(nd)
real*8 noux(nd),nouy(nd),nouz(nd)
real*8 cnoux(nd),cnouy(nd),cnouz(nd)
real*8 com
real*8 a,b,c,mix,miy,miz
real*8 covmat(3,3),eigenvec(3,3),eigenva(3)
real*8 prodesc
real*8 mind,maxd,angle,x,y,z

! find centroid
mix=0.0d0 ; miy=0.0d0 ; miz=0.0d0

do i=1,nd
  mix=mix+node(i)%x  
  miy=miy+node(i)%y  
  miz=miz+node(i)%z  
end do

a=1.0d0/nd

mix=mix*a ; miy=miy*a ; miz=miz*a

node(:nd)%x=node(:nd)%x-mix
node(:nd)%y=node(:nd)%y-miy
node(:nd)%z=node(:nd)%z-miz

! first we calculate the covariance matrix
! x vs y

!cov x x
b=0.0d0
do i=1,nd
 b=b+node(i)%x*node(i)%x
end do
b=b*a
covmat(1,1)=b

!cov x y
b=0.0d0
do i=1,nd
 b=b+node(i)%x*node(i)%y
end do
b=b*a
covmat(1,2)=b

!cov x z
b=0.0d0
do i=1,nd
 b=b+node(i)%x*node(i)%z
end do
b=b*a
covmat(1,3)=b

!cov y y
b=0.0d0
do i=1,nd
 b=b+node(i)%y*node(i)%y
end do
b=b*a
covmat(2,2)=b

!cov y z
b=0.0d0
do i=1,nd
 b=b+node(i)%y*node(i)%z
end do
b=b*a
covmat(2,3)=b

!cov z z
b=0.0d0
do i=1,nd
 b=b+node(i)%z*node(i)%z
end do
b=b*a
covmat(3,3)=b

covmat(2,1)=covmat(1,2)
covmat(3,1)=covmat(1,3)
covmat(3,2)=covmat(2,3)

print *,covmat

!call DSYEVJ3(covmat, eigenvec,eigenva)  

print *,eigenva,"eigenva"

print *,eigenvec(1:3,1),"eigenvec 1"
print *,eigenvec(1:3,2),"eigenvec 2"
print *,eigenvec(1:3,3),"eigenvec 3"

!now we check the distance in the direction of the first eigenvector
b=1.0d0/sqrt(eigenvec(1,1)**2+eigenvec(2,1)**2+eigenvec(3,1)**2)
do i=1,nd
  prodesc=(node(i)%x*eigenvec(1,1)+node(i)%y*eigenvec(2,1)+node(i)%z*eigenvec(3,1))*b
  if (prodesc>maxd) then ; maxd=prodesc ; imaxd=i ; end if
  if (prodesc<mind) then ; mind=prodesc ; imind=i ; end if
!  pro(i)=prodesc
!print *,i,prodesc,"prodesc"
end do

!print *,mind,maxd,"mind maxd",imaxd,imind

!we divide the whole thing in four quadrats with the first eigenvector as the axis
!the four quadrants need to be at least roughly equal for some kind of asymmetry to be possible
!for simplicity we assume a real radial symmetry (e.g.no pentaradial symmetry or alike)

c=0

do ii=1,100 !we change the angles

  b=2*3.1416/100.0d0
  angle=b*(ii-1)

  do i=1,nd
    noux(i)=cnoux(i)*cos(angle)-cnouy(i)*sin(angle)
    nouy(i)=cnouy(i)*cos(angle)-cnoux(i)*sin(angle)
  end do

!to determine the quadrant we simply make the dot product porjection between a position and the two other eigenvecs
b=1.0d0/sqrt(eigenvec(1,1)**2+eigenvec(2,1)**2+eigenvec(3,1)**2)
do i=1,nd
  nouz(i)=(node(i)%x*eigenvec(1,1)+node(i)%y*eigenvec(2,1)+node(i)%z*eigenvec(3,1))*b
end do
b=1.0d0/sqrt(eigenvec(1,2)**2+eigenvec(2,2)**2+eigenvec(3,2)**2)
do i=1,nd
  noux(i)=(node(i)%x*eigenvec(1,2)+node(i)%y*eigenvec(2,2)+node(i)%z*eigenvec(3,2))*b
end do
b=1.0d0/sqrt(eigenvec(1,3)**2+eigenvec(2,3)**2+eigenvec(3,3)**2)
do i=1,nd
  nouy(i)=(node(i)%x*eigenvec(1,3)+node(i)%y*eigenvec(2,3)+node(i)%z*eigenvec(3,3))*b
end do

!now we simply establish the quadrant based on those coordinates
do i=1,nd
  if (noux(i)>0.and.nouy(i)>0) qua(i)=1
  if (noux(i)>0.and.nouy(i)<0) qua(i)=2
  if (noux(i)<0.and.nouy(i)<0) qua(i)=3
  if (noux(i)<0.and.nouy(i)>0) qua(i)=4
node(i)%rep=qua(i)
end do

cnoux=noux ; cnouy=nouy ; cnouz=nouz


! EMD distance between between quadrant 1 and 2
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==1) then
    do j=1,nd
      if (qua(j)==2) then
        d=sqrt((-node(i)%x-node(j)%x)**2+(-node(i)%y-node(j)%y)**2+(-node(i)%z-node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 1 2"

! EMD distance between between quadrant 2 and 1
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==2) then
    do j=1,nd
      if (qua(j)==1) then
        d=sqrt((-node(i)%x-node(j)%x)**2+(-node(i)%y-node(j)%y)**2+(-node(i)%z-node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 1 2"


! EMD distance between between quadrant 1 and 3
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==1) then
    do j=1,nd
      if (qua(j)==3) then
        d=sqrt((-node(i)%x+node(j)%x)**2+(-node(i)%y+node(j)%y)**2+(-node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 1 3"

! EMD distance between between quadrant 3 and 1
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==3) then
    do j=1,nd
      if (qua(j)==1) then
        d=sqrt((-node(i)%x+node(j)%x)**2+(-node(i)%y+node(j)%y)**2+(-node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 3 1"


! EMD distance between between quadrant 1 and 4
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==1) then
    do j=1,nd
      if (qua(j)==4) then
        d=sqrt((node(i)%x+node(j)%x)**2+(node(i)%y+node(j)%y)**2+(node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 1 4"

! EMD distance between between quadrant 4 and 1
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==4) then
    do j=1,nd
      if (qua(j)==1) then
        d=sqrt((node(i)%x+node(j)%x)**2+(node(i)%y+node(j)%y)**2+(node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 4 1"


! EMD distance between between quadrant 2 and 3
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==2) then
    do j=1,nd
      if (qua(j)==3) then
        d=sqrt((node(i)%x+node(j)%x)**2+(node(i)%y+node(j)%y)**2+(node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 2 3"

! EMD distance between between quadrant 3 and 2
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==3) then
    do j=1,nd
      if (qua(j)==2) then
        d=sqrt((node(i)%x+node(j)%x)**2+(node(i)%y+node(j)%y)**2+(node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 2 3"

! EMD distance between between quadrant 2 and 4
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==2) then
    do j=1,nd
      if (qua(j)==4) then
        d=sqrt((-node(i)%x+node(j)%x)**2+(-node(i)%y+node(j)%y)**2+(-node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 2 4"

! EMD distance between between quadrant 4 and 2
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==4) then
    do j=1,nd
      if (qua(j)==2) then
        d=sqrt((-node(i)%x+node(j)%x)**2+(-node(i)%y+node(j)%y)**2+(-node(i)%z+node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 2 4"


! EMD distance between between quadrant 3 and 4
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==3) then
    do j=1,nd
      if (qua(j)==4) then
        d=sqrt((-node(i)%x-node(j)%x)**2+(-node(i)%y-node(j)%y)**2+(-node(i)%z-node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 3 4"

! EMD distance between between quadrant 3 and 4
a=0.0d0
b=0.0d0
do i=1,nd
  maxd=0.0d0
  if (qua(i)==4) then
    do j=1,nd
      if (qua(j)==3) then
        d=sqrt((-node(i)%x-node(j)%x)**2+(-node(i)%y-node(j)%y)**2+(-node(i)%z-node(j)%z)**2)
        if (d>maxd) maxd=d
      end if
    end do
    b=b+1
    a=a+maxd
  end if
end do
b=1.0d0/b
c=c+a*b
!print *,a*b,"qua 3 4"
end do


print *,c,"c"

end subroutine


!*********************************************************************************************************************
subroutine total_old(ind,com,npart) ! for each node, or one node, it calculates the distance and angle to all other nodes then make categories
                          ! of distance and calculates de variance of those angles per category and sums them
integer ind
real*8  com
integer npart

integer i,j,k,ii,jj,kk,iii,jjj,kkk
real*8  distaun(nd)
integer idistaun(nd)
real*8 ix,iy,iz  ! coordinates of i
real*8 cx,cy,cz  ! components of the cylinder vector
real*8 ccx,ccy,ccz  ! components to neighbor
real*8 dotpro,angle,moda,modb,a,b,c,bcont,cc
real*8 nangles(nd)
real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
integer dd
real*8 ma,mc,mb,tam_boxes

call neighbor_build

distaun=0.0d0
com=0.0d0

do i=1,nd

  if (node(i)%tipus/=1) cycle

  distaun=0.0d0
  idistaun=0
  dd=0

!do i=1,nd
!  if (node(i)%tipus==1) then
!    goto 14
!  end if
!end do

14 if (node(i)%tipus<3) then

    nangles=0.0d0
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    ii=node(i)%altre
    cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
    moda=sqrt(cx**2+cy**2+cz**2)

    do j=1,nd
      if (i==j) cycle
      if (node(j)%tipus/=node(i)%tipus) then ; dd=dd+1 ;cycle ; endif

      !distance
      distaun(j)=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+(node(i)%z-node(j)%z)**2)

      !angle
      ccx=node(j)%x-ix ; ccy=node(j)%y-iy ; ccz=node(j)%z-iz
      modb=sqrt(ccx**2+ccy**2+ccz**2)
      dotpro=cx*ccx+cy*ccy+cz*ccz
      a=dotpro/(moda*modb)
12    if (a>1) then ! this is because acos function domain is only from -1 to 1
        a=2-a       ! this is effectively a reflection
      else
        if (a<-1) then
          a=-2-a     ! this is also a reflection
        end if
      end if
      if (abs(a)>1) goto 12
      angle=acos(a)-pi*0.5d0
      nangles(j)=angle
    end do

    !now we want to order them
    call ordenarepe(distaun,idistaun,nd)

    tam_boxes=(nd-dd)/real(npart)
    coma=0
    do j=1,npart
      k=dd+tam_boxes*(j-1)+1
      kk=dd+tam_boxes*j
      a=0.0d0
      cc=0
      do jj=k,kk
        jjj=idistaun(jj)
        if (node(jjj)%tipus==node(i)%tipus) then
          cc=cc+1
          a=a+nangles(jjj)
        end if
      end do
      d=1.0d0/real(cc)
      a=a*d
      b=0.0d0
      do jj=k,kk
        jjj=idistaun(jj)
        if (node(jjj)%tipus==node(i)%tipus) then
          b=b+(nangles(jjj)-a)**2
        end if
      end do
      d=1.0d0/real(cc)
      b=b*d
      coma=(coma+b) !*d	!pfh, *d, = coma*d+b, d no * b (arriba, 213)
!print *,i,j,b,"variance among the ones at a given distance category",cc
    end do
  end if
  coma=coma/real(npart)
!  print *,i,coma,"coma",carg
  bcont=bcont+1
  com=com+coma
end do

com=com/bcont

print *,com,"com final"
end subroutine

!*********************************************************************************************************************
subroutine total(ind,com,npart) ! for each node, or one node, it calculates the distance and angle to all other nodes then make categories
                          ! of distance and calculates de variance of those angles per category and sums them
integer ind
real*8  com
integer npart

integer i,j,k,ii,jj,kk,iii,jjj,kkk
real*8  distaun(nd)
integer idistaun(nd)
real*8 ix,iy,iz  ! coordinates of i
real*8 cx,cy,cz  ! components of the cylinder vector
real*8 ccx,ccy,ccz  ! components to neighbor
real*8 dotpro,angle,moda,modb,a,b,c,bcont,cc
real*8 nangles(nd)
real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
integer dd
real*8 ma,mc,mb,tam_boxes,lomax,vain,vainmu

call neighbor_build

distaun=0.0d0
com=0.0d0

do i=1,nd

  if (node(i)%tipus/=1) cycle

  distaun=0.0d0
  idistaun=0
  dd=0

!do i=1,nd
!  if (node(i)%tipus==1) then
!    goto 14
!  end if
!end do

14 if (node(i)%tipus<3) then

    nangles=0.0d0
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    ii=node(i)%altre
    cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
    moda=sqrt(cx**2+cy**2+cz**2)

    do j=1,nd
      if (i==j) cycle
      if (node(j)%tipus/=node(i)%tipus) then ; dd=dd+1 ;cycle ; endif

      !distance
      distaun(j)=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+(node(i)%z-node(j)%z)**2)

      !angle
      ccx=node(j)%x-ix ; ccy=node(j)%y-iy ; ccz=node(j)%z-iz
      modb=sqrt(ccx**2+ccy**2+ccz**2)
      dotpro=cx*ccx+cy*ccy+cz*ccz
      a=dotpro/(moda*modb)
12    if (a>1) then ! this is because acos function domain is only from -1 to 1
        a=2-a       ! this is effectively a reflection
      else
        if (a<-1) then
          a=-2-a     ! this is also a reflection
        end if
      end if
      if (abs(a)>1) goto 12
      angle=acos(a)-pi*0.5d0
      nangles(j)=angle
    end do

    !now we want to order them
    call ordenarepe(distaun,idistaun,nd)

    ii=nd-dd+2
!    lomax=distaun(idistaun(nd))
!    lomax=maxval(node(:nd)%add) !*npart   !!>> HC 30-6-2020
lomax=0d0
do j=1,nd
  if(node(j)%tipus>2) cycle
  if(node(j)%add>lomax) lomax=node(j)%add  !!>> HC 30-6-2020
end do



    coma=0
    do j=3,npart          !we start at 5 so that small noise between neighbors nodes does not affect
      vain=real(j)*lomax
      vainmu=real(j-1)*lomax
      !k=dd+tam_boxes*(j-1)+1
      !kk=dd+tam_boxes*j
      a=0.0d0
      cc=0
!print *,vain,vainmu,lomax/npart
      do jj=1,nd     !here we calculate mean
        jjj=idistaun(jj)
        if (distaun(jjj)>vainmu.and.distaun(jjj)<vain) then
          if (node(jjj)%tipus==node(i)%tipus) then
            cc=cc+1
            a=a+nangles(jjj)
          end if
        end if
      end do
      d=1.0d0/real(cc)
      a=a*d
      b=0.0d0       !here we calcuate deviation
      do jj=1,nd
        jjj=idistaun(jj)
        if (distaun(jjj)>vainmu.and.distaun(jjj)<vain) then
          if (node(jjj)%tipus==node(i)%tipus) then
            b=b+(nangles(jjj)-a)**2
          end if
        end if
      end do
!print *,j,cc,"cc"
      if (cc>0.0) then
        d=1.0d0/real(cc)
        b=b*d
        coma=(coma+b) !*d	!pfh, *d, = coma*d+b, d no * b (arriba, 213)
      end if
!print *,i,j,b,"variance among the ones at a given distance category",cc
!read(*,*)
    end do
  end if
!  coma=coma/real(npart)
!  print *,i,coma,"coma",carg
  bcont=bcont+1
  com=com+coma
end do

com=com/bcont

!print *,com,"com final"
end subroutine


!*******************************************************************************************************
!*********************************************************************************************************************
subroutine total_tail(ind,com,npart) ! for each node, or one node, it calculates the distance and angle to all other nodes then make categories
                          ! of distance and calculates de variance of those angles per category and sums them
integer ind
real*8  com
integer npart

integer i,j,k,ii,jj,kk,iii,jjj,kkk
real*8  distaun(nd)
integer idistaun(nd)
real*8 ix,iy,iz  ! coordinates of i
real*8 cx,cy,cz  ! components of the cylinder vector
real*8 ccx,ccy,ccz  ! components to neighbor
real*8 dotpro,angle,moda,modb,a,b,c,bcont,cc
real*8 nangles(nd)
real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
integer dd
real*8 ma,mc,mb,tam_boxes,lomax,vain,vainmu

call neighbor_build

distaun=0.0d0
com=0.0d0

do i=1,nd

  if (node(i)%tipus/=1) cycle
  !ojo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !gex(:,:)  ! expression level: x is node id, y is the gene 
   ! if(gex(i,1)==1)cycle  !gen 1 en cola

  distaun=0.0d0
  idistaun=0
  dd=0

!do i=1,nd
!  if (node(i)%tipus==1) then
!    goto 14
!  end if
!end do

14 if (node(i)%tipus<3) then


   

    nangles=0.0d0
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    ii=node(i)%altre
    cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
    moda=sqrt(cx**2+cy**2+cz**2)

    do j=1,nd
      if (i==j) cycle
      if (node(j)%tipus/=node(i)%tipus) then ; dd=dd+1 ;cycle ; endif

      !distance
      distaun(j)=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+(node(i)%z-node(j)%z)**2)

      !angle
      ccx=node(j)%x-ix ; ccy=node(j)%y-iy ; ccz=node(j)%z-iz
      modb=sqrt(ccx**2+ccy**2+ccz**2)
      dotpro=cx*ccx+cy*ccy+cz*ccz
      a=dotpro/(moda*modb)
12    if (a>1) then ! this is because acos function domain is only from -1 to 1
        a=2-a       ! this is effectively a reflection
      else
        if (a<-1) then
          a=-2-a     ! this is also a reflection
        end if
      end if
      if (abs(a)>1) goto 12
      angle=acos(a)-pi*0.5d0
      nangles(j)=angle
    end do

    !now we want to order them
    call ordenarepe(distaun,idistaun,nd)

    ii=nd-dd+2
!    lomax=distaun(idistaun(nd))
!    lomax=maxval(node(:nd)%add) !*npart  !!>> HC 30-6-2020
lomax=0d0
do j=1,nd
  if(node(j)%tipus>2) cycle
  if(node(j)%add>lomax) lomax=node(j)%add !!>> HC 30-6-2020
end do



    coma=0
    do j=3,npart          !we start at 5 so that small noise between neighbors nodes does not affect
      vain=real(j)*lomax
      vainmu=real(j-1)*lomax
      !k=dd+tam_boxes*(j-1)+1
      !kk=dd+tam_boxes*j
      a=0.0d0
      cc=0
!print *,vain,vainmu,lomax/npart
      do jj=1,nd     !here we calculate mean
        jjj=idistaun(jj)
        if (distaun(jjj)>vainmu.and.distaun(jjj)<vain) then
          if (node(jjj)%tipus==node(i)%tipus) then
            cc=cc+1
            a=a+nangles(jjj)
          end if
        end if
      end do
      d=1.0d0/real(cc)
      a=a*d
      b=0.0d0       !here we calcuate deviation
      do jj=1,nd
        jjj=idistaun(jj)
        if (distaun(jjj)>vainmu.and.distaun(jjj)<vain) then
          if (node(jjj)%tipus==node(i)%tipus) then
            b=b+(nangles(jjj)-a)**2
          end if
        end if
      end do
!print *,j,cc,"cc"
      if (cc>0.0) then
        d=1.0d0/real(cc)
        b=b*d
        coma=(coma+b) !*d	!pfh, *d, = coma*d+b, d no * b (arriba, 213)
      end if
!print *,i,j,b,"variance among the ones at a given distance category",cc
!read(*,*)
    end do
  end if
!  coma=coma/real(npart)
!  print *,i,coma,"coma",carg
  bcont=bcont+1
  com=com+coma
end do

com=com/bcont

!print *,com,"com final"
end subroutine


!*******************************************************************************************************
!*********************************************************************************************************************

!*********************************************************************************************************************
subroutine totalpinta! for each node, or one node, it calculates the distance and angle to all other nodes then make categories
                          ! of distance and calculates de variance of those angles per category and sums them
integer ind
real*8  com
integer npart

integer i,j,k,ii,jj,kk,iii,jjj,kkk
real*8  distaun(nd)
integer idistaun(nd)
real*8 ix,iy,iz  ! coordinates of i
real*8 cx,cy,cz  ! components of the cylinder vector
real*8 ccx,ccy,ccz  ! components to neighbor
real*8 dotpro,angle,moda,modb,a,b,c,bcont,cc
real*8 nangles(nd)
real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
integer dd
real*8 ma,mc,mb,tam_boxes
com=0.0
npart=10.0
call neighbor_build

distaun=0.0d0
com=0.0d0

do i=1,nd

  if (node(i)%tipus/=1) cycle

  distaun=0.0d0
  idistaun=0
  dd=0

!do i=1,nd

!  if (node(i)%tipus==1) then
!    goto 14
!  end if
!end do

14 if (node(i)%tipus<3) then

    nangles=0.0d0
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    ii=node(i)%altre
    cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
    moda=sqrt(cx**2+cy**2+cz**2)

    do j=1,nd
      if (i==j) cycle
      if (node(j)%tipus/=node(i)%tipus) then ; dd=dd+1 ;cycle ; endif

      !distance
      distaun(j)=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+(node(i)%z-node(j)%z)**2)

      !angle
      ccx=node(j)%x-ix ; ccy=node(j)%y-iy ; ccz=node(j)%z-iz
      modb=sqrt(ccx**2+ccy**2+ccz**2)
      dotpro=cx*ccx+cy*ccy+cz*ccz
      a=dotpro/(moda*modb)
12    if (a>1) then ! this is because acos function domain is only from -1 to 1
        a=2-a       ! this is effectively a reflection
      else
        if (a<-1) then
          a=-2-a     ! this is also a reflection
        end if
      end if
      if (abs(a)>1) goto 12
      angle=acos(a)-pi*0.5d0
      nangles(j)=angle
    end do

    !now we want to order them
    call ordenarepe(distaun,idistaun,nd)

    tam_boxes=(nd-dd)/real(npart)
    coma=0
    do j=1,npart
      k=dd+tam_boxes*(j-1)+1
      kk=dd+tam_boxes*j
      a=0.0d0
      cc=0
      do jj=k,kk
        jjj=idistaun(jj)
        if (node(jjj)%tipus==node(i)%tipus) then
          cc=cc+1
          a=a+nangles(jjj)
        end if
      end do
      d=1.0d0/real(cc)
      a=a*d
      b=0.0d0
      do jj=k,kk
        jjj=idistaun(jj)
        if (node(jjj)%tipus==node(i)%tipus) then
          b=b+(nangles(jjj)-a)**2
        end if
      end do
      d=1.0d0/real(cc)
      b=b*d
      coma=coma+b
print *,i,j,b,"variance among the ones at a given distance category",cc
!read(*,*)
    end do
  end if
  coma=coma/real(npart)
  print *,i,coma,"coma",carg
  bcont=bcont+1
  com=com+coma
end do

com=com/bcont

print *,com,"com final"
end subroutine totalpinta

!*******************************************************************************************************

subroutine ordenarepe(ma,mt,rang)
  integer rang
  real*8 ma(rang)
  integer mt(rang)
  integer i,j,k
  real*8 a
    mt=0
el: do i=1,rang
      a=ma(i) ; k=1
      do j=1,rang ; if (a>ma(j)) k=k+1 ; end do 
      do j=k,rang ; if (mt(j)==0) then ; mt(j)=i ; cycle el ; end if ; end do
    end do el 
end subroutine ordenarepe

!*********************************************************************************************************************
subroutine curvature_ordred(ind,com)
  integer ind
  real*8 com

  integer i,j,k,ii,jj,kk,iii,jjj,kkk
  real*8 ix,iy,iz  ! coordinates of i
  real*8 cx,cy,cz  ! components of the cylinder vector
  real*8 ccx,ccy,ccz  ! components to neighbor
  real*8 dotpro,angle,moda,modb,a,b,c,bcont,cc
  real*8 nangles(nd)
  integer iik
  real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  real*8 ma,mc,mb,cl,tam_boxes

  call neighbor_build

  com=0.0
  bcont=0.0d0

  do i=1,nd
    if (node(i)%tipus==2) then

      nangles=0.0d0
      cc=0.0d0
      ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
      ii=node(i)%altre
      cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
      moda=sqrt(cx**2+cy**2+cz**2)

      do j=1,nneigh(i)
        k=neigh(i,j)
        if (node(k)%tipus/=2) cycle
        ccx=node(k)%x-ix ; ccy=node(k)%y-iy ; ccz=node(k)%z-iz
        modb=sqrt(ccx**2+ccy**2+ccz**2)
        dotpro=cx*ccx+cy*ccy+cz*ccz
        a=dotpro/(moda*modb)
!print *,a,dotpro,moda,modb
12      if (a>1) then ! this is because acos function domain is only from -1 to 1
          a=2-a       ! this is effectively a reflection
        else
          if (a<-1) then
            a=-2-a     ! this is also a reflection
          end if
        end if
        if (abs(a)>1) goto 12
!print *,a,"a"
        cc=cc+1
        angle=acos(a)-pi*0.5d0
!        com=com+abs(angle)
        nangles(int(cc))=angle
!print *,i,k,angle*360/(2*pi),"angle"        
      end do

      a=0.0d0
      b=0.0d0
      do j=1,int(cc)
        c=nangles(j)
        a=a+c
        b=b+c**2
      end do
      c=cc
      c=1d0/c
      a=a*c
      b=b*c
      com=com+b-a**2
      bcont=bcont+1
    end if
  end do
  com=com/bcont
print *,com,"scale 0"

  coma=0.0d0
  ! second order neighbors
  do i=1,nd
    if (node(i)%tipus==2) then

      nangles=0.0d0
      cc=0.0d0

      ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
      ii=node(i)%altre
      cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
      moda=sqrt(cx**2+cy**2+cz**2)

      do j=1,nneigh(i)
        k=neigh(i,j)
        if (node(k)%tipus/=2) cycle

er:     do jj=1,nneigh(k)
          kk=neigh(k,jj)
          if (node(kk)%tipus/=2) cycle
          if (kk==i) cycle
          do jjj=1,nneigh(i)
            if (neigh(i,jjj)==kk) cycle er
          end do

          ccx=node(kk)%x-ix ; ccy=node(kk)%y-iy ; ccz=node(kk)%z-iz
          modb=sqrt(ccx**2+ccy**2+ccz**2)
          dotpro=cx*ccx+cy*ccy+cz*ccz
          a=dotpro/(moda*modb)
!print *,a,dotpro,moda,modb
13        if (a>1) then ! this is because acos function domain is only from -1 to 1
            a=2-a       ! this is effectively a reflection
          else
            if (a<-1) then
              a=-2-a     ! this is also a reflection
            end if
          end if
          if (abs(a)>1) goto 13
!print *,a,"a"
          cc=cc+1
          angle=acos(a)-pi*0.5d0
!          coma=coma+abs(angle)
          nangles(int(cc))=angle
!print *,i,k,angle*360/(2*pi),"angle"        
        end do er
      end do

      a=0.0d0
      b=0.0d0
      do j=1,int(cc) !since we do not take all the neighbors of the neighbor
        c=nangles(j)
        a=a+c
        b=b+c**2
      end do
      c=cc
      c=1d0/c
      a=a*c
      b=b*c
      coma=coma+b-a**2
      bcont=bcont+1
    end if
  end do
  coma=coma/bcont
  com=com+coma
print *,coma,"scale 1"


end subroutine

!*********************************************************************************************************************
subroutine curvature_simple(ind,com)
  integer ind
  real*8 com

  integer i,j,k,ii,jj,kk,iii,jjj,kkk
  real*8 ix,iy,iz  ! coordinates of i
  real*8 cx,cy,cz  ! components of the cylinder vector
  real*8 ccx,ccy,ccz  ! components to neighbor
  real*8 dotpro,angle,moda,modb,a,b,c,bcont
  real*8 nangles(nd)
  integer iik
  real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  real*8 ma,mc,mb,cl,tam_boxes

  call neighbor_build

  com=0.0
  bcont=0.0d0

  do i=1,nd
    if (node(i)%tipus==2) then

      nangles=0.0d0

      ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
      ii=node(i)%altre
      cx=node(ii)%x-ix ; cy=node(ii)%y-iy ; cz=node(ii)%z-iz
      moda=sqrt(cx**2+cy**2+cz**2)

      do j=1,nneigh(i)
        k=neigh(i,j)
        if (node(k)%tipus/=2) cycle
        ccx=node(k)%x-ix ; ccy=node(k)%y-iy ; ccz=node(k)%z-iz
        modb=sqrt(ccx**2+ccy**2+ccz**2)
        dotpro=cx*ccx+cy*ccy+cz*ccz
        a=dotpro/(moda*modb)
print *,a,dotpro,moda,modb
12      if (a>1) then ! this is because acos function domain is only from -1 to 1
          a=2-a       ! this is effectively a reflection
        else
          if (a<-1) then
            a=-2-a     ! this is also a reflection
          end if
        end if
        if (abs(a)>1) goto 12
print *,a,"a"
        angle=acos(a)-pi*0.5d0
        com=com+abs(angle)
        nangles(j)=angle
print *,i,k,angle*360/(2*pi),"angle"        
      end do
      a=0.0d0
      b=0.0d0
      do j=1,nneigh(i)
        c=nangles(j)
        a=a+c
        b=b+c**2
      end do
      c=nneigh(i)
      c=1d0/c
      a=a*c
      b=b*c
      com=com+b-a**2
      bcont=bcont+1
    end if
  end do
  com=com/bcont
print *,com,"scale 0"

end subroutine

!*********************************************************************************************************************
subroutine curvature(ind,com)
  integer ind
  real*8 com

  integer i,j,k,ii,jj,kk,iii,jjj,kkk
  real*8 ix,iy,iz  ! coordinates of i
  real*8 cx,cy,cz  ! components of the cylinder vector
  real*8 ccx,ccy,ccz  ! components to neighbor
  real*8 dotpro,angle,moda,modb,a,b,c,bcont
  real*8 nangles(nd)
  integer iik
  real*8 sc,v,coma,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  real*8 ma,mc,mb,cl,tam_boxes
  integer inic,fi,nboo
  integer fr, ir, k2, kk2, kkk2
  
  a=0.0d0
  d=0.0d0
  do i=1,nd
!    if (outside(i)==1) cycle 
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    do j=1,nd
      ax=node(j)%x-ix ; ay=node(j)%y-iy ; az=node(j)%z-iz
      d=sqrt(ax**2+ay**2+az**2)
      if (d>a) then
        a=d
        bx=ax ; by=ay ; bz=az
        ii=i  ; jj=j
      end if
    end do
  end do

  maxd=a
  if (maxd<2.0d0) maxd=2.0d0  !THIS IS TO FORBIDE REALLY SMALL

  ! we calculate the angle between the vector of the longest distance and the x axis
  if (bz==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+bz**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
!   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  if (bx==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+by**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0.0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
!   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  !now we reposition the embryo so that it is all in positive coordinates
  ! ma mb mc
  ma=1d10 ; mb=1d10 ; mc=1d10
  do i=1,nd
!   if (outside(i)==1) cycle
    if (ma>node(i)%x) ma=node(i)%x
    if (mb>node(i)%y) mb=node(i)%y
    if (mc>node(i)%z) mc=node(i)%z
  end do

  do i=1,nd
!   if (outside(i)==1) cycle
    node(i)%x=node(i)%x-ma    
    node(i)%y=node(i)%y-mb    
    node(i)%z=node(i)%z-mc    
  end do

  !now we find the longest dimension assuming all start at 0,0,0
  a=0.0
  do i=1,nd
!   if (outside(i)==1) cycle
    if (a<node(i)%x) a=node(i)%x
    if (a<node(i)%y) a=node(i)%y
    if (a<node(i)%z) a=node(i)%z
  end do
  maxd=a

  !now we find the longest dimension assuming all start at 0,0,0
  maxdd=1.0d10
  do i=1,nd
!   if (outside(i)==1) cycle
   if (node(i)%add<maxdd) maxdd=node(i)%add !!>> HC 30-6-2020
  end do
  !maxdd=a/real(nd)

  print *,maxd,maxdd,"hu"

123 nbo=7 !nbo=2*maxd/maxdd ! the number of boxes
   print *,nbo,"nbo"
  if (allocated(bo)) deallocate(bo)
!  nboo=int(maxd/maxdd)
  nboo=nbo
  if (allocated(bo)) deallocate(bo)
  if (allocated(boave)) deallocate(boave)
  if (allocated(altre_boave)) deallocate(altre_boave)
  if (allocated(nboave)) deallocate(nboave)
  allocate(bo(-1:nboo+1,-1:nboo+1,-1:nboo+1))
  allocate(boave(-1:nboo+1,-1:nboo+1,-1:nboo+1,3))
  allocate(altre_boave(-1:nboo+1,-1:nboo+1,-1:nboo+1,3))
  allocate(nboave(-1:nboo+1,-1:nboo+1,-1:nboo+1))

  boave=0
  nboave=0
  altre_boave=0

  cl=0
  com=0.0d0
print *,maxd,"maxd"  
  do ir=1,nbo
    a=nbo-1
!    tam_boxes=maxdd+(maxd-maxdd)*(ir-1)/a
    tam_boxes=maxd/real(ir)
    bo=0
    fr=0
    coma=0
    v=0.0d0

    sc=1d0/tam_boxes

    do i=1,nd
!!ACHTUNG      if (outside(i)==1) cycle
      ii=int(node(i)%x*sc);jj=int(node(i)%y*sc);kk=int(node(i)%z*sc)
      if (node(i)%x==maxd) ii=ir-1
      if (node(i)%y==maxd) jj=ir-1
      if (node(i)%z==maxd) kk=ir-1
      if (node(i)%tipus==1) then
        bo(ii,jj,kk)=1
        nboave(ii,jj,kk)=nboave(ii,jj,kk)+1
        boave(ii,jj,kk,1)=node(i)%x
        boave(ii,jj,kk,2)=node(i)%y
        boave(ii,jj,kk,3)=node(i)%z
        iii=node(i)%altre
        altre_boave(ii,jj,kk,1)=node(iii)%x
        altre_boave(ii,jj,kk,2)=node(iii)%y
        altre_boave(ii,jj,kk,3)=node(iii)%z
      end if
    end do

    v=0
    do i=0,nboo-1
      do j=0,nboo-1
        do k=0,nboo-1
          if (bo(i,j,k)==1) then
            v=v+1
            a=nboave(i,j,k)
            a=1d0/a
            boave(i,j,k,:)=boave(i,j,k,:)*a  !average positions per box
            altre_boave(i,j,k,:)=altre_boave(i,j,k,:)*a
            altre_boave(i,j,k,:)=altre_boave(i,j,k,:)-boave(i,j,k,:)  ! the cylinder vector components for the box
          end if
        end do
      end do
    end do

    ! now we look which boxes are next to each other and the angles.
    coma=0.0d0
    do i=0,nboo-1
      do j=0,nboo-1
        do k=0,nboo-1
          if (bo(i,j,k)==1) then ! here we check in a box
            cc=0
            do ii=-1,1
              do jj=-1,1
                do kk=-1,1
                  if (ii==0.and.jj==0.and.kk==0) cycle
                  if (bo(i+ii,j+jj,k+kk)==1) then
                    call angleb(i,j,k,i+ii,j+jj,k+kk,angle) !THERE SHOULD BE A DISTANCE FILTER
                    cc=cc+1
!print *,ir,"ir",i,j,k,ii,jj,kk,cc,angle
                    nangles(int(cc))=angle
                  end if
                end do
              end do
            end do
            a=0.0d0
            b=0.0d0
            do jj=1,int(cc)
              c=nangles(jj)
              a=a+c
              b=b+c**2
            end do
            if (cc>0) then
              c=1d0/cc
              a=a*c
              b=b*c
              coma=coma+b-a**2
!print *,ir,coma,cc,"coma"
!read(*,*)
              bcont=bcont+1
            end if
          end if
        end do
      end do
    end do
!print *,coma,coma/bcont,bcont,"yuio"
    com=com+coma/bcont   !fitness is the sum of the complexity at each level
print *,ir,com,bcont,"ttt"
  end do  
  print *,com,"com"
end subroutine

!********************************************************************************************************************

subroutine angleb(cox,coy,coz,coox,cooy,cooz,angle)
  integer cox,coy,coz,coox,cooy,cooz
  integer i,j,k,ii,jj,kk,iii,jjj,kkk
  real*8  ccx,ccy,ccz,modb,moda,cx,cy,cz,a,b,c,cc,angle, dotpro

  cx=boave(cox,coy,coz,1)
  cy=boave(cox,coy,coz,2)
  cz=boave(cox,coy,coz,3)
  ccx=boave(coox,cooy,cooz,1)-cx
  ccy=boave(coox,cooy,cooz,2)-cy
  ccz=boave(coox,cooy,cooz,3)-cz
  moda=sqrt(altre_boave(cox,coy,coz,1)**2+altre_boave(cox,coy,coz,2)**2+altre_boave(cox,coy,coz,3)**2)
  modb=sqrt(ccx**2+ccy**2+ccz**2)
  dotpro=altre_boave(cox,coy,coz,1)*ccx+altre_boave(cox,coy,coz,2)*ccy+altre_boave(cox,coy,coz,2)*ccz
  a=dotpro/(moda*modb)
!print *,a,dotpro,moda,modb
12 if (a>1) then ! this is because acos function domain is only from -1 to 1
     a=2-a       ! this is effectively a reflection
   else
     if (a<-1) then
       a=-2-a     ! this is also a reflection
     end if
   end if
   if (abs(a)>1) goto 12
   angle=acos(a)-pi*0.5d0
end subroutine 

!*********************************************************************************************************************
subroutine surface_volume(ind,com)
  integer ind,i,j,k,ii,jj,kk,ir,iii,kkk! IMPORTANT, the minimal box size is the minimal da in the embryo and then it goes until all the embryo is included
  integer iik
  real*8 com,sc,v,a,b,coma,ix,iy,iz,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  integer up,do,an,po,le,ri,ce
  real*8 ma,mc,mb,cl,tam_boxes
  integer inic,fi,nboo
  integer, allocatable :: bo(:,:,:)
  real*8 fr(0:1,0:1,0:1,0:1,0:1,0:1)

  ! we calculculate the complexity as the join information based on the probability of each kind of neighborhood at radium one box

  ! it is independent of size because we split each embryo in the same number of boxes 

  ! we have to rotate everything so that it gets properly aligned to make the boxes along the longest axis of the embryo
  ! we first calculate the distance between the most distant nodes in the embryo
  a=0.0d0
  d=0.0d0
  do i=1,nd
    if (outside(i)==1) cycle 
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    do j=1,nd
      ax=node(j)%x-ix ; ay=node(j)%y-iy ; az=node(j)%z-iz
      d=sqrt(ax**2+ay**2+az**2)
      if (d>a) then
        a=d
        bx=ax ; by=ay ; bz=az
        ii=i  ; jj=j
      end if
    end do
  end do

  maxd=a
  if (maxd<2.0d0) maxd=2.0d0  !THIS IS TO FORBIDE REALLY SMALL


  ! we calculate the angle between the vector of the longest distance and the x axis
  if (bz==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+bz**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
!   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  if (bx==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+by**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0.0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
!   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  !now we reposition the embryo so that it is all in positive coordinates
  ! ma mb mc
  ma=1d10 ; mb=1d10 ; mc=1d10
  do i=1,nd
!   if (outside(i)==1) cycle
    if (ma>node(i)%x) ma=node(i)%x
    if (mb>node(i)%y) mb=node(i)%y
    if (mc>node(i)%z) mc=node(i)%z
  end do

  do i=1,nd
!   if (outside(i)==1) cycle
    node(i)%x=node(i)%x-ma    
    node(i)%y=node(i)%y-mb    
    node(i)%z=node(i)%z-mc    
  end do

  !now we find the longest dimension assuming all start at 0,0,0
  a=0.0
  do i=1,nd
!   if (outside(i)==1) cycle
    if (a<node(i)%x) a=node(i)%x
    if (a<node(i)%y) a=node(i)%y
    if (a<node(i)%z) a=node(i)%z
  end do
  maxd=a

  !now we find the longest dimension assuming all start at 0,0,0
  maxdd=1.0d10
  do i=1,nd
!   if (outside(i)==1) cycle
   if (node(i)%add<maxdd) maxdd=node(i)%add  !!>> HC 30-6-2020
  end do
  !maxdd=a/real(nd)

  print *,maxd,maxdd,"hu"

123 nbo=300 !nbo=2*maxd/maxdd ! the number of boxes
   print *,nbo,"nbo"
  if (allocated(bo)) deallocate(bo)
!  nboo=int(maxd/maxdd)
  nboo=nbo
  allocate(bo(-1:nboo+1,-1:nboo+1,-1:nboo+1))

  cl=0
  com=0.0d0
print *,maxd,"maxd"  
  !nbo=10
  do ir=1,nbo
    a=nbo-1
!    tam_boxes=maxdd+(maxd-maxdd)*(ir-1)/a
    tam_boxes=maxd/real(ir)
    bo=0
    fr=0
    coma=0
    v=0.0d0

    sc=1d0/tam_boxes

    do i=1,nd
!!ACHTUNG      if (outside(i)==1) cycle
      ii=int(node(i)%x*sc);jj=int(node(i)%y*sc);kk=int(node(i)%z*sc)
!print *,node(i)%x,node(i)%y,node(i)%z,int(node(i)%x*sc),int(node(i)%y*sc),int(node(i)%z*sc),maxd,maxdd,sc,tam_boxes
!read(*,*)
      if (node(i)%x==maxd) ii=ir-1
      if (node(i)%y==maxd) jj=ir-1
      if (node(i)%z==maxd) kk=ir-1
      if (node(i)%tipus==1) then
        bo(ii,jj,kk)=1
      end if
    end do

    v=0
    do i=0,nboo-1
      do j=0,nboo-1
        do k=0,nboo-1
          if (bo(i,j,k)==1) then
            v=v+1
          end if
        end do
      end do
    end do

    ! now the volume
    kkk=0
    iik=0
    do i=0,nboo-1
      do j=0,nboo-1
        fi=0
        inic=-1
        do k=0,nboo-1
          if (bo(i,j,k)>0) then
            if (inic==-1) then
              inic=k
!print *,inic,k,"ini"
!read(*,*)
            else
              fi=k
!print *,fi,k,"fi"
!read(*,*)
            end if
          end if
        end do
!print *,inic,fi,"ppp"
        if (inic>=fi) then
          iik=iik+1
        else
          if (inic>-1) then
            iik=iik+fi-inic+1
!print *,i,j,fi-inic+1
          end if
        end if
!        print *,inic,fi,iik,"rty"
      end do
    end do
    
    !surface area is:
    a=iik
    if (a>0) then
      coma=v/a
    else
      coma=0.0
    end if
print *,ir,v,a,coma,"err"

  com=com+coma   !fitness is the sum of the complexity at each level
  end do  

  cl=nbo
  com=com !/cl

111  print *,com,nbo,cl,maxd,maxdd,"fitness"

end subroutine

!*********************************************************************************************************************

subroutine pco_per_logpco(ind,com)   ! measures complexity as the mutual information between embryo and non-embryo boxes for different scales
  integer ind,i,j,k,ii,jj,kk,ir   ! NOT NORMALIZED BY SIZE, but with a minimal size being the minval
  real*8 com,scale,sc,v,a,b,coma
  integer, parameter :: nboo=10
!  integer, allocatable :: bo(:,:,:)
  real*8 fr(0:1,0:1)

  call extrem 

  ! this is the minimal scale 
  scale=0.0
  do i=1,nd
    if (node(i)%add>scale) then ; print *,i,node(i)%add ; scale=node(i)%add ; end if !!>> HC 30-6-2020
  end do
!  scale=maxval(node(:nd)%add) !!>> HC 30-6-2020
!print *,scale,"scale"
  nbo=extre/scale+2

  if (allocated(bo)) deallocate(bo)
  allocate(bo(-nbo:nbo,-nbo:nbo,-nbo:nbo))
  if (allocated(bom)) deallocate(bom)
  allocate(bom(-nbo:nbo,-nbo:nbo,-nbo:nbo,nbo))


  com=0.0d0

  do ir=1,nboo
    bo=0
    fr=0
    coma=0
    v=0.0d0
    sc=1.0d0/(scale*ir)

    do i=1,nd
      ii=nint(node(i)%x*sc);jj=nint(node(i)%y*sc);kk=nint(node(i)%z*sc)
      bo(ii,jj,kk)=1
    end do

    do i=-nbo+1,nbo-1
      do j=-nbo+1,nbo-1
        do k=-nbo+1,nbo-1
          v=v+bo(i,j,k)
        end do
      end do
    end do
!    print *,ir,v,v/nbo**3,"volume"

    ! now the mutual information
    b=0.0d0
    do i=-nbo+1,nbo-1
      do j=-nbo+1,nbo-1
        do k=-nbo+1,nbo-1
          kk=bo(i,j,k)
          kkk=bo(i+1,j,k)
          fr(kk,kkk)=fr(kk,kkk)+1
          kkk=bo(i,j+1,k)
          fr(kk,kkk)=fr(kk,kkk)+1
          kkk=bo(i,j,k+1)
          fr(kk,kkk)=fr(kk,kkk)+1
          b=b+3
        end do
      end do
    end do

    a=1.0d0/b
    fr=fr*a
    if (fr(0,0)/=0) coma=coma-fr(0,0)*log(fr(0,0))
    if (fr(1,0)+fr(0,1)/=0) coma=coma-(fr(1,0)+fr(0,1))*log(fr(1,0)+fr(0,1))
    if (fr(1,1)/=0) coma=coma-fr(1,1)*log(fr(1,1))

!    print *,fr(0,0),fr(0,1)+fr(1,0),fr(1,1),"frs"

    com=com+coma   !fitness is the sum of the complexity at each level

!print *,ir,coma,scale,"nivell com"
    if (v<4) exit
  end do



  print *,com,"fitness"

end subroutine

!****************************************************************************************************

subroutine ord4_pco_per_logpco_ns(ind,com)   ! it is like ord4_pco_per_logpco but its invariable with size
                                    ! measures complexity as the mutual information between embryo and non-embryo boxes for different scales
  integer ind,i,j,k,ii,jj,kk,ir,iii ! IMPORTANT, the minimal box size is the minimal da in the embryo and then it goes until all the embryo is included
  real*8 com,sc,v,a,b,coma,ix,iy,iz,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  integer up,do,an,po,le,ri,ce
  real*8 ma,mb,mc
  integer, parameter :: nboo=10
  integer, allocatable :: bo(:,:,:)
  real*8 fr(0:1,0:1,0:1,0:1,0:1,0:1)

  ! we calculculate the complexity as the join information based on the probability of each kind of neighborhood at radium one box

  ! it is independent of size because we split each embryo in the same number of boxes 

  ! we have to rotate everything so that it gets properly aligned to make the boxes along the longest axis of the embryo
  ! we first calculate the distance between the most distant nodes in the embryo
  a=0.0d0
  d=0.0d0
  do i=1,nd
    if (outside(i)==1) cycle 
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    do j=1,nd
      ax=node(j)%x-ix ; ay=node(j)%y-iy ; az=node(j)%z-iz
      d=sqrt(ax**2+ay**2+az**2)
      if (d>a) then
        a=d
        bx=ax ; by=ay ; bz=az
        ii=i  ; jj=j
      end if
    end do
  end do

  maxd=a
  if (maxd<2.0d0) maxd=2.0d0  !THIS IS TO FORBIDE REALLY SMALL


  ! we calculate the angle between the vector of the longest distance and the x axis
  if (bz==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+bz**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  if (bx==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+by**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0.0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
   if (outside(i)==1) cycle
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  !now we reposition the embryo so that it is all in positive coordinates
  ! ma mb mc
  ma=1d10 ; mb=1d10 ; mc=1d10
  do i=1,nd
   if (outside(i)==1) cycle
    if (ma>node(i)%x) ma=node(i)%x
    if (mb>node(i)%y) mb=node(i)%y
    if (mc>node(i)%z) mc=node(i)%z
  end do

  do i=1,nd
   if (outside(i)==1) cycle
    node(i)%x=node(i)%x-ma    
    node(i)%y=node(i)%y-mb    
    node(i)%z=node(i)%z-mc    
  end do

  !now we find the longest dimension assuming all start at 0,0,0
  a=0.0
  do i=1,nd
   if (outside(i)==1) cycle
    if (a<node(i)%x) a=node(i)%x
    if (a<node(i)%y) a=node(i)%y
    if (a<node(i)%z) a=node(i)%z
  end do
  maxd=a

  !now we find the longest dimension assuming all start at 0,0,0
  maxdd=1.0d10
  do i=1,nd
   if (outside(i)==1) cycle
   if (node(i)%add<maxdd) maxdd=node(i)%add !!>> HC 30-6-2020
  end do
  !maxdd=a/real(nd)

  print *,maxd,maxdd,"hu"

123  nbo=51 !nbo=2*maxd/maxdd ! the number of boxes
   print *,nbo,"nbo"
!  allocate(bo(-1:nbo,-1:nbo,-1:nbo))
  if (allocated(bo)) deallocate(bo)
  allocate(bo(-1:nbo,-1:nbo,-1:nbo))
  if (allocated(bom)) deallocate(bom)
  allocate(bom(-1:nbo,-1:nbo,-1:nbo,nbo))

  com=0.0d0
  
  !nbo=10
  do ir=1,nbo-1
!   do ir=3,3  
!  do ir=1,1
    bo=0
    fr=0
    coma=0
    v=0.0d0
!   sc=real(ir)/real(nbo)!maxd
    sc=real(ir)/maxd
!print *,ir,sc,"ir sc",maxd,maxdd,nbo

    if (sc<maxdd) cycle ! this is to avoid high complexity when nodes are unrealistically close
    do i=1,nd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ACHTUNG      if (outside(i)==1) cycle
!      ii=int(node(i)%x*sc+0.5d0);jj=int(node(i)%y*sc+0.5d0);kk=int(node(i)%z*sc+0.5d0)
      ii=int(node(i)%x*sc);jj=int(node(i)%y*sc);kk=int(node(i)%z*sc)
      if (node(i)%x==maxd) ii=ir-1
      if (node(i)%y==maxd) jj=ir-1
      if (node(i)%z==maxd) kk=ir-1
!print *,node(i)%x,node(i)%y,node(i)%z,maxd
!print *,ir,ii,jj,kk,"ii jj kk",sc
!print *,i,ir,maxd,sc,node(i)%x*sc,node(i)%y*sc,node(i)%z*sc
!read(*,*)
      bo(ii,jj,kk)=1
    end do

    bom(:,:,:,ir)=bo

    do i=0,nbo
      do j=0,nbo
        do k=0,nbo
          v=v+bo(i,j,k)
        end do
      end do
    end do
    kk=0
    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          if (bo(i,j,k)>0) then
!            print *,i,j,k,"kok"
            kk=kk+1
          end if
        end do
      end do
    end do

print *,kk,"kk"

!    print *,ir,v,v/nbo**3,nbo**3,"volum"

    ! now the mutual information
    b=0.0d0
    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          ce=bo(i,j,k)
          if (ce==1) then ! otherwise we are counting totally empty areas and this produces a bias based on the boxes and 
                                           ! the shape of the embryo
            ri=bo(i+1,j,k)
            le=bo(i-1,j,k)
            an=bo(i,j+1,k)
            po=bo(i,j-1,k)
            up=bo(i,j,k+1)
            do=bo(i,j,k-1)
            fr(ri,le,an,po,up,do)=fr(ri,le,an,po,up,do)+1
!print *,ri,le,an,po,up,do,fr(ri,le,an,po,up,do)
            b=b+1
          end if
        end do
      end do
    end do

    a=1.0d0/b
    fr=fr*a
!print *,a,"a"
      do j=0,1
        do k=0,1
          do ii=0,1
            do jj=0,1
              do kk=0,1
                do iii=0,1
                  if (fr(j,k,ii,jj,kk,iii)>0.0d0) then
                    coma=coma-fr(j,k,ii,jj,kk,iii)*log(fr(j,k,ii,jj,kk,iii))
!c=fr(j,k,ii,jj,kk,iii)
!cc=log(fr(j,k,ii,jj,kk,iii))
!print *,ir,j,k,ii,jj,kk,iii,c/a,c,cc,c*cc
!print *,ir,j,k,ii,jj,kk,iii,"aaa" !,fr(j,k,ii,jj,kk,iii),nd,"ww"
                    if (nint(fr(j,k,ii,jj,kk,iii)/a)==nd) then ; com=com+coma ;goto 111 ; end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      com=com+coma   !fitness is the sum of the complexity at each level
print *,ir,coma,"coma",v,nbo,"nbo",maxd,"maxd",com,"com"
print *,v,nd,"v nd"
    if (v>=nd) exit
  end do

111  print *,com,nbo,maxd,maxdd,"fitness"

end subroutine


!****************************************************************************************************

subroutine ord4_pco_per_logpco(ind,com)   ! measures complexity as the mutual information between embryo and non-embryo boxes for different scales
  integer ind,i,j,k,ii,jj,kk,ir,iii       ! NOT NORMALIZED BY SIZE, but with a minimal size being the minval
  real*8 com,sc,v,a,b,coma,ix,iy,iz,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd,maxdd
  integer up,do,an,po,le,ri,ce
  real*8 ma,mb,mc
  integer, parameter :: nboo=10
  integer, allocatable :: bo(:,:,:)
  real*8 fr(0:1,0:1,0:1,0:1,0:1,0:1)

  ! we calculculate the complexity as the join information based on the probability of each kind of neighborhood at radium one box

  ! it is independent of size because we split each embryo in the same number of boxes 

  ! we have to rotate everything so that it gets properly aligned to make the boxes along the longest axis of the embryo
  ! we first calculate the distance between the most distant nodes in the embryo
  a=0.0d0
  d=0.0d0
  do i=1,nd
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    do j=1,nd
      ax=node(j)%x-ix ; ay=node(j)%y-iy ; az=node(j)%z-iz
      d=sqrt(ax**2+ay**2+az**2)
      if (d>a) then
        a=d
        bx=ax ; by=ay ; bz=az
        ii=i  ; jj=j
      end if
    end do
  end do

  maxd=a
  if (maxd<2.0d0) maxd=2.0d0  !THIS IS TO FORBIDE REALLY SMALL


  ! we calculate the angle between the vector of the longest distance and the x axis
  if (bz==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+bz**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  if (bx==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+by**2)
    angleu=acos(abs(bx)/b)
  end if
angleu=0.0 !ACHTUNG
  !now we rotate the whole embryo for this angle 
  do i=1,nd
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  !now we reposition the embryo so that it is all in positive coordinates
  ! ma mb mc
  ma=1d10 ; mb=1d10 ; mc=1d10
  do i=1,nd
    if (ma>node(i)%x) ma=node(i)%x
    if (mb>node(i)%y) mb=node(i)%y
    if (mc>node(i)%z) mc=node(i)%z
  end do

  do i=1,nd
    node(i)%x=node(i)%x-ma    
    node(i)%y=node(i)%y-mb    
    node(i)%z=node(i)%z-mc    
  end do

  !now we find the longest dimension assuming all start at 0,0,0
  a=0.0
  do i=1,nd
    if (a<node(i)%x) a=node(i)%x
    if (a<node(i)%y) a=node(i)%y
    if (a<node(i)%z) a=node(i)%z
  end do
  maxd=a

  !now we find the longest dimension assuming all start at 0,0,0
  a=0.0d0
  do i=1,nd
    a=a+node(i)%add !!>> HC 30-6-2020
  end do
  maxdd=a/real(nd)

  print *,maxd,maxdd,"hu"

123  nbo=2*maxd/maxdd ! the number of boxes
   print *,nbo,"nbo"
!  allocate(bo(-1:nbo,-1:nbo,-1:nbo))
  if (allocated(bo)) deallocate(bo)
  allocate(bo(-1:nbo,-1:nbo,-1:nbo))
  if (allocated(bom)) deallocate(bom)
  allocate(bom(-1:nbo,-1:nbo,-1:nbo,nbo))

  com=0.0d0

  do ir=3,nbo !we make at least 3 partitions per dimension
    bo=0
    fr=0
    coma=0
    v=0.0d0
    sc=real(ir)/maxd
!print *,ir,sc,"ir sc",maxd
    do i=1,nd
      ii=int(node(i)%x*sc+0.5d0);jj=int(node(i)%y*sc+0.5d0);kk=int(node(i)%z*sc+0.5d0)
!print *,node(i)%x,node(i)%y,node(i)%z,"x y z"
!print *,ii,jj,kk,"ii jj kk",sc
!print *,i,ir,maxd,sc,node(i)%x*sc,node(i)%y*sc,node(i)%z*sc
      bo(ii,jj,kk)=1
    end do

    bom(:,:,:,ir)=bo

    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          v=v+bo(i,j,k)
        end do
      end do
    end do
    kk=0
    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          if (bo(i,j,k)>0) then
            print *,i,j,k,"kok"
            kk=kk+1
          end if
        end do
      end do
    end do


!    print *,ir,v,v/nbo**3,nbo**3,"volum"

    ! now the mutual information
    b=0.0d0
    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          ce=bo(i,j,k)
          if (ce==1) then ! otherwise we are counting totally empty areas and this produces a bias based on the boxes and 
                                           ! the shape of the embryo
            ri=bo(i+1,j,k)
            le=bo(i-1,j,k)
            an=bo(i,j+1,k)
            po=bo(i,j-1,k)
            up=bo(i,j,k+1)
            do=bo(i,j,k-1)
            fr(ri,le,an,po,up,do)=fr(ri,le,an,po,up,do)+1
!print *,ri,le,an,po,up,do,fr(ri,le,an,po,up,do)
            b=b+1
          end if
        end do
      end do
    end do

    a=1.0d0/b
    fr=fr*a
!print *,a,"a"
      do j=0,1
        do k=0,1
          do ii=0,1
            do jj=0,1
              do kk=0,1
                do iii=0,1
                  if (fr(j,k,ii,jj,kk,iii)>0.0d0) then
                    coma=coma-fr(j,k,ii,jj,kk,iii)*log(fr(j,k,ii,jj,kk,iii))
!print *,b,j,k,ii,jj,kk,iii,fr(j,k,ii,jj,kk,iii)/a,nd,fr(j,k,ii,jj,kk,iii),"estats"
!print *,ir,nint(fr(j,k,ii,jj,kk,iii)/a),fr(j,k,ii,jj,kk,iii),nd,"ww"
!                    if (nint(fr(j,k,ii,jj,kk,iii)/a)==nd) then ; com=com+coma ;goto 111 ; end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      com=com+coma   !fitness is the sum of the complexity at each level
print *,ir,coma,"coma",v,nbo,"nbo",maxd
!    if (v==nd) exit
  end do

111  print *,com,nbo,maxd,maxdd,"fitness"

end subroutine

!************************************************************************************

subroutine ord4_pco_per_logpco_forcells(ind,com)   ! measures complexity as the mutual information between embryo and non-embryo boxes for different scales
  integer ind,i,j,k,ii,jj,kk,ir,iii       ! NOT NORMALIZED BY SIZE, but with a minimal size being the minval
  real*8 com,sc,v,a,b,coma,ix,iy,iz,d,ax,ay,az,bx,by,bz,angleu,aa,bb,maxd
  integer up,do,an,po,le,ri,ce
  real*8 ma,mb,mc
  integer, parameter :: nboo=10
  integer, allocatable :: bo(:,:,:)
  real*8 fr(0:1,0:1,0:1,0:1,0:1,0:1)

  ! we calculculate the complexity as the join information based on the probability of each kind of neighborhood at radium one box

  ! it is independent of size because we split each embryo in the same number of boxes 

  ! we have to rotate everything so that it gets properly aligned to make the boxes along the longest axis of the embryo
  ! we first calculate the distance between the most distant nodes in the embryo
  a=0.0d0
  d=0.0d0
  do i=1,nd
    ix=node(i)%x ; iy=node(i)%y ; iz=node(i)%z
    do j=1,nd
      ax=node(j)%x-ix ; ay=node(j)%y-iy ; az=node(j)%z-iz
      d=sqrt(ax**2+ay**2+az**2)
      if (d>a) then
        a=d
        bx=ax ; by=ay ; bz=az
        ii=i  ; jj=j
      end if
    end do
  end do

  maxd=a
  if (maxd<1.0d0) maxd=2.0d0  !THIS IS TO FORBIDE REALLY SMALL
  ! we calculate the angle between the vector of the longest distance and the x axis
  if (bz==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+bz**2)
    angleu=acos(abs(bx)/b)
  end if

  !now we rotate the whole embryo for this angle 
  do i=1,nd
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  if (bx==0.0.or.a==0.0d0) then
    angleu=0.0d0
  else 
    b=sqrt(bx**2+by**2)
    angleu=acos(abs(bx)/b)
  end if

  !now we rotate the whole embryo for this angle 
  do i=1,nd
    aa=node(i)%x*cos(angleu)-node(i)%z*sin(angleu)
    bb=node(i)%x*sin(angleu)+node(i)%z*cos(angleu)
    node(i)%x=aa
    node(i)%z=bb
  end do

  !now we reposition the embryo so that it is all in positive coordinates
  ! ma mb mc
  ma=0.0d0 ; mb=0.0d0 ; mc=0.0d0
  do i=1,nd
    if (ma>node(i)%x) ma=node(i)%x
    if (mb>node(i)%y) mb=node(i)%y
    if (mc>node(i)%z) mc=node(i)%z
  end do

  do i=1,nd
    node(i)%x=node(i)%x-ma    
    node(i)%y=node(i)%y-mb    
    node(i)%z=node(i)%z-mc    
  end do

  !call extrem 

  nbo=7           ! that by *2+1 is the number of partitions per dimension
  allocate(bo(-1:nbo,-1:nbo,-1:nbo))

  com=0.0d0

  do ir=1,nbo
    bo=0
    fr=0
    coma=0
    v=0.0d0
    sc=real(ir)/maxd
    do i=1,nd
      if (node(i)%marge==0) then
        if (node(i)%tipus==2.or.node(i)%tipus==3) then
          ii=int(node(i)%x*sc);jj=int(node(i)%y*sc);kk=int(node(i)%z*sc)
          bo(ii,jj,kk)=1
        end if
      end if
    end do

    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          v=v+bo(i,j,k)
        end do
      end do
    end do
    print *,ir,v,v/nbo**3,"volum"

    ! now the mutual information
    b=0.0d0
    do i=0,nbo-1
      do j=0,nbo-1
        do k=0,nbo-1
          ce=bo(i,j,k)
          if (ce==1) then ! otherwise we are counting totally empty areas and this produces a bias based on the boxes and 
                                           ! the shape of the embryo
            ri=bo(i+1,j,k)
            le=bo(i-1,j,k)
            an=bo(i,j+1,k)
            po=bo(i,j-1,k)
            up=bo(i,j,k+1)
            do=bo(i,j,k-1)
            fr(ri,le,an,po,up,do)=fr(ri,le,an,po,up,do)+1
!print *,ri,le,an,po,up,do,fr(ri,le,an,po,up,do)
            b=b+1
          end if
        end do
      end do
    end do

    a=1.0d0/b
    fr=fr*a
      do j=0,1
        do k=0,1
          do ii=0,1
            do jj=0,1
              do kk=0,1
                do iii=0,1
                  if (fr(j,k,ii,jj,kk,iii)>0.0d0) then
                    coma=coma-fr(j,k,ii,jj,kk,iii)*log(fr(j,k,ii,jj,kk,iii))
!print *,b,j,k,ii,jj,kk,iii,fr(j,k,ii,jj,kk,iii)/a,fr(j,k,ii,jj,kk,iii),"estats"
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

    com=com+coma   !fitness is the sum of the complexity at each level
print *,ir,coma,"coma"
    if (v==ncels) exit
  end do

  print *,com,"fitness"
end subroutine

!****************************************************************************************************************

subroutine other_control(ind,com)  !this is to detect when the cells have fused
  integer ind
  real*8 com
  integer i,j,k,ii,jj,kk,iii,jjj,kkk
  real*8 a,b,c,aa,bb,cc,d,dd

  if (ncels<2) then ; com=1.0d0 ; return ; end if

  kkk=0
  iii=0

    !UPDATING CELL CENTROIDS
    do i=1,ncels
      a=0 ; b=0 ; c=0
      do j=1,cels(i)%nunodes
        k=cels(i)%node(j)
        if(node(k)%tipus==1.or.node(k)%tipus==3)then
          a=a+node(k)%x ; b=b+node(k)%y ; c=c+node(k)%z
        end if
      end do
      d=1d0/real(cels(i)%nunodes)
      if(node(cels(i)%node(1))%tipus<3) d=2d0*d !if it's epithelial; peta aqui, veces:3
      cels(i)%cex=a*d ; cels(i)%cey=b*d ; cels(i)%cez=c*d
    end do


  do i=1,ncels

    ! average distance between a cell center and its nodes
    d=0.0d0
    do k=1,cels(i)%nunodes
      ii=cels(i)%node(k)    
      d=d+sqrt((node(ii)%x-cels(i)%cex)**2+(node(ii)%y-cels(i)%cey)**2+(node(ii)%z-cels(i)%cez)**2)
    end do
    a=cels(i)%nunodes
    d=d/a
!print *,d,"d"
!print *,i,j,ncels,"ncels"
    do j=i+1,ncels
      
      ! average distance between the center of the previous cell and the nodes of another cell
      dd=0.0d0
      do k=1,cels(j)%nunodes
        ii=cels(j)%node(k)    
        dd=dd+sqrt((node(ii)%x-cels(i)%cex)**2+(node(ii)%y-cels(i)%cey)**2+(node(ii)%z-cels(i)%cez)**2)
      end do
      a=cels(j)%nunodes
      dd=dd/a

!      print *,d,dd,1-dd/d,"d dd 1-dd/d"
      if (dd/d<1.5) kkk=kkk+1
      iii=iii+1

    end do
  end do

  a=kkk*2 ; b=iii
  print *,kkk,iii  

  com=1-a/b  !this is the factor of which proportion of the cells overlap too much
!  print *,com,"com"
end subroutine


!*************************************************************************************************************************************
!EMD + gene per cell
subroutine functional_complexity(ind,EMD00,pernofi3)  
 real*8, allocatable  :: geneterri1a(:),geneterri2a(:),ngmax(:),prepoints0(:,:),prepoints(:,:)
 character*140 :: pernofi3
 integer :: ind
 real*8 :: EMD00,EMD01,scfit,aaa,maxdistance,start_time,stop_time,start_time0,stop_time0
 integer :: jjcom,jcom,icom,iicom,ngnumber,ndYY,ngYY,ndYYb,ndYYa,j_nd,i_nd,closest_nd
 real*8 :: closest_a,closest_b,closest_c
 type :: nodepatern
    integer :: pattern_id
    integer, allocatable :: genepattern(:),genepatterncopy(:)
 end type
 integer, allocatable :: ndExt(:),ndIntra(:)
 real*8, allocatable  :: aap(:),dIntra(:),dExt(:)
 integer :: id_pat,ab
 real*8 ::d1,d2,pi 

 real*8 :: a1, a2, b1, b2, c1, c2, abb, acc, genecon
 integer:: j_altre, metapatterns, typend
 
type(nodepatern),allocatable:: nodep(:)

if(allocated(nodep))deallocate(nodep)
allocate(nodep(nd))
do i=1,nd
 if(allocated(nodep(i)%genepattern))deallocate(nodep(i)%genepattern)
 allocate(nodep(i)%genepattern(ng))

 if(allocated(nodep(i)%genepatterncopy))deallocate(nodep(i)%genepatterncopy)
 allocate(nodep(i)%genepatterncopy(ng))
enddo

!sc=0.0
call cpu_time(start_time0)
call iniread
call readsnap(pernofi3)
EMD00=0d0;EMD01=0d0
aaa=0d0
ngnumber=0
ndYY=0;ngYY=0
!if (allocated(geneterri1a))deallocate(geneterri1a)
!allocate(geneterri1a(nd))
!if (allocated(geneterri2a))deallocate(geneterri2a)
!allocate(geneterri2a(nd))
if (allocated(ngmax))deallocate(ngmax)
allocate(ngmax(ng))
ngmax=0d0
maxdistance=0d0
do i=1,nd
  do ii=1,nd
     aaa=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+&
     &(node(i)%z-node(j)%z)**2)
     if(aaa>maxdistance)maxdistance=aaa
enddo;enddo

aaa=0
do i=1,ng 
  do j=1,nd
    if(ngmax(i)<gex(j,i))ngmax(i)=gex(j,i)!;print*,"ngmax",i,"--",ngmax(i);endif
  enddo
enddo
do i=1,ng 
    if(ngmax(i)>0d0)ngnumber=ngnumber+1!;print*,"ngmax",i,"--",ngmax(i);endif
enddo

do i=1,nd
   do ii=1,ng
    ab=1
    if(gex(i,ii).lt.1d-8)ab=0
    if(ngmax(ii)>0)then
       genecon=gex(i,ii)/ngmax(ii)*ab 
    else;genecon=0d0;endif 
    genecon=genecon*10 ; genecon=int(genecon)
    nodep(i)%genepattern(ii)=genecon    
    nodep(i)%genepatterncopy(ii)=genecon  
   enddo
enddo

if (allocated(aap))deallocate(aap)
allocate(aap(ng))

metapatterns=0
do i=1,nd
  aap=nodep(i)%genepatterncopy 
  nodep(i)%genepatterncopy=1000
  if(sum(aap)>999)cycle
  metapatterns=metapatterns+1;print*,"metapatterns",metapatterns
  do j=1,nd 
    if(all(aap==nodep(j)%genepattern))then
      nodep(j)%genepatterncopy=1000
      nodep(j)%pattern_id=metapatterns
    endif
  enddo
enddo

print*,"****"
do i=1,metapatterns
    print*,nodep(i)%genepattern
enddo
print*,"****"



!print*,nodep(:)%pattern_id

!stop

if (allocated(dIntra))deallocate(dIntra)
allocate(dIntra(metapatterns))

if (allocated(dExt))deallocate(dExt)
allocate(dExt(metapatterns))

if (allocated(ndExt))deallocate(ndExt)
allocate(ndExt(metapatterns))

if (allocated(ndIntra))deallocate(ndIntra)
allocate(ndIntra(metapatterns))
print*,"2689--fitmo"
dIntra=0;dExt=0
ndExt=0;ndIntra=0
!do metap=1,metapattern
 do i_nd=1,nd
  ! dIntra(i_nd)=0;dExt(i_nd)=0
  ! ndExt(i_nd)=0;ndIntra(i_nd)=0
   a1=node(i_nd)%x ; b1=node(i_nd)%y ; c1=node(i_nd)%z
   typend=node(i_nd)%tipus
   !aap=nodep(i_nd)%genepattern
   id_pat=nodep(i_nd)%pattern_id
   do j_nd=1,nd
     if(i_nd==j_nd)cycle           
     if(node(j_nd)%tipus.ne.typend)then
        if(node(j_nd)%tipus==1.or.node(j_nd)%tipus==2)then
          j_altre=node(j_nd)%altre
          a2=node(j_altre)%x ; b2=node(j_altre)%y ; c2=node(j_altre)%z 
        else 
          a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z 
        endif
     else 
          a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z   
     endif   
     if(id_pat==nodep(j_nd)%pattern_id)then  
       dIntra(id_pat)=((a2-a1)**2+(b2-b1)**2+(c2-c1)**2)/maxdistance+dIntra(id_pat)
       ndIntra(id_pat)=ndIntra(id_pat)+1
     else
       dExt(id_pat)=((a2-a1)**2+(b2-b1)**2+(c2-c1)**2)/maxdistance+dExt(id_pat)
       ndExt(id_pat)=ndExt(id_pat)+1
     endif
   enddo
 enddo
!enddo

print*,"2721--fitmo"
EMD00=0d0
do i=1,metapatterns
  print*,"metapatterns",i
  abb=real(ndIntra(i))!/nd
  print*,"nd intra",abb
!  print*,"d intra",dIntra(i)!/abb 
  
  acc=real(ndExt(i))!/nd
  print*,"nd ext",acc
!  print*,"d ext",dExt(i)!/acc
  
  if(dIntra(i)==0)then;dIntra(i)=0;abb=nd;endif

  d1=dIntra(i)/abb !average distance
  d2=dExt(i)/acc
  if(d1==0)d1=d2 !only one nd, no more intra
  if(acc==0)cycle;if(dExt(i)==0)cycle      !pattern everywhere
  d1=d1/d2 ; pi=1/abb 
  EMD00=-1*(d1*pi*log(pi))+EMD00
!  print*,"emd00",emd00
enddo
print*,"2728--fitmo"
!EMD00=aaa/ndYY
print*,"EMD00",EMD00,"*************************"
a=0
!print*,"mp",ndIntra
do i=1,metapatterns
   if(ndIntra(i)==0)then;a=a+1
   else;a=a+ndIntra(i);endif
enddo
print*,"all nd intra",a,"nd",nd


stop



!aaa=0d0
!do ii=1,ng
!   if(ngmax(ii)==0d0)cycle
!   do jj=1,ng
!     if(ii==jj)cycle
!     if(ngmax(jj)==0d0)cycle
!     do i_nd=1,nd
!        if(gex(i_nd,ii)<1e-8)cycle
!        a1=node(i_nd)%x ; b1=node(i_nd)%y ; c1=node(i_nd)%z
!        typend=node(i_nd)%tipus
!        aa=1000000
!        do j_nd=1,nd
!           if(i_nd==j_nd)cycle
!           if(gex(j_nd,jj)<1e-8)cycle
!           if(node(j_nd)%tipus.ne.typend)then
!             if(node(j_nd)%tipus==1.or.node(j_nd)%tipus==2)then
!               j_altre=node(j_nd)%altre
!               a2=node(j_altre)%x ; b2=node(j_altre)%y ; c2=node(j_altre)%z 
!             else 
!               a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z 
!             endif
!           else 
!             a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z   
!           endif           !
!
!           dd=(a2-a1)**2+(b2-b1)**2+(c2*2-c1)**2
!           if (dd<aa)then;aa=dd;closest_a=a2;closest_b=b2;closest_c=c2
!               closest_nd=j_nd
!           endif
!        enddo
!        ndYY=ndYY+1
!aaa=aaa+sqrt(((a1-node(closestnd)%x)**2)/maxdistance+((b2-node(closestnd)%y)**2)/maxdistance+&
!((c1-node(closestnd)%z)**2)/maxdistance)*sqrt(((gex(i_nd,ii)+ngmax(ii))/ngmax(ii))-(gex(closestnd,jj)+&
!ngmax(ii))/ngmax(jj))**2) 

!aaa=aaa+sqrt(((a1-closest_a)**2)/maxdistance+((b1-closest_b)**2)/maxdistance+&
!((c1-closest_c)**2)/maxdistance)*sqrt(((gex(i_nd,ii)+ngmax(ii))/ngmax(ii)-&
!(gex(closest_nd,jj)+ngmax(jj))/ngmax(jj))**2)        
!     enddo
!   enddo
!enddo

!EMD00=aaa/ndYY
!print*,"EMD00",EMD00,"*************************"
!stop



 ndYY=0
 ndYYa=0
 ndYYb=0
 do icom=1,ng-1 !ojo 1-->2!
    if(ngmax(icom)==0d0)cycle
    print*,"gene",icom,"compare to ***************"
    ndYYa=0
    do iicom=1,nd
     if(gex(iicom,icom)>1e-8)ndYYa=ndYYa+1
    enddo
    if(allocated(geneterri1a))deallocate(geneterri1a)
    allocate(geneterri1a(ndYYa))
    geneterri1a=0d0
    if(allocated(prepoints0))deallocate(prepoints0)
    allocate(prepoints0(ndYYa,3))
    prepoints0=0
    ndYY=0
    do iicom=1,nd
      if(gex(iicom,icom)>1e-8)then
        ndYY=ndYY+1
        geneterri1a(ndYY)=gex(iicom,icom)/ngmax(icom)
        prepoints0(ndYY,1)=node(iicom)%x
        prepoints0(ndYY,2)=node(iicom)%y
        prepoints0(ndYY,3)=node(iicom)%z
      endif
    enddo
!print*,prepoints0 ;stop
    do jcom=icom+1,ng
     if(ngmax(jcom)==0d0)cycle
     ndYYb=0
     do iicom=1,nd
      if(gex(iicom,jcom)>1e-8)ndYYb=ndYYb+1
     enddo
     if (allocated(geneterri2a))deallocate(geneterri2a)
     allocate(geneterri2a(ndYYb))
     geneterri2a=0d0
     if(allocated(prepoints))deallocate(prepoints)
     allocate(prepoints(ndYYb,3))
     prepoints=0
     print*,"gene",jcom," ************************"
     ndYY=0
     do jjcom=1,nd
      if(gex(jjcom,jcom)>1e-8)then
        ndYY=ndYY+1
        geneterri2a(ndYY)=gex(jjcom,jcom)/ngmax(jcom)
        prepoints(ndYY,1)=node(jjcom)%x
        prepoints(ndYY,2)=node(jjcom)%y
        prepoints(ndYY,3)=node(jjcom)%z
      endif
     enddo
!     call cpu_time(start_time)
     if(ndYYb>ndYYa)then
  !OJO    call conservativecomplexity(scfit,geneterri1a,geneterri2a,pernofi3,prepoints0,ndYYa,prepoints,ndYYb)
     else
  !OJO    call conservativecomplexity(scfit,geneterri2a,geneterri1a,pernofi3,prepoints,ndYYb,prepoints0,ndYYa)
     endif
!     call cpu_time(stop_time)  ; print*,"time********",(stop_time-start_time)/60;stop
     EMD01=EMD01+scfit
     ngYY=ngYY+1
!     print*,"EMD01",EMD01,"ngYY",ngYY
    enddo
 enddo
!print*,"ngYY",ngYY,"(ngnumber-1)*ngnumber",(ngnumber-1)*ngnumber
EMD01=EMD01/ngYY!((ngnumber-1)*ngnumber)
!print*,"EMD01 pre *1d3",EMD01
!EMD01=EMD01*1d3
print*,"active genes",ngnumber
print*,"EMD00",EMD00
print*,"EMD01",EMD01
EMD00=EMD00*EMD01
print*,"EMD final",EMD00,"**************"
call cpu_time(stop_time0)
print*,stop_time0-start_time0,"***time in fitness***"
!stop !****ojo*******
end subroutine functional_complexity


!*************************************************************************************************************************************
!EMD + gene per cell
subroutine functional_complexity2(ind,EMD00,pernofi3)  
 real*8, allocatable  :: geneterri1a(:),geneterri2a(:),ngmax(:),prepoints0(:,:),prepoints(:,:)
 character*140 :: pernofi3
 real*8 :: EMD00,EMD01,scfit,aaa,maxdistance,start_time,stop_time,start_time0,stop_time0
 integer :: ind, jjcom,jcom,icom,iicom,ngnumber,ndYY,ngYY,ndYYb,ndYYa,j_nd,i_nd,closest_nd,fhj
 real*8 :: closest_a,closest_b,closest_c
 type :: nodepatern
    integer :: pattern_id
    integer, allocatable :: genepattern(:),genepatterncopy(:)
 end type
 integer, allocatable :: ndExt(:),ndIntra(:)
 real*8, allocatable  :: aap(:),dIntra(:),dExt(:)
 integer :: id_pat,ab,yes
 real*8 ::d1,d2,pi 
 real*8 :: a1, a2, b1, b2, c1, c2, abb, acc, genecon
 integer:: j_altre, metapatterns, typend
type(nodepatern),allocatable:: nodep(:)

if(allocated(nodep))deallocate(nodep)
allocate(nodep(nd))
do i=1,nd
 if(allocated(nodep(i)%genepattern))deallocate(nodep(i)%genepattern)
 allocate(nodep(i)%genepattern(ng))

 if(allocated(nodep(i)%genepatterncopy))deallocate(nodep(i)%genepatterncopy)
 allocate(nodep(i)%genepatterncopy(ng))
enddo

!sc=0.0
call cpu_time(start_time0)
call iniread
call readsnap(pernofi3)
EMD00=0d0;EMD01=0d0
aaa=0d0
ngnumber=0
ndYY=0;ngYY=0

if (allocated(ngmax))deallocate(ngmax)
allocate(ngmax(ng))
ngmax=0d0
maxdistance=0d0
do i=1,nd
  do ii=1,nd
     aaa=sqrt((node(i)%x-node(j)%x)**2+(node(i)%y-node(j)%y)**2+&
     &(node(i)%z-node(j)%z)**2)
     if(aaa>maxdistance)maxdistance=aaa
enddo;enddo

aaa=0
do i=1,ng 
  do j=1,nd
    if(ngmax(i)<gex(j,i))ngmax(i)=gex(j,i)!;print*,"ngmax",i,"--",ngmax(i);endif
  enddo
enddo
do i=1,ng 
    if(ngmax(i)>0d0)ngnumber=ngnumber+1!;print*,"ngmax",i,"--",ngmax(i);endif
enddo

do i=1,nd

   do ii=1,ng
    ab=1
    if(gex(i,ii).lt.1d-8)ab=0
    if(ngmax(ii)>0)then
       genecon=gex(i,ii)/ngmax(ii)*ab!;print*,"gex",i,ii,gex(i,ii);stop
    else;genecon=0d0;endif 
    genecon=genecon*10 ; genecon=int(genecon)!;print*,"genecon",i,ii,genecon;stop
    nodep(i)%genepattern(ii)=genecon    
    nodep(i)%genepatterncopy(ii)=genecon  
   enddo
enddo

!do i=1,nd
! print*,i,nodep(i)%genepatterncopy 
!enddo

if (allocated(aap))deallocate(aap)
allocate(aap(ng))

metapatterns=0
do i=1,nd
  aap=nodep(i)%genepatterncopy 
  nodep(i)%genepatterncopy=1000
  if(sum(aap)>999)cycle
  metapatterns=metapatterns+1!;print*,"metapatterns",metapatterns,"--",aap
  fhj=0
  do j=1,nd 
    if(all(aap==nodep(j)%genepattern))then
      nodep(j)%genepatterncopy=1000
      nodep(j)%pattern_id=metapatterns ; fhj=fhj+1
    endif
  enddo
!  print*,"fhj",fhj

enddo



!do i=1,nd
! print*,i,nodep(i)%pattern_id
!enddo
!stop
if (allocated(dIntra))deallocate(dIntra)
allocate(dIntra(metapatterns))

if (allocated(dExt))deallocate(dExt)
allocate(dExt(metapatterns))

if (allocated(ndExt))deallocate(ndExt)
allocate(ndExt(metapatterns))

if (allocated(ndIntra))deallocate(ndIntra)
allocate(ndIntra(metapatterns))


!print*,"2689--fitmo"
dIntra=0;dExt=0
ndExt=0;ndIntra=0


do i=1,nd
   do ii=1,ng
    nodep(i)%genepatterncopy(ii)=nodep(i)%genepattern(ii)
enddo;enddo

metapatterns=0
do i=1,nd
  aap=nodep(i)%genepatterncopy 
  nodep(i)%genepatterncopy=1000
  if(sum(aap)>999)cycle
  metapatterns=metapatterns+1!;print*,"metapatterns",metapatterns,"--",aap
  fhj=0 ; yes=1
  do j=1,nd 
    if(all(aap==nodep(j)%genepattern))then
      nodep(j)%genepatterncopy=1000
      nodep(j)%pattern_id=metapatterns ; fhj=fhj+1
    endif
  enddo
  if(yes==1)then
!   print*,"fhj",fhj
   ndIntra(metapatterns)=fhj
   ndExt(metapatterns)=nd-fhj
   yes=0
  endif
enddo


!do i=1,metapatterns
! print*,ndIntra(i),ndExt(i)
!enddo


!do metap=1,metapattern
 do i_nd=1,nd
  ! dIntra(i_nd)=0;dExt(i_nd)=0
  ! ndExt(i_nd)=0;ndIntra(i_nd)=0
   a1=node(i_nd)%x ; b1=node(i_nd)%y ; c1=node(i_nd)%z
   typend=node(i_nd)%tipus
   !aap=nodep(i_nd)%genepattern
   id_pat=nodep(i_nd)%pattern_id
   do j_nd=1,nd
     if(i_nd==j_nd)cycle           
     if(node(j_nd)%tipus.ne.typend)then
        if(node(j_nd)%tipus==1.or.node(j_nd)%tipus==2)then
          j_altre=node(j_nd)%altre
          a2=node(j_altre)%x ; b2=node(j_altre)%y ; c2=node(j_altre)%z 
        else 
          a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z 
        endif
     else 
          a2=node(j_nd)%x ; b2=node(j_nd)%y ; c2=node(j_nd)%z   
     endif   
     if(id_pat==nodep(j_nd)%pattern_id)then  
       dIntra(id_pat)=((a2-a1)**2+(b2-b1)**2+(c2-c1)**2)/maxdistance+dIntra(id_pat)
     !  ndIntra(id_pat)=ndIntra(id_pat)+1
     else
       dExt(id_pat)=((a2-a1)**2+(b2-b1)**2+(c2-c1)**2)/maxdistance+dExt(id_pat)
     !  ndExt(id_pat)=ndExt(id_pat)+1
     endif
   enddo
 enddo
!enddo

EMD00=0d0
do i=1,metapatterns
 ! print*,"metapatterns",i
  abb=real(ndIntra(i))!/nd
!  print*,"nd intra",abb
!  print*,"d intra",dIntra(i)!/abb 
  
  acc=real(ndExt(i))!/nd
!  print*,"nd ext",acc
!  print*,"d ext",dExt(i)!/acc
  
  if(dIntra(i)==0)then;dIntra(i)=0;abb=nd;endif

  d1=dIntra(i)/abb !average distance
  d2=dExt(i)/acc
  if(d1==0)d1=d2 !only one nd, no more intra
  if(acc==0)cycle;if(dExt(i)==0)cycle      !pattern everywhere
  d1=d1/d2 ; pi=1/abb ;! print*,"d1",d1,"pi",pi,"pattern",i,-1*(d1*pi*log(pi))&
                       ! &,"intra",ndIntra(i),"ext",ndExt(i)
  EMD00=-1*(d1*pi*log(pi))+EMD00
!  print*,"emd00",emd00
enddo

call total(ind,emd01,10)

print*,"EMD00-gene distribution",EMD00
print*,"EMD01-angle variation",EMD01
EMD00=EMD00*EMD01
print*,"EMD final",EMD00,"**************"
call cpu_time(stop_time0)
print*,stop_time0-start_time0,"***time in fitness***"
!stop !****ojo*******
end subroutine functional_complexity2

end module
