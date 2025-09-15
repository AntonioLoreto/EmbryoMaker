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




module shell			!!!!!!!!!!!!!!!!!!!!!!!!!!!! miguel made module 22-5-13

use general

real*8,public,allocatable   :: shellp(:,:)
integer,public,allocatable  :: trs(:,:)
integer,public              :: pts
real*8,public               :: r1,r2,r11,r22,exc,vol,cte

contains

subroutine shellinicial
  !r1=1d0                     ! parameter 1 ; radius of the sphere. It determines the volume of the spheroid. 
  r1=(nd**(1d0/3d0))*maxval(node(:nd)%eqd)    !calculated as the sum of volumes of the "req" spheres for all nodes
  write(*,*)'r1',r1
  exc=0.75d0                 ! parameter 2 ; excentricity (r1=exc*r2)  
  ! exc<0   : as in exc>0
  ! exc=0   : straight line (0,0,Z)
  ! 0<exc<1 : elogated speroid (c.elegans shape)
  ! exc=1   : sphere
  ! exc>1   : discoidal spheroid (flattened shape)
  exc=1d0/exc                                            ! (because in our system r1<r2)
  vol=(4d0/3d0)*pi*(r1**3)                               ! sphere volume
  r2=((3d0*vol*(exc**2))/(4d0*pi))**(1d0/3d0)            ! second radius
  r1=r2/exc                                              ! recalculate the first radius so that volume is kept constant 

  if(r1.lt.1d0)then 					 ! r1 cannot be (for calculations) r1<1 
  cte=1d0/r1 ; r11=r1*cte ; r22=r2*cte  
  else 
  r11=r1 ; r22=r2 
  endif

  do i=1,nd
  movi=0
  call eggshell(i)
  if(movi.eq.1)then;write(*,*)'WARNING ! There are nodes outside the eggshell from the beginning!';end if
  end do

  call shellplot
end subroutine shellinicial

subroutine shellplot
integer              :: ptsc,trsc
real*8               :: pex,pey,pez,alfa,marg

  pts=20 ; ptsc=1 ; marg=maxval(node(:nd)%eqd)

  if(allocated(shellp))then;deallocate(shellp);endif
  allocate(shellp(2+int(pts*pts),3)) 
  if(allocated(trs))then;deallocate(trs);endif
  allocate(trs(int(pts*pts*2),3)) ; trs=0

  shellp(int(1+pts*pts),1)=0d0 ; shellp(int(1+pts*pts),2)=0d0 ; shellp(int(1+pts*pts),3)=-1d0*r2
  shellp(int(2+pts*pts),1)=0d0 ; shellp(int(2+pts*pts),2)=0d0 ; shellp(int(2+pts*pts),3)=r2

  do i=1,pts
      pez=-1d0*r2+((2d0*(r2))/real(pts))*real(i)
      pez=pez-(2d0*(r2))/(2d0*real(pts))
      pex=sqrt(((r1)**2)*(1d0-((pez**2)/((r2)**2))))    
      do j=1,pts
          alfa=(2d0*pi/real(pts))*real(j)
          shellp(ptsc,1)=sin(alfa)*pex 
          shellp(ptsc,2)=cos(alfa)*pex
          shellp(ptsc,3)=pez     
          ptsc=ptsc+1
      end do
  end do
  trsc=pts+1
  do i=1,pts
      trs(i,1)=int(1+pts*pts) ; trs(i,2)=i ; trs(i,3)=i+1 ; if(i.eq.pts)then; trs(i,3)=1 ; endif
      trs(int(pts*pts*2)-i+1,1)=int(2+pts*pts) ; trs(int(pts*pts*2)-i+1,2)=int(pts*pts)-i+1 
      trs(int(pts*pts*2)-i+1,3)=int(pts*pts)-i ; if(i.eq.pts)then; trs(int(pts*pts*2)-i+1,3)=int(pts*pts) ; endif
      if(i.ne.pts)then
          do j=1,pts
              trs(trsc,1)=trsc-pts ; trs(trsc,2)=(trsc-pts)+1 ; trs(trsc,3)=trsc           
              trs(trsc+(int((pts-1)*pts)),1)=(trsc-pts)+1          
              trs(trsc+(int((pts-1)*pts)),2)=trs(trsc+(int((pts-1)*pts)),1)+pts
              trs(trsc+(int((pts-1)*pts)),3)=trs(trsc+(int((pts-1)*pts)),1)+pts-1 
              if(mod(trsc-pts,pts).eq.0)then; trs(trsc,2)=trs(trsc,2)-pts
                  trs(trsc+(int((pts-1)*pts)),1)=(trsc-2*pts)+1 
                  trs(trsc+(int((pts-1)*pts)),2)=(trsc-pts)+1
                  trs(trsc+(int((pts-1)*pts)),3)=trs(trsc+(int((pts-1)*pts)),1)+pts-1+pts
              endif          
              trsc=trsc+1    
          end do
      end if
  end do
end subroutine shellplot

subroutine eggshell(nod)
integer :: nod
real*8  :: ds,dd,ie,ix,iy,iz,px,py,pz,b,d
real*8  :: sa,sb,sc,vx,vy,vz,ssx,ssy,ssz,ppx,ppy,ppz,iix,iiy,iiz,kkk

  ie=0d0 ; dd=0d0 ; d=10d0
  ix=abs(node(nod)%x)  ; iy=abs(node(nod)%y)  ; iz=abs(node(nod)%z)                                           
  if(r1.lt.1d0)then                                      ! r1 cannot be (for calculations) r1<1   
  ix=ix*cte ; iy=iy*cte ; iz=iz*cte
  endif

  dd=sqrt((ix**2)+(iy**2))                               ! square root
  ix=dd ; iy=iz ; iz=0d0                                 ! change in coordinate system (3d->2d)
  b=iy/ix                                                ! slope
  px=sqrt(((r11**2)*(r22**2))/((r22**2)+(b**2*r11**2)))  ! projected X point in spheroid surface
  py=px*b
  dd=sqrt((ix**2)+(iy**2))                               ! euclidean distance in the new coordinate system
  ds=sqrt((px**2)+(py**2))                               ! distance of  the point in the spheroid surface 
  if(dd.ge.ds)then                                       ! energy contribution               
      movi=1
  end if
end subroutine eggshell

subroutine eggshell_forces(nod,px,py,pz)
integer :: nod
real*8  :: dd,ds,ie,ix,iy,iz,px,py,pz,b,d,eqa,eqb
real*8  :: sa,sb,sc,vx,vy,vz,ssx,ssy,ssz,ppx,ppy,ppz,iix,iiy,iiz,kkk

  ie=0d0 ; dd=0d0 ; d=10d0 
  ix=abs(node(nod)%x) ; iy=abs(node(nod)%y) ; iz=abs(node(nod)%z)  
  
  if(r1.lt.1d0)then                                      ! r1 cannot be (for initial calculations) r1<1   
  ix=ix*cte ; iy=iy*cte ; iz=iz*cte
  endif

  dd=sqrt((ix**2)+(iy**2))                               ! square root
  ix=dd ; iy=iz ; iz=0d0                                 ! change in coordinate system (3d->2d)
  b=iy/ix                                                ! slope
  ppx=sqrt(((r11**2)*(r22**2))/((r22**2)+(b**2*r11**2))) ! projected X point in spheroid surface
  ppy=ppx*b					         ! It is simpler taking the equation of a line
  dd=sqrt((ix**2)+(iy**2))                               ! euclidean distance in the new coordinate system
  ds=sqrt((ppx**2)+(ppy**2))                             ! distance of  the point in the spheroid surface

  if(dd.ge.ds)then                                       ! the node is outside ....
    iix=node(nod)%x ; iiy=node(nod)%y ; iiz=node(nod)%z 
    vx=iix-px ; vy=iiy-py ; vz=iiz-pz                    ! displacement vector (coming from iterdiferencial)   
  
    sa=(((vx**2)+(vy**2))/((r1**2)*(vz**2)))+(1d0/(r2**2))                       ! second degree equation
    sb=(2d0/(vz*r1**2))*((px*vx)+(py*vy)-((pz*vx**2)/vz)-((pz*vy**2)/vz))        ! intersection between a line
    sc=(1d0/(r1**2))*((px**2)-((2d0*px*pz*vx)/vz)+(py**2)-((2d0*py*pz*vy)/vz)+&  ! (defined by a point and a vector)
    &(((pz**2)*(vx**2))/(vz**2))+(((pz**2)*(vy**2))/(vz**2))-(r1**2))            ! and the ellipsoid surface 

    !eqa=iix*vz-iiz*vx                                                             ! alternative form of the equation
    !eqb=iiy*vz-iiz*vy                                                             ! (simpler??)
    !sa=(r2**2)*((vx**2)+(vy**2))+(r1**2)*(vz**2)
    !sb=2d0*(r2**2)*((vx*eqa)+(vy*eqb))
    !sc=(r2**2)*(eqa**2+eqb**2)-(r1**2)*(r2**2)*(vz**2)
 
    if(((sb**2)-4d0*sa*sc).lt.0d0)then;write(*,*)(sb**2)-4d0*sa*sc,'rneg',nod;return;endif                
    ssz=(-sb+sqrt((sb**2)-4d0*sa*sc))/(2d0*sa)                       ! first possible solution to Z coordinate 

    if(((ssz.gt.pz).and.(ssz.gt.iiz)).or.((ssz.lt.pz).and.(ssz.lt.iiz)))then ! if intersection has intermediate values ...
        ssz=(-sb-sqrt((sb**2)-4d0*sa*sc))/(2d0*sa)                     ! second possible solution to Z coordinate
    end if

    if(((ssz.gt.pz).and.(ssz.gt.iiz)).or.((ssz.lt.pz).and.(ssz.lt.iiz)))then !sometimes ¿it happens? because floating point
       kkk=0.5d0*(iiz+pz)
       if(abs(ssz-kkk).gt.abs(((-sb+sqrt((sb**2)-4d0*sa*sc))/(2d0*sa))-kkk))then ! proximity criterium
       ssz=(-sb+sqrt((sb**2)-4d0*sa*sc))/(2d0*sa)                  
       end if
    end if
    
    ! taking the equation of a line:
    ssx=(vx*((ssz-pz)/vz))+px
    ssy=(vy*((ssz-pz)/vz))+py

    ! the node trajectory is intercepted by the shell
    if(ssx.lt.iix)then; node(nod)%x=ssx-0.9*(abs(ssx-px));else;node(nod)%x=ssx+0.9*(abs(ssx-px));end if
    if(ssy.lt.iiy)then; node(nod)%y=ssy-0.9*(abs(ssy-py));else;node(nod)%y=ssy+0.9*(abs(ssy-py));end if
    if(ssz.lt.iiz)then; node(nod)%z=ssz-0.9*(abs(ssz-pz));else;node(nod)%z=ssz+0.9*(abs(ssz-pz));end if
end if

end subroutine eggshell_forces
end module shell

