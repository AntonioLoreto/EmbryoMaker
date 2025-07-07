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




!***************************************************************************
!***************  MODUL ***************************************************
!***************************************************************************

module aleas
!del alea3d
integer, public                ::nvalo,nvaloq  !en el fons tenim nvalo**2
real*8 , public, allocatable   ::particions_esfera(:,:) !s'hauria de fer nomes de out
integer, public, parameter     ::tamstocas=100000   !es per fer la matriu de ran2
integer, public                ::nunualea,cotran !numero de numeros aleas k realment tenim
real*8 , public                ::stocas(tamstocas)

contains

!**************************************************************************

subroutine inialea3d(un)
  integer i,j,k
  real*8 a,b,c,e,x,y,z
  integer un
 ! un=30 
  nvalo=un
  nvaloq=un**2

  if (allocated(particions_esfera)) deallocate(particions_esfera)
  allocate(particions_esfera(nvaloq,3))
  e=3.141592
  j=nvaloq/2
  
  do i=1,j
      call random_number(x)
      call random_number(y)
      x=2d0*e*x
      y=(2d0*y)-1d0          
      particions_esfera(i,1)=sqrt(1-y**2)*cos(x)
      particions_esfera(i,2)=sqrt(1-y**2)*sin(x)
      particions_esfera(i,3)=y      
      particions_esfera(i+j,1)=-1d0*sqrt(1-y**2)*cos(x)
      particions_esfera(i+j,2)=-1d0*sqrt(1-y**2)*sin(x)
      particions_esfera(i+j,3)=-1d0*y   
  end do
end subroutine inialea3d


subroutine inialea3d_old(un)
  integer i,j,k
  real*8 a,b,c,e,x,y,z,aa
  integer un
  real*8 vec(3)
  nvalo=un
  nvaloq=un*un
  if (allocated(particions_esfera)) deallocate(particions_esfera)
  allocate(particions_esfera(nvaloq,3))
  e=2d0*3.141592
  b=nvalo
  vec=0
  do i=1,nvalo
    do j=1,nvalo
      a=i
      a=e*a/b    !this is the angle     
      a=dcos(a)
      aa=j
      aa=e*aa/b  !>>Miquel27-3-14
      !aa=acos(2*aa/b-1)    !this is the other angle     
      c=sqrt(1d0-a**2)
      x=c*dcos(aa)  !>>Miquel27-3-14
      y=c*dsin(aa)  !>>Miquel27-3-14
      z=a
      k=(i-1)*nvalo+j
      particions_esfera(k,1)=x
      particions_esfera(k,2)=y
      particions_esfera(k,3)=z
      vec=vec+particions_esfera(k,:)
      !print*,i,j,"xyz",dcos(aa),dsin(aa),z
    end do
  end do
  !print*,"vec",vec
end subroutine inialea3d_old

!**************************************************************************
subroutine llaleat
    integer ui,i,ii,ios
    open(10,file='alea.dat',status='old',iostat=ios)
    if (ios/=0) then
      print *,"";      print *,"";      print *,"";
      print *,"there is no alea.dat file with aleatory numbers"
      print *,"";      print *,"";      print *,"";
      print *,"I can't tolerate that: f**ck you I'm quitting"
      print *,"";      print *,"";      print *,"";
      stop 
    end if
    nunualea=0
    cotran=1
    do ui=1,tamstocas
      read (10,*,ERR=100,END=90) stocas(ui)
      nunualea=nunualea+1
    end do
    return
100 print *,"";      print *,"";      print *,"";
    print *,"error in the random file alea.dat"
    print *,"";      print *,"";      print *,"";
    print *," I can't tolerate that: f**ck you I'm quitting"
    print *,"";      print *,"";      print *,"";
    stop 
90  close(10)
    open(10,file='alea.dat')
	print*,"end alea.dat"
	rewind(10)
	ii=int(a*1000)+1
print *,a,ii,"rew"
	do i=1,ii-1
          read(10,*)
	end do
end subroutine

end module
