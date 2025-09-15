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




module pola  !>>>>>>>MIQUEL MADE MODULE 3-4-13
use general
use neighboring
use io
use aleas


contains

subroutine pol_physic(celd,vx,vy,vz)		!this subroutine applies a physical criterion to determine the cell polarisation for division to take place
											!(namely it will polarize in the cell's longer axis)
											!a linear regression is calculated using all the cell's nodes and the resulting vector will be normal to the plane of division
integer,intent(in)::celd
integer::num,tipi
real*8::unum
integer::i,j,k,ii,jj,kk
real*8 ::d
real*8::a,b,c,xm,ym,zm,sxx,sxy,syy,sxz,szz,syz,thet,cost,sint
real*8::k11,k22,k12,k10,k01,k00
real*8::c2,c1,c0
real*8::dm,rho,phi,p,q,r
real*8::ka,kb,ku,kv,kw
real*8,intent(out)::vx,vy,vz
real*8,allocatable::px(:),py(:),pz(:)


	tipi=node(cels(celd)%node(1))%tipus

	if(tipi<3)then	!if epithelial, we only calculate the regression for the top layer (we ignore cell depth,because it causes artifacts if we take all the epithelial nodes)
		num=cels(celd)%nunodes
		allocate(px(num),py(num),pz(num))
		ii=0
		do i=1,num				!we load the cell's nodes
			j=cels(celd)%node(i)
			if(node(j)%tipus==1)then
				ii=ii+1
				px(ii)=node(j)%x
				py(ii)=node(j)%y
				pz(ii)=node(j)%z
!				print*,"px",px(i),"py",py(i),"pz",pz(i)
			end if
		end do
		num=ii
	else
		num=cels(celd)%nunodes
		unum=1d0/num
		allocate(px(num),py(num),pz(num))
!		print*,"num",num
		do i=1,num				!we load the cell's nodes
			j=cels(celd)%node(i)
			px(i)=node(j)%x
			py(i)=node(j)%y
			pz(i)=node(j)%z
!			print*,"px",px(i),"py",py(i),"pz",pz(i)
		end do
	end if

	unum=1d0/num
!	print*,"num",num

	xm=0;ym=0;zm=0			!averages
	do i=1,num
		xm=xm+px(i)
		ym=ym+py(i)
		zm=zm+pz(i)
	end do

!	print*,"before xm",xm,"ym",ym,"zm",zm
	
	xm=xm*unum
	ym=ym*unum
	zm=zm*unum

!	print*,"xm",xm,"ym",ym,"zm",zm

	sxx=0;sxy=0;syy=0;sxz=0;szz=0;syz=0		!covariances
	do i=1,num
		sxx=sxx+px(i)**2
		sxy=sxy+px(i)*py(i)
		syy=syy+py(i)**2
		sxz=sxz+px(i)*pz(i)
		szz=szz+pz(i)**2
		syz=syz+py(i)*pz(i)
	end do
	sxx=sxx*unum-(xm**2)
	sxy=sxy*unum-xm*ym
	syy=syy*unum-(ym**2)
	sxz=sxz*unum-xm*zm
	szz=szz*unum-(zm**2)
	syz=syz*unum-ym*zm
	if(sxx<=epsilod) sxx=1d-5       !if one of the covariances are 0 (which means the surface is perfectly flat on one coordinate plane
    if(sxy<=epsilod) sxy=1d-5       !the algorithm crashes, so we do the trick of assigning a very small value
    if(syy<=epsilod) syy=1d-5
    if(sxz<=epsilod) sxz=1d-5
    if(syz<=epsilod) syz=1d-5
    if(szz<=epsilod) szz=1d-5


!	print*,"sxx",sxx,"sxy",sxy,"syy",syy,"sxz",sxz,"szz",szz,"syz",syz


	!from here on it's fucking black magic!
	thet=0.5d0*atan(2*sxy/(sxx-Syy))

	cost=cos(thet)
	sint=sin(thet)

!	print*,"thet",thet,"cost",cost,"sint",sint

	k11=(syy+szz)*cost**2+(sxx+szz)*sint**2-2*sxy*cost*sint
	k22=(syy+szz)*sint**2+(sxx+szz)*cost**2+2*sxy*cost*sint
	k12=-sxy*(cost**2-sint**2)+(sxx-syy)*cost*sint
	k10=sxz*cost+syz*sint
	k01=-sxz*sint+syz*cost
	k00=sxx+syy

!	print*,"k11",k11,"k22",k22,"k12",k12,"k10",k10,"k01",k01,"k00",k00

	c2=-k00-k11-k22
	c1=k00*k11+k00*k22+k11*k22-k01**2-k10**2
	c0=k01**2*k11+k10**2*k22-k00*k11*k22

!	print*,"c2",c2,"c1",c1,"c0",c0


	p=c1-(c2**2)/3d0
	q=2*(c2**3)/27d0-(c1*c2)/3d0+c0
	r=(q**2)/4d0+(p**3)/27d0


!	print*,"p",p,"q",q,"r",r

	if(r>0)then
		a=-c2/3d0
		b=(-0.5d0*q+sqrt(r))
		c=(-0.5d0*q-sqrt(r))
!		print*,"A",a,"B",b,"C",c
		if(b<0)then
			b=-(-b)**(1d0/3d0)
		end if
		if(c<0)then
			c=-(-c)**(1d0/3d0)
		end if
!		print*,"A",a,"B",b,"C",c
		dm=a+b+c
!		dm=-c2/3d0+(-0.5d0*q+sqrt(r))**(0.333333d0)+(-0.5d0*q-sqrt(r))**(0.33333d0)
!		print*,"dm",dm
	else
		rho=sqrt(-p**3/27d0)
		phi=acos(-q/(2*rho))

		a=-c2/3d0+2*(rho**(1d0/3d0))*cos(phi/3d0)
		b=-c2/3d0+2*(rho**(1d0/3d0))*cos((phi+2*pi)/3d0)
		c=-c2/3d0+2*(rho**(1d0/3d0))*cos((phi+4*pi)/3d0)
		dm=min(a,b,c)

!		print*,"rho",rho,"phi",phi
!		print*,"a",a,"b",b,"c",c
!		print*,"dm",dm
	end if




	ka=-k10*cost/(k11-dm)+k01*sint/(k22-dm)
	kb=-k10*sint/(k11-dm)-k01*cost/(k22-dm)

!	print*,"ka",ka,"kb",kb


	ku=((1+kb**2)*xm-ka*kb*ym+ka*zm)/(1+ka**2+kb**2)
	kv=(-ka*kb*xm+(1+ka**2)*ym+kb*zm)/(1+ka**2+kb**2)
	kw=(ka*xm+kb*ym+(ka**2+kb**2)*zm)/(1+ka**2+kb**2)


!	print*,"vectortip",ku,kv,kw
!	print*,"centerpoint",xm,ym,zm

	vx=xm-ku;vy=ym-kv;vz=zm-kw

        d=1d0/sqrt(vx**2+vy**2+vz**2)  ! >>> Is 21-6-14 

        vx=vx*d ; vy=vy*d ; vz=vz*d    ! >>> Is 21-6-14 

end subroutine pol_physic

end module pola
