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




module biomechanic_pola	!>>>>>>>> by Miquel 17-5-13

use general
use genetic
use neighboring
use polarization
use shell       ! miguel4-11-13
use io          !!>>HC 13-10-2020 We need this to put a filter into the number of neightbors
!use nexus       ! Is 3-1-14

contains

subroutine iterdiferencial_pola
integer::nodmo,i,j,k,ii
real*8::a,b,c
real*8::ox,oy,oz !miguel4-1-13

  !CALCULATING FORCES AND MOVEMENT VECTORS
  call forces_pola !>> TT 24-7-2020
  !print*,"start iteration",getot
  !print*,"check*** p1",px(1),py(1),pz(1)
  !print*,"check*** p2",px(2),py(2),pz(2)
  !print*,"check*** p3",px(3),py(3),pz(3)
  !print*,getot,"maxval fmeanl",maxval(fmeanl(1:nd))


  a=0.0d0
  do i=1,nd
    if (dex(i)>a) then ; a=dex(i) ; ii=i ; end if
  end do
  
  if (a<epsilod) then ; delta=deltamin; return ; end if  !>>> Is 29-8-13 !>>> Is 6-7-21
  if(ffu(7)==0)then
    delta=resmax/a      !delta is adjusted so the node that has the greatest force
                      !will always move the same absolute distance, and the rest will
                      !move proportionally to that force
  else
    delta=deltamin
!    delta=deltamax 
  end if
                      
  if (delta>deltamax) delta=deltamax
  if (delta<deltamin) delta=deltamin
end subroutine iterdiferencial_pola

!***************************************************************************************************

subroutine rungekutta4_pola(d)  ! Runge-Kutta fourth order integration
real*8 d,halfd,sixthd
real*8 ox(nd),oy(nd),oz(nd)
real*8 kux(nd),kuy(nd),kuz(nd)
real*8 kdx(nd),kdy(nd),kdz(nd)
real*8 ktx(nd),kty(nd),ktz(nd)
real*8 kqx(nd),kqy(nd),kqz(nd)

halfd=d*0.5d0
sixthd=d/6.0d0

ox=node(:nd)%x ; oy=node(:nd)%y ; oz=node(:nd)%z 

!k1
kux=px(:nd) ; kuy=py(:nd) ; kuz=pz(:nd)  !!>> HC>> 8-7-2021 We assume that we have already called to forces once in iterdiferencial

!k2
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)
if (ffu(27)==1)then                         !!>> HC 12-7-2021
   if(nd>1) call neighbor_build_pola             !!>> HC 12-7-2021 
   if(ffu(24)==0) call restore_neighbors    !!>> HC 12-7-2021  Algorithm that restores lost ADD neighbors
   call forces_pola                              !!>> HC 12-7-2021  This is the normal version of forces that takes distances from dneigh         
else                                        !!>> HC 12-7-2021
   call forces_calculating_distances_pola        !!>> HC 12-7-2021  In this version we calculate the distances in forces
endif                                       !!>> HC 12-7-2021

kdx=px(:nd) ; kdy=py(:nd) ; kdz=pz(:nd)

!k3
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)
if (ffu(27)==1)then                         !!>> HC 12-7-2021
   if(nd>1) call neighbor_build_pola             !!>> HC 12-7-2021
   if(ffu(24)==0) call restore_neighbors    !!>> HC 12-7-2021 Algorithm that restores lost ADD neighbors
   call forces_pola                              !!>> HC 12-7-2021 This is the normal version of forces that takes distances from dneigh             
else                                        !!>> HC 12-7-2021
   call forces_calculating_distances_pola        !!>> HC 12-7-2021 In this version we calculate the distances in forces
endif                                       !!>> HC 12-7-2021

ktx=px(:nd) ; kty=py(:nd) ; ktz=pz(:nd)

!k4 
node(:nd)%x=node(:nd)%x+d*px(:nd)
node(:nd)%y=node(:nd)%y+d*py(:nd)
node(:nd)%z=node(:nd)%z+d*pz(:nd)
if (ffu(27)==1)then                         !!>> HC 12-7-2021
   if(nd>1) call neighbor_build_pola             !!>> HC 12-7-2021
   if(ffu(24)==0) call restore_neighbors    !!>> HC 12-7-2021 Algorithm that restores lost ADD neighbors
   call forces_pola                              !!>> HC 12-7-2021 This is the normal version of forces that takes distances from dneigh             
else                                        !!>> HC 12-7-2021
   call forces_calculating_distances_pola        !!>> HC 12-7-2021 In this version we calculate the distances in forces
endif                                       !!>> HC 12-7-2021

kqx=px(:nd) ; kqy=py(:nd) ; kqz=pz(:nd)

!final
node(:nd)%x=ox+sixthd*(kux+2*kdx+2*ktx+kqx)
node(:nd)%y=oy+sixthd*(kuy+2*kdy+2*kty+kqy)
node(:nd)%z=oz+sixthd*(kuz+2*kdz+2*ktz+kqz)
end subroutine

!***************************************************************************************************

subroutine adaptive_rungekutta_pola
real*8 ox(nd),oy(nd),oz(nd)
real*8 aux(nd),auy(nd),auz(nd)
real*8 adx(nd),ady(nd),adz(nd)
real*8 r,halfdelta,invdelta,mdx,mdy,mdz,suggesteddelta

37 continue

halfdelta=0.5d0*delta
invdelta=1d0/delta

ox=node(:nd)%x ; oy=node(:nd)%y ; oz=node(:nd)%z

call rungekutta4_pola(delta)

aux=node(:nd)%x ; auy=node(:nd)%y ; auz=node(:nd)%z
node(:nd)%x=ox  ; node(:nd)%y=oy  ; node(:nd)%z=oz

call rungekutta4_pola(halfdelta)
call rungekutta4_pola(halfdelta)

adx=node(:nd)%x ; ady=node(:nd)%y ; adz=node(:nd)%z

mdx=maxval(abs(aux-adx)*invdelta)
mdy=maxval(abs(auy-ady)*invdelta)
mdz=maxval(abs(auz-adz)*invdelta)

!mdx=maxval(abs(aux-adx)/ox)
!mdy=maxval(abs(auy-ady)/oy)
!mdz=maxval(abs(auz-adz)/oz)

if (mdx>=mdy.and.mdx>=mdz) r=mdx
if (mdy>=mdz.and.mdy>=mdz) r=mdy
if (mdz>=mdy.and.mdz>=mdx) r=mdz

suggesteddelta=0.9d0*prec*delta/r

if (r>prec) then !the step is too large
  delta=suggesteddelta
!print *,"NO",delta,r,prec
  goto 37
else
!print *,"SI",delta,r,prec
  ! in theory we should make delta=suggested delta but we prefer delta to be decided based on resmax in each step
  node(:nd)%x=adx
  node(:nd)%y=ady
  node(:nd)%z=adz
end if

end subroutine




!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCES

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine forces_pola
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,ideqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda,posca
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::nuve,inuve
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::switch,twoep

integer  ::lateral,vertical !flags that tell if there is a lateral or vertical component to take into account !>>Miquel28-1-14

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::ad,fd  !>>Miquel28-1-14

real*8   ::upr(nd)
integer  ::iupr(nd)
real*8   ::r(nd),er(nd)

real*8   ::rcilx(nd),rcily(nd),rcilz(nd)
real*8   ::rtorx(nd),rtory(nd),rtorz(nd)
real*8   ::rstorx(nd),rstory(nd),rstorz(nd)
real*8   ::rsprx,rspry,rsprz
real*8   ::arcilxch(nd),arcilych(nd),arcilzch(nd) !!>> HC 15-6-2020
real*8   ::rrcilxch(nd),rrcilych(nd),rrcilzch(nd) !!>> HC 15-6-2020

real*8   ::ftch                !>> HC 14-5-2020 Total adh-rep force to compare with maximum
!integer ::nveins(nd)         !
integer  ::epinveins(nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone
integer  ::whichend           ! For filtering too many neighbors !!>>HC 17-11-2020
character*300 ::cxhc          ! For filtering too many neighbors !!>>HC 17-11-2020
integer :: itt, mitt

  if (aut==0) then
    vcilx=0 ; vcily=0 ; vcilz=0 ; vsprx=0     !force vectors storage matrices, for different components
    vtorx=0 ; vtory=0 ; vtorz=0 ; vspry=0     !of the resulting force
    vstorx=0 ; vstory=0 ; vstorz=0 ; vsprz=0
  end if

  fmeanl=0 ; fmeanv=0  !storage vector that makes the balance between compressive and tensile forces within a node    !>>>Miquel23-1-14

  rcilx=0d0  ; rcily=0d0  ; rcilz=0d0  !they store the force components for all the nodes !>>Miquel4-4-14
  rtorx=0d0  ; rtory=0d0  ; rtorz=0d0
  rstorx=0d0 ; rstory=0d0 ; rstorz=0d0
  arcilxch=0d0; arcilych=0d0; arcilzch=0d0 !!>> HC 15-6-2020
  rrcilxch=0d0; rrcilych=0d0; rrcilzch=0d0 !!>> HC 15-6-2020
  epinveins=0
  
  do nod=1,nd
    !lonely=0 !fossile? >>Miquel27-2-14

    if (node(nod)%fix==2) then                ! >>> Is 30-6-14
      dex(nod)=0  !module of the force vector ! >>> Is 30-6-14
      px(nod)=0 ; py(nod)=0 ; pz(nod)=0       ! >>> Is 30-6-14
      cycle                                   ! >>> Is 30-6-14
    end if                                    ! >>> Is 30-6-14

    rsprx=0.0d0  ; rspry=0.0d0; rsprz=0.0d0
    ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
    tipi=node(nod)%tipus ; celi=node(nod)%icel
    iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)	
    rvx=0d0    ; rvy=0d0    ; rvz=0d0
    nuve=0!nuve=nveins(nod) !>>>>>Miquel 7-8-13 : this is not always zero since we fill this matrix from its neighbours
    switch=0
    alone=0

    !SPRINGS
    if(tipi<3)then
      iv=node(nod)%altre
      ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      !if(ffu(1)==1)then
        uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !the unit vector
        ddd=dd-node(nod)%eqs-node(iv)%eqs  !>>Miquel5-2-14
        f=2*node(nod)%hoo*ddd    !the force
        rsprx=f*uvx ; rspry=f*uvy ; rsprz=f*uvz
        fmeanv(nod)=ddd !>>>Miquel5-2-14
      !end if
    else
      iv=0    !>>>>Miquel17-1-14
    end if
    younod=node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=node(nod)%rep                                            !
    adhnod=node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=node(nod)%rec                                      !
    tornod=node(nod)%erp                                            !
    stornod=node(nod)%est                                          !
    reqnod=node(nod)%eqd                                            !
    nodda=node(nod)%add

    !NODE's REPULSIONS AND ADHESIONS

    do i=1,nneigh(nod)
      ic=neigh(nod,i)
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if
      
      do itt=1, nneigh(ic)
        if(neigh(ic,itt) .eq. nod)then
          mitt=itt
          exit
        end if
      end do

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      twoep=0
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares
      if (tipi<3)then
        if(tipic<3)then
          ivv=node(ic)%altre
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 
          !iudd=1.0d0/sqrt(icx**2+icy**2+icz**2)   !ic's spring vector	
          posca=icx*cx+icy*cy+icz*cz
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
              dotp=mcx*ccx+mcy*ccy+mcz*ccz
              ddd=abs(dotp)*md  !vertical component
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-ADDe(nod,i)-ADDe(ic,mitt)>epsilod .and. ffu(23)==1)cycle
              pesco=dotp*md**2
              a=1.0d0/ad              
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins(nod)=epinveins(nod)+1
              epinveins(ic)=epinveins(ic)+1
            else
              ! apical/basal contact, two epithelia from the same side
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=mcx**2+mcy**2+mcz**2 ! this should be sqrt(md) but we can do it later to make it faster
                                                                            !we take as spring vector the sum of ic's and nod's spring vectors
              if (md<epsilod)then !the two cylinders are parallel, the vector used is the spring !>>Miquel20-8-14
                 dotp=cx*ccx+cy*ccy+cz*ccz
                 if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                    ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                    a=ADDe(nod,i)+ADDe(ic,mitt)
                    if (ddd-a<epsilod) then
                       ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                       ddd=sqrt(ddd)              !lateral component
                       if (ad-a>epsilod) cycle
                       uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
                       fd=ddd     !distance used, vertical component
                       twoep=2
                       goto 300
                    end if
                    cycle 
                 end if  
              else  !proper apical cylindric interface applied !>>Miquel20-8-14
                 md=1d0/sqrt(md)
                 dotp=mcx*ccx+mcy*ccy+mcz*ccz
                 ad=abs(dotp)*md  !vertical component
                 a=ADDe(nod,i)+ADDe(ic,mitt)
                 if (ad-a<epsilod)then
                    ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                    ddd=sqrt(ddd)              !lateral component
                    if (ddd-a>epsilod) cycle
                    pesco=dotp*md**2
                    a=1/ddd                    
                    uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
                    fd=ddd     !distance used, vertical component
                    twoep=2
                    goto 300
                 end if
                    cycle 
              end if  
            end if
          else
            if (posca<0.0d0) then           
              ! we are from contrary sides so we must have face adhesion
              ! this is adhesion and repulsion between epithelial sides
              dotp=cx*ccx+cy*ccy+cz*ccz                     
              if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                a=ADDe(nod,i)+ADDe(ic,mitt)
                if (ddd-a<epsilod) then
                  ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                  goto 300
                end if
              end if
            end if
            cycle
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=(cx*ccx+cy*ccy+cz*ccz)!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=ADDe(nod,i)+ADDe(ic,mitt)
            if (ddd-a<epsilod) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        end if
      else
        if(tipic<3)then          ! IC IS EPITHELIAL and nod is not
          ivv=node(ic)%altre     ! we check the distance to radial distance in the plane of the ic cylinder 
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz   
          idd=1d0/sqrt(icx**2+icy**2+icz**2)
          dotp=icx*ccx+icy*ccy+icz*ccz  !hi ha un vector del revés, per tant això està al revés també
          if (dotp>0.0) then
            ddd=dotp*idd        !vertical component
            a=ADDe(nod,i)+ADDe(ic,mitt)
            if (ddd-a<0.0) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              uvx=icx*idd ; uvy=icy*idd ; uvz=icz*idd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        else
          fd=d !BOTH NODES ARE NON-EPITHELIAL: we just consider the interactions between nodes as such
          if (fd-ADDe(nod,i)-ADDe(ic,mitt)>epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
        end if
      end if

300   nuve=nuve+1

      if (nuve>mnn) then; 
        print *,"PANIC!!! PANIC!!!PANIC!!!; this is going to crash because there is too many " ; 
        print *,"neighbors per a node to interact, and vuvx matrix and those are too small"
        if (ffu(21)==1)then                                       !!>> HC 17-11-2020 We filter individuals that have too many neighbors
           if (ffufi(2,1)==1.and.ffufi(2,3)>0)then                !!>> HC 20-2-2021
              if (mod(getot,ffufi(2,3)).eq.0)then                 !!>> HC 20-2-2021   
                 whichend=2                                       !!>> HC 17-11-2020
                 goto 171                                         !!>> HC 17-11-2020
              endif                                               !!>> HC 20-2-2021
           endif                                                  !!>> HC 17-11-2020
        endif                                                     !!>> HC 20-2-2021
      end if                                                      !!>> HC 20-2-2021
  
      !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, fd, NOW WE CALCULATE THE ACTUAL ENERGIES
        if(node(nod)%icel==node(ic)%icel)then    !>> HC 12-6-2020 forces back to normal
          deqe=(EQDe(nod,i)+EQDe(ic,mitt))             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
            f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
          else                                   !>> HC 12-6-2020
            f=(younod+node(ic)%you)*(fd-deqe)                                                  !>> HC 27-7-2020
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(EQDe(nod,i)+EQDe(ic,mitt))             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
            f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
          else                                   !>> HC 12-6-2020
            adhe=0.5d0*(adhnod+node(ic)%adh)     !>>>> Miquel 16-8-13
            if (npag(1)>0) then                  !>> HC 12-6-2020! we have adhesion molecules
               do j=1,npag(1)                     !>> HC 12-6-2020
                  k=whonpag(1,j)                   !>> HC 12-6-2020
                  do kjjj=1,npag(1)                !>> HC 12-6-2020
                     jkkk=whonpag(1,kjjj)           !>> HC 12-6-2020
                     if (gex(nod,k)>0.0d0.and.gex(ic,jkkk)>0.0d0) then  !>> HC 12--6-2020  
                        adhe=adhe+gex(nod,k)*gex(ic,jkkk)*kadh(int(gen(k)%e(1)),int(gen(jkkk)%e(1))) ! >>> Is 7-6-14    !this is specific adhesion             !>> HC 12-6-2020
                     end if                         !>> HC 12-6-2020
                  end do                           !>> HC 12-6-2020
               end do                             !>> HC 12-6-2020
            end if                               !>> HC 12-6-2020
            f=2*adhe*(fd-deqe)                   !>> HC 12-6-2020 !in fact that is the derivative of the energy            
          end if                                 !>> HC 12-6-2020
        end if !;print*,"f",f,"fd",fd        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
       
        if (ffu(20)==0)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
           if (f>0.0d0)then                     !!>> HC 15-6-2020 adhesion and repulsion 
              arcilxch(nod)=arcilxch(nod)+f*uvx !!>> HC 15-6-2020 separately to test later 
              arcilych(nod)=arcilych(nod)+f*uvy !!>> HC 15-6-2020 whether adhesion has or not exeded the 
              arcilzch(nod)=arcilzch(nod)+f*uvz !!>> HC 15-6-2020 limit (maxad)
              arcilxch(ic)=arcilxch(ic)-f*uvx   !!>> HC 16-6-2020
              arcilych(ic)=arcilych(ic)-f*uvy   !!>> HC 16-6-2020
              arcilzch(ic)=arcilzch(ic)-f*uvz   !!>> HC 16-6-2020
           else                                 !!>> HC 15-6-2020
              rrcilxch(ic)=rrcilxch(ic)-f*uvx   !!>> HC 16-6-2020
              rrcilych(ic)=rrcilych(ic)-f*uvy   !!>> HC 16-6-2020
              rrcilzch(ic)=rrcilzch(ic)-f*uvz   !!>> HC 16-6-2020
              rrcilxch(nod)=rrcilxch(nod)+f*uvx !!>> HC 15-6-2020
              rrcilych(nod)=rrcilych(nod)+f*uvy !!>> HC 15-6-2020
              rrcilzch(nod)=rrcilzch(nod)+f*uvz !!>> HC 15-6-2020
           endif                                !!>> HC 15-6-2020
        else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
           rcilx(nod)=rcilx(nod)+f*uvx ; rcily(nod)=rcily(nod)+f*uvy ; rcilz(nod)=rcilz(nod)+f*uvz  !>>>Miquel 7-8-13
           rcilx(ic)=rcilx(ic)-f*uvx ; rcily(ic)=rcily(ic)-f*uvy ; rcilz(ic)=rcilz(ic)-f*uvz  !>>>Miquel 4-4-14
        endif                                  !!>> HC 15-6-2020
                   
        if (ffu(9)==1) then  ! 4-3-2020
          fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        else
          fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        end if

      !TORSION
      if(ffu(3)==0 .and.twoep==1 ) then !it is only between epithelial nodes
        !surface tension-like torsion (original)
        ! this is already calculated mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
        ! this is already calculated dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
        dotp=dotp*md  ! this goes instead of the previous section
        !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
        if(abs(dotp)-angletor*d>epsilod)then
          uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
          f=(stornod+node(ic)%est)*dotp         !>>>>>Miquel 7-8-13                         !!>> HC 6-7-2021 THIS IS VERTICAL TORSION (EST)
          a=f*uvx ; b=f*uvy ; c=f*uvz
          rstorx(nod)=rstorx(nod)+a ; rstory(nod)=rstory(nod)+b ; rstorz(nod)=rstorz(nod)+c
          rstorx(ic)=rstorx(ic)-a   ; rstory(ic)=rstory(ic)-b   ; rstorz(ic)=rstorz(ic)-c

          !parallel springs torsion   !new  !>>Miquel14-3-14
          dotp=(cx*ccx+cy*ccy+cz*ccz)*udd !vertical projection, more stable than the angle
          f=(tornod+node(ic)%erp)*dotp
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
          rtorx(nod)=rtorx(nod)+f*uvx ; rtory(nod)=rtory(nod)+f*uvy ; rtorz(nod)=rtorz(nod)+f*uvz

          !this is the neighbor, this force is not simmetrical, but it's no problem
          iudd=1.0d0/sqrt(icx**2+icy**2+icz**2)   !ic's spring vector	          
          dotp=-(icx*ccx+icy*ccy+icz*ccz)*iudd !vertical projection, more stable than the angle
          f=(tornod+node(ic)%erp)*dotp
          rtorx(ic)=rtorx(ic)-f*uvx ; rtory(ic)=rtory(ic)-f*uvy ; rtorz(ic)=rtorz(ic)-f*uvz
        end if
      end if
324   continue
    end do
    
    if (ffu(20)==0) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
      ftch=0
      ftch=sqrt(arcilxch(nod)**2+arcilych(nod)**2+arcilzch(nod)**2) !!>> HC 15-6-2020 modul adhesion
      b=maxad/ftch
      if (ftch>maxad) then !>>HC 14-05-2020 When the adhesion vector is larger than the limit(maxad)
	arcilxch(nod)=arcilxch(nod)*b !!>> HC 15-6-2020 !>>HC 14-05-2020 We rescale the vector 
	arcilych(nod)=arcilych(nod)*b !!>> HC 15-6-2020
	arcilzch(nod)=arcilzch(nod)*b !!>> HC 15-6-2020
      endif 
      rcilx(nod)=arcilxch(nod)+rrcilxch(nod) !!>> HC 15-6-2020 We have to add the adhesion and repulsion
      rcily(nod)=arcilych(nod)+rrcilych(nod) !!>> HC 15-6-2020 components 
      rcilz(nod)=arcilzch(nod)+rrcilzch(nod) !!>> HC 15-6-2020
    endif


    if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
      !fmeanl(nod)=0
    !else
      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    end if

789 if (ffu(4)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if

    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
      vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
      vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
      a=epinveins(nod) ! >>> Is 7-6-14
      a=1.0d0/a        ! >>> Is 7-6-14
      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
    end if
    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%fix==1)then
      uvx=node(nod)%orix-ix
      uvy=node(nod)%oriy-iy
      uvz=node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%kfi
      rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if

    dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector !!>> HC 14-7-2021
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  
return !X! THIS IS THE END         !!>> HC 17-11-2020

! REASONS TO STOP A SIMULATION: WHICHEND VALUES          !!>> HC 18-9-2020
!  1 = OVER  4 HOURS RUNNING                             !!>> HC  8-1-2021
!  2 = TOO MANY NEIGHBORS                                !!>> HC  7-1-2021
!  3 = TOO LARGE                                         !!>> HC 18-9-2020
!  4 = TOO MANY NODES                                    !!>> HC  7-1-2021
!  5 = NOT ENOUGH NODES                                  !!>> HC  7-1-2021
!  6 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
!  7 = TIME HAS ENDED (GETOT)                            !!>> HC 18-9-2020
!  8 = TIME HAS ENDED (REALTIME)                         !!>> HC 18-9-2020
!  9 = SOME NODES ARE DISPROPORTIONALY LARGE             !!>> HC 14-1-2020
! 10 = THERE ARE TOO MANY NODES IN AVERAGE PER BOX       !!>> HC 28-1-2021
! 11 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
! 12 = BROKEN EPITHELIUM                                 !!>> HC 18-9-2020
! 13 = BLACK HOLES                                       !!>> HC 18-9-2020

171  print *,""                                                                                                 !!>> HC 20-11-2021
  if (ffu(22)==0)then                                                                                           !!>> HC 20-11-2021
      print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)          !!>> HC 20-11-2021
      print *,""                                                                                                !!>> HC 20-11-2021
      print*,"reason to end:",whichend                                                                          !!>> HC 20-11-2021
  endif                                                                                                         !!>> HC 20-11-2021
  open(616, file=trim(carg)//".loga")                                                                           !!>> HC 20-11-2021 We want to know
      write(616,*) "reason to end:",whichend                                                                    !!>> HC 20-11-2021 why we killed this indv
  close(616)                                                                                                    !!>> HC 20-11-2021
  cxhc="pwd >> "//trim(carg)//".logb"                                                                           !!>> HC 20-11-2021 Save the name of the individual (evolution)
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  cxhc="paste "//trim(carg)//".logb "//trim(carg)//".loga > "//trim(carg)//".log"                               !!>> HC 20-11-2021 paste the name and the reason to end
  call system(cxhc)                                                                                             !!>> HC 20-11-2021 Remove the temporal files
  cxhc="rm "//trim(carg)//".logb"                                                                               !!>> HC 20-11-2021
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  cxhc="rm "//trim(carg)//".loga"                                                                               !!>> HC 20-11-2021
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  if (ffufi(whichend,2)==1)then                                                                                 !!>> HC 20-11-2021  If the end is lethal 
      cxhc='echo "0.00" > '//trim(carg)//'fitness'                                                              !!>> HC 20-11-2021 We kill the embryo 
      call system(cxhc)                                                                                         !!>> HC 20-11-2021 and it has no fitness 
  endif                                                                                                        

  call writesnap                                                                                                !!>> HC 20-11-2021  
  
  stop                                                                                                          !!>> HC 20-11-2021
45 print *,"end of file error"                                                                                  !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021
46 print *,"error in writing",trim(carg)//"t"                                                                   !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021
48 print *,"other end"                                                                                          !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021

  
end subroutine 



!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCES_CALCULATING_DISTANCES

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine forces_calculating_distances_pola
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,ideqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda,posca
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::nuve,inuve
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::switch,twoep

integer  ::lateral,vertical !flags that tell if there is a lateral or vertical component to take into account !>>Miquel28-1-14

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::ad,fd  !>>Miquel28-1-14

real*8   ::upr(nd)
integer  ::iupr(nd)
real*8   ::r(nd),er(nd)

real*8   ::rcilx(nd),rcily(nd),rcilz(nd)
real*8   ::rtorx(nd),rtory(nd),rtorz(nd)
real*8   ::rstorx(nd),rstory(nd),rstorz(nd)
real*8   ::rsprx,rspry,rsprz
real*8   ::arcilxch(nd),arcilych(nd),arcilzch(nd) !!>> HC 15-6-2020
real*8   ::rrcilxch(nd),rrcilych(nd),rrcilzch(nd) !!>> HC 15-6-2020

real*8   ::ftch                !>> HC 14-5-2020 Total adh-rep force to compare with maximum
!integer ::nveins(nd)         !
integer  ::epinveins(nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone
integer  ::whichend           ! For filtering too many neighbors !!>>HC 17-11-2020
character*300 ::cxhc          ! For filtering too many neighbors !!>>HC 17-11-2020
integer :: itt, mitt

  if (aut==0) then
    vcilx=0 ; vcily=0 ; vcilz=0 ; vsprx=0     !force vectors storage matrices, for different components
    vtorx=0 ; vtory=0 ; vtorz=0 ; vspry=0     !of the resulting force
    vstorx=0 ; vstory=0 ; vstorz=0 ; vsprz=0
  end if

  fmeanl=0 ; fmeanv=0  !storage vector that makes the balance between compressive and tensile forces within a node    !>>>Miquel23-1-14

  rcilx=0d0  ; rcily=0d0  ; rcilz=0d0  !they store the force components for all the nodes !>>Miquel4-4-14
  rtorx=0d0  ; rtory=0d0  ; rtorz=0d0
  rstorx=0d0 ; rstory=0d0 ; rstorz=0d0
  arcilxch=0d0; arcilych=0d0; arcilzch=0d0 !!>> HC 15-6-2020
  rrcilxch=0d0; rrcilych=0d0; rrcilzch=0d0 !!>> HC 15-6-2020
  epinveins=0
  
  do nod=1,nd
    !lonely=0 !fossile? >>Miquel27-2-14

    if (node(nod)%fix==2) then                ! >>> Is 30-6-14
      dex(nod)=0  !module of the force vector ! >>> Is 30-6-14
      px(nod)=0 ; py(nod)=0 ; pz(nod)=0       ! >>> Is 30-6-14
      cycle                                   ! >>> Is 30-6-14
    end if                                    ! >>> Is 30-6-14

    rsprx=0.0d0  ; rspry=0.0d0; rsprz=0.0d0
    ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
    tipi=node(nod)%tipus ; celi=node(nod)%icel
    iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)	
    rvx=0d0    ; rvy=0d0    ; rvz=0d0
    nuve=0!nuve=nveins(nod) !>>>>>Miquel 7-8-13 : this is not always zero since we fill this matrix from its neighbours
    switch=0
    alone=0

    !SPRINGS
    if(tipi<3)then
      iv=node(nod)%altre
      ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      !if(ffu(1)==1)then
        uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !the unit vector
        ddd=dd-node(nod)%eqs-node(iv)%eqs  !>>Miquel5-2-14
        f=2*node(nod)%hoo*ddd    !the force
        rsprx=f*uvx ; rspry=f*uvy ; rsprz=f*uvz
        fmeanv(nod)=ddd !>>>Miquel5-2-14
      !end if
    else
      iv=0    !>>>>Miquel17-1-14
    end if
    younod=node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=node(nod)%rep                                            !
    adhnod=node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=node(nod)%rec                                      !
    tornod=node(nod)%erp                                            !
    stornod=node(nod)%est                                          !
    reqnod=node(nod)%eqd                                            !
    nodda=node(nod)%add

    !NODE's REPULSIONS AND ADHESIONS

    do i=1,nneigh(nod)
      ic=neigh(nod,i)
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      do itt=1, nneigh(ic)
        if(neigh(ic,itt) .eq. nod)then
          mitt=itt
          exit
        end if
      end do

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      d=sqrt(ccx**2+ccy**2+ccz**2)                                                  !!>>HC 12-7-2021 THIS IS THE ONLY DIFFERENCE WITH FORCES
      !d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      twoep=0
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares
      if (tipi<3)then
        if(tipic<3)then
          ivv=node(ic)%altre
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 
          !iudd=1.0d0/sqrt(icx**2+icy**2+icz**2)   !ic's spring vector	
          posca=icx*cx+icy*cy+icz*cz
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
              dotp=mcx*ccx+mcy*ccy+mcz*ccz
              ddd=abs(dotp)*md  !vertical component
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-ADDe(nod,i)-ADDe(ic,mitt)>epsilod)cycle  !!>> HC 1-10-2021 
              pesco=dotp*md**2
              a=1.0d0/ad              
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins(nod)=epinveins(nod)+1
              epinveins(ic)=epinveins(ic)+1
            else
              ! apical/basal contact, two epithelia from the same side
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=mcx**2+mcy**2+mcz**2 ! this should be sqrt(md) but we can do it later to make it faster
                                                                            !we take as spring vector the sum of ic's and nod's spring vectors
              if (md<epsilod)then !the two cylinders are parallel, the vector used is the spring !>>Miquel20-8-14
                 dotp=cx*ccx+cy*ccy+cz*ccz
                 if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                    ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                    a=ADDe(nod,i)+ADDe(ic,mitt)
                    if (ddd-a<epsilod) then
                       ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                       ddd=sqrt(ddd)              !lateral component
                       if (ad-a>epsilod) cycle
                       uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
                       fd=ddd     !distance used, vertical component
                       twoep=2
                       goto 300
                    end if
                    cycle 
                 end if  
              else  !proper apical cylindric interface applied !>>Miquel20-8-14
                 md=1d0/sqrt(md)
                 dotp=mcx*ccx+mcy*ccy+mcz*ccz
                 ad=abs(dotp)*md  !vertical component
                 a=ADDe(nod,i)+ADDe(ic,mitt)
                 if (ad-a<epsilod)then
                    ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                    ddd=sqrt(ddd)              !lateral component
                    if (ddd-a>epsilod) cycle
                    pesco=dotp*md**2
                    a=1/ddd                    
                    uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
                    fd=ddd     !distance used, vertical component
                    twoep=2
                    goto 300
                 end if
                    cycle 
              end if  
            end if
          else
            if (posca<0.0d0) then           
              ! we are from contrary sides so we must have face adhesion
              ! this is adhesion and repulsion between epithelial sides
              dotp=cx*ccx+cy*ccy+cz*ccz                     
              if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                a=ADDe(nod,i)+ADDe(ic,mitt)
                if (ddd-a<epsilod) then
                  ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                  goto 300
                end if
              end if
            end if
            cycle
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=(cx*ccx+cy*ccy+cz*ccz)!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=ADDe(nod,i)+ADDe(ic,mitt)
            if (ddd-a<epsilod) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        end if
      else
        if(tipic<3)then          ! IC IS EPITHELIAL and nod is not
          ivv=node(ic)%altre     ! we check the distance to radial distance in the plane of the ic cylinder 
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz   
          idd=1d0/sqrt(icx**2+icy**2+icz**2)
          dotp=icx*ccx+icy*ccy+icz*ccz  !hi ha un vector del revés, per tant això està al revés també
          if (dotp>0.0) then
            ddd=dotp*idd        !vertical component
            a=ADDe(nod,i)+ADDe(ic,mitt)
            if (ddd-a<0.0) then
              ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>epsilod) cycle
              uvx=icx*idd ; uvy=icy*idd ; uvz=icz*idd  !unit vector vertical
              fd=ddd     !distance used, vertical component
            else
              cycle
            end if
          else
            cycle
          end if
        else
          fd=d !BOTH NODES ARE NON-EPITHELIAL: we just consider the interactions between nodes as such
          if (fd-ADDe(nod,i)-ADDe(ic,mitt)>epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
        end if
      end if

300   nuve=nuve+1

      if (nuve>mnn) then; 
        print *,"PANIC!!! PANIC!!!PANIC!!!; this is going to crash because there is too many " ; 
        print *,"neighbors per a node to interact, and vuvx matrix and those are too small"
        if (ffu(21)==1)then                                       !!>> HC 17-11-2020 We filter individuals that have too many neighbors
           if (ffufi(2,1)==1.and.ffufi(2,3)>0)then                !!>> HC 20-2-2021
              if (mod(getot,ffufi(2,3)).eq.0)then                 !!>> HC 20-2-2021   
                 whichend=2                                       !!>> HC 17-11-2020
                 goto 171                                         !!>> HC 17-11-2020
              endif                                               !!>> HC 20-2-2021
           endif                                                  !!>> HC 17-11-2020
        endif                                                     !!>> HC 20-2-2021
      end if                                                      !!>> HC 20-2-2021
  
      !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, fd, NOW WE CALCULATE THE ACTUAL ENERGIES
        if(node(nod)%icel==node(ic)%icel)then    !>> HC 12-6-2020 forces back to normal
          deqe=(EQDe(nod,i)+EQDe(ic,mitt))             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
            f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
          else                                   !>> HC 12-6-2020
            f=(younod+node(ic)%you)*(fd-deqe)                                                  !>> HC 27-7-2020
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(EQDe(nod,i)+EQDe(ic,mitt))             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
            f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
          else                                   !>> HC 12-6-2020
            adhe=0.5d0*(adhnod+node(ic)%adh)     !>>>> Miquel 16-8-13
            if (npag(1)>0) then                  !>> HC 12-6-2020! we have adhesion molecules
               do j=1,npag(1)                     !>> HC 12-6-2020
                  k=whonpag(1,j)                   !>> HC 12-6-2020
                  do kjjj=1,npag(1)                !>> HC 12-6-2020
                     jkkk=whonpag(1,kjjj)           !>> HC 12-6-2020
                     if (gex(nod,k)>0.0d0.and.gex(ic,jkkk)>0.0d0) then  !>> HC 12--6-2020  
                        adhe=adhe+gex(nod,k)*gex(ic,jkkk)*kadh(int(gen(k)%e(1)),int(gen(jkkk)%e(1))) ! >>> Is 7-6-14    !this is specific adhesion             !>> HC 12-6-2020
                     end if                         !>> HC 12-6-2020
                  end do                           !>> HC 12-6-2020
               end do                             !>> HC 12-6-2020
            end if                               !>> HC 12-6-2020
            f=2*adhe*(fd-deqe)                   !>> HC 12-6-2020 !in fact that is the derivative of the energy            
          end if                                 !>> HC 12-6-2020
        end if !;print*,"f",f,"fd",fd        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
       
        if (ffu(20)==0)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
           if (f>0.0d0)then                     !!>> HC 15-6-2020 adhesion and repulsion 
              arcilxch(nod)=arcilxch(nod)+f*uvx !!>> HC 15-6-2020 separately to test later 
              arcilych(nod)=arcilych(nod)+f*uvy !!>> HC 15-6-2020 whether adhesion has or not exeded the 
              arcilzch(nod)=arcilzch(nod)+f*uvz !!>> HC 15-6-2020 limit (maxad)
              arcilxch(ic)=arcilxch(ic)-f*uvx   !!>> HC 16-6-2020
              arcilych(ic)=arcilych(ic)-f*uvy   !!>> HC 16-6-2020
              arcilzch(ic)=arcilzch(ic)-f*uvz   !!>> HC 16-6-2020
           else                                 !!>> HC 15-6-2020
              rrcilxch(ic)=rrcilxch(ic)-f*uvx   !!>> HC 16-6-2020
              rrcilych(ic)=rrcilych(ic)-f*uvy   !!>> HC 16-6-2020
              rrcilzch(ic)=rrcilzch(ic)-f*uvz   !!>> HC 16-6-2020
              rrcilxch(nod)=rrcilxch(nod)+f*uvx !!>> HC 15-6-2020
              rrcilych(nod)=rrcilych(nod)+f*uvy !!>> HC 15-6-2020
              rrcilzch(nod)=rrcilzch(nod)+f*uvz !!>> HC 15-6-2020
           endif                                !!>> HC 15-6-2020
        else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
           rcilx(nod)=rcilx(nod)+f*uvx ; rcily(nod)=rcily(nod)+f*uvy ; rcilz(nod)=rcilz(nod)+f*uvz  !>>>Miquel 7-8-13
           rcilx(ic)=rcilx(ic)-f*uvx ; rcily(ic)=rcily(ic)-f*uvy ; rcilz(ic)=rcilz(ic)-f*uvz  !>>>Miquel 4-4-14
        endif                                  !!>> HC 15-6-2020
                   
        if (ffu(9)==1) then  ! 4-3-2020
          fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        else
          fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        end if

      !TORSION
      if(ffu(3)==0 .and.twoep==1 ) then !it is only between epithelial nodes
        !surface tension-like torsion (original)
        ! this is already calculated mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
        ! this is already calculated dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
        dotp=dotp*md  ! this goes instead of the previous section
        !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
        if(abs(dotp)-angletor*d>epsilod)then
          uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
          f=(stornod+node(ic)%est)*dotp         !>>>>>Miquel 7-8-13                         !!>> HC 6-7-2021 THIS IS VERTICAL TORSION (EST)
          a=f*uvx ; b=f*uvy ; c=f*uvz
          rstorx(nod)=rstorx(nod)+a ; rstory(nod)=rstory(nod)+b ; rstorz(nod)=rstorz(nod)+c
          rstorx(ic)=rstorx(ic)-a   ; rstory(ic)=rstory(ic)-b   ; rstorz(ic)=rstorz(ic)-c

          !parallel springs torsion   !new  !>>Miquel14-3-14
          dotp=(cx*ccx+cy*ccy+cz*ccz)*udd !vertical projection, more stable than the angle
          f=(tornod+node(ic)%erp)*dotp
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
          rtorx(nod)=rtorx(nod)+f*uvx ; rtory(nod)=rtory(nod)+f*uvy ; rtorz(nod)=rtorz(nod)+f*uvz

          !this is the neighbor, this force is not simmetrical, but it's no problem
          iudd=1.0d0/sqrt(icx**2+icy**2+icz**2)   !ic's spring vector	          
          dotp=-(icx*ccx+icy*ccy+icz*ccz)*iudd !vertical projection, more stable than the angle
          f=(tornod+node(ic)%erp)*dotp
          rtorx(ic)=rtorx(ic)-f*uvx ; rtory(ic)=rtory(ic)-f*uvy ; rtorz(ic)=rtorz(ic)-f*uvz
        end if
      end if
324   continue
    end do
    
    if (ffu(20)==0) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
      ftch=0
      ftch=sqrt(arcilxch(nod)**2+arcilych(nod)**2+arcilzch(nod)**2) !!>> HC 15-6-2020 modul adhesion
      b=maxad/ftch
      if (ftch>maxad) then !>>HC 14-05-2020 When the adhesion vector is larger than the limit(maxad)
	arcilxch(nod)=arcilxch(nod)*b !!>> HC 15-6-2020 !>>HC 14-05-2020 We rescale the vector 
	arcilych(nod)=arcilych(nod)*b !!>> HC 15-6-2020
	arcilzch(nod)=arcilzch(nod)*b !!>> HC 15-6-2020
      endif 
      rcilx(nod)=arcilxch(nod)+rrcilxch(nod) !!>> HC 15-6-2020 We have to add the adhesion and repulsion
      rcily(nod)=arcilych(nod)+rrcilych(nod) !!>> HC 15-6-2020 components 
      rcilz(nod)=arcilzch(nod)+rrcilzch(nod) !!>> HC 15-6-2020
    endif


    if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
      !fmeanl(nod)=0
    !else
      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    end if

789 if (ffu(4)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if

    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
      vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
      vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
      a=epinveins(nod) ! >>> Is 7-6-14
      a=1.0d0/a        ! >>> Is 7-6-14
      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
    end if

    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%fix==1)then
      uvx=node(nod)%orix-ix
      uvy=node(nod)%oriy-iy
      uvz=node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%kfi
      rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if

    dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector !!>> HC 14-7-2021
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  
return !X! THIS IS THE END         !!>> HC 17-11-2020

! REASONS TO STOP A SIMULATION: WHICHEND VALUES          !!>> HC 18-9-2020
!  1 = OVER  4 HOURS RUNNING                             !!>> HC  8-1-2021
!  2 = TOO MANY NEIGHBORS                                !!>> HC  7-1-2021
!  3 = TOO LARGE                                         !!>> HC 18-9-2020
!  4 = TOO MANY NODES                                    !!>> HC  7-1-2021
!  5 = NOT ENOUGH NODES                                  !!>> HC  7-1-2021
!  6 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
!  7 = TIME HAS ENDED (GETOT)                            !!>> HC 18-9-2020
!  8 = TIME HAS ENDED (REALTIME)                         !!>> HC 18-9-2020
!  9 = SOME NODES ARE DISPROPORTIONALY LARGE             !!>> HC 14-1-2020
! 10 = THERE ARE TOO MANY NODES IN AVERAGE PER BOX       !!>> HC 28-1-2021
! 11 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
! 12 = BROKEN EPITHELIUM                                 !!>> HC 18-9-2020
! 13 = BLACK HOLES                                       !!>> HC 18-9-2020

171  print *,""                                                                                                 !!>> HC 20-11-2021
  if (ffu(22)==0)then                                                                                           !!>> HC 20-11-2021
      print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)          !!>> HC 20-11-2021
      print *,""                                                                                                !!>> HC 20-11-2021
      print*,"reason to end:",whichend                                                                          !!>> HC 20-11-2021
  endif                                                                                                         !!>> HC 20-11-2021
  open(616, file=trim(carg)//".loga")                                                                           !!>> HC 20-11-2021 We want to know
      write(616,*) "reason to end:",whichend                                                                    !!>> HC 20-11-2021 why we killed this indv
  close(616)                                                                                                    !!>> HC 20-11-2021
  cxhc="pwd >> "//trim(carg)//".logb"                                                                           !!>> HC 20-11-2021 Save the name of the individual (evolution)
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  cxhc="paste "//trim(carg)//".logb "//trim(carg)//".loga > "//trim(carg)//".log"                               !!>> HC 20-11-2021 paste the name and the reason to end
  call system(cxhc)                                                                                             !!>> HC 20-11-2021 Remove the temporal files
  cxhc="rm "//trim(carg)//".logb"                                                                               !!>> HC 20-11-2021
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  cxhc="rm "//trim(carg)//".loga"                                                                               !!>> HC 20-11-2021
  call system(cxhc)                                                                                             !!>> HC 20-11-2021
  if (ffufi(whichend,2)==1)then                                                                                 !!>> HC 20-11-2021  If the end is lethal 
      cxhc='echo "0.00" > '//trim(carg)//'fitness'                                                              !!>> HC 20-11-2021 We kill the embryo 
      call system(cxhc)                                                                                         !!>> HC 20-11-2021 and it has no fitness 
  endif                                                                                                        

  call writesnap                                                                                                !!>> HC 20-11-2021  
  
  stop                                                                                                          !!>> HC 20-11-2021
45 print *,"end of file error"                                                                                  !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021
46 print *,"error in writing",trim(carg)//"t"                                                                   !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021
48 print *,"other end"                                                                                          !!>> HC 20-11-2021
  stop                                                                                                          !!>> HC 20-11-2021

  
end subroutine  


subroutine ordenarepeq_pola(ma,mt,rang)
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
end subroutine ordenarepeq_pola


end module biomechanic_pola
