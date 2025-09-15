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




module biomechanic	!>>>>>>>> by Miquel 17-5-13

use general
use genetic
use neighboring
use shell       ! miguel4-11-13
use io          !!>>HC 13-10-2020 We need this to put a filter into the number of neightbors
!use nexus       ! Is 3-1-14

real*8, private, allocatable   :: ab_pol(:,:)      !apical-basal polarity vector for each epi cell
real*8, private, allocatable   :: normal_vec(:,:)  !vector normal to the surface in each cell, it is the average of the neighbors

contains

subroutine iterdiferencial
integer::nodmo,i,j,k,ii
real*8::a,b,c
real*8::ox,oy,oz !miguel4-1-13

  call forces 

  a=0.0d0
  do i=1,nd
    if (dex(i)>a) then ; a=dex(i) ; ii=i ; end if
  end do

  if (a<epsilod) then ; rtime=rtime+delta ; delta=deltamin; return ; end if  !>>> Is 29-8-13
  if(ffu(12)==0)then
    delta=resmax/a      !delta is adjusted so the node that has the greatest force
                      !will always move the same absolute distance, and the rest will
                      !move proportionally to that force
  else
    delta=deltamin
!    delta=deltamax 
  end if

  if (delta>deltamax) delta=deltamax
  if (delta<deltamin) delta=deltamin
end subroutine iterdiferencial

!***************************************************************************************************

subroutine rungekutta4(d)  ! Runge-Kutta fourth order integration
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
kux=px(:nd) ; kuy=py(:nd) ; kuz=pz(:nd)

!k2
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)

if(nd>1) call neighbor_build_2021
!if(ffu(28)==1) then; call forces_grid; else; call forces; endif !!HC 8-6-2020
call forces !>> TT 24-7-2020

kdx=px(:nd) ; kdy=py(:nd) ; kdz=pz(:nd)

!k3
node(:nd)%x=node(:nd)%x+halfd*px(:nd)
node(:nd)%y=node(:nd)%y+halfd*py(:nd)
node(:nd)%z=node(:nd)%z+halfd*pz(:nd)

call neighbor_build_2021  !3-3-21 IMPORTANT, POTSER NO CALDRIA CRIDAR AIXO AQUI I AIXI ESTALVIEM MOLT DE TEMPS  

!if(ffu(28)==1) then; call forces_grid; else; call forces; endif !!HC 8-6-2020
call forces !>> TT 24-7-2020

ktx=px(:nd) ; kty=py(:nd) ; ktz=pz(:nd)

!k4 
node(:nd)%x=node(:nd)%x+d*px(:nd)
node(:nd)%y=node(:nd)%y+d*py(:nd)
node(:nd)%z=node(:nd)%z+d*pz(:nd)

call neighbor_build_2021 !3-3-21 IMPORTANT, POTSER NO CALDRIA CRIDAR AIXO AQUI I AIXI ESTALVIEM MOLT DE TEMPS  

!if(ffu(28)==1) then; call forces_grid; else; call forces; endif !!HC 8-6-2020
call forces !>> TT 24-7-2020

kqx=px(:nd) ; kqy=py(:nd) ; kqz=pz(:nd)

!final
node(:nd)%x=ox+sixthd*(kux+2*kdx+2*ktx+kqx)
node(:nd)%y=oy+sixthd*(kuy+2*kdy+2*kty+kqy)
node(:nd)%z=oz+sixthd*(kuz+2*kdz+2*ktz+kqz)

end subroutine

!***************************************************************************************************

subroutine adaptive_rungekutta
real*8 ox(nd),oy(nd),oz(nd)
real*8 aux(nd),auy(nd),auz(nd)
real*8 adx(nd),ady(nd),adz(nd)
real*8 r,halfdelta,invdelta,mdx,mdy,mdz,suggesteddelta

37 continue

halfdelta=0.5d0*delta
invdelta=1d0/delta

ox=node(:nd)%x ; oy=node(:nd)%y ; oz=node(:nd)%z

call rungekutta4(delta)

aux=node(:nd)%x ; auy=node(:nd)%y ; auz=node(:nd)%z
node(:nd)%x=ox  ; node(:nd)%y=oy  ; node(:nd)%z=oz

call rungekutta4(halfdelta)
call rungekutta4(halfdelta)

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

subroutine forces
  if (ffu(15)==1) then
    call forces2021
  else
    if (ffu(15)==0) then
      call forces_april
    else
      call forces_old  
    end if
  end if
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCE_ERRORS

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine force_errors 
character*300 ::cxhc          ! For filtering too many neighbors !!>>HC 17-11-2020 

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

  print *,""                                                                                                 !!>> HC 17-11-2020
  if (ffu(32)==0)then                                                                                           !!>> HC 4-2-2021
      print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)          !!>> HC 17-11-2020
      print *,""                                                                                                !!>> HC 17-11-2020
      print*,"reason to end:",whichend                                                                          !!>> HC 17-11-2020
  endif                                                                                                         !!>> HC 4-2-2021
  open(616, file=trim(carg)//".log")                                                                            !!>> HC 20-2-2021 We want to know
      write(616,*) "reason to end:",whichend                                                                    !!>> HC 17-11-2020 why we killed this indv
  close(616)                                                                                                    !!>> HC 17-11-2020
  if (ffufi(whichend,3)==1)then                                                                                 !!>> HC 20-2-2021  If the end is lethal 
      cxhc='echo "0.00" > '//trim(carg)//'fitness'                                                              !!>> HC 17-11-2020 We kill the embryo 
      call system(cxhc)                                                                                         !!>> HC 17-11-2020 and it has no fitness 
  endif                                                                                                        

  call writesnap                                                                                                 !!>> HC 17-11-2020  
  

  if (len_trim(carg)/=0) then                                                                                   !!>> HC 17-11-2020
    print *,trim(carg)                                                                                          !!>> HC 17-11-2020
    print *,trim(carg)//trim(nofi),"kii"                                                                        !!>> HC 17-11-2020
    open(23,file=trim(carg)//"t",iostat=i)                                                                      !!>> HC 17-11-2020
    print *,"making...",trim(carg)//"t"                                                                         !!>> HC 17-11-2020
    write(23,*,ERR=46) trim(carg)//trim(nofi)                                                                   !!>> HC 17-11-2020
    flush(23)                                                                                                   !!>> HC 17-11-2020

    open(23,file=trim(carg)//"t",iostat=ii)                                                                     !!>> HC 17-11-2020
    print *,trim(carg)//"t",ii,"iostat"                                                                         !!>> HC 17-11-2020
    read(23,*,END=45,ERR=48) cxhc                                                                               !!>> HC 17-11-2020
print *,""                                                                                                      !!>> HC 17-11-2020
print *,cxhc,"cx HERE"                                                                                          !!>> HC 17-11-2020
print *,""                                                                                                      !!>> HC 17-11-2020
    !close(23)                                                                                                  !!>> HC 17-11-2020
    call flush(23)                                                                                              !!>> HC 17-11-2020
    print *,"done with iostat",ii                                                                               !!>> HC 17-11-2020
    cxhc="ls -alrt "//trim(carg)//"t"                                                                           !!>> HC 17-11-2020
    call system(cxhc)                                                                                           !!>> HC 17-11-2020
  end if                                                                                                        !!>> HC 17-11-2020
  stop                                                                                                          !!>> HC 17-11-2020
45 print *,"end of file error"                                                                                  !!>> HC 17-11-2020
  stop                                                                                                          !!>> HC 17-11-2020
46 print *,"error in writing",trim(carg)//"t"                                                                   !!>> HC 17-11-2020
  stop                                                                                                          !!>> HC 17-11-2020
48 print *,"other end"                                                                                          !!>> HC 17-11-2020
  stop                                                                                                          !!>> HC 17-11-2020
end subroutine 


!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCES2021

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine forces2021  ! CALCULATES FORCES BETWEEN NODES IN THE EPI_GRID, THERE IS ALWAYS ADHESION BETWEEN NODES CONNECTED IN THE GRID
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,ideqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::twoep

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::ad,fd  !>>Miquel28-1-14

real*8   ::rcilx(nd),rcily(nd),rcilz(nd)
real*8   ::rtorx(nd),rtory(nd),rtorz(nd)
real*8   ::rstorx(nd),rstory(nd),rstorz(nd)
real*8   ::rsprx,rspry,rsprz
real*8   ::arcilxch(nd),arcilych(nd),arcilzch(nd) !!>> HC 15-6-2020
real*8   ::rrcilxch(nd),rrcilych(nd),rrcilzch(nd) !!>> HC 15-6-2020

real*8   ::ftch                !>> HC 14-5-2020 Total adh-rep force to compare with maximum
integer  ::epinveins(nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone
integer  ::whichend           ! For filtering too many neighbors !!>>HC 17-11-2020

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

    if (node(nod)%fix==2) then
      dex(nod)=0  !module of the force vector ! >>> Is 30-6-14
      px(nod)=0 ; py(nod)=0 ; pz(nod)=0       ! >>> Is 30-6-14
      cycle                                   ! >>> Is 30-6-14
    end if                                    ! >>> Is 30-6-14

    rsprx=0.0d0  ; rspry=0.0d0; rsprz=0.0d0
    ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
    tipi=node(nod)%tipus ; celi=node(nod)%icel
    iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)	
    rvx=0d0    ; rvy=0d0    ; rvz=0d0
    switch=0
    alone=0

    !SPRINGS
    if(tipi<3)then
      iv=node(nod)%altre
      ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !the unit vector
      ddd=dd-node(nod)%eqs-node(iv)%eqs  !>>Miquel5-2-14
      f=2*node(nod)%hoo*ddd    !the force
      rsprx=f*uvx ; rspry=f*uvy ; rsprz=f*uvz
      fmeanv(nod)=ddd !>>>Miquel5-2-14
    end if
    younod=node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=node(nod)%rep                                            !
    adhnod=node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=node(nod)%rec                                      !
    !tornod=node(nod)%erp                                            !
    !stornod=node(nod)%est                                          !
    reqnod=node(nod)%eqd                                            !
    nodda=node(nod)%add

    !NODE's REPULSIONS AND ADHESIONS
    do i=1,epi_ngrid(nod)

      ic=epi_grid(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact

      if (nod>=ic) cycle
      !if (nod==iv) cycle
      !if (ic==iv) cycle
      if (node(ic)%fix/=0)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if
      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      d=sqrt(ccx**2+ccy**2+ccz**2)
      ud=1d0/d
      tipic=node(ic)%tipus
      fd=d 
      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the FORCE

      if(celi==node(ic)%icel)then   
        deqe=(reqnod+node(ic)%eqd)             
        if(fd-deqe<-epsilod)then 	        			
          f=(fd-deqe)*(repnod+node(ic)%rep)    
        else                                   
          f=(younod+node(ic)%you)*(fd-deqe) 
        end if                                 !>> HC 12-6-2020
      else                                     !>> HC 12-6-2020
        deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
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
        end if 
      end if

      if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
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

      ! REVISAR??  if (ffu(15)==3) then  ! 4-3-2020
      !    fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
      !  else
      !    fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
      !  end if
!print *,nod,ic,rcilx(nod),rcily(nod),rcilz(nod)          
    end do
goto 199
    do i=1,nneigh(nod)  ! this line the next and seven lines after are the only differences with the forces_epigrid
      ic=neigh(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact

      if (ic==node(ic)%altre) cycle
      if (nod>=ic) cycle
!      if(nod>ic.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
      if(node(ic)%fix/=0)then
        alone=1
        cycle
      end if
      
      if (co_griders(nod,ic)==1) cycle  ! these already interacted through the grid

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      fd=d 
      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the FORCE

      if (node(ic)%tipus<3) then  ! IMPORTANT WE ASSUME, FOR SIMPLICITY AND NUMERICAL EASINESS, THAT EPITHELIAL CELLS CAN 
                                  ! REPULSE EACH OTHER OUTSIDE THE GRID BUT NOT ADHERE (I.E. NO FACE ADHESION BETWEEN EPITHELIA
        if(celi==node(ic)%icel)then   
          deqe=(reqnod+node(ic)%eqd)             
          if(fd-deqe<-epsilod)then 	        			
            f=(fd-deqe)*(repnod+node(ic)%rep)    
          else                                   
            cycle 
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
            f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
          else                                   !>> HC 12-6-2020
            cycle
          end if 
        end if
      else
        if(celi==node(ic)%icel)then   
          deqe=(reqnod+node(ic)%eqd)             
          if(fd-deqe<-epsilod)then 	        			
            f=(fd-deqe)*(repnod+node(ic)%rep)    
          else                                   
            f=(younod+node(ic)%you)*(fd-deqe) 
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
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
          end if 
        end if
      end if

      if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
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
      endif

      ! REVISAR??  if (ffu(15)==3) then  ! 4-3-2020
      !    fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
      !  else
      !    fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
      !  end if
    end do

199 if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
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

    ! REVISAR if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
    !  !fmeanl(nod)=0
    !!else
!      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    !end if

    ! REVISAR if (ffu(8)==1) then   !this is kind of crappy in the sense that 
    !  if (alone==0) then
    !    node(nod)%talone=node(nod)%talone+1
    !  else
    !    node(nod)%talone=0
    !  end if
    !end if

    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
      !vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
      !vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

    !if (epinveins(nod)>0) then  ! >>> Is 7-6-14
    !  a=epinveins(nod) ! >>> Is 7-6-14
    !  a=1.0d0/a        ! >>> Is 7-6-14
    !  rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
    !  rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
    !  rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
    !else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
    !end if
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

    if (aut==0) dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  return
171 call force_errors    
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   REP_ADH_EPI_GRID2

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine rep_adh_epi_grid2  ! CALCULATES FORCES BETWEEN NODES IN THE EPI_GRID, THERE IS ALWAYS ADHESION BETWEEN NODES CONNECTED IN THE GRID
real*8   ::ix,iy,iz,dd,hdd
real*8   ::a,b,c,d,e,f,g,icd,uicd
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,deki,deda
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod,nodda    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz,pvx,pvy,pvz,ux,uy,uz
real*8   ::ud,u_long_c,uddd,dotpro,dist_forced
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco,rccx,rccy,rccz,rcx,rcy,rcz,ru_long_c
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::twoep

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::projection

real*8   ::ad,lat_dist,up_dist,dist_force,aup_dist

real*8   ::rcilx(nd),rcily(nd),rcilz(nd)
real*8   ::rtorx(nd),rtory(nd),rtorz(nd)
real*8   ::rsprx,rspry,rsprz
real*8   ::arcilxch(nd),arcilych(nd),arcilzch(nd) !!>> HC 15-6-2020
real*8   ::rrcilxch(nd),rrcilych(nd),rrcilzch(nd) !!>> HC 15-6-2020

real*8   ::ftch                !>> HC 14-5-2020 Total adh-rep force to compare with maximum
integer  ::epinveins(nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone
integer  ::whichend           ! For filtering too many neighbors !!>>HC 17-11-2020

  if (aut==0) then
    vcilx=0 ; vcily=0 ; vcilz=0 ; vsprx=0     !force vectors storage matrices, for different components
    vtorx=0 ; vtory=0 ; vtorz=0 ; vspry=0     !of the resulting force
    vstorx=0 ; vstory=0 ; vstorz=0 ; vsprz=0
  end if

  fmeanl=0 ; fmeanv=0  !storage vector that makes the balance between compressive and tensile forces within a node    !>>>Miquel23-1-14

  rcilx=0.0d0  ; rcily=0.0d0  ; rcilz=0.0d0  !they store the force components for all the nodes !>>Miquel4-4-14
  rtorx=0.0d0  ; rtory=0.0d0  ; rtorz=0.0d0
  arcilxch=0d0; arcilych=0d0; arcilzch=0d0 !!>> HC 15-6-2020
  rrcilxch=0d0; rrcilych=0d0; rrcilzch=0d0 !!>> HC 15-6-2020
  epinveins=0
  
  do nod=1,nd

    if (node(nod)%fix==2) then
      dex(nod)=0  !module of the force vector ! >>> Is 30-6-14
      px(nod)=0 ; py(nod)=0 ; pz(nod)=0       ! >>> Is 30-6-14
      cycle                                   ! >>> Is 30-6-14
    end if                                    ! >>> Is 30-6-14

    ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
    tipi=node(nod)%tipus ; celi=node(nod)%icel
    !rvx=0d0    ; rvy=0d0    ; rvz=0d0
    alone=0

    !SPRINGS
    if (tipi<3) then
      iv=node(nod)%altre
      ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx*cx+cy*cy+cz*cz)
      hdd=dd*0.5d0
      u_long_c=1d0/dd
      ddd=dd-node(nod)%eqs-node(iv)%eqs  !>>Miquel5-2-14
      f=2*node(nod)%hoo*ddd              !the force
      rsprx=f*cx ; rspry=f*cy ; rsprz=f*cz
      fmeanv(nod)=ddd !>>>Miquel5-2-14
      ux=cx*u_long_c ; uy=cy*u_long_c ; uz=cz*u_long_c ! we will need this later
    end if
    younod=node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=node(nod)%rep                                            !
    adhnod=node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=node(nod)%rec                                      !
    reqnod=node(nod)%eqd                                            !
    nodda=node(nod)%add

    !NODE's REPULSIONS AND ADHESIONS
    do i=1,epi_ngrid(nod)

      ic=epi_grid(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact

      if (nod==ic) cycle
      if (ic==iv)  cycle
!      if(nod>=ic) cycle    ! we force symmetry of the interactions fij=fji
      if (node(ic)%fix/=0)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      if(celi==node(ic)%icel)then   
        deki=repnod+node(ic)%rep
        adhe=younod+node(ic)%you
      else
        deki=repcelnod+node(ic)%rec
        adhe=0.5d0*(adhnod+node(ic)%adh)             !>>>> Miquel 16-8-13
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
      end if

      ! basic variables of ic
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      d=sqrt(ccx*ccx+ccy*ccy+ccz*ccz)
      ud=1d0/d
      tipic=node(ic)%tipus
      dist_force=d 
      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      ! the relevant node properties 
      deqe=reqnod+node(ic)%eqd

      ! the lateral and up distances
      dotpro=ccx*cx+ccy*cy+ccz*cz                           
      up_dist=dotpro*u_long_c                              ! distance normal on the face: is the projection of the cc vector 
                                                           ! into the cylinder vector c       
      lat_dist=sqrt(d*d-up_dist*up_dist)                   ! this is by pythagoras the magnitude of the lateral force                  
      deda=nodda+node(ic)%add
      if (dotpro>=0) then ! UPPER FACE: ic is in the upper face of nod
        e=deqe
        g=deda
      else      ! LOWER FACE: ic is down in the cylinder side of nod, so we can interact until half the cylinder length 
        e=hdd
        g=hdd
      end if  

      ! FOR repulsion ic needs to have both a lateral and up distance smaller than eqd1+eqd2 -> lateral repulsion, 
      ! if not adhesion
      ! FOR lateral adhesion both the lateral and up distance needs to be smaller than add1+add2  -> lateral adhesion,
      ! if not radial adhesion

      a=lat_dist-deqe                   
      if (lat_dist<deqe.and.up_dist<e) then ! lateral repulsion
        f=a*deki
      else
!        if (lat_dist<deda.and.up_dist<g) then ! lateral atraction ! A FER 24-3-21 que quan d=eqd la forca no canvii de cop
!          f=a*adhe
!        else        
          f=adhe*(d-deqe)  ! d is not a mistake 
          ud=1d0/d
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud          
print *,nod,ic,"fora",f        
          goto 33
!        end if
      end if                     

      b=1.0d0/lat_dist 
      uvx=(ccx-ux*up_dist)*b ; uvy=(ccy-uy*up_dist)*b ; uvz=(ccz-uz*up_dist)*b 
                                                       ! ux*up_dist is the projection vector of cc into c
                                                       ! and ccx-cx*c is the rest of these two vectors that is in
                                                       ! in fact the lateral component in the cylinder
                                                       ! then we multiply by b to make it unit vector
33    a=f*uvx ; b=f*uvy ; c=f*uvz
      if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
        if (f>0.0d0)then                     
          arcilxch(nod)=arcilxch(nod)+a ; arcilych(nod)=arcilych(nod)+b ; arcilzch(nod)=arcilzch(nod)+c 
!          arcilxch(ic)=arcilxch(ic)-a   ; arcilych(ic)=arcilych(ic)-b   ; arcilzch(ic)=arcilzch(ic)-c   !!>> HC 16-6-2020
        else                                
          rrcilxch(ic)=rrcilxch(ic)-a   ; rrcilych(ic)=rrcilych(ic)-b   ; rrcilzch(ic)=rrcilzch(ic)-c   !!>> HC 16-6-2020
!          rrcilxch(nod)=rrcilxch(nod)+a ; rrcilych(nod)=rrcilych(nod)+b ; rrcilzch(nod)=rrcilzch(nod)+c !!>> HC 15-6-2020
        endif                                
      else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
        rcilx(nod)=rcilx(nod)+a ; rcily(nod)=rcily(nod)+b ; rcilz(nod)=rcilz(nod)+c  !>>>Miquel 7-8-13
!        rcilx(ic)=rcilx(ic)-a   ; rcily(ic)=rcily(ic)-b   ; rcilz(ic)=rcilz(ic)-c  !>>>Miquel 4-4-14
      endif                                 

      ! REVISAR??  if (ffu(15)==3) then  ! 4-3-2020
      !    fmeanl(nod)=fmeanl(nod)+(dist_force-deqe) ; fmeanl(ic)=fmeanl(ic)+(dist_force-deqe) ! >>> Is 21-6-14 
      !  else
      !    fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
      !  end if
    end do

19  if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
      ftch=0
      ftch=sqrt(arcilxch(nod)**2+arcilych(nod)**2+arcilzch(nod)**2) !!>> HC 15-6-2020 modul adhesion
      b=maxad/ftch
      if (ftch>maxad) then !>>HC 14-05-2020 When the adhesion vector is larger than the limit(maxad)
	arcilxch(nod)=arcilxch(nod)*b ; arcilych(nod)=arcilych(nod)*b ; arcilzch(nod)=arcilzch(nod)*b !!>> HC 15-6-2020
      endif 
      rcilx(nod)=arcilxch(nod)+rrcilxch(nod) ; rcily(nod)=arcilych(nod)+rrcilych(nod) ; rcilz(nod)=arcilzch(nod)+rrcilzch(nod)
    endif

    ! REVISAR if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
    !  !fmeanl(nod)=0
    !!else
!      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    !end if

    ! REVISAR if (ffu(8)==1) then   !this is kind of crappy in the sense that 
    !  if (alone==0) then
    !    node(nod)%talone=node(nod)%talone+1
    !  else
    !    node(nod)%talone=0
    !  end if
    !end if

    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
      vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is display
      !vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

    !if (epinveins(nod)>0) then  ! >>> Is 7-6-14
    !  a=epinveins(nod) ! >>> Is 7-6-14
    !  a=1.0d0/a        ! >>> Is 7-6-14
    !  rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
    !  rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
    !  rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
    !else
print *,rcilx(nod),rcily(nod),rcilz(nod),rtorx(nod),rtory(nod),rtorz(nod),"ioo"    
      rvx=rsprx+rcilx(nod)+rtorx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)+rtory(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)+rtorz(nod)  ! >>> Is 7-6-14
    !end if
    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%fix==1)then
      uvx=node(nod)%orix-ix
      uvy=node(nod)%oriy-iy
      uvz=node(nod)%oriz-iz
      d=uvx*uvx+uvy*uvy+uvz*uvz ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%kfi
      rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if

    if (aut==0) dex(nod)=sqrt(rvx*rvx+rvy*rvy+rvz*rvz)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
end subroutine 


!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCES_OLD

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine forces_old
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe
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

      ic=neigh(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      if (tipi<3)then
        if(tipic<3)then
          ivv=node(ic)%altre
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 

          mcx=icx+cx; mcy=icy+cy; mcz=icz+cz ! this is the sum of the cylinder vectors of i and ic
          md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  ! we do that so that we can calculate what happens to i and ic and the same time
          dotp=mcx*ccx+mcy*ccy+mcz*ccz       ! dot product between mc and cc (the vector from i to ic)
          ddd=dotp*md  !vertical component   ! this is the projection of cc into mc  see Fig 2
                   
          posca=icx*cx+icy*cy+icz*cz      ! Is 3-3-21 THIS is FIG.1!!!!! 
          
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              !!! Is 3-3-21 FIG.2  
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-nodda-node(ic)%add>epsilod) cycle !>>HC ISAAC 14-5-2020 !>> HC 12-6-2020
              pesco=ddd*md   ! Is 3-3-21 IT WOULD BETTER BE pesco=ddd*md, it would be the same by faster
              a=1.0d0/ad              
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins(nod)=epinveins(nod)+1
              epinveins(ic)=epinveins(ic)+1
            else
              ! SO THIS IS APICAL-APICAL OR BASAL-BASAL INTERACTION but by cells that are not close in the plane of the epi
                                      ! FIG 3 IS 3-3-21 FIG 3
              ddd=abs(ddd)  !vertical component              
              a=nodda+node(ic)%add
              if(ddd-a<epsilod) cycle  ! ic is too far away in the vertical component
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component                              
              if (ad-a>epsilod) cycle  ! ic is too far away in the lateral component
              pesco=dotp*md*md  ! REVISAR
              a=1/ddd                    
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ddd     !distance used, vertical component
              twoep=2
            end if
          else
            if (posca<epsilod) then           
              ! here there is an apical interacting with a basal or vice versa
              ! this is adhesion and repulsion between epithelial sides
              if (dotp<epsilod) then         
                ddd=-dotp*udd                ! this is projection of vector cc into the cylinder vector of i (the direction UP face)
                a=nodda+node(ic)%add         ! this is still figure 3
                if (ddd-a<epsilod) then
                  ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                else
                  cycle
                end if
              end if
            end if
            cycle  ! Is 3-3-21 This never happens??? since it would be detected before with posca???
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=cx*ccx+cy*ccy+cz*ccz!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=nodda+node(ic)%add
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
            a=nodda+node(ic)%add
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
          if (fd-nodda-node(ic)%add>epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
        end if
      end if

300   nuve=nuve+1

      !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, fd, NOW WE CALCULATE THE ACTUAL ENERGIES
      !if (ffu(2)==1) then IS 25-12-13
        if(node(nod)%icel==node(ic)%icel)then    !>> HC 12-6-2020 forces back to normal
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
            f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
          else                                   !>> HC 12-6-2020
            if(ffu(28)==1 .and. node(nod)%tipus .le. 2 .and. node(ic)%tipus .le. 2) then !>> HC 27-7-2020
              ! We use the grid and both nodes are epithelial.
              if(any(ic .eq. epi_grid(nod, :epi_ngrid(nod)))) then !>> HC 27-7-2020
                f=(younod+node(ic)%you)*(fd-deqe)                                   !>> HC 27-7-2020
              else                                                 !>> HC 27-7-2020
                f=0.0d0                                            !>> HC 27-7-2020
              endif                                                !>> HC 27-7-2020
            else                                                   !>> HC 27-7-2020
              f=(younod+node(ic)%you)*(fd-deqe) 
            endif                                                  !>> HC 27-7-2020
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
            f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
          else                                   !>> HC 12-6-2020
            if(ffu(28)==1 .and. node(nod)%tipus .le. 2 .and. node(ic)%tipus .le. 2) then !>> TT 24-7-2020
              ! We use the grid and both nodes are epithelial.
              if(any(ic .eq. epi_grid(nod, :epi_ngrid(nod)))) then !>> TT 24-7-2020
                ! Nodes can adhere to each other only if they are connected in the grid
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
              else                                   !>> HC 27-7-2020
                f=0.0d0                              !>> HC 27-7-2020
              end if
            else !>> TT 24-7-2020 ffu(28)==0 or both nodes are not epithelial
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
            end if
          end if                                 !>> HC 12-6-2020
        end if !;print*,"f",f,"fd",fd        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
       
        if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
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
                   
        if (ffu(15)==3) then  ! 4-3-2020
          fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        else
          fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        end if

      !TORSION
      if(ffu(4)==0 .and.twoep==1 ) then !it is only between epithelial nodes
        !surface tension-like torsion (original)
        ! this is already calculated mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
        ! this is already calculated dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
        dotp=dotp*md  ! this goes instead of the previous section
        !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
        if(abs(dotp)-angletor*d>epsilod)then
          uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
          f=(stornod+node(ic)%est)*dotp ! NOTICE %EST AIXO ES UNA FORCA CAP AL PROMIG i en la direccio d'aquest >Miquel 7-8-13
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

    if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
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

789 if (ffu(8)==1) then   !this is kind of crappy in the sense that 
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

    if (aut==0) dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  return
171 call force_errors    
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   FORCES_OLDNEW

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine forces_april

  if (ffu(28)==1) then
    call rep_adh_epi_grid     !version using epi_grid  
    if (ffu(4)==0) then 
      call calculate_normals
       call torsion_classic_epi_grid
    end if    
  else
    call rep_adh_no_epi_grid   !version not using epi_grid but neigh 
    if (ffu(4)==0) then 
      call calculate_normals
      call torsion_classic_no_epi_grid
    end if    
  end if
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!
!                                      SUBROUTINE CALCULATE NORMALS
!
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
subroutine calculate_normals
integer i,iv,j,k,tipi
real*8 cx,cy,cz,dd,udd,a,b,c

if (allocated(ab_pol)) deallocate(ab_pol)
allocate(ab_pol(nd,3))     !apical-basal polarity vector for each epi cell
if (allocated(normal_vec)) deallocate(normal_vec)
allocate(normal_vec(nd,3))  !vector normal to the surface in each cell, it is the average of the neighbors

  !first we calculate the ab-pol vectors for each cell
  do i=1,nd
    if (node(i)%tipus<3) then ! ATENCIO, HAURIEM DE FER UNA LLISTA AMB TOTS ELS NODES QUE SON EPITELIALS, I ACTUALITZARLA 
      iv=node(i)%altre
      cx=node(iv)%x-node(i)%x        ; cy=node(iv)%y-node(i)%y         ; cz=node(iv)%z-node(i)%z
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      ab_pol(i,1)=cx*udd ; ab_pol(i,2)=cy*udd ; ab_pol(i,3)=cz*udd  !the unit vector            
    end if
  end do

  ! now we calculate the normal surface to each cell by summing the ab_pol of its epigrid neighbors
  do i=1,nd 
    a=0.0d0 ; b=0.0d0 ; c=0.0d0
    tipi=node(i)%tipus
    do j=1,nneigh(i)
      k=neigh(i,j)
      if (tipi/=node(k)%tipus) cycle
      a=a+ab_pol(k,1) ; b=b+ab_pol(k,2) ; c=c+ab_pol(k,3)
    end do
    dd=sqrt(a**2+b**2+c**2)
    udd=1d0/dd
    normal_vec(i,1)=a*udd ; normal_vec(i,2)=b*udd ; normal_vec(i,3)=c*udd  !the unit vector            
  end do

end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!
!                                      TORSION_classic_epi_grid
!
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine torsion_classic_epi_grid
integer i,j,k,ic,tipi
real*8 a,b,c,cx,cy,cz,ccx,ccy,ccz,f,dotp,tornod,stornod,ix,iy,iz,mcx,mcy,mcz,icx,icy,icz,md,uvx,uvy,uvz,udd,ud,iudd,d
real*8 aa,bb,cc
real*8 :: proych, michx, michy, michz, dotpich, uchx, uchy, uchz, modch

do i=1,nd
  tipi=node(i)%tipus
  if (tipi>2) cycle
  f=0.0d0
  tornod=node(i)%erp
  stornod=node(i)%est
  cx=ab_pol(i,1) ; cy=ab_pol(i,2) ; cz=ab_pol(i,3)  
  ix=node(i)%x   ; iy=node(i)%y   ; iz=node(i)%z  
  a=0.0d0  ; b=0.0d0  ; c=0.0d0
  aa=0.0d0 ; bb=0.0d0 ; cc=0.0d0
  do j=1,epi_ngrid(i)
    ic=epi_grid(i,j)
    if (node(ic)%tipus/=tipi) cycle   
    !surface tension-like torsion (original)
    ccx=node(ic)%x-ix ; ccy=node(ic)%y-iy ; ccz=node(ic)%z-iz
    icx=ab_pol(ic,1) ; icy=ab_pol(ic,2) ; icz=ab_pol(ic,3)
    mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
    
    !!!!!!!!!!!!!NEW TORSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                   !!>> HC 17-3-2021 
    !c  nod spring vector                                                                           !!>> HC 26-3-2021
    !ic ic spring vector                                                                            !!>> HC 26-3-2021
    !mc sum of spring vectors (ic+c)                                                                !!>> HC 26-3-2021
    !cc nod to ic vector                                                                            !!>> HC 26-3-2021
       
    dotp=cx*ccx+cy*ccy+cz*ccz    !!>> Is-13-4-2021
           
    !!VERTICAL TORSION                                                          !!>> HC 26-3-2021
    if (abs(dotp)>epsilod)then                                                  !!>> HC 26-3-2021
       mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)    !!>> HC 26-3-2021
       proych=(mcx*ccx+mcy*ccy+mcz*ccz)*md*md             !!>> HC 26-3-2021 proyection (scalar) of  the vector cc into mc 
       uvx=proych*mcx; uvy=proych*mcy; uvz=proych*mcz  !!>> HC 26-3-2021 now we multiply by each coordinate to obtain the direction of the torsion vector
       f=stornod+node(ic)%est                    !!>> HC 26-3-2021          
       a=a+f*ccx ; b=b+f*ccy ; c=c+f*ccz           
           
       dotpich=(cx*icx+cy*icy+cz*icz)                 !!>> HC 26-3-2021 dot product between the unit vectors of nod and ic cells
           
       !!LATERAL TORSION                                          !!>> HC 26-3-2021
       if (abs(dotpich)-angletor>epsilod)then                     !!>> HC 26-3-2021 if the dot product is not 0, cells are not in parallel
          f=tornod+node(ic)%erp                             !!>> HC 26-3-2021  angletor is the minimum angle to feel torsion usually=0  
          if (dotp>epsilod)then                                   !!>> HC 26-3-2021  if this dot product is positive, these nodes are too near one from each other               
             aa=aa+f*ccx ; bb=bb+f*ccy ; cc=cc+f*ccz  !!>> HC 26-3-2021 TORSION COMPONENTS for nod
          else                                                                                        !!>> HC 26-3-2021 if this dot product is negative, they are too far from each other
             aa=aa-f*ccx ; bb=bb-f*ccy ; cc=cc-f*ccz  !!>> HC 26-3-2021 TORSION COMPONENTS for nod                                     
          endif                                                                                       !!>> HC 26-3-2021
       endif                                                                                          !!>> HC 26-3-2021
    endif                                                                                             !!>> HC 26-3-2021
  end do
  px(i)=px(i)+a+aa ; py(i)=py(i)+b+bb ; pz(i)=pz(i)+c+cc
  if (aut==0)  then ; vstorx(i)=a ; vstory(i)=b ; vstorz(i)=c ; end if
  if (aut==0)  then ; vtorx(i)=aa ; vtory(i)=bb ; vtorz(i)=cc ; end if    
end do          

end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!
!                                      TORSION_classic_no_epi_grid
!
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine torsion_classic_no_epi_grid
integer i,j,k,ic,tipi
real*8 a,b,c,cx,cy,cz,ccx,ccy,ccz,f,dotp,tornod,stornod,ix,iy,iz,mcx,mcy,mcz,icx,icy,icz,md,uvx,uvy,uvz,udd,ud,iudd,d
real*8 aa,bb,cc

do i=1,nd
  tipi=node(i)%tipus
  if (tipi>2) cycle
  f=0.0d0
  tornod=node(i)%erp
  stornod=node(i)%est
  cx=ab_pol(i,1) ; cy=ab_pol(i,2) ; cz=ab_pol(i,3)  
  ix=node(i)%x   ; iy=node(i)%y   ; iz=node(i)%z  
  a=0.0d0  ; b=0.0d0  ; c=0.0d0
  aa=0.0d0 ; bb=0.0d0 ; cc=0.0d0
!print *,nneigh(i)  
  do j=1,nneigh(i)
    ic=neigh(i,j)
    if (node(ic)%tipus/=tipi) cycle 
    !surface tension-like torsion (original)
    ccx=node(ic)%x-ix ; ccy=node(ic)%y-iy ; ccz=node(ic)%z-iz
    d=sqrt(ccx**2+ccy**2+ccz**2)
    icx=ab_pol(ic,1) ; icy=ab_pol(ic,2) ; icz=ab_pol(ic,3)
    mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors

    dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
    dotp=dotp*md  ! this goes instead of the previous section
!print *,mcx,mcy,mcz,"m",dotp
    !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
    if(abs(dotp)-angletor*d>epsilod)then
      uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
      f=(stornod+node(ic)%est)*dotp ! NOTICE %EST AIXO ES UNA FORCA CAP AL PROMIG i en la direccio d'aquest >Miquel 7-8-13
!      a=a+f*uvx ; b=b+f*uvy ; c=c+f*uvz
      
!print *,a,b,c,"d"
      !parallel springs torsion   !new  !>>Miquel14-3-14
      dotp=(cx*ccx+cy*ccy+cz*ccz) !vertical projection, more stable than the angle
      f=(tornod+node(ic)%erp)*dotp
      ud=1.0d0/sqrt(ccx**2+ccy**2+ccz**2)
      uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
      aa=aa+f*uvx ; bb=bb+f*uvy ; cc=cc+f*uvz
!print *,f*uvx,f*uvy,f*uvz,"dd"
    end if
  end do
  px(i)=px(i)+a+aa ; py(i)=py(i)+b+bb ; pz(i)=pz(i)+c+cc  
  if (aut==0)  then ; vstorx(i)=a ; vstory(i)=b ; vstorz(i)=c ; end if
  if (aut==0)  then ; vtorx(i)=aa ; vtory(i)=bb ; vtorz(i)=cc ; end if  
  
end do          
          
end subroutine


!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!
!                                      TORSION
!
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine torsion_ap ! apico-basal tension
integer i,j,k,ic,tipi
real*8 a,b,c,cx,cy,cz,ccx,ccy,ccz,f,dotp,tornod,ix,iy,iz

do i=1,nd
  tipi=node(i)%tipus
  if (tipi>2) cycle
  !a=0.0d0 ; b=0.0d0 ; c=0.0d0
  f=0.0d0
  tornod=node(i)%erp
  cx=ab_pol(i,1) ; cy=ab_pol(i,2) ; cz=ab_pol(i,3)  
  ix=node(i)%x   ; iy=node(i)%y   ; iz=node(i)%z  
  do j=1,nneigh(i)
    ic=neigh(i,j)
    if (node(ic)%tipus/=tipi) cycle
    ccx=node(ic)%x-ix ; ccy=node(ic)%y-iy ; ccz=node(ic)%z-iz ! this could be saved, it is re-calculated          
    dotp=cx*ccx+cy*ccy+cz*ccz !vertical projection, cx is already normalized, ccx does not need to
    f=f+(tornod+node(ic)%erp)*dotp
    ! ATENCIO, tornod pot anar a fora del loop per anar mes rapid
!    if (i==1.or.i==8) print *,i,ic,f,dotp,cx,cy,cz
  end do
  if (aut==0)  then ; vstorx(i)=cx ; vstory(i)=cy ; vstorz(i)=cz ; end if
  if (aut==0)  then ; vtorx(i)=f*cx ; vtory(i)=f*cy ; vtorz(i)=f*cz ; end if  
  px(i)=px(i)+f*cx ; py(i)=py(i)+f*cy ; pz(i)=pz(i)+f*cz    
end do
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   REP_ADH

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine rep_adh
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe
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

      ic=neigh(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      if (tipi<3)then
        if(tipic<3)then
          ivv=node(ic)%altre
          icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 

          mcx=icx+cx; mcy=icy+cy; mcz=icz+cz ! this is the sum of the cylinder vectors of i and ic
          md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  ! we do that so that we can calculate what happens to i and ic and the same time
          dotp=mcx*ccx+mcy*ccy+mcz*ccz       ! dot product between mc and cc (the vector from i to ic)
          ddd=dotp*md  !vertical component   ! this is the projection of cc into mc  see Fig 2
                   
          posca=icx*cx+icy*cy+icz*cz      ! Is 3-3-21 THIS is FIG.1!!!!! 
          
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              !!! Is 3-3-21 FIG.2  
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-nodda-node(ic)%add>epsilod) cycle !>>HC ISAAC 14-5-2020 !>> HC 12-6-2020
              pesco=ddd*md   ! Is 3-3-21 IT WOULD BETTER BE pesco=ddd*md, it would be the same by faster
              a=1.0d0/ad              
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins(nod)=epinveins(nod)+1
              epinveins(ic)=epinveins(ic)+1
            else
              ! SO THIS IS APICAL-APICAL OR BASAL-BASAL INTERACTION but by cells that are not close in the plane of the epi
                                      ! FIG 3 IS 3-3-21 FIG 3
              ddd=abs(ddd)  !vertical component              
              a=nodda+node(ic)%add
              if(ddd-a<epsilod) cycle  ! ic is too far away in the vertical component
              ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
              ad=sqrt(ad)              !lateral component                              
              if (ad-a>epsilod) cycle  ! ic is too far away in the lateral component
              pesco=dotp*md*md  ! REVISAR
              a=1/ddd                    
              uvx=(ccx-mcx*pesco)*a ; uvy=(ccy-mcy*pesco)*a ; uvz=(ccz-mcz*pesco)*a  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              fd=ddd     !distance used, vertical component
              twoep=2
            end if
          else
            if (posca<epsilod) then           
              ! here there is an apical interacting with a basal or vice versa
              ! this is adhesion and repulsion between epithelial sides
              if (dotp<epsilod) then         
                ddd=-dotp*udd                ! this is projection of vector cc into the cylinder vector of i (the direction UP face)
                a=nodda+node(ic)%add         ! this is still figure 3
                if (ddd-a<epsilod) then
                  ad=d**2-ddd**2 ;if(ad<epsilod) ad=epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                else
                  cycle
                end if
              end if
            end if
            cycle  ! Is 3-3-21 This never happens??? since it would be detected before with posca???
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=cx*ccx+cy*ccy+cz*ccz!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=nodda+node(ic)%add
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
            a=nodda+node(ic)%add
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
          if (fd-nodda-node(ic)%add>epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
        end if
      end if

300   nuve=nuve+1

      !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, fd, NOW WE CALCULATE THE ACTUAL ENERGIES
      !if (ffu(2)==1) then IS 25-12-13
        if(node(nod)%icel==node(ic)%icel)then    !>> HC 12-6-2020 forces back to normal
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
            f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
          else                                   !>> HC 12-6-2020
            if(ffu(28)==1 .and. node(nod)%tipus .le. 2 .and. node(ic)%tipus .le. 2) then !>> HC 27-7-2020
              ! We use the grid and both nodes are epithelial.
              if(any(ic .eq. epi_grid(nod, :epi_ngrid(nod)))) then !>> HC 27-7-2020
                f=(younod+node(ic)%you)*(fd-deqe)                                   !>> HC 27-7-2020
              else                                                 !>> HC 27-7-2020
                f=0.0d0                                            !>> HC 27-7-2020
              endif                                                !>> HC 27-7-2020
            else                                                   !>> HC 27-7-2020
              f=(younod+node(ic)%you)*(fd-deqe) 
            endif                                                  !>> HC 27-7-2020
          end if                                 !>> HC 12-6-2020
        else                                     !>> HC 12-6-2020
          deqe=(reqnod+node(ic)%eqd)             !>>>> Miquel 16-8-13
          if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
            f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
          else                                   !>> HC 12-6-2020
            if(ffu(28)==1 .and. node(nod)%tipus .le. 2 .and. node(ic)%tipus .le. 2) then !>> TT 24-7-2020
              ! We use the grid and both nodes are epithelial.
              if(any(ic .eq. epi_grid(nod, :epi_ngrid(nod)))) then !>> TT 24-7-2020
                ! Nodes can adhere to each other only if they are connected in the grid
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
              else                                   !>> HC 27-7-2020
                f=0.0d0                              !>> HC 27-7-2020
              end if
            else !>> TT 24-7-2020 ffu(28)==0 or both nodes are not epithelial
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
            end if
          end if                                 !>> HC 12-6-2020
        end if !;print*,"f",f,"fd",fd        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
       
        a=f*uvx ; b=f*uvy ; c=f*uvz
        if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
          if (f>0.0d0)then                     
            arcilxch(nod)=arcilxch(nod)+a ; arcilych(nod)=arcilych(nod)+b ; arcilzch(nod)=arcilzch(nod)+c 
            arcilxch(ic)=arcilxch(ic)-a   ; arcilych(ic)=arcilych(ic)-b   ; arcilzch(ic)=arcilzch(ic)-c   !!>> HC 16-6-2020
          else                                
            rrcilxch(nod)=rrcilxch(nod)+a ; rrcilych(nod)=rrcilych(nod)+b ; rrcilzch(nod)=rrcilzch(nod)+c !!>> HC 15-6-2020        
            rrcilxch(ic)=rrcilxch(ic)-a   ; rrcilych(ic)=rrcilych(ic)-b   ; rrcilzch(ic)=rrcilzch(ic)-c   !!>> HC 16-6-2020
          endif                                
        else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
          rcilx(nod)=rcilx(nod)+a ; rcily(nod)=rcily(nod)+b ; rcilz(nod)=rcilz(nod)+c  !>>>Miquel 7-8-13
          rcilx(ic)=rcilx(ic)-a   ; rcily(ic)=rcily(ic)-b   ; rcilz(ic)=rcilz(ic)-c  !>>>Miquel 4-4-14
        endif                                 
                   
        !if (ffu(15)==3) then  ! 4-3-2020
        !  fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        !else
        !  fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        !end if
    end do

    if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
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

    if (ffu(8)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if
    
    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
!      vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
!      vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    end if

!    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
!      a=epinveins(nod) ! >>> Is 7-6-14
!      a=1.0d0/a        ! >>> Is 7-6-14
!      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
!      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
!      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
!    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
!    end if

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

    if (aut==0) dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  return
171 call force_errors    
end subroutine

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   REP_ADH_EPI_GRID

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine rep_adh_epi_grid
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda
real*8   :: lat_dist,up_dist
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

    do i=1,epi_ngrid(nod) !nneigh(nod)

!      ic=neigh(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact
      ic=epi_grid(nod,i)
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      d=sqrt(ccx**2+ccy**2+ccz**2)
      !d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      ivv=node(ic)%altre
      icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 

      mcx=icx+cx; mcy=icy+cy; mcz=icz+cz ! this is the sum of the cylinder vectors of i and ic
      md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  ! we do that so that we can calculate what happens to i and ic and the same time
      dotp=mcx*ccx+mcy*ccy+mcz*ccz       ! dot product between mc and cc (the vector from i to ic)
      up_dist=dotp*md  !vertical component   ! this is the projection of cc into mc  see Fig 2
      ad=d**2-up_dist**2 ; if(ad<epsilod) ad=epsilod !ATENCIO: aixo ultim potser no cal
      lat_dist=sqrt(ad)              !lateral component

      !WE MAKE THAT cells connected in the grid adhere to each other laterally independently of their distance                             
      if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
        !!! Is 3-3-21 FIG.2  
!        if (d<2*(node(i)%add+nodda)) then      !the 2 is arbotrary it should be a parameter
          up_dist=up_dist*md                 !this is the up component
          a=1.0d0/lat_dist              
          uvx=(ccx-mcx*up_dist)*a ; uvy=(ccy-mcy*up_dist)*a ; uvz=(ccz-mcz*up_dist)*a  
                                                               ! mcx*up_dist is the projection vector of cc into c
                                                               ! and ccx-cx*c is the rest of these two vectors that is in
                                                               ! in fact the lateral component in the cylinder
                                                               ! then we multiply by a to make it unit vector

           fd=lat_dist  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
!         else
!           ! if the distance is too large we make the force central
!print *,nod,ic,d,"d"           
!           uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the FORCE          
!           fd=d
!         end if
      else
        ! WARNING: DOES NOT WORK WELL WITH CROSS-LINKS
        ! SO THIS IS APICAL-APICAL OR BASAL-BASAL INTERACTION but by cells that are not close in the plane of the epi
        ! FIG 3 IS 3-3-21 FIG 3
        ad=d**2-up_dist**2 ; if(ad<epsilod) ad=epsilod
        lat_dist=sqrt(ad)              !lateral component                              
        up_dist=up_dist*md
        a=1/up_dist                    
        uvx=(ccx-mcx*up_dist)*a ; uvy=(ccy-mcy*up_dist)*a ; uvz=(ccz-mcz*up_dist)*a  
        fd=up_dist     !distance used, vertical component
      end if

      deqe=(reqnod+node(ic)%eqd)      
      if(node(nod)%icel==node(ic)%icel)then    ! ATENCIO: ES POT OPTIMITZAR PER AL CAS QUE JA SABEMQ QUE AIXO NO ES DONA
        if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
          f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
        else                                   !>> HC 12-6-2020
          f=(younod+node(ic)%you)*(fd-deqe)                                   !>> HC 27-7-2020
        end if
      else                                     !>> HC 12-6-2020
        if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
          f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
        else                                   !>> HC 12-6-2020
          adhe=0.5d0*(adhnod+node(ic)%adh)     !>>>> Miquel 16-8-13
          if (npag(1)>0) then                !ATENCIO, OPTIMITZABLE PQ NO VARIA DURANT UNA SIMULACIO
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
      end if
       
      a=f*uvx ; b=f*uvy ; c=f*uvz
        if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
          if (f>0.0d0)then                     
            arcilxch(nod)=arcilxch(nod)+a ; arcilych(nod)=arcilych(nod)+b ; arcilzch(nod)=arcilzch(nod)+c 
            arcilxch(ic)=arcilxch(ic)-a   ; arcilych(ic)=arcilych(ic)-b   ; arcilzch(ic)=arcilzch(ic)-c   !!>> HC 16-6-2020
          else                                
            rrcilxch(nod)=rrcilxch(nod)+a ; rrcilych(nod)=rrcilych(nod)+b ; rrcilzch(nod)=rrcilzch(nod)+c !!>> HC 15-6-2020        
            rrcilxch(ic)=rrcilxch(ic)-a   ; rrcilych(ic)=rrcilych(ic)-b   ; rrcilzch(ic)=rrcilzch(ic)-c   !!>> HC 16-6-2020
          endif                                
        else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
          rcilx(nod)=rcilx(nod)+a ; rcily(nod)=rcily(nod)+b ; rcilz(nod)=rcilz(nod)+c  !>>>Miquel 7-8-13
          rcilx(ic)=rcilx(ic)-a   ; rcily(ic)=rcily(ic)-b   ; rcilz(ic)=rcilz(ic)-c  !>>>Miquel 4-4-14
        endif                                 
                   
        !if (ffu(15)==3) then  ! 4-3-2020
        !  fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        !else
        !  fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        !end if
    end do

    if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
      ftch=0
      ftch=sqrt(arcilxch(nod)**2+arcilych(nod)**2+arcilzch(nod)**2) !!>> HC 15-6-2020 modul adhesion
      b=maxad/ftch
      if (ftch>maxad) then !>>HC 14-05-2020 When the adhesion vector is larger than the limit(maxad)
	arcilxch(nod)=arcilxch(nod)*b ; arcilych(nod)=arcilych(nod)*b ; arcilzch(nod)=arcilzch(nod)*b !!>> HC 15-6-2020
      endif 
      rcilx(nod)=arcilxch(nod)+rrcilxch(nod) ; rcily(nod)=arcilych(nod)+rrcilych(nod) ;rcilz(nod)=arcilzch(nod)+rrcilzch(nod) 
    endif
    if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
      !fmeanl(nod)=0
    !else
      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    end if

    if (ffu(8)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if
    
    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
    end if

!    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
!      a=epinveins(nod) ! >>> Is 7-6-14
!      a=1.0d0/a        ! >>> Is 7-6-14
!      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
!      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
!      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
!    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
!    end if

    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%fix==1)then
      uvx=node(nod)%orix-ix  ; uvy=node(nod)%oriy-iy ; uvz=node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%kfi
      rvx=rvx+uvx*a*ud ; rvy=rvy+uvy*a*ud  ; rvz=rvz+uvz*a*ud
    end if

    if (aut==0) dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  return
171 call force_errors    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

!                   REP_ADH_NO_EPI_GRID

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************

subroutine rep_adh_no_epi_grid
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda
real*8   :: lat_dist,up_dist
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
    if (node(nod)%tipus>2) cycle
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

      ic=neigh(nod,i)  ! ic THIS IS THE NODE WITH WHICH i may interact
!      ic=epi_grid(nod,i)
if (node(ic)%tipus/=tipi) cycle ! WARNING, THIS IS TEMPORARY
      if(ic<nod.and.node(ic)%fix/=2)then !we have calculated that interaction already !>>Miquel6-5-15
        alone=1
        cycle
      end if

      !so it turns out that the nod-ic interactions has not been calculated before
      bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=dneigh(nod,i)
      ud=1d0/d
      tipic=node(ic)%tipus
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares

      ivv=node(ic)%altre
      icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz 

      mcx=icx+cx; mcy=icy+cy; mcz=icz+cz ! this is the sum of the cylinder vectors of i and ic
      md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  ! we do that so that we can calculate what happens to i and ic and the same time
      dotp=mcx*ccx+mcy*ccy+mcz*ccz       ! dot product between mc and cc (the vector from i to ic)
      up_dist=dotp*md  !vertical component   ! this is the projection of cc into mc  see Fig 2
      ad=d**2-up_dist**2 ; if(ad<epsilod) ad=epsilod !ATENCIO: aixo ultim potser no cal
      lat_dist=sqrt(ad)              !lateral component

      !WE MAKE THAT cells connected in the grid adhere to each other laterally independently of their distance                             
      if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
        !!! Is 3-3-21 FIG.2  
!        if (d<2*(node(i)%add+nodda)) then      !the 2 is arbotrary it should be a parameter
          up_dist=up_dist*md                 !this is the up component
          a=1.0d0/lat_dist              
          uvx=(ccx-mcx*up_dist)*a ; uvy=(ccy-mcy*up_dist)*a ; uvz=(ccz-mcz*up_dist)*a  
                                                               ! mcx*up_dist is the projection vector of cc into c
                                                               ! and ccx-cx*c is the rest of these two vectors that is in
                                                               ! in fact the lateral component in the cylinder
                                                               ! then we multiply by a to make it unit vector

           fd=lat_dist  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
!         else
!           ! if the distance is too large we make the force central
!print *,nod,ic,lat_dist,"d"           
!           uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the FORCE          
!           fd=d
!         end if
      else
        ! WARNING: DOES NOT WORK WELL WITH CROSS-LINKS
        ! SO THIS IS APICAL-APICAL OR BASAL-BASAL INTERACTION but by cells that are not close in the plane of the epi
        ! FIG 3 IS 3-3-21 FIG 3
        ad=d**2-up_dist**2 ; if(ad<epsilod) ad=epsilod
        lat_dist=sqrt(ad)              !lateral component                              
        up_dist=up_dist*md
        a=1/up_dist                    
        uvx=(ccx-mcx*up_dist)*a ; uvy=(ccy-mcy*up_dist)*a ; uvz=(ccz-mcz*up_dist)*a  
        fd=up_dist     !distance used, vertical component
      end if

      deqe=(reqnod+node(ic)%eqd)      
      if(node(nod)%icel==node(ic)%icel)then    ! ATENCIO: ES POT OPTIMITZAR PER AL CAS QUE JA SABEMQ QUE AIXO NO ES DONA
        if(fd-deqe<-epsilod)then 	         !>> HC 12-6-2020 			
          f=(fd-deqe)*(repnod+node(ic)%rep)     !>>>> Miquel 16-8-13
        else                                   !>> HC 12-6-2020
          f=(younod+node(ic)%you)*(fd-deqe)                                   !>> HC 27-7-2020
        end if
      else                                     !>> HC 12-6-2020
        if(fd-deqe<-epsilod)then               !>> HC 12-6-2020
          f=(repcelnod+node(ic)%rec)*(fd-deqe)                  !>> HC 12-6-2020 !in fact that is the derivative of the energy
        else                                   !>> HC 12-6-2020
          adhe=0.5d0*(adhnod+node(ic)%adh)     !>>>> Miquel 16-8-13
          if (npag(1)>0) then                !ATENCIO, OPTIMITZABLE PQ NO VARIA DURANT UNA SIMULACIO
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
      end if
       
      a=f*uvx ; b=f*uvy ; c=f*uvz
        if (ffu(27)==1)then          !!>> HC 15-6-2020 If there is a limit to adhesion, we have to add
          if (f>0.0d0)then                     
            arcilxch(nod)=arcilxch(nod)+a ; arcilych(nod)=arcilych(nod)+b ; arcilzch(nod)=arcilzch(nod)+c 
            arcilxch(ic)=arcilxch(ic)-a   ; arcilych(ic)=arcilych(ic)-b   ; arcilzch(ic)=arcilzch(ic)-c   !!>> HC 16-6-2020
          else                                
            rrcilxch(nod)=rrcilxch(nod)+a ; rrcilych(nod)=rrcilych(nod)+b ; rrcilzch(nod)=rrcilzch(nod)+c !!>> HC 15-6-2020        
            rrcilxch(ic)=rrcilxch(ic)-a   ; rrcilych(ic)=rrcilych(ic)-b   ; rrcilzch(ic)=rrcilzch(ic)-c   !!>> HC 16-6-2020
          endif                                
        else                                    !!>> HC 15-6-2020 If there is no limin, we do it conventionally
          rcilx(nod)=rcilx(nod)+a ; rcily(nod)=rcily(nod)+b ; rcilz(nod)=rcilz(nod)+c  !>>>Miquel 7-8-13
          rcilx(ic)=rcilx(ic)-a   ; rcily(ic)=rcily(ic)-b   ; rcilz(ic)=rcilz(ic)-c  !>>>Miquel 4-4-14
        endif                                 
                   
        !if (ffu(15)==3) then  ! 4-3-2020
        !  fmeanl(nod)=fmeanl(nod)+(fd-deqe) ; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 
        !else
        !  fmeanl(nod)=fmeanl(nod)+f ; fmeanl(ic)=fmeanl(ic)+f ! >>> ! 4-3-2020 aixo es una mica mes correcte       
        !end if
    end do

    if (ffu(27)==1) then !>>HC 14-05-2020 This creates a limit to the adhesion force suffered by the node
      ftch=0
      ftch=sqrt(arcilxch(nod)**2+arcilych(nod)**2+arcilzch(nod)**2) !!>> HC 15-6-2020 modul adhesion
      b=maxad/ftch
      if (ftch>maxad) then !>>HC 14-05-2020 When the adhesion vector is larger than the limit(maxad)
	arcilxch(nod)=arcilxch(nod)*b ; arcilych(nod)=arcilych(nod)*b ; arcilzch(nod)=arcilzch(nod)*b !!>> HC 15-6-2020
      endif 
      rcilx(nod)=arcilxch(nod)+rrcilxch(nod) ; rcily(nod)=arcilych(nod)+rrcilych(nod) ;rcilz(nod)=arcilzch(nod)+rrcilzch(nod) 
    endif
    if(epinveins(nod)>0)then                  !>>>Miquel24-3-14
      !fmeanl(nod)=0
    !else
      fmeanl(nod)=fmeanl(nod)/epinveins(nod)
    end if

    if (ffu(8)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        node(nod)%talone=node(nod)%talone+1
      else
        node(nod)%talone=0
      end if
    end if
    
    if (aut==0) then ! if aut==0 there is not display   
      vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
      vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
    end if

!    if (epinveins(nod)>0) then  ! >>> Is 7-6-14
!      a=epinveins(nod) ! >>> Is 7-6-14
!      a=1.0d0/a        ! >>> Is 7-6-14
!      rvx=rsprx+rcilx(nod)+(rtorx(nod)+rstorx(nod))*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
!      rvy=rspry+rcily(nod)+(rtory(nod)+rstory(nod))*a  ! >>> Is 7-6-14
!      rvz=rsprz+rcilz(nod)+(rtorz(nod)+rstorz(nod))*a  ! >>> Is 7-6-14
!    else
      rvx=rsprx+rcilx(nod)  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily(nod)  ! >>> Is 7-6-14
      rvz=rsprz+rcilz(nod)  ! >>> Is 7-6-14
!    end if

    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(node(nod)%fix==1)then
      uvx=node(nod)%orix-ix  ; uvy=node(nod)%oriy-iy ; uvz=node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>epsilod)then ; ud=1d0/sqrt(d);else;d=0;end if
      a=node(nod)%kfi
      rvx=rvx+uvx*a*ud ; rvy=rvy+uvy*a*ud  ; rvz=rvz+uvz*a*ud
    end if

    if (aut==0) dex(nod)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    px(nod)=rvx ; py(nod)=rvy ; pz(nod)=rvz
  end do
  return
171 call force_errors    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine ordenarepeq(ma,mt,rang)
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
end subroutine ordenarepeq


end module biomechanic
