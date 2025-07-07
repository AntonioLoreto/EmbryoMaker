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




module energy
use general
use genetic

contains


!**************************************************************************

subroutine energia(nod) !>REMADE SUBSTANTIALLY Isaac 5-6-13
real*8   ::ix,iy,iz,dd
real*8   ::a,b,c,d,e,f,g
real*8   ::ie,posca
integer  ::nod
real*8   ::ene,ened,youe,repe,adhe,repcele,deqe,ideqe,adhed
real*8   ::reqnod,younod,repnod,adhnod,repcelnod,tornod,stornod
real*8   ::ax,ay,az,bx,by,bz
real*8   ::udd
real*8   ::cx,cy,cz,ccx,ccy,ccz,pesc
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md
real*8   ::nodda
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv
integer  ::kj,kjjj
integer  ::loi,loj,lok
integer  ::nuve,twoep
real*8   ::ener(nd),enea(nd),enet(nd),iener(nd),ienea(nd)		!we will store the energies calculated for each node, then calculate the shielding and apply to them
real*8   ::enes,inuve
integer  ::tipi,tipic	!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
real*8,  allocatable   ::upr(:)
integer, allocatable  ::iupr(:)
real*8   ::r(nd),er(nd)
real*8   ::dotp,ddd
real*8   ::ad,fd  !>>Miquel28-1-14

  if (node(nod)%fix==2) then  ! >>> Is 30-6-14
    node(nod)%e=huge(a)        ! >>> Is 30-6-14
    return                     ! >>> Is 30-6-14
  end if                       ! >>> Is 30-6-14

  nuve=0
  nodda=node(nod)%add
  ix=node(nod)%x      ; iy=node(nod)%y     ; iz=node(nod)%z
  tipi=node(nod)%tipus											!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
  !iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)
  ie=0.0
  ener=0d0 ; iener=0d0 ; ienea=0d0 ; enea=0d0 ; enet=0d0 ; enes=0d0 !energy storing matrices (es guarda l'energia aqui per aplicar després l'apantallament		!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
  younod=node(nod)%you																			
  repnod=node(nod)%rep																			
  adhnod=node(nod)%adh   !default inespecific adhesion of node nodmo																			
  repcelnod=node(nod)%rec
  tornod=node(nod)%erp                                            !
  stornod=node(nod)%est                                          !
  reqnod=node(nod)%eqd

  if(tipi<3)then	!Springs is only for epithelia										!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
   !SPRINGS calculem primer l'energia per a les molles
    iv=node(nod)%altre
    ax=node(iv)%x   ; ay=node(iv)%y    ; az=node(iv)%z
    cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
    dd=sqrt(cx**2+cy**2+cz**2)
    udd=1d0/dd
    !ENERGY COMING FROM SPRINGS BETWEEN EVERY LOWER AND UPPER NODE
    enes=node(nod)%hoo*abs(dd-node(nod)%eqs)**2
    ie=enes
  end if

  !NODE's REPULSIONS AND ADHESIONS
  !do ii=-1,1	!codi per explorar els cubs adjacents
    !do jj=-1,1
      !do kk=-1,1
      !print*,nod,"icel",node(nod)%icel,"neigh",neigh(nod,1:nneigh(nod))
  do ii=1,nneigh(nod)
    ic=neigh(nod,ii)
    !do while(ic.ne.0)
    !if(ic==nod)then;ic=list(ic);cycle;end if
    !if(ic==iv)then ;ic=list(ic);cycle;end if  !el vei el passem de llarg
    bx=node(ic)%x   ; by=node(ic)%y    ; bz=node(ic)%z
    ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
    ddd=0d0  !es necessari?
    d=sqrt(ccx**2+ccy**2+ccz**2)

    tipic=node(ic)%tipus		    
    twoep=0
    if(tipi<3)then
      if(tipic<3)then
        ivv=node(ic)%altre
        icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz
        idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
        posca=icx*cx+icy*cy+icz*cz
        if (tipi==tipic) then        ! we are equal so we must have lateral adhesion
          !we have now to check that they are contiguous cells, that we are not in two folds touching each other 
          if (posca>epsilod) then
            idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
            mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
            ddd=abs(mcx*ccx+mcy*ccy+mcz*ccz)*md
            ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
            ad=sqrt(ad)
            if (ad-nodda-node(ic)%add>epsilod) cycle
            fd=ad
            twoep=1
          else
            ! apical/basal contact, two epithelia from the same side
            mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
            if (md<epsilod)then !the two cylinders are parallel, the vector used is the spring !>>Miquel20-8-14
               dotp=(cx*ccx+cy*ccy+cz*ccz)
               if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                  ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                  a=nodda+node(ic)%add
                  if (ddd-a<epsilod) then
                     ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                     ddd=sqrt(ddd)              !lateral component
                     if (ad-a>epsilod) cycle
                     fd=ddd     !distance used, vertical component
                     twoep=2
                     goto 301
                  end if
                  cycle
               end if
            else  !proper apical cylindric interface applied !>>Miquel20-8-14
               md=1d0/md
               dotp=(mcx*ccx+mcy*ccy+mcz*ccz)
               ad=abs(dotp)*md  !vertical component
               a=nodda+node(ic)%add
               if (ad-a<epsilod)then
                  ddd=d**2-ad**2 ;if(ddd<epsilod) ddd=epsilod
                  ddd=sqrt(ddd)              !lateral component
                  if (ddd-a>epsilod) cycle
                  fd=ddd     !distance used, vertical component
                  twoep=2
                  goto 301
               end if
               cycle
            end if
          end if
        else
          if (posca<0.0d0) then           
            dotp=(cx*ccx+cy*ccy+cz*ccz)                     
            if (dotp<epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
              ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
              a=nodda+node(ic)%add
              if (ddd-a<epsilod) then
                ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
                ad=sqrt(ad)
                if (ad-a>epsilod) cycle
                fd=ddd
                twoep=2
                goto 301
              end if
            end if
          end if
          cycle
        end if
      else
        ! nod IS EPITHELIAL and ic is not
        ! we check the distance to radial distance in the plane of the ic cylinder 
        dotp=(cx*ccx+cy*ccy+cz*ccz) !hi ha un vector del revés, per tant això està al revés també
        if (dotp<0.0) then
          ddd=abs(dotp)*udd
          a=nodda+node(ic)%add
          if (ddd-a<epsilod) then
            ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
            ad=sqrt(ad)
            if (ad-a>epsilod) cycle
            fd=ddd
            !fd=d
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
        idd=1/sqrt(icx**2+icy**2+icz**2)
        dotp=(icx*ccx+icy*ccy+icz*ccz)              ! projection of the vector from nod to ic into the vector from ic to ivv
        if (dotp>0.0) then  ! >>> Is 29-6-14 there was a bug here, a serious typo
          ddd=dotp*idd
          a=nodda+node(ic)%add
          if (ddd-a<0.0) then
            ad=d**2-ddd**2 ; if(ad<epsilod) ad=epsilod
            ad=sqrt(ad)
            if (ad-a>epsilod) cycle
            fd=ddd
            !fd=d
          else
            cycle
          end if		
        else
          cycle
        end if
      else
        fd=d !BOTH NODES ARE NON-EPITHELIAL: we just consider the interactions between nodes as such
        if (fd-nodda-node(ic)%add>epsilod) cycle
      end if
    end if

301 nuve=nuve+1		!apantallaments
    r(nuve)=fd              

    !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, ddd, NOW WE CALCULATE THE ACTUAL ENERGIES
    if(node(nod)%icel==node(ic)%icel)then
      youe=0.5*(younod+node(ic)%you)
      repe=0.5*(repnod+node(ic)%rep)
      deqe=reqnod+node(ic)%eqd
      ideqe=((nodda+node(ic)%add-deqe)/deqe)**2
      if(fd-deqe<-epsilod)then 				
        ener(nuve)=ener(nuve)+repe*((fd-deqe)/deqe)**2-youe*ideqe !((node(nod)%add-deqe)/deqe)**2 !this is the repulsion energy for entering the cilinder
      else
        enea(nuve)=enea(nuve)+youe*((fd-deqe)/deqe)**2-youe*ideqe !((node(nod)%add-deqe)/deqe)**2 !this is the adhesion or you 
      end if
    else																							
      adhe=0.5d0*(adhnod+node(ic)%adh) !adhesion is necessarily symmetric
      repcele=0.5*(repcelnod+node(ic)%rec)
      deqe=reqnod+node(ic)%eqd
      ideqe=((nodda+node(ic)%add-deqe)/deqe)**2
      if (npag(1)>0) then ! we have adhesion molecules
        do j=1,npag(1)
          kj=whonpag(1,j)
          do kjjj=1,npag(1)
            kkkk=whonpag(1,kjjj)
            if (gex(nod,kj)>0.0d0.and.gex(ic,kkkk)>0.0d0) then     
     adhe=adhe+gex(nod,kj)*gex(ic,kkkk)*kadh(int(gen(kj)%e(1)),int(gen(kkkk)%e(1)))    !this is specific adhesion
            end if
          end do
        end do
      end if

      if(fd-deqe<-epsilod)then 
        iener(nuve)=iener(nuve)+repcele*((fd-deqe)/deqe)**2-adhe*ideqe !((node(nod)%add-deqe)/deqe)**2!this is the repulsion energy for entering the cilinder
                                                            !this is because energy needs to be smaller(negative) 
                                                            !this is from f(x)=ax**2-c we need to make it that 
                                                            !the function crosses the x axis at node(i)%add
                                                            !where the x axis is the distance between two nodes
                                                            !then c=a*da**2
      else
        ienea(nuve)=ienea(nuve)+adhe*((fd-deqe)/deqe)**2-adhe*ideqe !((node(nod)%add-deqe)/deqe)**2 !this is the adhesion or you 
                                                            !this is because energy needs to be smaller(negative) 
                                                            !this is from f(x)=ax**2-c we need to make it that 
                                                            !the function crosses the x axis at node(i)%add
                                                            !where the x axis is the distance between two nodes
                                                            !then c=a*da**2
      end if																							
    end if								
    er(nuve)=ener(nuve)+enea(nuve)+iener(nuve)+ienea(nuve)
  
    !TORSION																								
    if(ffu(3)==0 .and. twoep==1) then !it is only between epithelial nodes
      !ivv=node(ic)%altre
      !icx=node(ivv)%x-bx ; icy=node(ivv)%y-by ; icz=node(ivv)%z-bz; 
      !idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	


      mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors		!>>>>Miquel 30-4-13
      dotp=((mcx*ccx+mcy*ccy+mcz*ccz)/md)**2											!>>>>Miquel 30-4-13
      if(abs(dotp)-angletor*d>epsilod)then

        enet(nuve)=enet(nuve)+dotp*(stornod+node(ic)%est)


        dotp=((cx*ccx+cy*ccy+cz*ccz)*udd)**2 !vertical projection, more stable than the angle
        enet(nuve)=enet(nuve)+dotp*(tornod+node(ic)%erp) !tor torsion
      end if

      !dotp=((cx*ccx+cy*ccy+cz*ccz)*udd)**2 !vertical projection, more stable than the angle
      !enet(nuve)=enet(nuve)+(1-abs(cx*icx+cy*icy+cz*icz)*udd*iudd)*node(nod)%erp !making the springs in paralel
     
      !surface tension-like torsion (original)
      !mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors		!>>>>Miquel 30-4-13
      !dotp=((mcx*ccx+mcy*ccy+mcz*ccz)/md)**2											!>>>>Miquel 30-4-13
      !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
    end if
  end do
        !ic=list(ic)
        !end do
      !end do
    !end do
  !end do


    do i=1,nuve
      ie=ie+ener(i)+iener(i)+enea(i)+ienea(i)+enet(i)		!we may want to apply shielding only to some type of forces
    end do
  !end if
  node(nod)%e=ie
  ! end >>>>>>>>< Is 22-6-13

  erep(nod)=sum(ener)
  erepcel(nod)=sum(iener)
  eyou(nod)=sum(enea)
  eadh(nod)=sum(ienea)
  etor(nod)=sum(enet)
  espring(nod)=enes

end subroutine

end module
