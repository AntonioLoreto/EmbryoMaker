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




module ecm !by Miquel 4-6-13
use general
use genetic
use neighboring
use io
use aleas

contains

!********************************************************************************************

subroutine should_I_secrete
integer::i,j,k,ii,jj,kk,iii,jjj,kkk,l,iy
real*8::grac,c

  do iy=1,nd
    if (node(iy)%tipus==4) cycle
    if(node(iy)%fix==1) cycle !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
    grac=0
    c=1-node(iy)%dif
    !if(node(iy)%tipus<4)then
      do jj=1,npag(nparam_per_node+4)  !number of genes affecting secretion
        k=whonpag(nparam_per_node+4,jj)  !which are those genes
        grac=grac+gex(iy,k)*gen(k)%e(nparam_per_node+4)*c  ! wa in space-req units
      end do
      grac=grac*delta
      node(iy)%ecm=node(iy)%ecm+grac
      if(node(iy)%ecm>=ecmmax)then
        call ecm_secretion(iy)
        node(iy)%ecm=node(iy)%ecm-ecmmax
        if (node(iy)%ecm<0.0) node(iy)%ecm=0.0
      end if
    !end if
  end do
end subroutine

!********************************************************************************************

subroutine ecm_secretion(nods)
implicit none
integer::celd,tipi,i,j,k,ii,jj,kk,nnod,nods
real*8::ax,ay,az,cx,cy,cz,dd,d,pox,poy,poz,pesc,uvx,uvy,uvz,a,b,c,dix,diy,diz,bx,by,bz
integer,allocatable::clist(:)
type(nod),allocatable :: cnode(:),cpnode(:,:),cnodeo(:)
real*8,allocatable::cpx(:),cpy(:),cpz(:),cdex(:),cgex(:,:)               !>>>>>Miquel 17-6-13
real*8,allocatable::cvcil(:),cvtor(:),cvstor(:),cvspr(:)              !>>>>>Miquel 20-6-13
integer,allocatable::cneigh(:,:),cdneigh(:,:),cnneigh(:)
real*8,allocatable::cfmeanv(:),cfmeanl(:),cdidpol(:,:)
real*8 :: dch,rax,ray,raz,vxch,vych,vzch,axch,aych,azch,alfach  !!>>HC 9-10-2020
real*8 :: polxch,polych,polzch,cxch,cych,czch                   !!>>HC 9-10-2020
integer :: kch,jch,jjch,kkch,ic,icc                             !!>>HC 9-10-2020

  nd=nd+1 ; ndx=ndx+1

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
    
!    deallocate(agex) ; allocate(agex(nda,ng)) ; agex=0
    
    allocate(cvcil(nd),cvtor(nd),cvstor(nd),cvspr(nd))  !force components storing arrays !>>>>>>>>>>>>>> Miquel 20-6-13

    allocate(cneigh(nda,omnn))
    cneigh(1:nd,:)=neigh(1:nd,:)
    deallocate(neigh)
    allocate(neigh(nda,omnn))
    neigh(:nd,:)=cneigh(:nd,:)
    deallocate(cneigh)
    allocate(cdneigh(nd,omnn))
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
    
    !if(single==1 .and. npag(nparam_per_node+8)>0)then
    !  allocate(cdidpol(nd,npag(nparam_per_node+8)))
    !  cdidpol(1:nd,:)=didpol(1:nd,:)
    !  deallocate(didpol) ; allocate(didpol(nda,npag(nparam_per_node+8)))
    !  didpol(1:nd,:)=cdidpol(1:nd,:)
    !  deallocate(cdidpol)
    !end if

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

  celd=node(nods)%icel
  tipi=node(nods)%tipus

  if(tipi<3)then !when the cell is epithelial we put the ECM over the apical/basal surface, in the direction of the apical-basal polarity

    ii=node(nods)%altre
    ax=node(nods)%x ; ay=node(nods)%y ; az=node(nods)%z
    dd=1d0/sqrt((ax-node(ii)%x)**2+(ay-node(ii)%y)**2+(az-node(ii)%z)**2) ;
    cx=(node(nods)%x-node(ii)%x)*dd ; cy=(node(nods)%y-node(ii)%y)*dd ; cz=(node(nods)%z-node(ii)%z)*dd !apical-basal vector
    call random_number(a)
    d=node(nods)%eqd+a*0.1d0 !>>Miquel15-8-14
    node(nd)%x=ax+cx*d ; node(nd)%y=ay+cy*d ; node(nd)%z=az+cz*d

  else !when the cell is mesenchymal

    !the ECM node is placed on the line intersecting the centroid of the cell and the secreting node, and close to that node
    ii=node(nods)%icel

    
    cx=cels(celd)%cex ; cy=cels(celd)%cey ; cz=cels(celd)%cez
    dix=node(nods)%x-cx ; diy=node(nods)%y-cy ; diz=node(nods)%z-cz !>>>Miguel12-8-14

    if(ffu(19)==0.and.npag(nparam_per_node+11)>0)then                        !!>>HC 9-10-2020 POLARIZED ECM SECRETION	
      call random_number(a)                                                  !!>>HC 9-10-2020 RANDOM VECTOR (default)
      kch=int(a*nvaloq)+1          !                                         !!>>HC 9-10-2020
      dix=particions_esfera(kch,1) !                                         !!>>HC 9-10-2020
      diy=particions_esfera(kch,2) !                                         !!>>HC 9-10-2020  
      diz=particions_esfera(kch,3) !                                         !!>>HC 9-10-2020
      dch=node(nods)%eqd  !>>Miquel28-7-14                                   !!>>HC 9-10-2020
      rax=dix*dch ; ray=diy*dch ; raz=diz*dch                                !!>>HC 9-10-2020
      call random_number(a)                                                  !!>>HC 9-10-2020
      a=1.0d0+a*0.1d0                                                        !!>>HC 9-10-2020
      rax=rax*a ; ray=ray*a ; raz=raz*a                                      !!>>HC 9-10-2020 This is the default random vector

                                                                             !!>>HC 9-10-2020
      vxch=0 ; vych=0 ; vzch=0                                               !!>>HC 9-10-2020 GRADIENT VECTOR
      axch=node(nods)%x ; aych=node(nods)%y ; azch=node(nods)%z              !!>>HC 9-10-2020
      do jch=1,nneigh(nods)                                                  !!>>HC 9-10-2020 We first stablish the node polarity for ECM
        jjch=neigh(nods,jch)                                                 !!>>HC 9-10-2020
        ax=node(jjch)%x-axch ; ay=node(jjch)%y-aych ; az=node(jjch)%z-azch   !!>>HC 9-10-2020
        dch=1/sqrt(ax**2+ay**2+az**2)                                        !!>>HC 9-10-2020 wa 11 is used to establish ECM polarity when
        do kch=1,npag(nparam_per_node+11)                                    !!>>HC 9-10-2020 ffu(19)=1, we might want to create a new wa
          kkch=whonpag(nparam_per_node+11,kch)                               !!>>HC 9-10-2020
          if (gex(jjch,kkch)>0.0d0) then                                     !!>>HC 9-10-2020
            alfach=gex(jjch,kkch)*gen(kkch)%e(nparam_per_node+11)*dch        !!>>HC 9-10-2020
            vxch=vxch+ax*alfach ; vych=vych+ay*alfach ; vzch=vzch+az*alfach  !!>>HC 9-10-2020
          else                                                               !!>>HC 9-10-2020
            cycle                                                            !!>>HC 9-10-2020
          end if                                                             !!>>HC 9-10-2020
        end do                                                               !!>>HC 9-10-2020
      end do                                                                 !!>>HC 9-10-2020
      if(vxch==0.and.vych==0.and.vzch==0)then                                !!>>HC 9-10-2020
        dch=0                                                                !!>>HC 9-10-2020
      else                                                                   !!>>HC 9-10-2020
        dch=1/sqrt(vxch**2+vych**2+vzch**2)                                  !!>>HC 9-10-2020
      end if                                                                 !!>>HC 9-10-2020
      polxch=vxch*dch ; polych=vych*dch ; polzch=vzch*dch                    !!>>HC 9-10-2020 This is the gradient chemical vector

      dch=0.0d0                                                              !!>>HC 9-10-2020 HERE WE ASSESS THE IMPORTANCE    
      do jch=1,npag(nparam_per_node+11)                                      !!>>HC 9-10-2020 OF RANDOM VS CHEMICAL VECTORS 
         kch=whonpag(nparam_per_node+11,jch)                                 !!>>HC 9-10-2020
         dch=d+gex(nods,kch)*gen(kch)%e(nparam_per_node+11)                  !!>>HC 9-10-2020
      end do                                                                 !!>>HC 9-10-2020
      dch=1d0/(1d0+dch)  !ponderacion entre vector aleatortio i quimico      !!>>HC 9-10-2020
      !print*,"d ponderation",d,node(nd)%tipus                               !!>>HC 9-10-2020
      !d    ; dependence of the gradient vector                              !!>>HC 9-10-2020
      !dch=1  ; polarization vector is completely random (default mode)      !!>>HC 9-10-2020
      !dch=0  ; polarization vector comes only from the gradient             !!>>HC 9-10-2020
      !dch=0.5; polarization vector comes equally from both                  !!>>HC 9-10-2020    
      cxch=((1-dch)*polxch)+(dch*rax)                                        !!>>HC 9-10-2020 vector resultante para secreccion 
      cych=((1-dch)*polych)+(dch*ray)                                        !!>>HC 9-10-2020
      czch=((1-dch)*polzch)+(dch*raz)                                        !!>>HC 9-10-2020
      dch=1/sqrt(cxch**2+cych**2+czch**2)                                    !!>>HC 9-10-2020
      cxch=cxch*dch ; cych=cych*dch ; czch=czch*dch                          !!>>HC 9-10-2020
      node(nd)%x=node(nods)%x+desmax*cxch                                    !!>>HC 9-10-2020 ECM NODE COORDINATES
      node(nd)%y=node(nods)%y+desmax*cych                                    !!>>HC 9-10-2020
      node(nd)%z=node(nods)%z+desmax*czch                                    !!>>HC 9-10-2020


    else
       if(ffu(1)==1)then  !>>Miquel28-7-14
         b=0
         do i=1,cels(ii)%nunodes
           iii=cels(ii)%node(i)
           ax=node(iii)%x-cx ; ay=node(iii)%y-cy ; az=node(iii)%z-cz      
           a=dix*ax+diy*ay+diz*az
           if (a>b) then ; b=a ; bx=ax ; by=ay ; bz=az ; kk=iii ; end if
         end do            
         ax=bx ; ay=by ; az=bz             !>>>Miguel12-8-14
         !ax=bx*d ; ay=by*d ; az=bz*d      !>>>Miguel12-8-14 (comment)
       else                                !>>Miquel28-7-14
         call random_number(a)      !when mes. cell is 1 node the ECM node is secreted on a random direction !>>Miquel13-8-14
         k=int(a*nvaloq)+1          !
         dix=particions_esfera(k,1) !
         diy=particions_esfera(k,2) !
         diz=particions_esfera(k,3) !
         d=node(nods)%eqd  !>>Miquel28-7-14
         ax=dix*d ; ay=diy*d ; az=diz*d
       end if
       call random_number(a)
       a=1.0d0+a*0.1d0
       node(nd)%x=cx+ax*a ; node(nd)%y=cy+ay*a ; node(nd)%z=cz+az*a 
    endif      
  end if

  node(nd)%tipus=4

  node(nd)%you=0d0                !this has not effect
  node(nd)%eqd=ecmmax ; node(nd)%grd=ecmmax
!  node(nd)%codel=ecmmax
  node(nd)%rep=0d0
  node(nd)%eqs=0d0
  node(nd)%hoo=0d0
  node(nd)%erp=0d0
  node(nd)%est=0d0  !pfh 20-3-15
  node(nd)%altre=0
  node(nd)%icel=-nd
  node(nd)%dmo=desmax
  node(nd)%mov=temp
  node(nd)%marge=1
  node(nd)%orix=node(nd)%x
  node(nd)%oriy=node(nd)%y
  node(nd)%oriz=node(nd)%z
  node(nd)%ecm=0.0d0
  node(nd)%fix=0.0d0
  node(nd)%cod=0
  node(nd)%pld=0
  node(nd)%dif=0
  node(nd)%kfi=0
  node(nd)%pla=0
  node(nd)%kvol=0
  node(nd)%talone=0
  node(nd)%vod=0d0  !pfh 20-3-15
  node(nd)%adh=0d0   !pfh 20-3-15
  
  nneigh(nd)=0 !initialize the neighbors to not get errors in other parts
  neigh(nd,:)=0 !>>>Miquel14-3-14
  fmeanl(nd)=0d0
  
  a=0.0d0
  b=0.0d0 
  c=0.0d0
  do ii=1,npag(nparam_per_node+6)            !>>> Is 22-5-14
    i=whonpag(nparam_per_node+6,ii)          !>>> Is 22-5-14
    if (gex(nods,i)>0.0d0) then              !>>> Is 22-5-14
      if (gen(i)%e(nparam_per_node+6)>0) then              ! wa only 0-1, any value larger than 0 is like 1
        b=b+gex(nods,i)*gen(i)%e(nparam_per_node+6)!*0.5d0 ! wa in units of space-req
      end if                                 !>>> Is 22-5-14
    end if                                   !>>> Is 22-5-14
  end do
  !b=b!*delta !it's not a rate, so no delta !>>Miquel12-8-14

  node(nd)%rec=b !>>Miquel2-7-14
  !node(nods)%rec=node(nods)%rec-b  !>>> Is 22-5-15 otherwise repcel can be zero if no +5 or +6
  b=0.0d0
  do ii=1,npag(nparam_per_node+7)            !>>> Is 22-5-14
    i=whonpag(nparam_per_node+7,ii)          !>>> Is 22-5-14
    if (gex(nods,i)>0.0d0) then              !>>> Is 22-5-14
      if (gen(i)%e(nparam_per_node+7)>0) then              ! wa only 0-1, any value larger than 0 is like 1
        b=b+gex(nods,i)*gen(i)%e(nparam_per_node+7)!*0.5d0 ! wa in units of space-req
      end if                                 !>>> Is 22-5-14
    end if                                   !>>> Is 22-5-14
  end do                                     !>>> Is 22-5-14
  !b=b!*delta !it's not a rate, so no delta !>>Miquel12-8-14

  node(nd)%add=node(nd)%grd+b !>>Miquel2-7-14

  if (node(nd)%rec .le. 0) node(nd)%rec=50.0d0 !>>> HC 4-5-20 !>>> HC 7-10-20 >>> HC 18-2-21 default value
  !node(nd)%add=node(nd)%eqd+c 
  !node(nd)%adh=node(nd)%eqd+c
  !node(nd)%adh=node(nd)%eqd+c   !>>> Is 22-5-14 here we had that: node(nd)%adh=node(nd)%eqd+c that seemed a mistake
  !gene product composition
  agex(nd,:)=0d0 !just in case !>>Miquel7-8-14
  do i=1,npag(nparam_per_node+5)         !the ECM gene products are transfered to the new node
    k=whonpag(nparam_per_node+5,i)       
    a=gen(k)%e(nparam_per_node+5)                !>>> Is 22-5-14
    !if (gen(k)%e(nparam_per_node+5)>0.0) a=1.0d0 !>>> Is 22-5-14
    b=gex(nods,k)                                   !>>> Is 22-5-14 it think it was reversed !>>Miquel7-8-14 
    agex(nods,k)=b*(1.0d0-a) ! wa is in units of probability  !>>> Is 22-5-14
    agex(nd,k)=b*a                                        !>>> Is 22-5-14
  end do

  if(ffu(8)==0)then
    !put it into the boxes
    ic=boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))!position in boxes
    if(ic==0)then
      !the cube was empty
      boxes(nint(node(nd)%x*urv),nint(node(nd)%y*urv),nint(node(nd)%z*urv))=nd
      list(nd)=0
    else		
      !the cube was not empty so we search for the end																
      do
        icc=ic; ic=list(icc) ; if(ic==0)then ; list(icc)=nd ; list(nd)=0 ; exit ; end if
      end do
    end if
  end if

  nodeo(nd)=node(nd)

end subroutine ecm_secretion

end module ecm
