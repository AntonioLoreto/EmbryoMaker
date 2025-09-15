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
!***************  MODUL VEINATGE ********************************************
!***************************************************************************

module neighboring
use general

integer,public, allocatable  :: boxes(:,:,:),list(:)
integer,public               :: nboxes,iextre
integer,public, allocatable  :: borderface(:,:)  !>>Miquel3-2-14
integer,public               :: nborder
real*8 ,private              :: maxlen    ! >>> Is 7-2-15
real*8 ,parameter            :: arrel_de_dos=1.414213562d0 

contains

!**************************************************************************
subroutine extrem 
    ! troba quina es la cel que esta mes allunyada del centre del embryo
    extre=0.0d0 ; iextre=1 
    do i=1,nd
      a=sqrt((node(i)%x)**2+(node(i)%y)**2+(node(i)%z)**2)
      if (a>extre) then ; extre=a ; iextre=i ; end if ;
    end do
end subroutine extrem

!**************************************************************************
subroutine iniboxes
integer ic,ii,jj,kk
integer,allocatable:: cboxes(:,:,:)
   
    call extrem

    a=2*maxval(node(:nd)%add) !maximal interaction distance between mesenchymal nodes
    rv=a+1d-3
    urv=1.0d0/a

    nboxes=nint(extre*urv)+int(maxval(node(:nd)%dmo)+dmax)+1
    if (allocated(list)) deallocate(list)
    allocate(list(nda))	!>>Miquel 14-10-12
    list=0
    if (allocated(boxes)) deallocate(boxes)
    if (allocated(cboxes)) deallocate(cboxes)
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    allocate(cboxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    boxes=0
    cboxes=0
    do i=1,nd
      ii=nint(node(i)%x*urv);jj=nint(node(i)%y*urv);kk=nint(node(i)%z*urv)
      list(i)=boxes(ii,jj,kk)
      boxes(ii,jj,kk)=i
      cboxes(ii,jj,kk)=cboxes(ii,jj,kk)+1
    end do
    mnn_dyn=maxval(cboxes)
    deallocate(cboxes)  ! >>> IS 22-12-20
    
    mnn_dynam=mnn_dyn*27  !27 because it is the maximal number of boxes around a box (and then the maximal number of points can be larger than mnn_dyn by that  
    !if (mnn_dynam>600) then
    !  print *,"PANIC, some nodes have more than 20 neighbors, in fact they have",mnn_dynam
    !end if
end subroutine iniboxes

!**************************************************************************
subroutine iniboxes_p
integer ic,ii,jj,kk
integer,allocatable:: cboxes(:,:,:)

    nboxes=nint(extre*urv)+int(maxval(node(:nd)%dmo)+dmax)+2
    if (allocated(list)) deallocate(list)
    allocate(list(nda))	!>>Miquel 14-10-12
    list=0

    if (allocated(boxes)) deallocate(boxes)
    if (allocated(cboxes)) deallocate(cboxes)
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    allocate(cboxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    boxes=0
    cboxes=0
    do i=1,nd
      ii=nint(node(i)%x*urv);jj=nint(node(i)%y*urv);kk=nint(node(i)%z*urv)
      list(i)=boxes(ii,jj,kk)
      boxes(ii,jj,kk)=i
      cboxes(ii,jj,kk)=cboxes(ii,jj,kk)+1
    end do
    mnn_dyn=maxval(cboxes)
    deallocate(cboxes)  ! >>> IS 22-12-20  MIRAR SI CAL, POTSER HO PODEM ENDIVINAR SI ES SIGNIFICATIVAMENT MES RAPID
    
    mnn_dynam=mnn_dyn*27  !27 because it is the maximal number of boxes around a box (and then the maximal number of points can be larger than mnn_dyn by that  
    !if (mnn_dynam>600) then
    !  print *,"PANIC, some nodes have more than 20 neighbors, in fact they have",mnn_dynam
    !end if
    
end subroutine iniboxes_p

!**************************************************************************
subroutine iniboxesll
    call extrem

    nboxes=nint(extre*urv)+int(maxval(node(:nd)%dmo)+dmax)+1 !;if(nboxes<7) print*,"NBOXESll",nboxes,"extre",extre
    if (allocated(list)) deallocate(list)
    allocate(list(nda))	!>>Miquel 14-10-12
    list=0
    if (allocated(boxes)) deallocate(boxes)
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
    boxes=0

    do i=1,nd 
      ii=nint(node(i)%x*urv);jj=nint(node(i)%y*urv);kk=nint(node(i)%z*urv) !;print*,"ii",ii,"jj",jj,"kk",kk
      list(i)=boxes(ii,jj,kk)
      boxes(ii,jj,kk)=i                 
    end do
end subroutine iniboxesll

!************************************************

subroutine neighbor_build !this subroutine is called at each iteration to assess the neighbors of all the nodes !>>>Miquel24-2-14
  call extrem
  !maxlen=sqrt(2*maxval(node(:nd)%eqs))**2+(2*maxval(node(:nd)%add)**2) !this is the maximal interaction distance between epithelial nodes
  maxlen=2*maxval(node(:nd)%add)**2 !this is the maximal interaction distance between epithelial nodes
  rv=2*maxval(node(:nd)%add) !maximal interaction distance between mesenchymal nodes
  urv=1.0d0/rv

  !THIS IS JUST FOR CLARITY, TO SIMPLIFY EACH SUBROUTINE AND MAKE IT MORE READEABLE 
  if (ffu(9)==0) then              !!>> HC 14-1-2021
     if (ffu(8)==0) then !normal neighboring, extensive search of the boxes
        if (ffu(2)==0) then
           if (ffu(23)==1) then        !!>> HC 14-1-2021
              urv=1.d0/(arrel_de_dos*rv+1d-3) ! beware, that could produce troubles, somehow, I have not check it carefully
                                            ! COMMENT: that was mildly wrong (or just inefficient) in all previous versions
              call iniboxes_p
              if (ffu(24)==1) then
                 if (ffu(17)==0) then
                    call neighbor_build_complexneigh
                 else
                    call neighbor_build_complexneigh_inter_other_side_epi
                 end if
              else
                 print *,"The algorithm to recover neighbors is not implemented yet with extended neighbording"
                 stop
              end if
           else  
              urv=1.0d0/(rv+1d-3)
              call iniboxes_p
              if (ffu(17)==0) then
                 call neighbor_build_simpleneigh
              else
                 call neighbor_build_simpleneigh_inter_other_side_epi
              end if    
           end if
        else
          print *,"SORRY,..., GABRIEL IS NOT IMPLEMENTED WITH THOSE FFUs, BUT IT IS EASY TO DO IF YOU WANT"
          stop
        end if
     else
        urv=1.0d0/(rv+1d-3)
        call iniboxes_p
        if (ffu(2)==0) then !no GABRIEL
          call neighbor_build_triang_nogabriel    
        else
          call neighbor_build_triang_gabriel
        end if
     end if 
  else
     call neighbor_build_old
  end if
end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_node

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_node(i)
  integer i
  real*8 a

  call extrem
  !maxlen=sqrt((2*maxval(node(:nd)%eqs))**2+(2*maxval(node(:nd)%add)**2)) !this is the maximal interaction distance between epithelial nodes
  maxlen=2*maxval(node(:nd)%add)**2 !this is the maximal interaction distance between epithelial nodes
  a=2*maxval(node(:nd)%add) !maximal interaction distance between mesenchymal nodes
  rv=a
  urv=1.0d0/a

  !THIS IS JUST FOR CLARITY, TO SIMPLIFY EACH SUBROUTINE AND MAKE IT MORE READEABLE 
  if (ffu(9)==0) then        !!>> HC 14-1-2021
    if (ffu(8)==0) then !normal neighboring, extensive search of the boxes
      if (ffu(2)==0) then
        if (ffu(23)==1) then  !!>> HC 14-1-2021
          urv=1.d0/(arrel_de_dos*rv+1d-3)   ! COMMENT: that was mildly wrong (or just inefficient) in all previous versions   
          call iniboxes_p
          if (ffu(17)==0) then
            call neighbor_build_complexneigh_node(i)
          else
            call neighbor_build_complexneigh_inter_other_side_epi_node(i)
          end if
        else  
          urv=1.0d0/(rv+1d-3)      
          call iniboxes_p
          if (ffu(17)==0) then
            call neighbor_build_simpleneigh_node(i)
          else
            call neighbor_build_simpleneigh_inter_other_side_epi_node(i)
          end if    
        end if
      else
        print *,"SORRY,..., GABRIEL IS NOT IMPLEMENTED WITH THOSE FFUs, BUT IT IS EASY TO DO IF YOU WANT"
        stop
      end if
    else
      urv=1.0d0/(rv+1d-3)    
      call iniboxes_p
      if (ffu(2)==0) then !no GABRIEL
        call neighbor_build_node_old(i)  ! this may need to be changed if this is used often, that is probably not    
      else
        call neighbor_build_node_old(i)  ! this may need to be changed if this is used often, that is probably not    
      end if
    end if 
  else
    call neighbor_build_node_old(i)
  end if
end subroutine 



!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_COMPLEXNEIGH    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_complexneigh        ! COMPLEX NEIGHBORHOOD, EPITELIA PUSH NORMAL TO THEIR SURFACE, NO INTERACTION BETWEEN NODES FROM DIFFERENT SIDES OF EPITELIUM
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie
real*8:: ix,iy,iz,dist,udist

integer::tipi  !>>Miquel31-12-14
real*8::dai,sdai    !>>Miquel31-12-14
real*8::sqad(nd),sda(nd)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8 sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh,ti


    sqad(:nd)=node(:nd)%add*arrel_de_dos
    sda(:nd)=node(:nd)%add    
    omnn=0
    do i=1,nd
      ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
      ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
      ivv=node(i)%altre

      snneigh=0
      
      tipi=node(i)%tipus
      dai=node(i)%add
      sdai=sqad(i)
      ti=1
      
      if (tipi<3) then      !23-1-21 Isaac, l'ordre dels loops podria afecta la eficiciencia, mirar
        if (tipi==1) ti=2   
        do i1=-1,1,2 
          iii1=ii1+i1
          do i2=-1,1
            iii2=ii2+i2
            do i3=-1,1
              !iii3=ii3+i3
              ie=boxes(ii3+i3,iii2,iii1)
              do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
              end do
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        do i2=-1,1,2
          iii2=ii2+i2
          do i3=-1,1
            !iii3=ii3+i3
            ie=boxes(ii3+i3,iii2,ii1)
            do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
            
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        do i3=-1,1,2
          !iii3=ii3+i3
          ie=boxes(ii3+i3,ii2,ii1)
          do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)

          end do
        end do
        
        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        !i3=0
        !do i3=-1,1,2
        !iii3=ii3+i3
        
        ! this is the central box
        ie=boxes(ii3,ii2,ii1)
        do while(ie.ne.0)
          k=node(ie)%tipus
          if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
          if (ie==i) then ; ie=list(ie) ; cycle ; end if   
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)

        end do

      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i1=-1,1,2 
          iii1=ii1+i1
          do i2=-1,1
            iii2=ii2+i2
            do i3=-1,1
              !iii3=ii3+i3
              ie=boxes(ii3+i3,iii2,iii1)
              do while(ie.ne.0)
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
                  !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
                end if
                if (ie==i) then ; ie=list(ie) ; cycle ; end if
                 
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
              end do
            end do
          end do
        end do       

        !i1=0
        !iii1=ii1+i1
        do i2=-1,1,2
          iii2=ii2+i2
          do i3=-1,1
            !iii3=ii3+i3
            ie=boxes(ii3+i3,iii2,ii1)
            do while(ie.ne.0)
              dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
              if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
              else  !mesench/ECM vs epithelial
                if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
                !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
              end if
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
               
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              sdneigh(snneigh)=sqrt(dist)
              ie=list(ie)           
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        do i3=-1,1,2
          !iii3=ii3+i3
          ie=boxes(ii3+i3,ii2,ii1)
          do while(ie.ne.0)
            dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
            if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
              if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
            else  !mesench/ECM vs epithelial
              if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
              !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
            end if
            if (ie==i) then ; ie=list(ie) ; cycle ; end if
             
            snneigh=snneigh+1
            sneigh(snneigh)=ie
            sdneigh(snneigh)=sqrt(dist)
            ie=list(ie)          
          end do
        end do
        
        ie=boxes(ii3,ii2,ii1)
        do while(ie.ne.0)
          dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
          if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
            if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
          else  !mesench/ECM vs epithelial
            if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
            !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
          end if
          if (ie==i) then ; ie=list(ie) ; cycle ; end if
             
          snneigh=snneigh+1
          sneigh(snneigh)=ie
          sdneigh(snneigh)=sqrt(dist)
          ie=list(ie)                
        end do
      end if

      ! here we transfer the neighbors of i to a pre-neigh neighbor
    if(snneigh>omnn) omnn=snneigh
    nneigh(i)=snneigh
    trans_neigh(i,:omnn)=sneigh(:omnn)
    trans_dneigh(i,:omnn)=sdneigh(:omnn)     
  end do
    
  if(allocated(neigh)) deallocate(neigh)
  if(allocated(dneigh)) deallocate(dneigh)
  allocate(neigh(nda,omnn),dneigh(nda,omnn))
  neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
  dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
 end subroutine 

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_COMPLEXNEIGH_inter_other_side_epi    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_complexneigh_inter_other_side_epi      ! COMPLEX NEIGHBORHOOD AND NODES FROM DIFFERENT EPITELIAL SIDES DOOOO  INTERACT
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie
real*8:: ix,iy,iz,dist,udist

integer::tipi  !>>Miquel31-12-14
real*8::dai,sdai    !>>Miquel31-12-14
real*8::sqad(nd),sda(nd)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8 sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh,ti


    sqad(:nd)=node(:nd)%add*arrel_de_dos
    sda(:nd)=node(:nd)%add    
    omnn=0
    do i=1,nd
      ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
      ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
      ivv=node(i)%altre

      snneigh=0
      
      tipi=node(i)%tipus
      dai=node(i)%add
      sdai=sqad(i)
      ti=1
      
      if (tipi<3) then      !23-1-21 Isaac, l'ordre dels loops podria afecta la eficiciencia, mirar
        if (tipi==1) ti=2   
        do i1=-1,1,2 
          iii1=ii1+i1
          do i2=-1,1
            iii2=ii2+i2
            do i3=-1,1
              !iii3=ii3+i3
              ie=boxes(ii3+i3,iii2,iii1)
              do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+node(ie)%add)**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
              end do
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        do i2=-1,1,2
          iii2=ii2+i2
          do i3=-1,1
            !iii3=ii3+i3
            ie=boxes(ii3+i3,iii2,ii1)
            do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+node(ie)%add)**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
            
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        do i3=-1,1,2
          !iii3=ii3+i3
          ie=boxes(ii3+i3,ii2,ii1)
          do while(ie.ne.0)
                k=node(ie)%tipus
                if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+node(ie)%add)**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)

          end do
        end do
        
        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        !i3=0
        !do i3=-1,1,2
        !iii3=ii3+i3
        
        ! this is the central box
        ie=boxes(ii3,ii2,ii1)
        do while(ie.ne.0)
          k=node(ie)%tipus
          if (k==ti) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
          if (ie==i) then ; ie=list(ie) ; cycle ; end if   
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(tipi==k)then !same face epithelials
                  !b=(sdai+sqad(ie))**2 ! the square root of 2 comes because we may interact planarly (upwards from the round face of the cylinder                
                  if(dist>(sdai+sqad(ie))**2) then ; ie=list(ie) ; cycle ; end if
                else
                  !b=(sdai+node(ie)%add)**2 ! here we do not need the square root of 2                                
                  if(dist>(sdai+node(ie)%add)**2)then ; ie=list(ie) ; cycle ; end if    
                end if
                  
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)

        end do

      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i1=-1,1,2 
          iii1=ii1+i1
          do i2=-1,1
            iii2=ii2+i2
            do i3=-1,1
              !iii3=ii3+i3
              ie=boxes(ii3+i3,iii2,iii1)
              do while(ie.ne.0)
                dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
                  !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
                end if
                if (ie==i) then ; ie=list(ie) ; cycle ; end if
                 
                snneigh=snneigh+1
                sneigh(snneigh)=ie
                sdneigh(snneigh)=sqrt(dist)
                ie=list(ie)
              end do
            end do
          end do
        end do       

        !i1=0
        !iii1=ii1+i1
        do i2=-1,1,2
          iii2=ii2+i2
          do i3=-1,1
            !iii3=ii3+i3
            ie=boxes(ii3+i3,iii2,ii1)
            do while(ie.ne.0)
              dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
              if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
              else  !mesench/ECM vs epithelial
                if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
                !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
              end if
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
               
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              sdneigh(snneigh)=sqrt(dist)
              ie=list(ie)           
            end do
          end do
        end do

        !i1=0
        !iii1=ii1+i1
        !i2=0
        !iii2=ii2+i2
        do i3=-1,1,2
          !iii3=ii3+i3
          ie=boxes(ii3+i3,ii2,ii1)
          do while(ie.ne.0)
            dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
            if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
              if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
            else  !mesench/ECM vs epithelial
              if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
              !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
            end if
            if (ie==i) then ; ie=list(ie) ; cycle ; end if
             
            snneigh=snneigh+1
            sneigh(snneigh)=ie
            sdneigh(snneigh)=sqrt(dist)
            ie=list(ie)          
          end do
        end do
        
        ie=boxes(ii3,ii2,ii1)
        do while(ie.ne.0)
          dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
          if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
            if(dist>(dai+sda(ie))**2)then ; ie=list(ie) ; cycle ; end if
          else  !mesench/ECM vs epithelial
            if(dist>(dai+sqad(ie))**2)then ; ie=list(ie) ; cycle ; end if
            !if(dist>arrel_de_dos*a)then ; ie=list(ie) ; cycle ; end if
          end if
          if (ie==i) then ; ie=list(ie) ; cycle ; end if
             
          snneigh=snneigh+1
          sneigh(snneigh)=ie
          sdneigh(snneigh)=sqrt(dist)
          ie=list(ie)                
        end do
      end if

      ! here we transfer the neighbors of i to a pre-neigh neighbor
    if(snneigh>omnn) omnn=snneigh
    nneigh(i)=snneigh
    trans_neigh(i,:omnn)=sneigh(:omnn)
    trans_dneigh(i,:omnn)=sdneigh(:omnn)     
  end do
    
  if(allocated(neigh)) deallocate(neigh)
  if(allocated(dneigh)) deallocate(dneigh)
  allocate(neigh(nda,omnn),dneigh(nda,omnn))
  neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
  dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)

end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_SIMPLENEIGH    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_simpleneigh            ! ALL NEIGHBORS INTERACT IN THE SAME WAY, EXCEPT WE DO NOT ALLOW INTERACTION OF NODES FROM DIFFERENT SIDES OF THE EPITHELIUM

integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8 sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh

  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
!  neigh=0  ! 4-3-2020
!  nneigh=0 ! 4-3-2020
  omnn=0

  ! main node loop in this subroutine
  do i=1,nd
      
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    !if(node(i)%tipus<3)then
    !  snneigh=1
    !  sneigh(1)=ivv
    !  sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    !end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
!          iii3=ii3+i3
          ie=boxes(ii3+i3,iii2,iii1)
          do while(ie.ne.0)
            !if (ie==i) then ; ie=list(ie) ; cycle ; end if
            !if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
            if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
            if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
            dist=(node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2
            a=(dai+node(ie)%add)**2
            if(dist>a)then ; ie=list(ie) ; cycle ; end if             
            if (ie==i) then ; ie=list(ie) ; cycle ; end if                 

            snneigh=snneigh+1
            sneigh(snneigh)=ie
            sdneigh(snneigh)=sqrt(dist)
            ie=list(ie)
          end do
        end do
      end do
    end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if(snneigh>omnn) omnn=snneigh
     nneigh(i)=snneigh
     trans_neigh(i,:omnn)=sneigh(:omnn)
     trans_dneigh(i,:omnn)=sdneigh(:omnn)     
    end do  
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)

end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_SIMPLE_inter_other_side_epi    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
subroutine neighbor_build_simpleneigh_inter_other_side_epi            ! ALL NEIGHBORS INTERACT IN THE SAME WAY, EXCEPT WE DO NOT ALLOW INTERACTION OF NODES FROM DIFFERENT SIDES OF THE EPITHELIUM


integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8 sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh

  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
  neigh=0  ! 4-3-2020
  nneigh=0 ! 4-3-2020
  omnn=0

  ! main node loop in this subroutine
  do i=1,nd
      
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    !if(node(i)%tipus<3)then
    !  snneigh=1
    !  sneigh(1)=ivv
    !  sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    !end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
          iii3=ii3+i3
          ie=boxes(iii3,iii2,iii1)
          do while(ie.ne.0)
            !if (ie==i) then ; ie=list(ie) ; cycle ; end if
            !if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
            !if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
            !if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
            dist=((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
            a=(dai+node(ie)%add)**2
            if(dist>a)then ; ie=list(ie) ; cycle ; end if             
            if (ie==i) then ; ie=list(ie) ; cycle ; end if                 

            snneigh=snneigh+1
            sneigh(snneigh)=ie
            sdneigh(snneigh)=sqrt(dist)
            ie=list(ie)
          end do
        end do
      end do
    end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if(snneigh>omnn) omnn=snneigh
     nneigh(i)=snneigh
     trans_neigh(i,:omnn)=sneigh(:omnn)
     trans_dneigh(i,:omnn)=sdneigh(:omnn)     
    end do
    
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)

end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_triand_gabriel

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_triang_gabriel
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14


!!!triangulation variables
integer:: npt !number of points
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer,allocatable::oosneigh(:)
real*8,allocatable::oosdneigh(:)
integer::snneigh


!  a=maxval(node(:nd)%eqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%add)*dmax

  !omnn=0
  call iniboxes
  mnn_dynam=mnn_dyn*(2*nint(rdiffmax*urv)+1)**3
  allocate(trans_neigh(nd,mnn_dynam),trans_dneigh(nd,mnn_dynam))

  !*****3D triangulation neighooring***************************************************************
    
    npt=nd
!    sizht=3/2*npt
    sizht=3*npt
    maxfc=npt**2
    maxbf=npt**2
    allocate(vcl(3,npt),vm(npt))
    allocate(bf(1:3,maxbf),fc(1:7,maxfc))
    allocate(ht(0:sizht-1))
    
    
    !call extrem !not necessary, but sets some variables that later may be used in pinta and creix !>>Miquel21-3-14
    
    
    do i=1,npt  !building the input arrays for the triangualtion function
      vcl(1,i)=node(i)%x ; vcl(2,i)=node(i)%y ; vcl(3,i)=node(i)%z
      vm(i)=i
    end do !;print*,"pillat pre delau"
    !calling the triangulation subroutine
    call dtriw3(npt, sizht, maxbf, maxfc, vcl, vm, nbf, nfc, nface, ntetra, bf, fc, ht, ierr)
    !if(ierr/=0)then; print*,"error in the triangulation, avort or something",ierr,getot ; endif !call exit(24); end if
    !translating the ouptut into a neighbor matrix

    do i=1,nd !initialize neighbor matrix
!      neigh(i,1:nneigh(i))=0
      neigh(i,1:)=0  !IS 29-4-14
      nneigh(i)=0
    end do
    do j=1,nfc               ! it passes through all triangles 
      if((fc(1,j)>0).and.(fc(2,j)>0).and.(fc(3,j)>0))then ! it is a "valid" triangle
        ii=fc(1,j) ; jj=fc(2,j) ; kk=fc(3,j)
        iii=vm(ii);jjj=vm(jj);kkk=vm(kk) !; print*,"iii jjj kkk",iii,jjj,kkk
        tiii=node(iii)%tipus ; tjjj=node(jjj)%tipus ; tkkk=node(kkk)%tipus

        ivv=nneigh(iii)
        ij=0 ; ik=0 ; jk=0
        !connection iii-jjj*****
        if(ivv==0)then
          d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
          ivv=1
          trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
        else
          do i=1,nneigh(iii)
            if(trans_neigh(iii,i)==jjj)then; ij=1;exit;end if
          end do
          if(ij==0)then
            d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
            ivv=ivv+1
            trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
          end if
        end if
        if(ij==0)then
          !connection jjj-iii*******
          if(nneigh(jjj)==0)then ; nneigh(jjj)=1 ; trans_neigh(jjj,1)=iii ; trans_dneigh(jjj,1)=d
          else; iv=nneigh(jjj)+1 ; trans_neigh(jjj,iv)=iii ; trans_dneigh(jjj,iv)=d ; nneigh(jjj)=iv ; end if
        end if
        !connection iii-kkk*****
        do i=1,nneigh(iii)
          if(trans_neigh(iii,i)==kkk)then; ik=1;exit;end if
        end do
        if(ik==0)then
          d=sqrt((vcl(1,iii)-vcl(1,kkk))**2+(vcl(2,iii)-vcl(2,kkk))**2+(vcl(3,iii)-vcl(3,kkk))**2)
          ivv=ivv+1
          trans_neigh(iii,ivv)=kkk ; trans_dneigh(iii,ivv)=d
          !connection kkk-iii*******
          if(nneigh(kkk)==0)then ; nneigh(kkk)=1 ; trans_neigh(kkk,1)=iii ; trans_dneigh(kkk,1)=d
          else; iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=iii ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv ; end if
        end if
        nneigh(iii)=ivv
        !connection jjj-kkk*****
        do i=1,nneigh(jjj)
          if(trans_neigh(jjj,i)==kkk)then; jk=1;exit;end if
        end do
        if(jk==0)then
          d=sqrt((vcl(1,jjj)-vcl(1,kkk))**2+(vcl(2,jjj)-vcl(2,kkk))**2+(vcl(3,jjj)-vcl(3,kkk))**2)
          iv=nneigh(jjj)+1
          trans_neigh(jjj,iv)=kkk ; trans_dneigh(jjj,iv)=d ;nneigh(jjj)=iv
          !connection kkk-jjj*******
          iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=jjj ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv
        end if
      end if
    end do


    if(ffu(2)==1)then
    
      do i=1,nd
        !Screening by Gabriel graph !>>Miquel6-3-14
        !A neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
        !Ordering the neighors with increasing distance

        ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   


        !sorting algorithm by selection, it's ok
        do j=1,nneigh(i)-1  
          b=trans_dneigh(i,j)
          ii=0
          do k=j+1,nneigh(i)
            c=trans_dneigh(i,k)
             if(b>c)then
              ii=k ; b=trans_dneigh(i,k)
            end if
          end do
          if(ii/=0)then
            !jj=osneigh(j)
            kk=trans_neigh(i,ii) ; c=trans_dneigh(i,ii) !the swap
            trans_neigh(i,ii)=trans_neigh(i,j) ; trans_dneigh(i,ii)=trans_dneigh(i,j)
            trans_neigh(i,j)=kk ; trans_dneigh(i,j)=c
          end if
        end do
        
        !the screening
        ii=0 !the number of eliminated connections
        do j=nneigh(i),1,-1
          jj=trans_neigh(i,j)
          if(jj==ivv) cycle
          a=trans_dneigh(i,j)*0.5*screen_radius !the radius of the sphere
          jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
          cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
          do k=j-1,1,-1
            kk=trans_neigh(i,k)
            d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
            if(d<a)then !there is one node within the sphere, we must delete this connection
              do l=j,nneigh(i)-ii-1
                trans_neigh(i,l)=trans_neigh(i,l+1)
                trans_dneigh(i,l)=trans_dneigh(i,l+1)
              end do
              trans_neigh(i,nneigh(i)-ii)=0
              trans_dneigh(i,nneigh(i)-ii)=0
              ii=ii+1
              exit
            end if
          end do
        end do
        nneigh(i)=nneigh(i)-ii
      end do
    !else                                                 !>>> Is 23-4-15
    !  trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
    !  trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15      
    end if
    
    omnn=0
    do i=1,nd
      if(nneigh(i)>omnn) omnn=nneigh(i)
    end do
    !print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
!!do i=1,nd
!end do  



if (ffu(9)==0) then
  if(allocated(trans_neigh)) deallocate(trans_neigh)
  if(allocated(trans_dneigh)) deallocate(trans_dneigh)
  if(allocated(sdneigh)) deallocate(sdneigh)
  if(allocated(osdneigh)) deallocate(osdneigh)
  if(allocated(osneigh)) deallocate(osneigh)
  if(allocated(sneigh)) deallocate(sneigh)
  if(allocated(oosdneigh)) deallocate(oosdneigh)
  if(allocated(oosneigh)) deallocate(oosneigh)
  if(allocated(vcl)) deallocate(vcl)
  if(allocated(vm)) deallocate(vm)
  if(allocated(bf)) deallocate(bf)
  if(allocated(fc)) deallocate(fc)
  if(allocated(ht)) deallocate(ht)
end if
end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_triang_nogabriel

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_triang_nogabriel
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14


!!!triangulation variables
integer:: npt !number of points
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
!integer,dimension(:)::sneigh(mnn)
!real*8,dimension(:)::sdneigh(mnn)
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer,allocatable::oosneigh(:)
real*8,allocatable::oosdneigh(:)
integer::snneigh


!  a=maxval(node(:nd)%eqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%add)*dmax

  !omnn=0
  call iniboxes
  mnn_dynam=mnn_dyn*(2*nint(rdiffmax*urv)+1)**3
  allocate(trans_neigh(nd,mnn_dynam),trans_dneigh(nd,mnn_dynam))

  !*****3D triangulation neighooring***************************************************************
    
    npt=nd
!    sizht=3/2*npt
    sizht=3*npt
    maxfc=npt**2
    maxbf=npt**2
    allocate(vcl(3,npt),vm(npt))
    allocate(bf(1:3,maxbf),fc(1:7,maxfc))
    allocate(ht(0:sizht-1))
    
    
    !call extrem !not necessary, but sets some variables that later may be used in pinta and creix !>>Miquel21-3-14
    
    
    do i=1,npt  !building the input arrays for the triangualtion function
      vcl(1,i)=node(i)%x ; vcl(2,i)=node(i)%y ; vcl(3,i)=node(i)%z
      vm(i)=i
    end do !;print*,"pillat pre delau"
    !calling the triangulation subroutine
    call dtriw3(npt, sizht, maxbf, maxfc, vcl, vm, nbf, nfc, nface, ntetra, bf, fc, ht, ierr)
    !if(ierr/=0)then; print*,"error in the triangulation, avort or something",ierr,getot ; endif !call exit(24); end if
    !translating the ouptut into a neighbor matrix

    do i=1,nd !initialize neighbor matrix
!      neigh(i,1:nneigh(i))=0
      neigh(i,1:)=0  !IS 29-4-14
      nneigh(i)=0
    end do
    do j=1,nfc               ! it passes through all triangles 
      if((fc(1,j)>0).and.(fc(2,j)>0).and.(fc(3,j)>0))then ! it is a "valid" triangle
        ii=fc(1,j) ; jj=fc(2,j) ; kk=fc(3,j)
        iii=vm(ii);jjj=vm(jj);kkk=vm(kk) !; print*,"iii jjj kkk",iii,jjj,kkk
        tiii=node(iii)%tipus ; tjjj=node(jjj)%tipus ; tkkk=node(kkk)%tipus

        ivv=nneigh(iii)
        ij=0 ; ik=0 ; jk=0
        !connection iii-jjj*****
        if(ivv==0)then
          d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
          ivv=1
          trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
        else
          do i=1,nneigh(iii)
            if(trans_neigh(iii,i)==jjj)then; ij=1;exit;end if
          end do
          if(ij==0)then
            d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
            ivv=ivv+1
            trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
          end if
        end if
        if(ij==0)then
          !connection jjj-iii*******
          if(nneigh(jjj)==0)then ; nneigh(jjj)=1 ; trans_neigh(jjj,1)=iii ; trans_dneigh(jjj,1)=d
          else; iv=nneigh(jjj)+1 ; trans_neigh(jjj,iv)=iii ; trans_dneigh(jjj,iv)=d ; nneigh(jjj)=iv ; end if
        end if
        !connection iii-kkk*****
        do i=1,nneigh(iii)
          if(trans_neigh(iii,i)==kkk)then; ik=1;exit;end if
        end do
        if(ik==0)then
          d=sqrt((vcl(1,iii)-vcl(1,kkk))**2+(vcl(2,iii)-vcl(2,kkk))**2+(vcl(3,iii)-vcl(3,kkk))**2)
          ivv=ivv+1
          trans_neigh(iii,ivv)=kkk ; trans_dneigh(iii,ivv)=d
          !connection kkk-iii*******
          if(nneigh(kkk)==0)then ; nneigh(kkk)=1 ; trans_neigh(kkk,1)=iii ; trans_dneigh(kkk,1)=d
          else; iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=iii ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv ; end if
        end if
        nneigh(iii)=ivv
        !connection jjj-kkk*****
        do i=1,nneigh(jjj)
          if(trans_neigh(jjj,i)==kkk)then; jk=1;exit;end if
        end do
        if(jk==0)then
          d=sqrt((vcl(1,jjj)-vcl(1,kkk))**2+(vcl(2,jjj)-vcl(2,kkk))**2+(vcl(3,jjj)-vcl(3,kkk))**2)
          iv=nneigh(jjj)+1
          trans_neigh(jjj,iv)=kkk ; trans_dneigh(jjj,iv)=d ;nneigh(jjj)=iv
          !connection kkk-jjj*******
          iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=jjj ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv
        end if
      end if
    end do


    omnn=0
    do i=1,nd
      if(nneigh(i)>omnn) omnn=nneigh(i)
    end do
    !print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)

if (ffu(9)==0) then
  if(allocated(trans_neigh)) deallocate(trans_neigh)
  if(allocated(trans_dneigh)) deallocate(trans_dneigh)
  if(allocated(sdneigh)) deallocate(sdneigh)
  if(allocated(osdneigh)) deallocate(osdneigh)
  if(allocated(osneigh)) deallocate(osneigh)
  if(allocated(sneigh)) deallocate(sneigh)
  if(allocated(oosdneigh)) deallocate(oosdneigh)
  if(allocated(oosneigh)) deallocate(oosneigh)
  if(allocated(vcl)) deallocate(vcl)
  if(allocated(vm)) deallocate(vm)
  if(allocated(bf)) deallocate(bf)
  if(allocated(fc)) deallocate(fc)
  if(allocated(ht)) deallocate(ht)
end if
end subroutine



!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_node

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************


!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_COMPLEXNEIGH_node    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_complexneigh_node(i)         ! COMPLEX NEIGHBORHOOD, EPITELIA PUSH NORMAL TO THEIR SURFACE, NO INTERACTION BETWEEN NODES FROM DIFFERENT SIDES OF EPITELIUM
integer i
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8  sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh


  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
  ! main node loop in this subroutine
      
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
          iii3=ii3+i3
          ie=boxes(iii3,iii2,iii1)
           do while(ie.ne.0)
             if (ie==i) then ; ie=list(ie) ; cycle ; end if
             if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
             if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
             a=dai+node(ie)%add
             if(tipi>=3)then
               if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                 if(dist>a)then ; ie=list(ie) ; cycle ; end if
               else  !mesench/ECM vs epithelial
                 b=arrel_de_dos*a
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if
               end if
             else !epithelial vs epithelial
               b=arrel_de_dos*a
               if(tipi==node(ie)%tipus)then !same face epithelials
                 if(b<maxlen) b=maxlen
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if
               else
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if    
               end if
             end if
                  
             !>>Miquel31-12-14
             snneigh=snneigh+1
             sneigh(snneigh)=ie
             sdneigh(snneigh)=dist
             ie=list(ie)
           end do
         end do
       end do
     end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if (snneigh>omnn)then                                !!>> HC 15-1-2021 the node i has increased the maximum number of neighbors
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We save the original data (we need to change the dimension of neigh and dneigh)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        omnn=snneigh                                      !!>> HC 15-1-2021
        if(allocated(neigh)) deallocate(neigh)            !!>> HC 15-1-2021 We realocate the matrixes neigh and dneigh to change their dimension
        if(allocated(dneigh)) deallocate(dneigh)          !!>> HC 15-1-2021
        allocate(neigh(nda,omnn),dneigh(nda,omnn))        !!>> HC 15-1-2021
        neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We recover the original data
        dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 Here we add the new neighbors of the node i
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     else                                                 !!>> HC 15-1-2021 In this case the node i has the same or less neighbors than the rest of nodes
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 we add the new neighs to the matrixes neigh and dneigh directly
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     endif                                                !!>> HC 15-1-2021
  
end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_COMPLEXNEIGH_inter_other_side_epi_node    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_complexneigh_inter_other_side_epi_node(i)      ! COMPLEX NEIGHBORHOOD AND NODES FROM DIFFERENT EPITELIAL SIDES DOOOO  INTERACT
integer i
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8  sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh


  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
  ! main node loop in this subroutine
      
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
          iii3=ii3+i3
          ie=boxes(iii3,iii2,iii1)
           do while(ie.ne.0)
             if (ie==i) then ; ie=list(ie) ; cycle ; end if
             if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
             !if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             !if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
             a=dai+node(ie)%add
             if(tipi>=3)then
               if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                 if(dist>a)then ; ie=list(ie) ; cycle ; end if
               else  !mesench/ECM vs epithelial
                 b=arrel_de_dos*a
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if
               end if
             else !epithelial vs epithelial
               b=arrel_de_dos*a
               if(tipi==node(ie)%tipus)then !same face epithelials
                 if(b<maxlen) b=maxlen
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if
               else
                 if(dist>b)then ; ie=list(ie) ; cycle ; end if    
               end if
             end if
                  
             !>>Miquel31-12-14
             snneigh=snneigh+1
             sneigh(snneigh)=ie
             sdneigh(snneigh)=dist
             ie=list(ie)
           end do
         end do
       end do
     end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if (snneigh>omnn)then                                !!>> HC 15-1-2021 the node i has increased the maximum number of neighbors
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We save the original data (we need to change the dimension of neigh and dneigh)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        omnn=snneigh                                      !!>> HC 15-1-2021
        if(allocated(neigh)) deallocate(neigh)            !!>> HC 15-1-2021 We realocate the matrixes neigh and dneigh to change their dimension
        if(allocated(dneigh)) deallocate(dneigh)          !!>> HC 15-1-2021
        allocate(neigh(nda,omnn),dneigh(nda,omnn))        !!>> HC 15-1-2021
        neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We recover the original data
        dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 Here we add the new neighbors of the node i
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     else                                                 !!>> HC 15-1-2021 In this case the node i has the same or less neighbors than the rest of nodes
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 we add the new neighs to the matrixes neigh and dneigh directly
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     endif                                                !!>> HC 15-1-2021

end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_SIMPLENEIGH_node(i=    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

subroutine neighbor_build_simpleneigh_node(i)            ! ALL NEIGHBORS INTERACT IN THE SAME WAY, EXCEPT WE DO NOT ALLOW INTERACTION OF NODES FROM DIFFERENT SIDES OF THE EPITHELIUM
integer i
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8  sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh


  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
          iii3=ii3+i3
          ie=boxes(iii3,iii2,iii1)
           do while(ie.ne.0)
             if (ie==i) then ; ie=list(ie) ; cycle ; end if
             if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
             !if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             !if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
             a=dai+node(ie)%add
             if(dist>a)then ; ie=list(ie) ; cycle ; end if             
                  
             !>>Miquel31-12-14
             snneigh=snneigh+1
             sneigh(snneigh)=ie
             sdneigh(snneigh)=dist
             ie=list(ie)
           end do
         end do
       end do
     end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if (snneigh>omnn)then                                !!>> HC 15-1-2021 the node i has increased the maximum number of neighbors
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We save the original data (we need to change the dimension of neigh and dneigh)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        omnn=snneigh                                      !!>> HC 15-1-2021
        if(allocated(neigh)) deallocate(neigh)            !!>> HC 15-1-2021 We realocate the matrixes neigh and dneigh to change their dimension
        if(allocated(dneigh)) deallocate(dneigh)          !!>> HC 15-1-2021
        allocate(neigh(nda,omnn),dneigh(nda,omnn))        !!>> HC 15-1-2021
        neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We recover the original data
        dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 Here we add the new neighbors of the node i
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     else                                                 !!>> HC 15-1-2021 In this case the node i has the same or less neighbors than the rest of nodes
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 we add the new neighs to the matrixes neigh and dneigh directly
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     endif                                                !!>> HC 15-1-2021

end subroutine

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************

           !NEIGHBOR_BUILT_SIMPLE_inter_other_side_epi_node(i)    IS 22-12-20

!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
!*************************************************************************************************************************************************************************
subroutine neighbor_build_simpleneigh_inter_other_side_epi_node(i)            ! ALL NEIGHBORS INTERACT IN THE SAME WAY, EXCEPT WE DO NOT ALLOW INTERACTION OF NODES FROM DIFFERENT SIDES OF THE EPITHELIUM
integer i
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai    !>>Miquel31-12-14
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer sneigh(mnn_dynam),trans_neigh(nd,mnn_dynam)
real*8  sdneigh(mnn_dynam),trans_dneigh(nd,mnn_dynam)
integer::snneigh


  ! ASSUMPTION: WE ASSUME THAT BOXES ARE LARGE ENOUGH SO THAT DIFFUSION CAN ONLY HAPPEN TO NEIGHBOR BOXES (THAT IS REALISTIC IN MOST CASES AS LONG AS THERE ARE NO LARGE CAVITIIES
  
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre

    ! temporal neigh and dneigh
    snneigh=0
    sneigh=0 
    sdneigh=0d0

    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    tipi=node(i)%tipus
    dai=node(i)%add
    do i1=-1,1 
      iii1=ii1+i1
      do i2=-1,1
        iii2=ii2+i2
        do i3=-1,1
          iii3=ii3+i3
          ie=boxes(iii3,iii2,iii1)
           do while(ie.ne.0)
             if (ie==i) then ; ie=list(ie) ; cycle ; end if
             if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
             if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15 THAT LINE IS DIFFERENT WITH THE _inter_other_side_epi VERSION OF THE SUBROUTINE
             dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
             a=dai+node(ie)%add
             if(dist>a)then ; ie=list(ie) ; cycle ; end if             
                  
             !>>Miquel31-12-14
             snneigh=snneigh+1
             sneigh(snneigh)=ie
             sdneigh(snneigh)=dist
             ie=list(ie)
           end do
         end do
       end do
     end do

     ! here we transfer the neighbors of i to a pre-neigh neighbor
     if (snneigh>omnn)then                                !!>> HC 15-1-2021 the node i has increased the maximum number of neighbors
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We save the original data (we need to change the dimension of neigh and dneigh)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        omnn=snneigh                                      !!>> HC 15-1-2021
        if(allocated(neigh)) deallocate(neigh)            !!>> HC 15-1-2021 We realocate the matrixes neigh and dneigh to change their dimension
        if(allocated(dneigh)) deallocate(dneigh)          !!>> HC 15-1-2021
        allocate(neigh(nda,omnn),dneigh(nda,omnn))        !!>> HC 15-1-2021
        neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)       !!>> HC 15-1-2021 We recover the original data
        dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)     !!>> HC 15-1-2021
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 Here we add the new neighbors of the node i
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     else                                                 !!>> HC 15-1-2021 In this case the node i has the same or less neighbors than the rest of nodes
        nneigh(i)=snneigh                                 !!>> HC 15-1-2021 the new number of neighbors is registered in nneigh
        neigh(i,1:omnn)=sneigh(:omnn)                     !!>> HC 15-1-2021 we add the new neighs to the matrixes neigh and dneigh directly
        dneigh(i,1:omnn)=sdneigh(:omnn)                   !!>> HC 15-1-2021
     endif                                                !!>> HC 15-1-2021

end subroutine

!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************
!************************************************************************************************************************************************************

subroutine neighbor_build_node_old(i) !this calculates the neighbors for one node onlly, used for random noise !>>Miquel16-12-14
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh
integer::i

!!!triangulation variables
integer:: npt !number of points  
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
real*8::dai
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer::snneigh


  maxlen=sqrt((2*maxval(node(:nd)%eqs))**2+(2*maxval(node(:nd)%add)**2)) !this is the maximal interaction distance between epithelial nodes
  a=2*maxval(node(:nd)%add) !maximal interaction distance between mesenchymal nodes
  if(a>maxlen)then ; rv=a ; urv=1.0d0/a ; else ; rv=maxlen ; urv=1.0d0/maxlen ;end if


!  rv=2*maxval(node(:nd)%add);urv=1d0/rv !;print*,"RV",rv
!  a=maxval(node(:nd)%eqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%add)*dmax
  
  call iniboxes
 
    if (rdiffmax<2*node(i)%add)then
      nbh=2*node(i)%add  !the neighbor search range for diffusion is the same as for node interactions
    else
      nbh=rdiffmax
    end if
    ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
    ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
    ivv=node(i)%altre
    
    nbo=nint(nbh*urv) !;print*,"nbo",nbo,"rdiffmax",rdiffmax
    mnn_dynam=mnn_dyn*(2*nbo+1)**3 !;print*,"mnn_dyn def",mnn_dynam !calculating the alleged maximal width of the neigh matrix
    
    do j=1,nneigh(i)                            !!This to erase node i from the neighbor matrix, as it will be ovewritten later  !! >>>Miguel17-12-14
      k=neigh(i,j)                              ! "k" neighbor
      do ii=1,nneigh(k)                         ! it searches in the neighborhood of "k"         
        if(i.eq.neigh(k,ii))then                ! if it (i) appears    
          neigh(k,ii:nneigh(k)-1)=neigh(k,ii+1:nneigh(k))   ; neigh(k,nneigh(k))=0   !  neigh matrix is displaced
          dneigh(k,ii:nneigh(k)-1)=dneigh(k,ii+1:nneigh(k)) ; dneigh(k,nneigh(k))=0  !  dneigh matrix is diaplced
          nneigh(k)=nneigh(k)-1 ; exit                                               !  neighbor counter
        end if 
      end do        
    end do
    neigh(i,:)=0 ; nneigh(i)=0 ; dneigh(i,:)=0d0  !! >>>Miguel17-12-14
    
    allocate(sdneigh(mnn_dynam),sneigh(mnn_dynam))
    snneigh=0
    sneigh=0 ; sdneigh=0d0
    
    
    tipi=node(i)%tipus
    if(node(i)%tipus<3)then
      snneigh=1
      sneigh(1)=ivv 
      sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
    end if

    !>>> Is 18-4-15
    if (ffu(17)==0) then
      dai=node(i)%add
      do i1=-nbo,nbo 
        iii1=ii1+i1
        do i2=-nbo,nbo
          iii2=ii2+i2
          do i3=-nbo,nbo
            iii3=ii3+i3
            ie=boxes(iii3,iii2,iii1)
            do while(ie.ne.0)
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
              if (ie==ivv) then ; ie=list(ie) ; cycle ; end if   
              if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
              if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
              dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              !>>Miquel31-12-14
              a=dai+node(ie)%add
              if(tipi>=3)then
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>a)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  b=sqrt(2*(a**2))
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              else !epithelial vs epithelial
                b=sqrt(2*(a**2))
                if(tipi==node(ie)%tipus)then !same face epithelials
                  if(b<maxlen) b=maxlen
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                else
                  ie=list(ie) ; cycle !>>> Is 18-4-15
                end if
              end if
              !>>Miquel31-12-14
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              !dneigh(i,nneigh(i))=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              sdneigh(snneigh)=dist
              ie=list(ie)
            end do
          end do
        end do
      end do
      !<<< Is 18-4-15
    else
      dai=node(i)%add
      do i1=-nbo,nbo 
        iii1=ii1+i1
        do i2=-nbo,nbo
          iii2=ii2+i2
          do i3=-nbo,nbo
            iii3=ii3+i3
            ie=boxes(iii3,iii2,iii1)
            do while(ie.ne.0)
              if (ie==i) then ; ie=list(ie) ; cycle ; end if
              if (ie==ivv) then ; ie=list(ie) ; cycle ; end if   
              dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              !>>Miquel31-12-14
              a=dai+node(ie)%add
              if(tipi>=3)then
                if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                  if(dist>a)then ; ie=list(ie) ; cycle ; end if
                else  !mesench/ECM vs epithelial
                  b=sqrt(2*(a**2))
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              else !epithelial vs epithelial
                b=sqrt(2*(a**2))
                if(tipi==node(ie)%tipus)then !same face epithelials
                  if(b<maxlen) b=maxlen
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                else
                  if(dist>b)then ; ie=list(ie) ; cycle ; end if
                end if
              end if
              !>>Miquel31-12-14
              snneigh=snneigh+1
              sneigh(snneigh)=ie
              !dneigh(i,nneigh(i))=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
              sdneigh(snneigh)=dist
              ie=list(ie)
            end do
          end do
        end do
      end do
    end if

    allocate(osdneigh(snneigh),osneigh(snneigh))
    osdneigh(1:snneigh)=sdneigh(1:snneigh)
    osneigh(1:snneigh)=sneigh(1:snneigh)
    
    if((ffu(2)==1).or.(ffu(8)==1))then !! >>>Miguel17-12-14
      !screening by Gabriel graph !>>Miquel6-3-14
      !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
      do j=1,snneigh-1  !ordering the neighors with increasing distance
        b=osdneigh(j)
        ii=0
        do k=j+1,snneigh
          c=osdneigh(k)
           if(b>c)then
            ii=k ; b=osdneigh(k)
          end if
        end do
        if(ii/=0)then
          !jj=neigh(i,j)
          kk=osneigh(ii) ; c=osdneigh(ii) !the swap
          osneigh(ii)=osneigh(j) ; osdneigh(ii)=osdneigh(j)
          osneigh(j)=kk ; osdneigh(j)=c
        end if
      end do
      if(ffu(8).eq.1)then;b=0.75d0;else;b=screen_radius;endif ! >>>Miguel17-12-14
      !the screening
      ii=0 !the number of eliminated connections
      do j=snneigh,1,-1
        jj=osneigh(j)
        if(jj==ivv) cycle
        a=osdneigh(j)*0.5*b !the radius of the sphere !>>Miquel28-7-14   ! >>>Miguel17-12-14       
        jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
        cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
        do k=j-1,1,-1
          kk=osneigh(k)
          d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
          if(d-a<epsilod)then !there is one node within the sphere, we must delete this connection
            do l=j,snneigh-ii-1
              osneigh(l)=osneigh(l+1)
              osdneigh(l)=osdneigh(l+1)
            end do
            osneigh(snneigh-ii)=0
            osdneigh(snneigh-ii)=0
            ii=ii+1
            exit
          end if
        end do
      end do
      snneigh=snneigh-ii       
   ! else !>> HC 19-5-20202 trans_neigh has not been allocated yet this crashes
    !  trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
     ! trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15 
    end if
    
    if(snneigh>omnn)then !the whole matrix has to be reallocated
      ii=omnn
      allocate(trans_neigh(nd,omnn),trans_dneigh(nd,omnn))
      trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)
      trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)
      deallocate(neigh,dneigh)
      omnn=snneigh
      allocate(neigh(nda,omnn),dneigh(nda,omnn))
      neigh(1:nd,1:ii)=trans_neigh(1:nd,1:ii)
      dneigh(1:nd,1:ii)=trans_dneigh(1:nd,1:ii)
      deallocate(trans_neigh,trans_dneigh)
    end if
    
    nneigh(i)=snneigh
    neigh(i,1:snneigh)=osneigh(1:snneigh)
    dneigh(i,1:snneigh)=osdneigh(1:snneigh)
    
    do j=1,snneigh !here we put i on the reciprocal neighborhood of its neighbors
      k=osneigh(j)
      nneigh(k)=nneigh(k)+1
      if(nneigh(k)>omnn)then  !the whole matrix has to be reallocated
        ii=omnn
        allocate(trans_neigh(nd,omnn),trans_dneigh(nd,omnn))
        trans_neigh(1:nd,1:omnn)=neigh(1:nd,1:omnn)
        trans_dneigh(1:nd,1:omnn)=dneigh(1:nd,1:omnn)
        deallocate(neigh,dneigh)
        omnn=nneigh(k)
        allocate(neigh(nda,omnn),dneigh(nda,omnn))
        neigh(1:nd,1:ii)=trans_neigh(1:nd,1:ii)
        dneigh(1:nd,1:ii)=trans_dneigh(1:nd,1:ii)
        deallocate(trans_neigh,trans_dneigh)
      end if
      neigh(k,nneigh(k))=i
      dneigh(k,nneigh(k))=osdneigh(j)
    end do

if(allocated(vcl)) deallocate(vcl)
if(allocated(vm)) deallocate(vm)
if(allocated(ht)) deallocate(ht)
if(allocated(bf)) deallocate(bf)
if(allocated(fc)) deallocate(fc)
if(allocated(osneigh)) deallocate(osneigh)
if(allocated(osdneigh)) deallocate(osdneigh)
if(allocated(sneigh)) deallocate(sneigh)
if(allocated(sdneigh)) deallocate(sdneigh)
if(allocated(trans_neigh)) deallocate(trans_neigh)
if(allocated(trans_dneigh)) deallocate(trans_dneigh)


end subroutine

!************************************************************************************************************************************************
!************************************************************************************************************************************************
!************************************************************************************************************************************************


subroutine neighbor_build_old
integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh

integer::tipi  !>>Miquel31-12-14
real*8::dai,maxlen    !>>Miquel31-12-14

!!!triangulation variables
integer:: npt !number of points
integer:: sizht !size of hash table
integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
real*8,allocatable :: vcl(:,:) !point coordinates
integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14
integer,allocatable::osneigh(:),sneigh(:),trans_neigh(:,:)
real*8,allocatable::osdneigh(:),sdneigh(:),trans_dneigh(:,:)
integer,allocatable::oosneigh(:)
real*8,allocatable::oosdneigh(:)
integer::snneigh


!  a=maxval(node(:nd)%eqs)
!  if (rv<a) then ; rv=a;urv=1d0/rv ;end if !this is mostly to allow diffusion between the two faces of the epithelium
  rdiffmax=2*maxval(node(:nd)%add)*dmax

  !omnn=0
  call iniboxes
  mnn_dynam=mnn_dyn*(2*nint(rdiffmax*urv)+1)**3
  allocate(trans_neigh(nd,mnn_dynam),trans_dneigh(nd,mnn_dynam))

  if (ffu(8)==0)then !normal neighboring, extensive search of the boxes
    omnn=0
    neigh=0  ! 4-3-2020
    nneigh=0 ! 4-3-2020
    do i=1,nd
      if (rdiffmax<2*node(i)%add)then
        nbh=2*node(i)%add  !the neighbor search range for diffusion is the same as for node interactions
      else
        nbh=rdiffmax
      end if
      
      ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
      ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
      ivv=node(i)%altre

      nbo=nint(nbh*urv) !;print*,"nbo",nbo,"rdiffmax",rdiffmax
      mnn_dynam=mnn_dyn*(2*nbo+1)**3 !;print*,"mnn_dyn def",mnn_dynam !calculating the alleged maximal width of the neigh matrix

      allocate(sdneigh(mnn_dynam),sneigh(mnn_dynam))
      snneigh=0
      sneigh=0 ; sdneigh=0d0

      if(node(i)%tipus<3)then
        snneigh=1
        sneigh(1)=ivv
        sdneigh(1)=sqrt((node(ivv)%x-ix)**2+(node(ivv)%y-iy)**2+(node(ivv)%z-iz)**2)
      end if


      !>>> Is 18-4-15
      if (ffu(23)==1) then ! IS 22-12-20 WHETHER WE HAVE SIMPLIFIED NEIGHBORING !!>> HC 14-1-2021
      
        ! HERE WE DON'T
        if (ffu(17)==0) then  !if 0 epithelial nodes from one side do not consider as neighbors nodes from the other side
            tipi=node(i)%tipus
            dai=node(i)%add
            do i1=-nbo,nbo 
              iii1=ii1+i1
              do i2=-nbo,nbo
                iii2=ii2+i2
                do i3=-nbo,nbo
                  iii3=ii3+i3
                  ie=boxes(iii3,iii2,iii1)
                  do while(ie.ne.0)
                    if (ie==i) then ; ie=list(ie) ; cycle ; end if
                    if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                      if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                      if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                    dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
                    !>>Miquel31-12-14
                    a=dai+node(ie)%add
                    if(tipi>=3)then
                      if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                       if(dist>a)then ; ie=list(ie) ; cycle ; end if
                      else  !mesench/ECM vs epithelial
                        b=arrel_de_dos*a
                        if(dist>b)then ; ie=list(ie) ; cycle ; end if
                      end if
                    else !epithelial vs epithelial
                      b=arrel_de_dos*a
                      if(tipi==node(ie)%tipus)then !same face epithelials
                        if(b<maxlen) b=maxlen
                        if(dist>b)then ; ie=list(ie) ; cycle ; end if
                      else
                        if(dist>b)then ; ie=list(ie) ; cycle ; end if    
                      end if
                    end if
                  
                    !>>Miquel31-12-14
                    snneigh=snneigh+1
                    sneigh(snneigh)=ie
                    sdneigh(snneigh)=dist
                    ie=list(ie)
                  end do
                end do
              end do
            end do
       
          
        else  ! epithelial nodes from one side do consider as neighbors nodes from the other side
          tipi=node(i)%tipus
          dai=node(i)%add
          do i1=-nbo,nbo 
            iii1=ii1+i1
            do i2=-nbo,nbo
              iii2=ii2+i2
              do i3=-nbo,nbo
                iii3=ii3+i3
                ie=boxes(iii3,iii2,iii1)
                do while(ie.ne.0)
                  if (ie==i) then ; ie=list(ie) ; cycle ; end if
                  if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                  dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
                  !>>Miquel31-12-14
                  a=dai+node(ie)%add
                  if(tipi>=3)then
                    if(node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
                      if(dist>a)then ; ie=list(ie) ; cycle ; end if
                    else  !mesench/ECM vs epithelial
                      b=arrel_de_dos*a !sqrt(2*(a**2))
                      if(dist>b)then ; ie=list(ie) ; cycle ; end if
                    end if
                  else !epithelial vs epithelial
                    b=arrel_de_dos*a !sqrt(2*(a**2))
                    if(tipi==node(ie)%tipus)then !same face epithelials
                      if(b<maxlen) b=maxlen
                      if(dist>b)then ; ie=list(ie) ; cycle ; end if
                    else
                      if(dist>b)then ; ie=list(ie) ; cycle ; end if    
                    end if
                  end if
                  !>>Miquel31-12-14
                  snneigh=snneigh+1
                  sneigh(snneigh)=ie
                  sdneigh(snneigh)=dist
                  ie=list(ie)
                end do
              end do
            end do
          end do
        end if
        
      else          ! 

      ! IS 22-12-20 WITH SIMPLIFIED NEIGHBORING
        if (ffu(17)==0) then  !if 0 epithelial nodes from one side do not consider as neighbors nodes from the other side
            tipi=node(i)%tipus
            dai=node(i)%add
            do i1=-nbo,nbo 
              iii1=ii1+i1
              do i2=-nbo,nbo
                iii2=ii2+i2
                do i3=-nbo,nbo
                  iii3=ii3+i3
                  ie=boxes(iii3,iii2,iii1)
                  do while(ie.ne.0)
                    if (ie==i) then ; ie=list(ie) ; cycle ; end if
                    if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                      if (node(ie)%tipus==1.and.tipi==2) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                      if (node(ie)%tipus==2.and.tipi==1) then ; ie=list(ie) ; cycle ; end if !>>> Is 18-4-15
                    dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)                    
                    a=dai+node(ie)%add
                    if(dist>a)then ; ie=list(ie) ; cycle ; end if
                  
                    !>>Miquel31-12-14
                    snneigh=snneigh+1
                    sneigh(snneigh)=ie
                    sdneigh(snneigh)=dist
                    ie=list(ie)
                  end do
                end do
              end do
            end do
          
        else  ! epithelial nodes from one side do consider as neighbors nodes from the other side
        
          tipi=node(i)%tipus
          dai=node(i)%add
          do i1=-nbo,nbo 
            iii1=ii1+i1
            do i2=-nbo,nbo
              iii2=ii2+i2
              do i3=-nbo,nbo
                iii3=ii3+i3
                ie=boxes(iii3,iii2,iii1)
                do while(ie.ne.0)
                  if (ie==i) then ; ie=list(ie) ; cycle ; end if
                  if (ie==ivv) then ; ie=list(ie) ; cycle ; end if
                  dist=sqrt((node(ie)%x-ix)**2+(node(ie)%y-iy)**2+(node(ie)%z-iz)**2)
                  a=dai+node(ie)%add
                  if(dist>a)then ; ie=list(ie) ; cycle ; end if
                  
                  !>>Miquel31-12-14
                  snneigh=snneigh+1
                  sneigh(snneigh)=ie
                  sdneigh(snneigh)=dist
                  ie=list(ie)
                end do
              end do
            end do
          end do
        end if
      
      end if

      allocate(osdneigh(snneigh),osneigh(snneigh))
      osdneigh(1:snneigh)=sdneigh(1:snneigh)
      osneigh(1:snneigh)=sneigh(1:snneigh)

   !print*,"sneigh",sneigh(1:nneigh(i))
      if(ffu(2)==1)then
        !screening by Gabriel graph !>>Miquel6-3-14
        !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
        !ordering the neighors with increasing distance

        allocate(oosdneigh(snneigh),oosneigh(snneigh))

        !sorting algorithm by selection, it's ok
        do j=1,snneigh-1  
          b=osdneigh(j)
          ii=0
          do k=j+1,snneigh
            c=osdneigh(k)
             if(b>c)then
              ii=k ; b=osdneigh(k)
            end if
          end do
          if(ii/=0)then
            !jj=osneigh(j)
            kk=osneigh(ii) ; c=osdneigh(ii) !the swap
            osneigh(ii)=osneigh(j) ; osdneigh(ii)=osdneigh(j)
            osneigh(j)=kk ; osdneigh(j)=c
          end if
        end do

        !the screening
        !ii=0 !the number of eliminated connections
        do j=snneigh,1,-1
          jj=osneigh(j)
          if(jj==ivv) cycle
          a=osdneigh(j)*0.5*screen_radius !the radius of the sphere !>>Miquel28-7-14
          jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
          cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
          do k=j-1,1,-1
            kk=osneigh(k)
            d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
            if(d-a<epsilod)then !there is one node within the sphere, we must delete this connection
              osneigh(j)=0
              exit
            end if
          end do
        end do
        ii=0
        do j=1,snneigh
          jj=osneigh(j)
          if(jj/=0)then
            ii=ii+1
            oosneigh(ii)=jj ; oosdneigh(ii)=osdneigh(j)
          end if
        end do
        snneigh=ii
        !snneigh=snneigh-ii
        !if(i==1) print*,nneigh(1),"neigh1",neigh(1,1:nneigh(1))

        trans_neigh(i,1:snneigh)=oosneigh(1:snneigh)
        trans_dneigh(i,1:snneigh)=oosdneigh(1:snneigh)
      else                                                 !>>> Is 23-4-15
        trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
        trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15
      end if                                            
      
      if(snneigh>omnn) omnn=snneigh
      nneigh(i)=snneigh
      deallocate(osdneigh,osneigh,sdneigh,sneigh)
      if (ffu(2)==1) deallocate(oosneigh,oosdneigh)

    end do
    !omnn=maxval(nneigh,dim=1,mask=nneigh<=nd) ; print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
   
  else   !*****3D triangulation neighooring***************************************************************
    
    npt=nd
!    sizht=3/2*npt
    sizht=3*npt
    maxfc=npt**2
    maxbf=npt**2
    allocate(vcl(3,npt),vm(npt))
    allocate(bf(1:3,maxbf),fc(1:7,maxfc))
    allocate(ht(0:sizht-1))
    
    
    !call extrem !not necessary, but sets some variables that later may be used in pinta and creix !>>Miquel21-3-14
    
    
    do i=1,npt  !building the input arrays for the triangualtion function
      vcl(1,i)=node(i)%x ; vcl(2,i)=node(i)%y ; vcl(3,i)=node(i)%z
      vm(i)=i
    end do !;print*,"pillat pre delau"
    !calling the triangulation subroutine
    call dtriw3(npt, sizht, maxbf, maxfc, vcl, vm, nbf, nfc, nface, ntetra, bf, fc, ht, ierr)
    !if(ierr/=0)then; print*,"error in the triangulation, avort or something",ierr,getot ; endif !call exit(24); end if
    !translating the ouptut into a neighbor matrix

    do i=1,nd !initialize neighbor matrix
!      neigh(i,1:nneigh(i))=0
      neigh(i,1:)=0  !IS 29-4-14
      nneigh(i)=0
    end do
    do j=1,nfc               ! it passes through all triangles 
      if((fc(1,j)>0).and.(fc(2,j)>0).and.(fc(3,j)>0))then ! it is a "valid" triangle
        ii=fc(1,j) ; jj=fc(2,j) ; kk=fc(3,j)
        iii=vm(ii);jjj=vm(jj);kkk=vm(kk) !; print*,"iii jjj kkk",iii,jjj,kkk
        tiii=node(iii)%tipus ; tjjj=node(jjj)%tipus ; tkkk=node(kkk)%tipus

        ivv=nneigh(iii)
        ij=0 ; ik=0 ; jk=0
        !connection iii-jjj*****
        if(ivv==0)then
          d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
          ivv=1
          trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
        else
          do i=1,nneigh(iii)
            if(trans_neigh(iii,i)==jjj)then; ij=1;exit;end if
          end do
          if(ij==0)then
            d=sqrt((vcl(1,iii)-vcl(1,jjj))**2+(vcl(2,iii)-vcl(2,jjj))**2+(vcl(3,iii)-vcl(3,jjj))**2)
            ivv=ivv+1
            trans_neigh(iii,ivv)=jjj ; trans_dneigh(iii,ivv)=d
          end if
        end if
        if(ij==0)then
          !connection jjj-iii*******
          if(nneigh(jjj)==0)then ; nneigh(jjj)=1 ; trans_neigh(jjj,1)=iii ; trans_dneigh(jjj,1)=d
          else; iv=nneigh(jjj)+1 ; trans_neigh(jjj,iv)=iii ; trans_dneigh(jjj,iv)=d ; nneigh(jjj)=iv ; end if
        end if
        !connection iii-kkk*****
        do i=1,nneigh(iii)
          if(trans_neigh(iii,i)==kkk)then; ik=1;exit;end if
        end do
        if(ik==0)then
          d=sqrt((vcl(1,iii)-vcl(1,kkk))**2+(vcl(2,iii)-vcl(2,kkk))**2+(vcl(3,iii)-vcl(3,kkk))**2)
          ivv=ivv+1
          trans_neigh(iii,ivv)=kkk ; trans_dneigh(iii,ivv)=d
          !connection kkk-iii*******
          if(nneigh(kkk)==0)then ; nneigh(kkk)=1 ; trans_neigh(kkk,1)=iii ; trans_dneigh(kkk,1)=d
          else; iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=iii ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv ; end if
        end if
        nneigh(iii)=ivv
        !connection jjj-kkk*****
        do i=1,nneigh(jjj)
          if(trans_neigh(jjj,i)==kkk)then; jk=1;exit;end if
        end do
        if(jk==0)then
          d=sqrt((vcl(1,jjj)-vcl(1,kkk))**2+(vcl(2,jjj)-vcl(2,kkk))**2+(vcl(3,jjj)-vcl(3,kkk))**2)
          iv=nneigh(jjj)+1
          trans_neigh(jjj,iv)=kkk ; trans_dneigh(jjj,iv)=d ;nneigh(jjj)=iv
          !connection kkk-jjj*******
          iv=nneigh(kkk)+1 ; trans_neigh(kkk,iv)=jjj ; trans_dneigh(kkk,iv)=d ; nneigh(kkk)=iv
        end if
      end if
    end do


    if(ffu(2)==1)then
    
      do i=1,nd
        !Screening by Gabriel graph !>>Miquel6-3-14
        !A neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
        !Ordering the neighors with increasing distance

        if (ffu(9)==0) then   !!!! Is 4-3-2020 APARENTLY THAT WAS A BUT IN THE 2016 VERSION, BUT WE KEEP FOR COMPATIBILITY WITH THE EXAMPLES IN THE OLD VERSION
          ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
        end if

        !sorting algorithm by selection, it's ok
        do j=1,nneigh(i)-1  
          b=trans_dneigh(i,j)
          ii=0
          do k=j+1,nneigh(i)
            c=trans_dneigh(i,k)
             if(b>c)then
              ii=k ; b=trans_dneigh(i,k)
            end if
          end do
          if(ii/=0)then
            !jj=osneigh(j)
            kk=trans_neigh(i,ii) ; c=trans_dneigh(i,ii) !the swap
            trans_neigh(i,ii)=trans_neigh(i,j) ; trans_dneigh(i,ii)=trans_dneigh(i,j)
            trans_neigh(i,j)=kk ; trans_dneigh(i,j)=c
          end if
        end do
        
        !the screening
        ii=0 !the number of eliminated connections
        do j=nneigh(i),1,-1
          jj=trans_neigh(i,j)
          if(jj==ivv) cycle
          a=trans_dneigh(i,j)*0.5*screen_radius !the radius of the sphere
          jx=node(jj)%x ; jy=node(jj)%y ; jz=node(jj)%z
          cx=(ix+jx)*0.5 ; cy=(iy+jy)*0.5 ; cz=(iz+jz)*0.5 !the midpoint
          do k=j-1,1,-1
            kk=trans_neigh(i,k)
            d=sqrt((cx-node(kk)%x)**2+(cy-node(kk)%y)**2+(cz-node(kk)%z)**2)
            if(d<a)then !there is one node within the sphere, we must delete this connection
              do l=j,nneigh(i)-ii-1
                trans_neigh(i,l)=trans_neigh(i,l+1)
                trans_dneigh(i,l)=trans_dneigh(i,l+1)
              end do
              trans_neigh(i,nneigh(i)-ii)=0
              trans_dneigh(i,nneigh(i)-ii)=0
              ii=ii+1
              exit
            end if
          end do
        end do
        nneigh(i)=nneigh(i)-ii
      end do
    !else                                                 !>>> Is 23-4-15
    !  trans_neigh(i,1:snneigh)=osneigh(1:snneigh)        !>>> Is 23-4-15
    !  trans_dneigh(i,1:snneigh)=osdneigh(1:snneigh)      !>>> Is 23-4-15      
    end if
    
    omnn=0
    do i=1,nd
      if(nneigh(i)>omnn) omnn=nneigh(i)
    end do
    !print*,"omnn",omnn
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)
    allocate(neigh(nda,omnn),dneigh(nda,omnn))
    neigh(1:nd,1:omnn)=trans_neigh(1:nd,1:omnn)
    dneigh(1:nd,1:omnn)=trans_dneigh(1:nd,1:omnn)
  end if
!!do i=1,nd
!end do  



if (ffu(9)==0) then
  if(allocated(trans_neigh)) deallocate(trans_neigh)
  if(allocated(trans_dneigh)) deallocate(trans_dneigh)
  if(allocated(sdneigh)) deallocate(sdneigh)
  if(allocated(osdneigh)) deallocate(osdneigh)
  if(allocated(osneigh)) deallocate(osneigh)
  if(allocated(sneigh)) deallocate(sneigh)
  if(allocated(oosdneigh)) deallocate(oosdneigh)
  if(allocated(oosneigh)) deallocate(oosneigh)
  if(allocated(vcl)) deallocate(vcl)
  if(allocated(vm)) deallocate(vm)
  if(allocated(bf)) deallocate(bf)
  if(allocated(fc)) deallocate(fc)
  if(allocated(ht)) deallocate(ht)
end if
end subroutine neighbor_build_old

!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************RECOVERING NEIGHBORS ALGORITHM******************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!

subroutine fill_co_griders                                                    !!>> HC 10-6-2021 This subroutine fills the co_griders matrix with the neigh matrix
implicit none                                                                 !!>> HC 10-6-2021 
integer :: ich, jch, nmaxneigh, tipch, neich, ntipch, ord                     !!>> HC 10-6-2021
real*8 :: deich,ichx,ichy,ichz, jchx, jchz, jchy                              !!>> HC 10-6-2021
ndch=nd                                                                       !!>> HC 15-6-2021 This is a public variable storing the number of nodes registered in the cogrid_matrix

if (allocated(oneigh)) deallocate(oneigh)                                     !!>> HC 29-6-2021
allocate(oneigh(1:nd,1:maxval(nneigh(1:nd))))                                 !!>> HC 29-6-2021
if (allocated(onneigh)) deallocate(onneigh)                                   !!>> HC 29-6-2021
allocate(onneigh(1:nd))                                                       !!>> HC 29-6-2021
oneigh=0; onneigh=0                                                           !!>> HC 29-6-2021

do ich=1,nd                                                                   !!>> HC 10-6-2021 search all the epithelial nodes
   if (node(ich)%tipus.ge.3) cycle                                            !!>> HC 10-6-2021
   do jch=1,nneigh(ich)                                                       !!>> HC 10-6-2021 go over their neigh matrix
      neich=neigh(ich,jch)                                                    !!>> HC 10-6-2021
      if(neich<ich)cycle                                                      !!>> HC 30-6-2021
      if (node(neich)%tipus.ge.3)cycle                                        !!>> HC 10-6-2021 find epithelial neighs
      onneigh(ich)=onneigh(ich)+1                                             !!>> HC 29-6-2021
      oneigh(ich,onneigh(ich))=neich                                          !!>> HC 29-6-2021
      onneigh(neich)=onneigh(neich)+1                                         !!>> HC 30-6-2021
      oneigh(neich,onneigh(neich))=ich                                        !!>> HC 30-6-2021
   end do                                                                     !!>> HC 10-6-2021
end do                                                                        !!>> HC 10-6-2021

end subroutine                                                                !!>> HC 10-6-2021


!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!
!************************************************************************************************************************************************************************!



subroutine restore_neighbors                                            !!>> HC 15-6-2021 This subroutine compares neigh and oneigh and restores the links that have been lost
implicit none                                                           !!>> HC 15-6-2021
integer :: ich, jch, newich, newjch, oldch                              !!>> HC 15-6-2021
integer :: neich, leich, peich, lch, kch, pch, recover                  !!>> HC 15-6-2021
integer, allocatable, dimension(:,:) :: trans_neich                     !!>> HC 15-6-2021
real, allocatable, dimension(:,:) :: trans_deich                        !!>> HC 15-6-2021
real :: deich, ichx, ichz, ichy, jchx, jchz, jchy, ichadd,jchadd        !!>> HC 15-6-2021
integer :: addneigh,isaneigh                                            !!>> HC 15-6-2021
integer:: mch,meich,gch,geich,ijch                                      !!>> HC 15-6-2021


do ich=1,ndch                                                           !!>> HC 15-6-2021 This goes through ndch nodes that were there before any division happened 
   if (node(ich)%tipus>2)cycle                                          !!>> HC 15-6-2021 only epithelial cells
   ichadd=node(ich)%add                                                 !!>> HC 28-6-2021
   do ijch=1,onneigh(ich)                                               !!>> HC 15-6-2021 Loop though all the interactions between nodes 
      jch=oneigh(ich,ijch)                                              !!>> HC 29-6-2021
      if (jch<ich)cycle                                                 !!>> HC 29-6-2021 This neighbor has already been checked
      if (jch>ndch)cycle                                                !!>> HC 29-6-2021 This neighbor has just appeared by division 
      if (node(jch)%tipus>2)cycle                                       !!>> HC 15-6-2021 only epithelial nodes
      jchadd=node(ich)%add                                              !!>> HC 28-6-2021
      addneigh=0                                                        !!>> HC 17-6-2021
      do kch=1, nneigh(ich)                                             !!>> HC 17-6-2021 we want to know if ich and jch are still connected in the add range
         neich=neigh(ich,kch)                                           !!>> HC 17-6-2021 look over the neighs of ich
         if (neich==jch)then                                            !!>> HC 17-6-2021 Here we find that jch is still in neigh, BUT it still can be a neighbor outside the ADD range
            addneigh=1                                               !!>> HC 17-6-2021 if addneigh is 1 it means that jch is still and ADD neighbor of ich and we do not need to restore it
            exit
         endif                                                          !!>> HC 17-6-2021 if it is 0 it means that either it is not a neigh at all or it is a neigh outside the ADD range
      enddo                                                             !!>> HC 17-6-2021
      if (addneigh==1) cycle                                            !!>> HC 15-6-2021 only nodes that are not connected in the ADD range anymore       
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 15-6-2021 SHOULD I STAY OR SHOULD I GO ALGORITHM (determining if we have to recover or not this interaction)
      recover=1                                                            !!>> HC 15-6-2021 Let us assume that the interaction is going to be recovered
      do kch=1,nneigh(ich)                                                 !!>> HC 29-6-2021 look the neighbors of node ich
         neich=neigh(ich,kch)                                              !!>> HC 29-6-2021
!         if (dneigh(ich,kch)-node(neich)%add-ichadd>epsilod)cycle          !!>> HC 29-6-2021 neich has to be an ADD neighbor of ich
         do lch=1,nneigh(jch)                                              !!>> HC 29-6-2021 look the neighbors of node jch
            leich=neigh(jch,lch)                                           !!>> HC 29-6-2021
            if (leich.ne.neich)cycle                                       !!>> HC 29-6-2021 
!            if (dneigh(jch,lch)-node(neich)%add-jchadd>epsilod)exit        !!>> HC 29-6-2021 neich has to be an ADD neighbor of jch too
            do pch=1,nneigh(ich)                                           !!>> HC 29-6-2021 look again all the neighbors of ich to obtain the second node of the pair (peich)
               peich=neigh(ich,pch)                                        !!>> HC 29-6-2021
               if (peich==neich)cycle                                      !!>> HC 29-6-2021 peich and neich cannot be the same node
!               if (dneigh(ich,pch)-node(peich)%add-ichadd>epsilod)cycle    !!>> HC 29-6-2021 peich has to be an ADD neighbor of ich
               do mch=1,nneigh(jch)                                        !!>> HC 29-6-2021 look the neighbors of node jch
                  meich=neigh(jch,mch)                                     !!>> HC 29-6-2021
                  if (meich.ne.peich)cycle                                 !!>> HC 29-6-2021
!                  if (dneigh(jch,mch)-node(peich)%add-jchadd>epsilod)exit  !!>> HC 29-6-2021 peich has to be an ADD neighbor of jch too
                  do gch=1,nneigh(peich)                                   !!>> HC 29-6-2021 look the neighbors of peich to know whether peich and neich are neighbors
                     geich=neigh(peich,gch)                                !!>> HC 29-6-2021
                     if (geich.ne.neich)cycle                              !!>> HC 29-6-2021  peich and neich are neighbors
!                     if (dneigh(peich,gch)-node(peich)%add-node(geich)%add>epsilod)exit  !!>> HC 29-6-2021 peich and neich are neighbors in the ADD range
                     recover=0                                             !!>> HC 29-6-2021  WE DO NOT RECOVER THIS INTERACTION
!                     print*, "connection not to be recovered", ich, jch, "connected", peich, neich 
                     go to 667                                             !!>> HC 29-6-2021 This is a dirty trick that saves time. when we find the first couple of connected peich and neich
                  enddo                                                    !!>> HC 29-6-2021 we exit the loop because we already know that this connection is not going to be recovered
               enddo                                                       !!>> HC 29-6-2021
            enddo                                                          !!>> HC 29-6-2021
         enddo                                                             !!>> HC 29-6-2021
      enddo                                                                !!>> HC 29-6-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 15-6-2021
      
      
      
667   if (recover==0)cycle                                     !!>> HC 17-6-2021 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 15-6-2021 WE ARE GOING TO RECOVER THE LINK BETWEEN NODE ICH AND JCH      
      newich=nneigh(ich)+1                                                           !!>> HC 15-6-2021
      newjch=nneigh(jch)+1                                                           !!>> HC 15-6-2021
      oldch=size(neigh,2)                                                            !!>> HC 15-6-2021
      if (newich>oldch .or. newjch>oldch)then                                        !!>> HC 15-6-2021 We have to increase the size of the neigh matrix
         if (newjch>newich) newich=newjch                                            !!>> HC 15-6-2021 New size
         if (allocated(trans_neich)) deallocate(trans_neich)                         !!>> HC 10-6-2021 
         allocate(trans_neich(1:nd,1:oldch))                                         !!>> HC 15-6-2021
         trans_neich=0                                                               !!>> HC 15-6-2021
         trans_neich(1:nd,1:oldch)=neigh(1:nd,1:oldch)                               !!>> HC 15-6-2021 Save data
         if (allocated(neigh)) deallocate(neigh)                                     !!>> HC 10-6-2021 
         allocate(neigh(1:nda,1:newich))                                             !!>> HC 15-6-2021
         neigh=0                                                                     !!>> HC 15-6-2021
         neigh(1:nd,1:oldch)=trans_neich(1:nd,1:oldch)                               !!>> HC 15-6-2021 recover data
             
         if (allocated(trans_deich)) deallocate(trans_deich)                         !!>> HC 10-6-2021 we have to do the same with the dneigh matrix 
         allocate(trans_deich(1:nd,1:oldch))                                         !!>> HC 15-6-2021
         trans_deich=0.0d0                                                           !!>> HC 15-6-2021
         trans_deich(1:nd,1:oldch)=dneigh(1:nd,1:oldch)                              !!>> HC 15-6-2021 save data
         if (allocated(dneigh)) deallocate(dneigh)                                   !!>> HC 10-6-2021 
         allocate(dneigh(1:nda,1:newich))                                            !!>> HC 15-6-2021
         dneigh=0.0d0                                                                !!>> HC 15-6-2021
         dneigh(1:nd,1:oldch)=trans_deich(1:nd,1:oldch)                              !!>> HC 15-6-2021 recover data
         omnn=newich                                                                 !!>> HC 15-6-2021 
      endif                                                                          !!>> HC 15-6-2021 
      nneigh(ich)=nneigh(ich)+1                                                      !!>> HC 15-6-2021 Here we recover the conection
      nneigh(jch)=nneigh(jch)+1                                                      !!>> HC 15-6-2021 
      neigh(ich,nneigh(ich))=jch                                                     !!>> HC 15-6-2021 
      neigh(jch,nneigh(jch))=ich                                                     !!>> HC 15-6-2021 
        
      ichx=node(ich)%x; ichy=node(ich)%y; ichz=node(ich)%z                           !!>> HC 15-6-2021 
      jchx=node(jch)%x; jchy=node(jch)%y; jchz=node(jch)%z                           !!>> HC 15-6-2021 
         
      deich=sqrt((ichx-jchx)**2+(ichy-jchy)**2+(ichz-jchz)**2)                       !!>> HC 15-6-2021 We also have to recover the distance
      dneigh(ich,nneigh(ich))=deich                                                  !!>> HC 15-6-2021 
      dneigh(jch,nneigh(jch))=deich                                                  !!>> HC 15-6-2021 


   enddo                                                                             !!>> HC 15-6-2021 
enddo                                                                                !!>> HC 15-6-2021 

if(allocated(trans_deich)) deallocate(trans_deich)                                   !!>> HC 15-6-2021  trash allocated space
if(allocated(trans_neich)) deallocate(trans_neich)                                   !!>> HC 15-6-2021 

call fill_co_griders                                                                 !!>> HC 15-6-2021 fill the co_grides matrix for the next comparison

end subroutine

end module neighboring
