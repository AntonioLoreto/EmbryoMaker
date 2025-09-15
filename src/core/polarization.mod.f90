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

module polarization
use general
use neighboring
use genetic

real*8, public, allocatable  :: EQDe(:,:)
real*8, public, allocatable  :: ADDe(:,:)


contains

subroutine fill_pol_distances
  implicit none
  real*8 :: dotp, crossp, crosspx, crosspy, crosspz
  real*8 :: polx, poly, polz
  real*8 :: projx, projy, projz, toneighx, toneighy, toneighz
  real*8 :: dproj, dpol, dtoneigh
  real*8 :: higheqd, highadd, loweqd, lowadd, myfactor, diffeqd
  integer :: ich,jch,kch,lch, myneigh, dimneigh, kk
  
  dimneigh = maxval(nneigh(1:nd))+5
  
  if(allocated(EQDe))deallocate(EQDe)
  allocate(EQDe(1:nd,1:dimneigh)); EQDe=0d0
  if(allocated(ADDe))deallocate(ADDe)
  allocate(ADDe(1:nd,1:dimneigh)); ADDe=0d0
  
  do ich=1,nd
    if(npag(nparam_per_node+17) > 0) then
      ! Compute the factor by which ADD and EQD are modified
      myfactor=1.0d0
      do jch=1, npag(nparam_per_node+17)
        kk=whonpag(nparam_per_node+17,jch)
        if (gex(ich,kk)>0.0d0) then
          myfactor=myfactor+gex(ich,kk)*gen(kk)%e(nparam_per_node+17)
        end if
      end do
      if(myfactor > maxelong)myfactor=maxelong
      ! The polarization vector for the node is the one of its cell.
      if (node(ich)%tipus==4)then                                      !!>> HC 2-11-2021 ECM nodes do not elongate
         ADDe(ich,1:nneigh(ich))=node(ich)%add                         !!>> HC 2-11-2021 because they are not cells
         EQDe(ich,1:nneigh(ich))=node(ich)%eqd                         !!>> HC 2-11-2021
      else                                                             !!>> HC 2-11-2021
         polx=cels(node(ich)%icel)%polx
         poly=cels(node(ich)%icel)%poly
         polz=cels(node(ich)%icel)%polz
         dpol=sqrt(polx**2+poly**2+polz**2)
         if(dpol .eq. 0)then
           ADDe(ich,1:nneigh(ich))=node(ich)%add
           EQDe(ich,1:nneigh(ich))=node(ich)%eqd
         else
           ! High and low EQD and ADD values:
           diffeqd=node(ich)%add-node(ich)%eqd
           loweqd=node(ich)%eqd/myfactor; lowadd=loweqd + diffeqd
           higheqd=myfactor*node(ich)%eqd; highadd=higheqd + diffeqd
           do jch=1, nneigh(ich)
              myneigh=neigh(ich,jch)
              ! Compute dot product / cosinus
              toneighx=node(myneigh)%x-node(ich)%x
              toneighy=node(myneigh)%y-node(ich)%y
              toneighz=node(myneigh)%z-node(ich)%z
              dtoneigh=sqrt(toneighx**2+toneighy**2+toneighz**2)
              dotp = (polx*toneighx+poly*toneighy+polz*toneighz)/(dpol*dtoneigh)
              ! Compute cross product / sinus
              crossp=sqrt(1-dotp**2)
              ! ADDe and EQDe values depend on major/minor radii and cosinus/sinus
              ADDe(ich,jch)=sqrt((highadd*dotp)**2 + (lowadd*crossp)**2)
              EQDe(ich,jch)=sqrt((higheqd*dotp)**2 + (loweqd*crossp)**2)
           end do
         endif                                                        !!>> HC 2-11-2021
      end if
    else
      ADDe(ich,1:nneigh(ich))=node(ich)%add
      EQDe(ich,1:nneigh(ich))=node(ich)%eqd
    end if
  end do

end subroutine fill_pol_distances

!**************************************************************************

subroutine fill_pol_distances_node(i)
  implicit none
  real*8 :: dotp, crossp, crosspx, crosspy, crosspz
  real*8 :: polx, poly, polz
  real*8 :: projx, projy, projz, toneighx, toneighy, toneighz
  real*8 :: dproj, dpol, dtoneigh
  real*8 :: higheqd, highadd, loweqd, lowadd, diffeqd, myfactor
  integer :: ich,jch,kch,lch, myneigh, i, kk

  do ich=1,nd
    if(ich .ne. i .and. all(ich .ne. neigh(i,1:nneigh(i))))cycle
    if(npag(nparam_per_node+17) > 0) then
      ! Compute the factor by which ADD and EQD are modified
      myfactor=1.0d0
      do jch=1, npag(nparam_per_node+17)
        kk=whonpag(nparam_per_node+17,jch)
        if (gex(ich,kk)>0.0d0) then
          myfactor=myfactor+gex(ich,kk)*gen(kk)%e(nparam_per_node+17)
        end if
      end do
      if(myfactor > 2d0)myfactor=2d0
      ! The polarization vector for the node is the one of its cell.
      polx=cels(node(ich)%icel)%polx
      poly=cels(node(ich)%icel)%poly
      polz=cels(node(ich)%icel)%polz
      dpol=sqrt(polx**2+poly**2+polz**2)
      if(dpol .eq. 0)then
        ADDe(ich,1:nneigh(ich))=node(ich)%add
        EQDe(ich,1:nneigh(ich))=node(ich)%eqd
      else
        ! High and low EQD and ADD values:
        diffeqd=node(ich)%add-node(ich)%eqd
        loweqd=node(ich)%eqd/myfactor; lowadd=loweqd + diffeqd
        higheqd=myfactor*node(ich)%eqd; highadd=higheqd + diffeqd
        do jch=1, nneigh(ich)
          myneigh=neigh(ich,jch)
          ! Compute dot product / cosinus
          toneighx=node(myneigh)%x-node(ich)%x
          toneighy=node(myneigh)%y-node(ich)%y
          toneighz=node(myneigh)%z-node(ich)%z
          dtoneigh=sqrt(toneighx**2+toneighy**2+toneighz**2)
          dotp = (polx*toneighx+poly*toneighy+polz*toneighz)/(dpol*dtoneigh)
          ! Compute cross product / sinus
          crossp=sqrt(1-dotp**2)
          ! ADDe and EQDe values depend on major/minor radii and cosinus/sinus
          ADDe(ich,jch)=sqrt((highadd*dotp)**2 + (lowadd*crossp)**2)
          EQDe(ich,jch)=sqrt((higheqd*dotp)**2 + (loweqd*crossp)**2)
        end do
      end if
    else
      ADDe(ich,1:nneigh(ich))=node(ich)%add
      EQDe(ich,1:nneigh(ich))=node(ich)%eqd
    end if
  end do

end subroutine fill_pol_distances_node

!**************************************************************************
subroutine iniboxes_pola
integer ic,ii,jj,kk
integer,allocatable:: cboxes(:,:,:)
real*8 :: mymaxval
   
    call extrem
    mymaxval=0d0
    do i=1,nd
      if(maxval(ADDe(i,:)) > mymaxval)mymaxval=maxval(ADDe(i,:))
    end do
    a=2*mymaxval !maximal interaction distance between mesenchymal nodes
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
end subroutine iniboxes_pola

!**************************************************************************
subroutine iniboxes_p_pola
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
    
end subroutine iniboxes_p_pola

!**************************************************************************
subroutine iniboxesll_pola
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
end subroutine iniboxesll_pola

!************************************************


subroutine neighbor_build_pola !this subroutine is called at each iteration to assess the neighbors of all the nodes !>>>Miquel24-2-14
real*8 :: mymaxval

  call neighbor_build
  call fill_pol_distances
  !maxlen=sqrt(2*maxval(node(:nd)%eqs))**2+(2*maxval(node(:nd)%add)**2) !this is the maximal interaction distance between epithelial nodes
  mymaxval=0d0
  do i=1,nd
    if(maxval(ADDe(i,:)) > mymaxval)mymaxval=maxval(ADDe(i,:))
  end do
  maxlen=2*mymaxval**2 !this is the maximal interaction distance between epithelial nodes
  rv=2*mymaxval !maximal interaction distance between mesenchymal nodes
  urv=1.0d0/rv

  !THIS IS JUST FOR CLARITY, TO SIMPLIFY EACH SUBROUTINE AND MAKE IT MORE READEABLE 
  if (ffu(9)==0) then              !!>> HC 14-1-2021
     if (ffu(8)==0) then !normal neighboring, extensive search of the boxes
        if (ffu(2)==0) then
           if (ffu(23)==1) then        !!>> HC 14-1-2021
              urv=1.d0/(arrel_de_dos*rv+1d-3) ! beware, that could produce troubles, somehow, I have not check it carefully
                                            ! COMMENT: that was mildly wrong (or just inefficient) in all previous versions
              call iniboxes_p_pola
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
              call iniboxes_p_pola
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
        call iniboxes_p_pola
        if (ffu(2)==0) then !no GABRIEL
          call neighbor_build_triang_nogabriel    
        else
          call neighbor_build_triang_gabriel
        end if
     end if 
  else
     call neighbor_build_old
  end if
  call fill_pol_distances
end subroutine

!*************************************************************************************************************************************************************************

subroutine neighbor_build_node_pola(i)
  integer i
  real*8 a
  real*8 :: mymaxval

  call extrem
  call fill_pol_distances_node(i)
  !maxlen=sqrt((2*maxval(node(:nd)%eqs))**2+(2*maxval(node(:nd)%add)**2)) !this is the maximal interaction distance between epithelial nodes
  mymaxval=0d0
  do i=1,nd
    if(maxval(ADDe(i,:)) > mymaxval)mymaxval=maxval(ADDe(i,:))
  end do
  maxlen=2*mymaxval**2 !this is the maximal interaction distance between epithelial nodes
  a=2*mymaxval !maximal interaction distance between mesenchymal nodes
  rv=a
  urv=1.0d0/a

  !THIS IS JUST FOR CLARITY, TO SIMPLIFY EACH SUBROUTINE AND MAKE IT MORE READEABLE 
  if (ffu(9)==0) then        !!>> HC 14-1-2021
    if (ffu(8)==0) then !normal neighboring, extensive search of the boxes
      if (ffu(2)==0) then
        if (ffu(23)==1) then  !!>> HC 14-1-2021
          urv=1.d0/(arrel_de_dos*rv+1d-3)   ! COMMENT: that was mildly wrong (or just inefficient) in all previous versions   
          call iniboxes_p_pola
          if (ffu(17)==0) then
            call neighbor_build_complexneigh_node(i)
          else
            call neighbor_build_complexneigh_inter_other_side_epi_node(i)
          end if
        else  
          urv=1.0d0/(rv+1d-3)      
          call iniboxes_p_pola
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
      call iniboxes_p_pola
      if (ffu(2)==0) then !no GABRIEL
        call neighbor_build_node_old(i)  ! this may need to be changed if this is used often, that is probably not    
      else
        call neighbor_build_node_old(i)  ! this may need to be changed if this is used often, that is probably not    
      end if
    end if 
  else
    call neighbor_build_node_old(i)
  end if
  
  call fill_pol_distances_node(i)
  
end subroutine 



end module polarization
