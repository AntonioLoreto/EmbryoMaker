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



module death    ! ROLAND

use general
use neighboring
use io

integer,public:: doomsday,deemed,toosmall ! moved toggle: Roland 25-9-13

contains

subroutine should_I_die !>>>>>>>> Miquel11-6-13 it's a copy of growth_gradual, but the opposite
integer::i,j,k,ii,jj,kk,iii,jjj,kkk,l,nnod ! RZ 17-11-14 added nnod
real*8::grate,c
real*8::minreq,maxreq !these should be implementation parameters, it's the minimum size of the node at birth and death, and max size for growth
integer,dimension(:)::discel(ncels+1)
integer::tipi, naltre !>>> TT 14-05-2020 added tipi to check cell type 

  minreq=0.0001d0;   !this is an arbitrary value that should not matter much
  maxreq=mmae
  iii=0 !counter of how many cells are disappearing
  discel=0

  do i=1,ncels
    if(node(cels(i)%node(1))%fix==1) cycle !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
    nnod=cels(i)%nunodes
    jjj=0
    tipi=cels(i)%ctipus !>>> TT 14-05-2020
    do ii=1,nnod
      j=cels(i)%node(ii)
      grate=0d0
      c=1-node(j)%dif
      do jj=1,npag(nparam_per_node+3)  !number of genes affecting req
        k=whonpag(nparam_per_node+3,jj)  !which are those genes
        grate=grate+gex(j,k)*gen(k)%e(nparam_per_node+3)*c !this is the differential of death for the node !wa in space-req units, kind of
      end do
      grate=grate*delta
      if(grate==0) cycle
      if(node(j)%grd>=maxreq)then  !>>> Is 17-3-14
        !grac=grac+grate
      else
        a=node(j)%grd ; !b=node(j)%codel
        !c=node(j)%add-a
        node(j)%grd=a-grate
        if(node(j)%grd<minreq)then   !this is the maximal req when growing (the same as for invagination: recycling)
          d=minreq-node(j)%grd
          !grac=grac+d   !if the node is "full" we save the growth to apply it to another node
          if(tipi<3) naltre = node(j)%altre  !!>> HC 11-6-2020 !>>> TT 14-05-2020 
          call apoptosis(j)
          if (ffu(24)==0)then                                       !!>> HC 6-7-2021
             if (tipi<3) then                                       !!>> HC 30-6-2021
                call neighbor_build                                 !!>> HC 6-10-2021
                call fill_co_griders                                !!>> HC 30-6-2021
             else                                                   !!>> HC 30-6-2021
                call neighbor_build                                 !!>> HC 6-10-2021
                call fill_co_griders                                !!>> HC 6-10-2021
             end if                                                 !!>> HC 30-6-2021
          endif                                                     !!>> HC 6-7-2021
          if (nd<0) then 
            print *,"THIS IS THE END MY FRIEND: NO NODES LEFT"; 
            open(23,file=trim(carg)//"e")
            print *,"making...",trim(carg)//"e"
            write(23,*) trim(carg)//trim(nofi)
            close(23)
            call exit(status) 
          end if
          if(cels(i)%nunodes==0.and.discel(i)==0) then  !>>> Is 10-5-14
            discel(i)=i                                 !>>> Is 10-5-14
            iii=iii+1                                     !>>> Is 10-5-14
          end if ! this means we have to delete this cell i 
          cycle
        end if
        !node(j)%add=node(j)%eqd+c
        nodeo(j)%grd=node(j)%grd       ! when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
        !nodeo(j)%add=nodeo(j)%add-grate     !values have to be updated, otherwise it will create a conflict with nexus
        jjj=jjj+1
        if(node(j)%tipus<3)then
          k=node(j)%altre
          a=node(k)%grd
          !c=node(k)%add-a
          node(k)%grd=a-grate
          !node(k)%add=node(k)%eqd+c
          nodeo(k)%grd=node(k)%grd  !when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
          !nodeo(k)%add=nodeo(k)%add-grate    !values have to be updated, otherwise it will create a conflict with nexus
        end if
!        nodsmall(iii)=j
      end if
    end do
    if(jjj==0)then
      call random_number(a)
      k=int(a*nnod)+1
      j=cels(i)%node(k)
      a=node(j)%grd 
      !c=node(j)%add-a
      node(j)%grd=a-grate
      !node(j)%add=node(j)%eqd+c
      nodeo(j)%grd=node(j)%grd      !when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
      !nodeo(j)%add=nodeo(j)%add-grate        !values have to be updated, otherwise it will create a conflict with nexus
    end if
  end do

  ! here we eliminate the cells that have lost all their nodes, there are iii of them
  if (iii>0) then
    do k=1,ncels
23    if (discel(k)==0) cycle !>>> Is 11-4-14
      do ii=k,ncels-1         !>>> Is 11-4-14
        cels(ii)%nunodes=cels(ii+1)%nunodes
        cels(ii)%cex=cels(ii+1)%cex ; cels(ii)%cey=cels(ii+1)%cey ; cels(ii)%cez=cels(ii+1)%cez
        cels(ii)%polx=cels(ii+1)%polx ; cels(ii)%poly=cels(ii+1)%poly ; cels(ii)%polz=cels(ii+1)%polz
        cels(ii)%ctipus=cels(ii+1)%ctipus
        discel(ii)=discel(ii+1)
        if(cels(ii)%nodela<cels(ii+1)%nodela)then
          deallocate(cels(ii)%node)
          cels(ii)%nodela=cels(ii+1)%nodela
          allocate(cels(ii)%node(cels(ii)%nodela))
          cels(ii)%node(1:cels(ii)%nunodes)=cels(ii+1)%node(1:cels(ii)%nunodes) !!! ii+1 may be undefined
        else
          cels(ii)%node(1:cels(ii)%nunodes)=cels(ii+1)%node(1:cels(ii)%nunodes)
        end if
      end do
      cels(ncels)%nunodes=0
      cels(ncels)%cex=0 ; cels(ncels)%cey=0 ; cels(ncels)%cez=0
      cels(ncels)%polx=0 ; cels(ncels)%poly=0 ; cels(ncels)%polz=0
      cels(ncels)%ctipus=0
      cels(ncels)%nodela=0
      deallocate(cels(ncels)%node)
      ncels=ncels-1
      if (ncels==0) then; print *,"THIS IS THE END MY FRIEND, NO CELLS LEFT" ; 
        open(23,file=trim(carg)//"e")
        print *,"making...",trim(carg)//"e"
        write(23,*) trim(carg)//trim(nofi)
        close(23)
        call exit(status) ; 
      end if
      goto 23 
    end do

    do k=1,ncels
      do ii=1,cels(k)%nunodes
        node(cels(k)%node(ii))%icel=k
      end do
    end do
  end if


  if(ndx>0)then  !apoptosis for ECM nodes (proteolisis) !>>>>Miquel16-12-13
    do j=1,nd    !***************************this could be highly optimizable if we kept a list of the ECM nodes**************************
      if(node(j)%tipus==4)then
        grate=0d0
        c=1-node(j)%dif
        do jj=1,npag(nparam_per_node+14)  !number of genes affecting req
          k=whonpag(nparam_per_node+14,jj)  !which are those genes
          grate=grate+gex(j,k)*gen(k)%e(nparam_per_node+14)*c !this is the differential of death for the node
        end do                                                 ! wa in units of space-req
        grate=grate*delta
        if(grate==0) cycle
        a=node(j)%grd ;
        !c=node(j)%add-a
        node(j)%grd=a-grate
        if(node(j)%grd<minreq)then   !!>> HC 17-6-2020 its grd not eqd here !this is the maximal req when growing (the same as for invagination: recycling)
          d=minreq-node(j)%grd
          !grac=grac+d   !if the node is "full" we save the growth to apply it to another node
          call apoptosis(j)
          if (ffu(24)==0)then                                    !!>> HC 6-7-2021
             call neighbor_build                                 !!>> HC 6-10-2021
             call fill_co_griders                                !!>> HC 6-10-2021         
          endif                                                  !!>> HC 6-10-2021
          cycle
        end if
        !node(j)%add=node(j)%eqd+c

        nodeo(j)%grd=node(j)%grd   !when req and da are irreversively modified by growth or apoptosis, the nodeo  !>>>>Miquel16-12-13
        !nodeo(j)%add=nodeo(j)%add-grate    !values have to be updated, otherwise it will create a conflict with nexus
      end if
    end do
  end if
end subroutine

!***************************************************************************************************

subroutine eliminate_cell(cell)
  integer cell
  integer i,ii,iii,j,jj,jjj,k,kk,kkk
    do ii=cell,ncels-1
      cels(ii)%nunodes=cels(ii+1)%nunodes
      cels(ii)%cex=cels(ii+1)%cex ; cels(ii)%cey=cels(ii+1)%cey ; cels(ii)%cez=cels(ii+1)%cez
      cels(ii)%polx=cels(ii+1)%polx ; cels(ii)%poly=cels(ii+1)%poly ; cels(ii)%polz=cels(ii+1)%polz
      cels(ii)%ctipus=cels(ii+1)%ctipus
      if(cels(ii)%nodela<cels(ii+1)%nodela)then
        deallocate(cels(ii)%node)
        cels(ii)%nodela=cels(ii+1)%nodela
        allocate(cels(ii)%node(cels(ii)%nodela))
        cels(ii)%node(1:cels(ii)%nunodes)=cels(ii+1)%node(1:cels(ii)%nunodes) !!! ii+1 may be undefined
      else
        cels(ii)%node(1:cels(ii)%nunodes)=cels(ii+1)%node(1:cels(ii)%nunodes)
      end if
    end do
    cels(ncels)%nunodes=0
    cels(ncels)%cex=0 ; cels(ncels)%cey=0 ; cels(ncels)%cez=0
    cels(ncels)%polx=0 ; cels(ncels)%poly=0 ; cels(ncels)%polz=0
    cels(ncels)%ctipus=0
    cels(ncels)%nodela=0
    deallocate(cels(ncels)%node)
    ncels=ncels-1

    !!!!! Now, the updated part: RZ 30-07-13
    do j=1,ncels
      do l=1,cels(j)%nunodes
        if(j.ge.cell)then
          node(cels(j)%node(l))%icel=node(cels(j)%node(l))%icel-1
        endif
      end do
    end do
    !!!!!!

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine apoptosis(i)

integer:: j,jj,j1,j2,ic,icc,toggle,togglea,i1,i2,i3 ! toggle: epithelium/mesenchyma distinction
integer:: plych                                     !!>> HC 17-6-2020
    if (nd==0) then 
      print *,"no nodes left: I quit"  
      open(23,file=trim(carg)//"t")
      print *,"making...",trim(carg)//"t"
      write(23,*) trim(carg)//trim(nofi)
      close(23)
      call exit(status)
    end if

    !if (node(i)%tipus==4) return !!>> HC 17-6-2020 Proteolysis implemented

    if(node(i)%tipus<3) then;
      toggle=2 !epi
      plych=0  !!>> HC 18-6-2020 Proteolysis flag OFF
    else
      toggle=1 !non-epi
      if (node(i)%tipus==4)then; plych=1; else; plych=0; endif !!>> HC 17-6-2020 Proteolysis flag
    endif 

    ! >>> Is 11-6-14 
    if (ffu(8)==0) then
      do ic=1,nd
        if(list(ic).eq.i)then
          list(ic)=list(i);
          cycle
        end if
      end do

      if(toggle.eq.2) then
        do ic=1,nd
          if(list(ic).eq.node(i)%altre)then
            list(ic)=list(node(i)%altre);
            cycle;
          end if
        end do
      end if
    end if ! >>> Is 11-6-14

    nd=nd-toggle
    k=node(i)%icel
    j1=i;j2=node(i)%altre

    if(toggle.eq.2)then

      if(node(i)%altre.lt.i)then;j1=node(i)%altre;j2=i;endif

      do j=j2,nd+2
        node(j)=node(j+1)    !>>> Is 4-1-14 
        gex(j,:)=gex(j+1,:)  !>>>>>>> Miquel18-7-13
        agex(j,:)=agex(j+1,:)!>>> Is 13-5-14
      end do

      do j=1,nd+2
        if (node(j)%altre>j2-1) node(j)%altre=node(j)%altre-1 ! Is 4-1-4
      end do

      if (ffu(8)==0) then ! >>> Is 15-6-14
        do j=1,nd+1    
          if(j.lt.j2)then
            if(list(j).lt.j2)then;list(j)=list(j)
            elseif(list(j).eq.j2)then
              if(list(j2).lt.j2)then;list(j)=list(j2)
              else;if(list(j2).ne.0)then;list(j)=list(j2)-1;endif
              end if
            else;if(list(j).ne.0)then;list(j)=list(j)-1;endif
            end if
          else
           if(list(j+1).lt.j2)then;
             list(j)=list(j+1)
           elseif(list(j+1).eq.j2)then
             if(list(j2).lt.j2)then;
               list(j)=list(j2)
             else;
               if(list(j2).ne.0)then;list(j)=list(j2)-1;endif
            end if
          else; if(list(j+1).ne.0)then;list(j)=list(j+1)-1;endif
          end if
        end if
      end do
      end if ! >>> Is 15-6-14
    end if

    do j=j1,nd
      node(j)=node(j+1)    !>>> Is 4-1-14
      gex(j,:)=gex(j+1,:)  !>>>>>>> Miquel18-7-13
      agex(j,:)=agex(j+1,:)!>>> Is 13-5-14
    end do

    do j=1,nd
      if (node(j)%altre>j1-1) node(j)%altre=node(j)%altre-1  ! Is 4-1-4
    end do
 
    if (ffu(8)==0) then ! >>> Is 15-6-14
    do j=1,nd+1 
      if(j.lt.j1)then
        if(list(j).lt.j1)then;list(j)=list(j)
        elseif(list(j).eq.j1)then
          if(list(j1).lt.j1)then;list(j)=list(j1)
          else;if(list(j1).ne.0)then;list(j)=list(j1)-1;endif
          end if
        else;if(list(j).ne.0)then; list(j)=list(j)-1;endif
        end if 
      else
        if(list(j+1).lt.j1)then;list(j)=list(j+1)
        elseif(list(j+1).eq.j1)then
          if(list(j1).lt.j1)then;list(j)=list(j1)
          else;if(list(j1).ne.0)then;list(j)=list(j1)-1;endif
          end if
        else; if(list(j+1).ne.0)then;list(j)=list(j+1)-1;endif
        end if
      end if
    end do
    end if

! here was the major problem. The number of boxes actualized wasn't dynamical, so i1,2,3 ran only from -2 to 2.

    if(ffu(8)==0)then
      do i1=-nboxes,nboxes,1 ! the boxes matrix Roland 25-9-13
        do i2=-nboxes,nboxes,1 ! Roland 25-9-13
          do i3=-nboxes,nboxes,1 ! Roland 25-9-13
            ic=boxes(i1,i2,i3)
            if(ic.eq.0)then;cycle;endif  ! Roland 25-9-13
            if(ic.eq.j1)then
              boxes(i1,i2,i3)=list(ic); if(list(ic).ne.0)then; ic=boxes(i1,i2,i3); endif  ! Roland 25-9-13
            end if
            if(ic.eq.j2)then
              boxes(i1,i2,i3)=list(ic); if(list(ic).ne.0)then; ic=boxes(i1,i2,i3); endif  ! Roland 25-9-13
            end if
            if(ic.eq.j1)then
              boxes(i1,i2,i3)=list(ic); if(list(ic).ne.0)then; ic=boxes(i1,i2,i3); endif  ! Roland 25-9-13
            end if
            if(ic.gt.j1)then
              boxes(i1,i2,i3)=ic-1; if(list(ic).ne.0)then; ic=boxes(i1,i2,i3); endif  ! Roland 25-9-13
            end if
            if(ic.gt.j2)then
              boxes(i1,i2,i3)=ic-1; if(list(ic).ne.0)then; ic=boxes(i1,i2,i3); endif  ! Roland 25-9-13
            end if
          end do
        end do
      end do
    end if

    if (plych==0)then !!>> HC 17-6-2020 We do not do this if the node is ecm because it is not in a cell
       ii=k      ! >>> just the comment :: Is 14-4-14 k=node(i)%icel
       togglea=0
       do jj=1,cels(ii)%nunodes-1
         if(cels(ii)%node(jj).eq.j1)then;togglea=1;endif
         if(togglea.eq.0)then
           if(cels(ii)%node(jj).gt.j1)then
             cels(ii)%node(jj)=cels(ii)%node(jj)-1
           endif
         else
           if(cels(ii)%node(jj+1).gt.j1)then
             cels(ii)%node(jj)=cels(ii)%node(jj+1)-1
           else
             cels(ii)%node(jj)=cels(ii)%node(jj+1)
           endif
         endif
       enddo
    endif             !!>> HC 17-6-2020

    do jj=1,toggle
      j=nd+jj
      node(j)%eqs=0
      node(j)%hoo=0
      node(j)%you=0
      node(j)%rep=0
      node(j)%eqd=0
      node(j)%adh=0
      node(j)%rec=0
      node(j)%add=0.25
      node(j)%erp=0
      node(j)%est=0				
      node(j)%tipus=0
      node(j)%altre=j  ! this one may be troublesome, as it defines a value of %altre in unpaired nodes

      gex(j,:)=0   !>>>>Miquel18-7-13
      agex(j,:)=0  !>>> Is 13-5-14

      if (ffu(8)==0)then

        call iniboxes

!        ic=boxes(nint(node(j)%x*urv),nint(node(j)%y*urv),nint(node(j)%z*urv))
!        boxes(nint(node(j)%x*urv),nint(node(j)%y*urv),nint(node(j)%z*urv))=0
!        list(j)=0
      end if
    end do

    if (plych==0)then !!>> HC 17-6-2020 We do not do this if the node is ecm because it is not in a cell
       togglea=0
       if(toggle.eq.2)then
         do jj=1,cels(ii)%nunodes-1
           if(cels(ii)%node(jj).eq.j2-1)then;togglea=1;endif
           if(togglea.eq.0)then
             if(cels(ii)%node(jj).gt.j2-1)then
               cels(ii)%node(jj)=cels(ii)%node(jj)-1
             endif
           else
             if(cels(ii)%node(jj+1).gt.j2-1)then
               cels(ii)%node(jj)=cels(ii)%node(jj+1)-1
             else
               cels(ii)%node(jj)=cels(ii)%node(jj+1)
             endif
           endif
         enddo
       end if  ! >>> Is 14-4-14
 
       if (cels(k)%nunodes>0) then
         cels(k)%nunodes=cels(k)%nunodes-toggle
         do ii=1,ncels
           if (ii.ne.k) then
             do jj=1,cels(ii)%nunodes
               if(cels(ii)%node(jj).ge.j1) then
                 cels(ii)%node(jj)=cels(ii)%node(jj)-1
               end if
               if(cels(ii)%node(jj).ge.j2)then
                 cels(ii)%node(jj)=cels(ii)%node(jj)-(toggle-1)
               end if
             end do
           end if
         end do
       end if
    endif       !!>> HC 17-6-2020
    !!!!end if  ! <<< Is 14-4-14
end subroutine apoptosis

end module death


! TO BE CONSIDERED
! make doomsday a function incorporating a threshold of some apoptosis factor
! save some space by reallocation of the matrices; it's not necessary



