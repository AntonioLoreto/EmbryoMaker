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
module model
use general
use aleas
use neighboring
use io
!use creix	!>>Miquel 14-10-12
use shell
use nexus
use biomechanic
use genetic
use energy ! >>> Is 5-6-14
use polarization
use biomechanic_pola

public:: iteracio

integer, public :: passed
integer, public :: getott !>>>Miguel 8-10-14
integer,public  :: lock   !>>>Miguel29-10-14
real*8,  public, allocatable  :: gext1(:,:,:),gext2(:,:,:),ejex1(:),ejex2(:) ! to store gex(:,:) in time lapses ...!>>>Miguel 8-10-14


contains

!**************************************************************************
subroutine eneinicial
integer::i

  node(:)%e=0
  do i=1,nd
    call energia(i)
  end do
end subroutine eneinicial

!**************************************************************************
subroutine prints
!print *,"getot",getot,"real time",rtime,"nd",nd,"delta",delta,"deltamin",deltamin,node(nd)%x,minval(node(:nd)%dif)," file ",carg,maxval(gex(:,5))
print *,"getot",getot,"real time",rtime,"nd",nd!," file ",carg,delta,deltamin,cels(1)%fase
!do i=1,ng
!  print *,i,maxval(gex(:,i))
!end do

end subroutine prints

!**************************************************************************


!**************************************************************************

subroutine iteracio(tf)
integer tf,itf,il,status
real*8::rtf,realtf  !>>>>>>>>>>>Miquel 17-6-13
integer::cont       !>>>>>>>>>>>Miquel 17-6-13
real*8::ox,oy,oz !miguel4-1-13
real*8::rad,nrad,ax,ay,az
character*200 cxhc  !!>> HC 20-2-2021

  realtf=real(tf)   !>>>>>>>>>>>Miquel 17-6-13
  rtf=0d0           !>>>>>>>>>>>Miquel 17-6-13
  cont=0            !>>>>>>>>>>>Miquel 17-6-13
  geu=0
  itf=0


  if(getot.eq.0)then                                        !>>>Miguel29-10-14
    lock=0 ! matrices gext are filled by default            !>>>Miguel29-10-14
    if((aut.ne.1).and.(aut.ne.5))then;call printgex1;endif  !>>>Miguel 8-10-14
  end if                                                    !>>>Miguel29-10-14


  if (itviactual<itvi) then
    call go_iteration_forth(tf)  !recovers remembered iterations before runing new ones
  end if

  do  !*******this loop does the number of iterations equat to tf

    getot=getot+1
    geu=geu+1
    itf=itf+1

    !UPDATING CELL CENTROIDS
    do i=1,ncels
      a=0 ; b=0 ; c=0
      do j=1,cels(i)%nunodes
        k=cels(i)%node(j)
        if(node(k)%tipus==1.or.node(k)%tipus==3)then
          a=a+node(k)%x ; b=b+node(k)%y ; c=c+node(k)%z
        end if
      end do
      d=1d0/real(cels(i)%nunodes)
      if(node(cels(i)%node(1))%tipus<3) d=2d0*d !if it's epithelial
      cels(i)%cex=a*d ; cels(i)%cey=b*d ; cels(i)%cez=c*d
    end do

    rdiffmax=2*maxval(node(:nd)%add)*dmax !>>Miquel27-2-14

    !calculating neighbors
    if(nd>1)then
      if(ffu(28)==1 .or. npag(nparam_per_node+17) .eq. 0)then
        call neighbor_build
      else
        call neighbor_build_pola
      end if
    end if
    if(ffu(24)==0)  call restore_neighbors
    


    if(ffu(16)==0)then
      if(ffu(28)==1 .or. npag(nparam_per_node+17) .eq. 0)then
        call iterdiferencial
      else
        call iterdiferencial_pola
      end if
    else !forces is disabled
      delta=deltamin
    end if


    call gene_stuff    !ffu ensemble


    if (ffu(26)==1)then        !!>> HC 15-2-2021 
       call nexe               !!>> HC 15-2-2021  !nexe should be first >>> Is 13-2-14
    else                       !!>> HC 15-2-2021 
       call nexe_gradual       !!>> HC 15-2-2021 Gradual implementation of genetic control of cell properties 
    endif                      !!>> HC 15-2-2021 (it comes from nexe 2 with some 2021 changes)
  
    

    !UPDATINGS********
    !node positions from forces
    if (ffu(12)==0) then
      if (ffu(5)==1) then !euler numerical integration
        do i=1,nd                                       ! miguel4-11-13
          node(i)%x=node(i)%x+delta*px(i)               ! miguel4-11-13
          node(i)%y=node(i)%y+delta*py(i)               ! miguel4-11-13
          node(i)%z=node(i)%z+delta*pz(i)               ! miguel4-11-13
        end do
      else  ! Runge-Kutta order 4 numerical integration
        if(ffu(28)==1 .or. npag(nparam_per_node+17) .eq. 0)then
          call rungekutta4(delta)
        else
          call rungekutta4_pola(delta)
        end if
      end if
    else
      if(ffu(28)==1 .or. npag(nparam_per_node+17) .eq. 0)then
        call adaptive_rungekutta
      else
        call adaptive_rungekutta_pola
      end if
    end if

    !gene expression
    gex(1:nd,1:ng) = agex(1:nd,1:ng)
    where(gex.lt.0) gex=0.0d0


    !RANDOM NOISE
    if(ffu(16)==0)then
      c=nd*prop_noise*delta/deltamin !now the proportion of nodes is still dynamic, but equal to prop_noise on default  !>>Miquel28-7-14
     !print*,"proportion c",c
      if (c>1) then
        do il=1,int(c)
          call itera                  !****we add noise, though the behaviour of the system it'
        end do
      else
        call random_number(a)
        if (a<c) then
          call itera
        end if
      end if
    else  !forces disabled
      do il=1,nd
        call itera
      end do
    end if  
    
    ! WHICHEND 4: TOO MANY NODES 
   ! whichend=0                                                             !!>> HC 26-2-2021 This is already initialized in elli
    if ((ffu(21)==1.and.ffufi(4,1)==1)) then                                !!>> HC 12-7-2021
       if (mod(getot,ffufi(4,3)).eq.0)then                                  !!>> HC 20-2-2021
          if (nd>ndmax) then                                                !!>> HC 20-2-2021
             out_of_control=1                                               !!>> HC 20-2-2021
             status=10 !implies a 20560 exit status  Is 1-10-14             !!>> HC 20-2-2021
             if (ffu(22)==0)then                                            !!>> HC 7-1-2021 FILTER OF MAXIMUM NUMBER OF NODES
                print *,"too many cells, I quit"                            !!>> HC 7-1-2021 Silent mode
                write(0,*) nd,"nd",ndmax,"ndmax too many cells, I quit"     !!>> HC 7-1-2021
             endif                                                          !!>> HC 7-1-2021
             whichend=4                                                     !!>> HC 7-1-2021
          endif                                                             !!>> HC 20-2-2021
       endif                                                                !!>> HC 20-2-2021
    endif                                                                   !!>> HC 20-2-2021

    ! WHICHEND 5: NO NODES   
    if (ffu(21)==1)then                                                     !!>> HC 20-2-2021
       if (ffufi(5,1)==1.and.ffufi(5,3)>0)then                              !!>> HC 20-2-2021
          if (mod(getot,ffufi(5,3)).eq.0)then                               !!>> HC 20-2-2021
             if (nd<1.or.ncels<1)then                                       !!>> HC 20-2-2021
                out_of_control=1                                            !!>> HC 20-2-2021
                status=10 !implies a 20560 exit status  Is 1-10-14          !!>> HC 20-2-2021
                if (ffu(22)==0)then                                         !!>> HC 7-1-2021 FILTER OF NO NODES
                   print *,"too few nodes, they all die"                    !!>> HC 7-1-2021 Silent mode
                   write(0,*) nd,"nd",ndmax,"too few cells, they all die"   !!>> HC 7-1-2021
                endif                                                       !!>> HC 7-1-2021
                whichend=5                                                  !!>> HC 7-1-2021
             endif                                                          !!>> HC 20-2-2021
          endif                                                             !!>> HC 20-2-2021
       endif                                                                !!>> HC 20-2-2021
    end if                                                                  !!>> HC 20-2-2021
   
    if (whichend==4 .or. whichend==5)then                                   !!>> HC 20-2-2021        
       open(616, file=trim(carg)//".loga")                                   !!>> HC 20-2-2021 We want to know
       write(616,*) "reason to end:",whichend                               !!>> HC 17-11-2020 why we killed this indv
       close(616)                                                           !!>> HC 17-11-2020
       cxhc="pwd >> "//trim(carg)//".logb"                                  !!>> HC 29-9-2021 Save the name of the individual (evolution)
       call system(cxhc)                                                    !!>> HC 29-9-2021
       cxhc="paste "//trim(carg)//".logb "//trim(carg)//".loga > "//trim(carg)//".log"   !!>> HC 29-9-2021 paste the name and the reason to end
       call system(cxhc)                                                    !!>> HC 29-9-2021 Remove the temporal files
       cxhc="rm "//trim(carg)//".logb"                                      !!>> HC 29-9-2021
       call system(cxhc)                                                    !!>> HC 29-9-2021
       cxhc="rm "//trim(carg)//".loga"                                      !!>> HC 29-9-2021
       call system(cxhc)                                                    !!>> HC 29-9-2021
       if (ffufi(whichend,2)==1)then                                        !!>> HC 20-2-2021  If the end is lethal 
           cxhc='echo "0.00" > '//trim(carg)//'fitness'                     !!>> HC 17-11-2020 We kill the embryo 
           call system(cxhc)                                                !!>> HC 17-11-2020 and it has no fitness
       endif                                                                !!>> HC 17-11-2020

       call writesnap                                                       !!>> HC 7-1-2021
       
       call exit(status)                                                    !!>> HC 20-2-2021
       stop                                                                 !!>> HC 20-2-2021
    endif                                                                   !!>> HC 20-2-2021

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rtf=rtf+delta
    rtime=rtime+delta

    call put_param_to_matrix(param)

    if (mod(getot,fprint)==0 .and. ffu(22)==0) call prints  !!>> HC 30-11-2020 I have added a silent mode 

    if (mod(getot,freqsnap).eq.0) then ; if (fsnap==1) then ;call writesnap; end if ; end if

    if (mod(getot,100000).eq.0) print *,getot,nd,"getot nd",node(nd)%x

    if (ffu(7)==0) then
      if (ffu(14)==0) then
!        print *,"the actual time is",rtf,"delta is",delta
        if (rtf>=tf) then                   
          if (ffu(22)==0)then              !!>> HC 30-11-2020 I have added a silent mode
             print *,"the real time is ",rtime,"we now runned",rtf,"the number of iterations runned are",itf,"and",getot,&
             "in total : we have",nd,"nodes and",ncels,"cells" 
          endif                            !!>> HC 30-11-2020 
          exit 
        end if
      else
!        print *,delta,"delta",nd,"nd"
        if (itf==tf) then                   
           if (ffu(22)==0)then             !!>> HC 30-11-2020 I have added a silent mode
              print *,"real time is",rtf,"to run a time of ",tf,"delta",delta,"total number of iterations is",getot,&
              "nd",nd,"ncels",ncels
           endif                           !!>> HC 30-11-2020 
           exit
        end if
      end if
    else 
      if (itf==tf) then                    
        if (ffu(22)==0)then                !!>> HC 30-11-2020 I have added a silent mode 
           print *,"the number of iterations is",getot," and we just runned",itf,"we have",nd,"nodes and",ncels,"cells" 
        endif                              !!>> HC 30-11-2020 
        exit
      end if
    end if

  end do

  if((aut.ne.1).and.(aut.ne.5).and.(lock.eq.0))then  !>>>Miguel8-10-14   !>>>Miguel29-10-14
    call printgex1                   !>>>Miguel8-10-14
  end if                             !>>>Miguel8-10-14



!  print*,"rtime",rtime,"RTIME",delta,"delta",rtime,"nd",nd,"ncels",ncels

  if (itviactual==itvi) then
    itvi=itvi+1
    if (itvi>mamax) itvi=1
    pnode(itvi,:)=node
    call put_param_to_matrix(param)
    pparam(itvi,:)=param
    pvarglobal_out(itvi,:)=varglobal_out
    itviactual=itvi
  end if

  passed=1

end subroutine

!**************************************************************************
subroutine itera
integer   ::i,accepta,val,ic,t,icc,iccc,nodmo,celi,nnod
real*8    ::ox,oy,oz,oe,de,ae,ocx,ocy,ocz,ax,ay,az,kl
integer   ::nr  !>>> Is 4-3-14

    accepta=0
    recub=0
    movi=0             !>>>> miguel23-7-13

    !escull un node
    !cotran=cotran+1 ;if (cotran>nunualea) call llaleat; a=stocas(cotran)       
    !nodmo=int(a*nd)+1
    nr=0
10  call random_number(a)
    nr=nr+1   !>>> Is 4-3-14
    nodmo=int(a*nd)+1 !ORIGNAL

    if (nr>nd*2) return             !>>> Is 4-3-14
    if(node(nodmo)%fix==2) goto 10 !>>> Is 4-3-14

    if(ffu(15)==1) call energia(nodmo)	!energia abans de moure'l !only when noise by energy is activated !>>Miquel28-7-14

    oe=node(nodmo)%e
    !mou-lo
    ox=node(nodmo)%x
    oy=node(nodmo)%y
    oz=node(nodmo)%z

    if (node(nodmo)%tipus<4) then
      celi=node(nodmo)%icel 
      nnod=cels(celi)%nunodes
      ocx=cels(celi)%cex		!the old centroid
      ocy=cels(celi)%cey
      ocz=cels(celi)%cez
    end if

    !cotran=cotran+1 ;if (cotran>nunualea) call llaleat; a=stocas(cotran)       
    !desplacament=a*node(nodmo)%dmo
    !a=ran2(idum);
    call random_number(a)
    if (ffu(9)==1) then                  ! 4-3-2020
      desplacament=a*node(nodmo)%dmo       ! 4-3-2020
    else                                   ! 4-3-2020
      desplacament=a*node(nodmo)%dmo*delta ! 4-3-2020 THIS ONE IS MORE CORRECT    
    end if
!if (node(nodmo)%tipus==3) print *,nodmo,node(nodmo)%dmo,a,desplacament,"ds"
    if(desplacament<epsilod) return !no need to run all the energies if the movement is going to be 0 !>>Miquel27-8-14
  !insert biased noise here
    if(npag(nparam_per_node+16)>0.and.node(nodmo)%tipus<4) then !>>> IS 10-5-14
      a=0d0
      do k=1,npag(nparam_per_node+16)
        kk=whonpag(nparam_per_node+16,k)
        if (gex(nodmo,kk)>0.0d0) then
          a=a+gex(nodmo,kk)*gen(kk)%e(nparam_per_node+16)  !wa in units of probability
        end if
      end do
      call random_number(b)
      k=int(b*nvaloq)+1
      ax=particions_esfera(k,1)+a*cels(celi)%polx  !we need to use a unit vector !>>Miquel28-7-14
      ay=particions_esfera(k,2)+a*cels(celi)%poly
      az=particions_esfera(k,3)+a*cels(celi)%polz
      d=1d0/sqrt(ax**2+ay**2+az**2)
      node(nodmo)%x=node(nodmo)%x+ax*d*desplacament
      node(nodmo)%y=node(nodmo)%y+ay*d*desplacament
      node(nodmo)%z=node(nodmo)%z+az*d*desplacament 
    else
      call random_number(a)
      k=int(a*nvaloq)+1
      node(nodmo)%x=node(nodmo)%x+particions_esfera(k,1)*desplacament
      node(nodmo)%y=node(nodmo)%y+particions_esfera(k,2)*desplacament
      node(nodmo)%z=node(nodmo)%z+particions_esfera(k,3)*desplacament 
    end if

    if(node(nodmo)%tipus==1)then   !>>>>Miquel 11-6-13 faster way to recalculate centroid
      a=2d0/real(nnod)
      cels(celi)%cex=cels(celi)%cex+(node(nodmo)%x-ox)*a	!recalculate the centroid according to nod's new position
      cels(celi)%cey=cels(celi)%cey+(node(nodmo)%y-oy)*a
      cels(celi)%cez=cels(celi)%cez+(node(nodmo)%z-oz)*a
    else if(node(nodmo)%tipus==3)then
      a=1d0/real(nnod)
      cels(celi)%cex=cels(celi)%cex+(node(nodmo)%x-ox)*a	!recalculate the centroid according to nod's new position
      cels(celi)%cey=cels(celi)%cey+(node(nodmo)%y-oy)*a
      cels(celi)%cez=cels(celi)%cez+(node(nodmo)%z-oz)*a
    end if

    if(ffu(15)==0.and.ffu(16)==0) return  !if unbiased nise is activated, we skip all the energy part, the node is simply moved randomly !>>Miquel28-7-14
    !Calculate energy for new position
    if(nd>1) call neighbor_build_node(nodmo)   ! >>> Is 29-6-14  !>> HC 27-5-2020 !>> HC 16-9-2020
    if(ffu(24)==0) call restore_neighbors
    !call neighbor_build  ! >>> Is 29-6-14  ! THIS COULD BE OPTIMIZED

    call energia(nodmo)
    !acceptes?
    ae=node(nodmo)%e-oe
    if(movi.eq.1)then; goto 432 ; endif   !miguel4-11-13 
    if(ae<-epsilod)then
      accepta=1
      itacc=itacc+1
    else
      !cotran=cotran+1 ;if (cotran>nunualea) call llaleat; g=stocas(cotran)
      call random_number(a)
      !a=ran2(idum)
      kl=temp+node(nodmo)%mov  ! >>> Is 10-10-14
      if (kl<0) kl=epsilod    ! >>> Is 10-10-14
!if (node(nodmo)%tipus==3) print *,nodmo,node(nodmo)%mov,a,desplacament,"dsr",nue**(-ae/kl)
      if(a<nue**(-ae/kl)) then ! >>> Is 10-10-14         
        !if(movi.eq.0)then    !miguel22-7-13           ! miguel4-11-13
          accepta=1
          itacc=itacc+1
        !end if                              !miguel22-7-13
      else
432     node(nodmo)%x=ox                    !miguel4-11-13
        node(nodmo)%y=oy
        node(nodmo)%z=oz
        node(nodmo)%e=oe
        if (node(nodmo)%tipus<4) then
          cels(celi)%cex=ocx		!reestablish the old centroid
          cels(celi)%cey=ocy
          cels(celi)%cez=ocz
        end if
        !call neighbor_build  ! >>> Is 29-6-14 !THIS COULD BE OPTIMIZED
        if(nd>1) then
          if(ffu(28)==1 .or. npag(nparam_per_node+17) .eq. 0)then
            call neighbor_build_node(nodmo)   ! >>> Is 29-6-14  !>> HC 27-5-2020 !>> HC 16-9-2020
          else
            call neighbor_build_node_pola(nodmo)
          end if
        end if
        if(ffu(24)==0) call restore_neighbors

        goto 150
      end if
    end if
    !if(ffu(8)==0)then !only necessary when classic neighboring !this is not necessary anymore !>>Miquel7-3-14
    !  if (accepta==1) then
    !    if(abs(node(nodmo)%x*urv)>=nboxes-1.or.abs(node(nodmo)%y*urv)>=nboxes-1.or.abs(node(nodmo)%z*urv)>=nboxes-1) then		!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
    !      call iniboxesll  ;print*,"we remake the boxes"!we got out of the grid and have to remake the grid
    !    else
    !      icc=boxes(nint(ox*urv),nint(oy*urv),nint(oz*urv))	                            !old position boxes		
    !      iccc=boxes(nint(node(nodmo)%x*urv),nint(node(nodmo)%y*urv),nint(node(nodmo)%z*urv))!new position boxes
    !      if (icc/=iccc) then !Return the node to old box
    !        ic=icc
    !        ! we take away the node from its old cube
    !        if (ic==nodmo) then ;print*,"first in the box"
    !          boxes(nint(ox*urv),nint(oy*urv),nint(oz*urv))=list(nodmo)
    !        else
    !          print*,"ic comença",ic
    !          do ; i=ic ; ic=list(i) ; if(ic==nodmo) exit ; end do
    !          list(i)=list(ic)														
    !        end if																		
    !        ! we put it in the new cub															
    !        ic=iccc																		
    !        if (ic==0) then   !the cube was empty
    !          boxes(nint(node(nodmo)%x*urv),nint(node(nodmo)%y*urv),nint(node(nodmo)%z*urv))=nodmo
    !          list(nodmo)=0
    !        else		
    !          !the cube was not empty so we search for the end																
    !          do
    !            i=ic ; ic=list(i) ; if(ic==0)then ; list(i)=nodmo ; list(nodmo)=0 ; exit ; end if
    !          end do
    !        end if
    !      end if
    !    end if
    !  end if
    !end if
    
150 continue

end subroutine itera

!*******************************************SUBROUTINE********************************************************
subroutine go_iteration_back(it)  !FUNKY MODULE POSITION IT IS CALLED FROM PINTA AND DOES NOT RUN SIMULATION
  integer it
  itviactual=itviactual-it
  if (itviactual<1) then ; itviactual=1 ; print *,"you've gone too far" ; endif
  call move_iteration(it)
end subroutine

!*******************************************SUBROUTINE********************************************************
subroutine go_iteration_forth(it)  !FUNKY MODULE POSITION IT IS CALLED FROM PINTA AND DOES NOT RUN SIMULATION
  integer it
  itviactual=itviactual+it
  if (itviactual>itvi) then ; itviactual=itvi; print *,"now we are were we were" ; endif
  call move_iteration(it)
end subroutine

!*******************************************SUBROUTINE********************************************************
subroutine move_iteration(it)
  node(:nd)=pnode(itviactual,:nd)
  param=pparam(itviactual,:)
  varglobal_out=pvarglobal_out(itviactual,:)
  call get_param_from_matrix(param)
  print *,"read ",getot,"getot",itviactual,"itviactual",itvi,"itvi"
end subroutine

!*************************************************************************************************************
subroutine ordenarepe(ma,mt,rang)
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
end subroutine ordenarepe

subroutine printgex1 !>>>Miguel 8-10-14 made subroutine
integer :: faktor,col1,col2
faktor=10

if(getot.eq.0)then   ; getott=1 
else;getott=getott+1 ; endif         

if(getott.eq.1)then 
   allocate(gext1(nd,ng,faktor))  
   allocate(ejex1(faktor))  
end if

col1=size(gext1(1,1,:)) ; col2=size(gext1(:,1,1))

if((getott.gt.col1).or.(nd.gt.col2))then
  allocate(gext2(nd,ng,col1),ejex2(col1))   
  gext2=gext1 ; ejex2=ejex1
  deallocate(gext1,ejex1) 
  allocate(gext1(nd,ng,getott+faktor),ejex1(getott+faktor))
  gext1(1:col2,1:ng,1:col1)=gext2(1:col2,1:ng,1:col1)
  ejex1(1:col1)=ejex2(1:col1)
  deallocate(gext2,ejex2)
end if
gext1(1:nd,1:ng,getott)=gex(1:nd,1:ng)
ejex1(getott)=getot

end subroutine printgex1

end module model
