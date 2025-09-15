

!LIMIT ARBITRAR_Y A DIFFUSIO


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




!created by Is 29-4-13

module genetic

use general
!use neighboring

real*8,  public, allocatable  :: gex(:,:)  ! expression level: x is node id, y is the gene
				           ! we put it in a matrix to be faster
real*8,  public, allocatable :: dgex(:,:) ! the increment: only for calculations
integer, public, allocatable :: checkn(:)
real*8 , public, allocatable  :: agex(:,:) !gex matrix used to store updated results !>>Miquel11-3-14
real*8,  public, allocatable  :: kadh(:,:) ! adhesions between adhesion molecules
integer, public               :: ntipusadh ! it is the same than npag(1), it is just an alias for the number of different adhesion molecules

integer              ::nkindof(8)   ! number of each kindof gene
integer,allocatable  ::wkindof(:,:) ! which are those for each kindof


!real*8, public, allocatable  :: didpol(:,:) !this one stores the diffusion differentials of the polarising genes (to be used only in single mode,
!                                            !when there is polarisation)

type, public :: genes

  !of genes on genes
  real*8, allocatable  ::t(:)      ! interaction strenghts of the genes regulating the transcription of gene i
  integer              ::nww       ! number of reactions catalyzed by i
  real*8,allocatable   ::r(:,:)   ! which reactions i regulates, dim1 index of the reaction(as in nww), dim2(1) the pre form, dim2(2) the post reaction, dim2(3) the w strength
  real*8               ::diffu     ! extracellular diffusivity of gene proteins or protein forms
  real*8               ::idiffu    ! intracellular diffusivity of gene proteins or protein forms
  real*8               ::mu        ! degradation of gene proteins or protein forms
  real*8               ::mich      ! Constant of Michaelis-Menten HC 31-3-20
  integer              ::npre      ! the number of pre form that gives rise to this form 
  integer,allocatable  ::pre(:)    ! the pre forms of a gene: the forms that give rise to it, NOTICE the reverse reaction only occurs if the post of a gene is its pre
                                   ! if that does not occur it means that the reaction is effectivelly irreversible
  integer              ::npost     ! number of post forms of a gene
  integer,allocatable  ::post(:)   ! the forms that arises from this form
  real*8               ::kindof    ! 1 transcribed/translated protein or RNA without further modification; 
                                   ! 2 transcribed/translated protein or RNA without that can be modified postraductionally
                                   ! 3 comes from a previous form   
                                   ! 4 extracelular signal: NOTICE THE MUST NOT HAVE A PRE: since they get secreted immediately and thus is like an irreversible reaction
                                   ! 5 form that is intrinsically transported to the apical side of an epithelial cell (through kinesins)  !in these two cases, transport depends only
                                   ! 6 form that is intrinsically transported to the basal side of an epithelial cell (through dineins)    !on the absolute concentration on one node,
                                   ! 7 membrane-bound signal/receptor: like notch-delta ; it requires a pre since 7 is only the activated form of the receptor
                                   ! 8 bound receptor of extracellular signal, it has to have one (and only one) pre, and promote its own change to the 
                                                                         !pre form (unbounding of the ligand)
                                   ! notice:                                                                                               !so the %diffus probably should be smaller.
				   !          extracellular signals are automatically secreted and can only affect genes in a node through receptors
				   !

  character*40        ::label      !this is intended to be used as a label to identify and track the gene (trhough evolution for example) !>>Miquel8-8-14 !!>> HC 30-6-2020


  !of genes on cell behaviours or properties
  real*8, allocatable  ::e(:)    ! effect of a gene on each node parameter. These are ordered by number; after the number of the node parameter
                                  ! come each of the cell behaviours
                                  ! nparam_per_node+1 = cell growth
                                  ! nparam_per_node+2 = cell cycle increase, when 1 the cell can divide if it has the right size
                                  ! nparam_per_node+3 = apoptosis
                                  ! nparam_per_node+4 = this gene promotes secretion at this rate 
                                  ! nparam_per_node+5 = 1 means this gene is secretable 
                                  !                     there is no parameter for adh because that comes from the genes expressed there
                                  ! nparam_per_node+6 = repcel of the secreted node
                                  ! nparam_per_node+7 = da (as prportion of reqcel) of the secreted node
                                  ! nparam_per_node+8 = this gene is polarizing the cell [affecting pol vector of cel]
                                  ! nparam_per_node+9 = this gene is telling the cells that they should grow in the polarized direction with
                                  !                     a noise that is 1-(this value)
                                  ! nparam_per_node+10= this affects the cell property minsize_for_div, the number of nodes required for a cell
                                  !                     to divide
                                  ! nparam_per_node+11= dependence of the plane of division relative to chemical gradients over the Hertwig vector
                                  ! nparam_per_node+12= activates assymetric division
                                  ! nparam_per_node+13= activates epitelial-to-mesenchymal transition (EMT)
                                  ! nparam_per_node+14= ECM proteolisis (like apoptosis for ECM nodes)
                                  ! nparam_per_node+15= changes maxsize_for_div the max number of nodes that allows division
                                  ! nparam_per_node+16= noise is biased by polarization vector
                                  ! nparam_per_node+17= polarization creates ellipsoid deformation of cell
end type
                                             !these two variables are to make things faster since it can be done only once
integer, public, allocatable :: npag(:)      !number of genes that have an effect in that cellular parameter
integer, public, allocatable :: whonpag(:,:) !list of the indexes of the genes affecting each cellular parameter

type(genes),public,allocatable :: gen(:)


! derived matrices that do not change and are just for making the code run faster
real*8 , public, allocatable :: nw(:),w(:,:)
integer, public, allocatable :: nwpre(:,:)    ! number of forms catalyzing the reaction pre j of i 
real*8 , public, allocatable :: wpre(:,:,:,:) ! dim1 the form i, dim2 index of the reactions leading to i from its pres (as in npre)
                                               ! dim3(1) the list of enzymes, dim4(1) their gene indices dim4(2) their ws  
integer, public, allocatable :: nwpost(:,:)   ! number of forms catalyzing the reaction post j of i 
real*8 , public, allocatable :: wpost(:,:,:,:)! dim1 the form i, dim2 index of the reactions leading to i from its posts (as in npost)
                                               ! dim3(1) the list of enzymes, dim4(1) their gene indices dim4(2) their ws

real*8 , private :: nbh
real*8 , public, allocatable :: geoch(:) !!>>HC 18-9-2020 This is to store total gene expression, used in nexe for filtering no changes in gene expression.

contains

!*************************************************************************************************
subroutine initiate_gene  ! allocates the matrices and fills npag and whonpag
  integer i

  if (allocated(gex)) deallocate(gex)
  allocate(gex(nda,ng))
  gex=0.0d0

  if (allocated(agex)) deallocate(agex)
  allocate(agex(nda,ng))
  agex=0

  if (allocated(gen)) deallocate(gen)
  allocate(gen(ng))

  do i=1,ng

    !w
    if (allocated(gen(i)%t)) deallocate(gen(i)%t)
    allocate(gen(i)%t(ng))
    gen(i)%t=0.0d0

    gen(i)%nww=0

    !ww 
    if (allocated(gen(i)%r)) deallocate(gen(i)%r)
    allocate(gen(i)%r(ng*ng,3))
    gen(i)%r=0.0d0


    !wa
    if (allocated(gen(i)%e)) deallocate(gen(i)%e)
    allocate(gen(i)%e(nga))
    gen(i)%e=0.0d0

  end do

  if (allocated(npag)) deallocate(npag)
  allocate(npag(nga))
  npag=0

  if (allocated(whonpag)) deallocate(whonpag)  
  allocate(whonpag(nga,ng))
  whonpag=0

  if (allocated(whonpag)) deallocate(whonpag)  
  allocate(whonpag(nga,ng))
  whonpag=0


  gen(:)%npre=0
  gen(:)%npost=0
  gen(:)%diffu=0.0d0
  gen(:)%idiffu=0.0d0
  gen(:)%mich = 1.0d0 !!!!! Michaelis Menten constant initialize HC 3-31-20
  gen(:)%kindof=1
  gen(:)%mu=0.0d0

  call update_npag
  
end subroutine

!*************************************************************************************************

subroutine update_npag   ! fils the npag and whonpag once wa is filled: this needs to be done only once per simulation (if there is no mutation)
  integer i,j,k

  npag=0
  whonpag=0
  do i=1,nga
    npag(i)=0
    do j=1,ng
      if (gen(j)%e(i)/=0.0d0) then
        npag(i)=npag(i)+1
        whonpag(i,npag(i))=j
      end if
    end do
  end do

  nkindof=0   !we make a matrix recording how many of each kindof gene there is and which are those
  if (allocated(wkindof)) deallocate(wkindof)
  allocate(wkindof(8,ng))
  do i=1,ng
    if (gen(i)%kindof==0) gen(i)%kindof=1  !every gene should have a kindof, if not we give the 1 value, for transcripts
    nkindof(int(gen(i)%kindof))=nkindof(int(gen(i)%kindof))+1   ! >>> Is 7-6-14   
    wkindof(int(gen(i)%kindof),nkindof(int(gen(i)%kindof)))=i   ! >>> Is 7-6-14 
  end do

  if (allocated(wpre)) deallocate(wpre)! number of non-zero w elements in each gene 
  k=maxval(gen(:)%npre)
  allocate(wpre(ng,k,ng,2))
  if (allocated(nwpre)) deallocate(nwpre)! number of non-zero w elements in each gene 
  allocate(nwpre(ng,k))

  if (allocated(wpost)) deallocate(wpost)! number of non-zero w elements in each gene 
  k=maxval(gen(:)%npost)
  allocate(wpost(ng,k,ng,2))
  if (allocated(nwpost)) deallocate(nwpost)! number of non-zero w elements in each gene 
  allocate(nwpost(ng,k))

  nwpre=0
  do i=1,ng
    do j=1,gen(i)%nww
rt:   do k=1,ng
        do kk=1,gen(k)%npre
          if (gen(i)%r(j,1)==gen(k)%pre(kk).and.gen(i)%r(j,2)==k) then ! these are the reaction from the pre to k catalyzed by i
            nwpre(k,kk)=nwpre(k,kk)+1                ! this is the number of enzymes catalyzing the reaction from k to pre(kk)
            wpre(k,kk,nwpre(k,kk),1)=i               ! this is the enzyme catalyzing that reaction
            wpre(k,kk,nwpre(k,kk),2)=gen(i)%r(j,3)  ! this is the activity 
          end if
        end do
      end do rt
    end do
  end do

  nwpost=0
  do i=1,ng
    do j=1,gen(i)%nww
rtd:  do k=1,ng
        do kk=1,gen(k)%npost
          if (gen(i)%r(j,1)==k.and.gen(i)%r(j,2)==gen(k)%post(kk)) then  ! those are the reactions from k to its post catalyzed by i
            nwpost(k,kk)=nwpost(k,kk)+1              ! this is the number of enzymes catalyzing the reaction from k to post(kk)
            wpost(k,kk,nwpost(k,kk),1)=i             ! this is the enzyme catalyzing that reaction
            wpost(k,kk,nwpost(k,kk),2)=gen(i)%r(j,3)! this is the activity
          end if
        end do
      end do rtd
    end do
  end do

  ! make the w and nw matrices
  if (allocated(nw)) deallocate(nw)! number of non-zero w elements un each gene
  allocate(nw(ng))  
  if (allocated(w)) deallocate(w)! which are those (indices)
  allocate(w(ng,ng))

  nw=0
  w=0

  do i=1,ng
    do j=1,ng
      if (gen(i)%t(j)/=0.0d0) then
        nw(i)=nw(i)+1
        w(i,int(nw(i)))=j  ! >>> Is 7-6-14
      end if
    end do
  end do


end subroutine update_npag

!***********************************************************************************************************************

subroutine gene_stuff ! gene interactions and gene diffusion


  call genestep

  if (ffu(5)/=2.and.ffu(12)/=2) then
    agex(:nd,1:ng) = gex(:nd,1:ng) + delta*dgex(:nd,1:ng) ;        
  else
    call rungekutta4gex
  end if

  do i=1,nd
    do j=1,ng
       if (abs(agex(i,j))<epsilod) gex(i,j)=0.0d0
    end do
  end do

  do i=1,nd        
    if (checkn(i)==0) then ! disruption of lateral signalling
      do jj=1,nkindof(7)
        k=wkindof(7,jj)
        if (gen(k)%npre>0) then ! >>> Is 9-5-14 
          kkk=gen(k)%pre(1)   !each notch molecule should have one and only one pre to make anything
          agex(i,kkk)=agex(i,kkk)+gex(i,k)! we lose all the bound forms !!!! RZ 4-3-14
          agex(i,k)=0.0d0     ! we lose all the bound forms
        end if  ! >>> Is 9-5-14
         ! mind that we do not consider partial disruptions here
      end do
    end if
  end do
end subroutine gene_stuff

!***********************************************************************************************************************
subroutine genestep
integer :: i, k, kk ! RZ 17-11-14 ! added k, kk
integer :: ii1,i1,ii2,i2,ii3,i3,iii1,iii2,iii3
integer :: kki,ivv
real*8  :: sm,smd,smp
integer :: ie,tipi,celi
real*8  :: ix,iy,iz
real*8  :: gext(ng)!,agex(nda,ng)
real*8  :: dist,nbh,udist,udistn
integer :: nv,nbo
real*8  :: did(ng)
!integer :: checkn(nd)
!integer :: c_nneigh(nd)
real*8  ::krev,kbound

       if (allocated(dgex)) deallocate(dgex)
       allocate(dgex(nd,ng))
       dgex=0.0d0
       if (allocated(checkn)) deallocate(checkn)
       allocate(checkn(nd))
       checkn=0.0d0

!      rdiffmax=2*maxval(node(:nd)%add)*dmax !now it's in the iteracio subroutine !>>Miquel27-2-14
      do i=1,nd
        did=0.0d0
        tipi=node(i)%tipus
        celi=node(i)%icel !>>Miquel31-3-14
        !if (rdiffmax<2*node(i)%add)   then ;print *,"RRRRRRRRRRRRR"; nbh=2*node(i)%add   ; else ; nbh=rdiffmax ; end if  ! this may make nbh different for different nodes 
        !!! and then non-symmetric diffusion 
        nbh=rdiffmax          
        nv=0
        ix=node(i)%x     ; iy=node(i)%y     ; iz=node(i)%z   
        ii1=nint(iz*urv) ; ii2=nint(iy*urv) ; ii3=nint(ix*urv)
        ivv=node(i)%altre

        nbo=nint(nbh*urv)
        !nneigh(i)=0 
        checkn(i)=0 ! RZ 4-3-14
        gext=0d0  !>>Miquel22-8-14
        ! ************************* STANDARD MINIMAL LOOP (1,2,3,4) **********************************
        ! DIFFUSION
        ! DIFFUSION WITH ALTRE                                                       !!>> TT 21-11-25
        if(tipi < 3) then                                                            !!>> TT 21-11-25
          ie = ivv                                                                   !!>> TT 21-11-25
          dist = sqrt((ix-node(ie)%x)**2+(iy-node(ie)%y)**2+(iz-node(ie)%z)**2)      !!>> TT 21-11-25
          udist=1.0d0/(1d0+dist)                                                     !!>> TT 21-11-25
          do jj=1,3                                                                  !!>> TT 21-11-25
            if (nkindof(jj)>0) then                                                  !!>> TT 21-11-25
              do jjj=1,nkindof(jj)                                                   !!>> TT 21-11-25
                j=wkindof(jj,jjj)                                                    !!>> TT 21-11-25
                did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist                !!>> TT 21-11-25
              end do                                                                 !!>> TT 21-11-25
            end if                                                                   !!>> TT 21-11-25
          end do                                                                     !!>> TT 21-11-25
          ! kindof 4                                                                 !!>> TT 21-11-25
          if (nkindof(4)>0) then                                                     !!>> TT 21-11-25
!            if (dist<nbh) then                                                      !!>> HC 26-11-2021 we do not check the diffusion limit for one cell and its altre
              do jjj=1,nkindof(4) !!>> TT 21-11-25 kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really) 
                j=wkindof(4,jjj)                                                     !!>> TT 21-11-25
                did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist                !!>> TT 21-11-25
              end do                                                                 !!>> TT 21-11-25
!            end if                                                                  !!>> HC 26-11-2021 
          end if                                                                     !!>> TT 21-11-25
          ! kindof 5                                                                 !!>> TT 21-11-25
          if (nkindof(5)>0) then                                                     !!>> TT 21-11-25
            do jjj=1,nkindof(5)                                                      !!>> TT 21-11-25 active apical-basal diffusion due to microtubule transport
              j=wkindof(5,jjj)                                                       !!>> TT 21-11-25
              if (tipi==2) then                                                      !!>> TT 21-11-25
                did(j)=did(j)-gen(j)%diffu*gex(i,j)*udist  !!>> TT 21-11-25 loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
              else if(tipi==1)then                         !!>> TT 21-11-25
                did(j)=did(j)+gen(j)%diffu*gex(ie,j)*udist !!>> TT 21-11-25 gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
              end if                                                                 !!>> TT 21-11-25
            end do                                                                   !!>> TT 21-11-25
          end if                                                                     !!>> TT 21-11-25
          if (nkindof(6)>0) then                                                     !!>> TT 21-11-25
            do jjj=1,nkindof(6)                                                      !!>> TT 21-11-25 active apical-basal diffusion due to microtubule transport
              j=wkindof(6,jjj)                                                       !!>> TT 21-11-25
              if(tipi==2)then                                                        !!>> TT 21-11-25
                did(j)=did(j)-gen(j)%diffu*gex(i,j)*udist  !!>> TT 21-11-25 loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
              else if(tipi==1)then                         !!>> TT 21-11-25
                did(j)=did(j)+gen(j)%diffu*gex(ie,j)*udist !!>> TT 21-11-25 gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
              end if                                                                 !!>> TT 21-11-25
            end do                                                                   !!>> TT 21-11-25
          end if                                                                     !!>> TT 21-11-25
        end if                                                                       !!>> TT 21-11-25
        ! DIFFUSION WITH NEIGHBORS  !!>> TT 21-11-25
        do i1=1,nneigh(i)
          ie=neigh(i,i1)
          dist=dneigh(i,i1)
          udist=1.0d0/(1d0+dist) ! >>> Is 25-5-14
          ! now in case we have kindof 5,6 and 7: NOTICE that these are membrane molecules so they do not diffuse within the cell
          !     then in order to reach the membrane they have to come through a previous form of a lower kindof 
          do jj=1,3
            if (nkindof(jj)>0) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                if(tipi<3)then
                  do jjj=1,nkindof(jj)
                    j=wkindof(jj,jjj) 
                    did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                  end do
                  cycle
                end if
              elseif (celi==node(ie)%icel) then
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  do jjj=1,nkindof(jj)
                    j=wkindof(jj,jjj) !; print*,jj,"j diffu",j,"gex",gex(i,j)
                    did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                  end do
                end if
              end if
            end if
          end do

          if (nkindof(4)>0) then
            if (dist<nbh) then
              do jjj=1,nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
                j=wkindof(4,jjj)
                if ((node(ie)%tipus==1 .and. tipi==2).or.(node(ie)%tipus==2 .and. tipi==1)) then 
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))/(dist+1) no diffusion across epitelia for signaling molecules
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist 
                else 
                  did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist !;print*,"i,ie",i,ie,"udist",udist
                end if
              end do
            end if
          end if
  
          ! kindof 5
          if (nkindof(5)>0) then
            if (tipi<3) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                do jjj=1,nkindof(5)  !active apical-basal diffusion due to microtubule transport
                  j=wkindof(5,jjj)
                  if (tipi==2) then
                    did(j)=did(j)-gen(j)%diffu*gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  else if(tipi==1)then
                    did(j)=did(j)+gen(j)%diffu*gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  end if
                end do
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  if(tipi==node(ie)%tipus)then
                    do jjj=1,nkindof(5)
                      j=wkindof(5,jjj)
                      did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (nkindof(6)>0) then
            if (tipi<3) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                do jjj=1,nkindof(6)  !active apical-basal diffusion due to microtubule transport
                  j=wkindof(6,jjj)
                  if(tipi==2)then
                    did(j)=did(j)-gen(j)%diffu*gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  else if(tipi==1)then
                    did(j)=did(j)+gen(j)%diffu*gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  end if
                end do
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
                  if(tipi==node(ie)%tipus)then
                    do jjj=1,nkindof(6)
                      j=wkindof(6,jjj)
                      did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (nkindof(7)>0) then !membrane-bound notch-delta kind of interaction >>> Is 9-2-14 
            if(tipi<4) then ! there is no reason to exclude mesenchymal cells RZ 4-3-14
              if (node(ie)%icel.ne.celi) then   ! ONLY from different cells
                if(tipi==node(ie)%tipus)then !>>Miquel2-10-14
                  if(dist.le.(node(i)%add+node(ie)%add)) then     !new kinetics !>>Miquel22-8-14
                    udistn=(node(i)%add+node(ie)%add)/(node(i)%add+node(ie)%add+dist) !notice this is the inverse Is
                    checkn(i)=1 !nodes close enough for binding
                    do jj=1,nkindof(7)
                      j=wkindof(7,jj)
                      if (gen(j)%npost>0) then !the dissociation reaction
                        do kki=1,gen(j)%npost   
                          !smd=0.0d0
                          kkk=gen(j)%post(kki) !the post form (either ligand or receptor)
                          do k=1,nwpost(j,kki)
                            kk=gen(int(wpost(j,kki,k,1)))%post(1) !the other bound form  ! RZ 17-11-14 INT
                            krev= wpost(j,kki,k,2) !dissociation constant
                          end do                           
                          b=krev*gex(i,j)*gex(ie,kk)   !loss of bound form by dissociation
                          gext(kkk) = gext(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
                          gext(j)   = gext(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)
                        end do
                      end if
                      if (gen(j)%npre>0) then !the binding reaction
                        do kki=1,gen(j)%npre   
                          kkk=gen(j)%pre(kki) !the pre form (either ligand or receptor)
                          do k=1,nwpre(j,kki)
                            kk=wpre(j,kki,k,1) !catalist (either ligand or receptor) and also the one
                            kbound=wpre(j,kki,k,2) !binding constant
                            smd=gex(i,kkk)*gex(ie,kk) !productory of the concentrations of ligand and receptor (kk and kkk)
                          end do
                          b=kbound*smd         ! gain of bound form by binding of receptor and ligand
                          gext(kkk) = gext(kkk) - b !this is the loss of dissociated form due to binding
                          gext(j)   = gext(j)   + b !this is the gain of bound form
                        end do
                      end if
                    end do
                  end if
                end if
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
                  if(tipi==node(ie)%tipus)then
                    do jjj=1,nkindof(7)
                      j=wkindof(7,jjj)
                      did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (nkindof(8)>0) then  ! receptors of extracellular diffusible molecules ONLY THE ACTIVE FORM, the inactive is kindof2 or 3
            if (celi==node(ie)%icel) then
              if (tipi==node(ie)%tipus) then ! the active form diffuses only in the membrane (the inactive everywhere)
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  do jjj=1,nkindof(8)
                    j=wkindof(8,jjj)
                    did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist
                  end do
                end if
              end if
            end if
          end if
        end do




        !LIMIT TO DIFFUSION, so that diffusion does not get to infinite distance with small amounts
        !do j=1,ng
        !  if (abs(did(j))<ldi) did(j)=0.0d0
        !end do

        !REACTION 
        if (node(i)%marge==0.and.node(i)%tipus<4) then ! IN THE NUCLEUS (THERE WILL BE TRANSCRIPTION)
          do jjj=1,nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=wkindof(1,jjj)     
            a=gex(i,j)
            sm=0.0d0
            do k=1,int(nw(j))
              kk=w(j,k)
              sm = sm + gen(j)%t(kk)*gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            gext(j) = sm/(gen(j)%mich+sm) - gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13 !!!!! Michaelis Menten constant can be different in each gene HC 31-3-20
          end do 
          do jjj=1,nkindof(4)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=wkindof(4,jjj)     
            a=gex(i,j)
            sm=0.0d0
            do k=1,int(nw(j))
              kk=w(j,k)
              sm = sm + gen(j)%t(kk)*gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            gext(j) = gext(j) + sm/(gen(j)%mich+sm) ! - gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13 !!!!! Michaelis Menten constant can be different in each gene HC 31-3-20 !!! Diffusion must be added later
          end do 
          do jjj=1,nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
            j=wkindof(2,jjj)
            a=gex(i,j)
            sm=0.0d0
            do k=1,int(nw(j))
              kk=w(j,k)
              sm = sm + gen(j)%t(kk)*gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            if (gen(j)%npost>0) then 
              do kki=1,gen(j)%npost   
                kkk=gen(j)%post(kki); if(gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                smd=0.0d0
                do k=1,nwpost(j,kki)
                  kk=int(wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                  smd=smd + wpost(j,kki,k,2)*gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                end do
                if (smd<0.0) smd=0.0   ! Is 9-2-14
                b=smd*a/(gen(j)%mich+a)   ! Michaelis-Menten !!Made variable between genes HC 31-3-20
                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
              end do
              gext(j) = gext(j) + sm/(gen(j)%mich+sm) - gen(j)%mu*a + did(j) !>>> Is 9-2-14 !!Made variable between genes HC 31-3-20
            else
              gext(j) = sm/(gen(j)%mich+sm) - gen(j)%mu*a + did(j)  !>>> Miquel20-12-13 !!Made variable between genes HC 31-3-20
            end if
          end do
        else  !****************** NOT IN THE NUCLEUS (NO TRANSCRIPTION)

          do jjj=1,nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=wkindof(1,jjj)  
            gext(j) = - gen(j)%mu*gex(i,j) + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          end do 
          do jjj=1,nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
            j=wkindof(2,jjj)     
            a=gex(i,j)
            if (gen(j)%npost>0) then !
              do kki=1,gen(j)%npost
                kkk=gen(j)%post(kki); if(gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                smd=0.0d0
                do k=1,nwpost(j,kki)
                  kk=int(wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                  smd=smd + wpost(j,kki,k,2)*gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                end do                           
                if (smd<0.0) smd=0.0   ! Is 9-2-14
                b=smd*a/(gen(j)%mich+a)   ! Michaelis-Menten !!Made variable between genes HC 31-3-20
                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
              end do
              gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
            else
              gext(j) = - gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
            end if
          end do
        end if
        
        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        do jjj=1,nkindof(3)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
          j=wkindof(3,jjj)
          a=gex(i,j)
          if (gen(j)%npost>0) then 
            do kki=1,gen(j)%npost
              kkk=gen(j)%post(kki); if(gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
              smd=0.0d0
              do k=1,nwpost(j,kki)
                kk=int(wpost(j,kki,k,1))  ! RZ 17-11-14 INT
                smd=smd + wpost(j,kki,k,2)*gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
              end do                           
              if (smd<0.0) smd=0.0   ! Is 9-2-14
              b=smd*a/(gen(j)%mich+a)   ! Michaelis-Menten !!Made variable between genes HC 31-3-20
              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
            end do
            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
          end if
          gext(j) = gext(j)  - gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do 

        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        do jjj=1,nkindof(4)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
          j=wkindof(4,jjj)     
          a=gex(i,j)
          if (gen(j)%npost>0) then 
            do kki=1,gen(j)%npost   
              kkk=gen(j)%post(kki); if(gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
              smd=0.0d0
              do k=1,nwpost(j,kki)
                kk=int(wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                smd=smd + wpost(j,kki,k,2)*gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
              end do                           
              if (smd<0.0) smd=0.0   ! Is 9-2-14
              b=smd*a/(gen(j)%mich+a)   ! Michaelis-Menten !!Made variable between genes HC 31-3-20
              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
            end do
            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
          end if

          ! NOTICE: growth factors do not return any thing to pre forms, there is not back reaction, because they get secreted as they get produced

          gext(j) = gext(j)  - gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do 

        do jj=5,7 !********PRODUCTION OF OTHER SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        !  if (jj==7) cycle   !>>> Is 1-3-14 ! outcommented RZ 4-3-14
          if (jj==7) then ! >>> RZ 4-3-14
            do jjj=1,nkindof(jj)
              j=wkindof(jj,jjj)     
              a=gex(i,j)
              gext(j) = gext(j)  - gen(j)%mu*a + did(j)  ! Is 9-2-14
            end do
          else   ! <<< RZ 4-3-14
            do jjj=1,nkindof(jj)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
              j=wkindof(jj,jjj)     
              a=gex(i,j)
              if (gen(j)%npost>0) then 
                do kki=1,gen(j)%npost   
                  kkk=gen(j)%post(kki); if(gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                  smd=0.0d0
                  do k=1,nwpost(j,kki)
                    kk=int(wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                    smd=smd + wpost(j,kki,k,2)*gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                  end do                           
                  if (smd<0.0) smd=0.0   ! Is 9-2-14
                  b=smd*a/(gen(j)%mich+a)   ! Michaelis-Menten !!Made variable between genes HC 31-3-20
                  gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                  gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
                end do
                !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
              end if
              gext(j) = gext(j)  - gen(j)%mu*a + did(j)  ! Is 9-2-14
            end do 
          end if  ! RZ 4-3-14
        end do

        !bound receptors !>>Miquel18-8-14
        do jjj=1,nkindof(8)
          j=wkindof(8,jjj)
          if (gen(j)%npost<1) cycle !they should though ! Is 1-10-14
          a=gex(i,j)
          !the dissociation reaction
          krev=0d0
          kkk=gen(j)%post(1) !the post form (either ligand or receptor)
          iii=gen(j)%post(2) !the post form (either ligand or receptor)
          krev= wpost(j,1,1,2) !dissociation constant
          b=krev*a       !loss of bound form by dissociation
          gext(kkk) = gext(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
          gext(iii) = gext(iii) + b !this is the gain of ligand or receptor by dissociation of bound form

          gext(j)   = gext(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)

          !the binding reaction
          if (gen(j)%npre<1) cycle !they should though !>>Miquel1-10-14
          kbound=0d0
          kkk=gen(j)%pre(1) !the pre form (either ligand or receptor)
          iii=gen(j)%pre(2) !the pre form (either ligand or receptor)
          kbound=wpre(j,1,1,2) !binding constant
          smd=gex(i,iii)*gex(i,kkk) !productory of the concentrations of ligand and receptor
          b=kbound*smd         ! gain of bound form by binding of receptor and ligand

          gext(kkk) = gext(kkk) - b !this is the loss of dissociated form due to binding
          gext(iii) = gext(iii) - b !this is the loss of dissociated form due to binding

          gext(j)   = gext(j)   + b !this is the gain of bound form due to binding
          gext(j) = gext(j)  - gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do
         
        dgex(i,:)=gext  ! this is the increment in the gene this step
      end do

end subroutine

!**********************************************************************************************

subroutine rungekutta4gex

real*8 d,halfd,sixthd
real*8 ogex(nd,ng)
real*8 kux(nd,ng)
real*8 kdx(nd,ng)
real*8 ktx(nd,ng)
real*8 kqx(nd,ng)

d=delta
halfd=d*0.5d0
sixthd=d/6.0d0

ogex=gex(:nd,:ng)

!k1
kux=dgex(:nd,:ng) 

!k2
gex(:nd,:ng)=gex(:nd,:ng)+halfd*dgex(:nd,:ng)

call genestep

kdx=dgex(:nd,:ng) 

!k3
gex(:nd,:ng)=gex(:nd,:ng)+halfd*dgex(:nd,:ng)

call genestep

ktx=dgex(:nd,:ng) 

!k4 
gex(:nd,:ng)=gex(:nd,:ng)+d*dgex(:nd,:ng)

call genestep

kqx=dgex(:nd,:ng)

!final
agex(:nd,:ng)=ogex+sixthd*(kux+2*kdx+2*ktx+kqx)

end subroutine

end module genetic
