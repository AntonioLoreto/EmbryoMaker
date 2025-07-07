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
module mutation
 
use general
use genetic
use io      !!>>HC 6-10-2021

implicit none

real*8 ismurate         ! IS mutation rate
real*8 mag_ismurate     ! magnitude of the IS mutation rate in % of the original value
real*8 tmuratew_g       !   T mutation that change w: gain a new interaction
real*8 tmuratew_l       !   T mutation that change w: lose an existing interaction 
real*8 tmurateww        !   T mutation that changes ww
real*8 tmuratewa        !   T mutation that change wa
real*8 tmurateadh       !   T mutation that change kadh
real*8 mag_tmurate      ! magnitude of the T mutation rate in % of the original value
real*8 duprate          ! gene duplication rate (also of deletion)
real*8 new_form_mu_rate ! mutation rate to get a new form, or to lose an existing form
integer,allocatable    :: dupi(:)        ! 1 if the gene has been duplicated in this round, 0 otherwise
integer,allocatable    :: father_son(:)  ! which gene duplicated to produce which
integer,allocatable    :: todel(:)
integer dupcount        ! number of duplications from a given gen
integer ng_before_dup   
integer ang,ij,gnr !gnr is variable for do while loop for dupdels

integer it_has_mutated  ! 1 if there has been a mutation, 0 otherwise
real*8:: myran

real*8 :: b_max,b_min,w_max,w_min,d_max,d_min
real*8 :: m_max,m_min,a_max,a_min,e_max,e_min
real*8 :: ap_max,ap_min,fp_max,fp_min,em_max,em_min
real*8 :: ma_max,ma_min,dv_max,dv_min
real*8 :: mich_max, mich_min   !!>> HC 7-7-2020

integer ::free, empty, full !for T mutations of w per gene
integer, allocatable ::freespot(:), fullspot(:) !to store positions of free and full spots in w matrix

contains

! IMPORTANT:::::::
! IT ASSUMES THAT ALL THE FORMS DERIVED FROM A GENE CAN NOT HAVE A PRE COMING FROM A DIFFERENT GENE (IT TAKES SENSE)
! THE kindof CAN NOT BE CHANGED SO TO GET NEW GROWTH FACTORS, notchs or alike ONE OF THEM HAS TO BE DUPLCATED
! RECEPTORS OF GROWTH FACTORS CAN ARISE SIMPLY BY A NEW w GETTING TO A GENE FROM A GROWTH FACTOR, BUT IT HAS TO BE A GENE THAT IS A FORM (A POST
! OF SOME OTHER FORM (AND THIS FORM HAS TO REVERT TO ITS PRE BY its own w).
! WE DO NOT ENFORCE THE LATTER BUT THE FORMER: SO A NEW W FROM A GROWTH FACTOR CAN ONLY CONNECT TO A GENE THAT HAS A PRE FORM AND THAT IS KINDOF 8.

! FOR KINDOF 4, THAT IS GROWTH FACTORS, WE CAN NOT HAVE MUTATIONS IN WA BECAUSE WA IS AFFECTED INTRACELLULARLY AND WE CAN HAVE T-MUTATIONS IN W BECAUSE
! THAT SIMPLY MEANS MORE PRODUCTION OF IT (OR INHIBITION OF THAT PRODUCTION)

! FOR KINDOF 8, BOTH KINDOF 4 AND OTHER CAN AFFECT ITS TRANSITION FROM A KINDOF 3

subroutine muta(limit, inviable)
  real*8 intu,intd
  integer :: limit,nwa, inviable
  real*8, allocatable :: copkadh(:,:)
  real*8 :: probmod,xx !the modifier used to bias T mutations of W towards interaction loss, and the temp var to store the ratio between empty and filled w positions (for readability only)
  integer genadh 
  it_has_mutated=0
 ! call random_number(a)
 ! if(a>0.001)return		!one of every thousand mutates pfh
  
!!  print*,"muta:: start mutating..."
  
  ! Set limits for the values that should not be crossed
  b_max=5.d0*10000 ; b_min=0d0!0.0053d0/10 !max used to be 5.33d0*100
  w_max=1.0d0*100 ; w_min=0d0 !used to be 33.6434d0 w_max is absolute value
  d_max=2.d0*100 ; d_min=0d0!0.00021d0/10
  m_max=32d0 ; m_min=0d0!1d0
  a_max=5d0*10000 ; a_min=0d0!0.000375d0/10
  e_max=5d0*10000 ; e_min=0d0!0.083d0/10
  ap_max=(0.033d0/deltamax)*100 ; ap_min=0d0
  fp_max=75d0*100 ; fp_min=0d0!0.0075d0/10
  em_max=5d0*1000 ; em_min=0d0  !max used to be 1d0*100
  ma_max=(100./deltamax)*100 ; ma_min=0d0
  dv_max=5d0*10000 ; dv_min=0d0
  mich_max=3.0d0; mich_min=0.0d0 !!>> HC 7-7-2020
 
  if (allocated(dupi)) deallocate(dupi)
  if (allocated(father_son)) deallocate(father_son)
  if (allocated(todel)) deallocate(todel)
  allocate(dupi(ng*ng))       !since each gene can no duplicate more than ng times per mutation round
  allocate(father_son(ng*ng)) !since each gene can no duplicate more than ng times per mutation round
  allocate(todel(ng*ng*10))
  dupi=0
  father_son=0
  todel=0
  nwa=28
!call
  !mutations
  ! the p mutation is per site so we can just pass over and mutate as it goes

  ! IS MUTATIONS
444  do i=1,ng

    genadh=int(gen(i)%e(1)) !handy for the kadh arrays later >>> Renske !!>> HC 30-6-2020
    
     do j=1,ng
      if (gen(i)%t(j)/=0.0d0) then !modify existing interaction (will also modify if almost zero though...) !!>> HC 30-6-2020
        call random_number(a)    
        if (a<ismurate) then
          call random_number(a) 
!!print *,i,j,"IS W",gen(i)%t(j),gen(i)%t(j)+(1-2*a)*gen(i)%t(j)*mag_ismurate !!>> HC 30-6-2020
          gen(i)%t(j)=gen(i)%t(j)+(1-2*a)*gen(i)%t(j)*mag_ismurate  !W !!>> HC 30-6-2020
          it_has_mutated=1
        end if
      end if
    end do

    do j=1,gen(i)%nww !mutate reactions catalyzed by i (dim2->3 is the strength of the reaction)
      call random_number(a)    
      if (a<ismurate) then
        call random_number(a) 
!!print *,i,j,"IS WW",gen(i)%r(j,3),gen(i)%r(j,3)+(1-2*a)*gen(i)%r(j,3)*mag_ismurate  !!>> HC 30-6-2020
        gen(i)%r(j,3)=gen(i)%r(j,3)+(1-2*a)*gen(i)%r(j,3)*mag_ismurate  !WW  !!>> HC 30-6-2020
        it_has_mutated=1
      end if
    end do
    call random_number(a) 
    if (a<ismurate) then !mutate gene product degradation rate
      call random_number(a)
      if (gen(i)%mu==0.0d0) then
        if(i/=5) gen(i)%mu=m_max*a!10**(log10(m_min)+a*log10(m_max/m_min))!a*mag_tmurate    !MU !pfh***11-09-15 What is special about gene 5? >>Renske
      else
        gen(i)%mu=gen(i)%mu+(1-2*a)*gen(i)%mu*mag_ismurate    !MU
      end if
      if (gen(i)%mu<0.0d0) gen(i)%mu=0.0d0
!!print *,i,"IS MU",gen(i)%mu,gen(i)%mu+(1-2*a)*gen(i)%mu*mag_ismurate
      it_has_mutated=1
    end if
    call random_number(a) 
    if (a<ismurate) then !mutate extracellular diffusivity of gene proteins or protein forms
      call random_number(a)
      if (gen(i)%diffu==0) then
        gen(i)%diffu=d_max*a!10**(log10(d_min)+a*log10(d_max/d_min))!a*mag_tmurate    !DIF				!pfh***11-09-15
      else
        gen(i)%diffu=gen(i)%diffu+(1-2*a)*gen(i)%diffu*mag_ismurate !DIF
      end if
      if (gen(i)%diffu<0.0d0) gen(i)%diffu=0.0d0
!!print *,i,"IS DIF",gen(i)%diffu,gen(i)%diffu+(1-2*a)*gen(i)%diffu*mag_ismurate !DIF 
      it_has_mutated=1
    end if
    call random_number(a) 
    if (a<ismurate) then !mutate the constant of Michaelis-Menten                                   !!>> HC 7-7-2020
      call random_number(a)                                                                         !!>> HC 7-7-2020
      if (gen(i)%mich==0) then                                                                      !!>> HC 7-7-2020
        gen(i)%mich=mich_max*a!10**(log10(d_min)+a*log10(d_max/d_min))!a*mag_tmurate    !DIF	    !!>> HC 7-7-2020	!pfh***11-09-15
      else                                                                                          !!>> HC 7-7-2020
        gen(i)%mich=gen(i)%mich+(1-2*a)*gen(i)%mich*mag_ismurate !DIF                               !!>> HC 7-7-2020
      end if                                                                                        !!>> HC 7-7-2020
      if (gen(i)%mich<0.0d0) gen(i)%mich=0.0d0                                                      !!>> HC 7-7-2020
!!print *,i,"IS DIF",gen(i)%diffu,gen(i)%diffu+(1-2*a)*gen(i)%diffu*mag_ismurate !DIF               !!>> HC 7-7-2020
      it_has_mutated=1                                                                              !!>> HC 7-7-2020
    end if                                                                                          !!>> HC 7-7-2020
    do j=2,nga  !mutates effect of a gene on each node parameter. 
      if (gen(i)%e(j)/=0.0) then ! IS mutation of wa  !!>> HC 30-6-2020
        call random_number(a) 
        if (a<ismurate) then
          call random_number(a)
!!print *,"gen",i,"wa",j,"IS wa",gen(i)%e(j),"new val.",gen(i)%e(j)+(1-2*a)*gen(i)%e(j)*mag_ismurate !!>> HC 30-6-2020
          gen(i)%e(j)=gen(i)%e(j)+(1-2*a)*gen(i)%e(j)*mag_ismurate !WA !!>> HC 30-6-2020
          if (gen(i)%e(j)<0.0d0) gen(i)%e(j)=0.0d0 !!>> HC 30-6-2020
          it_has_mutated=1
        end if
      end if
    end do
    if (genadh/=0) then ! the wa=1, can not be changed since it is the ntipusadh index of the gene if it is adhesive
      do j=1,ntipusadh
        call random_number(a)
        if (a<ismurate) then
          if (kadh(genadh,j)/=0.0d0) then                   !KADH adhesions between adhesion molecules
            call random_number(a)
!!print *,i,j,"IS kadh",kadh(gen(i)%e(1),j),kadh(gen(i)%e(1),j)+(1-2*a)*kadh(gen(i)%e(1),j)*mag_ismurate !!>> HC 30-6-2020
            kadh(genadh,j)=kadh(genadh,j)+(1-2*a)*kadh(genadh,j)*mag_ismurate
            it_has_mutated=1
          end if
        end if
      end do
    end if
  end do !end of IS mutations loop over all the genes

!!T MUTATIONS
  do i=1,ng
    !!!!instead of going per gene, we will look per interaction. You are more likely to lose an existing interaction than to gain a new one (seems more realistic) >> renske 08-05-18
!     do j=1, ng                      
!       call random_number(a) ! T mutation of w: gain interaction
!       if (a<tmuratew_g .and. gen(i)%t(j)==0.0d0) then !gain interaction !!>> HC 30-6-2020
!         it_has_mutated=1
!         call random_number(a)      
!         call random_number(b)
!         if(b>0.5)then; b=-1; else; b=1; endif
! !!print *,i,j,"++T",b*a*w_max!(1.0d0-2.0d0*a)*mag_tmurate
!         gen(i)%t(j)=b*a*w_max!(1.0d0-2.0d0*a)*mag_tmurate		!pfh 11-09-15 !!>> HC 30-6-2020
!         if (gen(i)%e(j)<0.0d0) gen(i)%e(j)=0.0d0 !why is this here?  !!>> HC 30-6-2020
!       else if (a<tmuratew_l .and. gen(i)%t(j)/=0.0d0) then	!lose interaction !!>> HC 30-6-2020
! !!print *,i,j,"--T"
!         gen(i)%t(j)=0.0d0 !!>> HC 30-6-2020
!         it_has_mutated=1  
!       end if    
!     end do 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! per-gene fixed probability of a gaining an interaction in the upstream region >> Renske 13-12-18
    ! per-interaction probability of losing that interaction >>> Renske 13-12-18
    
    empty=0
    full=0
    if (allocated(freespot)) deallocate(freespot)
    allocate(freespot(ng))
    !if (allocated(fullspot)) deallocate(fullspot)
    !allocate(fullspot(ng))
    do j=1, ng 
      if (gen(i)%t(j)==0.0d0) then !where are the empty spots? !!>> HC 30-6-2020
        empty=empty+1
        freespot(empty)=j
      else
        !full=full+1
        !fullspot(full)=j
        !now, do we delete this TFBS? 
        call random_number(a)
        if(a<tmuratew_l) then 
          it_has_mutated=1
          gen(i)%t(j)=0.0d0 !!>> HC 30-6-2020
        endif  
      endif
    end do 
    !1-tmuratew_g is the prob of nothing happening. Whether you gain or lose depends on the ratio between the two.

    ! gain max 1 interaction per gene. If more positions are open, the probability is higher, to max tmuratew_g
    call random_number(a)
    probmod=real(empty,8)/real(ng,8)
    if(a<probmod*tmuratew_g) then
      it_has_mutated=1
      call random_number(b)
      call random_number(c) 
      j=int(c*empty)+1
      if(b>0.5)then; b=-1; else; b=1; endif
      call random_number(c)
      gen(i)%t(freespot(j))=b*c*w_max  !(1.0d0-2.0d0*a)*mag_tmurate !!>> HC 30-6-2020
    endif
    !old per-gene probability of gaining or losing an interaction
    ! we want the probability of losing an interaction to be slightly higher than gaining one
!     xx=real(empty,8)/real(ng,8)
!     probmod=xx**tmuratew_l !The larger tmuratew_l, the larger the bias to loss (>1). Set probmod equal to xx if you want equal prob of gain and loss, or tmuratew_l equal to 1
!     call random_number(a)
!     if (a<tmuratew_g*probmod) then !you gain an interaction
!       it_has_mutated=1
!       call random_number(b)
!       call random_number(c) 
!       j=int(c*empty)+1
!       if(b>0.5)then; b=-1; else; b=1; endif
!       call random_number(c)
!       gen(i)%t(freespot(j))=b*c*w_max  !(1.0d0-2.0d0*a)*mag_tmurate !!>> HC 30-6-2020
!     elseif (a<tmuratew_g) then !you lose an interaction
!       it_has_mutated=1
!       call random_number(b)
!       j=int(b*full)+1      
!       gen(i)%t(fullspot(j))=0.0d0 !(1.0d0-2.0d0*a)*mag_tmurate !!>> HC 30-6-2020
!     endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    call random_number(a) ! T mutation of ww
    if (a<tmurateww.and.gen(i)%nww<ng) then 
      call random_number(a)
      k=int(a*ng)+1
      kk=-1
      do j=k,ng
        if (gen(j)%npost>0) then ; kk=j ; goto 51 ; end if
      end do
      do j=1,k-1
        if (gen(j)%npost>0) then ; kk=j ; goto 51 ; end if
      end do
51    if (kk>0) then !if there are no posttranslational forms, skip this (might want to check this sooner, before drawing a bunch of random nrs)
        call random_number(a)          
        k=int(gen(kk)%npost*a)+1
        do j=1,gen(i)%nww
          if (gen(i)%r(j,1)==kk.and.gen(i)%r(j,2)==k) then  ! i was already catalyzing the reaction from kk to k !!>> HC 30-6-2020
            call random_number(a) 
            gen(i)%r(j,3)=gen(i)%r(j,3)+(1-2*a)*gen(i)%r(j,3)*mag_ismurate  !WW so it is like a IS mutation !!>> HC 30-6-2020
!!print *,i,j,"T WW",gen(i)%r(j,3),gen(i)%r(j,3)+(1-2*a)*gen(i)%r(j,3)*mag_ismurate !!>> HC 30-6-2020
            goto 344
          end if
        end do      
        gen(i)%nww=gen(i)%nww+1
        gen(i)%r(gen(i)%nww,1)=kk !!>> HC 30-6-2020
        gen(i)%r(gen(i)%nww,2)=k  !!>> HC 30-6-2020
        call random_number(a)    
        gen(i)%r(gen(i)%nww,3)=(1.0d0-2.0d0*a)*mag_tmurate  !!>> HC 30-6-2020
!!print *,i,gen(i)%nww,kk,k,"T WW"
344     continue
      end if !posttranslational forms
    end if 
    call random_number(a)  !Delete a posttranslational interaction    
    if (a<tmurateww) then  ! delete T
      if (gen(i)%nww>0) then
        call random_number(a)          
        k=int(a*ng)+1
        do j=k,gen(i)%nww-1
          gen(i)%r(j,:)=gen(i)%r(j+1,:) !!>> HC 30-6-2020
        end do
        gen(i)%nww=gen(i)%nww-1
!!print *,i,gen(i)%nww,k,"-T WW"
      end if
      it_has_mutated=1 
    end if  
34  if (gen(i)%kindof<9) then  ! Is >>> 23-1-15 !mutate the effect of gene on nodes
      call random_number(a)
      if (a<tmuratewa) then
!        call random_number(a)
!        j=int(a*(nga-1))+2
        call random_number(a)
        if (a<0.5d0) then !add an interaction
!           print *,i,"T wa"
        

!          call random_number(a)  !wa% here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1/*/*/*/*/*//*////*/---down
!          print *,"gen",i,"wa",j,"T"!,a*mag_tmurate
         ! gen(i)%e(j)=(1.0d0-2.0d0*a)*mag_tmurate  !!>> HC 30-6-2020
 
!          do i=1,ng  					!pfh***11-09-15
446        call random_number(a)  ! RZ which one? 
           call random_number(b)  ! RZ how strong?
           call random_number(c)  ! RZ + or - ?
            if(c<0.5)then;c=1;else;c=-1;endif
          j=0;jj=0
          if(a<1d0/nwa)then ! new adh
            if(gen(i)%e(1)==0)then !!>> HC 30-6-2020
             ntipusadh=ntipusadh+1 
             gen(i)%e(1)=ntipusadh !!>> HC 30-6-2020
             if(allocated(copkadh))deallocate(copkadh)
             allocate(copkadh(ntipusadh,ntipusadh))
             copkadh(:ntipusadh-1,:ntipusadh-1)=kadh(:ntipusadh-1,:ntipusadh-1)
             deallocate(kadh)
             allocate(kadh(ntipusadh,ntipusadh))
             kadh=0.0d0
             kadh(:ntipusadh-1,:ntipusadh-1)=copkadh(:ntipusadh-1,:ntipusadh-1)
              
             call random_number(a)      
             jj=int(ntipusadh*a)+1 !pick random gene to adhere to?
             kadh(ntipusadh,jj)=b_max*c*b!c*10**(log10(b_min)+b*log10(b_max/b_min))
             kadh(jj,ntipusadh)=b_max*c*b!c*10**(log10(b_min)+b*log10(b_max/b_min))
             j=1
            else;goto 446 ! Why would you want to keep trying until you pass the if statement? To find an empty one?
            endif
          ! all remnant wa's according to the limits of their activation strengths
          elseif(a<2d0/nwa)then
            gen(i)%e(5)=e_max*c*b;j=5!c*10**(log10(e_min)+b*log10(e_max/e_min));j=5 !!>> HC 30-6-2020
          elseif(a<4d0/nwa)then
            gen(i)%e(int(a*nwa+4))=a_max*c*b;j=int(a*nwa+4)!c*10**(log10(a_min)+b*log10(a_max/a_min));j=int(a*nwa+4) ! RZ !!>> HC 30-6-2020
          elseif(a<12d0/nwa)then
            gen(i)%e(int(a*nwa+4))=b_max*c*b;j=int(a*nwa+4)!c*10**(log10(b_min)+b*log10(b_max/b_min));j=int(a*nwa+4) ! RZ !!>> HC 30-6-2020
          elseif(a<13d0/nwa)then
            gen(i)%e(16)=fp_max*c*b;j=16!c*10**(log10(fp_min)+b*log10(fp_max/fp_min));j=16 ! RZ !!>> HC 30-6-2020
          elseif(a<18d0/nwa)then
            gen(i)%e(int(a*nwa+7))=e_max*c*b;j=int(a*nwa+7)!c*10**(log10(e_min)+b*log10(e_max/e_min));j=int(a*nwa+7) ! RZ ! e.g.contraction 21 !!>> HC 30-6-2020
          elseif(a<21d0/nwa)then
            gen(i)%e(int(a*nwa+8))=b_max*c*b;j=int(a*nwa+8)!c*10**(log10(b_min)+b*log10(b_max/b_min));j=int(a*nwa+8) ! RZ *100 !!>> HC 30-6-2020
          elseif(a<22d0/nwa)then  ! growth
             a=b*dv_max!dv_min+10**(b*(log10(dv_max)))
             gen(i)%e(nparam_per_node+2)=a;j=nparam_per_node+2 !!>> HC 30-6-2020
             gen(i)%e(nparam_per_node+1)=a ! THIS IS IRRELEVANT FOR SINGLE_NODE, instead of dv !!>> HC 30-6-2020
            ! do k=1,ncels
            !  call random_number(c)
            !  cels(k)%fase=0d0 ! change this if you want cells to divide asynchronically
            !end do
             487 call random_number(b) ! if polarized
             if(b.gt.0.5)then ! 50% of the divisons shall be polarized.
              call random_number(b)
              gen(i)%e(nparam_per_node+9)=b*b_max;jj=nparam_per_node+9 !10**(log10(b_min)+b*log10(b_max/b_min)) !RZ !!>> HC 30-6-2020
            endif
          elseif(a<23d0/nwa)then
            a=b*ap_max!ap_min+10**(b*(log10(ap_max)))
            gen(i)%e(nparam_per_node+3)=a;j=nparam_per_node+3 !!>> HC 30-6-2020
          elseif(a<24d0/nwa)then
            a=ma_max*b!ma_min+10**(b*(log10(ma_max)))
            gen(i)%e(nparam_per_node+4)=a;j=nparam_per_node+4 !!>> HC 30-6-2020
            call random_number(b)
            gen(i)%e(nparam_per_node+6)=b;j=nparam_per_node+6 ! RZ !!>> HC 30-6-2020
            call random_number(b)
            gen(i)%e(nparam_per_node+7)=b;j=nparam_per_node+7 ! RZ !!>> HC 30-6-2020
          elseif(a<25d0/nwa)then
            gen(i)%e(nparam_per_node+8)=b_max*c*b;j=nparam_per_node+8!10**(log10(b_min)+b*log10(b_max/b_min));j=nparam_per_node+8!*10**(7+2*c) ! RZ ! abit arbitrary, but b_max seems ok !!>> HC 30-6-2020
          elseif(a<26d0/nwa)then
            a=b*em_max!em_min+10**(b*(log10(em_max)))
            gen(i)%e(nparam_per_node+13)=a;j=nparam_per_node+13 !!>> HC 30-6-2020
          elseif(a<27d0/nwa)then
            gen(i)%e(nparam_per_node+14)=m_max*b;j=nparam_per_node+14!10**(log10(m_min)+b*log10(m_max/m_min));j=nparam_per_node+14 ! RZ !!>> HC 30-6-2020
          else
            gen(i)%e(nparam_per_node+16)=a*1000;j=nparam_per_node+16 !arbitrary  !!>> HC 30-6-2020
          end if
!!print*,"gen",i,"wa",j,"other",jj,"+T wa"
           
      

        else  !	(a>0.5d0): delete an interaction	 
!          call random_number(a)
!          j=int(a*(nga-1))+2
           j=0;jj=0		!pfh***11-09-15

447        call random_number(a)  ! RZ which one? 
           call random_number(b)  ! RZ how strong?
           call random_number(c)  ! RZ + or - ?
            if(c<0.5)then;c=1;else;c=-1;endif

          if(a<1d0/nwa)then ! adh
            if(gen(i)%e(1)/=0)then !!>> HC 30-6-2020
              if(maxval(kadh(ntipusadh,:))>0)then
336             call random_number(a)      
                jj=int(ntipusadh*a)+1
                if(kadh(ntipusadh,jj)/=0)then
                  kadh(ntipusadh,jj)=0
                  kadh(jj,ntipusadh)=0
                  j=1
                else;goto 336;endif
              else;goto 447;endif      
            else;goto 447;endif 
           
          ! all remnant wa's according to the limits of their activation strengths
          elseif(a<2d0/nwa)then
            gen(i)%e(5)=0;j=5 !!>> HC 30-6-2020

          elseif(a<4d0/nwa)then
            gen(i)%e(int(a*nwa+4))=0;j=int(a*nwa+4) !!>> HC 30-6-2020
          elseif(a<12d0/nwa)then
            gen(i)%e(int(a*nwa+4))=0;j=int(a*nwa+4) !!>> HC 30-6-2020
          elseif(a<13d0/nwa)then
            gen(i)%e(16)=0;j=16                     !!>> HC 30-6-2020
          elseif(a<18d0/nwa)then
            gen(i)%e(int(a*nwa+7))=0;j=int(a*nwa+7) !!>> HC 30-6-2020
          elseif(a<21d0/nwa)then
            gen(i)%e(int(a*nwa+8))=0;j=int(a*nwa+8) !!>> HC 30-6-2020
          elseif(a<22d0/nwa)then  ! growth
             gen(i)%e(nparam_per_node+2)=0;j=nparam_per_node+1 !!>> HC 30-6-2020
             gen(i)%e(nparam_per_node+1)=0                     !!>> HC 30-6-2020
             490 call random_number(b) ! if polarized
             if(b.gt.0.5)then ! 50% of the divisons shall be polarized.
              call random_number(b)
              gen(i)%e(nparam_per_node+9)=0 ; jj=nparam_per_node+9 !!>> HC 30-6-2020
            endif
          elseif(a<23d0/nwa)then
            gen(i)%e(nparam_per_node+3)=0;j=nparam_per_node+3      !!>> HC 30-6-2020
          elseif(a<24d0/nwa)then
            gen(i)%e(nparam_per_node+4)=0;j=nparam_per_node+4      !!>> HC 30-6-2020
            call random_number(b)
            gen(i)%e(nparam_per_node+6)=0;j=nparam_per_node+6      !!>> HC 30-6-2020
            call random_number(b)
            gen(i)%e(nparam_per_node+7)=0;j=nparam_per_node+7      !!>> HC 30-6-2020
          elseif(a<25d0/nwa)then
            gen(i)%e(nparam_per_node+8)=0;j=nparam_per_node+8      !!>> HC 30-6-2020
          elseif(a<26d0/nwa)then
            gen(i)%e(nparam_per_node+13)=0;j=nparam_per_node+13    !!>> HC 30-6-2020
          elseif(a<27d0/nwa)then
            gen(i)%e(nparam_per_node+14)=0;j=nparam_per_node+14    !!>> HC 30-6-2020
          end if 

!!print *,"gen",i,"wa",j,"other",jj,"-T wa"
         
        end if       
        it_has_mutated=1
      end if
!      if (gen(i)%e(1)/=0) then !existing adh  !!>> HC 30-6-2020   !!>> HC 18-11-2020  This was commented out because the mutations in kadh
!        call random_number(a)                                     !!>> HC 18-11-2020  matrix have already been done above
!        if (a<tmurateadh) then                                    !!>> HC 18-11-2020  this chunk increased the probabilitu of mutations in
!          call random_number(a)                                   !!>> HC 18-11-2020  adhesion molecules.
!          j=int(ntipusadh*a)+1                                    !!>> HC 18-11-2020
!          call random_number(a)                                   !!>> HC 18-11-2020
!          if (a<0.5d0) then                                       !!>> HC 18-11-2020
!            call random_number(a) ! T mutation of kadh            !!>> HC 18-11-2020
!print *,i,j,"T kadh",a*mag_tmurate                                !!>> HC 18-11-2020
!            genadh=int(gen(i)%e(1)) !better to cast it >>> Renske !!>> HC 18-11-2020 !!>> HC 30-6-2020
!            kadh(genadh,j)=(1.0d0-2.0d0*a)*mag_tmurate            !!>> HC 18-11-2020
!          else                                                    !!>> HC 18-11-2020
!            call random_number(a) ! T mutation to take out a kadh !!>> HC 18-11-2020
!!print *,i,j,"-T kadh"                                            !!>> HC 18-11-2020
!            genadh=int(gen(i)%e(1))   !!>> HC 30-6-2020           !!>> HC 18-11-2020
!            kadh(genadh,j)=0.0d0                                  !!>> HC 18-11-2020
!          end if                                                  !!>> HC 18-11-2020
!          it_has_mutated=1                                        !!>> HC 18-11-2020
!        end if                                                    !!>> HC 18-11-2020
!      end if                                                      !!>> HC 18-11-2020
    end if

  end do !probably the end of the loop over all genes for T mutations

  !duplications and deletions modified from earlier code that did not seem to run >>Renske
  ! any gene in the old genome either gets duplicated or deleted, not both; and that only once.
  ! Note: ang is the original nr of genes before any duplications happened. (so we don't duplicate/delete new genes)
  ! we subtract 1 from ang for deletions because those reduce the nr of leftover original genes.
 ang=ng
 gnr=2 ! we don't duplicate or delete the first gene
 do while (gnr<=ang)
   if (ng<ang) exit  
   if (gen(gnr)%npre==0) then ! only primary forms get duplicated as such
     call random_number(a)
     if (a<duprate) then
       !call random_number(myran)
       !if(myran<0.5000) then 
         !print *,gnr,"del",ng
         if (ng>2) then
!!print *, "gene label is "//trim(gen(gnr)%label) 
           call deletion(gnr)
           ang=ang-1 !so that we don't start duplicating and deleting newly duplicated genes
           it_has_mutated=1
         end if
       else if (a<duprate*2) then
         !print *,gnr,"dup"
         call duplication(gnr)
         it_has_mutated=1
       end if
     !end if
   end if
   gnr=gnr+1
 end do

 
  ! a new form arises (or gets deleted): it has initially all the w's and was equal to its father but there is no 
  ! gene affecting it (no w leading to it), and it is a duplication of a non-primary form
!   ang=ng
!   do i=1,ang
!     call random_number(a)
!     if (a<new_form_mu_rate) then 
! print *,i,"add form",ng
!       call dupli(i)
!       
!       !the duplicate is simply the post of the original
!       call add_form(i)
! 
!       !we make that sombody is catalyzing this new form
!       !call random_number(a)        
!       !k=int(a*ng)+1
!       !if (gen(k)%nww<ng) then
!       !  gen(k)%nww=gen(k)%nww+1
!       !  gen(k)%r(gen(k)%nww,1)=gen(ng)%pre(1)  !!>> HC 30-6-2020
!       !  gen(k)%r(gen(k)%nww,2)=ng              !!>> HC 30-6-2020
!       !  call random_number(a)        
!       !  gen(k)%r(gen(k)%nww,3)=a*mag_tmurate   !!>> HC 30-6-2020
!       !end if
!       it_has_mutated=1
!     end if
!   end do
!   ang=ng
!   do i=1,ang
!     !DELETION OF FORMS IS ACTUALLY A GENETIC CHANGE IN THE STRUCTURE OF THE PROTEIN THAT MAKES THAT SOME PHOSPORYLATION OR SO CAN NOT TAKE PLACE ANY MORE
!     call random_number(a)
!     if (a<new_form_mu_rate) then 
!       if (i<ng) then
! print *,"deletion of form",i,"ng",ng
!         call deletion(i)
! print *,"form deleted; new ng is;",ng
!         it_has_mutated=1
!       end if
!     end if
!   end do


  if (allocated(npag)) deallocate(npag)
  allocate(npag(nga))
  npag=0

  if (allocated(whonpag)) deallocate(whonpag)  
  allocate(whonpag(nga,ng))
  whonpag=0

  !if (allocated(whonpag)) deallocate(whonpag)  !was here double
  !allocate(whonpag(nga,ng))
  !whonpag=0

!do i=1,ng
!  if (gen(i)%nww>0) then
!    do j=1,gen(i)%nww
!      print *,i,j,gen(i)%nww,gen(i)%r(j,:) !!>> HC 30-6-2020
!    end do 
!  else
!   print *,nww,0
!  end if
!end do

  call update_npag

!!print *,"mutation finished"

!if(it_has_mutated==0)goto 444   !I don't want this: a mutation prob determines if and how often I get mutations

!  remutate=remutate+1
!  idum=preidum+idum+ind+ind*(gene-1)+remutate
!  print*,"preidum",preidum,"idum",idum,"ind",ind,remutate
!  print*,"idum",idum
!  call random_seed(put=idum)
!  print*,"*****************/*/*/*/*/*/*777777777777777777777*/*/*/*/*"
!  goto 444
!endif
call limits(limit)

do i=1,ng
  do j=1,gen(i)%npre
    if (gen(i)%pre(j)==0) then
      print *,i,j,"aqui hi ha un zero a pre" 
      stop
    end if
  end do

  do j=1,gen(i)%npost
    if (gen(i)%post(j)==0) then
      print *,i,j,"aqui hi ha un zero a post" 
      stop
    end if
  end do
end do

inviable=0
if (ng<1) then
  print *,"inviable individual because there are no genes left"
  inviable=1
else
  do j=1,ng
    if (gen(j)%kindof<3) then
      inviable=0
      goto 37
    end if 
  end do
  print *,"no gene is transcribable so this is inviable"
  inviable=1
end if
37 continue

end subroutine muta

!***********************************************************************************************************'




subroutine suremuta(limit, inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue,rangfile)
!!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file
!!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
!!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
!!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
!!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
!!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
  implicit none
  integer :: limit,inviable,mind,tind
  real*8, allocatable :: copkadh(:,:)
  real*8 :: probmod,xx !the modifier used to bias T mutations of W towards interaction loss, and the temp var to store the ratio between empty and filled w positions (for readability only)
  integer genadh 
  integer :: props,nume,numt,numr,numk !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
  integer :: luckyprop, luckylink, luckyg  !!>> HC 24-11-2020 selected gene, interaction or property
  integer :: ich,jch,lch,ord,ord2  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
  integer,allocatable,dimension(:) :: plinki,plinkj !!>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
  integer,allocatable,dimension(:) :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM
  integer,allocatable,dimension(:) :: isadh !!>> HC 24-11-2020 if 1 the interaction is in the kadh matrix (TM)
  integer :: full,empty !!>> HC 24-11-2020 full an empty spaces for T mutations
  real*8, dimension(1:5) :: probs !!>> HC 24-11-2020 relative probabilities of different T mutations
  real*8, dimension(1:7) :: probis !!>> HC 24-11-2020 relative probabilities of different IS mutations
  real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
  integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  real*8 :: cum !!>> HC 24-11-2020 funny position to store cumulative Fi
  integer :: stilladi, stilladj,isad, mnadh, mxadh !!>> HC 24-11-2020 pointers related to T mutations in the Kadh matrix
  integer :: signo  !!>> HC 24-11-2020 sign of the mutation  
  real*8 :: inch, newval, anch !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
  real*8, dimension(1:5) :: bhc !!>>HC 10-9-2021 Patch to a problem with biased random numbers
  integer :: mutacode,mutageni,mutagenj  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
  real*8 :: prevalue, newvalue
  character*400 :: rangfile !!>> 6-10-2021 file where the ranges are stored
  real*8, dimension (1:5) :: min_glim, max_glim                                  !!>>HC 6-10-2021
  
  !!!!!!!!!!!!!!!!!!!!!!! LIMITS AND RANGES !!!!!!!!!!!!!!!!!!!!!!! 
                   
  limit=1                                               !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
  rembeh=0; max_elim=0.0d0; min_elim=0.0d0              !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
  max_glim=0.0d0; min_glim=0.0d0                        !!>> HC 6-10-2020 Same for other gene properties
  
  call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
  
  d_max=max_glim(1)    ; d_min=min_glim(1)                       !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020
  m_max=max_glim(2)    ; m_min=min_glim(2)                       !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
  mich_max=max_glim(3) ; mich_min=min_glim(3)                    !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
  w_max=max_glim(4)    ; w_min=min_glim(4)                       !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
  b_max=max_glim(5)    ; b_min=min_glim(5)                       !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  if (allocated(dupi)) deallocate(dupi)
  if (allocated(father_son)) deallocate(father_son)
  if (allocated(todel)) deallocate(todel)
  allocate(dupi(ng*ng))       !since each gene can no duplicate more than ng times per mutation round
  allocate(father_son(ng*ng)) !since each gene can no duplicate more than ng times per mutation round
  allocate(todel(ng*ng*10))
  dupi=0
  father_son=0
  todel=0
  empty=0; full=0             !!>> HC 24-11-2020
  inch=0.0d0; anch=0.0d0      !!>> HC 28-11-2020
  mutacode=0                  !!>> HC 16-9-2021  
  mutageni=0                  !!>> HC 29-11-2021
  mutagenj=0                  !!>> HC 29-11-2021
  prevalue=0.0d0              !!>> HC 29-11-2021
  newvalue=0.0d0              !!>> HC 29-11-2021
  call random_number(bhc)     !!>> HC 25-11-2020 THESE ARE BIASED I DO NOT KNOW WHY

  if (mind==1)then                                        !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!! IS MUTATIONS!!!!!!!!!!!!!!!!!!!!!
     nume=0; numt=0; numk=0; numr=0                       !!>> HC 25-11-2020 
     do ich=1,ng                                          !!>> HC 25-11-2020 
        do jch=1,ng                                       !!>> HC 25-11-2020 
           if(gen(ich)%t(jch)>0.0d0) numt=numt+1          !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
        enddo                                             !!>> HC 25-11-2020 
        do jch=1,gen(ich)%nww                             !!>> HC 25-11-2020 
           if(gen(ich)%r(jch,3)>0.0d0) numr=numr+1        !!>> HC 25-11-2020 Number of existing links in the r matrix (post-transcription reactions)
        enddo                                             !!>> HC 25-11-2020 
        do jch=2,nga                                      !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
           if(gen(ich)%e(jch)>0.0d0) nume=nume+1          !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
        enddo                                             !!>> HC 25-11-2020 
     enddo                                                !!>> HC 25-11-2020 
     
     if (ntipusadh/=0)then                                !!>> HC 25-11-2020 
        do ich=1, ntipusadh                               !!>> HC 25-11-2020 
           do jch=1, ntipusadh                            !!>> HC 25-11-2020 
             if(kadh(ich,jch)>0.0d0) numk=numk+1          !!>> HC 25-11-2020 Number of links in the kadh (B) matrix (specific adhesion)
           enddo                                          !!>> HC 25-11-2020 
        enddo                                             !!>> HC 25-11-2020 
     endif                                                !!>> HC 25-11-2020 
     
     probis=0.0d0;props=0                                 !!>> HC 25-11-2020 props is the total number of parameters that can change
     props=3*ng+numt+nume+numr+numk                       !!>> HC 25-11-2020 3*ng are diff, mu and mich of all genes
     probis(1)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mu
     probis(2)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability diff
     probis(3)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mich
     probis(4)=real(numt)/real(props)                     !!>> HC 25-11-2020 probability t matrix
     probis(5)=real(nume)/real(props)                     !!>> HC 25-11-2020 probability e matrix
     probis(6)=real(numr)/real(props)                     !!>> HC 25-11-2020 probability nww matrix
     probis(7)=real(numk)/real(props)                     !!>> HC 25-11-2020 probability kadh matrix
     cum=0.0d0                                            !!>> HC 25-11-2020
     call random_number(a)                                !!>> HC 25-11-2020
     do ich=1,size(probis)                                !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
        cum=cum+probis(ich)                               !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
        if (a<cum)then                                    !!>> HC 25-11-2020
           luckyprop=ich                                  !!>> HC 25-11-2020
           exit                                           !!>> HC 25-11-2020
        endif                                             !!>> HC 25-11-2020
     enddo                                                !!>> HC 25-11-2020
     mutacode=luckyprop                                   !!>> HC 16-9-2021  Store mutation code for output
     if(luckyprop==7) mutacode=mutacode-1                 !!>> HC 17-9-2021 !!WARNING!! We are not using reactions, so number 6 in mutacode if for kadh 
     
     if (luckyprop==1)then                                                !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DEGRADATION RATE (mu)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020
        if (i==0) i=1                                                     !!>> HC 25-11-2020 we randomly choose a gene i
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%mu                                                !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%mu==0.0d0) then                                        !!>> HC 25-11-2020
            gen(i)%mu=m_max*a                                             !!>> HC 25-11-2020 if mu=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
            inch=(1-2*a)*gen(i)%mu*mag_ismurate                           !!>> HC 27-11-2020 inch saves the magnitude of the change
            gen(i)%mu=gen(i)%mu+inch                                      !!>> HC 27-11-2020 more commonly, we increase or decrease mu randomly
            if (gen(i)%mu<0.0d0) gen(i)%mu=0.0d0                          !!>> HC 27-11-2020 negative values of mu are imposible
            if (gen(i)%mu>m_max)then                                      !!>> HC 27-11-2020 maximum degradation rate surpased, we apply the limit
               gen(i)%mu=gen(i)%mu-inch                                   !!>> HC 27-11-2020
               limit=0                                                    !!>> HC 27-11-2020
            endif                                                         !!>> HC 27-11-2020
        end if                                                            !!>> HC 27-11-2020
        newvalue=gen(i)%mu                                                !!>> HC 29-11-2021  we store the new value
     elseif (luckyprop==2)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
        if (i==0) i=1                                                     !!>> HC 25-11-2020
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%diffu                                             !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%diffu==0.0d0) then                                     !!>> HC 25-11-2020
           gen(i)%diffu=d_max*a                                           !!>> HC 25-11-2020 if diffu=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
           inch=(1-2*a)*gen(i)%diffu*mag_ismurate                         !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
           gen(i)%diffu=gen(i)%diffu+inch                                 !!>> HC 27-11-2020 more commonly, we increase or decrease diffu randomly
           if (gen(i)%diffu<0.0d0) gen(i)%diffu=0.0d0                     !!>> HC 27-11-2020 negative values of diffu are imposible
           if (gen(i)%diffu>d_max)then                                    !!>> HC 27-11-2020 maximum value of diffusion surpased, we apply the limit
              gen(i)%diffu=gen(i)%diffu-inch                              !!>> HC 27-11-2020
              limit=0                                                     !!>> HC 27-11-2020
           endif                                                          !!>> HC 27-11-2020
        end if                                                            !!>> HC 27-11-2020
        newvalue=gen(i)%diffu                                             !!>> HC 29-11-2021  we store the new value        
     elseif (luckyprop==3)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
        if (i==0) i=1                                                     !!>> HC 25-11-2020
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%mich                                              !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%mich==0.0d0) then                                      !!>> HC 25-11-2020
           gen(i)%mich=mich_max*a                                         !!>> HC 25-11-2020 if mich=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
           inch=(1-2*a)*gen(i)%mich*mag_ismurate                          !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
           gen(i)%mich=gen(i)%mich+inch                                   !!>> HC 27-11-2020 more commonly, we increase or decrease mich randomly  
           if (gen(i)%mich<0.0d0) gen(i)%mich=0.0d0                       !!>> HC 27-11-2020 negative values of mich are imposible
           if (gen(i)%mich>mich_max)then                                  !!>> HC 27-11-2020 maximum michaelis-menten constant surpased, we apply the limit
              gen(i)%mich=gen(i)%mich-inch                                !!>> HC 27-11-2020
              limit=0                                                     !!>> HC 27-11-2020
           endif                                                          !!>> HC 27-11-2020
        end if                                                            !!>> HC 25-11-2020 
        newvalue=gen(i)%mich                                              !!>> HC 29-11-2021  we store the new value  
     elseif (luckyprop==4)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
        if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020  
        allocate(plinki(1:numt))                                          !!>> HC 25-11-2020
        if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
        allocate(plinkj(1:numt))                                          !!>> HC 25-11-2020
        ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020
        do ich=1,ng                                                       !!>> HC 25-11-2020
           do jch=1, ng                                                   !!>> HC 25-11-2020
              if (gen(ich)%t(jch)>0.0d0)then                              !!>> HC 25-11-2020
                 ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (from 1 to the number of iterations=numt)
                 plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that are transcription factors
                 plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of genes j whose transcription is regulated
              endif                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
        enddo                                                             !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        luckylink=ceiling(a*numt)                                         !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
        if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        i=plinki(luckylink)                                               !!>> HC 25-11-2020 gene i is the randomly selected transcription factor
        j=plinkj(luckylink)                                               !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
        mutageni=i                                                        !!>> HC 29-11-2021    
        mutagenj=j                                                        !!>> HC 29-11-2021
        prevalue=gen(i)%t(j)                                              !!>> HC 29-11-2021 store the previous value
        inch=(1-2*a)*gen(i)%t(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
        gen(i)%t(j)=gen(i)%t(j)+inch                                      !!>> HC 27-11-2020 we apply a mutation of magnitude inch in the link between i and j
        if (abs(gen(i)%t(j))>w_max)then                                   !!>> HC 27-11-2020 if the limit for transcription interactions is surpased
           gen(i)%t(j)=gen(i)%t(j)-inch                                   !!>> HC 27-11-2020 we go back an activate the flag (limit) to make a new mutation
           limit=0                                                        !!>> HC 27-11-2020
        endif                                                             !!>> HC 27-11-2020
        newvalue=gen(i)%t(j)                                              !!>> HC 29-11-2021
     elseif (luckyprop==5)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
        if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020
        allocate(plinki(1:nume))                                          !!>> HC 25-11-2020
        if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020
        allocate(plinkj(1:nume))                                          !!>> HC 25-11-2020
        ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
        do ich=1,ng                                                       !!>> HC 25-11-2020
           do jch=2,nga                                                   !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
              if (gen(ich)%e(jch)>0.0d0)then                              !!>> HC 25-11-2020
                 ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (1 to nume)
                 plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that regulate cell behaviors
                 plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of regulated cell behaviors
              endif                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
        enddo                                                             !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        luckylink=ceiling(a*nume)                                         !!>> HC 25-11-2020 Randomly selecting 
        if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        i=plinki(luckylink)                                               !!>> HC 25-11-2020 Randomly selected gene i that regulates
        j=plinkj(luckylink)                                               !!>> HC 25-11-2020 Randomly selected behavior/property j
        mutageni=i                                                        !!>> HC 29-11-2021 store for output    
        mutagenj=j                                                        !!>> HC 29-11-2021 store for output
        prevalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output
        inch=(1-2*a)*gen(i)%e(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
        gen(i)%e(j)=gen(i)%e(j)+inch                                      !!>> HC 27-11-2020 we apply the mutation
        newval=gen(i)%e(j)                                                !!>> HC 27-11-2020 saving the new value for comparison
        if (newval>max_elim(j) .or. newval<min_elim(j))then               !!>> HC 27-11-2020 We compare with the limits of the cell behavior/property j
           gen(i)%e(j)=gen(i)%e(j)-inch                                   !!>> HC 27-11-2020 if the limit is surpased, we undo the change
           limit=0                                                        !!>> HC 27-11-2020 and activate the flag to make the IS mutation again
        endif                                                             !!>> HC 27-11-2020
        newvalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output  
     elseif (luckyprop==6)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
        if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020
        allocate(plinki(1:numr))                                          !!>> HC 25-11-2020                 (NOT USED RIGHT NOW)
        if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020
        allocate(plinkj(1:numr))                                          !!>> HC 25-11-2020
        ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020
        do ich=1,ng                                                       !!>> HC 25-11-2020
           do jch=1,ng                                                    !!>> HC 25-11-2020
              if (gen(ich)%e(jch)>0.0d0)then                              !!>> HC 25-11-2020
                 ord=ord+1                                                !!>> HC 25-11-2020
                 plinki(ord)=ich                                          !!>> HC 25-11-2020
                 plinkj(ord)=jch                                          !!>> HC 25-11-2020
              endif                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
        enddo                                                             !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        luckylink=ceiling(a*numr)                                         !!>> HC 25-11-2020
        if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        i=plinki(luckylink)                                               !!>> HC 25-11-2020
        j=plinkj(luckylink)                                               !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        gen(i)%r(j,3)=gen(i)%r(j,3)+(1-2*a)*gen(i)%r(j,3)*mag_ismurate    !!>> HC 25-11-2020
     elseif (luckyprop==7)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION SPECIFIC ADHESION (kadh matrix AKA B matrix)!!!!!!
            if(allocated(plinki)) deallocate(plinki)                      !!>> HC 25-11-2020
            allocate(plinki(1:numk))                                      !!>> HC 25-11-2020
            if(allocated(plinkj)) deallocate(plinkj)                      !!>> HC 25-11-2020
            allocate(plinkj(1:numk))                                      !!>> HC 25-11-2020 These vectors will store the existing links in kadh matrix
            ord=0;plinki=0;plinkj=0                                       !!>> HC 25-11-2020
            do ich=1, ntipusadh                                           !!>> HC 25-11-2020
               do jch=1, ntipusadh                                        !!>> HC 25-11-2020
                  if (kadh(ich,jch)>0.0d0)then                            !!>> HC 25-11-2020
                     ord=ord+1                                            !!>> HC 25-11-2020 order of the vectors (1 to numk)
                     plinki(ord)=ich                                      !!>> HC 25-11-2020 adhesion molecule i that specifically adheres to
                     plinkj(ord)=jch                                      !!>> HC 25-11-2020 adhesion molecule j (indexes in kadh no gene index!!)
                  endif                                                   !!>> HC 25-11-2020
               enddo                                                      !!>> HC 25-11-2020
            enddo                                                         !!>> HC 25-11-2020
            call random_number(a)                                         !!>> HC 25-11-2020
            luckylink=ceiling(a*numk)                                     !!>> HC 25-11-2020 Randomly picking
            if(luckylink==0)luckylink=1                                   !!>> HC 25-11-2020
            call random_number(a)                                         !!>> HC 25-11-2020
            i=plinki(luckylink)                                           !!>> HC 25-11-2020 This is the randomly picked adhesion molecule i 
            j=plinkj(luckylink)                                           !!>> HC 25-11-2020 that specifically adheres to adhesion molecule j 
            prevalue=kadh(i,j)                                            !!>> HC 29-11-2021 store for output
            inch=(1-2*a)*kadh(i,j)*mag_ismurate                           !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
            kadh(i,j)=kadh(i,j)+inch                                      !!>> HC 27-11-2020 and we mutate the interaction                
            if (kadh(i,j)>b_max .or. kadh(i,j)<b_min)then                 !!>> HC 27-11-2020 whe apply the limits b_max and b_min to the new value                
               kadh(i,j)=kadh(i,j)-inch                                   !!>> HC 27-11-2020 if the limits are surpased, we undo the mutation
               limit=0                                                    !!>> HC 27-11-2020 and activate the flag to mutate this individual again
            else                                                          !!>> HC 27-11-2020 if the mutation is accepted
               if (j .ne. i) kadh(j,i)=kadh(i,j)                          !!>> HC 27-11-2020 we have to be sure that is applied to both sides of the matrix
            endif                                                         !!>> HC 27-11-2020
            do ich=1,ng                                                   !!>> HC 16-9-2021
               if(gen(ich)%e(1)==i) mutageni=ich                          !!>> HC 16-9-2021 Store positions affected
               if(gen(ich)%e(1)==j) mutagenj=ich                          !!>> HC 16-9-2021
            enddo                                                         !!>> HC 16-9-2021
            newvalue=kadh(i,j)                                            !!>> HC 29-11-2021 store for output
     endif                                                                !!>> HC 25-11-2020
     
  else                      !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!!
     probs(1)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
     probs(2)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
     probs(3)=0d0           !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now)
     probs(4)=32d-2         !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
     probs(5)=2d-2          !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals
     call random_number(a)  !!>> HC 25-11-2020
     cum=0.0d0              !!>> HC 25-11-2020
     do ich=1,size(probs)   !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
        cum=cum+probs(ich)  !!>> HC 25-11-2020 taking into account probabilites in probs
        if (a<cum)then      !!>> HC 25-11-2020
           tind=ich         !!>> HC 25-11-2020
           exit             !!>> HC 25-11-2020
        endif               !!>> HC 25-11-2020
     enddo                  !!>> HC 25-11-2020
    if (tind .le. 2)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
       mutacode=4                                                         !!>> HC 29-11-2021 Store the matrix that will mutate
       do ich=1,ng                                                        !!>> HC 25-11-2020
         do jch=1,ng                                                      !!>> HC 25-11-2020
            if (gen(ich)%t(jch)==0) empty=empty+1                         !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
         enddo                                                            !!>> HC 25-11-2020 and the number of already existing links in t matrix
       enddo                                                              !!>> HC 25-11-2020
       full=ng*ng-empty                                                   !!>> HC 25-11-2020 
       if (tind==1 .and. full>0)then                                      !!>> HC 25-11-2020 !!!!!! DELECTION OF A LINK IN THE T MATRIX!!!!!!
          if (allocated(fullspoti)) deallocate(fullspoti)                 !!>> HC 25-11-2020
          allocate(fullspoti(1:full))                                     !!>> HC 25-11-2020
          if (allocated(fullspotj)) deallocate(fullspotj)                 !!>> HC 25-11-2020
          allocate(fullspotj(1:full))                                     !!>> HC 25-11-2020
          ord=0;fullspoti=0;fullspotj=0                                   !!>> HC 25-11-2020 These vectors store the already existing links
          do ich=1,ng                                                     !!>> HC 25-11-2020
            do jch=1,ng                                                   !!>> HC 25-11-2020
               if (gen(ich)%t(jch) .ne. 0.0d0)then                        !!>> HC 25-11-2020
                  ord=ord+1                                               !!>> HC 25-11-2020 order of vectors (1 to full)
                  fullspoti(ord)=ich                                      !!>> HC 25-11-2020 this gene i regulates the expression of
                  fullspotj(ord)=jch                                      !!>> HC 25-11-2020 gene j
               endif                                                      !!>> HC 25-11-2020
            enddo                                                         !!>> HC 25-11-2020
          enddo                                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          luckyg=ceiling(a*full)                                          !!>> HC 25-11-2020 Randomly pickin an interaction to delete
          if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          i=fullspoti(luckyg)                                             !!>> HC 25-11-2020 Removing interaction between gene i and j
          j=fullspotj(luckyg)                                             !!>> HC 25-11-2020
          mutageni=i                                                      !!>> HC 29-11-2021 Store for output
          mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
          prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
          gen(i)%t(j)=0.0d0                                               !!>> HC 25-11-2020
          newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 
       else                                                               !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
          if (allocated(freespoti)) deallocate(freespoti)                 !!>> HC 25-11-2020
          allocate(freespoti(1:empty))                                    !!>> HC 25-11-2020
          if (allocated(freespotj)) deallocate(freespotj)                 !!>> HC 25-11-2020
          allocate(freespotj(1:empty))                                    !!>> HC 25-11-2020
          ord=0;freespoti=0;freespotj=0                                   !!>> HC 25-11-2020 These vectors store the free spots
          do ich=1,ng                                                     !!>> HC 25-11-2020
            do jch=1,ng                                                   !!>> HC 25-11-2020
               if (gen(ich)%t(jch)==0)then                                !!>> HC 25-11-2020
                  ord=ord+1                                               !!>> HC 25-11-2020 order of the vectors (1 to empty)
                  freespoti(ord)=ich                                      !!>> HC 25-11-2020 gene i has NOT and interaction
                  freespotj(ord)=jch                                      !!>> HC 25-11-2020 with gene j
               endif                                                      !!>> HC 25-11-2020
            enddo                                                         !!>> HC 25-11-2020
          enddo                                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          luckyg=ceiling(a*empty)                                         !!>> HC 25-11-2020 Randomly picking a spot
          if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          call random_number(b)                                           !!>> HC 25-11-2020
          if(b<0.50d0)then;signo=1;else;signo=-1;endif                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
          i=freespoti(luckyg)                                             !!>> HC 25-11-2020
          j=freespotj(luckyg)                                             !!>> HC 25-11-2020 Making new interaction
          mutageni=i                                                      !!>> HC 29-11-2021 Store for output
          mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
          prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
          gen(i)%t(j)=a*signo*w_max                                       !!>> HC 25-11-2020
          newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 
       endif                                                              !!>> HC 25-11-2020
    elseif (tind==3)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
    !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW        !!>> HC 25-11-2020
    elseif (tind==4)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e and kadh matrixes)!!!!!!
           empty=0                                                        !!>> HC 25-11-2020
           do ich=1,ng                                                    !!>> HC 25-11-2020
              do jch=1,nga                                                !!>> HC 25-11-2020
                 if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                 if (gen(ich)%e(jch)==0.0d0) then                         !!>> HC 27-11-2020
                     empty=empty+1                                        !!>> HC 27-11-2020 Number of empty spots in e Matrix
                 else                                                     !!>> HC 27-11-2020
                    full=full+1                                           !!>> HC 27-11-2020 Number of existing interactions in e matrix
                 endif                                                    !!>> HC 27-11-2020
              enddo                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
           if (ntipusadh/=0)then                                          !!>> HC 25-11-2020
               do ich=1, ntipusadh                                        !!>> HC 25-11-2020
                  do jch=1, ntipusadh                                     !!>> HC 25-11-2020
                     if(kadh(ich,jch)==0.0d0)then                         !!>> HC 27-11-2020
                        empty=empty+1                                     !!>> HC 27-11-2020 Number of empty spots in kadh Matrix (AKA B Matrix)
                     else                                                 !!>> HC 27-11-2020
                        full=full+1                                       !!>> HC 27-11-2020 Number of existing interaction in kadh Matrix
                     endif                                                !!>> HC 27-11-2020
                  enddo                                                   !!>> HC 25-11-2020
               enddo                                                      !!>> HC 25-11-2020
           endif                                                          !!>> HC 25-11-2020
           call random_number(a)                                          !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
           if (a<0.50d0 .and. full>0.or.empty==0)then                     !!>> HC 10-5-2022 !!!REMOVE INTERACTION           
              if (allocated(fullspoti)) deallocate(fullspoti)             !!>> HC 25-11-2020
              allocate(fullspoti(1:full))                                 !!>> HC 25-11-2020
              if (allocated(fullspotj)) deallocate(fullspotj)             !!>> HC 25-11-2020
              allocate(fullspotj(1:full))                                 !!>> HC 25-11-2020
              if (allocated(isadh)) deallocate(isadh)                     !!>> HC 25-11-2020
              allocate(isadh(1:full))                                     !!>> HC 25-11-2020
              ord=0;fullspoti=0;fullspotj=0;isadh=0                       !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
              do ich=1,ng                                                 !!>> HC 25-11-2020
                 do jch=1,nga                                             !!>> HC 25-11-2020 we start in 2 because to stop being an adhesion mol you need to will every adh interaction
                    if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
                    if (gen(ich)%e(jch).ne. 0.0d0)then                    !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                       ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                       fullspoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                       fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
                       isadh(ord)=0                                       !!>> HC 25-11-2020 this is not in kadh
                    endif                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              enddo                                                       !!>> HC 25-11-2020
              if (ntipusadh/=0)then                                       !!>> HC 25-11-2020
                 do ich=1, ntipusadh                                      !!>> HC 25-11-2020
                    do jch=1, ntipusadh                                   !!>> HC 25-11-2020
                       if (kadh(ich,jch).ne. 0.0d0)then                   !!>> HC 25-11-2020
                          ord=ord+1                                       !!>> HC 25-11-2020
                          fullspoti(ord)=ich                              !!>> HC 25-11-2020 adhesion molecule i 
                          fullspotj(ord)=jch                              !!>> HC 25-11-2020 adheres to adhesion molecule j
                          isadh(ord)=1                                    !!>> HC 25-11-2020 this interaction is in kadh (adhesion mols no genes)
                       endif                                              !!>> HC 25-11-2020
                    enddo                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              endif                                                       !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              luckyg=ceiling(a*real(full))                                !!>> HC 25-11-2020 Randomly picking an interaction to kill
              i=fullspoti(luckyg)                                         !!>> HC 25-11-2020
              j=fullspotj(luckyg)                                         !!>> HC 25-11-2020
              isad=isadh(luckyg)                                          !!>> HC 25-11-2020     
              if (isad==0)then                                            !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
                 mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
                 mutageni=i                                               !!>> HC 29-11-2021 Store positions affected
                 mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
                 prevalue=gen(i)%e(j)                                     !!>> HC 29-11-2021 save for output the prev value
              else                                                        !!>> HC 29-11-2021 
                 mutacode=6                                               !!>> HC 29-11-2021  The mutation is in the B matrix (kadh)
                 do ich=1,ng                                              !!>> HC 29-11-2021
                    if(gen(ich)%e(1)==i) mutageni=ich                     !!>> HC 29-11-2021 Store positions affected
                    if(gen(ich)%e(1)==j) mutagenj=ich                     !!>> HC 29-11-2021
                 enddo                                                    !!>> HC 29-11-2021 
                 prevalue=kadh(i,j)                                       !!>> HC 29-11-2021  save for output the prev value
              endif                                                       !!>> HC 29-11-2021  --END-- 
              newvalue=0.0d0                                              !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
              if (isad==0)then                                            !!>> HC 25-11-2020 if the interaction is not in kadh
                 if (j==1)then                                            !!>> HC 22-11-2021 DELECTION OF AN ADHESION MOLECULE
                    if (ntipusadh.le.1)then                               !!>> HC 22-11-2021 This is the only adhesion molecule left
                       if (allocated(kadh)) deallocate(kadh)              !!>> HC 22-11-2021 we remove all the kadh matrix
                       ntipusadh=0                                        !!>> HC 22-11-2021 and there are no more adhesion molecules
                    else                                                  !!>> HC 22-11-2021 There are some other molecules left
                       ntipusadh=ntipusadh-1                              !!>> HC 22-11-2021
                       if(allocated(copkadh))deallocate(copkadh)          !!>> HC 22-11-2021
                       allocate(copkadh(ntipusadh,ntipusadh))             !!>> HC 22-11-2021 Copy the data in a smaller matrix (one col and row less)
                       copkadh=0.0d0;ord=0;ord2=0                         !!>> HC 22-11-2021
                       do ich=1,ntipusadh+1                               !!>> HC 22-11-2021 Copy data
                          if (ich==int(gen(i)%e(j)))cycle                 !!>> HC 22-11-2021
                          ord=ord+1                                       !!>> HC 22-11-2021
                          ord2=0                                          !!>> HC 20-12-2021 
                          do jch=1,ntipusadh+1                            !!>> HC 22-11-2021
                             if (jch==int(gen(i)%e(j)))cycle              !!>> HC 22-11-2021
                             ord2=ord2+1                                  !!>> HC 22-11-2021
                             copkadh(ord,ord2)=kadh(ich,jch)              !!>> HC 22-11-2021
                           enddo                                          !!>> HC 22-11-2021
                       enddo                                              !!>> HC 22-11-2021
                       if(allocated(kadh))deallocate(kadh)                !!>> HC 25-11-2020 reallocate kadh with one molecule less
                       allocate(kadh(ntipusadh,ntipusadh))                !!>> HC 22-11-2021
                       kadh=0.0d0                                         !!>> HC 22-11-2021
                       kadh(1:ntipusadh,1:ntipusadh)=copkadh(1:ntipusadh,1:ntipusadh)  !!>> HC 22-11-2021 Save previous data
                       do ich=1,ng                                        !!>> HC 22-11-2021
                          if (gen(ich)%e(1).le.gen(i)%e(j))cycle          !!>> HC 22-11-2021 Update the indexes of the resting adhesion molecules
                          gen(ich)%e(1)=gen(ich)%e(1)-1                   !!>> HC 22-11-2021
                       enddo                                              !!>> HC 22-11-2021
                    endif                                                 !!>> HC 22-11-2021
                 endif                                                    !!>> HC 22-11-2021
                 gen(i)%e(j)=0.0d0                                        !!>> HC 25-11-2020 we kill it easily
              else                                                        !!>> HC 25-11-2020 But if it is in kadh, shit comes...
                 kadh(i,j)=0.0d0                                          !!>> HC 25-11-2020
                 kadh(j,i)=0.0d0                                          !!>> HC 25-11-2020 we kill the interaction
                 stilladi=0; stilladj=0                                   !!>> HC 25-11-2020
                 do ich=1,ntipusadh                                       !!>> HC 25-11-2020 we check if the adhesion molecules i and j have some 
                    if (kadh(i,ich) .ne. 0.0d0) stilladi=1                !!>> HC 25-11-2020 interactions left in kadh
                    if (kadh(j,ich) .ne. 0.0d0) stilladj=1                !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
                 if (stilladi==0 .and. stilladj==1)then                   !!>> HC 25-11-2020 adhesion molecule i has no interactions left but j has some
                    ntipusadh=ntipusadh-1                                 !!>> HC 25-11-2020
                    if(allocated(copkadh))deallocate(copkadh)             !!>> HC 25-11-2020 copy kadh in a smaller matrix
                    allocate(copkadh(ntipusadh,ntipusadh))                !!>> HC 25-11-2020
                    copkadh=0.0d0;ord=0                                   !!>> HC 25-11-2020
                    do ich=1,ntipusadh+1                                  !!>> HC 25-11-2020 We delete the row i and column i of kadh matrix 
                       if (ich==i)cycle                                   !!>> HC 25-11-2020
                       ord=ord+1                                          !!>> HC 25-11-2020
                       ord2=0                                             !!>> HC 25-11-2020
                       do jch=1,ntipusadh+1                               !!>> HC 25-11-2020
                          if(jch==i) cycle                                !!>> HC 25-11-2020
                          ord2=ord2+1                                     !!>> HC 25-11-2020
                          copkadh(ord,ord2)=kadh(ich,jch)                 !!>> HC 25-11-2020
                       enddo                                              !!>> HC 25-11-2020
                    enddo                                                 !!>> HC 25-11-2020
                    if(allocated(kadh))deallocate(kadh)                   !!>> HC 25-11-2020 realocate kadh and recovering saved data
                    allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 25-11-2020
                    kadh=0.0d0; kadh=copkadh                              !!>> HC 25-11-2020
                    do ich=1,ng                                           !!>> HC 25-11-2020 Remove e(1) in the e matrix of the deleted adhesion molecule i
                       if (gen(ich)%e(1)==i) gen(ich)%e(1)=0              !!>> HC 25-11-2020 and adjust the indexes of all the other adhesion molecules
                       if (gen(ich)%e(1)>i) gen(ich)%e(1)=gen(ich)%e(1)-1 !!>> HC 25-11-2020
                    enddo                                                 !!>> HC 25-11-2020
                 endif                                                    !!>> HC 25-11-2020                          
                 if (stilladi==1 .and. stilladj==0)then                   !!>> HC 25-11-2020  adhesion molecule j has no interactions left but i has some 
                    ntipusadh=ntipusadh-1                                 !!>> HC 25-11-2020
                    if(allocated(copkadh))deallocate(copkadh)             !!>> HC 25-11-2020
                    allocate(copkadh(ntipusadh,ntipusadh))                !!>> HC 25-11-2020 copying kadh in a smaller matrix
                    copkadh=0.0d0;ord=0                                   !!>> HC 25-11-2020
                    do ich=1,ntipusadh+1                                  !!>> HC 25-11-2020  We delete the column j and row j of kadh matrix 
                       if (ich==j)cycle                                   !!>> HC 25-11-2020
                       ord=ord+1                                          !!>> HC 25-11-2020
                       ord2=0                                             !!>> HC 25-11-2020
                       do jch=1,ntipusadh+1                               !!>> HC 25-11-2020
                          if(jch==j) cycle                                !!>> HC 25-11-2020
                          ord2=ord2+1                                     !!>> HC 25-11-2020
                          copkadh(ord2,ord)=kadh(ich,jch)                 !!>> HC 25-11-2020
                       enddo                                              !!>> HC 25-11-2020
                    enddo                                                 !!>> HC 25-11-2020
                    if(allocated(kadh))deallocate(kadh)                   !!>> HC 25-11-2020 reallocating kadh and recovering saved data
                    allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 25-11-2020
                    kadh=0.0d0; kadh=copkadh                              !!>> HC 25-11-2020
                    do ich=1,ng                                           !!>> HC 25-11-2020 Remove e(1) in the e matrix of the deleted adhesion molecule j
                       if (gen(ich)%e(1)==j) gen(ich)%e(1)=0              !!>> HC 25-11-2020 and adjust the indexes of all the other adhesion molecules
                       if (gen(ich)%e(1)>j) gen(ich)%e(1)=gen(ich)%e(1)-1 !!>> HC 25-11-2020 
                    enddo                                                 !!>> HC 25-11-2020
                 endif                                                    !!>> HC 25-11-2020
                 if (stilladi==0 .and. stilladj==0)then                   !!>> HC 25-11-2020 it was the last adhesion interaction for both i and j molecules
                    if (i==j)then                                            !!>> HC 17-9-2021  i and j molecules are actually the SAME
                       ntipusadh=ntipusadh-1                                 !!>> HC 17-9-2021  then there is only one molecule to remove
                       if(allocated(copkadh))deallocate(copkadh)             !!>> HC 17-9-2021
                       allocate(copkadh(ntipusadh,ntipusadh))                !!>> HC 17-9-2021
                       copkadh=0.0d0;ord=0                                   !!>> HC 17-9-2021 Removing  column i==j and row i==j of kadh matrix
                       do ich=1,ntipusadh+1                                  !!>> HC 17-9-2021
                          if (ich==i)cycle                                   !!>> HC 17-9-2021
                          ord=ord+1                                          !!>> HC 17-9-2021
                          ord2=0                                             !!>> HC 17-9-2021
                          do jch=1,ntipusadh+1                               !!>> HC 17-9-2021
                             if(jch==j) cycle                                !!>> HC 17-9-2021
                             ord2=ord2+1                                     !!>> HC 17-9-2021
                             copkadh(ord,ord2)=kadh(ich,jch)                 !!>> HC 17-9-2021
                          enddo                                              !!>> HC 17-9-2021
                       enddo                                                 !!>> HC 17-9-2021
                       if(allocated(kadh))deallocate(kadh)                   !!>> HC 17-9-2021 Reallocate kadh and recover saved data
                       allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 17-9-2021
                       kadh=0.0d0; kadh=copkadh                              !!>> HC 17-9-2021
                       do ich=1,ng                                           !!>> HC 17-9-2021 Remove the e(1) of adhesion molecule i
                          if (gen(ich)%e(1)==i) gen(ich)%e(1)=0.0d0          !!>> HC 17-9-2021 and adjust the indexes of all the other adhesion molecules
                          if (gen(ich)%e(1)>i) gen(ich)%e(1)=gen(ich)%e(1)-1 !!>> HC 17-9-2021
                       enddo                                                 !!>> HC 17-9-2021
                    else
                       ntipusadh=ntipusadh-2                                 !!>> HC 25-11-2020
                       if(allocated(copkadh))deallocate(copkadh)             !!>> HC 25-11-2020
                       allocate(copkadh(ntipusadh,ntipusadh))                !!>> HC 25-11-2020 Copy the data in a smaller matrix (two cols and rows less)
                       copkadh=0.0d0;ord=0                                   !!>> HC 25-11-2020 Removing  column i, row i, column j and row j of kadh matrix
                       do ich=1,ntipusadh+2                                  !!>> HC 25-11-2020
                          if (ich==i .or. ich==j)cycle                       !!>> HC 25-11-2020
                          ord=ord+1                                          !!>> HC 25-11-2020
                          ord2=0                                             !!>> HC 25-11-2020
                          do jch=1,ntipusadh+2                               !!>> HC 25-11-2020
                             if(jch==j .or. jch==i) cycle                    !!>> HC 25-11-2020
                             ord2=ord2+1                                     !!>> HC 25-11-2020
                             copkadh(ord,ord2)=kadh(ich,jch)                 !!>> HC 25-11-2020
                          enddo                                              !!>> HC 25-11-2020
                       enddo                                                 !!>> HC 25-11-2020
                       if(allocated(kadh))deallocate(kadh)                   !!>> HC 25-11-2020 Reallocate kadh and recover saved data
                       allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 25-11-2020
                       kadh=0.0d0; kadh=copkadh                              !!>> HC 25-11-2020
                       if(i>j)then;mnadh=j;mxadh=i;else;mnadh=i;mxadh=j;endif!!>> HC 25-11-2020
                       do ich=1,ng                                           !!>> HC 25-11-2020 Remove the e(1) of adhesion molecules i and j and
                          if (gen(ich)%e(1)==i) gen(ich)%e(1)=0.0d0          !!>> HC 25-11-2020 and adjust the indexes of all the other adhesion molecules
                          if (gen(ich)%e(1)==j) gen(ich)%e(1)=0.0d0          !!>> HC 25-11-2020
                          if (gen(ich)%e(1)>mnadh .and. gen(ich)%e(1)<mxadh) gen(ich)%e(1)=gen(ich)%e(1)-1  !!>> HC 25-11-2020
                          if (gen(ich)%e(1)>mxadh) gen(ich)%e(1)=gen(ich)%e(1)-2                            !!>> HC 25-11-2020
                       enddo                                                                                !!>> HC 25-11-2020
                    endif
                 endif                                                                                   !!>> HC 25-11-2020
              endif              
           else                                                           !!>> HC 25-11-2020 !ADD AN INTERACTION   
              if (allocated(freespoti)) deallocate(freespoti)             !!>> HC 25-11-2020
              allocate(freespoti(1:empty))                                !!>> HC 25-11-2020
              if (allocated(freespotj)) deallocate(freespotj)             !!>> HC 25-11-2020
              allocate(freespotj(1:empty))                                !!>> HC 25-11-2020
              if (allocated(isadh)) deallocate(isadh)                     !!>> HC 25-11-2020 These vectors will store the free spots in e and kadh matrixes
              allocate(isadh(1:empty))                                    !!>> HC 25-11-2020 and whether the spot is in kadh matrix or not
              ord=0;freespoti=0;freespotj=0; isadh=0                      !!>> HC 25-11-2020 
              do ich=1,ng                                                 !!>> HC 25-11-2020
                 do jch=1,nga                                             !!>> HC 25-11-2020
                    if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                    if (gen(ich)%e(jch)==0.0d0)then                       !!>> HC 25-11-2020
                       ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                       freespoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                       freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
                       isadh(ord)=0                                       !!>> HC 25-11-2020 and this is not in kadh matrix
                    endif                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              enddo                                                       !!>> HC 25-11-2020
              if (ntipusadh/=0)then                                       !!>> HC 25-11-2020
                 do ich=1, ntipusadh                                      !!>> HC 25-11-2020
                    do jch=1, ntipusadh                                   !!>> HC 25-11-2020
                       if (kadh(ich,jch)==0.0d0)then                      !!>> HC 25-11-2020
                          ord=ord+1                                       !!>> HC 25-11-2020 order of the vector (1 to nume)
                          freespoti(ord)=ich                              !!>> HC 25-11-2020 adhesion molecule i adheres t
                          freespotj(ord)=jch                              !!>> HC 25-11-2020 adhesion molecule j
                          isadh(ord)=1                                    !!>> HC 25-11-2020 and this is in kadh matrix
                       endif                                              !!>> HC 25-11-2020
                    enddo                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              endif                                                       !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              luckyg=ceiling(a*empty)                                     !!>> HC 25-11-2020 Randomly picking an interaction
              if(luckyg==0)luckyg=1                                       !!>> HC 25-11-2020
              i=freespoti(luckyg)                                         !!>> HC 25-11-2020
              j=freespotj(luckyg)                                         !!>> HC 25-11-2020
              isad=isadh(luckyg)                                          !!>> HC 25-11-2020
              if (isad==0.and.j>1)then                                    !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
                 mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
                 mutageni=i                                               !!>> HC 29-11-2021 Store positions affected
                 mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
                 prevalue=gen(i)%e(j)                                     !!>> HC 29-11-2021 save for output the prev value
              else                                                        !!>> HC 29-11-2021 
                 mutacode=6                                               !!>> HC 29-11-2021  The mutation is in the B matrix (kadh)
                 do ich=1,ng                                              !!>> HC 29-11-2021
                    if(gen(ich)%e(1)==i) mutageni=ich                     !!>> HC 29-11-2021 Store positions affected
                    if(gen(ich)%e(1)==j) mutagenj=ich                     !!>> HC 29-11-2021
                 enddo                                                    !!>> HC 29-11-2021 
                 prevalue=0.0d0                                           !!>> HC 29-11-2021  save for output the prev value
              endif                                                       !!>> HC 29-11-2021  --END--
              if (j==1 .and. isad .ne. 1 )then                            !!>> HC 25-11-2020 !!ADD AN ADHESION MOLECULE
                 ntipusadh=ntipusadh+1                                    !!>> HC 25-11-2020
                 gen(i)%e(1)=ntipusadh                                    !!>> HC 25-11-2020 
                 if (ntipusadh>1)then                                     !!>> HC 25-11-2020 there are already other adhesion molecules
                    if(allocated(copkadh))deallocate(copkadh)             !!>> HC 25-11-2020
                    allocate(copkadh(ntipusadh,ntipusadh))                !!>> HC 25-11-2020
                    copkadh(:ntipusadh-1,:ntipusadh-1)=kadh(:ntipusadh-1,:ntipusadh-1)      !!>> HC 25-11-2020
                    deallocate(kadh)                                      !!>> HC 25-11-2020
                    allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 25-11-2020
                    kadh=0.0d0                                            !!>> HC 25-11-2020
                    kadh(:ntipusadh-1,:ntipusadh-1)=copkadh(:ntipusadh-1,:ntipusadh-1)      !!>> HC 25-11-2020  We copy the data, realocate kadh and recover the data               
                    call random_number(a)                                 !!>> HC 25-11-2020 adhesion molecule 
                    call random_number(c)                                 !!>> HC 25-11-2020 magnitude of interaction
                        !!>> HC 25-11-2020
                    jj=int(ntipusadh*a)+1                                 !!>> HC 25-11-2020 pick random gene to adhere to
                    if (jj>ntipusadh) jj=ntipusadh                        !!>> HC 25-11-2020
                    anch=b_max-b_min                                      !!>> HC 28-11-2020 anch saves the interval
                    kadh(ntipusadh,jj)=anch*c+b_min                       !!>> HC 28-11-2020 This introduces the new interaction with a value in the range of the limits
                    kadh(jj,ntipusadh)=anch*c+b_min                       !!>> HC 28-11-2020
                    do ich=1,ng                                           !!>> HC 16-9-2021
                       if(gen(ich)%e(1)==i)  mutageni=ich                 !!>> HC 29-11-2021 Store positions affected
                       if(gen(ich)%e(1)==jj) mutagenj=ich                 !!>> HC 29-11-2021 Store positions affected positions
                    enddo                                                 !!>> HC 29-11-2021 
                    newvalue=kadh(ntipusadh,jj)                           !!>> HC 29-11-2021 store the new value
                 else                                                     !!>> HC 25-11-2020 This is the first adhesion molecule
                    if(allocated(kadh)) deallocate(kadh)                  !!>> HC 25-11-2020
                    allocate(kadh(ntipusadh,ntipusadh))                   !!>> HC 25-11-2020 no need of saving and recovering data but we have to allocate kadh
                    call random_number(c)                                 !!>> HC 25-11-2020 Magnitude of the interaction
                    anch=b_max-b_min                                      !!>> HC 28-11-2020 anch saves the interval
                    kadh(ntipusadh,ntipusadh)=anch*c+b_min                !!>> HC 28-11-2020 This introduces the new interaction with a value in the range of the limits
                    newvalue=kadh(ntipusadh,ntipusadh)                    !!>> HC 29-11-2021 store the new value
                    do ich=1,ng                                           !!>> HC 23-3-2022 
                       if (gen(ich)%e(1)>0.0d0)then                       !!>> HC 23-3-2022 Store positions affected
                          mutageni=ich                                    !!>> HC 23-3-2022 
                          mutagenj=ich                                    !!>> HC 23-3-2022
                       endif                                              !!>> HC 23-3-2022
                    enddo                                                 !!>> HC 23-3-2022
                 endif                                                    !!>> HC 25-11-2020
              elseif (j .ge. 2 .and. isad .ne. 1 )then                    !!>> HC 25-11-2020 !!ADD A NEW REGULATOR OF CELL BEHAVIORS/PROPERTIES (e matrix)
                     call random_number(c)                                !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
                     anch=max_elim(j)-min_elim(j)                         !!>> HC 28-11-2020 we introduce an interaction in the range of the limits
                     gen(i)%e(j)=anch*c+min_elim(j)                       !!>> HC 28-11-2020 
                                                                          !!>> HC 12-11-2020 !!!ATENTION!! this is a trick to make PCP more likely to appear by mutation
                     if (j==nparam_per_node+8)then                        !!>> HC 7-12-2020
                        gen(3)%e(nparam_per_node+17)=0.20d0               !!>> HC 12-11-2020 !!!ATENTION!! gene 3 is homogeneously expressed in the blastula ics 
                        gen(6)%e(nparam_per_node+17)=0.20d0               !!>> HC 7-12-2020
                     endif                                                !!>> HC 7-12-2020
                                                                          !!>> HC 12-11-2020 !!!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
                     newvalue=gen(i)%e(j)                                 !!>> HC 29-11-2020  Store the new value                 
              elseif (isad==1)then                                        !!>> HC 25-11-2020 !!ADD A NEW INTERACTION OF SPECIFIC ADHESION (kadh matrix)
                     call random_number(c)                                !!>> HC 25-11-2020 Magnitude
                     anch=b_max-b_min                                     !!>> HC 28-11-2020 anch saves the lenght of the valid interval
                     kadh(i,j)=anch*c+b_min                               !!>> HC 28-11-2020 we introduce an interaction in the range of the limits
                     kadh(j,i)=anch*c+b_min                               !!>> HC 28-11-2020
                     newvalue=kadh(i,j)                                   !!>> HC 29-11-2020  Store the new value  
              endif                                                       !!>> HC 25-11-2020
           endif                                                          !!>> HC 25-11-2020
    elseif (tind==5)then                                                  !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
           mutacode=7                                                     !!>> HC 29-11-2021 Store that this is a gene duplication/delection
           call random_number(a)                                          !!>> HC 10-9-2021
           i=ceiling(a*(ng))                                              !!>> HC 10-9-2021
           if (i==0) i=1                                                  !!>> HC 25-11-2020
           mutageni=i                                                     !!>> HC 16-9-2021 save the chosen gene for the output
           if (ng>2)then                                                  !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              if (a<0.50d0)then                                           !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
                 mutagenj=-1                                              !!>> HC 29-11-2021 store that this is a delection
                 call deletion(i)                                         !!>> HC 25-11-2020
              else                                                        !!>> HC 25-11-2020
                 mutagenj=1                                               !!>> HC 29-11-2021 store that this is a duplication
                 call duplication(i)                                      !!>> HC 25-11-2020
              endif                                                       !!>> HC 25-11-2020
           else                                                           !!>> HC 25-11-2020
              mutagenj=1                                                  !!>> HC 29-11-2021 store that this is a duplication
              call duplication(i)                                         !!>> HC 25-11-2020
           endif                                                          !!>> HC 25-11-2020
    endif                                                                 !!>> HC 25-11-2020
  endif                                                                   !!>> HC 25-11-2020
     
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 
!!>> HC 16-9-2021 !! ATENTION !! We are not using pre/post forms right now, this should be changed if you want them!!      !!>> HC 16-9-2021
!  if (allocated(npag)) deallocate(npag)                             !!>> HC 25-11-2020  This last chunk of code was inherited from muta
!  allocate(npag(nga))                                               !!>> HC 25-11-2020  it realocates matrixes 
!  npag=0                                                            !!>> HC 25-11-2020

!  if (allocated(whonpag)) deallocate(whonpag)                       !!>> HC 25-11-2020
!  allocate(whonpag(nga,ng))                                         !!>> HC 25-11-2020
!  whonpag=0                                                         !!>> HC 25-11-2020

!  call update_npag                                                  !!>> HC 25-11-2020

!  do i=1,ng                                                         !!>> HC 25-11-2020
!     do j=1,gen(i)%npre                                             !!>> HC 25-11-2020
!        if (gen(i)%pre(j)==0) then                                  !!>> HC 25-11-2020
!           print *,i,j,"aqui hi ha un zero a pre"                   !!>> HC 25-11-2020
!           stop                                                     !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020

!     do j=1,gen(i)%npost                                            !!>> HC 25-11-2020
!        if (gen(i)%post(j)==0) then                                 !!>> HC 25-11-2020
!        print *,i,j,"aqui hi ha un zero a post"                     !!>> HC 25-11-2020
!        stop                                                        !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020
!  end do                                                            !!>> HC 25-11-2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 

  inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
  if (ng<1) then                                                    !!>> HC 25-11-2020
     print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
     inviable=1                                                     !!>> HC 25-11-2020
  else                                                              !!>> HC 25-11-2020
     do j=1,ng                                                      !!>> HC 25-11-2020
        if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
           inviable=0                                               !!>> HC 25-11-2020
           goto 37                                                  !!>> HC 25-11-2020
        end if                                                      !!>> HC 25-11-2020
     end do                                                         !!>> HC 25-11-2020
     print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
     inviable=1                                                     !!>> HC 25-11-2020
  end if                                                            !!>> HC 25-11-2020
37 continue                                                         !!>> HC 25-11-2020
       
end subroutine suremuta

!***********************************************************************************************************'

!***********************************************************************************************************'

subroutine suremuta_no_kadh(limit, inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue,rangfile)
!!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file 
!!>> HC 14-4-2023 This version is the same but ignoring specific adhesion molecules
!!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
!!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
!!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
!!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
!!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
  implicit none
  integer :: limit,inviable,mind,tind
  real*8, allocatable :: copkadh(:,:)
  real*8 :: probmod,xx !the modifier used to bias T mutations of W towards interaction loss, and the temp var to store the ratio between empty and filled w positions (for readability only)
  integer genadh 
  integer :: props,nume,numt,numr,numk !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
  integer :: luckyprop, luckylink, luckyg  !!>> HC 24-11-2020 selected gene, interaction or property
  integer :: ich,jch,lch,ord,ord2  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
  integer,allocatable,dimension(:) :: plinki,plinkj !!>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
  integer,allocatable,dimension(:) :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM
  integer,allocatable,dimension(:) :: isadh !!>> HC 24-11-2020 if 1 the interaction is in the kadh matrix (TM)
  integer :: full,empty !!>> HC 24-11-2020 full an empty spaces for T mutations
  real*8, dimension(1:5) :: probs !!>> HC 24-11-2020 relative probabilities of different T mutations
  real*8, dimension(1:7) :: probis !!>> HC 24-11-2020 relative probabilities of different IS mutations
  real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 

  integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  real*8 :: cum !!>> HC 24-11-2020 funny position to store cumulative Fi
  integer :: stilladi, stilladj,isad, mnadh, mxadh !!>> HC 24-11-2020 pointers related to T mutations in the Kadh matrix
  integer :: signo  !!>> HC 24-11-2020 sign of the mutation  
  real*8 :: inch, newval, anch !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
  real*8, dimension(1:5) :: bhc !!>>HC 10-9-2021 Patch to a problem with biased random numbers
  integer :: mutacode,mutageni,mutagenj  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
  real*8 :: prevalue, newvalue
  character*400 :: rangfile !!>> 6-10-2021 file where the ranges are stored
  real*8, dimension (1:5) :: min_glim, max_glim                                  !!>>HC 6-10-2021
  real*8 :: t_mag !!>> HC 31-5-2023 Magnitude of the initial value in T mutations (should be a free parameter)
  
  !************************** LIMITS AND RANGES **************************!        
  limit=1                                               !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
  rembeh=0; max_elim=0.0d0; min_elim=0.0d0              !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
  max_glim=0.0d0; min_glim=0.0d0                        !!>> HC 6-10-2020 Same for other gene properties 
  call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
  
  d_max=max_glim(1)    ; d_min=min_glim(1)                       !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020
  m_max=max_glim(2)    ; m_min=min_glim(2)                       !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
  mich_max=max_glim(3) ; mich_min=min_glim(3)                    !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
  w_max=max_glim(4)    ; w_min=min_glim(4)                       !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
  b_max=max_glim(5)    ; b_min=min_glim(5)                       !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

  !************************************************************************!
  if (allocated(dupi)) deallocate(dupi)
  if (allocated(father_son)) deallocate(father_son)
  if (allocated(todel)) deallocate(todel)
  allocate(dupi(ng*ng))       !since each gene can no duplicate more than ng times per mutation round
  allocate(father_son(ng*ng)) !since each gene can no duplicate more than ng times per mutation round
  allocate(todel(ng*ng*10))   !!>>AL 9-4-24: why not only ng*ng..?
  dupi=0
  father_son=0
  todel=0
  empty=0; full=0             !!>> HC 24-11-2020
  inch=0.0d0; anch=0.0d0      !!>> HC 28-11-2020
  mutacode=0                  !!>> HC 16-9-2021  
  mutageni=0                  !!>> HC 29-11-2021
  mutagenj=0                  !!>> HC 29-11-2021
  prevalue=0.0d0              !!>> HC 29-11-2021
  newvalue=0.0d0              !!>> HC 29-11-2021
  t_mag=0.10d0                !!>> HC 31-5-2023
  call random_number(bhc)     !!>> HC 25-11-2020 THESE ARE BIASED I DO NOT KNOW WHY !!>>AL 5-4-24 chechar si estan sesgados realmente o no
  if (mind==1)then                                        !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!! IS MUTATIONS!!!!!!!!!!!!!!!!!!!!!
     nume=0; numt=0; numk=0; numr=0                       !!>> HC 25-11-2020 
     do ich=1,ng                                          !!>> HC 25-11-2020 
        do jch=1,ng                                       !!>> HC 25-11-2020 
           if(gen(ich)%t(jch).ne.0.0d0) numt=numt+1          !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
        enddo                                             !!>> HC 25-11-2020 
        do jch=1,gen(ich)%nww                             !!>> HC 25-11-2020 
           if(gen(ich)%r(jch,3)>0.0d0) numr=numr+1        !!>> HC 25-11-2020 Number of existing links in the r matrix (post-transcription reactions)
        enddo                                             !!>> HC 25-11-2020 
        do jch=2,nga                                      !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
           if(gen(ich)%e(jch).ne.0.0d0) nume=nume+1          !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
        enddo                                             !!>> HC 25-11-2020 
     enddo                                                !!>> HC 25-11-2020 
     
     probis=0.0d0;props=0                                 !!>> HC 25-11-2020 props is the total number of parameters that can change
     props=3*ng+numt+nume+numr+numk                       !!>> HC 25-11-2020 3*ng are diff, mu and mich of all genes
     probis(1)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mu
     probis(2)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability diff
     probis(3)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mich
     probis(4)=real(numt)/real(props)                     !!>> HC 25-11-2020 probability t matrix
     probis(5)=real(nume)/real(props)                     !!>> HC 25-11-2020 probability e matrix
     probis(6)=real(numr)/real(props)                     !!>> HC 25-11-2020 probability nww matrix == 0
     probis(7)=real(numk)/real(props)                     !!>> HC 25-11-2020 probability kadh matrix == 0
     cum=0.0d0                                            !!>> HC 25-11-2020
     call random_number(a)                                !!>> HC 25-11-2020
     do ich=1,size(probis)                                !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
        cum=cum+probis(ich)                               !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
        if (a<cum)then                                    !!>> HC 25-11-2020
           luckyprop=ich                                  !!>> HC 25-11-2020
           exit                                           !!>> HC 25-11-2020
        endif                                             !!>> HC 25-11-2020
     enddo                                                !!>> HC 25-11-2020
     mutacode=luckyprop                                   !!>> HC 16-9-2021  Store mutation code for output   
     if (luckyprop==1)then                                                !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DEGRADATION RATE (mu)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020
        if (i==0) i=1                                                     !!>> HC 25-11-2020 we randomly choose a gene i
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%mu                                                !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%mu==0.0d0) then                                        !!>> HC 25-11-2020
            gen(i)%mu=m_max*a                                             !!>> HC 25-11-2020 if mu=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
            inch=(1-2*a)*gen(i)%mu*mag_ismurate                           !!>> HC 27-11-2020 inch saves the magnitude of the change
            gen(i)%mu=gen(i)%mu+inch                                      !!>> HC 27-11-2020 more commonly, we increase or decrease mu randomly
            if (gen(i)%mu<0.0d0) gen(i)%mu=0.0d0                          !!>> HC 27-11-2020 negative values of mu are imposible
            if (gen(i)%mu>m_max)then                                      !!>> HC 27-11-2020 maximum degradation rate surpased, we apply the limit
               gen(i)%mu=gen(i)%mu-inch                                   !!>> HC 27-11-2020
               limit=0                                                    !!>> HC 27-11-2020
            endif                                                         !!>> HC 27-11-2020
        end if                                                            !!>> HC 27-11-2020
        newvalue=gen(i)%mu                                                !!>> HC 29-11-2021  we store the new value
     elseif (luckyprop==2)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
        if (i==0) i=1                                                     !!>> HC 25-11-2020
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%diffu                                             !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%diffu==0.0d0) then                                     !!>> HC 25-11-2020
           gen(i)%diffu=d_max*a                                           !!>> HC 25-11-2020 if diffu=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
           inch=(1-2*a)*gen(i)%diffu*mag_ismurate                         !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
           gen(i)%diffu=gen(i)%diffu+inch                                 !!>> HC 27-11-2020 more commonly, we increase or decrease diffu randomly
           if (gen(i)%diffu<0.0d0) gen(i)%diffu=0.0d0                     !!>> HC 27-11-2020 negative values of diffu are imposible
           if (gen(i)%diffu>d_max)then                                    !!>> HC 27-11-2020 maximum value of diffusion surpased, we apply the limit
              gen(i)%diffu=gen(i)%diffu-inch                              !!>> HC 27-11-2020
              limit=0                                                     !!>> HC 27-11-2020
           endif                                                          !!>> HC 27-11-2020
        end if                                                            !!>> HC 27-11-2020
        newvalue=gen(i)%diffu                                             !!>> HC 29-11-2021  we store the new value        
     elseif (luckyprop==3)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!
        call random_number(a)                                             !!>> HC 25-11-2020
        i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
        if (i==0) i=1                                                     !!>> HC 25-11-2020
        mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
        prevalue=gen(i)%mich                                              !!>> HC 29-11-2021 We store the old value
        call random_number(a)                                             !!>> HC 25-11-2020
        if (gen(i)%mich==0.0d0) then                                      !!>> HC 25-11-2020
           gen(i)%mich=mich_max*a                                         !!>> HC 25-11-2020 if mich=0.0d0, we change for a random proportion of maximum value
        else                                                              !!>> HC 25-11-2020
           inch=(1-2*a)*gen(i)%mich*mag_ismurate                          !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
           gen(i)%mich=gen(i)%mich+inch                                   !!>> HC 27-11-2020 more commonly, we increase or decrease mich randomly  
           if (gen(i)%mich<0.0d0) gen(i)%mich=0.0d0                       !!>> HC 27-11-2020 negative values of mich are imposible
           if (gen(i)%mich>mich_max)then                                  !!>> HC 27-11-2020 maximum michaelis-menten constant surpased, we apply the limit
              gen(i)%mich=gen(i)%mich-inch                                !!>> HC 27-11-2020
              limit=0                                                     !!>> HC 27-11-2020
           endif                                                          !!>> HC 27-11-2020
        end if                                                            !!>> HC 25-11-2020 
        newvalue=gen(i)%mich                                              !!>> HC 29-11-2021  we store the new value  
     elseif (luckyprop==4)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
        if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020  
        allocate(plinki(1:numt))                                          !!>> HC 25-11-2020
        if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
        allocate(plinkj(1:numt))                                          !!>> HC 25-11-2020
        ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020
        do ich=1,ng                                                       !!>> HC 25-11-2020
           do jch=1, ng                                                   !!>> HC 25-11-2020
              if (gen(ich)%t(jch).ne.0.0d0)then                              !!>> HC 25-11-2020
                 ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (from 1 to the number of iterations=numt)
                 plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that are transcription factors
                 plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of genes j whose transcription is regulated
              endif                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
        enddo                                                             !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        luckylink=ceiling(a*numt)                                         !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
        if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        i=plinki(luckylink)                                               !!>> HC 25-11-2020 gene i is the randomly selected transcription factor
        j=plinkj(luckylink)                                               !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
        mutageni=i                                                        !!>> HC 29-11-2021    
        mutagenj=j                                                        !!>> HC 29-11-2021
        prevalue=gen(i)%t(j)                                              !!>> HC 29-11-2021 store the previous value
        inch=(1-2*a)*gen(i)%t(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
        gen(i)%t(j)=gen(i)%t(j)+inch                                      !!>> HC 27-11-2020 we apply a mutation of magnitude inch in the link between i and j
        if (abs(gen(i)%t(j))>w_max)then                                   !!>> HC 27-11-2020 if the limit for transcription interactions is surpased
           gen(i)%t(j)=gen(i)%t(j)-inch                                   !!>> HC 27-11-2020 we go back an activate the flag (limit) to make a new mutation
           limit=0                                                        !!>> HC 27-11-2020
        endif                                                             !!>> HC 27-11-2020
        newvalue=gen(i)%t(j)                                              !!>> HC 29-11-2021
     elseif (luckyprop==5)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
        if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020
        allocate(plinki(1:nume))                                          !!>> HC 25-11-2020
        if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020
        allocate(plinkj(1:nume))                                          !!>> HC 25-11-2020
        ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
        do ich=1,ng                                                       !!>> HC 25-11-2020
           do jch=2,nga                                                   !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
              if (gen(ich)%e(jch).ne.0.0d0)then                           !!>> HC 25-11-2020
                 ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (1 to nume)
                 plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that regulate cell behaviors
                 plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of regulated cell behaviors
              endif                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
        enddo                                                             !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        luckylink=ceiling(a*nume)                                         !!>> HC 25-11-2020 Randomly selecting 
        if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        call random_number(a)                                             !!>> HC 25-11-2020
        i=plinki(luckylink)                                               !!>> HC 25-11-2020 Randomly selected gene i that regulates
        j=plinkj(luckylink)                                               !!>> HC 25-11-2020 Randomly selected behavior/property j
        mutageni=i                                                        !!>> HC 29-11-2021 store for output    
        mutagenj=j                                                        !!>> HC 29-11-2021 store for output
        prevalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output
        inch=(1-2*a)*gen(i)%e(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
        gen(i)%e(j)=gen(i)%e(j)+inch                                      !!>> HC 27-11-2020 we apply the mutation
        newval=gen(i)%e(j)                                                !!>> HC 27-11-2020 saving the new value for comparison
        if (newval>max_elim(j) .or. newval<min_elim(j))then               !!>> HC 27-11-2020 We compare with the limits of the cell behavior/property j
           gen(i)%e(j)=gen(i)%e(j)-inch                                   !!>> HC 27-11-2020 if the limit is surpased, we undo the change
           limit=0                                                        !!>> HC 27-11-2020 and activate the flag to make the IS mutation again
        endif                                                             !!>> HC 27-11-2020
        newvalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output  
     endif                                                                !!>> HC 25-11-2020
     
  else                      !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!!
     probs(1)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
     probs(2)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
     probs(3)=0d0           !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now)
     probs(4)=32d-2         !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
     probs(5)=2d-2          !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals
     call random_number(a)  !!>> HC 25-11-2020
     cum=0.0d0              !!>> HC 25-11-2020
     do ich=1,size(probs)   !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
        cum=cum+probs(ich)  !!>> HC 25-11-2020 taking into account probabilites in probs
        if (a<cum)then      !!>> HC 25-11-2020
           tind=ich         !!>> HC 25-11-2020
           exit             !!>> HC 25-11-2020
        endif               !!>> HC 25-11-2020
     enddo                  !!>> HC 25-11-2020
    if (tind .le. 2)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
       mutacode=4                                                         !!>> HC 29-11-2021 Store the matrix that will mutate
       do ich=1,ng                                                        !!>> HC 25-11-2020
         do jch=1,ng                                                      !!>> HC 25-11-2020
            if (gen(ich)%t(jch)==0) empty=empty+1                         !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
         enddo                                                            !!>> HC 25-11-2020 and the number of already existing links in t matrix
       enddo                                                              !!>> HC 25-11-2020
       full=ng*ng-empty                                                   !!>> HC 25-11-2020 
       if (tind==1 .and. full>0)then                                      !!>> HC 25-11-2020 !!!!!! DELECTION OF A LINK IN THE T MATRIX!!!!!!
          if (allocated(fullspoti)) deallocate(fullspoti)                 !!>> HC 25-11-2020
          allocate(fullspoti(1:full))                                     !!>> HC 25-11-2020
          if (allocated(fullspotj)) deallocate(fullspotj)                 !!>> HC 25-11-2020
          allocate(fullspotj(1:full))                                     !!>> HC 25-11-2020
          ord=0;fullspoti=0;fullspotj=0                                   !!>> HC 25-11-2020 These vectors store the already existing links
          do ich=1,ng                                                     !!>> HC 25-11-2020
            do jch=1,ng                                                   !!>> HC 25-11-2020
               if (gen(ich)%t(jch) .ne. 0.0d0)then                        !!>> HC 25-11-2020
                  ord=ord+1                                               !!>> HC 25-11-2020 order of vectors (1 to full)
                  fullspoti(ord)=ich                                      !!>> HC 25-11-2020 this gene i regulates the expression of
                  fullspotj(ord)=jch                                      !!>> HC 25-11-2020 gene j
               endif                                                      !!>> HC 25-11-2020
            enddo                                                         !!>> HC 25-11-2020
          enddo                                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          luckyg=ceiling(a*full)                                          !!>> HC 25-11-2020 Randomly pickin an interaction to delete
          if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          i=fullspoti(luckyg)                                             !!>> HC 25-11-2020 Removing interaction between gene i and j
          j=fullspotj(luckyg)                                             !!>> HC 25-11-2020
          mutageni=i                                                      !!>> HC 29-11-2021 Store for output
          mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
          prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
          gen(i)%t(j)=0.0d0                                               !!>> HC 25-11-2020
          newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 
       else                                                               !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
          if (allocated(freespoti)) deallocate(freespoti)                 !!>> HC 25-11-2020
          allocate(freespoti(1:empty))                                    !!>> HC 25-11-2020
          if (allocated(freespotj)) deallocate(freespotj)                 !!>> HC 25-11-2020
          allocate(freespotj(1:empty))                                    !!>> HC 25-11-2020
          ord=0;freespoti=0;freespotj=0                                   !!>> HC 25-11-2020 These vectors store the free spots
          do ich=1,ng                                                     !!>> HC 25-11-2020
            do jch=1,ng                                                   !!>> HC 25-11-2020
               if (gen(ich)%t(jch)==0)then                                !!>> HC 25-11-2020
                  ord=ord+1                                               !!>> HC 25-11-2020 order of the vectors (1 to empty)
                  freespoti(ord)=ich                                      !!>> HC 25-11-2020 gene i has NOT and interaction
                  freespotj(ord)=jch                                      !!>> HC 25-11-2020 with gene j
               endif                                                      !!>> HC 25-11-2020
            enddo                                                         !!>> HC 25-11-2020
          enddo                                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          luckyg=ceiling(a*empty)                                         !!>> HC 25-11-2020 Randomly picking a spot
          if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020
          call random_number(a)                                           !!>> HC 25-11-2020
          call random_number(b)                                           !!>> HC 25-11-2020
          if(b<0.50d0)then;signo=1;else;signo=-1;endif                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
          i=freespoti(luckyg)                                             !!>> HC 25-11-2020
          j=freespotj(luckyg)                                             !!>> HC 25-11-2020 Making new interaction
          mutageni=i                                                      !!>> HC 29-11-2021 Store for output
          mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
          prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
          gen(i)%t(j)=a*signo*w_max*t_mag                                 !!>> HC 31-5-2023
          newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 
       endif                                                              !!>> HC 25-11-2020
    elseif (tind==3)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
    !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW        !!>> HC 25-11-2020
    elseif (tind==4)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e matrix)!!!!!!
           empty=0                                                        !!>> HC 25-11-2020
           do ich=1,ng                                                    !!>> HC 25-11-2020
              do jch=1,nga                                                !!>> HC 25-11-2020
                 if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                 if (gen(ich)%e(jch)==0.0d0) then                         !!>> HC 27-11-2020
                     empty=empty+1                                        !!>> HC 27-11-2020 Number of empty spots in e Matrix
                 else                                                     !!>> HC 27-11-2020
                    full=full+1                                           !!>> HC 27-11-2020 Number of existing interactions in e matrix
                 endif                                                    !!>> HC 27-11-2020
              enddo                                                       !!>> HC 25-11-2020
           enddo                                                          !!>> HC 25-11-2020
           call random_number(a)                                          !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
           if (a<0.50d0 .and. full>0.or.empty==0)then                     !!>> HC 10-5-2022 !!!REMOVE INTERACTION           
              if (allocated(fullspoti)) deallocate(fullspoti)             !!>> HC 25-11-2020
              allocate(fullspoti(1:full))                                 !!>> HC 25-11-2020
              if (allocated(fullspotj)) deallocate(fullspotj)             !!>> HC 25-11-2020
              allocate(fullspotj(1:full))                                 !!>> HC 25-11-2020
              ord=0;fullspoti=0;fullspotj=0                               !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
              do ich=1,ng                                                 !!>> HC 25-11-2020
                 do jch=2,nga                                             !!>> HC 25-11-2020 we start in 2 because 1 is being an adhesion mol
                    if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
                    if (gen(ich)%e(jch).ne. 0.0d0)then                    !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                       ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                       fullspoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                       fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
                    endif                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              enddo                                                       !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              luckyg=ceiling(a*real(full))                                !!>> HC 25-11-2020 Randomly picking an interaction to kill
              i=fullspoti(luckyg)                                         !!>> HC 25-11-2020
              j=fullspotj(luckyg)                                         !!>> HC 25-11-2020
                                                                          !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
              mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
              mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
              mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected
              prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
              newvalue=0.0d0                                              !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
              gen(i)%e(j)=0.0d0                                           !!>> HC 25-11-2020 we kill it easily
                          
           else                                                           !!>> HC 25-11-2020 !ADD AN INTERACTION   
              if (allocated(freespoti)) deallocate(freespoti)             !!>> HC 25-11-2020
              allocate(freespoti(1:empty))                                !!>> HC 25-11-2020
              if (allocated(freespotj)) deallocate(freespotj)             !!>> HC 25-11-2020
              allocate(freespotj(1:empty))                                !!>> HC 25-11-2020
              ord=0;freespoti=0;freespotj=0;                              !!>> HC 25-11-2020 
              do ich=1,ng                                                 !!>> HC 25-11-2020
                 do jch=1,nga                                             !!>> HC 25-11-2020
                    if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                    if (gen(ich)%e(jch)==0.0d0)then                       !!>> HC 25-11-2020
                       ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                       freespoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                       freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
                    endif                                                 !!>> HC 25-11-2020
                 enddo                                                    !!>> HC 25-11-2020
              enddo                                                       !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              luckyg=ceiling(a*empty)                                     !!>> HC 25-11-2020 Randomly picking an interaction
              if(luckyg==0)luckyg=1                                       !!>> HC 25-11-2020
              i=freespoti(luckyg)                                         !!>> HC 25-11-2020
              j=freespotj(luckyg)                                         !!>> HC 25-11-2020 SAVE THE INFORMATION IN THE OUTPUT --START--
              mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
              mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
              mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected
              prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
              call random_number(c)                                       !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
              anch=max_elim(j)-min_elim(j)                                !!>> HC 28-11-2020 we introduce an interaction in the range of the limits            
              gen(i)%e(j)=anch*c+min_elim(j)                              !!>> HC 31-5-2023 
              gen(i)%e(j)=gen(i)%e(j)*t_mag                               !!>> HC 31-5-2023  Correct the magnitude of the new interaction
              if (j==nparam_per_node+8)then                               !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation
                 gen(ng)%e(nparam_per_node+17)=0.20d0                     !!>> HC 12-11-2020 !!!ATENTION!! gene ng is homogeneously expressed in the blastula ics 
              endif                                                       !!>> HC 7-12-2020 !!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
              newvalue=gen(i)%e(j)                                        !!>> HC 29-11-2020  Store the new value                 
           endif                                                          !!>> HC 25-11-2020
    elseif (tind==5)then                                                  !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
           mutacode=7                                                     !!>> HC 29-11-2021 Store that this is a gene duplication/delection
           call random_number(a)                                          !!>> HC 10-9-2021
           i=ceiling(a*(ng))                                              !!>> HC 10-9-2021
           if (i==0) i=1                                                  !!>> HC 25-11-2020
           mutageni=i                                                     !!>> HC 16-9-2021 save the chosen gene for the output
           if (ng>2)then                                                  !!>> HC 25-11-2020
              call random_number(a)                                       !!>> HC 25-11-2020
              if (a<0.50d0)then                                           !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
                 mutagenj=-1                                              !!>> HC 29-11-2021 store that this is a delection
                 call deletion(i)                                         !!>> HC 25-11-2020
              else                                                        !!>> HC 25-11-2020
                 mutagenj=1                                               !!>> HC 29-11-2021 store that this is a duplication
                 call duplication(i)                                      !!>> HC 25-11-2020
              endif                                                       !!>> HC 25-11-2020
           else                                                           !!>> HC 25-11-2020
              mutagenj=1                                                  !!>> HC 29-11-2021 store that this is a duplication
              call duplication(i)                                         !!>> HC 25-11-2020
           endif                                                          !!>> HC 25-11-2020
    endif                                                                 !!>> HC 25-11-2020
  endif                                                                   !!>> HC 25-11-2020
     
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 
!!>> HC 16-9-2021 !! ATENTION !! We are not using pre/post forms right now, this should be changed if you want them!!      !!>> HC 16-9-2021
!  if (allocated(npag)) deallocate(npag)                             !!>> HC 25-11-2020  This last chunk of code was inherited from muta
!  allocate(npag(nga))                                               !!>> HC 25-11-2020  it realocates matrixes 
!  npag=0                                                            !!>> HC 25-11-2020

!  if (allocated(whonpag)) deallocate(whonpag)                       !!>> HC 25-11-2020
!  allocate(whonpag(nga,ng))                                         !!>> HC 25-11-2020
!  whonpag=0                                                         !!>> HC 25-11-2020

!  call update_npag                                                  !!>> HC 25-11-2020

!  do i=1,ng                                                         !!>> HC 25-11-2020
!     do j=1,gen(i)%npre                                             !!>> HC 25-11-2020
!        if (gen(i)%pre(j)==0) then                                  !!>> HC 25-11-2020
!           print *,i,j,"aqui hi ha un zero a pre"                   !!>> HC 25-11-2020
!           stop                                                     !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020

!     do j=1,gen(i)%npost                                            !!>> HC 25-11-2020
!        if (gen(i)%post(j)==0) then                                 !!>> HC 25-11-2020
!        print *,i,j,"aqui hi ha un zero a post"                     !!>> HC 25-11-2020
!        stop                                                        !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020
!  end do                                                            !!>> HC 25-11-2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 

  inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
  if (ng<1) then                                                    !!>> HC 25-11-2020
     print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
     inviable=1                                                     !!>> HC 25-11-2020
  else                                                              !!>> HC 25-11-2020
     do j=1,ng                                                      !!>> HC 25-11-2020
        if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
           inviable=0                                               !!>> HC 25-11-2020
           goto 37                                                  !!>> HC 25-11-2020
        end if                                                      !!>> HC 25-11-2020
     end do                                                         !!>> HC 25-11-2020
     print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
     inviable=1                                                     !!>> HC 25-11-2020
  end if                                                            !!>> HC 25-11-2020
37 continue                                                         !!>> HC 25-11-2020
       
end subroutine suremuta_no_kadh

!***********************************************************************************************************!
!add argument to function so if no functional part of subroutine is mutated then do not run development
subroutine suremuta_no_kadh_functional_net(limit, inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue,rangfile,runornot,&
pcpval) !!>> AL 29-11-24
   !!>> AL 9-4-24: this function determines whether the mutation alters the functional part of the gen network of an individual or not.
   !!>> AL 9-4-24: if not then development of mutated individual is not runned and fitness is inherited from parent.
   !!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file 
   !!>> HC 14-4-2023 This version is the same but ignoring specific adhesion molecules
   !!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
   !!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
   !!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
   !!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
   !!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
   implicit none

   integer             :: limit,inviable,mind,tind
   real*8, allocatable :: copkadh(:,:)
   real*8              :: probmod,xx !the modifier used to bias T mutations of W towards interaction loss, and the temp var to store the ratio between empty and filled w positions (for readability only)
   integer             ::genadh 
   integer             :: props,nume,numt,numr,numk !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
   integer             :: luckyprop, luckylink, luckyg  !!>> HC 24-11-2020 selected gene, interaction or property
   integer             :: ich,jch,lch,ord,ord2  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
   integer                     :: full,empty !!>> HC 24-11-2020 full an empty spaces for T mutations
   real*8, dimension(1:5)      :: probs !!>> HC 24-11-2020 relative probabilities of different T mutations
   !real*8, dimension(1:7)      :: probis !!>> HC 24-11-2020 relative probabilities of different IS mutations
   real*8, dimension(5)        :: probis  !AL 2-12-24 this is 5 since there are no reactions 
   real*8, dimension(1:nga)    :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 

   integer, dimension(1:nga)   :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
   real*8                      :: cum !!>> HC 24-11-2020 funny position to store cumulative Fi
   integer                     :: stilladi, stilladj,isad, mnadh, mxadh !!>> HC 24-11-2020 pointers related to T mutations in the Kadh matrix
   integer                     :: signo  !!>> HC 24-11-2020 sign of the mutation  
   real*8                      :: inch, newval, anch,pcpval !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
   real*8, dimension(1:5)      :: bhc !!>>HC 10-9-2021 Patch to a problem with biased random numbers
   integer                     :: mutacode,mutageni,mutagenj  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
   real*8                      :: prevalue, newvalue
   character*400               :: rangfile !!>> 6-10-2021 file where the ranges are stored
   real*8, dimension (1:5)     :: min_glim, max_glim                                  !!>>HC 6-10-2021
   real*8                      :: t_mag,c !!>> HC 31-5-2023 Magnitude of the initial value in T mutations (should be a free parameter)
   integer,allocatable,dimension(:)  :: plinki,plinkj !!>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
   integer,allocatable,dimension(:)  :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM
   integer,allocatable,dimension(:)  :: isadh !!>> HC 24-11-2020 if 1 the interaction is in the kadh matrix (TM)
   logical :: file_exists
   !*******functional subnetwork*********! !!>>AL 9-4-24
   integer :: runornot                          !!>>AL 10-4-24: if == 1 you run development, if 0 you do not          
   integer, allocatable, dimension(:) :: good   !!>>AL 10-4-24: this vector has elements = 1 when genes are functional
   if(allocated(good)) deallocate(good)         !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
   allocate(good(1:ng))                 
   !! tal vez necesitamos vector despues de mutaciones porque puede ser que red functional cambie.
   !*************************************! !!>>AL 9-4-24
   
   !**************************limits and ranges **************************!        
   limit=1                                               !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
   rembeh=0; max_elim=0.0d0; min_elim=0.0d0              !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
   max_glim=0.0d0; min_glim=0.0d0                        !!>> HC 6-10-2020 Same for other gene properties 
   call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
   d_max=max_glim(1)    ; d_min=min_glim(1)              !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020 !>>AL 28-1-25 this is redundant 
   m_max=max_glim(2)    ; m_min=min_glim(2)              !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
   mich_max=max_glim(3) ; mich_min=min_glim(3)           !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
   w_max=max_glim(4)    ; w_min=min_glim(4)              !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
   b_max=max_glim(5)    ; b_min=min_glim(5)              !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

   call functional_net(good)                             !!>>AL 29-11-24

   if (allocated(dupi)) deallocate(dupi)
   allocate(dupi(ng*ng))       !since each gene can no duplicate more than ng times per mutation round

   if (allocated(father_son)) deallocate(father_son)
   allocate(father_son(ng*ng)) !since each gene can no duplicate more than ng times per mutation round

   if (allocated(todel)) deallocate(todel)
   allocate(todel(ng*ng*10))   !!>>AL 9-4-24: why not only ng*ng..?

   dupi=0
   father_son=0
   todel=0
   empty=0; full=0             !!>> HC 24-11-2020
   inch=0.0d0; anch=0.0d0      !!>> HC 28-11-2020
   mutacode=0                  !!>> HC 16-9-2021  
   mutageni=0                  !!>> HC 29-11-2021
   mutagenj=0                  !!>> HC 29-11-2021
   prevalue=0.0d0              !!>> HC 29-11-2021
   newvalue=0.0d0              !!>> HC 29-11-2021
   t_mag=0.20d0                !!>> HC 31-5-2023
   pcpval = 0.0d0                               !!>>AL 28-12-24:

   !call random_number(bhc)     !!>> HC 25-11-2020 THESE ARE BIASED I DO NOT KNOW WHY !!>>AL 10-4-24 bhc ni si quiera se usa... borrar(?)
   
   if (mind==1)then                                        !!>> HC 25-11-2020 !******************* IS MUTATIONS***************************!
      nume=0; numt=0; numk=0; numr=0                       !!>> HC 25-11-2020 
      do ich=1,ng                                          !!>> HC 25-11-2020 
         do jch=1,ng                                       !!>> HC 25-11-2020 
            if(gen(ich)%t(jch).ne.0.0d0) numt=numt+1       !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
         enddo                                             !!>> HC 25-11-2020 
         ! do jch=1,gen(ich)%nww                             !!>> HC 25-11-2020 
         !    if(gen(ich)%r(jch,3)>0.0d0) numr=numr+1        !!>> HC 25-11-2020 Number of existing links in the r matrix (post-transcription reactions)
         ! enddo                                             !!>> HC 25-11-2020 
         do jch=2,nga                                      !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
            if(gen(ich)%e(jch).ne.0.0d0) nume=nume+1       !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
         enddo                                             !!>> HC 25-11-2020 
      enddo                                                !!>> HC 25-11-2020 
      
      probis=0.0d0;props=0                                 !!>> HC 25-11-2020 props is the total number of parameters that can change
      
      !props=3*ng+numt+nume+numr+numk                       !!>> HC 25-11-2020 3*ng are diff, mu and mich of all genes
      
      props=3*ng+numt+nume                           !AL 2-12-24 

      probis(1)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mu
      probis(2)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability diff
      probis(3)=real(ng)/real(props)                       !!>> HC 25-11-2020 probability mich
      probis(4)=real(numt)/real(props)                     !!>> HC 25-11-2020 probability t matrix
      probis(5)=real(nume)/real(props)                     !!>> HC 25-11-2020 probability e matrix

      !probis(6)=real(numr)/real(props)                     !!>> HC 25-11-2020 probability nww matrix == 0
      !probis(7)=real(numk)/real(props)                     !!>> HC 25-11-2020 probability kadh matrix == 0
     
      cum=0.0d0                                            !!>> HC 25-11-2020
      call random_number(a)                                !!>> HC 25-11-2020
      do ich=1,size(probis)                                !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
         cum=cum+probis(ich)                               !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
         if (a<cum)then                                    !!>> HC 25-11-2020
            luckyprop=ich                                  !!>> HC 25-11-2020
            exit                                           !!>> HC 25-11-2020
         endif                                             !!>> HC 25-11-2020
      enddo                                                !!>> HC 25-11-2020

      !luckyprop = 5
      mutacode=luckyprop                                   !!>> HC 16-9-2021  Store mutation code for output   
      
      if (luckyprop==1)then                                                !!>> HC 25-11-2020 !***********IS MUTATION IN DEGRADATION RATE (mu)**************!
         
         call random_number(a)                                             !!>> HC 25-11-2020
         i=ceiling(a*ng)                                                   !!>> HC 25-11-2020
         if (i==0) i=1                                                     !!>> HC 25-11-2020 we randomly choose a gene i
         
         mutageni=i                                                        !!>> HC 29-11-2021 We store the gene for the output

         if(good(mutageni) == 1)then; runornot = 1; end if                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%mu                                                !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         
         if (gen(i)%mu==0.0d0) then                                        !!>> HC 25-11-2020
            inch=m_max*a 
         else 
            inch=(1-2*a)*gen(i)%mu*mag_ismurate                            !!>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if ((gen(i)%mu+inch) < 0.0d0)then                              !!>> ALI 2-12-24 we always want mutations to change the selected parameters
            if(gen(i)%mu == m_min)then
               inch = -inch                                         
               gen(i)%mu = gen(i)%mu+inch  
            else
               gen(i)%mu = m_min 
            end if
         else if ((gen(i)%mu+inch) > m_max) then 
            if(gen(i)%mu == m_max)then 
               inch = -inch
               gen(i)%mu=gen(i)%mu+inch
            else 
               gen(i)%mu = m_max
            end if
         else 
            gen(i)%mu=gen(i)%mu+inch
         end if

         newvalue=gen(i)%mu                                               
      
      elseif (luckyprop==2)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!
         
         call random_number(a)                                             !!>> HC 25-11-2020
         i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
         if (i==0) i=1                                                     !!>> HC 25-11-2020
         
         mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
         
         if(good(mutageni) == 1)then; runornot = 1; end if                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%diffu                                             !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(i)%diffu==0.0d0) then                                        !!>> HC 25-11-2020
            inch=d_max*a 
         else 
            inch=(1-2*a)*gen(i)%diffu*mag_ismurate                            !!>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if ((gen(i)%diffu+inch) < 0.0d0)then                              !!>> AL 2-12-24
            if(gen(i)%diffu == d_min)then
               inch = -inch                                         
               gen(i)%diffu = gen(i)%diffu+inch  
            else
               gen(i)%diffu = d_min 
            end if
         else if ((gen(i)%diffu+inch) > d_max) then 
            if(gen(i)%diffu == d_max)then 
               inch = -inch
               gen(i)%diffu=gen(i)%diffu+inch
            else 
               gen(i)%diffu = d_max
            end if
         else 
            gen(i)%diffu = gen(i)%diffu + inch
         end if
         
         newvalue=gen(i)%diffu 

      elseif (luckyprop==3)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!
         
         call random_number(a)                                             !!>> HC 25-11-2020
         i=ceiling(a*ng)                                                   !!>> HC 25-11-2020 we randomly choose a gene i
         if (i==0) i=1                                                     !!>> HC 25-11-2020
         
         mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
         
         if(good(mutageni) == 1)then; runornot = 1; end if                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%mich                                              !!>> HC 29-11-2021 We store the old value
         
        call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(i)%mich==0.0d0) then                                     !!>> AL 3-12-24 
            inch=mich_max*a 
         else 
            inch=(1-2*a)*gen(i)%mich*mag_ismurate                         
         end if
         
         if ((gen(i)%mich+inch) < 0.0d0)then                              !!>> AL 2-12-24
            if(gen(i)%mich == mich_min)then
               inch = -inch                                         
               gen(i)%mich = gen(i)%mich + inch  
            else
               gen(i)%mich = mich_min 
            end if
         else if ((gen(i)%mich + inch) > mich_max) then 
            if(gen(i)%mich == mich_max)then 
               inch = -inch
               gen(i)%mich = gen(i)%mich + inch
            else 
               gen(i)%mich = mich_max
            end if
         else 
            gen(i)%mich = gen(i)%mich + inch
         end if

         newvalue=gen(i)%mich                                              !!>> HC 29-11-2021  we store the new value  
      
      elseif (luckyprop==4)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020  
         allocate(plinki(1:numt))                                          !!>> HC 25-11-2020
         if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
         allocate(plinkj(1:numt))                                          !!>> HC 25-11-2020
         ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020
         
         do ich=1,ng                                                       !!>> HC 25-11-2020
            do jch=1,ng                                                    !!>> HC 25-11-2020
               if (gen(ich)%t(jch).ne.0.0d0)then                           !!>> HC 25-11-2020
                  ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (from 1 to the number of iterations=numt)
                  plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that are transcription factors
                  plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of genes j whose transcription is regulated
               endif                                                       !!>> HC 25-11-2020
            enddo                                                          !!>> HC 25-11-2020
         enddo                                                             !!>> HC 25-11-2020
         
         call random_number(a)                                             !!>> HC 25-11-2020
         luckylink=ceiling(a*numt)                                         !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
         if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        
         if(good(plinki(luckylink)) == 1)then; runornot = 1; end if        !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
                                                                           !!>> AL 10-4-24: if plinki(luckylink) then plinkj is functional (by definition)
        
         i=plinki(luckylink)                                               !!>> HC 25-11-2020 gene i is the randomly selected transcription factor
         j=plinkj(luckylink)                                               !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
         mutageni=i                                                        !!>> HC 29-11-2021    
         mutagenj=j                                                        !!>> HC 29-11-2021
                  
         prevalue=gen(i)%t(j)                                              !!>> HC 29-11-2021 store the previous value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations
        
         inch=(1-2*a)*gen(i)%t(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
         if ((gen(i)%t(j) + inch) < w_min)then                              !!>> AL 2-12-24
            if(gen(i)%t(j) == w_min)then
               inch = -inch                                         
               gen(i)%t(j) = gen(i)%t(j) + inch  
            else
               gen(i)%t(j) = w_min
            end if
         else if ((gen(i)%t(j) + inch) > w_max) then 
            if(gen(i)%t(j) == w_max)then
               inch = -inch                                         
               gen(i)%t(j) = gen(i)%t(j) + inch  
            else
               gen(i)%t(j) = w_max
            end if
         else 
            gen(i)%t(j) = gen(i)%t(j) + inch
         end if
         
         newvalue=gen(i)%t(j)                                              !!>> HC 29-11-2021

      elseif (luckyprop==5)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
         
         if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020
         allocate(plinki(1:nume))                                          !!>> HC 25-11-2020
         if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020
         allocate(plinkj(1:nume))                                          !!>> HC 25-11-2020
         
         ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
         
         do ich=1,ng                                                       !!>> HC 25-11-2020
            do jch=2,nga                                                   !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
               if (gen(ich)%e(jch).ne.0.0d0)then                           !!>> HC 25-11-2020
                  ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (1 to nume)
                  plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that regulate cell behaviors
                  plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of regulated cell behaviors
               endif                                                       !!>> HC 25-11-2020
            enddo                                                          !!>> HC 25-11-2020
         enddo                                                             !!>> HC 25-11-2020

         call random_number(a)                                             !!>> HC 25-11-2020
         luckylink=ceiling(a*nume)                                         !!>> HC 25-11-2020 Randomly selecting 
         if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        
         runornot = 1
         
         i=plinki(luckylink)                                               !!>> HC 25-11-2020 Randomly selected gene i that regulates
         j=plinkj(luckylink)                                               !!>> HC 25-11-2020 Randomly selected behavior/property j
         mutageni=i                                                        !!>> HC 29-11-2021 store for output    
         mutagenj=j                                                        !!>> HC 29-11-2021 store for output
         prevalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         inch = (1-2*a)*gen(i)%e(j)*mag_ismurate                           !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
         
         if ((gen(i)%e(j) + inch) < min_elim(j))then                       !!>> AL 2-12-24
            if(gen(i)%e(j) == min_elim(j))then
               inch = -inch                                         
               gen(i)%e(j) = gen(i)%e(j) + inch  
            else
               gen(i)%e(j) = min_elim(j)
            end if
         else if ((gen(i)%e(j) + inch) > max_elim(j)) then 
            if(gen(i)%e(j) == max_elim(j))then 
               inch = -inch
               gen(i)%e(j) = gen(i)%e(j) + inch  
            else 
               gen(i)%e(j) = max_elim(j)
            end if
         else 
            gen(i)%e(j) = gen(i)%e(j) + inch
         end if

         newvalue=gen(i)%e(j)            
      endif                                                                !!>> HC 25-11-2020

   else                      !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!!
      probs(1)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
      probs(2)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
      probs(3)=0d0           !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now)
      probs(4)=32d-2         !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
      probs(5)=2d-2          !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals
      
      call random_number(a)  !!>> HC 25-11-2020
      cum=0.0d0              !!>> HC 25-11-2020
      do ich=1,size(probs)   !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
         cum=cum+probs(ich)  !!>> HC 25-11-2020 taking into account probabilites in probs
         if (a<cum)then      !!>> HC 25-11-2020
            tind=ich         !!>> HC 25-11-2020
            exit             !!>> HC 25-11-2020
         endif               !!>> HC 25-11-2020
      enddo                  !!>> HC 25-11-2020
      
      if (tind .le. 2)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         mutacode=4                                                         !!>> HC 29-11-2021 Store the matrix that will mutate
         
         do ich=1,ng                                                        !!>> HC 25-11-2020
            do jch=1,ng                                                      !!>> HC 25-11-2020
               if (gen(ich)%t(jch)==0) empty=empty+1                         !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
            enddo                                                            !!>> HC 25-11-2020 and the number of already existing links in t matrix
         enddo                                                              !!>> HC 25-11-2020
         
         full=ng*ng-empty                                                   !!>> HC 25-11-2020 
         
         if (tind==1.and.full>0)then                                      !!>> HC 25-11-2020 !!!!!! DELETION OF A LINK IN THE T MATRIX!!!!!!
            
            if (allocated(fullspoti)) deallocate(fullspoti)                 !!>> HC 25-11-2020
            allocate(fullspoti(1:full))                                     !!>> HC 25-11-2020
            if (allocated(fullspotj)) deallocate(fullspotj)                 !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                                     !!>> HC 25-11-2020
            
            ord=0;fullspoti=0;fullspotj=0                                   !!>> HC 25-11-2020 These vectors store the already existing links
           
            do ich=1,ng                                                     !!>> HC 25-11-2020
               do jch=1,ng                                                   !!>> HC 25-11-2020
                  if (gen(ich)%t(jch).ne.0.0d0)then                        !!>> HC 25-11-2020
                     ord=ord+1                                               !!>> HC 25-11-2020 order of vectors (1 to full)
                     fullspoti(ord)=ich                                      !!>> HC 25-11-2020 this gene i regulates the expression of
                     fullspotj(ord)=jch                                      !!>> HC 25-11-2020 gene j
                  endif                                                      !!>> HC 25-11-2020
               enddo                                                         !!>> HC 25-11-2020
            enddo                                                           !!>> HC 25-11-2020
            
            call random_number(a)                                           !!>> HC 25-11-2020
            luckyg=ceiling(a*full)                                          !!>> HC 25-11-2020 Randomly pickin an interaction to delete
            if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020

            if(good(fullspoti(luckyg)) == 1)then; runornot = 1; end if        !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

            i=fullspoti(luckyg)                                             !!>> HC 25-11-2020 Removing interaction between gene i and j
            j=fullspotj(luckyg)                                             !!>> HC 25-11-2020
            mutageni=i                                                      !!>> HC 29-11-2021 Store for output
            mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
            
            prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
            gen(i)%t(j)=0.0d0                                               !!>> HC 25-11-2020
            newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 

         else                                                               !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
            if (allocated(freespoti)) deallocate(freespoti)                 !!>> HC 25-11-2020
            allocate(freespoti(1:empty))                                    !!>> HC 25-11-2020
            if (allocated(freespotj)) deallocate(freespotj)                 !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                                    !!>> HC 25-11-2020
            
            ord=0;freespoti=0;freespotj=0                                   !!>> HC 25-11-2020 These vectors store the free spots
            
            do ich=1,ng                                                     !!>> HC 25-11-2020
               do jch=1,ng                                                   !!>> HC 25-11-2020
                  if (gen(ich)%t(jch)==0)then                                !!>> HC 25-11-2020
                     ord=ord+1                                               !!>> HC 25-11-2020 order of the vectors (1 to empty)
                     freespoti(ord)=ich                                      !!>> HC 25-11-2020 gene i has NOT and interaction
                     freespotj(ord)=jch                                      !!>> HC 25-11-2020 with gene j
                  endif                                                      !!>> HC 25-11-2020
               enddo                                                         !!>> HC 25-11-2020
            enddo                                                           !!>> HC 25-11-2020
            
            call random_number(a)                                           !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                         !!>> HC 25-11-2020 Randomly picking a spot
            if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020

             call random_number(a)                                           !!>> HC 25-11-2020
            if(a < 1.0d-8) a = 0.01                                         !!>> AL 2-12-24 we don't want neutral mutations
            call random_number(b)                                           !!>> HC 25-11-2020
            if(b<0.50d0)then;signo=1;else;signo=-1;endif                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
        
            i=freespoti(luckyg)                                             !!>> HC 25-11-2020
            j=freespotj(luckyg)                                             !!>> HC 25-11-2020 Making new interaction
            mutageni=i                                                      !!>> HC 29-11-2021 Store for output
            mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
            
            prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
            gen(i)%t(j)=a*signo*w_max*t_mag                                 !!>> HC 31-5-2023
            newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 
            
            call functional_net(good)                                       !AL 29-11-24 important
            if(good(freespoti(luckyg)) == 1) runornot = 1                   !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

         endif                                                              !!>> HC 25-11-2020
      
      elseif (tind==3)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
      !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW        !!>> HC 25-11-2020
      
      elseif (tind==4)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e matrix)!!!!!!
         empty=0                                                        !!>> HC 25-11-2020
         
         do ich=1,ng                                                    !!>> HC 25-11-2020
            do jch=2,nga                                                !!>> AL 2-12-24
               if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
               if (gen(ich)%e(jch)==0.0d0) then                         !!>> HC 27-11-2020
                  empty=empty+1                                        !!>> HC 27-11-2020 Number of empty spots in e Matrix
               else                                                     !!>> HC 27-11-2020
                  full=full+1                                           !!>> HC 27-11-2020 Number of existing interactions in e matrix
               endif                                                    !!>> HC 27-11-2020
            enddo                                                       !!>> HC 25-11-2020
         enddo                                                          !!>> HC 25-11-2020
         
         call random_number(a)                                          !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
         if ((a<0.50d0 .and. full>1) .or. empty==0)then                 !!>> HC 10-5-2022 !!!REMOVE INTERACTION !!>> AL 28-1-25 important: we dont do negative T mutation if we only have 1 cell beh. regualted.      
            runornot = 1                                                !!>> AL 10-4-24: modifying a cell behaviour always means altering functional newtwork 
                                                                        !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(fullspoti)) deallocate(fullspoti)             !!>> HC 25-11-2020
            allocate(fullspoti(1:full))                                 !!>> HC 25-11-2020
            if (allocated(fullspotj)) deallocate(fullspotj)             !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                                 !!>> HC 25-11-2020
            
            ord=0;fullspoti=0;fullspotj=0                               !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
            
            do ich=1,ng                                                 !!>> HC 25-11-2020
               do jch=2,nga                                             !!>> HC 25-11-2020 we start in 2 because 1 is being an adhesion mol
                  if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
                  if (gen(ich)%e(jch).ne. 0.0d0)then                    !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                     ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                     fullspoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                     fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
                  endif                                                 !!>> HC 25-11-2020
               enddo                                                    !!>> HC 25-11-2020
            enddo                                                       !!>> HC 25-11-2020
            
            call random_number(a)                                       !!>> HC 25-11-2020
            luckyg=ceiling(a*real(full))                                !!>> HC 25-11-2020 Randomly picking an interaction to kill
            if(luckyg==0)luckyg=1                                       !!>> AL 3-12-2024   
            
            i=fullspoti(luckyg)                                         !!>> HC 25-11-2020
            j=fullspotj(luckyg)                                         !!>> HC 25-11-2020
            !if(i==24 .and. j==37)then                                   !!>> AL 28-1-25 important: we always have cell division consitutively being regulated by gene 24 !important: this is wrong
            !   limit = 0
            !end if 
                                                                        !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected

            prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
            newvalue=0.0d0                                              !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
            gen(i)%e(j)=0.0d0                                           !!>> HC 25-11-2020 we kill it easily
                        
         else                                                           !!>> HC 25-11-2020 !ADD AN INTERACTION   
            runornot = 1                                                !!>> AL 10-4-24: modifying a cell behaviour always means altering functional newtwork 
                                                                        !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(freespoti)) deallocate(freespoti)             !!>> HC 25-11-2020
            allocate(freespoti(1:empty))                                !!>> HC 25-11-2020
            if (allocated(freespotj)) deallocate(freespotj)             !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                                !!>> HC 25-11-2020
            
            ord=0;freespoti=0;freespotj=0;                              !!>> HC 25-11-2020 
            
            do ich=1,ng                                                 !!>> HC 25-11-2020
               do jch=2,nga                                             !!>> HC 25-11-2020
                  if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                  if (gen(ich)%e(jch)==0.0d0)then                       !!>> HC 25-11-2020
                     ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                     freespoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                     freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
                  endif                                                 !!>> HC 25-11-2020
               enddo                                                    !!>> HC 25-11-2020
            enddo                                                       !!>> HC 25-11-2020
            
            call random_number(a)                                       !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                     !!>> HC 25-11-2020 Randomly picking an interaction
            if(luckyg==0)luckyg=1                                       !!>> HC 25-11-2020
            
            i=freespoti(luckyg)                                         !!>> HC 25-11-2020
            j=freespotj(luckyg)                                         !!>> HC 25-11-2020 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
            
            call random_number(c)                                       !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
            if(c < 1.0d-8) c = 0.01            
            
            anch=max_elim(j)-min_elim(j)                                !!>> HC 28-11-2020 we introduce an interaction in the range of the limits            
            
            gen(i)%e(j)=anch*c+min_elim(j)                              !!>> HC 31-5-2023 
            gen(i)%e(j)=gen(i)%e(j)*t_mag                               !!>> HC 31-5-2023  Correct the magnitude of the new interaction
            
            !if (j==nparam_per_node+8)then                               !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation
            !   gen(ng)%e(nparam_per_node+17)=0.20d0                     !!>> HC 12-11-2020 !!!ATENTION!! gene ng is homogeneously expressed in the blastula ics 
            !endif                                                       !!>> HC 7-12-2020 !!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
            
            if (j==nparam_per_node+8)then                               !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation !>> AL 2-12-24 nparam_per_node = 35 PCP is the papers called PCO for whatever reason 
               call random_number(c)                                    !!>> AL 2-12-24 
               if(c < 1.0d-8) c = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations
              
               gen(ng)%e(nparam_per_node+17) = 0.20d0*c                 !!>> AL 2-12-24 we assign a random magnitud with respecto to max val (0.20d0)
            endif              
            pcpval = gen(ng)%e(nparam_per_node+17)                      !!>> AL 10-1-25 
            newvalue=gen(i)%e(j)                                        !!>> HC 29-11-2020  Store the new value          
         endif                                                          !!>> HC 25-11-2020
      
      elseif (tind==5)then                                                  !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
            
            mutacode=7                                                     !!>> HC 29-11-2021 Store that this is a gene duplication/delection
            call random_number(a)                                          !!>> HC 10-9-2021
            i=ceiling(a*(ng))                                              !!>> HC 10-9-2021
            if (i==0) i=1                                                  !!>> HC 25-11-2020
            mutageni=i                                                     !!>> HC 16-9-2021 save the chosen gene for the output

            if(good(i) == 1)then; runornot = 1; end if        !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

            if (ng>2)then                                                  !!>> HC 25-11-2020
               call random_number(a)                                       !!>> HC 25-11-2020
               if (a<0.50d0)then                                           !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
                  mutagenj=-1                                              !!>> HC 29-11-2021 store that this is a delection
                  call deletion(i)                                         !!>> HC 25-11-2020
               else                                                        !!>> HC 25-11-2020
                  mutagenj=1                                               !!>> HC 29-11-2021 store that this is a duplication
                  call duplication(i)                                      !!>> HC 25-11-2020
               endif                                                       !!>> HC 25-11-2020
            else                                                           !!>> HC 25-11-2020
               mutagenj=1                                                  !!>> HC 29-11-2021 store that this is a duplication
               call duplication(i)                                         !!>> HC 25-11-2020
            endif                                                          !!>> HC 25-11-2020
      endif                                                                 !!>> HC 25-11-2020
   endif                                                                   !!>> HC 25-11-2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 
!!>> HC 16-9-2021 !! ATENTION !! We are not using pre/post forms right now, this should be changed if you want them!!      !!>> HC 16-9-2021
!  if (allocated(npag)) deallocate(npag)                             !!>> HC 25-11-2020  This last chunk of code was inherited from muta
!  allocate(npag(nga))                                               !!>> HC 25-11-2020  it realocates matrixes 
!  npag=0                                                            !!>> HC 25-11-2020

!  if (allocated(whonpag)) deallocate(whonpag)                       !!>> HC 25-11-2020
!  allocate(whonpag(nga,ng))                                         !!>> HC 25-11-2020
!  whonpag=0                                                         !!>> HC 25-11-2020

!  call update_npag                                                  !!>> HC 25-11-2020

!  do i=1,ng                                                         !!>> HC 25-11-2020
!     do j=1,gen(i)%npre                                             !!>> HC 25-11-2020
!        if (gen(i)%pre(j)==0) then                                  !!>> HC 25-11-2020
!           print *,i,j,"aqui hi ha un zero a pre"                   !!>> HC 25-11-2020
!           stop                                                     !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020

!     do j=1,gen(i)%npost                                            !!>> HC 25-11-2020
!        if (gen(i)%post(j)==0) then                                 !!>> HC 25-11-2020
!        print *,i,j,"aqui hi ha un zero a post"                     !!>> HC 25-11-2020
!        stop                                                        !!>> HC 25-11-2020
!        end if                                                      !!>> HC 25-11-2020
!     end do                                                         !!>> HC 25-11-2020
!  end do                                                            !!>> HC 25-11-2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 

   inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
   if (ng<1) then                                                    !!>> HC 25-11-2020
      print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   else                                                              !!>> HC 25-11-2020
      do j=1,ng                                                      !!>> HC 25-11-2020
         if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
            inviable=0                                               !!>> HC 25-11-2020
            goto 37                                                  !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   end if                                                            !!>> HC 25-11-2020
37 continue                                                         !!>> HC 25-11-2020
         
  end subroutine suremuta_no_kadh_functional_net

!************************************************************************************************************!

subroutine suremuta_no_kadh_functional_net_rec(limit,inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue,rangfile,runornot,&
pcpval) !!>> AL 29-11-24
   !!>> AL 9-4-24: this function determines whether the mutation alters the functional part of the gen network of an individual or not.
   !!>> AL 9-4-24: if not then development of mutated individual is not runned and fitness is inherited from parent.
   !!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file 
   !!>> HC 14-4-2023 This version is the same but ignoring specific adhesion molecules
   !!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
   !!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
   !!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
   !!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
   !!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
   implicit none

   integer                     :: limit,inviable,mind,tind,isitreal,realgenes
   integer                     :: props,nume,numt,numr,numk !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
   integer                     :: luckyprop, luckylink, luckyg  !!>> HC 24-11-2020 selected gene, interaction or property
   integer                     :: ich,jch,lch,ord  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
   integer                     :: full,empty !!>> HC 24-11-2020 full an empty spaces for T mutations
   real*8, dimension(1:5)      :: probs !!>> HC 24-11-2020 relative probabilities of different T mutations
   !real*8, dimension(1:7)      :: probis !!>> HC 24-11-2020 relative probabilities of different IS mutations
   real*8, dimension(5)        :: probis  !AL 2-12-24 this is 5 since there are no reactions 
   real*8, dimension(1:nga)    :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 

   integer, dimension(1:nga)   :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
   real*8                      :: cum !!>> HC 24-11-2020 funny position to store cumulative Fi
   integer                     :: signo  !!>> HC 24-11-2020 sign of the mutation  
   real*8                      :: inch, newval, anch !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
   integer                     :: mutacode,mutageni,mutagenj  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
   real*8                      :: prevalue, newvalue
   character*400               :: rangfile !!>> 6-10-2021 file where the ranges are stored
   real*8, dimension (1:5)     :: min_glim, max_glim                                  !!>>HC 6-10-2021
   real*8                      :: t_mag,pcpval   !!>> HC 31-5-2023 Magnitude of the initial value in T mutations (should be a free parameter) !!>> AL 28-12-2024 pcpval = magnitude of planar cell contraction 
   integer                     :: errror     !!>>AL 17-1-25
   real(8)                     :: tolerance  !!>>AL 21-1-25   
   integer,allocatable,dimension(:)  :: plinki,plinkj !!>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
   integer,allocatable,dimension(:)  :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM+
   
   !functional subnetwork                       !!>>AL 9-4-24
   integer :: runornot                          !!>>AL 10-4-24: if == 1 you run development, if 0 you do not          
   integer, allocatable, dimension(:) :: good   !!>>AL 10-4-24: this vector has elements = 1 when genes are functional
   if(allocated(good)) deallocate(good)         !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
   allocate(good(1:ng))                 
   pcpval = 0.0d0                               !!>>AL 28-12-24:
   errror = 0                                   !!>>AL 17-1-25:
   
   !limits and ranges        
   limit=1                                               !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
   rembeh=0; max_elim=0.0d0; min_elim=0.0d0              !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
   max_glim=0.0d0; min_glim=0.0d0                        !!>> HC 6-10-2020 Same for other gene properties 
   call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
   d_max=max_glim(1)    ; d_min=min_glim(1)              !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020
   m_max=max_glim(2)    ; m_min=min_glim(2)              !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
   mich_max=max_glim(3) ; mich_min=min_glim(3)           !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
   w_max=max_glim(4)    ; w_min=min_glim(4)              !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
   b_max=max_glim(5)    ; b_min=min_glim(5)              !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

   call functional_net(good)                             !!>>AL 29-11-24
   empty=0; full=0                                          !!>> HC 24-11-2020
   inch=0.0d0; anch=0.0d0                                   !!>> HC 28-11-2020
   mutacode=0                                               !!>> HC 16-9-2021  
   mutageni=0                                               !!>> HC 29-11-2021
   mutagenj=0                                               !!>> HC 29-11-2021
   prevalue=0.0d0                                           !!>> HC 29-11-2021
   newvalue=0.0d0                                           !!>> HC 29-11-2021
   t_mag=0.20d0
   realgenes = 0                                            !! AL 2-12-24
   
   if (mind==1)then                                        !!>> HC 25-11-2020 !******************* IS MUTATIONS***************************!
      nume=0; numt=0; numk=0; numr=0                       !!>> HC 25-11-2020 
      do ich=1,ng                                          !!>> HC 25-11-2020 
         do jch=1,ng                                       !!>> HC 25-11-2020 
            if(gen(ich)%t(jch).ne.0.0d0) numt=numt+1       !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
         enddo                                             !!>> HC 25-11-2020 
         
         do jch=2,nga                                      !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
            if(gen(ich)%e(jch).ne.0.0d0) nume=nume+1       !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
         enddo                                             !!>> HC 25-11-2020 

         if(int(gen(ich)%kindof) .ne. 9) realgenes = realgenes + 1 !AL 21-1-25 finding number of real genes
      end do                                                
      
      probis=0.0d0;props=0                                 !!>> HC 25-11-2020 props is the total number of parameters that can change
      
      props=3*realgenes+numt+nume                           !AL 2-12-24 

      probis(1)=real(realgenes)/real(props)                  !!>> HC 25-11-2020 probability mu     
      probis(2)=real(realgenes)/real(props)                  !!>> HC 25-11-2020 probability diff
      probis(3)=real(realgenes)/real(props)                  !!>> HC 25-11-2020 probability mich
      probis(4)=real(numt)/real(props)                     !!>> HC 25-11-2020 probability t matrix
      probis(5)=real(nume)/real(props)                     !!>> HC 25-11-2020 probability e matrix
      !probis(6)=real(numr)/real(props)                     !!>> HC 25-11-2020 probability nww matrix == 0  !AL 2-12-24 this could be ereased 
      !probis(7)=real(numk)/real(props)                     !!>> HC 25-11-2020 probability kadh matrix == 0 !AL 2-12-24 this could be ereased 
      
      cum=0.0d0                                            !!>> HC 25-11-2020
      call random_number(a)                                !!>> HC 25-11-2020
      do ich=1,size(probis)                                !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
         cum=cum+probis(ich)                               !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
         if (a<cum)then                                    !!>> HC 25-11-2020
            luckyprop=ich                                  !!>> HC 25-11-2020
            exit                                           !!>> HC 25-11-2020
         endif                                             !!>> HC 25-11-2020
      enddo                                                !!>> HC 25-11-2020
      mutacode=luckyprop                                   !!>> HC 16-9-2021  Store mutation code for output   
      
      if (luckyprop==1)then                                !!>> HC 25-11-2020 !***********IS MUTATION IN DEGRADATION RATE (mu)**************!

         isitreal=0  
         do while(isitreal .eq. 0)                         !!>> AL 29-11-24 randomly choose a gene to mutate and check if is real or a ghost                  
            call random_number(a)                                         
            i=ceiling(a*ng)     
            if (i==0) i=1             
            if(int(gen(i)%kindof) .ne. 9) then                !!>> AL 29-11-24 if kindof == -1 then is a ghost gene 
               isitreal=1                             
               exit
            end if
         end do

         mutageni=i                                                        !!>> HC 29-11-2021 We store the gene for the output !AL 3-12-24 we could just directly use mutageni instead of i

         if(good(mutageni) == 1) runornot = 1                              !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%mu                                                !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(i)%mu==0.0d0) then                                        !!>> HC 25-11-2020
            inch=m_max*a 
         else 
            inch=(1-2*a)*gen(i)%mu*mag_ismurate                            !!>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if ((gen(i)%mu+inch) < 0.0d0)then                              !!>> ALI 2-12-24 we always want mutations to change the selected parameters
            if(gen(i)%mu == m_min)then
               inch = -inch                                         
               gen(i)%mu = gen(i)%mu+inch  
            else
               gen(i)%mu = m_min 
            end if
         else if ((gen(i)%mu+inch) > m_max) then 
            if(gen(i)%mu == m_max)then 
               inch = -inch
               gen(i)%mu=gen(i)%mu+inch
            else 
               gen(i)%mu = m_max
            end if
         else 
            gen(i)%mu=gen(i)%mu+inch
         end if

         newvalue=gen(i)%mu                                                !!>> HC 29-11-2021  we store the new value
         
      elseif (luckyprop==2)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!

         isitreal=0
         do while(isitreal .eq. 0)                        !!>> AL 29-11-24 check if the selected gene is real or a ghost                  
            call random_number(a)                                             
            i=ceiling(a*ng)                                                  
            if (i==0) i=1             
            if(int(gen(i)%kindof) .ne. 9) then !!>> AL 29-11-24 if kindof == 9 then is a ghost gene
               isitreal=1
               exit
            end if
         end do
         
         mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
         
         if(good(mutageni) == 1) runornot = 1                              !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%diffu                                             !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(i)%diffu==0.0d0) then                                        !!>> HC 25-11-2020
            inch=d_max*a 
         else 
            inch=(1-2*a)*gen(i)%diffu*mag_ismurate                            !!>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if ((gen(i)%diffu+inch) < 0.0d0)then                              !!>> AL 2-12-24
            if(gen(i)%diffu == d_min)then
               inch = -inch                                         
               gen(i)%diffu = gen(i)%diffu+inch  
            else
               gen(i)%diffu = d_min 
            end if
         else if ((gen(i)%diffu+inch) > d_max) then 
            if(gen(i)%diffu == d_max)then 
               inch = -inch
               gen(i)%diffu=gen(i)%diffu+inch
            else 
               gen(i)%diffu = d_max
            end if
         else 
            gen(i)%diffu = gen(i)%diffu + inch
         end if
         
         newvalue=gen(i)%diffu                                             !!>> HC 29-11-2021  we store the new value      


      elseif (luckyprop==3)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!

         isitreal=0
         do while(isitreal .eq. 0)                        !!>> AL 29-11-24 check if the selected gene is real or a ghost                  
            call random_number(a)                                             
            i=ceiling(a*ng)                                                  
            if (i==0) i=1             
            if(int(gen(i)%kindof) .ne. 9) then !!>> AL 29-11-24 if kindof == 9 then is a ghost gene
               isitreal=1
               exit
            end if
         end do

         mutageni=i                                                        !!>> HC 29-11-2021  We store the gene for the output
         
         if(good(mutageni) == 1) runornot = 1                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         
         prevalue=gen(i)%mich                                              !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(i)%mich==0.0d0) then                                     !!>> AL 3-12-24 
            inch=mich_max*a 
         else 
            inch=(1-2*a)*gen(i)%mich*mag_ismurate                         
         end if
         
         if ((gen(i)%mich+inch) < 0.0d0)then                              !!>> AL 2-12-24
            if(gen(i)%mich == mich_min)then
               inch = -inch                                         
               gen(i)%mich = gen(i)%mich + inch  
            else
               gen(i)%mich = mich_min 
            end if
         else if ((gen(i)%mich + inch) > mich_max) then 
            if(gen(i)%mich == mich_max)then 
               inch = -inch
               gen(i)%mich = gen(i)%mich + inch
            else 
               gen(i)%mich = mich_max
            end if
         else 
            gen(i)%mich = gen(i)%mich + inch
         end if

         newvalue=gen(i)%mich                                              !!>> HC 29-11-2021  we store the new value  
      

      elseif (luckyprop==4)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!

         if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020  
         allocate(plinki(1:numt))                                          !!>> HC 25-11-2020
         if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
         allocate(plinkj(1:numt))                                          !!>> HC 25-11-2020
         ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020
        
         do ich=1,ng                                                       !!>> HC 25-11-2020

            if(int(gen(ich)%kindof) .eq. 9) cycle              !!>> AL 29-11-24 ghost genes are discriminated
            
            do jch=1,ng                                                    !!>> HC 25-11-2020
               if(int(gen(jch)%kindof) .eq. 9) cycle              !!>> AL 29-11-24 ghost genes are discriminated
               if (gen(ich)%t(jch).ne.0.0d0)then                           !!>> HC 25-11-2020
                  ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (from 1 to the number of iterations=numt)
                  plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that are transcription factors
                  plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of genes j whose transcription is regulated
               endif                                                       !!>> HC 25-11-2020
            enddo                                                          !!>> HC 25-11-2020
         enddo                                                             !!>> HC 25-11-2020
         
         call random_number(a)                                             !!>> HC 25-11-2020
         luckylink=ceiling(a*numt)                                         !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
         if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
        
         if(good(plinki(luckylink)) == 1) runornot = 1                     !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
                                                                           !!>> AL 10-4-24: if plinki(luckylink) then plinkj is functional (by definition)
         i=plinki(luckylink)                                               !!>> HC 25-11-2020 gene i is the randomly selected transcription factor
         j=plinkj(luckylink)                                               !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
         mutageni=i                                                        !!>> HC 29-11-2021    
         mutagenj=j                                                        !!>> HC 29-11-2021
         
         prevalue=gen(i)%t(j)                                              !!>> HC 29-11-2021 store the previous value

         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                           !!>> AL 2-12-24 we don't want neutral mutations
        
         inch=(1-2*a)*gen(i)%t(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
         if ((gen(i)%t(j) + inch) < w_min)then                              !!>> AL 2-12-24
            if(gen(i)%t(j) == w_min)then
               inch = -inch                                         
               gen(i)%t(j) = gen(i)%t(j) + inch  
            else
               gen(i)%t(j) = w_min
            end if
         else if ((gen(i)%t(j) + inch) > w_max) then 
            if(gen(i)%t(j) == w_max)then
               inch = -inch                                         
               gen(i)%t(j) = gen(i)%t(j) + inch  
            else
               gen(i)%t(j) = w_max
            end if
         else 
            gen(i)%t(j) = gen(i)%t(j) + inch
         end if
         
         newvalue=gen(i)%t(j)                                              !!>> HC 29-11-2021

      elseif (luckyprop==5)then                                            !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
      
         if(allocated(plinki)) deallocate(plinki)                          !!>> HC 25-11-2020
         allocate(plinki(1:nume))                                          !!>> HC 25-11-2020
         if(allocated(plinkj)) deallocate(plinkj)                          !!>> HC 25-11-2020
         allocate(plinkj(1:nume))                                          !!>> HC 25-11-2020
      
         ord=0;plinki=0;plinkj=0                                           !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
         
         do ich=1,ng                                                       !!>> HC 25-11-2020
            
            if(int(gen(ich)%kindof) .eq. 9) cycle                                  !!>> AL 29-11-24 ghost genes are discriminated
            do jch=2,nga                                                   !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
               if (gen(ich)%e(jch) .ne. 0.0d0)then                         !!>> HC 25-11-2020
                  ord=ord+1                                                !!>> HC 25-11-2020 order of vectors (1 to nume)
                  plinki(ord)=ich                                          !!>> HC 25-11-2020 indexes of genes i that regulate cell behaviors
                  plinkj(ord)=jch                                          !!>> HC 25-11-2020 indexes of regulated cell behaviors
               endif                                                       !!>> HC 25-11-2020
            enddo                                                          !!>> HC 25-11-2020
         enddo                                                             !!>> HC 25-11-2020
      
         call random_number(a)                                             !!>> HC 25-11-2020
         luckylink=ceiling(a*nume)                                         !!>> HC 25-11-2020 Randomly selecting 
         if(luckylink==0)luckylink=1                                       !!>> HC 25-11-2020
         runornot = 1 
         !if(good(plinki(luckylink)) == 1) runornot = 1                     !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
 
         i=plinki(luckylink)                                               !!>> HC 25-11-2020 Randomly selected gene i that regulates
         j=plinkj(luckylink)                                               !!>> HC 25-11-2020 Randomly selected behavior/property j
         mutageni=i                                                        !!>> HC 29-11-2021 store for output    
         mutagenj=j                                                        !!>> HC 29-11-2021 store for output
         
         prevalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output
         
         call random_number(a)                                             !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.02                                           !!>> AL 2-12-24 we don't want neutral mutations

         inch=(1-2*a)*gen(i)%e(j)*mag_ismurate                             !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation

         if ((gen(i)%e(j) + inch) < min_elim(j))then                       !!>> AL 2-12-24
            if(gen(i)%e(j) == min_elim(j))then
               inch = -inch                                         
               gen(i)%e(j) = gen(i)%e(j) + inch  
            else
               gen(i)%e(j) = min_elim(j)
            end if
         else if ((gen(i)%e(j) + inch) > max_elim(j)) then 
            if(gen(i)%e(j) == max_elim(j))then 
               inch = -inch
               gen(i)%e(j) = gen(i)%e(j) + inch  
            else 
               gen(i)%e(j) = max_elim(j)
            end if
         else 
            gen(i)%e(j) = gen(i)%e(j) + inch
         end if
         
         newvalue=gen(i)%e(j)                                              !!>> HC 29-11-2021 store for output  

      endif                                                                !!>> HC 25-11-2020

   else                      !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!! mind == 2
      probs(1)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
      probs(2)=33d-2         !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
      probs(3)=0d0           !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now) !AL 2-12-24 this could be ereased
      probs(4)=32d-2         !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
      probs(5)=2d-2          !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals
      
      call random_number(a)  !!>> HC 25-11-2020
      cum=0.0d0              !!>> HC 25-11-2020
      do ich=1,size(probs)   !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
         cum=cum+probs(ich)  !!>> HC 25-11-2020 taking into account probabilites in probs
         if (a<cum)then      !!>> HC 25-11-2020
            tind=ich         !!>> HC 25-11-2020
            exit             !!>> HC 25-11-2020
         endif               !!>> HC 25-11-2020
      enddo                  !!>> HC 25-11-2020
      
      !we force duplication for testing 

      if (tind .le. 2)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         mutacode=4                                                         !!>> HC 29-11-2021 Store the matrix that will mutate
         realgenes = 0
         
         do ich=1,ng                                                        !!>> HC 25-11-2020
            
            if(int(gen(ich)%kindof) .ne. 9) realgenes = realgenes + 1               !AL 21-1-25 finding number of real genes
            
            if(int(gen(ich)%kindof) .eq. 9) cycle                !!>> AL 29-11-24 ghost genes are discriminated
            do jch=1,ng                                                     !!>> HC 25-11-2020
               if(int(gen(jch)%kindof) .eq. 9) cycle                             !!>> AL 29-11-24 ghost genes are discriminated
               if(gen(ich)%t(jch)==0) empty=empty+1                         !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
            enddo                                                           !!>> HC 25-11-2020 and the number of already existing links in t matrix
         enddo                                                              !!>> HC 25-11-2020
         
         full = realgenes*realgenes - empty                                !AL 2-12-24 important

         if (tind==1.and.full>0)then                                      !!>> HC 25-11-2020 !!!!!! DELETION OF A LINK IN THE T MATRIX!!!!!! !!>> AL 29-11-24 the full>0 shouldnt be necessary. Check: if its needed, fix, if not, erease
            
            if (allocated(fullspoti)) deallocate(fullspoti)                 !!>> HC 25-11-2020
            allocate(fullspoti(1:full))                                     !!>> HC 25-11-2020
            if (allocated(fullspotj)) deallocate(fullspotj)                 !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                                     !!>> HC 25-11-2020
            
            ord=0;fullspoti=0;fullspotj=0                                   !!>> HC 25-11-2020 These vectors store the already existing links
            
            do ich=1,ng                                                     !!>> HC 25-11-2020
               if(int(gen(ich)%kindof) .eq. 9) cycle              !!>> AL 29-11-24 ghost genes are discriminated
               do jch=1,ng                                                  !!>> HC 25-11-2020
                  if(int(gen(jch)%kindof) .eq. 9) cycle                        !!>> AL 29-11-24 ghost genes are discriminated
                  if (gen(ich)%t(jch).ne.0.0d0)then                         !!>> HC 25-11-2020
                     ord=ord+1                                              !!>> HC 25-11-2020 order of vectors (1 to full)
                     fullspoti(ord)=ich                                     !!>> HC 25-11-2020 this gene i regulates the expression of
                     fullspotj(ord)=jch                                     !!>> HC 25-11-2020 gene j
                  endif                                                     !!>> HC 25-11-2020
               enddo                                                        !!>> HC 25-11-2020
            enddo                                                           !!>> HC 25-11-2020
            
            call random_number(a)                                           !!>> HC 25-11-2020
            luckyg=ceiling(a*full)                                          !!>> HC 25-11-2020 Randomly pickin an interaction to delete
            if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020

            if(good(fullspoti(luckyg)) == 1) runornot = 1                   !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

            i=fullspoti(luckyg)                                             !!>> HC 25-11-2020 Removing interaction between gene i and j
            j=fullspotj(luckyg)                                             !!>> HC 25-11-2020
            mutageni=i                                                      !!>> HC 29-11-2021 Store for output
            mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
           
            prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output
            gen(i)%t(j)=0.0d0                                               !!>> HC 25-11-2020

            newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 

         else                                                               !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
            if (allocated(freespoti)) deallocate(freespoti)                 !!>> HC 25-11-2020
            allocate(freespoti(1:empty))                                    !!>> HC 25-11-2020
            if (allocated(freespotj)) deallocate(freespotj)                 !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                                    !!>> HC 25-11-2020
            
            ord=0;freespoti=0;freespotj=0                                   !!>> HC 25-11-2020 These vectors store the free spots
            
            do ich=1,ng                                                      !!>> HC 25-11-2020
               if(int(gen(ich)%kindof) .eq. 9) cycle              !!>> AL 29-11-24 ghost genes are discriminated
               do jch=1,ng                                                   !!>> HC 25-11-2020
                  if(int(gen(jch)%kindof) .eq. 9) cycle                         !!>> AL 29-11-24 ghost genes are discriminated
                  if (gen(ich)%t(jch)==0)then                                !!>> HC 25-11-2020
                     ord=ord+1                                               !!>> HC 25-11-2020 order of the vectors (1 to empty)
                     freespoti(ord)=ich                                      !!>> HC 25-11-2020 gene i has NOT and interaction !!AL 2-12-24 wtf is NOT?
                     freespotj(ord)=jch                                      !!>> HC 25-11-2020 with gene j
                  endif                                                      !!>> HC 25-11-2020
               enddo                                                         !!>> HC 25-11-2020
            enddo                                                            !!>> HC 25-11-2020
            
            call random_number(a)                                           !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                         !!>> HC 25-11-2020 Randomly picking a spot
            if(luckyg==0)luckyg=1                                           !!>> HC 25-11-2020

            call random_number(a)                                           !!>> HC 25-11-2020
            if(a < 1.0d-8) a = 0.01                                         !!>> AL 2-12-24 we don't want neutral mutations
            call random_number(b)                                           !!>> HC 25-11-2020
            if(b<0.50d0)then;signo=1;else;signo=-1;endif                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
            
            i=freespoti(luckyg)                                             !!>> HC 25-11-2020
            j=freespotj(luckyg)                                             !!>> HC 25-11-2020 Making new interaction
            mutageni=i                                                      !!>> HC 29-11-2021 Store for output
            mutagenj=j                                                      !!>> HC 29-11-2021 Store for output
            
            prevalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output

            gen(i)%t(j)=a*signo*w_max*t_mag                                !!>> HC 31-5-2023 

            newvalue=gen(i)%t(j)                                            !!>> HC 29-11-2021 Store for output 

            call functional_net(good)                                       !AL 29-11-24 important
            if(good(freespoti(luckyg)) == 1) runornot = 1                   !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

         endif                                                              !!>> HC 25-11-2020
     
      elseif (tind==3)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
      !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW        !!>> HC 25-11-2020 
      
      elseif (tind==4)then                                                  !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e matrix)!!!!!!
         empty=0                                                        !!>> HC 25-11-2020
         
         do ich=1,ng                                                    !!>> HC 25-11-2020  
            if(int(gen(ich)%kindof) .eq. 9) cycle                               !!>> AL 29-11-24 ghost genes are discriminated            
            do jch=2,nga                                                !!>> AL 2-12-24
               if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch !!>> AL 2-12-24 not sure what HC wanted to communicate but if rembeh= 1 that behaviour is off in ranges.dat so we dont consider it
               if (gen(ich)%e(jch)==0.0d0) then                         !!>> HC 27-11-2020
                  empty=empty+1                                         !!>> HC 27-11-2020 Number of empty spots in e Matrix
               else                                                     !!>> HC 27-11-2020
                  full=full+1                                           !!>> HC 27-11-2020 Number of existing interactions in e matrix
               endif                                                    !!>> HC 27-11-2020
            enddo                                                       !!>> HC 25-11-2020
         enddo                                                          !!>> HC 25-11-2020

         call random_number(a)                                          !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
         if ((a<0.50d0 .and. full>1) .or. empty==0)then                 !!>> HC 10-5-2022 !!!REMOVE INTERACTION !!>> AL 28-1-25: note that full>1, otherwhise we can remove the only cell behaviour and end up without a developmental mechanism 

            runornot = 1                                                !!>> AL 10-4-24: modifying a cell behaviour always means altering functional network 
                                                                        !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(fullspoti)) deallocate(fullspoti)             !!>> HC 25-11-2020
            allocate(fullspoti(1:full))                                 !!>> HC 25-11-2020
            if (allocated(fullspotj)) deallocate(fullspotj)             !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                                 !!>> HC 25-11-2020
            ord=0;fullspoti=0;fullspotj=0                               !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
            
            do ich=1,ng                                                 !!>> HC 25-11-2020
               if(int(gen(ich)%kindof) .eq. 9) cycle                            !!>> AL 29-11-24 ghost genes are discriminated
               do jch=2,nga                                             !!>> HC 25-11-2020 we start in 2 because 1 is being an adhesion mol
                  if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
                  if (gen(ich)%e(jch).ne. 0.0d0)then                    !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                     ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                     fullspoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                     fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
                  endif                                                 !!>> HC 25-11-2020
               enddo                                                    !!>> HC 25-11-2020
            enddo                                                       !!>> HC 25-11-2020
            
            call random_number(a)                                       !!>> HC 25-11-2020
            luckyg=ceiling(a*real(full))                                !!>> HC 25-11-2020 Randomly picking an interaction to kill
            if(luckyg==0)luckyg=1                                       !!>> AL 3-12-2024   

            i=fullspoti(luckyg)                                         !!>> HC 25-11-2020
            j=fullspotj(luckyg)                                         !!>> HC 25-11-2020
            
            !if(i==24 .and. j==37)then                                   !!>> AL 28-1-25 important: we always have cell division consitutively being regulated by gene 24
            !   limit = 0
            !end if 
                                                                        !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
            newvalue=0.0d0                                              !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
            gen(i)%e(j)=0.0d0                                           !!>> HC 25-11-2020 we kill it easily

         else                                                           !!>> HC 25-11-2020 !ADD AN INTERACTION   
            runornot = 1                                                !!>> AL 10-4-24: modifying a cell behaviour always means altering functional newtwork 
                                                                        !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(freespoti)) deallocate(freespoti)             !!>> HC 25-11-2020
            allocate(freespoti(1:empty))                                !!>> HC 25-11-2020
            if (allocated(freespotj)) deallocate(freespotj)             !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                                !!>> HC 25-11-2020
            
            ord=0;freespoti=0;freespotj=0;                              !!>> HC 25-11-2020 
            
            do ich=1,ng                                                 !!>> HC 25-11-2020
               if(int(gen(ich)%kindof) .eq. 9) cycle                            !!>> AL 29-11-24 ghost genes are discriminated 
               do jch=2,nga                                             !!>> HC 25-11-2020
                  if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
                  if (gen(ich)%e(jch)==0.0d0)then                       !!>> HC 25-11-2020
                     ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                     freespoti(ord)=ich                                 !!>> HC 25-11-2020 gene i regulates
                     freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
                  endif                                                 !!>> HC 25-11-2020
               enddo                                                    !!>> HC 25-11-2020
            enddo                                                       !!>> HC 25-11-2020
            
            call random_number(a)                                       !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                     !!>> HC 25-11-2020 Randomly picking an interaction
            if(luckyg==0)luckyg=1                                       !!>> HC 25-11-2020
            
            i=freespoti(luckyg)                                         !!>> HC 25-11-2020
            j=freespotj(luckyg)                                         !!>> HC 25-11-2020 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                                  !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=i                                                  !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                                  !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(i)%e(j)                                        !!>> HC 29-11-2021 save for output the prev value
            
            call random_number(c)                                       !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
            if(c < 1.0d-8) c = 0.01                                      !!>> AL 2-12-24 we don't want neutral mutations

            anch=max_elim(j)-min_elim(j)                                !!>> HC 28-11-2020 we introduce an interaction in the range of the limits            

            gen(i)%e(j)=anch*c+min_elim(j)                              !!>> AL 3-12-24
            gen(i)%e(j)=gen(i)%e(j)*t_mag                               !!>> HC 31-5-2023  Correct the magnitude of the new interaction 

            if (j==nparam_per_node+8)then                               !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation !>> AL 2-12-24 nparam_per_node = 35 PCP in the papers called PCO for whatever reason 
               call random_number(c)                                    !!>> AL 2-12-24 
               if(c < 1.0d-8) c = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations
               
               gen(24)%e(nparam_per_node+17) = 0.20d0*c                 !!>> AL 2-12-24 we assign a random magnitud with respecto to max val (0.20d0)!important: we should find which gene regulates ceel division, not assume that its the 24th
               
               !pcpval = 0.20d0*c 
               !gen(ng)%e(nparam_per_node+17) =0.20d0                   !!>> HC 12-11-2020 !!!ATENTION!! gene ng is homogeneously expressed in the blastula ics !>> AL 2-12-24 important: this shouldnt start at the max possible value..
            endif                                                       !!>> HC 7-12-2020 !!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
            pcpval = gen(24)%e(nparam_per_node+17)                      !!>> AL 10-1-25 
            newvalue=gen(i)%e(j)                                        !!>> HC 29-11-2020  Store the new value 
         endif                                                          !!>> HC 25-11-2020
      
      elseif (tind==5)then                                              !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
   
         mutacode=7                                                   !!>> HC 29-11-2021 Store that this is a gene duplication/delection
         call random_number(a)                                        !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
         if (a<0.50d0) then                                           !!>> AL 17-1-25 Deletion   
            if(ng<2)then 
               go to 12345                                            !!>> AL 17-1-25 If the deletion will erease the genome we do a duplication instead 
            end if 
            lch = 0
            isitreal = 0  
            do while(isitreal .eq. 0)                                !!>> AL 29-11-24 check if the selected gene is real or a ghost  !!>> AL 27-1-25 just identify ghost genes and chose from them for for fuck's sake
               call random_number(a)                                             
               i=ceiling(a*ng)                                                  
               if (i==0) i=1             
               if(int(gen(i)%kindof) .ne. 9)then 
                  isitreal=1
                  exit
               end if
               lch = lch+1
               if(lch .eq. 5*ng) then                                    !!>> AL 17-1-25 Deletion 
                  errror = 1
                  exit
               end if
            end do
            
            if(errror .eq. 1)then                                      !!>> AL 17-1-25 Deletion 
               mutageni=0 
               mutagenj=-1
            else 
               mutageni=i 
               mutagenj=-1                                              !!>> HC 29-11-2021 store that this is a delection
               if(good(i) == 1) runornot = 1                            !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
               call deletion_rec(i)                                     !!>> AL 29-11-24
            end if
         else                                                           !!>> HC 25-11-2020 !!>> AL 17-1-25 Duplication
      12345 continue
            lch = 0
            do jch = 1,ng                                               !!>> AL 17-1-25 we have to check if there are ghost genes available
               if(int(gen(jch)%kindof) .eq. 9) then
                  exit
               end if
               lch = lch + 1
            end do 

            if(lch .eq. ng)then                                          !!>> AL 17-1-25 there are no ghost genes left so no gen duplications are possible      
               errror = 2
               go to 757
            end if

            lch      = 0
            isitreal = 0
            do while(isitreal .eq. 0)                                   !!>> AL 29-11-24 check if the selected gene is real or a ghost                  
               call random_number(a)                                             
               i=ceiling(a*ng)                                                  
               if (i==0) i=1             
               if(int(gen(i)%kindof) .ne. 9)then
                  isitreal=1
                  exit
               end if

               lch = lch + 1
               if(lch .eq. 5*ng) then                                     
                  errror = 2
                  exit
               end if
            end do

            if(errror .eq. 2)then                                      
               mutageni=0 
               mutagenj=1
            else 
               mutageni=i 
               mutagenj=1                                              !!>> HC 29-11-2021 store that this is a duplication
               if(good(i) == 1) runornot = 1                                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
               call duplication_rec(i)                                     !!>> AL 29-11-24
            end if
         endif                                                       !!>> HC 25-11-2020
      end if 
   endif                                                                   !!>> HC 25-11-2020
 
   if (errror == 1 .or. errror == 2) then
757   continue   
      inviable=1  
   end if  

   inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
   if (ng<1) then                                                    !!>> HC 25-11-2020 
      print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   else                                                              !!>> HC 25-11-2020
      do j=1,ng                                                      !!>> HC 25-11-2020 !!>> AL 2-12-24 important, this might not do what it should
         if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
            inviable=0                                               !!>> HC 25-11-2020
            goto 37                                                  !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   end if                                                            !!>> HC 25-11-2020

37 continue      
end subroutine suremuta_no_kadh_functional_net_rec

subroutine suremuta_no_kadh_functional_net_rec_2(whichgene,limit,inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue &
                                                ,rangfile,runornot,pcpval) !!>> AL 14-2-25
   !!>> AL 14-2-25 in this function you already chose wich gene to mutate
   
   !!>> AL 9-4-24: this function determines whether the mutation alters the functional part of the gen network of an individual or not.
   !!>> AL 9-4-24: if not then development of mutated individual is not runned and fitness is inherited from parent.
   !!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file 
   !!>> HC 14-4-2023 This version is the same but ignoring specific adhesion molecules
   !!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
   !!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
   !!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
   !!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
   !!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
   implicit none

   integer                     :: limit,inviable,mind,tind,isitreal,realgenes
   integer                     :: props,nume,numt,numr,numk          !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
   integer                     :: luckyprop, luckylink, luckyg       !!>> HC 24-11-2020 selected gene, interaction or property
   integer                     :: ich,jch,lch,ord                    !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
   integer                     :: full,empty                         !!>> HC 24-11-2020 full an empty spaces for T mutations
   real*8, dimension(1:5)      :: probs                              !!>> HC 24-11-2020 relative probabilities of different T mutations
   !real*8, dimension(1:7)      :: probis                            !!>> HC 24-11-2020 relative probabilities of different IS mutations
   real*8, dimension(5)        :: probis                             !AL 2-12-24 this is 5 since there are no reactions 
   real*8, dimension(1:nga)    :: max_elim, min_elim                 !!>> HC 27-11-2020 These vectors store the limits of e matrix 

   integer, dimension(1:nga)   :: rembeh                             !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
   real*8                      :: cum                                !!>> HC 24-11-2020 funny position to store cumulative Fi
   integer                     :: signo                              !!>> HC 24-11-2020 sign of the mutation  
   real*8                      :: inch, newval, anch                 !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
   integer                     :: mutacode,mutageni,mutagenj,whichgene  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
   real*8                      :: prevalue, newvalue
   character*400               :: rangfile                           !!>> 6-10-2021 file where the ranges are stored
   real*8, dimension (1:5)     :: min_glim, max_glim                 !!>>HC 6-10-2021
   real*8                      :: t_mag,pcpval                       !!>> HC 31-5-2023 Magnitude of the initial value in T mutations (should be a free parameter) !!>> AL 28-12-2024 pcpval = magnitude of planar cell contraction 
   integer                     :: errror                             !!>>AL 17-1-25
   real(8)                     :: tolerance                          !!>>AL 21-1-25   
   integer,allocatable,dimension(:)  :: plinki,plinkj                !>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
   integer,allocatable,dimension(:)  :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM+
   
   !functional subnetwork                                            !!>>AL 9-4-24
   integer :: runornot                                               !!>>AL 10-4-24: if == 1 you run development, if 0 you do not          
   integer, allocatable, dimension(:) :: good                        !!>>AL 10-4-24: this vector has elements = 1 when genes are functional

   if(allocated(good)) deallocate(good)                              !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
   allocate(good(1:ng))                                     

   pcpval = 0.0d0                                                    !!>>AL 28-12-24:
   errror = 0                                                        !!>>AL 17-1-25:
   
   !limits and ranges        
   limit=1                                                           !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
   rembeh=0          
   max_elim=0.0d0 ; min_elim=0.0d0                                   !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
   max_glim=0.0d0 ; min_glim=0.0d0                                   !!>> HC 6-10-2020 Same for other gene properties 
   
   call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
   
   d_max=max_glim(1)    ; d_min=min_glim(1)                          !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020
   m_max=max_glim(2)    ; m_min=min_glim(2)                          !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
   mich_max=max_glim(3) ; mich_min=min_glim(3)                       !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
   w_max=max_glim(4)    ; w_min=min_glim(4)                          !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
   b_max=max_glim(5)    ; b_min=min_glim(5)                          !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

   call functional_net(good)                                         !!>>AL 29-11-24

   empty    = 0
   full     = 0                                                      !!>> HC 24-11-2020
   mutacode = 0                                                      !!>> HC 16-9-2021  
   mutageni = 0                                                      !!>> HC 29-11-2021
   mutagenj = 0                                                      !!>> HC 29-11-2021
   realgenes= 0                                                      !! AL 2-12-24
   inch     = 0.0d0           
   anch     = 0.0d0                                                  !!>> HC 28-11-2020
   prevalue = 0.0d0                                                  !!>> HC 29-11-2021
   newvalue = 0.0d0                                                  !!>> HC 29-11-2021
   t_mag    = 0.20d0                                                 !!AL 2-12-24 important: this should be defined manually but in paraopc.par
   
   if (mind == 1) then                                               !!>> HC 25-11-2020 !******************* IS MUTATIONS***************************!
      
      nume=0; numt=0; numk=0; numr=0                       

      do jch=1,ng                                       
         if(gen(whichgene)%t(jch) .ne. 0.0d0) numt = numt + 1        !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
      enddo                                              
      
      do jch=2,nga                                                   !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
         if(gen(whichgene)%e(jch) .ne. 0.0d0) nume = nume + 1        !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
      enddo                                             

      do jch=1,ng                                       
         if(int(gen(jch)%kindof) .ne. 9) realgenes = realgenes + 1   !AL 21-1-25 finding number of real genes
      enddo        

      probis = 0.0d0                                                 !!>> HC 25-11-2020 props is the total number of parameters that can change

      props = 3 + numt + nume                                        !AL 2-12-24 

      probis(1) = 1/real(props)                                      !!>> AL: 14-2-25 probability mu     
      probis(2) = 1/real(props)                                      !!>> AL: 14-2-25 probability diff
      probis(3) = 1/real(props)                                      !!>> AL: 14-2-25 probability mich
      probis(4) = real(numt)/real(props)                             !!>> HC 25-11-2020 probability t matrix
      probis(5) = real(nume)/real(props)                             !!>> HC 25-11-2020 probability e matrix
      
      cum = 0.0d0                                                    !!>> HC 25-11-2020
      call random_number(a)                                          !!>> HC 25-11-2020

      do ich = 1, size(probis)                                       !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
         cum = cum + probis(ich)                                     !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
         if (a < cum) then                                           !!>> HC 25-11-2020
            luckyprop = ich                                          !!>> HC 25-11-2020
            exit                                                     !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      
      mutacode = luckyprop                                           !!>> HC 16-9-2021  Store mutation code for output   

      mutageni = whichgene

      if (good(whichgene) == 1) runornot = 1                         !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

      if (luckyprop == 1)then                                        !!>> HC 25-11-2020 !***********IS MUTATION IN DEGRADATION RATE (mu)**************!
         
         prevalue = gen(whichgene)%mu                                !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                       !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                     !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(whichgene)%mu == 0.0d0) then                        !!>> HC 25-11-2020
            !inch = m_max*a
            inch = (m_max*a)*mag_ismurate                            !!>> AL 28-2-25 we assume a simple GPM                               
         else 
            inch = (1 - 2*a)*gen(whichgene)%mu*mag_ismurate          !!>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if((gen(whichgene)%mu + inch) < m_min)then                  !!>> AL 2-12-24 we always want mutations to change the selected parameters
            if(gen(whichgene)%mu == m_min)then
               inch      = -inch                                         
               gen(whichgene)%mu = gen(whichgene)%mu+inch  
            else
               gen(whichgene)%mu = m_min 
            end if
         else if ((gen(whichgene)%mu+inch) > m_max) then 
            if(gen(whichgene)%mu == m_max)then 
               inch      = -inch
               gen(whichgene)%mu = gen(whichgene)%mu + inch
            else 
               gen(whichgene)%mu = m_max
            end if
         else 
            gen(whichgene)%mu = gen(whichgene)%mu + inch
         end if

         newvalue = gen(whichgene)%mu                                !!>> HC 29-11-2021  we store the new value
         
      elseif (luckyprop==2)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!

         prevalue = gen(whichgene)%diffu                             !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                       !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                     !!>> AL 2-12-24 we don't want neutral mutations

         if (gen(whichgene)%diffu == 0.0d0) then                     !!>> HC 25-11-2020
            inch = (d_max*a)*mag_ismurate                            !!>> AL 28-2-25 we assume a simple sequence-structure map  
         else 
            inch = (1-2*a)*gen(whichgene)%diffu*mag_ismurate         !>> HC 27-11-2020 inch saves the magnitude of the change
         end if

         if((gen(whichgene)%diffu + inch) < d_min)then               !!>> AL 2-12-24
            if(gen(whichgene)%diffu == d_min)then
               inch = -inch                                         
               gen(whichgene)%diffu = gen(whichgene)%diffu+inch  
            else
               gen(whichgene)%diffu = d_min 
            end if
         else if((gen(whichgene)%diffu + inch) > d_max)then 
            if(gen(whichgene)%diffu == d_max)then 
               inch         = -inch
               gen(whichgene)%diffu = gen(whichgene)%diffu + inch
            else 
               gen(whichgene)%diffu = d_max
            end if
         else 
            gen(whichgene)%diffu = gen(whichgene)%diffu + inch
         end if
         
         newvalue = gen(whichgene)%diffu                             !!>> HC 29-11-2021  we store the new value      

      elseif (luckyprop==3)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!

         prevalue = gen(whichgene)%mich                              !!>> HC 29-11-2021 We store the old value
         
         call random_number(a)                                       !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                     !!>> AL 2-12-24 we don't want neutral mutations
         if (gen(whichgene)%mich == 0.0d0) then                      !!>> AL 3-12-24 
            inch = (mich_max*a)*mag_ismurate                         !!>> AL 28-2-25 we assume a simple sequence-structure map   
         else 
            inch = (1-2*a)*gen(whichgene)%mich*mag_ismurate                         
         end if
         
         if ((gen(whichgene)%mich + inch) < mich_min)then            !!>> AL 2-12-24
            if(gen(whichgene)%mich == mich_min)then
               inch = -inch                                         
               gen(whichgene)%mich = gen(whichgene)%mich + inch  
            else
               gen(whichgene)%mich = mich_min 
            end if
         else if ( (gen(whichgene)%mich + inch) > mich_max) then 
            if(gen(whichgene)%mich == mich_max)then 
               inch = -inch
               gen(whichgene)%mich = gen(whichgene)%mich + inch
            else 
               gen(whichgene)%mich = mich_max
            end if
         else 
            gen(whichgene)%mich = gen(whichgene)%mich + inch
         end if

         newvalue = gen(whichgene)%mich                              !!>> HC 29-11-2021  we store the new value  
      
      elseif (luckyprop==4)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         if(allocated(plinkj)) deallocate(plinkj)                    !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
         allocate(plinkj(1:numt))                                    !!>> HC 25-11-2020
         
         ord   = 0
         plinkj= 0                                           
        
         do jch=1,ng                                                    
            if(int(gen(jch)%kindof) .eq. 9) cycle                    !!>> AL 29-11-24 ghost genes are discriminated
            if (gen(whichgene)%t(jch) .ne. 0.0d0) then                          
               ord = ord + 1                                                
               plinkj(ord) = jch                                         
            end if                                                       
         end do                                                      !!>> HC 25-11-2020
                 
         call random_number(a)                                       !!>> HC 25-11-2020
         
         luckylink = ceiling(a*numt)                                 !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
         if (luckylink == 0) luckylink = 1                           !!>> HC 25-11-2020
        
         j = plinkj(luckylink)                                       !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
         mutagenj = j                                                !!>> HC 29-11-2021
         
         prevalue = gen(whichgene)%t(j)                              !!>> HC 29-11-2021 store the previous value

         call random_number(a)                                       !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.01                                     !!>> AL 2-12-24 we don't want neutral mutations
        
         inch = (1 - 2*a)*gen(whichgene)%t(j)*mag_ismurate           !!>> HC 27-11-2020 inch saves the magnitude of the IS mutation
         if ((gen(whichgene)%t(j) + inch) < w_min)then               !!>> AL 2-12-24
            if(gen(i)%t(j) == w_min)then
               inch = -inch                                         
               gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch  
            else
               gen(whichgene)%t(j) = w_min
            end if

         else if ((gen(whichgene)%t(j) + inch) > w_max) then 
            if(gen(whichgene)%t(j) == w_max)then
               inch = -inch                                         
               gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch  
            else
               gen(whichgene)%t(j) = w_max
            end if
         else 
            gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch
         end if
         
         newvalue=gen(whichgene)%t(j)                                !!>> HC 29-11-2021

      elseif (luckyprop==5)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
      
         if(allocated(plinkj)) deallocate(plinkj)                    !!>> HC 25-11-2020
         allocate(plinkj(1:nume))                                    !!>> HC 25-11-2020
      
         ord=0
         plinkj=0                                                    !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
                     
         do jch = 2,nga                                              !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
            if (gen(whichgene)%e(jch) .ne. 0.0d0) then               !!>> HC 25-11-2020
               ord = ord + 1                                         !!>> HC 25-11-2020 order of vectors (1 to nume)
               plinkj(ord) = jch                                     !!>> HC 25-11-2020 indexes of regulated cell behaviors
            endif                                                    !!>> HC 25-11-2020
         enddo                                                       !!>> HC 25-11-2020
      
         call random_number(a)                                       !!>> HC 25-11-2020
         luckylink=ceiling(a*nume)                                   !!>> HC 25-11-2020 Randomly selecting 
         if (luckylink == 0) luckylink=1                             !!>> HC 25-11-2020
         
         runornot = 1 
         j = plinkj(luckylink)                                       !!>> HC 25-11-2020 Randomly selected behavior/property j
         mutagenj = j                                                        
         prevalue = gen(whichgene)%e(j)                                              
         
         call random_number(a)                                       !!>> HC 25-11-2020
         if(a < 1.0d-8) a = 0.02                                     !!>> AL 2-12-24 we don't want neutral mutations

         inch = (1-2*a)*gen(whichgene)%e(j)*mag_ismurate                             
         if ((gen(whichgene)%e(j) + inch) < min_elim(j)) then        !!>> AL 2-12-24
            if(gen(whichgene)%e(j) == min_elim(j)) then
               inch = -inch                                         
               gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch  
            else
               gen(whichgene)%e(j) = min_elim(j)
            end if
         else if ( (gen(whichgene)%e(j) + inch) > max_elim(j)) then 
            if(gen(whichgene)%e(j) == max_elim(j))then 
               inch = -inch
               gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch  
            else 
               gen(whichgene)%e(j) = max_elim(j)
            end if
         else 
            gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch
         end if
         newvalue = gen(whichgene)%e(j)                                               
      endif                             
   else                                                              !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!! mind == 2
     
      probs(1)=33d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
      probs(2)=33d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
      probs(3)=0d0                                                   !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now) !AL 2-12-24 this could be ereased
      probs(4)=32d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
      probs(5)=2d-2                                                  !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals

      call random_number(a)                                          !!>> HC 25-11-2020
      cum=0.0d0                                                      !!>> HC 25-11-2020
      do ich=1,size(probs)                                           !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
         cum=cum+probs(ich)                                          !!>> HC 25-11-2020 taking into account probabilites in probs
         if (a<cum)then                                              !!>> HC 25-11-2020
            tind=ich                                                 !!>> HC 25-11-2020
            exit                                                     !!>> HC 25-11-2020
         endif                                                       !!>> HC 25-11-2020
      enddo                                                          !!>> HC 25-11-2020
      
      if (tind .le. 2)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         mutacode  = 4                                               !!>> HC 29-11-2021 Store the matrix that will mutate
         realgenes = 0
         mutageni = whichgene                                        !!>> HC 29-11-2021 Store for output

         do jch = 1,ng                                               !!>> HC 25-11-2020
            if (int(gen(jch)%kindof) .eq. 9) cycle                   !!>> AL 29-11-24 ghost genes are discriminated
            realgenes = realgenes + 1
            if (gen(whichgene)%t(jch) .eq. 0.0d0) empty = empty + 1  !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
         end do                                                      !!>> HC 25-11-2020 and the number of already existing links in t matrix

         full = realgenes - empty                          !AL 2-12-24 important

         if ((tind == 1 .and. full > 0) .or. empty == 0) then        !!>> HC 25-11-2020 !!!!!! DELETION OF A LINK IN THE T MATRIX!!!!!! !!>> AL 29-11-24 the full>0 shouldnt be necessary. Check: if its needed, fix, if not, erease
            if (allocated(fullspotj)) deallocate(fullspotj)          !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                              !!>> HC 25-11-2020
            
            ord=0
            fullspotj=0                                              !!>> HC 25-11-2020 These vectors store the already existing links
            
            do jch = 1,ng                                            !!>> HC 25-11-2020
               if (int(gen(jch)%kindof) .eq. 9) cycle                !!>> AL 29-11-24 ghost genes are discriminated
               if (gen(whichgene)%t(jch) .ne. 0.0d0) then            !!>> HC 25-11-2020
                  ord = ord + 1                                      !!>> HC 25-11-2020 order of vectors (1 to full)
                  fullspotj(ord) = jch                               !!>> HC 25-11-2020 gene j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg = ceiling(a*ord)                                  !!>> HC 25-11-2020 Randomly pickin an interaction to delete
            if(luckyg == 0) luckyg = 1                               !!>> HC 25-11-2020

            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

            j = fullspotj(luckyg)                                    !!>> HC 25-11-2020
            mutagenj = j                                             !!>> HC 29-11-2021 Store for output
            prevalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output
            gen(whichgene)%t(j) = 0.0d0                              !!>> HC 25-11-2020
            newvalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output 

         else if(empty .ne. 0) then                                  !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
            if (allocated(freespotj)) deallocate(freespotj)          !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                             !!>> HC 25-11-2020
            
            ord=0;
            freespotj=0                                              !!>> HC 25-11-2020 These vectors store the free spots
            do jch = 1,ng                                            !!>> HC 25-11-2020
               if (int(gen(jch)%kindof) .eq. 9) cycle                !!>> AL 29-11-24 ghost genes are discriminated
               if (gen(whichgene)%t(jch) == 0.0d0) then              !!>> HC 25-11-2020
                  ord = ord + 1                                      !!>> HC 25-11-2020 order of the vectors (1 to empty)
                  freespotj(ord) = jch                               !!>> HC 25-11-2020 with gene j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg = ceiling(a*empty)                                !!>> HC 25-11-2020 Randomly picking a spot
            if(luckyg==0) luckyg=1                                   !!>> HC 25-11-2020
            call random_number(a)                                    !!>> HC 25-11-2020
            if(a < 1.0d-8) a = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations
            call random_number(b)                                    !!>> HC 25-11-2020
            if (b < 0.50d0) then
               signo = 1
            else
               signo = -1
            endif                                                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
            
            j = freespotj(luckyg)                                    !!>> HC 25-11-2020 Making new interaction
            mutagenj = j                                             !!>> HC 29-11-2021 Store for output
            
            prevalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output
            gen(whichgene)%t(j) = a*signo*w_max*t_mag                !!>> HC 31-5-2023 
            newvalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output 
            call functional_net(good)                                !AL 29-11-24 important
            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         else 
            print*,"We cannot add interactions since there are no positions left"
            go to 757
         endif                                                       !!>> HC 25-11-2020
     
      elseif (tind==3)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
      !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW !!>> HC 25-11-2020 
      !important
      elseif (tind==4)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e matrix)!!!!!!
         empty=0                                                     !!>> HC 25-11-2020
   
         do jch=2,nga                                                !!>> AL 2-12-24
            if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch !!>> AL 2-12-24 not sure what HC wanted to communicate but if rembeh= 1 that behaviour is off in ranges.dat so we dont consider it
            if (gen(whichgene)%e(jch)==0.0d0) then                   !!>> HC 27-11-2020
               empty = empty + 1                                     !!>> HC 27-11-2020 Number of empty spots in e Matrix
            else                                                     !!>> HC 27-11-2020
               full = full + 1                                       !!>> HC 27-11-2020 Number of existing interactions in e matrix
            end if                                                   !!>> HC 27-11-2020
         end do                                                      !!>> HC 25-11-2020

         call random_number(a)                                       !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
        
         if ((a < 0.50d0 .and. full > 1) .or. empty == 0)then        !!>> HC 10-5-2022 !!!REMOVE INTERACTION !!>> AL 28-1-25: note that full>1, otherwhise we can remove the only cell behaviour and end up without a developmental mechanism 
            runornot = 1                                             !!>> AL 10-4-24: modifying a cell behaviour always means altering functional network 
                                                                     !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(fullspotj)) deallocate(fullspotj)          !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                              !!>> HC 25-11-2020
            
            ord=0
            fullspotj=0                                              !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
            
            do jch=2,nga                                             !!>> HC 25-11-2020 we start in 2 because 1 is being an adhesion mol
               if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
               if (gen(whichgene)%e(jch).ne. 0.0d0)then              !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                  ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                  fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
         
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg=ceiling(a*real(full))                             !!>> HC 25-11-2020 Randomly picking an interaction to kill
            if(luckyg==0)luckyg=1                                    !!>> AL 3-12-2024   

            j=fullspotj(luckyg)                                      !!>> HC 25-11-2020
            
            !if(i==24 .and. j==37)then                               !!>> AL 28-1-25 important: we always have cell division consitutively being regulated by gene 24
            !   limit = 0
            !end if 
                                                                     !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=whichgene                                       !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2021 save for output the prev value
            newvalue=0.0d0                                           !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
            gen(whichgene)%e(j)=0.0d0                                !!>> HC 25-11-2020 we kill it easily

         else                                                        !!>> HC 25-11-2020 !ADD AN INTERACTION   
            runornot = 1                                             !!>> AL 10-4-24: modifying a cell behaviour always means altering functional newtwork 
                                                                     !!>> AL 10-4-24: so we run development to see how it affects this pal

            if (allocated(freespotj)) deallocate(freespotj)          !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                             !!>> HC 25-11-2020
            
            ord=0
            freespotj=0                                              !!>> HC 25-11-2020 
            do jch = 2,nga                                           !!>> HC 25-11-2020
               if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
               if (gen(whichgene)%e(jch) == 0.0d0)then               !!>> HC 25-11-2020
                  ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                  freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                  !!>> HC 25-11-2020 Randomly picking an interaction
            if(luckyg==0)luckyg=1                                    !!>> HC 25-11-2020
            j=freespotj(luckyg)                                      !!>> HC 25-11-2020 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=whichgene                                       !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2021 save for output the prev value
            
            call random_number(c)                                    !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
            if(c < 1.0d-8) c = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations

            anch=max_elim(j)-min_elim(j)                             !!>> HC 28-11-2020 we introduce an interaction in the range of the limits            
            gen(whichgene)%e(j)=(anch*c+min_elim(j))*t_mag           !!>> AL 3-12-24
            !gen(whichgene)%e(j)=gen(whichgene)%e(j)*t_mag           !!>> HC 31-5-2023  Correct the magnitude of the new interaction 

            if (j==nparam_per_node+8)then                            !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation !>> AL 2-12-24 nparam_per_node = 35 PCP in the papers called PCO for whatever reason 
               call random_number(c)                                 !!>> AL 2-12-24 
               if(c < 1.0d-8) c = 0.01                               !!>> AL 2-12-24 we don't want neutral mutations
               
               gen(24)%e(nparam_per_node+17) = (0.20d0*c)*t_mag      !!>> AL 2-12-24 we assign a random magnitud with respecto to max val (0.20d0)
               
               !pcpval = 0.20d0*c 
               !gen(ng)%e(nparam_per_node+17) =0.20d0                !!>> HC 12-11-2020 !!!ATENTION!! gene ng is homogeneously expressed in the blastula ics !>> AL 2-12-24 important: this shouldnt start at the max possible value..
            endif                                                    !!>> HC 7-12-2020 !!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
            
            pcpval = gen(24)%e(nparam_per_node+17)                   !!>> AL 10-1-25 
            newvalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2020  Store the new value
         endif                                                       !!>> HC 25-11-2020
      
      elseif (tind==5)then                                           !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
   
         mutacode=7                                                  !!>> HC 29-11-2021 Store that this is a gene duplication/delection
         call random_number(a)                                       !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
         if (a < 0.50d0) then                                        !!>> AL 17-1-25 Deletion   
            if(ng < 2)then 
               go to 12345                                           !!>> AL 17-1-25 If the deletion will erease the genome we do a duplication instead 
            end if 

            mutageni = whichgene 
            mutagenj = -1                                            !!>> HC 29-11-2021 store that this is a delection
            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
            call deletion_rec(whichgene)                             !!>> AL 29-11-24

      else                                                           !!>> HC 25-11-2020 !!>> AL 17-1-25 Duplication
      
      12345 continue
            lch = 0
            do jch = 1,ng                                            !!>> AL 17-1-25 we have to check if there are ghost genes available
               if (int(gen(jch)%kindof) .eq. 9) then 
                  exit
               end if
               lch = lch + 1
            end do 

            if(lch .eq. ng)then                                      !!>> AL 17-1-25 there are no ghost genes left so no gen duplications are possible      
               errror = 2
               go to 757
            end if

            if(errror .eq. 2)then   
               mutageni=0 
               mutagenj=1
            else 
               mutageni=whichgene 
               mutagenj=1                                            !!>> HC 29-11-2021 store that this is a duplication
               if(good(whichgene) == 1) runornot = 1                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
               call duplication_rec(whichgene)                       !!>> AL 29-11-24
            end if
         endif                                                       !!>> HC 25-11-2020
      end if 
   endif                                                             !!>> HC 25-11-2020
 
   if (errror == 1 .or. errror == 2) then
757   continue   
      inviable=1  
   end if  

   inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
   if (ng<1) then                                                    !!>> HC 25-11-2020 
      print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   else                                                              !!>> HC 25-11-2020
      do j=1,ng                                                      !!>> HC 25-11-2020 !!>> AL 2-12-24 important, this might not do what it should
         if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
            inviable=0                                               !!>> HC 25-11-2020
            goto 37                                                  !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   end if                                                            !!>> HC 25-11-2020

37 continue      

   pcpval = gen(24)%e(nparam_per_node+17)                            !!>> AL 3-3-25

end subroutine suremuta_no_kadh_functional_net_rec_2

subroutine suremuta_no_kadh_functional_net_rec_3(whichgene,limit,inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue &
                                                ,rangfile,runornot,pcpval) !!>> AL 14-2-25
   !!>> AL 14-2-25 in this function you already chose wich gene to mutate
   
   !!>> AL 9-4-24: this function determines whether the mutation alters the functional part of the gen network of an individual or not.
   !!>> AL 9-4-24: if not then development of mutated individual is not runned and fitness is inherited from parent.
   !!>> HC 25-11-2020 This subroutine was created by Hugo Cano to make sure mutations in an imput file 
   !!>> HC 14-4-2023 This version is the same but ignoring specific adhesion molecules
   !!>> HC 25-11-2020 if mind=1 we do a IS mutation and if mind=2 we do a T mutation
   !!>> HC 25-11-2020 when doing and IS mutation, all the parameters have the same probability to mutate
   !!>> HC 25-11-2020 but when we do a T mutation each kind of mutation has a different probability to occur
   !!>> HC 25-11-2020 (i.e. gain of function of transcription factor, loss of function of transcription factor,
   !!>> HC 25-11-2020  gain/loss function regulation cell activites and gene delection7duplication)
   
   implicit none

   integer                     :: limit,inviable,mind,tind,isitreal,realgenes
   integer                     :: props,nume,numt,numr,numk          !!>> HC 24-11-2020 number of e,t,nww,kadh interaction for IS mutations and total of possible IS mutations
   integer                     :: luckyprop, luckylink, luckyg       !!>> HC 24-11-2020 selected gene, interaction or property
   integer                     :: ich,jch,lch,ord                    !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
   integer                     :: full,empty                         !!>> HC 24-11-2020 full an empty spaces for T mutations
   real*8, dimension(1:5)      :: probs                              !!>> HC 24-11-2020 relative probabilities of different T mutations
   !real*8, dimension(1:7)      :: probis                            !!>> HC 24-11-2020 relative probabilities of different IS mutations
   real*8, dimension(5)        :: probis                             !AL 2-12-24 this is 5 since there are no reactions 
   real*8, dimension(1:nga)    :: max_elim, min_elim                 !!>> HC 27-11-2020 These vectors store the limits of e matrix 

   integer, dimension(1:nga)   :: rembeh                             !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
   real*8                      :: cum                                !!>> HC 24-11-2020 funny position to store cumulative Fi
   integer                     :: signo                              !!>> HC 24-11-2020 sign of the mutation  
   real*8                      :: inch, newval, anch                 !!>> HC 27-11-2020 actual magnitude of the IS mutation and new value of the parameter
   integer                     :: mutacode,mutageni,mutagenj,whichgene  !!>> HC 16-9-2021 This stores the kind of mutation that we had 
   real*8                      :: prevalue, newvalue
   character*400               :: rangfile                           !!>> 6-10-2021 file where the ranges are stored
   real*8, dimension (1:5)     :: min_glim, max_glim                 !!>>HC 6-10-2021
   real*8                      :: t_mag,pcpval                       !!>> HC 31-5-2023 Magnitude of the initial value in T mutations (should be a free parameter) !!>> AL 28-12-2024 pcpval = magnitude of planar cell contraction 
   integer                     :: errror                             !!>>AL 17-1-25
   real*8                      :: mean_norm,stddev_norm,beta               !!>>AL 13-3-25   
   integer,allocatable,dimension(:)  :: plinki,plinkj                !>> HC 24-11-2020 Vectors containing present interactions (e,t,nww or kadh) for IS
   integer,allocatable,dimension(:)  :: freespoti,freespotj,fullspoti,fullspotj !!>> HC 24-11-2020 Vectors containing present and absent interactions (e,t,nww or kadh) fot TM+
   
   !functional subnetwork                                            !!>>AL 9-4-24
   integer :: runornot                                               !!>>AL 10-4-24: if == 1 you run development, if 0 you do not          
   integer, allocatable, dimension(:) :: good                        !!>>AL 10-4-24: this vector has elements = 1 when genes are functional

   if(allocated(good)) deallocate(good)                              !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
   allocate(good(1:ng))                                     

   pcpval      = 0.0d0                                               !!>>AL 28-12-24:
   errror      = 0                                                   !!>>AL 17-1-25:
   mean_norm   = 0                                                   !! AL 13-3-25
   stddev_norm = sqrt(0.1)                                           !! AL 13-3-25
   
   !limits and ranges        
   limit=1                                                           !!>> HC 27-11-2020 if limit stays in 1 no limits have been surpased if =0, we will repeat mutation
   rembeh=0          
   max_elim=0.0d0 ; min_elim=0.0d0                                   !!>> HC 6-10-2020 These vectors will store the max and min values of activation of cell properties and behaviors
   max_glim=0.0d0 ; min_glim=0.0d0                                   !!>> HC 6-10-2020 Same for other gene properties 
   
   call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2020 the maximum values are read from an expernal file
   
   d_max=max_glim(1)    ; d_min=min_glim(1)                          !!>> HC 27-11-2020 Limits for diffusion rate (diffu) maximum value from Hagolani et al. 2020
   m_max=max_glim(2)    ; m_min=min_glim(2)                          !!>> HC 27-11-2020 Limits for degradation rate (mu) from Hagolani et al. 2020
   mich_max=max_glim(3) ; mich_min=min_glim(3)                       !!>> HC 27-11-2020 Limits for Michaelis-Menten constant (mich)
   w_max=max_glim(4)    ; w_min=min_glim(4)                          !!>> HC 27-11-2020 Limits for transcription factors from Hagolani et al. 2020
   b_max=max_glim(5)    ; b_min=min_glim(5)                          !!>> HC 27-11-2020 Limits for specific adhesion (kadh AKA b matrix)

   call functional_net(good)                                         !!>>AL 29-11-24

   empty    = 0
   full     = 0                                                      !!>> HC 24-11-2020
   mutacode = 0                                                      !!>> HC 16-9-2021  
   mutageni = 0                                                      !!>> HC 29-11-2021
   mutagenj = 0                                                      !!>> HC 29-11-2021
   realgenes= 0                                                      !! AL 2-12-24
   inch     = 0.0d0           
   anch     = 0.0d0                                                  !!>> HC 28-11-2020
   prevalue = 0.0d0                                                  !!>> HC 29-11-2021
   newvalue = 0.0d0                                                  !!>> HC 29-11-2021
   t_mag    = 0.20d0                                                 !!AL 2-12-24 important: this shouldnt be defined manually but in paraopc.par
   
   if (mind == 1) then                                               !!>> HC 25-11-2020 !******************* IS MUTATIONS***************************!
      
      nume=0; numt=0; numk=0; numr=0                       

      do jch=1,ng                                       
         if(gen(whichgene)%t(jch) .ne. 0.0d0) numt = numt + 1        !!>> HC 25-11-2020 Number of existing links in the t matrix (transcription factors)
      enddo                                              
      
      do jch=2,nga                                                   !!>> HC 25-11-2020 !e(1), can not be changed since it is the index in the kadh matrix
         if(gen(whichgene)%e(jch) .ne. 0.0d0) nume = nume + 1        !!>> HC 25-11-2020 Number of existing links in the e matrix (cell behaviors)
      enddo                                             

      do jch=1,ng                                       
         if(int(gen(jch)%kindof) .ne. 9) realgenes = realgenes + 1   !AL 21-1-25 finding number of real genes
      enddo        

      probis = 0.0d0                                                 !!>> HC 25-11-2020 props is the total number of parameters that can change

      props = 3 + numt + nume                                        !AL 2-12-24 

      probis(1) = 1/real(props)                                      !!>> AL: 14-2-25 probability mu     
      if(int(gen(whichgene)%kindof) .ne. 4)then                     !!>> AL: 16-4-25 diff. coef. only changes if its an extrecellular signal
         probis(2) = 0
      else 
         probis(2) = 1/real(props)                                   !!>> AL: 14-2-25 probability diff
      end if
      probis(3) = 1/real(props)                                      !!>> AL: 14-2-25 probability mich
      probis(4) = real(numt)/real(props)                             !!>> HC 25-11-2020 probability t matrix
      probis(5) = real(nume)/real(props)                             !!>> HC 25-11-2020 probability e matrix
      
      cum = 0.0d0                                                    !!>> HC 25-11-2020
      call random_number(a)                                          !!>> HC 25-11-2020

      do ich = 1, size(probis)                                       !!>> HC 25-11-2020 This decides randomly what matrix is going to mutate
         cum = cum + probis(ich)                                     !!>> HC 25-11-2020 taking into acount that the number of parameters in each matrix
         if (a < cum) then                                           !!>> HC 25-11-2020
            luckyprop = ich                                          !!>> HC 25-11-2020
            exit                                                     !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      
      mutacode = luckyprop                                           !!>> HC 16-9-2021  Store mutation code for output   

      mutageni = whichgene

      if (good(whichgene) == 1) runornot = 1                         !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

      if (luckyprop == 1)then                                        !!>> HC 25-11-2020 !***********IS MUTATION IN DEGRADATION RATE (mu)**************!
         
         prevalue = gen(whichgene)%mu                                !!>> HC 29-11-2021 We store the old value
         
         beta = gaussian_random(mean_norm,stddev_norm)               !! AL 13-3-25 generate a randnum from N(0,sqrt(0.1)). -1 <= beta <= +1
         if(beta < -1) beta = -1  
         if(beta > +1) beta = +1

         inch = (m_max - m_min)*mag_ismurate*beta                    !! AL 13-3-25

         if((gen(whichgene)%mu + inch) < m_min)then                 !!>> AL 2-12-24 we always want mutations to change the selected parameters
            if(gen(whichgene)%mu == m_min)then
               inch = -inch                                         
               gen(whichgene)%mu = gen(whichgene)%mu+inch  
            else
               gen(whichgene)%mu = m_min 
            end if
         else if ((gen(whichgene)%mu+inch) > m_max) then 
            if(gen(whichgene)%mu == m_max)then 
               inch = -inch
               gen(whichgene)%mu = gen(whichgene)%mu + inch
            else 
               gen(whichgene)%mu = m_max
            end if
         else 
            gen(whichgene)%mu = gen(whichgene)%mu + inch
         end if

         newvalue = gen(whichgene)%mu                                !!>> HC 29-11-2021  we store the new value
         
      elseif (luckyprop==2)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN DIFFUSION COEFFICIENT (diffu)!!!!!!

         prevalue = gen(whichgene)%diffu                             !!>> HC 29-11-2021 We store the old value
               
         beta = gaussian_random(mean_norm,stddev_norm)             !! AL 13-3-25 generate a randnum from N(0,sqrt(0.1)). -1 <= beta <=+1
         if(beta < -1) beta = -1  
         if(beta > +1) beta = +1

         inch = (d_max - d_min)*mag_ismurate*beta                 !! AL 13-3-25

         if((gen(whichgene)%diffu + inch) < d_min)then               !!>> AL 2-12-24
            if(gen(whichgene)%diffu == d_min)then
               inch = -inch                                         
               gen(whichgene)%diffu = gen(whichgene)%diffu+inch  
            else
               gen(whichgene)%diffu = d_min 
            end if
         else if((gen(whichgene)%diffu + inch) > d_max)then 
            if(gen(whichgene)%diffu == d_max)then 
               inch = -inch
               gen(whichgene)%diffu = gen(whichgene)%diffu + inch
            else 
               gen(whichgene)%diffu = d_max
            end if
         else 
            gen(whichgene)%diffu = gen(whichgene)%diffu + inch
         end if
         
         newvalue = gen(whichgene)%diffu                             !!>> HC 29-11-2021  we store the new value      

      elseif (luckyprop==3)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN MICHAELIS-MENTEN CONSTANT (mich)!!!!!!

         prevalue = gen(whichgene)%mich                              !!>> HC 29-11-2021 We store the old value
      
         beta = gaussian_random(mean_norm,stddev_norm)             !! AL 13-3-25 generate a randnum from N(0,sqrt(0.1)). -1 <= beta <=+1
         if(beta < -1) beta = -1  
         if(beta > +1) beta = +1

         inch = (mich_max - mich_min)*mag_ismurate*beta           !! AL 13-3-25
           
         if ((gen(whichgene)%mich + inch) < mich_min)then            !!>> AL 2-12-24
            if(gen(whichgene)%mich == mich_min)then
               inch = -inch                                         
               gen(whichgene)%mich = gen(whichgene)%mich + inch  
            else
               gen(whichgene)%mich = mich_min 
            end if
         else if ( (gen(whichgene)%mich + inch) > mich_max) then 
            if(gen(whichgene)%mich == mich_max)then 
               inch = -inch
               gen(whichgene)%mich = gen(whichgene)%mich + inch
            else 
               gen(whichgene)%mich = mich_max
            end if
         else 
            gen(whichgene)%mich = gen(whichgene)%mich + inch
         end if

         newvalue = gen(whichgene)%mich                              !!>> HC 29-11-2021  we store the new value  
      
      elseif (luckyprop==4)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         if(allocated(plinkj)) deallocate(plinkj)                    !!>> HC 25-11-2020 These vectors store the number of links in the t matrix
         allocate(plinkj(1:numt))                                    !!>> HC 25-11-2020
         
         ord   = 0
         plinkj= 0                                           
        
         do jch=1,ng                                                    
            if(int(gen(jch)%kindof) .eq. 9) cycle                    !!>> AL 29-11-24 ghost genes are discriminated
            if (gen(whichgene)%t(jch) .ne. 0.0d0) then                          
               ord = ord + 1                                                
               plinkj(ord) = jch                                         
            end if                                                       
         end do                                                      !!>> HC 25-11-2020
                 
         call random_number(a)                                       !!>> HC 25-11-2020
         luckylink = ceiling(a*numt)                                 !!>> HC 25-11-2020 We choose randomly one transcription interaction in matrix T
         if (luckylink == 0) luckylink = 1                           !!>> HC 25-11-2020
        
         j = plinkj(luckylink)                                       !!>> HC 25-11-2020 gene j is the randomly selectied regulated gene
         mutagenj = j                                                !!>> HC 29-11-2021
         
         prevalue = gen(whichgene)%t(j)                              !!>> HC 29-11-2021 store the previous value
  
         beta = gaussian_random(mean_norm,stddev_norm)                !! AL 13-3-25 generate a randnum from N(0,sqrt(0.1)). -1 <= beta <=+1
         if(beta < -1) beta = -1  
         if(beta > +1) beta = +1

         inch = (w_max - w_min)*mag_ismurate*beta                    !! AL 13-3-25   

         if ((gen(whichgene)%t(j) + inch) < w_min)then               !!>> AL 2-12-24
            if(gen(i)%t(j) == w_min)then
               inch = -inch                                         
               gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch  
            else
               gen(whichgene)%t(j) = w_min
            end if

         else if ((gen(whichgene)%t(j) + inch) > w_max) then 
            if(gen(whichgene)%t(j) == w_max)then
               inch = -inch                                         
               gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch  
            else
               gen(whichgene)%t(j) = w_max
            end if
         else 
            gen(whichgene)%t(j) = gen(whichgene)%t(j) + inch
         end if
         
         newvalue = gen(whichgene)%t(j)                                !!>> HC 29-11-2021

      elseif (luckyprop==5)then                                      !!>> HC 25-11-2020 !!!!!!IS MUTATION REGULATION BEHAVIORS/PROPERTIES (e matrix)!!!!!!
      
         if(allocated(plinkj)) deallocate(plinkj)                    !!>> HC 25-11-2020
         allocate(plinkj(1:nume))                                    !!>> HC 25-11-2020
      
         ord=0
         plinkj=0                                                    !!>> HC 25-11-2020 These vectors store the number of links in the e matrix
                     
         do jch = 2,nga                                              !!>> HC 25-11-2020 There can be no IS mutations in being an adhesion molecule e(1)
            if (gen(whichgene)%e(jch) .ne. 0.0d0) then               !!>> HC 25-11-2020
               ord = ord + 1                                         !!>> HC 25-11-2020 order of vectors (1 to nume)
               plinkj(ord) = jch                                     !!>> HC 25-11-2020 indexes of regulated cell behaviors
            endif                                                    !!>> HC 25-11-2020
         enddo                                                       !!>> HC 25-11-2020
      
         call random_number(a)                                       !!>> HC 25-11-2020
         luckylink=ceiling(a*nume)                                   !!>> HC 25-11-2020 Randomly selecting 
         if (luckylink == 0) luckylink=1                             !!>> HC 25-11-2020
         
         runornot = 1 
         j = plinkj(luckylink)                                       !!>> HC 25-11-2020 Randomly selected behavior/property j
         mutagenj = j                                                        
         prevalue = gen(whichgene)%e(j)                                              

         beta = gaussian_random(mean_norm,stddev_norm)                !! AL 13-3-25 generate a randnum from N(0,sqrt(0.1)). -1 <= beta <=+1
         if(beta < -1) beta = -1  
         if(beta > +1) beta = +1
         inch = (max_elim(j) - min_elim(j))*mag_ismurate*beta        !! AL 13-3-25

         if ((gen(whichgene)%e(j) + inch) < min_elim(j)) then         !!>> AL 2-12-24
            if(gen(whichgene)%e(j) == min_elim(j)) then
               inch = -inch                                         
               gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch  
            else
               gen(whichgene)%e(j) = min_elim(j)
            end if
         else if ( (gen(whichgene)%e(j) + inch) > max_elim(j)) then 
            if(gen(whichgene)%e(j) == max_elim(j))then 
               inch = -inch
               gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch  
            else 
               gen(whichgene)%e(j) = max_elim(j)
            end if
         else 
            gen(whichgene)%e(j) = gen(whichgene)%e(j) + inch
         end if
         newvalue = gen(whichgene)%e(j)                                               
      endif                             
   else                                                              !!>> HC 25-11-2020 !!!!!!!!!!!!!!!!!!!!!T MUTATIONS!!!!!!!!!!!!!!!!!!!!! mind == 2
     
      probs(1)=33d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations with loss of function in the mutated individuals
      probs(2)=33d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations with gain of function in the mutated individuals
      probs(3)=0d0                                                   !!>> HC 23-11-2020 Proportion of TM mutations in post transcriptional reactions (we do not have these right now) !AL 2-12-24 this could be ereased
      probs(4)=32d-2                                                 !!>> HC 20-11-2020 Proportion of TM mutations in regulation of cell behaviors in the mutated individuals
      probs(5)=2d-2                                                  !!>> HC 20-11-2020 Proportion of duplications/delections in the mutated individuals

      call random_number(a)                                          !!>> HC 25-11-2020
      cum=0.0d0                                                      !!>> HC 25-11-2020
      do ich=1,size(probs)                                           !!>> HC 25-11-2020 Randomly deciding what kind of T mutation we are going to have
         cum=cum+probs(ich)                                          !!>> HC 25-11-2020 taking into account probabilites in probs
         if (a<cum)then                                              !!>> HC 25-11-2020
            tind=ich                                                 !!>> HC 25-11-2020
            exit                                                     !!>> HC 25-11-2020
         endif                                                       !!>> HC 25-11-2020
      enddo                                                          !!>> HC 25-11-2020
      
      if (tind .le. 2)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION IN TRANSCRIPTION FACTORS (t matrix)!!!!!!
         
         mutacode  = 4                                               !!>> HC 29-11-2021 Store the matrix that will mutate
         realgenes = 0
         mutageni = whichgene                                        !!>> HC 29-11-2021 Store for output

         do jch = 1,ng                                               !!>> HC 25-11-2020
            if (int(gen(jch)%kindof) .eq. 9) cycle                   !!>> AL 29-11-24 ghost genes are discriminated
            realgenes = realgenes + 1
            if (gen(whichgene)%t(jch) .eq. 0.0d0) empty = empty + 1  !!>> HC 25-11-2020 We find the number of empty spots in the t matrix
         end do                                                      !!>> HC 25-11-2020 and the number of already existing links in t matrix

         full = realgenes - empty                                    !AL 2-12-24 important

         if ((tind == 1 .and. full > 0) .or. empty == 0) then        !!>> HC 25-11-2020 !!!!!! DELETION OF A LINK IN THE T MATRIX!!!!!! !!>> AL 29-11-24 the full>0 shouldnt be necessary. Check: if its needed, fix, if not, erease
            if (allocated(fullspotj)) deallocate(fullspotj)          !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                              !!>> HC 25-11-2020
            
            ord=0
            fullspotj=0                                              !!>> HC 25-11-2020 These vectors store the already existing links
            
            do jch = 1,ng                                            !!>> HC 25-11-2020
               if (int(gen(jch)%kindof) .eq. 9) cycle                !!>> AL 29-11-24 ghost genes are discriminated
               if (gen(whichgene)%t(jch) .ne. 0.0d0) then            !!>> HC 25-11-2020
                  ord = ord + 1                                      !!>> HC 25-11-2020 order of vectors (1 to full)
                  fullspotj(ord) = jch                               !!>> HC 25-11-2020 gene j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg = ceiling(a*ord)                                  !!>> HC 25-11-2020 Randomly pickin an interaction to delete
            if(luckyg == 0) luckyg = 1                               !!>> HC 25-11-2020

            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal

            j = fullspotj(luckyg)                                    !!>> HC 25-11-2020
            mutagenj = j                                             !!>> HC 29-11-2021 Store for output
            prevalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output
            gen(whichgene)%t(j) = 0.0d0                              !!>> HC 25-11-2020
            newvalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output 

         else if(empty .ne. 0) then                                  !!>> HC 25-11-2020 !!!!!! ADDING OF A LINK IN THE T MATRIX!!!!!!
            if (allocated(freespotj)) deallocate(freespotj)          !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                             !!>> HC 25-11-2020
            
            ord=0;
            freespotj=0                                              !!>> HC 25-11-2020 These vectors store the free spots
            do jch = 1,ng                                            !!>> HC 25-11-2020
               if (int(gen(jch)%kindof) .eq. 9) cycle                !!>> AL 29-11-24 ghost genes are discriminated
               if (gen(whichgene)%t(jch) == 0.0d0) then              !!>> HC 25-11-2020
                  ord = ord + 1                                      !!>> HC 25-11-2020 order of the vectors (1 to empty)
                  freespotj(ord) = jch                               !!>> HC 25-11-2020 with gene j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg = ceiling(a*empty)                                !!>> HC 25-11-2020 Randomly picking a spot
            if(luckyg==0) luckyg=1                                   !!>> HC 25-11-2020
            call random_number(a)                                    !!>> HC 25-11-2020
            if(a < 1.0d-8) a = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations
            call random_number(b)                                    !!>> HC 25-11-2020
            if (b < 0.50d0) then
               signo = 1
            else
               signo = -1
            endif                                                    !!>> HC 25-11-2020 Deciding randomly the sign of the interaction
            
            j = freespotj(luckyg)                                    !!>> HC 25-11-2020 Making new interaction
            mutagenj = j                                             !!>> HC 29-11-2021 Store for output
            
            prevalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output
            gen(whichgene)%t(j) = a*signo*w_max*t_mag                !!>> HC 31-5-2023 
            newvalue = gen(whichgene)%t(j)                           !!>> HC 29-11-2021 Store for output 
            call functional_net(good)                                !AL 29-11-24 important
            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
         else 
            print*,"We cannot add interactions since there are no positions left"
            go to 757
         endif                                                       !!>> HC 25-11-2020
     
      elseif (tind==3)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION POST TRANSCRIPTIONAL REACTIONS (r matrix)!!!!!!
      !POST TRANSCRIPTIONAL REACTIONS, NO T MUTATIONS HERE RIGHT NOW !!>> HC 25-11-2020 
      !important
      elseif (tind==4)then                                           !!>> HC 25-11-2020 !!!!!! T MUTATION CELL BEHAVIORS/PROPERTIES (e matrix)!!!!!!
         empty=0                                                     !!>> HC 25-11-2020
   
         do jch=2,nga                                                !!>> AL 2-12-24
            if (rembeh(jch)==1)cycle                                 !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch !!>> AL 2-12-24 not sure what HC wanted to communicate but if rembeh= 1 that behaviour is off in ranges.dat so we dont consider it
            if (gen(whichgene)%e(jch)==0.0d0) then                   !!>> HC 27-11-2020
               empty = empty + 1                                     !!>> HC 27-11-2020 Number of empty spots in e Matrix
            else                                                     !!>> HC 27-11-2020
               full = full + 1                                       !!>> HC 27-11-2020 Number of existing interactions in e matrix
            end if                                                   !!>> HC 27-11-2020
         end do                                                      !!>> HC 25-11-2020

         call random_number(a)                                       !!>> HC 25-11-2020 50% chances we remove or 50% we add an interaction
        
         if ((a < 0.50d0 .and. full > 1) .or. empty == 0)then        !!>> HC 10-5-2022 !!!REMOVE INTERACTION !!>> AL 28-1-25: note that full>1, otherwhise we can remove the only cell behaviour and end up without a developmental mechanism 
            runornot = 1                                             !!>> AL 10-4-24: modifying a cell behaviour always means altering functional network 
                                                                     !!>> AL 10-4-24: so we run development to see how it affects this pal
            if (allocated(fullspotj)) deallocate(fullspotj)          !!>> HC 25-11-2020
            allocate(fullspotj(1:full))                              !!>> HC 25-11-2020
            
            ord=0
            fullspotj=0                                              !!>> HC 25-11-2020 These vectors will store the existing interactions and whether they are not in kadh
            
            do jch=2,nga                                             !!>> HC 25-11-2020 we start in 2 because 1 is being an adhesion mol
               if (rembeh(jch)==1)cycle                              !!>> HC 10-5-2022
               if (gen(whichgene)%e(jch).ne. 0.0d0)then              !!>> HC 25-11-2020 Here we do not need filtering unused cell behaviors/properties because we assume they are not
                  ord=ord+1                                          !!>> HC 25-11-2020 in the original file
                  fullspotj(ord)=jch                                 !!>> HC 25-11-2020 property j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
         
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg=ceiling(a*real(full))                             !!>> HC 25-11-2020 Randomly picking an interaction to kill
            if(luckyg==0)luckyg=1                                    !!>> AL 3-12-2024   

            j=fullspotj(luckyg)                                      !!>> HC 25-11-2020
            
            !if(i==24 .and. j==37)then                               !!>> AL 28-1-25 important: we always have cell division consitutively being regulated by gene 24
            !   limit = 0
            !end if 
                                                                     !!>> HC 29-11-2021 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=whichgene                                       !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2021 save for output the prev value
            newvalue=0.0d0                                           !!>> HC 29-11-2021 store that the new value should be 0 (delection)    
            gen(whichgene)%e(j)=0.0d0                                !!>> HC 25-11-2020 we kill it easily

         else                                                        !!>> HC 25-11-2020 !ADD AN INTERACTION   
            runornot = 1                                             !!>> AL 10-4-24: modifying a cell behaviour always means altering functional newtwork 
                                                                     !!>> AL 10-4-24: so we run development to see how it affects this pal

            if (allocated(freespotj)) deallocate(freespotj)          !!>> HC 25-11-2020
            allocate(freespotj(1:empty))                             !!>> HC 25-11-2020
            
            ord=0
            freespotj=0                                              !!>> HC 25-11-2020 
            do jch = 2,nga                                           !!>> HC 25-11-2020
               if (rembeh(jch)==1)cycle                              !!>> HC 6-10-2021  if this is 1 we are considering mutations in the prop/beh jch
               if (gen(whichgene)%e(jch) == 0.0d0)then               !!>> HC 25-11-2020
                  ord=ord+1                                          !!>> HC 25-11-2020 order of the vector (1 to nume)
                  freespotj(ord)=jch                                 !!>> HC 25-11-2020 cell behafior/property j
               end if                                                !!>> HC 25-11-2020
            end do                                                   !!>> HC 25-11-2020
            
            call random_number(a)                                    !!>> HC 25-11-2020
            luckyg=ceiling(a*empty)                                  !!>> HC 25-11-2020 Randomly picking an interaction
            if(luckyg==0)luckyg=1                                    !!>> HC 25-11-2020
            j=freespotj(luckyg)                                      !!>> HC 25-11-2020 SAVE THE INFORMATION IN THE OUTPUT --START--
            mutacode=5                                               !!>> HC 29-11-2021 The mutation is in the E matrix
            mutageni=whichgene                                       !!>> HC 29-11-2021 Store positions affected
            mutagenj=j                                               !!>> HC 29-11-2021 Store positions affected
            
            prevalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2021 save for output the prev value
            
            call random_number(c)                                    !!>> HC 25-11-2020 Magnitude, anch saves the lenght of the valid interval
            if(c < 1.0d-8) c = 0.01                                  !!>> AL 2-12-24 we don't want neutral mutations

            anch=max_elim(j)-min_elim(j)                             !!>> HC 28-11-2020 we introduce an interaction in the range of the limits            
            gen(whichgene)%e(j)=(anch*c+min_elim(j))*t_mag           !!>> AL 3-12-24
            !gen(whichgene)%e(j)=gen(whichgene)%e(j)*t_mag           !!>> HC 31-5-2023  Correct the magnitude of the new interaction 

            if (j==nparam_per_node+8)then                            !!>> HC 7-12-2020 !!ATENTION!! this is a trick to make PCP more likely to appear by mutation !>> AL 2-12-24 nparam_per_node = 35 PCP in the papers called PCO for whatever reason 
               call random_number(c)                                 !!>> AL 2-12-24 
               if(c < 1.0d-8) c = 0.01                               !!>> AL 2-12-24 we don't want neutral mutations
               
               gen(24)%e(nparam_per_node+17) = (0.20d0*c)*t_mag      !!>> AL 2-12-24 we assign a random magnitud with respecto to max val (0.20d0)
               
               !pcpval = 0.20d0*c 
               !gen(ng)%e(nparam_per_node+17) =0.20d0                !!>> HC 12-11-2020 !!!ATENTION!! gene ng is homogeneously expressed in the blastula ics !>> AL 2-12-24 important: this shouldnt start at the max possible value..
            endif                                                    !!>> HC 7-12-2020 !!ATENTION!! nparam_+er node +8 and +17 are needed for PCP to happen
            
            pcpval = gen(24)%e(nparam_per_node+17)                   !!>> AL 10-1-25 
            newvalue=gen(whichgene)%e(j)                             !!>> HC 29-11-2020  Store the new value
         endif                                                       !!>> HC 25-11-2020
      
      elseif (tind==5)then                                           !!>> HC 25-11-2020 GENE DELECTION/DUPLICATION
   
         mutacode=7                                                  !!>> HC 29-11-2021 Store that this is a gene duplication/delection
         call random_number(a)                                       !!>> HC 25-11-2020 50% chances we remove the gene 50% we duplicate it
         if (a < 0.50d0) then                                        !!>> AL 17-1-25 Deletion   
            if(ng < 2)then 
               go to 12345                                           !!>> AL 17-1-25 If the deletion will erease the genome we do a duplication instead 
            end if 

            mutageni = whichgene 
            mutagenj = -1                                            !!>> HC 29-11-2021 store that this is a delection
            if(good(whichgene) == 1) runornot = 1                    !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
            call deletion_rec(whichgene)                             !!>> AL 29-11-24

      else                                                           !!>> HC 25-11-2020 !!>> AL 17-1-25 Duplication
      
      12345 continue
            lch = 0
            do jch = 1,ng                                            !!>> AL 17-1-25 we have to check if there are ghost genes available
               if (int(gen(jch)%kindof) .eq. 9) then 
                  exit
               end if
               lch = lch + 1
            end do 

            if(lch .eq. ng)then                                      !!>> AL 17-1-25 there are no ghost genes left so no gen duplications are possible      
               errror = 2
               go to 757
            end if

            if(errror .eq. 2)then   
               mutageni=0 
               mutagenj=1
            else 
               mutageni=whichgene 
               mutagenj=1                                            !!>> HC 29-11-2021 store that this is a duplication
               if(good(whichgene) == 1) runornot = 1                 !!>> AL 10-4-24: mutation in functional part, run development to see how it affects this pal
               call duplication_rec(whichgene)                       !!>> AL 29-11-24
            end if
         endif                                                       !!>> HC 25-11-2020
      end if 
   endif                                                             !!>> HC 25-11-2020
 
   if (errror == 1 .or. errror == 2) then
757   continue   
      inviable=1  
   end if  

   inviable=0                                                        !!>> HC 25-11-2020 evaluates if the individual is inviable
   if (ng<1) then                                                    !!>> HC 25-11-2020 
      print *,"inviable individual because there are no genes left"  !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   else                                                              !!>> HC 25-11-2020
      do j=1,ng                                                      !!>> HC 25-11-2020 !!>> AL 2-12-24 important, this might not do what it should
         if (gen(j)%kindof<3) then                                   !!>> HC 25-11-2020
            inviable=0                                               !!>> HC 25-11-2020
            goto 37                                                  !!>> HC 25-11-2020
         end if                                                      !!>> HC 25-11-2020
      end do                                                         !!>> HC 25-11-2020
      print *,"no gene is transcribable so this is inviable"         !!>> HC 25-11-2020
      inviable=1                                                     !!>> HC 25-11-2020
   end if                                                            !!>> HC 25-11-2020

37 continue      

   pcpval = gen(24)%e(nparam_per_node+17)                            !!>> AL 3-3-25

end subroutine suremuta_no_kadh_functional_net_rec_3

! Function to generate a normally distributed random number using the Box-Muller transform.
function gaussian_random(mu, sigma) result(r)
   implicit none
   real*8, intent(in) :: mu, sigma
   real*8 :: r, u1, u2, z

   call random_number(u1)
   call random_number(u2)

   ! Apply the Box-Muller transform:
   z = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
   r = mu + sigma * z
end function gaussian_random


subroutine duplication_rec(gendup)             !!>> AL 28-11-2020
!important: We are not using reactions right now
   implicit none
   integer, intent(in)  :: gendup
   integer              :: ihc,jhc,first_empty_gen
   real*8               :: dup_mu,dup_mich,dup_diffu,dup_kindof
   real*8, allocatable, dimension(:) :: dup_gen_info,dup_beh

   first_empty_gen = -1

   do ihc=1,ng
      if (int(gen(ihc)%kindof) .eq. 9) then 
         first_empty_gen = ihc
         exit
      end if
   end do 

   !if(first_empty_gen == -1) then      !!There are no positions left, in theory this cannot happen
   !   go to 757
   !end if

   !duplication
   gen(first_empty_gen)%diffu    = gen(gendup)%diffu 
   gen(first_empty_gen)%mu       = gen(gendup)%mu 
   gen(first_empty_gen)%mich     = gen(gendup)%mich 
   gen(first_empty_gen)%kindof   = gen(gendup)%kindof 
   gen(first_empty_gen)%label    = gen(gendup)%label
   gen(first_empty_gen)%e(:)     = gen(gendup)%e(:)

   gen(first_empty_gen)%t(:)     = gen(gendup)%t(:)                    !who regulates gendup
   gen(first_empty_gen)%t(first_empty_gen) = gen(gendup)%t(gendup)     !self regulation

   do ihc=1,ng
      gen(ihc)%t(first_empty_gen) = gen(ihc)%t(gendup)                 !who does gendup regulate
   end do 

   gex(:,first_empty_gen)=gex(:,gendup)                                !we also copy gene expression throughout the tissue

end subroutine duplication_rec

subroutine deletion_rec(gendel)
!important: We are not using reactions right now
   implicit none
   integer, intent(in)  :: gendel
   integer              :: ihc

   !deletion
   gen(gendel)%diffu    = 0.0d0
   gen(gendel)%mu       = 0.0d0
   gen(gendel)%mich     = 0.0d0
   gen(gendel)%kindof   = 9                !AL 3-12-24 ghost gene identifier
   gen(gendel)%e(:)     = 0.0d0
   gen(gendel)%t(:)     = 0.0d0                   !who regulates gendup

   do ihc=1,ng
      gen(ihc)%t(gendel) = 0.0d0              !who does gendup regulates
   end do 

   gex(:,gendel) = 0.0d0                      !we also deletion gene expression throughout the tissue
end subroutine deletion_rec 

! AL 1-7-25: this was slightly wrong. 
! subroutine functional_net(good) 
!    integer, allocatable, dimension(:) :: good   !!>>AL 10-4-24: this vector has elements = 1 when genes are functional
!    integer :: ich,jch,lch
   
!    if(allocated(good)) deallocate(good)         !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
!    allocate(good(1:ng)) 
!    !determine functional subnetwork
!    !search the genes that activate any cell property or behaviour
!    good=0
!    do ich=1,nga !nga = this is the number of genetically affectable node and cell parameters 
!       do jch=1,ng
!          if (gen(jch)%kindof .eq. 9) cycle
!          if (gen(jch)%e(ich) .eq. 0.0d0) cycle
!          good(jch)=1                           
!       enddo
!    enddo
   
!    !! Search any gene upstream a gene that activate cell properties/behaviours
!    lch=1
!    do while(lch>0)                              !! do until we cannot add any new gene
!       lch=0                                     !! newly added genes
!       do ich=1,ng                               !! go through the genes
!          if(int(gen(ich)%kindof) .eq. 9) cycle
!          if(good(ich).ne.1) cycle                !! That are upstream of cell properties and behaviours
!          do jch=1,ng                            !! Go over the rest of the genes
!             if(gen(jch)%kindof .eq. 9) cycle
!             if(good(jch)==1) cycle               !! test whether is already in our list
!             if (gen(ich)%t(jch) .ne. 0.0d0)then   !! genes ich and jhc are connected
!                good(jch)=1                      !! then this gene is upstream
!                lch=lch+1                        !! we added one gene to our list
!             endif
!          enddo
!       enddo
!    enddo
! end subroutine functional_net

subroutine functional_net(good)

integer, allocatable, dimension(:) :: up_beh,dw_const,good
integer :: ihc,jhc,lhc,khc         

if(allocated(up_beh)) deallocate(up_beh)                  !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
allocate(up_beh(1:ng))                           
if(allocated(dw_const)) deallocate(dw_const)                  !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
allocate(dw_const(1:ng))      
if(allocated(good)) deallocate(good)                  !!>> HC 14-2-2024 Here we save the genes that are upstream of a cell property/behaviour
allocate(good(1:ng))

good = 0
dw_const = 0
up_beh=0

!! Search any gene downstream of cosntitutivly expressed gene [1,2,23,24]
dw_const(1)  = 1
dw_const(2)  = 1
dw_const(23) = 1
dw_const(24) = 1

khc=1
do while(khc>0)               !! do until we cannot add any new gene
   khc=0                      !! newly added genes
   do ihc=1,ng                !! go through the genes
      if(dw_const(ihc) .ne. 1) cycle !! That are upstream of cell properties and behaviours
      do jhc=1,ng             !! Go over the rest of the genes
         if(dw_const(jhc) == 1)cycle  !! test whether is already in our list
         if(gen(jhc)%t(ihc) .gt. 0.0d0)then  !! genes ihc and jhc are connected
            dw_const(jhc) = 1                     !! then this gene is upstream
            khc=khc+1                       !! we added one gene to our list
         end if
      end do
   end do
end do

!! Search any gene upstream of cell properties/behaviours
!! search the genes that activate any cell property or behaviour
do ihc=1,nga
   do jhc=1,ng
      if (gen(jhc)%e(ihc)==0.0d0)cycle
      up_beh(jhc)=1                           
   enddo
enddo

khc=1
do while(khc>0)               !! do until we cannot add any new gene
   khc=0                      !! newly added genes
   do ihc=1,ng                !! go through the genes
      if(up_beh(ihc).ne.1)cycle !! That are upstream of cell properties and behaviours
      do jhc=1,ng             !! Go over the rest of the genes
         if(up_beh(jhc)==1)cycle  !! test whether is already in our list
         if(gen(ihc)%t(jhc) .ne. 0.0d0)then  !! genes ihc and jhc are connected
            up_beh(jhc)=1                     !! then this gene is upstream
            khc=khc+1                       !! we added one gene to our list
         endif
      end do
   end do
end do

!look for set intersection betwween up_beh and dw_const
do ihc = 1,ng
   if(up_beh(ihc) == 1 .and. dw_const(ihc) == 1)then
      good(ihc) = 1
   end if 
end do 

end subroutine functional_net

subroutine duplication(i)
  integer, intent(in) :: i
  !integer j,ii
  !if (allocated(father_son)) deallocate(father_son)
  !allocate(father_son(ng*ng)) !since each gene can no duplicate more than ng times per mutation round
  !print*,"father_son allocated"
  !ng_before_dup=ng
  !dupcount=0
!print *,i,"the gene being duplicated, and it has a npost of",gen(i)%npost,"that is",gen(i)%post
  call dupli(i)
  !father_son(i)=ng
  !print*,'Estamos ss'
  !if (allocated(dupi)) deallocate(dupi)
  !allocate(dupi(ng*ng))       !since each gene can no duplicate more than ng times per mutation round
  !print*,"dupi allocated"
  !dupi(i)=1
  !dupcount=dupcount+1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 
!!>> HC 16-9-2021 !! ATENTION !! We are not using pre/post forms right now, this should be changed if you want them!!      !!>> HC 16-9-2021
!  call follow_to_the_end(i,ng) !all the post forms of the gene also get duplicated and they get posts and pre within them !!>> HC 16-9-2021
!  do j=ng_before_dup+1,ng_before_dup+dupcount                                                                             !!>> HC 16-9-2021
!    if (gen(j)%npre>0) then                                                                                               !!>> HC 16-9-2021
!      do ii=1,gen(j)%npre                                                                                                 !!>> HC 16-9-2021
!        gen(j)%pre(ii)=father_son(gen(j)%pre(ii)) !ACHTUNG ERROR HERE                                                     !!>> HC 16-9-2021
!      end do                                                                                                              !!>> HC 16-9-2021
!    end if                                                                                                                !!>> HC 16-9-2021
!    if (gen(j)%npost>0) then                                                                                              !!>> HC 16-9-2021
!      do ii=1,gen(j)%npost                                                                                                !!>> HC 16-9-2021
!        gen(j)%post(ii)=father_son(gen(j)%post(ii))                                                                       !!>> HC 16-9-2021
!      end do                                                                                                              !!>> HC 16-9-2021
!    end if                                                                                                                !!>> HC 16-9-2021
!  end do                                                                                                                  !!>> HC 16-9-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021 

end subroutine

!**********************************************************************************************************
subroutine dupli(gendup)
!  type(genes) :: cgen(ng+1) !increase the size of the list of genes by one
  integer cpost(ng)
  integer cpre(ng)
  real*8  w(ng,ng) !stores the w array for all genes
  real*8  wa(ng,nga) !stores the wa array for all genes
  real*8  ckadh(ntipusadh,ntipusadh)
  real*8  gdi(ng) !stores the diffusion rates
  real*8  gmu(ng) !stores the genes' degradation rates
  real*8  gmich(ng) !!>> HC 7-7-2020 Stores the constant of michaelis-menten
  character*300 glab(ng) !stores the genes' labels
  real*8  cgex(nda,ng),cww(ng,ng,3)
  integer gkindof(ng) !stores the genes' types (functions)
  integer mpost(ng,ng)
  integer mpre(ng,ng)
  integer, intent(in) ::  gendup
  integer i,j, check
  integer genadhori, genadhnew
  character*300 templabel
   
  do i=1,ng              !we deallocate the matrices within the type gene
    gdi(i)=gen(i)%diffu  
    gmu(i)=gen(i)%mu
    gmich(i)=gen(i)%mich  !!>> HC 7-7-2020
    gkindof(i)=int(gen(i)%kindof) !added explicit cast. Kindof should have been an int. >>Renske oct 18
    glab(i)=trim(gen(i)%label)
    !print *, "copied label is "//trim(glab(i))
    w(i,:)=gen(i)%t(:)   !!>> HC 30-6-2020
    wa(i,:)=gen(i)%e(:)  !!>> HC 30-6-2020
    deallocate(gen(i)%t) !!>> HC 30-6-2020
    deallocate(gen(i)%e) !!>> HC 30-6-2020
    cpre(i)=gen(i)%npre
    if (gen(i)%npre>0) then
      mpre(i,:gen(i)%npre)=gen(i)%pre(:gen(i)%npre)  !we save the values in pre and post since we dont know their sizes
      deallocate(gen(i)%pre)
    end if
    cpost(i)=gen(i)%npost
    if (gen(i)%npost>0) then
      mpost(i,:gen(i)%npost)=gen(i)%post(:gen(i)%npost)  !we save the values in pre and post since we dont know their sizes
      deallocate(gen(i)%post)
    end if
    cww(i,:ng,:)=gen(i)%r(:ng,:)  !!>> HC 30-6-2020
    deallocate(gen(i)%r)          !!>> HC 30-6-2020
  end do
  
  deallocate(gen) !deallocate the entire array, now that we stored all info elsewhere, and re-allocate a bigger array
  ng=ng+1
  allocate(gen(ng))
  do i=1,ng-1
    gen(i)%nww=0 !!>> HC 16-9-2021 !! ATENTION !! We are not using reactions right now, this should be changed if you want them!!
    gen(i)%diffu=gdi(i)
    gen(i)%mu=gmu(i)
    gen(i)%mich=gmich(i)
    gen(i)%kindof=gkindof(i)
!    gen(i)%label=trim(glab(i))
    !print *, "recopied label is "//trim(gen(i)%label)
    allocate(gen(i)%t(ng))     !!>> HC 30-6-2020
    gen(i)%t(:ng-1)=w(i,:ng-1) !!>> HC 30-6-2020
    allocate(gen(i)%e(nga))    !!>> HC 30-6-2020
    gen(i)%e(:)=wa(i,:)        !!>> HC 30-6-2020
    gen(i)%npre=0              !!>> HC 11-11-2021 !! ATENTION  we are not using pre or post forms right now, so 
    gen(i)%npost=0             !!>> HC 11-11-2021 !! ATENTION  this is set to 0 to be sure they are not activated
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021  ATENTION!!    
!    gen(i)%npre=cpre(i)                                     !!>> HC 11-11-2021  PRE and POST forms and REACTIONS
!    if (gen(i)%npre>0) then                                 !!>> HC 11-11-2021  not in use right now, uncomment to activate them
!      allocate(gen(i)%pre(gen(i)%npre))                     !!>> HC 11-11-2021
!      gen(i)%pre(:gen(i)%npre)=mpre(i,:gen(i)%npre)         !!>> HC 11-11-2021
!    end if                                                  !!>> HC 11-11-2021
!    gen(i)%npost=cpost(i)                                   !!>> HC 11-11-2021
!    if (gen(i)%npost>0) then                                !!>> HC 11-11-2021
!      allocate(gen(i)%post(gen(i)%npost))                   !!>> HC 11-11-2021
!      gen(i)%post(:gen(i)%npost)=mpost(i,:gen(i)%npost)     !!>> HC 11-11-2021
!    end if                                                  !!>> HC 11-11-2021
!    allocate(gen(i)%r(ng,3))         !!>> HC 30-6-2020      !!>> HC 11-11-2021
    !gen(i)%r=cww(i,:,:)             !!>> HC 30-6-2020       !!>> HC 11-11-2021
!    gen(i)%r(:ng-1,:)=cww(i,:ng-1,:) !!>> HC 30-6-2020      !!>> HC 11-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021
  end do 
  
  !now for the new gene
  gen(ng)%diffu=gdi(gendup)
  gen(ng)%mu=gmu(gendup)
  gen(ng)%mich=gmich(gendup)
  gen(ng)%kindof=gkindof(gendup)
  !now we make the label. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021  ATENTION !!
!  write(templabel, '(I3)') ng                                             !!>> HC 11-11-2021 Labels not in use right now, uncomment to activate
!  !Check if the label is not going to be too big                          !!>> HC 11-11-2021
!  if (len(gen(gendup)%label)-len_trim(gen(gendup)%label)>4) then          !!>> HC 11-11-2021
!    gen(ng)%label=trim(gen(gendup)%label)//","//trim(templabel)           !!>> HC 11-11-2021
!  else                                                                    !!>> HC 11-11-2021
!    do i = 1, len_trim(gen(gendup)%label)                                 !!>> HC 11-11-2021
!      if (gen(gendup)%label(i:i).eq.",") then                             !!>> HC 11-11-2021
!        check = i+1                                                       !!>> HC 11-11-2021
!        exit                                                              !!>> HC 11-11-2021
!      end if                                                              !!>> HC 11-11-2021
!    end do                                                                !!>> HC 11-11-2021
!    gen(ng)%label=trim(gen(gendup)%label(check:))//","//trim(templabel)   !!>> HC 11-11-2021
!  end if                                                                  !!>> HC 11-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021 
  
  allocate(gen(ng)%t(ng))          !!>> HC 30-6-2020
  gen(ng)%t(:ng-1)=w(gendup,:ng-1) !!>> HC 30-6-2020
  gen(ng)%t(ng)=w(gendup,gendup) !the regulation of this gene by itself is the same as that of the parent gene towards itself !!>> HC 30-6-2020

    gen(ng)%npre=0
    gen(ng)%npost=0    
!  allocate(gen(ng)%r(ng,3))                 !!>> HC 30-6-2020  !!>> HC 10-11-2021
!  gen(ng)%r(:ng-1,:)=gen(gendup)%r(:ng-1,:) !!>> HC 30-6-2020  !!>> HC 10-11-2021
!  gen(ng)%r=gen(gendup)%r                   !!>> HC 30-6-2020  !!>> HC 10-11-2021
  gen(ng)%nww=gen(gendup)%nww
  allocate(gen(ng)%e(nga))  !!>> HC 30-6-2020
  gen(ng)%e(:)=wa(gendup,:) !!>> HC 30-6-2020
  gen(ng)%npre=0            !!>> HC 11-11-2021 NPRE and NPOST not used right now
  gen(ng)%npost=0           !!>> HC 11-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021  ATENTION!!
!  gen(ng)%npre=cpre(gendup)                                 !!>> HC 11-11-2021  NPRE and NPOST not used right now
!  if (gen(ng)%npre>0) then                                  !!>> HC 11-11-2021  not in use right now, uncomment to activate them
!    allocate(gen(ng)%pre(gen(gendup)%npre))                 !!>> HC 11-11-2021
!    gen(ng)%pre(:gen(ng)%npre)=mpre(gendup,:gen(ng)%npre)   !!>> HC 11-11-2021
!  end if                                                    !!>> HC 11-11-2021
!  gen(ng)%npost=cpost(gendup)                               !!>> HC 11-11-2021
!  if (gen(ng)%npost>0) then                                 !!>> HC 11-11-2021
!    allocate(gen(ng)%post(gen(ng)%npost))                   !!>> HC 11-11-2021
!    gen(ng)%post(:gen(ng)%npost)=mpost(gendup,:gen(ng)%npost) !!>> HC 11-11-2021
!  end if                                                    !!>> HC 11-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 11-11-2021  ATENTION!!
  !we have to re-index the w 
  do i=1,ng-1
    gen(i)%t(ng)=gen(i)%t(gendup)     !!>> HC 30-6-2020
  end do

  !we have to add the self reference
  gen(ng)%t(ng)=gen(gendup)%t(gendup) !!>> HC 30-6-2020

  !we have to re-make the ww
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021
!!>> HC 16-9-2021 WE ARE NOT USING PRE_POS AT THE MOMENT        !!>> HC 16-9-2021
  !we have to re-index the pre and post matrices                !!>> HC 16-9-2021
  !post                                                         !!>> HC 16-9-2021
!  do i=1,ng-1                                                  !!>> HC 16-9-2021
!    if (gen(i)%npost/=0) then                                  !!>> HC 16-9-2021
!      do j=1,gen(i)%npost                                      !!>> HC 16-9-2021
!        if (gen(i)%post(j)==gendup) then                       !!>> HC 16-9-2021
!          cpost=0                                              !!>> HC 16-9-2021
!          cpost(1:gen(i)%npost)=gen(i)%post(gen(i)%npost)      !!>> HC 16-9-2021
!          deallocate(gen(i)%post)                              !!>> HC 16-9-2021
!          gen(i)%npost=gen(i)%npost+1                          !!>> HC 16-9-2021
!          allocate(gen(i)%post(gen(i)%npost))                  !!>> HC 16-9-2021
!          gen(i)%post(1:gen(i)%npost)=cpost(1:gen(i)%npost)    !!>> HC 16-9-2021
!          gen(i)%post(gen(i)%npost)=ng                         !!>> HC 16-9-2021
!        end if                                                 !!>> HC 16-9-2021
!      end do                                                   !!>> HC 16-9-2021
!    end if                                                     !!>> HC 16-9-2021
!  end do                                                       !!>> HC 16-9-2021

  !pre                                                          !!>> HC 16-9-2021
!  do i=1,ng-1                                                  !!>> HC 16-9-2021
!    if (gen(i)%npre/=0) then                                   !!>> HC 16-9-2021
!      do j=1,gen(i)%npre                                       !!>> HC 16-9-2021
!        if (gen(i)%pre(j)==gendup) then                        !!>> HC 16-9-2021
!          cpre=0                                               !!>> HC 16-9-2021
!          cpre(1:gen(i)%npre)=gen(i)%pre(gen(i)%npre)          !!>> HC 16-9-2021
!          deallocate(gen(i)%pre)                               !!>> HC 16-9-2021
!          gen(i)%npre=gen(i)%npre+1                            !!>> HC 16-9-2021
!          allocate(gen(i)%pre(gen(i)%npre))                    !!>> HC 16-9-2021
!          gen(i)%pre(1:gen(i)%npre)=cpre(1:gen(i)%npre)        !!>> HC 16-9-2021
!          gen(i)%pre(gen(i)%npre)=ng                           !!>> HC 16-9-2021
!        end if                                                 !!>> HC 16-9-2021
!      end do                                                   !!>> HC 16-9-2021
!    end if                                                     !!>> HC 16-9-2021
!  end do                                                       !!>> HC 16-9-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 16-9-2021
  ! now if the form is an adhesion molecule we have to re-size kadh and renumber wa(1)
  if (gen(gendup)%e(1)/=0) then  !!>> HC 30-6-2020
    ntipusadh=ntipusadh+1
    gen(ng)%e(1)=ntipusadh       !!>> HC 30-6-2020
    ckadh(:ntipusadh-1,:ntipusadh-1)=kadh(:ntipusadh-1,:ntipusadh-1)
    deallocate(kadh)
    allocate(kadh(ntipusadh,ntipusadh))
    kadh=0.0d0
    kadh(:ntipusadh-1,:ntipusadh-1)=ckadh(:ntipusadh-1,:ntipusadh-1)
    genadhnew=int(gen(ng)%e(1))     !!>> HC 30-6-2020
    genadhori=int(gen(gendup)%e(1)) !!>> HC 30-6-2020
    kadh(genadhnew,:)=kadh(genadhori,:)
    kadh(genadhnew,genadhnew)=kadh(genadhori,genadhori)
  end if
  cgex=gex
  deallocate(gex)
  allocate(gex(nda,ng))
  gex(:,:ng-1)=cgex
  gex(:,ng)=gex(:,gendup)
 

end subroutine

!*****************************************************************

subroutine add_form(gendup)
  integer cpost(ng)
  integer gendup
  integer i,j,jj

  !we have to delete any connection that ng may get or give to any gene
  do i=1,ng
    do j=1,gen(i)%npre
      if (gen(i)%pre(j)==ng) then !this has to be deleted
        do jj=j,gen(i)%npre-1
          gen(i)%pre(jj)=gen(i)%pre(jj-1)
        end do
        gen(i)%npre=gen(i)%npre-1
        exit
      end if
    end do
    do j=1,gen(i)%npost
      if (gen(i)%post(j)==ng) then !this has to be deleted
        do jj=j,gen(i)%npost-1
          gen(i)%post(jj)=gen(i)%post(jj-1)
        end do
        gen(i)%npost=gen(i)%npost-1
        exit
      end if
    end do
  end do

  !post; so the gendup gets a new form, thus a new post
  cpost=0
  if (gen(gendup)%npost>0) then
    cpost(:gen(gendup)%npost)=gen(gendup)%post(:gen(gendup)%npost)
    deallocate(gen(gendup)%post)
  end if
  gen(gendup)%npost=gen(gendup)%npost+1
  allocate(gen(gendup)%post(gen(gendup)%npost))
  gen(gendup)%post(:gen(gendup)%npost)=cpost(:gen(gendup)%npost)
  gen(gendup)%post(gen(gendup)%npost)=ng

  ! the kindof of gendup may change
  gen(ng)%kindof=gen(gendup)%kindof       !we copy that
  if (gen(ng)%kindof==2) gen(ng)%kindof=3 !because now it has a pre too
  if (gen(ng)%kindof==1) gen(ng)%kindof=3

  !pre the only pre of the new gene is gendup
  gen(ng)%npre=1
  if (allocated(gen(ng)%pre)) deallocate(gen(ng)%pre)
  allocate(gen(ng)%pre(1))
  gen(ng)%pre(1)=gendup
  gen(ng)%npost=0

end subroutine

!****************************************************************

recursive subroutine follow_to_the_end(father,son)
  integer, intent(in) :: father,son
  integer last,sono
  integer i,j,k,ii,kk,jp


  if (gen(father)%npost>0) then
!    if (allocated(gen(son)%pre)) deallocate(gen(son)%pre)
!    if (allocated(gen(son)%post)) deallocate(gen(son)%post)

!    gen(son)%npre=gen(father)%npre
!    gen(son)%npost=gen(father)%npost

!    allocate(gen(son)%pre(gen(son)%npre))
!    allocate(gen(son)%post(gen(son)%npost))

    sono=son

    do jp=1,gen(father)%npost
      last=gen(father)%post(jp)
      if (dupi(last)==0) then
        call dupli(last)
        father_son(last)=ng
        dupcount=dupcount+1
        dupi(last)=1
      end if
!      gen(sono)%post(jp)=ng
!      if (allocated(gen(ng)%pre)) deallocate(gen(ng)%pre)
!      if (allocated(gen(ng)%post)) deallocate(gen(ng)%post)
!      gen(ng)%npre=gen(last)%npre 
!      gen(ng)%npost=gen(last)%npost
!      allocate(gen(ng)%pre(gen(ng)%npre))
!      allocate(gen(ng)%post(gen(ng)%npost))
!      gen(ng)%pre(1)=sono
      call follow_to_the_end(last,ng)
    end do  
  end if
end subroutine

!*******************************************************************

subroutine del(delgen)

  !note that all the forms derived from a gene would also disappear
  type(genes) :: cgen(ng-1)
  integer, intent(in) :: delgen
  integer pres(ng),posts(ng)
  real*8 cgex(nda,ng-1)
  real*8 cw(ng-1),ckadh(ntipusadh,ntipusadh),cww(ng-1,3),ccww(ng*ng,3)
  integer i,j,moldadh

  pres=0 ; posts=0

  !for the contents of WW
  do i=1,ng
    k=gen(i)%nww
    
    ccww(:ng,:)=gen(i)%r(:ng,:)   !!>> HC 30-6-2020

    kk=0
    do j=1,k
      if (gen(i)%r(j,1)==delgen.or.gen(i)%r(j,2)==delgen) then !!>> HC 30-6-2020
        kk=kk+1
        do jj=j,k-1
          ccww(jj,:)=ccww(jj+1,:)
        end do
        ccww(k,:)=0.0d0
      else
        if (gen(i)%r(j,1)>delgen) then !!>> HC 30-6-2020
          ccww(j,1)=ccww(j,1)-1
        end if
        if (gen(i)%r(j,2)>delgen) then !!>> HC 30-6-2020
          ccww(j,2)=ccww(j,2)-1
        end if
      end if
    end do
    gen(i)%nww=gen(i)%nww-kk
    gen(i)%r(:ng,:)=ccww(:ng,:)  !!>> HC 30-6-2020
  end do

  ! if delgen is an adhesion molecule we have to re-size kadh
  ! now if the form is an adhesion molecule we have to re-size kadh and renumber wa(1)
  if (gen(delgen)%e(1)/=0) then   !!>> HC 30-6-2020
    moldadh=int(gen(delgen)%e(1)) !Added explicit cast >> Renske oct 18 !!>> HC 30-6-2020
    do j=1,ng
      if (gen(j)%e(1)>moldadh) gen(j)%e(1)=gen(j)%e(1)-1 !!>> HC 30-6-2020
    end do    
    ntipusadh=ntipusadh-1
    if (ntipusadh>0) then
      ckadh(:ntipusadh+1,:ntipusadh+1)=kadh(:ntipusadh+1,:ntipusadh+1)
      deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0.0d0
      kadh(:moldadh-1,:ntipusadh)=ckadh(:moldadh-1,:ntipusadh)
      kadh(moldadh:ntipusadh,:ntipusadh)=ckadh(moldadh+1:ntipusadh+1,:ntipusadh)
      kadh(:ntipusadh,:moldadh-1)=ckadh(:ntipusadh,:moldadh-1)
      kadh(:ntipusadh,moldadh:ntipusadh)=ckadh(:ntipusadh,moldadh+1:ntipusadh+1)
    else
      deallocate(kadh)
    end if
  end if

  ! we take away from pre and post the references to the deleted gene
  do i=1,ng-1
    do j=1,gen(i)%npre
      if (gen(i)%pre(j)==delgen) then
        pres(:j-1)=gen(i)%pre(:j-1)
        pres(j:gen(i)%npre-1)=gen(i)%pre(j+1:gen(i)%npre)
        deallocate(gen(i)%pre)
        gen(i)%npre=gen(i)%npre-1
        allocate(gen(i)%pre(gen(i)%npre))
        gen(i)%pre(:gen(i)%npre)=pres(:gen(i)%npre)
        exit
      end if
    end do
    do j=1,gen(i)%npost
      if (gen(i)%post(j)==delgen) then
        posts(:j-1)=gen(i)%post(:j-1)
        posts(j:gen(i)%npost-1)=gen(i)%post(j+1:gen(i)%npost)
        deallocate(gen(i)%post)
        gen(i)%npost=gen(i)%npost-1
        allocate(gen(i)%post(gen(i)%npost))
        gen(i)%post(:gen(i)%npost)=posts(:gen(i)%npost)
        exit
      end if
    end do
  end do

  ! the indices in pre and post need to be -1 for the genes after delgen
  do i=1,ng
    do j=1,gen(i)%npre
      if (gen(i)%pre(j)>=delgen) gen(i)%pre(j)=gen(i)%pre(j)-1
    end do
    do j=1,gen(i)%npost
      if (gen(i)%post(j)>=delgen) gen(i)%post(j)=gen(i)%post(j)-1
    end do
  end do

  ! re-make W
!  do i=1,ng-1
!    cw(:delgen-1)=gen(i)%t(:delgen-1)     !!>> HC 30-6-2020
!    cw(delgen:ng-1)=gen(i)%t(delgen+1:ng) !!>> HC 30-6-2020
!    deallocate(gen(i)%t)                  !!>> HC 30-6-2020
!    allocate(gen(i)%t(ng-1))              !!>> HC 30-6-2020
!  end do
  cgex(:,:delgen-1)=gex(:,:delgen-1)
  cgex(:,delgen:ng-1)=gex(:,delgen+1:ng)

  deallocate(gex)
  allocate(gex(nda,ng-1))

  gex=cgex

  !do the actual deletion and modify the w and ww matrices of the other genes to exclude the deleted gene  
  do i=1,delgen-1
    cgen(i)=gen(i)

    deallocate(cgen(i)%t)                     !!>> HC 30-6-2020
    allocate(cgen(i)%t(ng-1))                 !!>> HC 30-6-2020
    cgen(i)%t(:delgen-1)=gen(i)%t(:delgen-1)  !!>> HC 30-6-2020
    if (delgen<ng) then
      cgen(i)%t(delgen:ng-1)=gen(i)%t(delgen+1:ng) !!>> HC 30-6-2020
    end if

    deallocate(cgen(i)%r)                !!>> HC 30-6-2020
    allocate(cgen(i)%r(ng-1,3))          !!>> HC 30-6-2020
    cgen(i)%r(:ng-1,:)=gen(i)%r(:ng-1,:) !!>> HC 30-6-2020

    deallocate(cgen(i)%e)    !!>> HC 30-6-2020
    allocate(cgen(i)%e(nga)) !!>> HC 30-6-2020
    cgen(i)%e=gen(i)%e       !!>> HC 30-6-2020

    if (gen(i)%npre>0) then
      if (allocated(cgen(i)%pre)) deallocate(cgen(i)%pre)
      allocate(cgen(i)%pre(gen(i)%npre))
      cgen(i)%pre=gen(i)%pre
    end if
    if (gen(i)%npost>0) then
      if (allocated(cgen(i)%post)) deallocate(cgen(i)%post)
      allocate(cgen(i)%post(gen(i)%npost))
      cgen(i)%post=gen(i)%post
    end if
  end do

  do i=delgen,ng-1

    cgen(i)=gen(i+1)
    deallocate(cgen(i)%t)                            !!>> HC 30-6-2020
    allocate(cgen(i)%t(ng-1))                        !!>> HC 30-6-2020
    cgen(i)%t(:delgen-1)=gen(i+1)%t(:delgen-1)       !!>> HC 30-6-2020
    if (delgen<ng) then
      cgen(i)%t(delgen:ng-1)=gen(i+1)%t(delgen+1:ng) !!>> HC 30-6-2020
    end if

    deallocate(cgen(i)%r)                   !!>> HC 30-6-2020
    allocate(cgen(i)%r(ng-1,3))             !!>> HC 30-6-2020
    cgen(i)%r(:ng-1,:)=gen(i+1)%r(:ng-1,:)  !!>> HC 30-6-2020

    deallocate(cgen(i)%e)    !!>> HC 30-6-2020
    allocate(cgen(i)%e(nga)) !!>> HC 30-6-2020
    cgen(i)%e=gen(i+1)%e     !!>> HC 30-6-2020

    if (gen(i+1)%npre>0) then
      if (allocated(cgen(i)%pre)) deallocate(cgen(i)%pre)
      allocate(cgen(i)%pre(gen(i+1)%npre))
      cgen(i)%pre=gen(i+1)%pre
    end if
    if (gen(i+1)%npost>0) then
      if (allocated(cgen(i)%post)) deallocate(cgen(i)%post)
      allocate(cgen(i)%post(gen(i+1)%npost))
      cgen(i)%post=gen(i+1)%post
    end if
  end do
  deallocate(gen)
  allocate(gen(ng-1))

  do i=1,ng-1
    gen(i)=cgen(i)
    if (allocated(gen(i)%t)) deallocate(gen(i)%t) !!>> HC 30-6-2020
    allocate(gen(i)%t(ng-1))                      !!>> HC 30-6-2020
    gen(i)%t=cgen(i)%t                            !!>> HC 30-6-2020
    if (allocated(gen(i)%r)) deallocate(gen(i)%r) !!>> HC 30-6-2020
    allocate(gen(i)%r(ng-1,3))                    !!>> HC 30-6-2020
    gen(i)%r=cgen(i)%r                            !!>> HC 30-6-2020
    if (allocated(gen(i)%e)) deallocate(gen(i)%e) !!>> HC 30-6-2020
    allocate(gen(i)%e(nga))                       !!>> HC 30-6-2020
    gen(i)%e=cgen(i)%e                            !!>> HC 30-6-2020

    if (gen(i)%npre>0) then
      if (allocated(gen(i)%pre)) deallocate(gen(i)%pre)
      allocate(gen(i)%pre(cgen(i)%npre))
      gen(i)%pre=cgen(i)%pre
    end if
    if (gen(i)%npost>0) then
      if (allocated(gen(i)%post)) deallocate(gen(i)%post)
      allocate(gen(i)%post(cgen(i)%npost))
      gen(i)%post=cgen(i)%post
    end if
  end do

  ng=ng-1

end subroutine

!********************************************************************************
subroutine deletion(i)
  integer, intent(in) :: i
  integer j,ii,jj,jjj
  integer last

  last=i
  todel=0
  todel(i)=1
  call follow_to_the_end_del(last)
  ang=ng    
  do j=ng,1,-1
    if (todel(j)==1) then ; if (ng>1) then ; call del(j) ; else ; return ;  end if ; end if
  end do
  do j=1,gen(ng)%npre
    if (gen(ng)%pre(j)>=ng) then
      do jj=1,ng
        do jjj=1,gen(jj)%npost
          if (gen(jj)%post(jjj)==ng) gen(ng)%pre(j)=jj
        end do
      end do
    end if
  end do
end subroutine

!****************************************************************

recursive subroutine follow_to_the_end_del(father)
  integer, intent(in) :: father
  integer last
  integer i,j,k,ii,kk,jp,jpp

  if (gen(father)%npost>0) then
ro: do jp=1,gen(father)%npost
      last=gen(father)%post(jp)
      todel(last)=1
      if (gen(father)%npre>0) then
        do jpp=1,gen(father)%npre
          if (last==gen(father)%pre(jpp)) exit ro
        end do
      end if
      call follow_to_the_end_del(last)
    end do ro 
  end if
end subroutine

!*******************************************************************

subroutine write_evaparam(pop_size,s,faname)
  integer pop_size
  character(len=*), intent(in) :: faname
  real*8  s

  open(11,file=trim(faname)//"params_eva")

  write (11,*) "this file contains all the evolution parameters used in the simulation from the:"
  write (11,*) faname,"this is the founder father"

  write (11,*) pop_size,"populational size"
  write (11,*) s,"selection coefficient: simply multiplies the fitness of each individual" 
  write (11,*) ismurate,"IS mutation rate"
  write (11,*) mag_ismurate,"magnitude of the IS mutation rate in % of the original value"
  write (11,*) tmuratew_g,"T mutation that change w: gain new interaction"
  write (11,*) tmuratew_l,"T mutation that change w: lose existing interaction"
  write (11,*) tmurateww,"T mutation that changes ww"
  write (11,*) tmuratewa,"T mutation that change wa"
  write (11,*) tmurateadh,"T mutation that change kadh"
  write (11,*) mag_tmurate,"magnitude of the T mutation rate in % of the original value"
  write (11,*) duprate,"gene duplication rate (also of deletion)"
  write (11,*) new_form_mu_rate,"mutation rate to get a new form, or to lose an existing form"
print*,"muta1152"
!close(11)
end subroutine

!***********************************************************************

subroutine read_evaparam(pop_size,s,maxrtf, faname)     
  integer pop_size, s
  character(len=*), intent(in) :: faname
  real*8  maxrtf

  print *,len(faname)

  open(11,file=faname,iostat=i)
  print *,"this is zero if the file was opened successfully",i
  read (11,*) 
  read (11,*) pop_size  !,ERR=373,END=212 The error messages are triggered by default, even if reading went fine
  print *,pop_size,"psize"
  read (11, *) maxrtf
  print *, maxrtf, "max real time duration of development (if run with founding father, should match its first param)"
  read (11,*) s
  print *,s,"the random seed" 
  read (11,*) ismurate
  print *,ismurate,"IS mutation rate"
  read (11,*) mag_ismurate
  print *,mag_ismurate,"magnitude of the IS mutation rate in % of the original value"
  read (11,*) tmuratew_g
  print *,tmuratew_g,"T mutation that change w: gain new interaction" 
  read (11,*) tmuratew_l
  print *,tmuratew_l,"T mutation that change w: lose existing interaction"
  read (11,*) tmurateww
  print *,tmurateww,"T mutation that changes ww"
  read (11,*) tmuratewa
  print *,tmuratewa,"T mutation that change wa"
  read (11,*) tmurateadh
  print *,tmurateadh,"T mutation that change kadh"
  read (11,*) mag_tmurate
  print *,mag_tmurate,"magnitude of the T mutation rate in % of the original value"
  read (11,*) duprate
  print *,duprate,"gene duplication rate (also of deletion)"
  read (11,*) new_form_mu_rate
  print *,new_form_mu_rate,"mutation rate to get a new form, or to lose an existing form"
  read (11,*) distfitmag
  print *,distfitmag,"scales the height of the hill function for fitness due to morph distance"
  read (11,*) distfitscale
  print *,distfitscale,"scales the slope of the hill function for fitness due to morph distance"
close(11)
!373 print *,"error, but can safely be ignored"
!212 print *,"end of file, but not really"

end subroutine

subroutine read_elli_evaparam(faname, effu)     
  character(len=*), intent(in) :: faname
  real*8  maxrtf
  integer :: ich, neffu
  integer, allocatable, dimension(:) :: effu

  !print *,len(faname)

  open(11,file=faname,iostat=i)
  !print *,"this is zero if the file was opened successfully",i
  read (11,*) 
  read (11,*) ismurate
   !print *,ismurate,"IS mutation rate"
  read (11,*) mag_ismurate
   !print *,mag_ismurate,"magnitude of the IS mutation rate in % of the original value"
  read (11,*) tmuratew_g
   !print *,tmuratew_g,"T mutation that change w: gain new interaction" 
  read (11,*) tmuratew_l
   !print *,tmuratew_l,"T mutation that change w: lose existing interaction"
  read (11,*) tmurateww
   !print *,tmurateww,"T mutation that changes ww"
  read (11,*) tmuratewa
   !print *,tmuratewa,"T mutation that change wa"
  read (11,*) tmurateadh
   !print *,tmurateadh,"T mutation that change kadh"
  read (11,*) mag_tmurate
   !print *,mag_tmurate,"magnitude of the T mutation rate in % of the original value"
  read (11,*) duprate
   !print *,duprate,"gene duplication rate (also of deletion)"
  read (11,*) new_form_mu_rate
   !print *,new_form_mu_rate,"mutation rate to get a new form, or to lose an existing form"
  read (11,*) distfitmag
   !print *,distfitmag,"scales the height of the hill function for fitness due to morph distance"
  read (11,*) distfitscale
  read(11, *) tarfitmorphfile
  read(11, *) maxnodenr
  read(11, *)
  read(11, *) neffu
  if(allocated(effu)) deallocate(effu)
  allocate(effu(1:neffu))
  effu=0
  do ich=1, neffu
     read(11, *) effu(ich)
  enddo
  
   !print *,distfitscale,"scales the slope of the hill function for fitness due to morph distance"
close(11)
!373 print *,"error, but can safely be ignored"
!212 print *,"end of file, but not really"

end subroutine

!***************************************************
subroutine limits(limit2) !pfh***11-09-15

integer ii,jj,limit2
!real*8 :: aa

!print*,"limits"
limit2=1  !if limit =1 no limit has been surpased, else, =0
if (ntipusadh>0)then
if((abs(maxval(kadh))>b_max).or.(abs(minval(kadh))<b_min))then
  limit2=0;print*,"kadh",maxval(kadh);return;endif;endif
do ii=1,ng;!print*,"ii",ii
 do jj=1,ng
  if(abs(gen(ii)%t(jj))>w_max)then;limit2=0;print*,"larger w than allowed: ",ii,jj,gen(ii)%t(jj);return;endif !1, transcription !!>> HC 30-6-2020
  if(abs(gen(ii)%t(jj))<w_min)then;limit2=0;print*,"smaller w than allowed: ",ii,jj,gen(ii)%t(jj);return;endif  !!>> HC 30-6-2020
 enddo

  if((abs(gen(ii)%diffu)>d_max).or.abs(gen(ii)%diffu)<d_min)then;limit2=0;print*,"gen(ii)%diffu",gen(ii)%diffu;return;endif !2 diffusion

  if((abs(gen(ii)%mich)>mich_max).or.abs(gen(ii)%mich)<mich_min)then;limit2=0;print*,"gen(ii)%mich",gen(ii)%mich;return;endif !2 diffusion


  if((abs(gen(ii)%mu)>m_max).or.abs(gen(ii)%mu)<m_min)then
    if(ii/=5)then;limit2=0;print*,"gen(ii)%mu",gen(ii)%mu;return;endif;endif !3 decay
  if((abs(gen(ii)%e(nparam_per_node+14))>m_max).or.abs(gen(ii)%e(nparam_per_node+14))<m_min)then  !ecm decay  !!>> HC 30-6-2020
  limit2=0;print*,"gen(ii)%e(nparam_per_node+14",gen(ii)  %e(nparam_per_node+14);return;endif	!!>> HC 30-6-2020							 
  
  if((abs(maxval(gen(ii)%e(8:15)))>b_max).or.(abs(minval(gen(ii)%e(8:15)))<b_min))then !4 adhesion related and others !!>> HC 30-6-2020
  limit2=0;print*,"(maxval(gen(ii)%e(8:15))",maxval(gen(ii)%e(8:15));return;endif  !!>> HC 30-6-2020
 
  if((abs(maxval(gen(ii)%e(26:28)))>b_max).or.(abs(minval(gen(ii)%e(26:28)))<b_min))then !!>> HC 30-6-2020
  limit2=0;print*,"wa 26:28 over limit: ",maxval(gen(ii)%e(26:28));return;endif !!>> HC 30-6-2020

  if((abs(gen(ii)%e(nparam_per_node+9))>b_max).or.(abs(gen(ii)%e(nparam_per_node+9))<b_min))then !!>> HC 30-6-2020
  limit2=0;print*,"+9",gen(ii)%e(nparam_per_node+9);return;endif !!>> HC 30-6-2020
 
  if((abs(gen(ii)%e(nparam_per_node+8))>b_max).or.(abs(gen(ii)%e(nparam_per_node+8)))<b_min)then !!>> HC 30-6-2020
  limit2=0;print*,"+8",gen(ii)%e(nparam_per_node+8);return;endif !!>> HC 30-6-2020

  if((abs(maxval(gen(ii)%e(6:7)))>a_max).or.(abs(minval(gen(ii)%e(6:7)))<a_min))then !5 da, you !!>> HC 30-6-2020
  limit2=0;print*,"wa 6:7 over limit: ",maxval(gen(ii)%e(6:7));return;endif  !!>> HC 30-6-2020
  
  if((abs(maxval(gen(ii)%e(21:24)))>e_max).or.(abs(minval(gen(ii)%e(21:24)))<e_min))then !!>> HC 30-6-2020
  limit2=0;print*,"wa 21:24 over limit: ",maxval(gen(ii)%e(21:24));return;endif !6 Req wa's !!>> HC 30-6-2020

  if((abs(gen(ii)%e(nparam_per_node+3))>ap_max).or.(abs(gen(ii)%e(nparam_per_node+3))<ap_min))then !7 Apoptosis !!>> HC 30-6-2020
  limit2=0;print*,"+3",gen(ii)%e(nparam_per_node+3);return;endif !!>> HC 30-6-2020

  if((abs(gen(ii)%e(16))>fp_max).or.(abs(gen(ii)%e(16))<fp_min))then !8 dmo !!>> HC 30-6-2020
  limit2=0;print*,"wa 16 over limit: ", gen(ii)%e(16);return;endif !!>> HC 30-6-2020

  if((abs(gen(ii)%e(nparam_per_node+13))>em_max).or.(abs(gen(ii)%e(nparam_per_node+13))<em_min))then !9 EMT !!>> HC 30-6-2020
  limit2=0;print*,"+13",nparam_per_node+13;return;endif 
 
  if((abs(gen(ii)%e(nparam_per_node+4))>ma_max).or.(abs(gen(ii)%e(nparam_per_node+4))<ma_min))then !10 ECM segretion !!>> HC 30-6-2020
  limit2=0;print*,"+4",gen(ii)%e(nparam_per_node+4);return;endif !!>> HC 30-6-2020

  if((abs(gen(ii)%e(nparam_per_node+2))>dv_max).or.(abs(gen(ii)%e(nparam_per_node+2)))<dv_min)then !!>> HC 30-6-2020
  limit2=0;print*,"+2",gen(ii)%e(nparam_per_node+2);return;endif !11 cell division !!>> HC 30-6-2020
 
  if((abs(gen(ii)%e(nparam_per_node+16))>1000).or.(abs(gen(ii)%e(nparam_per_node+16)))<0)then !!>> HC 30-6-2020
  limit2=0;print*,"+16",gen(ii)%e(nparam_per_node+16);return;endif ! noise is biased by polarization vector, arbitrary limits !!>> HC 30-6-2020

enddo

end subroutine limits

end module
