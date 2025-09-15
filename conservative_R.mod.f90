!    GNOMO software (General Node Model)
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

! This module was constructed by Renske Vroomans in March 2018
! It contains simple comparative measures of similarity between an 
! individual and a "target" morphology. 
! This could either be used for conservative evolution or for evolution towards a target morphology
! EMD and D2 measures courtesy of PF H

module conservative_R

  use general
  use neighboring
  use io
  use genetic
  
  implicit none
  
  !variables
  integer :: tarindalloc=0 !flag for whether target individual has been allocated
  real*8, allocatable :: compind(:,:),tarind(:,:) !the individual to compare and the target individual.
  real*8 :: tarcentroidx=0d0, tarcentroidy=0d0, tarcentroidz=0d0 !target individual centroid
  integer :: tar_epind ! nr of epithelial nodes of target individual
  !***********************************************************
  
  contains
  
  !!We assume the individual to compare to has already been read and stored. 
  subroutine simple_EMD(center, EMDval, pernofi2, tarpernofi) !whether to use centering, storage var, the file with the individual to compare, target individual file (optional)
  integer::center
  real*8::EMDval
  character*140 :: pernofi2, tarpernofi
  real*8:: centroidx, centroidy, centroidz, req1, tarreq 
  integer:: epi_nd, epind, ndp 
  real*8::a,b,c, storenorm
  
  print *, "starting comparison..."
  
  if (len(trim(tarpernofi))==0 .and.tarindalloc==0) then
    print *, "conservative_R error -> simple_EMD: target individual not specified; set distance very large"
    EMDval=1d18
    return
  else if (tarindalloc==0) then !we have to read the file with the target individual now
   ! print*,pernofi2
    call read_targetind(center, tarpernofi)
  
  else if (len(trim(tarpernofi))/=0 .and.tarindalloc==1) then
    ! print *, "conservative_R warning -> simple_EMD: Target individual already specified, overriding stored morphology" !!>> HC 30-11-2020 Optimizing less prints
    call read_targetind(center, tarpernofi)
  endif
  
  !!! Read the individual to be compared to the target
  call readsnap(pernofi2)
  
  epi_nd=0
  do i=1,nd
     if(node(i)%tipus>2)cycle  
     epi_nd=epi_nd+1
  enddo
  
  if (nd>0) then 
  !now we transfer the coordinates of one individual to matrix compind
    if (allocated(compind)) deallocate(compind)
    allocate(compind(epi_nd,4))
         
    epind=0
               
    do i=1,nd
      if(node(i)%tipus>2)cycle  
      epind=epind+1
      compind(epind,1)=node(i)%x
      compind(epind,2)=node(i)%y
      compind(epind,3)=node(i)%z
      compind(epind,4)=node(i)%eqd  !!>>HC 30-6-2020
    end do
        
    if (center/=0) then  
      if (tarcentroidx==0d0) print *, "Conservative_R warning -> simple_EMD: tarind is not centered but compind is"
      
      centroidx=sum(compind(:,1))/epind     !centroid
      centroidy=sum(compind(:,2))/epind
      centroidz=sum(compind(:,3))/epind 
  
      do i=1,epind
        compind(i,1)=compind(i,1)-centroidx   !centering 
        compind(i,2)=compind(i,2)-centroidy
        compind(i,3)=compind(i,3)-centroidz
       
!         cox=sqrt(compind(i,1)**2)+cox !centroid size
!         coy=sqrt(compind(i,2)**2)+coy
!         coz=sqrt(compind(i,3)**2)+coz
        
      end do
      
      
    else if (tarcentroidx/=0d0) then
      print *, "Conservative_R warning -> simple_EMD: tarind is centered but compind is not"  
      
    end if
!       cox=cox/epind		!centroid size
!       coy=coy/epind
!       coz=coz/epind       
!    
!       centroid_size=(cox+tarcox+coy+tarcoy+coz+tarcoz)/6   !centroid of both morphologies 
  
  else 
    print *, "Conservative_R error -> simple_EMD: given individual has no nodes"
    stop
 
  end if
 
print*,'epind',epind

!now EMD
   d=0.0d0
   do i=1,epi_nd
        a=compind(i,1) ; b=compind(i,2) ; c=compind(i,3) ; req1=compind(i,4)
        aa=1000000
        do j=1,tar_epind
            dd=sqrt((tarind(j,1)-a)**2+(tarind(j,2)-b)**2+(tarind(j,3)-c)**2)
            if (dd<aa) then ; aa=dd ; storenorm=(tarind(j,4)+req1)/2 ; end if
        end do           
        d=d+aa/storenorm
   end do
   do i=1,tar_epind
        a=tarind(i,1) ; b=tarind(i,2) ; c=tarind(i,3) ; tarreq=tarind(i,4)
        aa=1000000
        do j=1,epi_nd
            dd=sqrt((compind(j,1)-a)**2+(compind(j,2)-b)**2+(compind(j,3)-c)**2)
            if (dd<aa) then; aa=dd ; storenorm=(compind(j,4)+tarreq)/2 ; end if
        end do          
        d=d+aa/storenorm
   enddo
   print*,'d ',d,aa,storenorm
   EMDval=d/real(epi_nd+tar_epind)
   !print*,'EMDval',EMDval
   !if (center) d2=d/centroid_size
   !print *, "I find a val of ", EMDval
  
  
end subroutine simple_EMD  


!!This function reads the nodes and req of the target individual and stores its values in global vars for 
!!later use by the function simple_EMD (or other comparison programs)
subroutine read_targetind(center, tarpernofi)
  integer:: center
  character*140 :: tarpernofi
  integer::epind
  
  !print *, "entering readtargetind"
  
  call readsnap(tarpernofi)
  print *, tarpernofi
  !how many epithelial nodes are there?
  tar_epind=0
  do i=1,nd
     if(node(i)%tipus>2)cycle  
     tar_epind=tar_epind+1
  enddo
  
  if (nd>0) then 
  !now we transfer the coordinates of the individual to matrix tarind
    if (allocated(tarind)) deallocate(tarind)
    allocate(tarind(tar_epind,4))
         
    epind=0
    
    do i=1,nd
      if(node(i)%tipus>2)cycle  
      epind=epind+1
      tarind(epind,1)=node(i)%x
      tarind(epind,2)=node(i)%y
      tarind(epind,3)=node(i)%z
      tarind(epind,4)=node(i)%eqd  !!>>HC 30-6-2020
    end do
        
    if (center/=0) then                    !we bring the tissue back to the center   
      tarcentroidx=sum(tarind(:,1))/epind     !centroid
      tarcentroidy=sum(tarind(:,2))/epind
      tarcentroidz=sum(tarind(:,3))/epind 
  
      do i=1,epind
        tarind(i,1)=tarind(i,1)-tarcentroidx   !centering 
        tarind(i,2)=tarind(i,2)-tarcentroidy
        tarind(i,3)=tarind(i,3)-tarcentroidz
       
      end do
    end if
    
    tarindalloc=1 !to signal that we stored the values 

  else
    print *, "Conservative_R error -> read_targetind: given individual has no nodes"
    stop
  
  end if

end subroutine read_targetind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!We assume the individual to compare to has already been read and stored. 
  subroutine EMD_only_with_original_cells(center, EMDval, pernofi2, tarpernofi)                                      !!>> HC 24-9-2021
  !whether to use centering, storage var, the file with the individual to compare, target individual file (optional) !!>> HC 24-9-2021
  integer::center, oripich                                                                                           !!>> HC 24-9-2021
  real*8::EMDval                                                                                                     !!>> HC 24-9-2021
  character*140 :: pernofi2, tarpernofi                                                                              !!>> HC 24-9-2021
  real*8:: centroidx, centroidy, centroidz, req1, tarreq                                                             !!>> HC 24-9-2021
  integer:: epind, ndp                                                                                               !!>> HC 24-9-2021
  integer:: ich, maxcelch                                                                                            !!>> HC 24-9-2021
  real*8::a,b,c, storenorm                                                                                           !!>> HC 24-9-2021
  
 ! print *, "starting comparison..."
    
  if (len(trim(tarpernofi))==0 .and.tarindalloc==0) then                                                        !!>> HC 24-9-2021
    print *, "conservative_R error -> simple_EMD: target individual not specified; set distance very large"     !!>> HC 24-9-2021
    EMDval=1d18                                                                                                 !!>> HC 24-9-2021
    return                                                                                                      !!>> HC 24-9-2021
  else if (tarindalloc==0) then !we have to read the file with the target individual now                        !!>> HC 24-9-2021
  
    call read_targetind_only_with_original_cells(center, tarpernofi)                                            !!>> HC 24-9-2021
  
  else if (len(trim(tarpernofi))/=0 .and.tarindalloc==1) then                                                   !!>> HC 24-9-2021
    ! print *, "conservative_R warning -> simple_EMD: Target individual already specified, overriding stored morphology" !!>> HC 30-11-2020 Optimizing less prints
    call read_targetind_only_with_original_cells(center, tarpernofi)                                            !!>> HC 24-9-2021
  endif                                                                                                         !!>> HC 24-9-2021
  
  !!! Read the individual to be compared to the target                                                          !!>> HC 24-9-2021
  
  call readsnap(pernofi2)                                                                                       !!>> HC 24-9-2021
  
  oripich=0                                                  !!>> HC 24-9-2021 Here we calculate the original number of epithelial nodes
  do ich=1,nd                                                !!>> HC 24-9-2021
     if (nodeo(ich)%tipus.ge.3)cycle                         !!>> HC 24-9-2021 WARNING 1 we are assuming that the epitelial cells are listed BEFORE the mesenchyme
     if (nodeo(ich)%icel>oripich) oripich=nodeo(ich)%icel    !!>> HC 24-9-2021 the original number of epithelial nodes is then the maximum number of epitelial cells
  enddo                                                      !!>> HC 24-9-2021 in the nodeo list (WARNING 2 assuming single_node cells)
  
  if (nd>0) then                                                            !!>> HC 24-9-2021
  !now we transfer the coordinates of one individual to matrix compind      !!>> HC 24-9-2021
    if (allocated(compind)) deallocate(compind)                             !!>> HC 24-9-2021
    allocate(compind(oripich,4))                                            !!>> HC 24-9-2021
         
    epind=0                                                                 !!>> HC 24-9-2021
               
    do i=1,nd                                                               !!>> HC 24-9-2021
      if(node(i)%tipus>2)cycle                                              !!>> HC 24-9-2021
      epind=epind+1                                                         !!>> HC 24-9-2021
      if (epind>oripich)exit                                                !!>> HC 24-9-2021
      compind(epind,1)=node(i)%x                                            !!>> HC 24-9-2021
      compind(epind,2)=node(i)%y                                            !!>> HC 24-9-2021
      compind(epind,3)=node(i)%z                                            !!>> HC 24-9-2021
      compind(epind,4)=node(i)%eqd  !!>>HC 30-6-2020                        !!>> HC 24-9-2021
    end do
        
    if (center/=0) then  
      if (tarcentroidx==0d0) print *, "Conservative_R warning -> simple_EMD: tarind is not centered but compind is" !!>> HC 24-9-2021
      centroidx=sum(compind(:,1))/epind                                     !!>> HC 24-9-2021 !centroid
      centroidy=sum(compind(:,2))/epind                                     !!>> HC 24-9-2021
      centroidz=sum(compind(:,3))/epind                                     !!>> HC 24-9-2021
      do i=1,epind                                                          !!>> HC 24-9-2021
        compind(i,1)=compind(i,1)-centroidx                                 !!>> HC 24-9-2021 !centering 
        compind(i,2)=compind(i,2)-centroidy                                 !!>> HC 24-9-2021
        compind(i,3)=compind(i,3)-centroidz                                 !!>> HC 24-9-2021
      end do                                                                !!>> HC 24-9-2021
      
    else if (tarcentroidx/=0d0) then                                        !!>> HC 24-9-2021
      print *, "Conservative_R warning -> simple_EMD: tarind is centered but compind is not" !!>> HC 24-9-2021  
    end if                                                                  !!>> HC 24-9-2021
  else                                                                      !!>> HC 24-9-2021
    print *, "Conservative_R error -> simple_EMD: given individual has no nodes"   !!>> HC 24-9-2021
    stop                                                                    !!>> HC 24-9-2021
  end if                                                                    !!>> HC 24-9-2021
 
!now EMD                                                                        !!>> HC 24-9-2021
   d=0.0d0                                                                      !!>> HC 24-9-2021
   do i=1,oripich                                                               !!>> HC 24-9-2021
        a=compind(i,1) ; b=compind(i,2) ; c=compind(i,3) ; req1=compind(i,4)    !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do j=1,oripich                                                          !!>> HC 24-9-2021
            dd=sqrt((tarind(j,1)-a)**2+(tarind(j,2)-b)**2+(tarind(j,3)-c)**2)   !!>> HC 24-9-2021
            if (dd<aa) then ; aa=dd ; storenorm=(tarind(j,4)+req1)/2 ; end if   !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   end do                                                                       !!>> HC 24-9-2021
   do i=1,oripich                                                               !!>> HC 24-9-2021
        a=tarind(i,1) ; b=tarind(i,2) ; c=tarind(i,3) ; tarreq=tarind(i,4)      !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do j=1,oripich                                                          !!>> HC 24-9-2021
            dd=sqrt((compind(j,1)-a)**2+(compind(j,2)-b)**2+(compind(j,3)-c)**2)!!>> HC 24-9-2021
            if (dd<aa) then; aa=dd ; storenorm=(compind(j,4)+tarreq)/2 ; end if !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   enddo                                                                        !!>> HC 24-9-2021
   EMDval=d/real(oripich+oripich)                                               !!>> HC 24-9-2021
  
end subroutine EMD_only_with_original_cells                                     !!>> HC 24-9-2021

!!This function reads the nodes and req of the target individual and stores its values in global vars for 
!!later use by the function simple_EMD (or other comparison programs)
subroutine read_targetind_only_with_original_cells(center, tarpernofi) !!>> HC 24-9-2021
  integer:: center, oripich, ich                                       !!>> HC 24-9-2021
  character*140 :: tarpernofi                                          !!>> HC 24-9-2021
  integer::epind                                                       !!>> HC 24-9-2021
  
  call readsnap(tarpernofi)                                            !!>> HC 24-9-2021
  
  
  if (nd>0) then                                                       !!>> HC 24-9-2021
  !now we transfer the coordinates of the individual to matrix tarind  !!>> HC 24-9-2021

     oripich=0                                                  !!>> HC 24-9-2021 Here we calculate the original number of epithelial nodes
     do ich=1,nd                                                !!>> HC 24-9-2021
        if (nodeo(ich)%tipus.ge.3)cycle                         !!>> HC 24-9-2021 WARNING 1 we are assuming that the epitelial cells are listed BEFORE the mesenchyma
        if (nodeo(ich)%icel>oripich) oripich=nodeo(ich)%icel    !!>> HC 24-9-2021 the original number of epithelial nodes is then the maximum number of epitelial cells
     enddo                                                      !!>> HC 24-9-2021 in the nodeo list (WARNING 2 assuming single_node cells)

    if (allocated(tarind)) deallocate(tarind)                   !!>> HC 24-9-2021
    allocate(tarind(oripich,4))                                 !!>> HC 24-9-2021
         
    epind=0                                                     !!>> HC 24-9-2021
    
    do i=1,nd                                                   !!>> HC 24-9-2021
      if(node(i)%tipus>2)cycle                                  !!>> HC 24-9-2021
      epind=epind+1                                             !!>> HC 24-9-2021
      if (epind>oripich)exit                                    !!>> HC 24-9-2021
      tarind(epind,1)=node(i)%x                                 !!>> HC 24-9-2021
      tarind(epind,2)=node(i)%y                                 !!>> HC 24-9-2021
      tarind(epind,3)=node(i)%z                                 !!>> HC 24-9-2021
      tarind(epind,4)=node(i)%eqd  !!>>HC 30-6-2020             !!>> HC 24-9-2021
    end do
        
    if (center/=0) then                                         !!>> HC 24-9-2021 we bring the tissue back to the center   
      tarcentroidx=sum(tarind(:,1))/epind                       !!>> HC 24-9-2021 !centroid
      tarcentroidy=sum(tarind(:,2))/epind                       !!>> HC 24-9-2021
      tarcentroidz=sum(tarind(:,3))/epind                       !!>> HC 24-9-2021
      do i=1,epind                                              !!>> HC 24-9-2021
        tarind(i,1)=tarind(i,1)-tarcentroidx                    !!>> HC 24-9-2021 !centering 
        tarind(i,2)=tarind(i,2)-tarcentroidy                    !!>> HC 24-9-2021
        tarind(i,3)=tarind(i,3)-tarcentroidz                    !!>> HC 24-9-2021
      end do                                                    !!>> HC 24-9-2021
    end if                                                      !!>> HC 24-9-2021
    
    tarindalloc=1                                               !!>> HC 24-9-2021 to signal that we stored the values 

  else                                                          !!>> HC 24-9-2021
    print *, "Conservative_R error -> read_targetind: given individual has no nodes" !!>> HC 24-9-2021
    stop                                                        !!>> HC 24-9-2021
  
  end if                                                        !!>> HC 24-9-2021

end subroutine read_targetind_only_with_original_cells          !!>> HC 24-9-2021

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!We assume the individual to compare to has already been read and stored. 
subroutine symmetry_hEMD(EMDval)                                      !!>> HC 24-9-2021
Implicit none
integer  ich, jch, kch, lch, oripich, half_oripich, epind
real*8, allocatable, dimension(:,:) ::compind, tarind, orinodes
real*8 :: EMDval, centroidx, centroidy, centroidz, req1, storenorm, tarreq

oripich=0                                                  !!>> HC 24-9-2021 Here we calculate the original number of epithelial nodes
do ich=1,nd                                                !!>> HC 24-9-2021
   if (nodeo(ich)%tipus.ge.3)cycle                         !!>> HC 24-9-2021 WARNING 1 we are assuming that the epitelial cells are listed BEFORE the mesenchyme
   if (nodeo(ich)%icel>oripich) oripich=nodeo(ich)%icel    !!>> HC 24-9-2021 the original number of epithelial nodes is then the maximum number of epitelial cells
enddo                                                      !!>> HC 24-9-2021 in the nodeo list (WARNING 2 assuming single_node cells)

half_oripich=ceiling(real(oripich)/2.0d0)
  
if (allocated(orinodes)) deallocate(orinodes)                !!>> HC 24-9-2021
allocate(orinodes(nd,4))                               !!>> HC 24-9-2021


if (allocated(compind)) deallocate(compind)                !!>> HC 24-9-2021
allocate(compind(oripich,4))                               !!>> HC 24-9-2021

if (allocated(tarind)) deallocate(tarind)                   !!>> HC 24-9-2021
allocate(tarind(oripich,4))                                 !!>> HC 24-9-2021

compind=0.0d0; tarind=0.0d0
lch=0; kch=0; epind=0; EMDval=0.0d0; orinodes=0.0d0

do ich=1,nd
   if(node(ich)%tipus.ne.2)cycle
   epind=epind+1
   if (epind>oripich)exit
   orinodes(ich,1)=node(ich)%x                                            !!>> HC 24-9-2021
   orinodes(ich,2)=node(ich)%y                                            !!>> HC 24-9-2021
   orinodes(ich,3)=node(ich)%z                                            !!>> HC 24-9-2021
   orinodes(ich,4)=node(ich)%eqd  !!>>HC 30-6-2020                        !!>> HC 24-9-2021
enddo
       
 centroidx=sum(orinodes(:,1))/epind                                  !!>> HC 24-9-2021 !centroid
 centroidy=sum(orinodes(:,2))/epind                                  !!>> HC 24-9-2021
 centroidz=sum(orinodes(:,3))/epind                                  !!>> HC 24-9-2021
 epind=0
 do ich=1,nd                                                         !!>> HC 24-9-2021
    if(node(ich)%tipus.ne.2)cycle
    epind=epind+1
    if (epind>oripich)exit
    orinodes(ich,1)=orinodes(ich,1)-centroidx                        !!>> HC 24-9-2021 !centering 
    orinodes(ich,2)=orinodes(ich,2)-centroidy                        !!>> HC 24-9-2021
    orinodes(ich,3)=orinodes(ich,3)-centroidz                        !!>> HC 24-9-2021
 end do                                                              !!>> HC 24-9-2021    

epind=0
do ich=1,nd
   if(node(ich)%tipus.ne.2)cycle
   epind=epind+1
   if (epind>oripich)exit
   if (nodeo(ich)%x>0.0d0)then
      lch=lch+1
      compind(lch,1)=orinodes(ich,1)                                            !!>> HC 24-9-2021
      compind(lch,2)=orinodes(ich,2)                                           !!>> HC 24-9-2021
      compind(lch,3)=orinodes(ich,3)                                            !!>> HC 24-9-2021
      compind(lch,4)=orinodes(ich,4) !!>>HC 30-6-2020                        !!>> HC 24-9-2021
   else
      kch=kch+1
      tarind(kch,1)=-orinodes(ich,1)                                            !!>> HC 24-9-2021
      tarind(kch,2)=orinodes(ich,2)                                            !!>> HC 24-9-2021
      tarind(kch,3)=orinodes(ich,3)                                            !!>> HC 24-9-2021
      tarind(kch,4)=orinodes(ich,4)  !!>>HC 30-6-2020                        !!>> HC 24-9-2021
   endif
enddo

!now EMD                                                                        !!>> HC 24-9-2021
   d=0.0d0                                                                      !!>> HC 24-9-2021
   do ich=1,lch                                                          !!>> HC 24-9-2021
        a=compind(ich,1) ; b=compind(ich,2) ; c=compind(ich,3) ; req1=compind(ich,4)    !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do jch=1,kch                                                     !!>> HC 24-9-2021
            dd=sqrt((tarind(jch,1)-a)**2+(tarind(jch,2)-b)**2+(tarind(jch,3)-c)**2)   !!>> HC 24-9-2021
            if (dd<aa) then ; aa=dd ; storenorm=(tarind(jch,4)+req1)/2 ; end if   !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   end do                                                                       !!>> HC 24-9-2021
   do ich=1,kch                                                          !!>> HC 24-9-2021
        a=tarind(ich,1) ; b=tarind(ich,2) ; c=tarind(ich,3) ; tarreq=tarind(ich,4)      !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do jch=1,lch                                                     !!>> HC 24-9-2021
            dd=sqrt((compind(jch,1)-a)**2+(compind(jch,2)-b)**2+(compind(jch,3)-c)**2)!!>> HC 24-9-2021
            if (dd<aa) then; aa=dd ; storenorm=(compind(jch,4)+tarreq)/2 ; end if !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   enddo                                                                        !!>> HC 24-9-2021
   EMDval=d/real(oripich)                                               !!>> HC 24-9-2021

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!We assume the individual to compare to has already been read and stored. 
subroutine symmetry_EMD(EMDval)                                      !!>> HC 24-9-2021
Implicit none
integer ::  ich, jch, kch, lch, epind
real*8, allocatable, dimension(:,:) ::compind,  tarind, orinodes
real*8 :: EMDval, centroidx, centroidy, centroidz, req1, storenorm, tarreq

if (allocated(orinodes)) deallocate(orinodes)                !!>> HC 24-9-2021
allocate(orinodes(nd,4))                               !!>> HC 24-9-2021


if (allocated(compind)) deallocate(compind)                !!>> HC 24-9-2021
allocate(compind(nd,4))                               !!>> HC 24-9-2021

if (allocated(tarind)) deallocate(tarind)                   !!>> HC 24-9-2021
allocate(tarind(nd,4))                                 !!>> HC 24-9-2021


compind=0.0d0; tarind=0.0d0
lch=0; kch=0; epind=0; EMDval=0.0d0; orinodes=0.0d0

do ich=1,nd
   if(node(ich)%tipus.ne.2)cycle
   epind=epind+1
   orinodes(ich,1)=node(ich)%x                                            !!>> HC 24-9-2021
   orinodes(ich,2)=node(ich)%y                                            !!>> HC 24-9-2021
   orinodes(ich,3)=node(ich)%z                                            !!>> HC 24-9-2021
   orinodes(ich,4)=node(ich)%eqd  !!>>HC 30-6-2020                        !!>> HC 24-9-2021
enddo
       
 centroidx=sum(orinodes(:,1))/epind                                  !!>> HC 24-9-2021 !centroid
 centroidy=sum(orinodes(:,2))/epind                                  !!>> HC 24-9-2021
 centroidz=sum(orinodes(:,3))/epind                                  !!>> HC 24-9-2021
 epind=0
 do ich=1,nd                                                         !!>> HC 24-9-2021
    if(node(ich)%tipus.ne.2)cycle
    orinodes(ich,1)=orinodes(ich,1)-centroidx                        !!>> HC 24-9-2021 !centering 
    orinodes(ich,2)=orinodes(ich,2)-centroidy                        !!>> HC 24-9-2021
    orinodes(ich,3)=orinodes(ich,3)-centroidz                        !!>> HC 24-9-2021
 end do                                                              !!>> HC 24-9-2021    

epind=0
do ich=1,nd
   if(node(ich)%tipus.ne.2)cycle
   epind=epind+1
   if (orinodes(ich,1)>0.0d0)then
      lch=lch+1
      compind(lch,1)=orinodes(ich,1)                                            !!>> HC 24-9-2021
      compind(lch,2)=orinodes(ich,2)                                           !!>> HC 24-9-2021
      compind(lch,3)=orinodes(ich,3)                                            !!>> HC 24-9-2021
      compind(lch,4)=orinodes(ich,4) !!>>HC 30-6-2020                        !!>> HC 24-9-2021
   else
      kch=kch+1
      tarind(kch,1)=-orinodes(ich,1)                                            !!>> HC 24-9-2021
      tarind(kch,2)=orinodes(ich,2)                                            !!>> HC 24-9-2021
      tarind(kch,3)=orinodes(ich,3)                                            !!>> HC 24-9-2021
      tarind(kch,4)=orinodes(ich,4)  !!>>HC 30-6-2020                        !!>> HC 24-9-2021
   endif
enddo

 

!now EMD                                                                        !!>> HC 24-9-2021
   d=0.0d0                                                                      !!>> HC 24-9-2021
   do ich=1,lch                                                          !!>> HC 24-9-2021
        a=compind(ich,1) ; b=compind(ich,2) ; c=compind(ich,3) ; req1=compind(ich,4)    !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do jch=1,kch                                                     !!>> HC 24-9-2021
            dd=sqrt((tarind(jch,1)-a)**2+(tarind(jch,2)-b)**2+(tarind(jch,3)-c)**2)   !!>> HC 24-9-2021
            if (dd<aa) then ; aa=dd ; storenorm=(tarind(jch,4)+req1)/2 ; end if   !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   end do                                                                       !!>> HC 24-9-2021
   do ich=1,kch                                                          !!>> HC 24-9-2021
        a=tarind(ich,1) ; b=tarind(ich,2) ; c=tarind(ich,3) ; tarreq=tarind(ich,4)      !!>> HC 24-9-2021
        aa=1000000                                                              !!>> HC 24-9-2021
        do jch=1,lch                                                     !!>> HC 24-9-2021
            dd=sqrt((compind(jch,1)-a)**2+(compind(jch,2)-b)**2+(compind(jch,3)-c)**2)!!>> HC 24-9-2021
            if (dd<aa) then; aa=dd ; storenorm=(compind(jch,4)+tarreq)/2 ; end if !!>> HC 24-9-2021
        end do                                                                  !!>> HC 24-9-2021
        d=d+aa/storenorm                                                        !!>> HC 24-9-2021
   enddo                                                                        !!>> HC 24-9-2021
   EMDval=d/real(epind)                                               !!>> HC 24-9-2021
   
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ranges(rangfile,onetwin)
!!>> HC 28-9-2021 This subroutine was created by Hugo Cano to make sure that the node properties in the final embryos
!!>> HC 28-9-2021 are still in the ranges of property activation that are biologicaly realistic
  implicit none
  integer :: limit
  integer :: ich,jch,lch,ord,ord2,unused  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
  real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
  integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  real :: wrongval, maxl, minl,ratio !!>> HC 28-9-2021
  character*400 cx, rangfile,onetwin !!>> HC 29-9-2021
  real*8, dimension (1:5) :: min_glim, max_glim   !!>> HC 6-10-2021
  
  call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2021 Read ranges
             
  limit=0                                               !!>> HC 27-11-2020 if limit stays in 0 no limits have been surpased if =0, we will the embryo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  do ich=1,nd                                                                                !!>> HC 28-9-2021 Here we check whether any node has surpased the limits
     !! 6 ADD                                                                                !!>> HC 28-9-2021 of node properties that correspond to the maximum/minimum activation
     ratio=node(ich)%add/node(ich)%eqd                                                       !!>> HC 28-9-2021 in a gene with concentration = 1
     if (ratio>2.3)then; limit=1; lch=6; wrongval=ratio;  exit; endif                        !!>> HC 28-9-2021
     if (ratio<1.2)then; limit=2; lch=6; wrongval=ratio;  exit; endif                        !!>> HC 28-9-2021 In ADD however, we check the ADD/EHD ration
                                                                                             !!>> HC 28-9-2021
     !! 8 ADH                                                                                !!>> HC 28-9-2021
     maxl=max_elim(8)+nodeo(ich)%adh                                                         !!>> HC 28-9-2021 This is the maximum realistic value of adh
     minl=min_elim(8)+nodeo(ich)%adh                                                         !!>> HC 28-9-2021 This is the minimum realistic value of adh
     if (node(ich)%adh>maxl)then; limit=1; lch=8; wrongval=node(ich)%adh;  exit; endif       !!>> HC 28-9-2021
     if (node(ich)%adh<minl)then; limit=2; lch=8; wrongval=node(ich)%adh;  exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 10 REC                                                                               !!>> HC 28-9-2021 And the same for all the properties
     maxl=max_elim(10)+nodeo(ich)%rec                                                        !!>> HC 28-9-2021
     minl=min_elim(10)+nodeo(ich)%rec                                                        !!>> HC 28-9-2021
     if (node(ich)%rec>maxl)then; limit=1; lch=10; wrongval=node(ich)%rec; exit; endif       !!>> HC 28-9-2021
     if (node(ich)%rec<minl)then; limit=2; lch=10; wrongval=node(ich)%rec; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 11 ERP                                                                               !!>> HC 28-9-2021
     maxl=max_elim(11)+nodeo(ich)%erp                                                        !!>> HC 28-9-2021
     minl=min_elim(11)+nodeo(ich)%erp                                                        !!>> HC 28-9-2021
     if (node(ich)%erp>maxl)then; limit=1; lch=11; wrongval=node(ich)%erp; exit; endif       !!>> HC 28-9-2021
     if (node(ich)%erp<minl)then; limit=2; lch=11; wrongval=node(ich)%erp; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 12 EST                                                                               !!>> HC 28-9-2021
     maxl=max_elim(12)+nodeo(ich)%est                                                        !!>> HC 28-9-2021
     minl=min_elim(12)+nodeo(ich)%est                                                        !!>> HC 28-9-2021
     if (node(ich)%est>maxl)then; limit=1; lch=12; wrongval=node(ich)%est; exit; endif       !!>> HC 28-9-2021
     if (node(ich)%est<minl)then; limit=2; lch=12; wrongval=node(ich)%est; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 13 EQS                                                                               !!>> HC 28-9-2021
     maxl=max_elim(13)+nodeo(ich)%eqs                                                        !!>> HC 28-9-2021
     minl=min_elim(13)+nodeo(ich)%eqs                                                        !!>> HC 28-9-2021
     if (node(ich)%eqs>maxl)then; limit=1; lch=13; wrongval=node(ich)%eqs; exit; endif       !!>> HC 28-9-2021
     if (node(ich)%eqs<minl)then; limit=2; lch=13; wrongval=node(ich)%eqs; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 14 HOO                                                                               !!>> HC 28-9-2021
     maxl=max_elim(14)+nodeo(ich)%hoo                                                        !!>> HC 28-9-2021
     minl=min_elim(14)+nodeo(ich)%hoo                                                        !!>> HC 28-9-2021
     if (node(ich)%hoo>maxl)then; limit=1; lch=14; wrongval=node(ich)%hoo; exit; endif       !!>> HC 28-9-2021
     if (node(ich)%hoo<minl)then; limit=2; lch=14; wrongval=node(ich)%hoo; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 21 COD                                                                               !!>> HC 28-9-2021
     maxl=max_elim(21)+nodeo(ich)%cod                                                        !!>> HC 28-9-2021
     minl=min_elim(21)+nodeo(ich)%cod                                                        !!>> HC 28-9-2021
     if (node(ich)%cod>maxl)then; limit=1; lch=21; wrongval=node(ich)%cod; exit;  endif      !!>> HC 28-9-2021
     if (node(ich)%cod<minl)then; limit=2; lch=21; wrongval=node(ich)%cod; exit;  endif      !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
     !! 22 GRD                                                                               !!>> HC 28-9-2021
     maxl=max_elim(22)+nodeo(ich)%grd                                                        !!>> HC 28-9-2021
     minl=min_elim(22)+nodeo(ich)%grd                                                        !!>> HC 28-9-2021
     if (node(ich)%grd>maxl)then; limit=1; lch=22; wrongval=node(ich)%grd; exit;  endif      !!>> HC 28-9-2021
     if (node(ich)%grd<minl)then; limit=2; lch=22; wrongval=node(ich)%grd; exit;  endif      !!>> HC 28-9-2021
  enddo                                                                                      !!>> HC 28-9-2021
  
  
  if (limit>0)then                            !!>> HC 28-9-2021 If limits have been surpased, we kill this individual
     open(819,file="individual.datfitness")   !!>> HC 28-9-2021
     write(819,*) "0.00"                      !!>> HC 28-9-2021
     close(819)                               !!>> HC 28-9-2021
  endif                                       !!>> HC 28-9-2021
  
  if (limit==1)then                                                                                         !!>> HC 28-9-2021 If limits have been surpased,
     open(818,file=trim(onetwin)//".loga" )                                                                 !!>> HC 28-9-2021 we write a logfile
     write(818,*) "Surpased the positive activation of the property/behaviour", lch, "value", wrongval      !!>> HC 28-9-2021
     close(818)                                                                                             !!>> HC 28-9-2021
     cx="pwd >> "//trim(onetwin)//".logb"                                                                   !!>> HC 6-10-2021 Save the name of the individual (evolution)
     call system(cx)                                                                                        !!>> HC 29-9-2021
     cx="paste "//trim(onetwin)//".logb "//trim(onetwin)//".loga > "//trim(onetwin)//".log"                 !!>> HC 6-10-2021 paste the name and the reason to end
     call system(cx)                                                                                        !!>> HC 29-9-2021 Remove the temporal files
     cx="rm "//trim(onetwin)//".logb"                                                                       !!>> HC 6-10-2021
     call system(cx)                                                                                        !!>> HC 29-9-2021
     cx="rm "//trim(onetwin)//".loga"                                                                       !!>> HC 6-10-2021
     call system(cx)                                                                                        !!>> HC 29-9-2021
  elseif(limit==2)then                                                                                      !!>> HC 28-9-2021
     open(818,file=trim(onetwin)//".loga")                                                                  !!>> HC 6-10-2021
     write(818,*) "Surpased the negative activation of the property/behaviour", lch, "value", wrongval      !!>> HC 28-9-2021
     close(818)                                                                                             !!>> HC 28-9-2021
     cx="pwd >> "//trim(onetwin)//".logb"                                                                   !!>> HC 6-10-2021 Save the name of the individual (evolution)
     call system(cx)                                                                                        !!>> HC 29-9-2021
     cx="paste "//trim(onetwin)//".logb "//trim(onetwin)//".loga > "//trim(onetwin)//".log"                 !!>> HC 6-10-2021 paste the name and the reason to end
     call system(cx)                                                                                        !!>> HC 29-9-2021 Remove the temporal files
     cx="rm "//trim(onetwin)//".logb"                                                                       !!>> HC 6-10-2021
     call system(cx)                                                                                        !!>> HC 29-9-2021
     cx="rm "//trim(onetwin)//".loga"                                                                       !!>> HC 6-10-2021
     call system(cx)                                                                                        !!>> HC 29-9-2021
  endif                                                                                                     !!>> HC 28-9-2021
end subroutine

subroutine assymetry(dieornot)                                                                !!>> AL 19-4-2024 
  integer :: dieornot
  real*8 :: fithc, ifithc, mefithc, emd                                                       !!>> HC 3-7-2023
  real*8 :: osfithc, otfithc, offithc, stfithc, sffithc, tffithc, symch                       !!>> HC 26-11-2020
  character*140 :: onetwin, comando                                                           !!>> AL 19-4-2024    
  integer ::  ihc, jhc, khc, nnodes, naltech, surface, volume, outside, total, inside, differ, centrow
  integer ::  iihc,jjhc,kkhc, ord1, ord2, ord3, lhc, prev, changes, mhc, nhc, neich, neichi, newdots, vhc, phc
  real*8, dimension(1:3) ::  u
  real*8 :: upv, sumd,  modu, ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, pershared
  real*8 :: maxx, minx, maxy, miny, maxz, minz, maxadd, minadd, anchx, anchy, anchz, anchadd, totmax, totmin, totanch
  real*8, allocatable, dimension(:,:) :: ncoords
  integer, allocatable, dimension(:,:,:) :: bfillz
  logical :: existes                                                                          !!>> AL 19-4-24
  fithc=0.0d0; ifithc=0.0d0 ; mefithc=0.0d0                                                   !!>> HC 26-11-2020
  osfithc=0.0d0; otfithc=0.0d0; offithc=0.0d0; stfithc=0.0d0; sffithc=0.0d0; tffithc=0.0d0    !!>> HC 26-11-2020
  distfitscale=2.250d0; distfitmag=1.0d0
  !call getarg(1,onetwin)

  !print*, onetwin
  
  !call iniread
  !call readsnap(onetwin)
  call iniboxes
  call neighbor_build
  
maxx=0.0d0; minx=10000; maxy=0.0d0; miny=10000; maxz=0.0d0; minz=10000                    !!>> HC 4-3-2024
anchx=0.0d0; anchy=0.0d0; anchz=0.0d0; totmax=0.0d0; totmin=0.0d0; totanch=0.0d0          !!>> HC 4-3-2024
maxadd=0.0d0                                                                              !!>> HC 4-3-2024

newdots=0                                                                                 !!>> HC 4-3-2024 NUMBER OF NEW POINTS
do ihc=1,nd                                                                               !!>> HC 4-3-2024 placed randomly to fill in gaps in the epithelium
   if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
   newdots=newdots+1                                                                      !!>> HC 4-3-2024
   do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
      neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
      if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
      do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
         neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
         if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
         if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
            newdots=newdots+10                                                            !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

if(allocated(ncoords))deallocate(ncoords)                                                 !!>> HC 4-3-2024 This stores the coordinates of the morphology
allocate(ncoords(1:newdots,1:3))                                                          !!>> HC 4-3-2024  and the newly added random points
ncoords=0.0d0                                                                             !!>> HC 4-3-2024

ord1=0                                                                                    !!>> HC 4-3-2024 STORE THE COORDINATES OF NODES
do ihc=1,nd                                                                               !!>> HC 4-3-2024
   if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
   ord1=ord1+1                                                                            !!>> HC 4-3-2024
   ncoords(ord1,1)=node(ihc)%x                                                            !!>> HC 4-3-2024
   ncoords(ord1,2)=node(ihc)%y                                                            !!>> HC 4-3-2024
   ncoords(ord1,3)=node(ihc)%z                                                            !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

do ihc=1,nd                                                                               !!>> HC 4-3-2024 ADD THE RANDOM POINTS
   if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
   do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
      neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
      if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
      do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
         neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
         if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
         if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
            do phc=1,10                                                                   !!>> HC 4-3-2024
                ahc = node(neich)%x - node(ihc)%x                                         !!>> HC 4-3-2024
                bhc = node(neich)%y - node(ihc)%y                                         !!>> HC 4-3-2024
                chc = node(neich)%z - node(ihc)%z                                         !!>> HC 4-3-2024
                
                ahc2 = node(neichi)%x - node(ihc)%x                                       !!>> HC 4-3-2024
                bhc2 = node(neichi)%y - node(ihc)%y                                       !!>> HC 4-3-2024
                chc2 = node(neichi)%z - node(ihc)%z                                       !!>> HC 4-3-2024
                
                call random_number(u1); call random_number(u2);                           !!>> HC 4-3-2024
                if ( (u1+u2) > 1.0d0)then                                                 !!>> HC 4-3-2024
                   u1=1-u1                                                                !!>> HC 4-3-2024
                   u2=1-u2                                                                !!>> HC 4-3-2024
                endif                                                                     !!>> HC 4-3-2024
                
                ahc3 = u1*ahc + u2*ahc2                                                   !!>> HC 4-3-2024
                bhc3 = u1*bhc + u2*bhc2                                                   !!>> HC 4-3-2024
                chc3 = u1*chc + u2*chc2                                                   !!>> HC 4-3-2024
                ord1=ord1+1                                                               !!>> HC 4-3-2024
                ncoords(ord1,1)= node(ihc)%x + ahc3                                       !!>> HC 4-3-2024
                ncoords(ord1,2)= node(ihc)%y + bhc3                                       !!>> HC 4-3-2024
                ncoords(ord1,3)= node(ihc)%z + chc3                                       !!>> HC 4-3-2024
            enddo                                                                         !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

do ihc=1,nd                                                                               !!>> HC 4-3-2024 CALCULATE THE BIGGEST COORDINATE IN THE MORPH
   if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)                                        !!>> HC 4-3-2024
   if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)                                        !!>> HC 4-3-2024
   if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)                                        !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

totmax=maxx                                                                               !!>> HC 4-3-2024
if (maxy>totmax) totmax=maxy                                                              !!>> HC 4-3-2024
if (maxz>totmax) totmax=maxz                                                              !!>> HC 4-3-2024 The biggest coordinate

maxadd=nodeo(1)%add*2                                                                     !!>> HC 4-3-2024
urv=1/maxadd                                                                              !!>> HC 4-3-2024
nboxes=(urv*totmax) +2                                                                    !!>> HC 4-3-2024  THE NUMBER OF BOXES

if(allocated(bfillz))deallocate(bfillz)                                                   !!>> HC 4-3-2024 THIS ARRAY STORES THE BOXES
allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))                            !!>> HC 4-3-2024
bfillz=0                                                                                  !!>> HC 4-3-2024

do ihc=1,ord1                                                                             !!>> HC 4-3-2024 INTRODUCE THE COORDINATES IN THEIR CORRESPONDING BOXES
   jhc=nint(ncoords(ihc,1)*urv)                                                           !!>> HC 4-3-2024 
   khc=nint(ncoords(ihc,2)*urv)                                                           !!>> HC 4-3-2024
   nhc=nint(ncoords(ihc,3)*urv)                                                           !!>> HC 4-3-2024
   bfillz(jhc,khc,nhc)=1                                                                  !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

bfillz(nboxes,nboxes,nboxes)=2                                                            !!>> HC 4-3-2024 THE FILLING ALGORITHM
changes=11; ord1=0                                                                        !!>> HC 4-3-2024 Fills with either 2 or 3 the boxes OUTSIDE the morphology
do while(changes>0)                                                                       !!>> HC 4-3-2024 2-3= OUT
   changes=0                                                                              !!>> HC 4-3-2024 0= IN
   do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024 1= SURFACE
      do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
         do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
            if( bfillz(ihc,jhc,khc).ne.2)cycle                                            !!>> HC 4-3-2024
            do nhc=-1,1                                                                   !!>> HC 4-3-2024
               ord1=ihc+nhc                                                               !!>> HC 4-3-2024
               if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
               if(bfillz(ord1,jhc,khc).ne.0)cycle                                         !!>> HC 4-3-2024
               bfillz(ord1,jhc,khc)=2                                                     !!>> HC 4-3-2024
               changes=changes+1                                                          !!>> HC 4-3-2024
            enddo                                                                         !!>> HC 4-3-2024
            do nhc=-1,1                                                                   !!>> HC 4-3-2024
               ord1=jhc+nhc                                                               !!>> HC 4-3-2024
               if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
               if(bfillz(ihc,ord1,khc).ne.0)cycle                                         !!>> HC 4-3-2024
               bfillz(ihc,ord1,khc)=2                                                     !!>> HC 4-3-2024
               changes=changes+1                                                          !!>> HC 4-3-2024
            enddo                                                                         !!>> HC 4-3-2024
            do nhc=-1,1                                                                   !!>> HC 4-3-2024
               ord1=khc+nhc                                                               !!>> HC 4-3-2024
               if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
               if(bfillz(ihc,jhc,ord1).ne.0)cycle                                         !!>> HC 4-3-2024
               bfillz(ihc,jhc,ord1)=2                                                     !!>> HC 4-3-2024
               changes=changes+1                                                          !!>> HC 4-3-2024
            enddo                                                                         !!>> HC 4-3-2024
            bfillz(ihc,jhc,khc)=3                                                         !!>> HC 4-3-2024
         enddo                                                                            !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

surface=0; volume=0; outside=0; total=0                                                   !!>> HC 4-3-2024  Counting outside and inside boxes
do ihc=-nboxes,nboxes                                                                     !!>> HC 4-3-2024
   do jhc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
      do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
         if (bfillz(ihc,jhc,khc)==1)then                                                  !!>> HC 4-3-2024
            surface=surface+1                                                             !!>> HC 4-3-2024
         elseif(bfillz(ihc,jhc,khc)==0)then                                               !!>> HC 4-3-2024
            volume=volume+1                                                               !!>> HC 4-3-2024
         else                                                                             !!>> HC 4-3-2024
            outside=outside+1                                                             !!>> HC 4-3-2024
         endif                                                                            !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024
enddo                                                                                     !!>> HC 4-3-2024

if (volume>0)then                                                                         !!>> HC 4-3-2024 If the morphology is not broken
   total=surface+volume                                                                   !!>> HC 4-3-2024 CALCULATE THE ASYMMETRY BY SHARED VOLUME
   do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
      do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
         do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024 Transform all the boxes into 1=IN 0=out
            if ( bfillz(ihc,jhc,khc).le.1)then                                            !!>> HC 4-3-2024
               bfillz(ihc,jhc,khc)=1                                                      !!>> HC 4-3-2024
            else                                                                          !!>> HC 4-3-2024
               bfillz(ihc,jhc,khc)=0                                                      !!>> HC 4-3-2024
            endif                                                                         !!>> HC 4-3-2024
         enddo                                                                            !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024

   centrow=0                                                                              !!>> HC 4-3-2024 NON-SHARED VOLUME
   do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
      do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
      if ( bfillz(ihc,0,jhc)==0)cycle                                                     !!>> HC 4-3-2024
      centrow=centrow+1                                                                   !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024

   differ=0                                                                               !!>> HC 4-3-2024
   do ihc=1,nboxes                                                                        !!>> HC 4-3-2024
      jhc=-ihc                                                                            !!>> HC 4-3-2024
      do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
         do nhc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
            if( bfillz(ihc,khc,nhc) == bfillz(jhc,khc,nhc))cycle                          !!>> HC 4-3-2024
            differ=differ+1                                                               !!>> HC 4-3-2024
         enddo                                                                            !!>> HC 4-3-2024
      enddo                                                                               !!>> HC 4-3-2024
   enddo                                                                                  !!>> HC 4-3-2024
   symch=real(differ)/real(total-centrow)                                                 !!>> HC 4-3-2024
else                                                                                      !!>> HC 4-3-2024 if the morphology is broken we will discard it
   symch=666.0d0                                                                          !!>> HC 4-3-2024
endif                                                                                     !!>> HC 4-3-2024
  
if (symch<0.0500d0) then                                                                  !!>> HC 20-12-2020 We only want robust individuals
    dieornot=0                                                                            !!>> HC 20-12-2020 The threshold comes my observations
else                                                                                      !!>> HC 20-12-2020  
    dieornot=1                                                                            !!>> HC 20-12-2020 Avoid extintion
endif                                                                                     !!>> HC 20-12-2020
  
end subroutine assymetry
 
subroutine assymetry_save(dieornot)                                                                !!>> AL 19-4-2024 
   integer :: dieornot
   real*8 :: fithc, ifithc, mefithc, emd                                                       !!>> HC 3-7-2023
   real*8 :: osfithc, otfithc, offithc, stfithc, sffithc, tffithc, symch                       !!>> HC 26-11-2020
   character*140 :: onetwin, comando                                                           !!>> AL 19-4-2024    
   integer ::  ihc, jhc, khc, nnodes, naltech, surface, volume, outside, total, inside, differ, centrow
   integer ::  iihc,jjhc,kkhc, ord1, ord2, ord3, lhc, prev, changes, mhc, nhc, neich, neichi, newdots, vhc, phc
   real*8, dimension(1:3) ::  u
   real*8 :: upv, sumd,  modu, ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, pershared
   real*8 :: maxx, minx, maxy, miny, maxz, minz, maxadd, minadd, anchx, anchy, anchz, anchadd, totmax, totmin, totanch
   real*8, allocatable, dimension(:,:) :: ncoords
   integer, allocatable, dimension(:,:,:) :: bfillz
   logical :: existes                                                                          !!>> AL 19-4-24
   fithc=0.0d0; ifithc=0.0d0 ; mefithc=0.0d0                                                   !!>> HC 26-11-2020
   osfithc=0.0d0; otfithc=0.0d0; offithc=0.0d0; stfithc=0.0d0; sffithc=0.0d0; tffithc=0.0d0    !!>> HC 26-11-2020
   distfitscale=2.250d0; distfitmag=1.0d0

   call iniboxes
   call neighbor_build
   
 maxx=0.0d0; minx=10000; maxy=0.0d0; miny=10000; maxz=0.0d0; minz=10000                    !!>> HC 4-3-2024
 anchx=0.0d0; anchy=0.0d0; anchz=0.0d0; totmax=0.0d0; totmin=0.0d0; totanch=0.0d0          !!>> HC 4-3-2024
 maxadd=0.0d0                                                                              !!>> HC 4-3-2024
 
 newdots=0                                                                                 !!>> HC 4-3-2024 NUMBER OF NEW POINTS
 do ihc=1,nd                                                                               !!>> HC 4-3-2024 placed randomly to fill in gaps in the epithelium
    if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
    newdots=newdots+1                                                                      !!>> HC 4-3-2024
    do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
       neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
       if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
       do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
          neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
          if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
          if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
             newdots=newdots+10                                                            !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 if(allocated(ncoords))deallocate(ncoords)                                                 !!>> HC 4-3-2024 This stores the coordinates of the morphology
 allocate(ncoords(1:newdots,1:3))                                                          !!>> HC 4-3-2024  and the newly added random points
 ncoords=0.0d0                                                                             !!>> HC 4-3-2024
 
 ord1=0                                                                                    !!>> HC 4-3-2024 STORE THE COORDINATES OF NODES
 do ihc=1,nd                                                                               !!>> HC 4-3-2024
    if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
    ord1=ord1+1                                                                            !!>> HC 4-3-2024
    ncoords(ord1,1)=node(ihc)%x                                                            !!>> HC 4-3-2024
    ncoords(ord1,2)=node(ihc)%y                                                            !!>> HC 4-3-2024
    ncoords(ord1,3)=node(ihc)%z                                                            !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 do ihc=1,nd                                                                               !!>> HC 4-3-2024 ADD THE RANDOM POINTS
    if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
    do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
       neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
       if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
       do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
          neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
          if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
          if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
             do phc=1,10                                                                   !!>> HC 4-3-2024
                 ahc = node(neich)%x - node(ihc)%x                                         !!>> HC 4-3-2024
                 bhc = node(neich)%y - node(ihc)%y                                         !!>> HC 4-3-2024
                 chc = node(neich)%z - node(ihc)%z                                         !!>> HC 4-3-2024
                 
                 ahc2 = node(neichi)%x - node(ihc)%x                                       !!>> HC 4-3-2024
                 bhc2 = node(neichi)%y - node(ihc)%y                                       !!>> HC 4-3-2024
                 chc2 = node(neichi)%z - node(ihc)%z                                       !!>> HC 4-3-2024
                 
                 call random_number(u1); call random_number(u2);                           !!>> HC 4-3-2024
                 if ( (u1+u2) > 1.0d0)then                                                 !!>> HC 4-3-2024
                    u1=1-u1                                                                !!>> HC 4-3-2024
                    u2=1-u2                                                                !!>> HC 4-3-2024
                 endif                                                                     !!>> HC 4-3-2024
                 
                 ahc3 = u1*ahc + u2*ahc2                                                   !!>> HC 4-3-2024
                 bhc3 = u1*bhc + u2*bhc2                                                   !!>> HC 4-3-2024
                 chc3 = u1*chc + u2*chc2                                                   !!>> HC 4-3-2024
                 ord1=ord1+1                                                               !!>> HC 4-3-2024
                 ncoords(ord1,1)= node(ihc)%x + ahc3                                       !!>> HC 4-3-2024
                 ncoords(ord1,2)= node(ihc)%y + bhc3                                       !!>> HC 4-3-2024
                 ncoords(ord1,3)= node(ihc)%z + chc3                                       !!>> HC 4-3-2024
             enddo                                                                         !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 do ihc=1,nd                                                                               !!>> HC 4-3-2024 CALCULATE THE BIGGEST COORDINATE IN THE MORPH
    if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)                                        !!>> HC 4-3-2024
    if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)                                        !!>> HC 4-3-2024
    if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)                                        !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 totmax=maxx                                                                               !!>> HC 4-3-2024
 if (maxy>totmax) totmax=maxy                                                              !!>> HC 4-3-2024
 if (maxz>totmax) totmax=maxz                                                              !!>> HC 4-3-2024 The biggest coordinate
 
 maxadd=nodeo(1)%add*2                                                                     !!>> HC 4-3-2024
 urv=1/maxadd                                                                              !!>> HC 4-3-2024
 nboxes=(urv*totmax) +2                                                                    !!>> HC 4-3-2024  THE NUMBER OF BOXES
 
 if(allocated(bfillz))deallocate(bfillz)                                                   !!>> HC 4-3-2024 THIS ARRAY STORES THE BOXES
 allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))                            !!>> HC 4-3-2024
 bfillz=0                                                                                  !!>> HC 4-3-2024
 
 do ihc=1,ord1                                                                             !!>> HC 4-3-2024 INTRODUCE THE COORDINATES IN THEIR CORRESPONDING BOXES
    jhc=nint(ncoords(ihc,1)*urv)                                                           !!>> HC 4-3-2024 
    khc=nint(ncoords(ihc,2)*urv)                                                           !!>> HC 4-3-2024
    nhc=nint(ncoords(ihc,3)*urv)                                                           !!>> HC 4-3-2024
    bfillz(jhc,khc,nhc)=1                                                                  !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 bfillz(nboxes,nboxes,nboxes)=2                                                            !!>> HC 4-3-2024 THE FILLING ALGORITHM
 changes=11; ord1=0                                                                        !!>> HC 4-3-2024 Fills with either 2 or 3 the boxes OUTSIDE the morphology
 do while(changes>0)                                                                       !!>> HC 4-3-2024 2-3= OUT
    changes=0                                                                              !!>> HC 4-3-2024 0= IN
    do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024 1= SURFACE
       do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
          do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
             if( bfillz(ihc,jhc,khc).ne.2)cycle                                            !!>> HC 4-3-2024
             do nhc=-1,1                                                                   !!>> HC 4-3-2024
                ord1=ihc+nhc                                                               !!>> HC 4-3-2024
                if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                if(bfillz(ord1,jhc,khc).ne.0)cycle                                         !!>> HC 4-3-2024
                bfillz(ord1,jhc,khc)=2                                                     !!>> HC 4-3-2024
                changes=changes+1                                                          !!>> HC 4-3-2024
             enddo                                                                         !!>> HC 4-3-2024
             do nhc=-1,1                                                                   !!>> HC 4-3-2024
                ord1=jhc+nhc                                                               !!>> HC 4-3-2024
                if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                if(bfillz(ihc,ord1,khc).ne.0)cycle                                         !!>> HC 4-3-2024
                bfillz(ihc,ord1,khc)=2                                                     !!>> HC 4-3-2024
                changes=changes+1                                                          !!>> HC 4-3-2024
             enddo                                                                         !!>> HC 4-3-2024
             do nhc=-1,1                                                                   !!>> HC 4-3-2024
                ord1=khc+nhc                                                               !!>> HC 4-3-2024
                if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                if(bfillz(ihc,jhc,ord1).ne.0)cycle                                         !!>> HC 4-3-2024
                bfillz(ihc,jhc,ord1)=2                                                     !!>> HC 4-3-2024
                changes=changes+1                                                          !!>> HC 4-3-2024
             enddo                                                                         !!>> HC 4-3-2024
             bfillz(ihc,jhc,khc)=3                                                         !!>> HC 4-3-2024
          enddo                                                                            !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 surface=0; volume=0; outside=0; total=0                                                   !!>> HC 4-3-2024  Counting outside and inside boxes
 do ihc=-nboxes,nboxes                                                                     !!>> HC 4-3-2024
    do jhc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
       do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
          if (bfillz(ihc,jhc,khc)==1)then                                                  !!>> HC 4-3-2024
             surface=surface+1                                                             !!>> HC 4-3-2024
          elseif(bfillz(ihc,jhc,khc)==0)then                                               !!>> HC 4-3-2024
             volume=volume+1                                                               !!>> HC 4-3-2024
          else                                                                             !!>> HC 4-3-2024
             outside=outside+1                                                             !!>> HC 4-3-2024
          endif                                                                            !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 enddo                                                                                     !!>> HC 4-3-2024
 
 if (volume>0)then                                                                         !!>> HC 4-3-2024 If the morphology is not broken
    total=surface+volume                                                                   !!>> HC 4-3-2024 CALCULATE THE ASYMMETRY BY SHARED VOLUME
    do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
       do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
          do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024 Transform all the boxes into 1=IN 0=out
             if ( bfillz(ihc,jhc,khc).le.1)then                                            !!>> HC 4-3-2024
                bfillz(ihc,jhc,khc)=1                                                      !!>> HC 4-3-2024
             else                                                                          !!>> HC 4-3-2024
                bfillz(ihc,jhc,khc)=0                                                      !!>> HC 4-3-2024
             endif                                                                         !!>> HC 4-3-2024
          enddo                                                                            !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 
    centrow=0                                                                              !!>> HC 4-3-2024 NON-SHARED VOLUME
    do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
       do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
       if ( bfillz(ihc,0,jhc)==0)cycle                                                     !!>> HC 4-3-2024
       centrow=centrow+1                                                                   !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
 
    differ=0                                                                               !!>> HC 4-3-2024
    do ihc=1,nboxes                                                                        !!>> HC 4-3-2024
       jhc=-ihc                                                                            !!>> HC 4-3-2024
       do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
          do nhc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
             if( bfillz(ihc,khc,nhc) == bfillz(jhc,khc,nhc))cycle                          !!>> HC 4-3-2024
             differ=differ+1                                                               !!>> HC 4-3-2024
          enddo                                                                            !!>> HC 4-3-2024
       enddo                                                                               !!>> HC 4-3-2024
    enddo                                                                                  !!>> HC 4-3-2024
    symch=real(differ)/real(total-centrow)                                                 !!>> HC 4-3-2024
 else                                                                                      !!>> HC 4-3-2024 if the morphology is broken we will discard it
    symch=666.0d0                                                                          !!>> HC 4-3-2024
 endif                                                                                     !!>> HC 4-3-2024

   inquire(file="assy_devtime.txt", exist=existes)                                         !!>> AL 17-7-2024
   if(existes)then
      open(6021, file="assy_devtime.txt", action='write',position='append')
         write(6021, *) symch
      close(6021)
   else
      open(6021, file="assy_devtime.txt", action="write")
         write(6021, *) symch
      close(6021)
   end if  

 if (symch<0.0500d0) then                                                                  !!>> HC 20-12-2020 We only want robust individuals
     dieornot=0                                                                            !!>> HC 20-12-2020 The threshold comes my observations
 else                                                                                      !!>> HC 20-12-2020  
     dieornot=1                                                                            !!>> HC 20-12-2020 Avoid extintion
 endif                                                                                     !!>> HC 20-12-2020
   
 end subroutine assymetry_save

end module conservative_R