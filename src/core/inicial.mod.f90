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
!***************  MODUL INICIAL ********************************************
!***************************************************************************
module inicial
use general
use genetic    !>>>>>>>>>>>> Is 29-4-13
use aleas
use neighboring
use ic
use io
use nexus !>>Miquel24-9-14

character*8, public :: cdate

contains

!**************************************************************************
subroutine initials
character*200 kk
  itvi=1
  itviactual=1
  call getarg(1,carg)
  if (len_trim(carg)==0.or.carg=="0") then
    !call default_param_values
    call default_ic
    erep=0.0d0
    erepcel=0.0d0
    eyou=0.0d0
    eadh=0.0d0
    etor=0.0d0
    espring=0.0d0
    !<<< Is 14-3-15
    if (allocated(nodeo)) deallocate(nodeo)
    allocate(nodeo(nda))
    print *,nda,"nda"
    nodeo=node  
    !<<< Is 14-3-15
  else
    call iniread
    call readsnap(carg)
    if(nd>1) call neighbor_build  !>>Miquel24-9-14
    !call nexe  !>>Miquel24-9-14
    if (errorlec==1.and.eva==0) stop  !>>> Is 22-1-14
  end if
  call iniio              ! this is just to allocate and inicializate the matrices for the variable names and stuff
  call inialea3d(nparti)  ! this is to inicialize the partition of random numbers in a sphere
  !call llaleat            ! this is to read a bunch of random numbers
  if(ffu(8)==0)then
    call iniboxes           ! this to inicialize the boxes
  end if
  if(nd>1) call neighbor_build
  if (ffu(24)==0)then                              !!>>HC 6-7-2021
     if (allocated(cogrid_r)) deallocate(cogrid_r) !!>>HC 17-6-2021
     allocate(cogrid_r(1:nd,1:nd)); cogrid_r=0     !!>>HC 17-6-2021
     call fill_co_griders                          !!>>HC 17-6-2021
     call restore_neighbors                        !!>>HC 17-6-2021
  endif                                            !!>>HC 6-7-2021
  call put_param_to_matrix(param)
  paramo=param


  if(aut>0)  call writesnapini    !here we write the output file with the initial conditions              !>>>>Miquel18-11-13

end subroutine

!**************************************************************************
subroutine default_ic

   call default_values          !IS 2-1-14 this is just some default values that get overrided by the ic
                                !but I put them here in case you forget
!  call epi_sphere              !OK

!  call mes_shell               !doesn't crash, but doesn't work
!  call mes_ecm                 !WTF it's not even mesenchyme, it's epithelium, DELETE?


!  call grn_test                !doesn't crash, but doesn't seem to work

!print*,"aut",aut
  if(aut/=1 .and. aut/=5)then
    !IC for the mechanisms
    print*,"ATENTION,*************"
    print*,"you didn't load any input file,"
    print*,"1 Cell death mechanism"
    print*,"2 Differential adhesion mechanism"
    print*,"3 Cell contraction mechanism"
    print*,"4 Polarized cell growth and division mechanism"
    print*,"5 Extracellular Matrix secretion"
    print*,"6 Cell migration"
    read(*,*) ic_load
  end if
!print*,"aut",aut,"icload",ic_load
!ic_load=1
  select case(ic_load)
    case(1);  call epi_apoptosis            !OK ; CHECKED 25-8-14
    case(2);  call mes_cell_sorting         !OK ; CHECKED 25-8-14
    case(3);  call invagination             !OK ; CHECKED 25-8-14
    case(4);  call epi_polar_growth_single         !OK ; CHECKED 25-8-14
    case(5);  call epi_mes_ecm         !OK ; CHECKED 25-8-14
    case(6);  call migration              !OK ; CHECKED 25-8-14
    case(7);  call mes_polar_growth                !OK ; CHECKED 25-8-14
    case(8);  call blastula_icosahedron 
    case(9);  call blastula_fibonacci
    case(10);  call father_inva
  end select


!COOL SUBROUTINES

!  call polarized                !OK ; clean ; no ; CHECKED 2-7-14
!  call epi_growth_and_division(2,2)  !OK ; clean ; on ; CHECKED 7-7-14
! PROBLEM in mitosis call epi_mes_growth_and_division  !OK ; clean ; no ; CHECKED 2-7-14
!  call epi_mes                  !OK ; clean ; no ; CHECKED 2-7-14
!  call mes_ecm                  !OK ; clean ; si ; CHECKED 7-7-14
!  call neg_adh
!  call differential_growth_epi  !OK ; clean ; si ; CHECKED 7-7-14
!  call diffusion_test           !OK ; clean ; no
!  call mesenchyme               !OK ; clean ; no
!  call mesenchyme_apop          !OK ; clean ; si ; CHECKED 7-7-14
!  call teloblast                 !                ; CHECKED 1-7-14 NO !growth does weird things, adds nodes in a single row and outside the cell
!  call differential_growth_mes  !OK ; clean ; si ; CHECKED 4-7-14
!  call mes_polar_growth         !OK ; clean ; si  ; CHECKED 4-7-14
!  call directed_mit_mes         !OK ; clean ; si ; CHECKED 7-7-14 !kinda works,supposed to do assym. div. , but the gradient it's not strong enough I guess... 
!  call mes_shell                ! miguel4-11-13
!  call ic_emt                   !OK ; clean ;si  ; CHECKED 4-7-14
!  call twoeps
!  call diffusion_test
!  call epi_active_transport
!  call mes_ecm_degradation
!  call epi_mes_primordium
!  call mes_primordium
!  call epi_mes_bud
!  call hair_placode
!  call feather_placode
!  call epi_mes_bud_ingrowth
!  call tooth_bud
!  call delta_notch

!call tooth40a
!call blastuloid

!call founding_father

!call epi_growth_and_division

!call invaginacio_diff

!call blastula

!call blastula_ensemble

!call father

end subroutine default_ic
!si si si
end module inicial

