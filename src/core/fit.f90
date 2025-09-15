

program fitness

  use general
  use genetic
  use neighboring
  use io
  use analysis
  use conservative_R ! module created by Renske, contains simple EMD and PF H measures of similarity to ancestor
  use fitmo
  use mutation


  implicit none
  character*400 onetwin, parfile, rangfile !!>> HC 6-10-2021
  character*300 cx
  real*8 fitelli, prevfit, emd  !!>> HC 30-6-2023
  integer :: ich
  integer, allocatable, dimension(:) :: effu
  call getarg(1,onetwin)
  call getarg(2,parfile)
  call getarg(3,rangfile) !!>> HC 6-10-2021
  call iniread
  call readsnap(onetwin)           !!>> HC 17-11-2020 we read the individual
  call neighbor_build              !!>> HC 17-11-2020
  call read_elli_evaparam(parfile, effu) !!>> HC 17-11-2020 we read the parfile to know whether there is a target or not
  
  
 
  if (len_trim(tarfitmorphfile)>1)then                       !!>> HC 17-11-2020  If we have stated some Target file
      call fit(0,fitelli,emd,onetwin,tarfitmorphfile)        !!>> HC 30-6-2023   we calculate the EMD distance to that target
      open(202,file="individual.datfitness")                 !!>> HC 17-11-2020 
      write(202,*)fitelli, emd                               !!>> HC 30-6-2023 
      close(202)                                             !!>> HC 17-11-2020 
      open(23,file=trim(carg)//"t",iostat=i)                 !!>> HC 17-11-2020 
      write(23,*,ERR=46) trim(carg)//trim(nofi)              !!>> HC 17-11-2020 
      open(24,file=trim(carg)//"x",iostat=i)                 !!>> HC 17-11-2020 
      write(24,*,ERR=46) 0                                   !!>> HC 17-11-2020 
  
      flush(24)                                              !!>> HC 17-11-2020 
      flush(23)                                              !!>> HC 17-11-2020 
 
      close(23)  ! The read statement further down gives an error without this! >>>Renske 05-03-18 !!>> HC 17-11-2020 
      open(23,file=trim(carg)//"t",iostat=ii)                !!>> HC 17-11-2020 

      read(23,*,END=45,ERR=48) cx                            !!>> HC 17-11-2020 
      call flush(23)                                         !!>> HC 17-11-2020 
     ! print *,"done with iostat",ii                          !!>> HC 17-11-2020 
     ! print *, "almost end of auto"                          !!>> HC 17-11-2020 
      call exit(231)                                         !!>> HC 17-11-2020 
    45 print *,"end of file error 111"                       !!>> HC 17-11-2020 
      call exit(231)                                         !!>> HC 17-11-2020 
    46 print *,"error in writing",trim(carg)//"t"            !!>> HC 17-11-2020 
      call exit(231)                                         !!>> HC 17-11-2020 
    48 print *,"other end 222"                               !!>> HC 17-11-2020 
      call exit(231)                                         !!>> HC 17-11-2020 
     
  else                                                                          !!>> HC 11-11-2021 If there is no target file we check which effu is == 1
      if (effu(1)==1)then; call fill_node_arrays; call OPCval; endif            !!>> HC 11-11-2021 OPC
      if (effu(2)==1) call joint_entropy_fourth_original_neighbors(fitelli)     !!>> HC 11-11-2021 GLOBAL COMPLEXITY based on curvature
      if (effu(3)==1)then                                                       !!>> HC 21-4-2022  LOCAL COMPLEXITY  based on curvature
         call local_joint_entropy_second(fitelli)                               !!>> HC 21-4-2022
         prevfit=fitelli                                                        !!>> HC 21-4-2022
         call local_joint_entropy_fourth(fitelli)                               !!>> HC 21-4-2022
         open(206,file="complexity_vals.dat")                                   !!>> HC 21-4-2022
             write(206,*) fitelli, prevfit                                      !!>> HC 21-4-2022
         close(206)                                                             !!>> HC 21-4-2022
         fitelli=fitelli+prevfit                                                !!>> HC 21-4-2022
      endif                                                                     !!>> HC 21-4-2022
      if (effu(4)==1)then                                                       !!>> HC 11-11-2021  COMBINED COMPLEXITY (global and local) based on curvature
         call joint_entropy_fourth_original_neighbors(fitelli)                  !!>> HC 11-11-2021
         prevfit=fitelli                                                        !!>> HC 11-11-2021
         call local_joint_entropy_second(fitelli)                                      !!>> HC 11-11-2021
         open(206,file="complexity_vals.dat")                                   !!>> HC 11-11-2021
             write(206,*) fitelli, prevfit                                      !!>> HC 11-11-2021
         close(206)                                                             !!>> HC 11-11-2021
         fitelli=fitelli*prevfit                                                !!>> HC 11-11-2021
      endif                                                                     !!>> HC 11-11-2021
      if (effu(5)==1) call flux_resistance_decrease(fitelli)                    !!>> HC 11-11-2021 Decreased flux resistance
      if (effu(6)==1) call flux_resistance_increase(fitelli)                    !!>> HC 11-11-2021 Increase flux resistance
      if (effu(7)==1) call surface_volume_ratio(fitelli)                        !!>> HC 11-11-2021 Surface volume (untested)
      if (effu(8)==1) call flux_half_increase_half_decrease(fitelli)            !!>> HC 12-5-2022  half decrease halv increase flux resistance (untested)
      if (effu(9)==1) call curvature_patch_count(fitelli)                       !!>> HC 12-5-2022  Patches of cells with the same curvature (untested)
      if (effu(10)==1) call directional_selection(fitelli)                      !!>> HC 18-2-2024  Directional selection (works)
      
          !!!OUTPUT!!!
      if (effu(1)==0)then
         open(202,file="individual.datfitness")                              !!>> HC 11-11-2021 WRITING OUTPUT
         write(202,*) fitelli                                                !!>> HC 11-11-2021
         close(202)                                                          !!>> HC 11-11-2021
         open(203,file="OPC.val")                                            !!>> HC 11-11-2021
         write(203,*) fitelli                                                !!>> HC 11-11-2021
         close(203)   
      endif
  endif



end program fitness
