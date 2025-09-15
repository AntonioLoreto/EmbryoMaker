

!!This is a program that makes the random seeds and random fitnesses           !!>> HC 20-12-2021
!!for all the individuals in the generation "00000" of an evolution experiment !!>> HC 20-12-2021
!!This script is called by bootstrap.sh                                        !!>> HC 20-12-2021
!!This script was written by Hugo Cano-Fernandez on 20-12-2021

program  seeds                                                   !!>> HC 20-12-2021
   use general                                                   !!>> HC 20-12-2021
   use genetic                                                   !!>> HC 20-12-2021
   use io                                                        !!>> HC 20-12-2021
   implicit none                                                 !!>> HC 20-12-2021
   integer :: ich, jch, numids, npop                             !!>> HC 20-12-2021
   integer, allocatable, dimension(:) :: childidum               !!>> HC 20-12-2021
   real :: fit, numb                                             !!>> HC 20-12-2021
   character*700 :: workdir, childir, allfather, rseed, datfit   !!>> HC 20-12-2021
   character*5 :: ind, child                                     !!>> HC 20-12-2021
   real, allocatable, dimension(:) :: fits
    
  call getarg(1,workdir)                                         !!>> HC 20-12-2021 The general working directory
  call getarg(2,allfather)                                       !!>> HC 20-12-2021 The location of the starting file
  call getarg(3,ind)                                             !!>> HC 20-12-2021 number of individuals running
  read (ind,*) numids                                            !!>> HC 20-12-2021
  call getarg(4,ind)                                             !!>> HC 20-12-2021 number of individuals running
  read (ind,*) npop                                              !!>> HC 20-12-2021
   
  call iniread                                                   !!>> HC 20-12-2021
  call readsnap(allfather)                                       !!>> HC 20-12-2021 Read the starting file
  
  if (allocated(childidum)) deallocate(childidum)                !!>> HC 20-12-2021 Allocate the vector that will store the new
  allocate(childidum(1:size(idum)))                              !!>> HC 20-12-2021 random seed
  childidum=0                                                    !!>> HC 20-12-2021
 
  if (allocated(fits)) deallocate(fits)                          !!>> HC 20-1-2022 
  allocate(fits(1:npop))                                         !!>> HC 20-1-2022 
  fits=0                                                         !!>> HC 20-1-2022

  datfit=trim(workdir)//"/population/population.datfitness"      !!>> HC 20-1-2022  path to the file that will store the random initial fitness  
  call random_number(fits)                                       !!>> HC 20-1-2022  Random initial fitness
  fits=fits/10.0d0                                               !!>> HC 20-1-2022
  open(836,file=datfit)                                          !!>> HC 20-1-2022 
      do ich=1,npop                                              !!>> HC 20-1-2022 
         write(836,*) fits(ich)                                  !!>> HC 20-1-2022 
      enddo                                                      !!>> HC 20-1-2022 
  close(836)                                                     !!>> HC 20-1-2022   
  
  fit=0.0d0; numb=0.0d0                                          !!>> HC 20-1-2022 
  
  
  do ich=1, numids                                               !!>> HC 20-12-2021 GO over all the individuals in the population
     write(child, '(I5.5)') ich                                  !!>> HC 20-12-2021 This is the number of the individual (0 to N-1)
     childir=trim(workdir)//"/running/"//trim(child)//"/"        !!>> HC 20-1-2022  This is the directory where the individual files are stored
     rseed=trim(childir)//"rseed.dat"                            !!>> HC 20-12-2021 path to the file that will store the random seed for mutations
     
     
      
     do jch=1,size(childidum)                                    !!>> HC 20-12-2021 Random seed for mutations
        call random_number(numb)                                 !!>> HC 20-12-2021 (it will be changed by muta.f90 in each generation)
        childidum(jch)=int(numb*1000)                            !!>> HC 20-12-2021
     enddo                                                       !!>> HC 20-12-2021
     
     open(836,file=rseed)                                        !!>> HC 20-12-2021 write the random seed
        write(836,*) size(childidum)                             !!>> HC 20-12-2021
        write(836,*) childidum                                   !!>> HC 20-12-2021
     close(836)                                                  !!>> HC 20-12-2021
  enddo                                                          !!>> HC 20-12-2021
  
  rseed=trim(workdir)//"/finished/rseed.dat"                     !!>> HC 21-12-2021 
  
  do jch=1,size(childidum)                                       !!>> HC 21-12-2021 Random seed for deciding who mutates in each generation
     call random_number(numb)                                    !!>> HC 21-12-2021 (it will be changed by substitute.f90 in each generation)
     childidum(jch)=int(numb*1000)                               !!>> HC 21-12-2021
  enddo                                                          !!>> HC 21-12-2021
     
  open(836,file=rseed)                                           !!>> HC 21-12-2021 write the random seed
      write(836,*) size(childidum)                               !!>> HC 21-12-2021
      write(836,*) childidum                                     !!>> HC 21-12-2021
  close(836)                                                     !!>> HC 21-12-2021
  

end program                                                      !!>> HC 20-12-2021
