!!This is a program that makes the random seeds and random fitnesses           !!>> HC 20-12-2021
!!for all the individuals in the generation "00000" of an evolution experiment !!>> HC 20-12-2021
!!This script is called by bootstrap.sh                                        !!>> HC 20-12-2021
!!This script was written by Hugo Cano-Fernandez on 20-12-2021
program  seeds                                                   !!>> HC 20-12-2021
   use general                                                   !!>> HC 20-12-2021
   use genetic                                                   !!>> HC 20-12-2021
   use io                                                        !!>> HC 20-12-2021
  
   implicit none                                                 !!>> HC 20-12-2021
   
   character*700  :: workdir, childir, allfather, rseed, datfit  !!>> HC 20-12-2021
   character*500  :: pathrandseed,rsid                           !!>> AL 5-4-2024
   character*5    :: ind,time
   character*6    :: child
   integer        :: ich,jch,numids,numids2,npop,i_cont       !>> HC 20-12-2021
   integer        :: randseed,enene,clock                                    !!>> AL 5-4-2024
   real           :: fit, numb                                   !!>> HC 20-12-2021
   integer, allocatable, dimension(:)  :: childidum              !!>> HC 20-12-2021
   real, dimension(:), allocatable     :: fits
   integer, dimension(:), allocatable  :: seed_evo, seed_subs    !!>> AL 1-4-2024

   call getarg(1,workdir)                                         !!>> HC 20-12-2021 The general working directory
   call getarg(2,allfather)                                       !!>> HC 20-12-2021 The location of the starting file
   call getarg(3,ind)                                             !!>> HC 20-12-2021 number of individuals running (cpus)
   read (ind,*) numids                                            !!>> HC 20-12-2021
   call getarg(4,ind)                                             !!>> HC 20-12-2021 number of individuals running
   read (ind,*) npop                                              !!>> HC 20-12-2021
   
   call getarg(5,rsid)                                            !!>> AL 1-4-2024 
   read (rsid,*) randseed                                         !!>> AL 1-4-2024
   if(randseed == 1) then                                         !!>> AL 1-4-2024
      print*, 'pathrandseed='
      call getarg(6,pathrandseed)                                 !!>> AL 1-4-2024
      print*, pathrandseed    
   else 
      print*,'randseed=',randseed
      print*,'random seed will be initialize to system clock...'
   end if                                                         !!>> AL 1-4-2024

   call iniread                                                   !!>> HC 20-12-2021
   call readsnap(allfather)                                       !!>> HC 20-12-2021 Read the starting file !!>>AL 3-4-24 here in theory you initiliaze random_seed() with random seed in initial condicionts file
   
!!!This part was written by Antonio Loreto 1-4-2024. It produces a random_seed array using
!!!the system time and it saves such array in a text file. Using this specific random_seed array
!!!enables you to repeat the sequence of all random_numbers generated during this evolution experiment->!!3-4-24  Not sure if this apply to EMaker development process..
if (randseed == 0) then
   call random_seed(size=enene)
   if(allocated(seed_evo)) deallocate(seed_evo)                
   allocate(seed_evo(enene))
   call system_clock(count=clock)
   seed_evo = clock + 31 * (/ (i_cont - 1, i_cont = 1, enene) /) 
   call random_seed(put=seed_evo)
   !write(time,*) clock
   rseed=trim(workdir)//"/rseed_master_evo_exp.dat"                  
   open(836,file=rseed)                                        
      write(836,*) numids        !!evolution experiments are repeatable if you use the same number of cpus (cpus=numids) !!AL 28-2-25 this is not true in a HPC
      write(836,*) size(seed_evo)                             
      write(836,*) seed_evo                                  
   close(836)                                                                                                        
   deallocate(seed_evo)
else 
   open(333,file=trim(pathrandseed))                                        
      read(333,*) numids2        !!evolution experiments are repeatable if you use the same number of cpus (cpus=numids)
      read(333,*) enene                             
      if(allocated(seed_evo)) deallocate(seed_evo)                
      allocate(seed_evo(enene))
      read(333,*) seed_evo                                  
   close(333)  
   print*,'We will use the following random seed:'
   print*,numids2
   print*,enene
   print*,seed_evo
   call random_seed(put=seed_evo)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (allocated(childidum)) deallocate(childidum)                !!>> HC 20-12-2021 Allocate the vector that will store the new
   allocate(childidum(enene))                                     !!>> HC 20-12-2021 random seed 
   childidum=0                                                    !!>> HC 20-12-2021
   
!    if (allocated(fits)) deallocate(fits)                          !!>> HC 20-1-2022 
!    allocate(fits(1:npop))                                         !!>> HC 20-1-2022 
!    fits=0                                                         !!>> HC 20-1-2022

!   datfit=trim(workdir)//"/population/population.datfitness"      !!>> HC 20-1-2022  path to the file that will store the random initial fitness  
!   call random_number(fits)                                       !!>> HC 20-1-2022  Random initial fitness !!AL 26-11-24 I dont get why random initial fintess (erease maybe)
!   !fits=fits/10.0d0                                              !!>> HC 20-1-2022 
!   open(836,file=datfit)                                          !!>> HC 20-1-2022 
!       do ich=1,npop                                              !!>> HC 20-1-2022 
!          write(836,*) 1001                                       !!>> HC 20-1-2022 
!       enddo                                                      !!>> HC 20-1-2022 
!   close(836)                                                     !!>> HC 20-1-2022   
  
  fit=0.0d0; numb=0.0d0                                          !!>> HC 20-1-2022 

  do ich=1, numids                                               !!>> HC 20-12-2021 GO over all the individuals in the population
     write(child, '(I6.6)') ich                                  !!>> HC 20-12-2021 This is the number of the individual (0 to N-1)
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
  end do                                                          !!>> HC 20-12-2021

  rseed=trim(workdir)//"/finished/rseed.dat"                     !!>> HC 21-12-2021 
  do jch=1,size(childidum)                                       !!>> HC 21-12-2021 Random seed for deciding who mutates in each generation
     call random_number(numb)                                    !!>> HC 21-12-2021 (it will be changed by substitute.f90 in each generation)
     childidum(jch)=int(numb*1000)                               !!>> HC 21-12-2021
  enddo                                                          !!>> HC 21-12-2021
  open(836,file=rseed)                                           !!>> HC 21-12-2021 write the random seed
      write(836,*) size(childidum)                               !!>> HC 21-12-2021
      write(836,*) childidum                                     !!>> HC 21-12-2021
  close(836)                                                     !!>> HC 21-12-2021

  if (allocated(seed_subs)) deallocate(seed_subs)                !!>> HC 20-12-2021 Allocate the vector that will store the new
  allocate(seed_subs(enene))                                     !!>> HC 20-12-2021 random seed 
  seed_subs=0                                                    !!>> HC 20-12-2021

  do jch=1,enene                                                 !!>> AL 5-4-24  
     call random_number(numb)                                    !!>> AL 5-4-24  
     seed_subs(jch)=int(numb*1000)                               !!>> AL 5-4-24 
  enddo                                                          !!>> AL 5-4-24 
  rseed=trim(workdir)//"/rseedsubs.dat"                          !!>> AL 5-4-24 random seed to be used in cpu that makes substitutions
  open(836,file=rseed)                                           !!>> AL 5-4-24  write the random seed
      write(836,*) size(seed_subs)                               !!>> AL 5-4-24 
      write(836,*) seed_subs                                     !!>> AL 5-4-24 
  close(836)

end program                                                      
