
!!>> HC 20-11-2020  Created by Hugo Cano Fernandez 20-11-2020
!!>> HC 20-11-2020  This program creates a list where we assign mutations of different types to a fixed number of individuals
!!>> HC 20-11-2020  the arguments are the population size and the direction to write the output
!!>> HC 20-11-2020  the proportions of the different mutation types must also be specified 
!!>> HC 20-11-2020  by default we mutate the half of the population

program whomut

implicit none
integer :: N, ich, nonmutated, mutated, IS, TM      !!>> HC 20-11-2020 
integer :: counter, lucky, mutations                !!>> HC 20-11-2020 
real*8 :: a, prop_IS, prop_TM                       !!>> HC 24-11-2020 
character*200 :: direction, output, rseed           !!>> HC 20-11-2020 
character*12 :: pop                                 !!>> HC 20-11-2020
integer, allocatable, dimension(:) :: whomutates    !!>> HC 20-11-2020 
integer :: nseed
integer, allocatable, dimension (:) :: idum, newidum
character*700 :: order

call getarg(1,direction)                        !!>> HC 20-11-2020 chiddir direction to write who mutates
output=trim(direction)//"/whomut.dat"           !!>> HC 20-11-2020 path to whomutates.dat
call getarg(2,pop)                              !!>> HC 20-11-2020 number of individuals read as character
read (pop,*) N                                  !!>> HC 20-11-2020 number of individuals integer

rseed=trim(direction)//"/rseed.dat"             !!>> HC 21-11-2020 READING THE RANDOM SEED
open(386,file=trim(rseed))                      !!>> HC 21-11-2020  it was saved here by seeding.f90
    read(386,*)nseed                            !!>> HC 21-11-2020  and it is calculated from the master random seed
    if (allocated(idum)) deallocate(idum)       !!>> HC 21-11-2020  (the one in the starting i/o Emaker file)
    allocate(idum(nseed))                       !!>> HC 21-11-2020
    idum=0                                      !!>> HC 21-11-2020
    read(386,*) idum                            !!>> HC 21-11-2020
 close(386)                                     !!>> HC 21-11-2020

if (allocated(newidum)) deallocate(newidum)     !!>> HC 21-11-2020
allocate(newidum(nseed))                        !!>> HC 21-11-2020
newidum=0                                       !!>> HC 21-11-2020
    
call random_seed(put=idum)                      !!>> HC 21-11-2020 
do ich=1,nseed                                  !!>> HC 21-11-2020 RANDOM SEED FOR THE NEXT GENERATION
   call random_number(a)                        !!>> HC 21-11-2020
   newidum(ich)=int(a*1000)                     !!>> HC 21-11-2020
enddo                                           !!>> HC 21-11-2020

order="cp "//trim(rseed)//" "//trim(direction)//"/rseed_original.dat"  !!>> HC 21-11-2020 
call system(order)                              !!>> HC 21-11-2020
 open(386,file=trim(rseed))                     !!>> HC 21-11-2020 
     write(386,*) nseed                         !!>> HC 21-11-2020
     write(386,*) newidum                       !!>> HC 21-11-2020
 close(386)                                     !!>> HC 21-11-2020
 
 
 
 
if(allocated(whomutates)) deallocate(whomutates)    !!>> HC 20-11-2020 This list will store who mutates and 
allocate(whomutates(1:N))                           !!>> HC 20-11-2020 with what type of mutation
whomutates=0                                        !!>> HC 20-11-2020
nonmutated=N                                        !!>> HC 20-11-2020 

prop_IS=80d-2        !!>> HC 18-11-2021 Proportion of IS mutations in the mutated individuals
prop_TM=20d-2        !!>> HC 18-11-2021 Proportion of T mutations 

mutated=N/2       !!>> HC 20-11-2020 We want to set artificially that the half of the population mutates

IS=nint(prop_IS*mutated)      !!>> HC 24-11-2020 Number of mutated individuals with IS mutations
TM=nint(prop_TM*mutated)       !!>> HC 24-11-2020 Number of mutated individuals with T mutations 
mutations=IS+TM
!print*, "population size", N, "mutated inds", mutated, "output dir", output
!print*, "IS", IS, "TM", TM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!>> HC 20-11-2020 ASSIGN MUTATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do while(mutations>0)                      !!>> HC 20-11-2020 We are going to assign mutations to the individuals randomly
   call random_number(a)                   !!>> HC 20-11-2020
   lucky=ceiling(a*nonmutated)             !!>> HC 20-11-2020 This is the lucky one that is going to be mutated (number 1 to nonmutated)
   if (lucky==0) lucky=1                   !!>> HC 20-11-2020 in case a=0.0d0
!print*, "mutations left", mutations, "nonmutated inds", nonmutated, "going to mutate", lucky
   counter=0                               !!>> HC 20-11-2020
   do ich=1,N                              !!>> HC 20-11-2020
      if(whomutates(ich)>0) cycle          !!>> HC 20-11-2020 If the individual has to be nonmutated
      counter=counter+1                    !!>> HC 20-11-2020 This is the counter of nonmutated individuals (goes 1 to nonmutated)
      if (counter==lucky)then              !!>> HC 20-11-2020 if true, this is the lucky one!
         if (IS>0)then                     !!>> HC 20-11-2020 We assing IS mutations first
            whomutates(ich)=1              !!>> HC 20-11-2020 Assign the mutation
            IS=IS-1                        !!>> HC 20-11-2020 update the number of IS mutations left
            nonmutated=nonmutated-1        !!>> HC 20-11-2020 update the number of nonmutated individuals
            mutations=mutations-1          !!>> HC 20-11-2020 update the number of mutations left
            exit                           !!>> HC 20-11-2020 exit the loop and go to the do while
         elseif (TM>0)then                 !!>> HC 24-11-2020 We do the same for the T mutations
                whomutates(ich)=2          !!>> HC 20-11-2020
                TM=TM-1                    !!>> HC 24-11-2020
                nonmutated=nonmutated-1    !!>> HC 20-11-2020
                mutations=mutations-1      !!>> HC 20-11-2020
                exit                       !!>> HC 20-11-2020
         endif                             !!>> HC 20-11-2020
      endif                                !!>> HC 20-11-2020
   enddo                                   !!>> HC 20-11-2020
enddo                                      !!>> HC 20-11-2020

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!>> HC 20-11-2020 OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 open(1,file=output)                    !!>> HC 20-11-2020
    do ich=1,N                          !!>> HC 24-11-2020
       write(1,*) whomutates(ich)       !!>> HC 24-11-2020
    enddo                               !!>> HC 24-11-2020
 close(1)                               !!>> HC 20-11-2020



end program whomut
