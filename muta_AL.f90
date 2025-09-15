!!>> HC 2-9-2020 This programme was written by Hugo Cano. it:
! Reads and embryo and makes the mutations 
program mutar

   use io
   use mutation
   use general
   use genetic
   use neighboring
   use inicial

   implicit none
   character*140 embryo, evaparams, whomut, embryoi
   character*400 moving, rangfile                         !!>> HC 6-10-2021
   character*80 cu
   character*80 ind
   integer :: remutate, inviable, limit1, individual
   integer :: ich,mind, wrsize
   integer :: runornot                                    !!>> AL 10-4-24: if = 1 the mutation affects functional subnetwork so development of individual is runned. 
                                                         !!>> AL 10-4-24: if = 0 no development is needed, this individual inherits parent fitness
   integer, allocatable, dimension(:) :: wridum, svidum
   integer, allocatable, dimension(:) :: whomutates, effu !!>>HC 20-11-2020 Here we store what individuals are going to mutate
   integer :: mutacode, mutageni, mutagenj                !!>> HC 16-9-2021 This stores the kind of mutation that we had 
   real*8 :: newvalue, prevalue, rval                     !!>> HC 29-11-2021 Output information of mutations
   call getarg(1,embryo)
   call getarg(2,evaparams)
   call getarg(3,rangfile)   !!>> HC 6-10-2021 file with the ranges of parameters

   aut=2
   carg=embryo                              !!>> HC 2-9-2020 We read the embryo 
   call random_seed(size = nseed)                          !!>> HC 2-9-2020
   call getarg(2,cu)                                       !!>> HC 2-9-2020

   if(allocated(idum))deallocate(idum)                   !!>> HC 2-9-2020 random seed length
   allocate(idum(nseed))                                 !!>> HC 2-9-2020
   if(allocated(idumoriginal))deallocate(idumoriginal)   !!>> HC 2-9-2020
   allocate(idumoriginal(nseed))                         !!>> HC 2-9-2020
   call read_elli_evaparam(evaparams, effu)                      !!>> HC 2-9-2020 Read parameterfile
   call iniread
   call readsnap(embryo)
      
   if (allocated(svidum)) deallocate(svidum)               !!>> HC 20-12-2021 allocate
   allocate(svidum(size(idum)))                            !!>> HC 20-12-2021 
   svidum=0                                                !!>> HC 20-12-2021
   svidum=idum                                             !!>> HC 20-12-2021
      
   idum=0; wrsize=0                                        !!>> HC 20-12-2021 READ THE RANDOM SEED TO MAKE MUTATIONS
   open(837,file="rseed.dat")                              !!>> HC 20-12-2021 the seed should be stored in the individual directory
      read(837,*) wrsize                                  !!>> HC 20-12-2021 read the sixe
      if (allocated(wridum)) deallocate(wridum)           !!>> HC 20-12-2021 allocate
      allocate(wridum(wrsize))                            !!>> HC 20-12-2021 
      read(837,*) wridum                                  !!>> HC 20-12-2021 read the seed
   close(837)                                              !!>> HC 20-12-2021
   do ich=1,size(idum)                                     !!>> HC 20-12-2021 Change the seed
      if (ich>wrsize)cycle                                 !!>> HC 20-12-2021
      idum(ich)=wridum(ich)                                !!>> HC 20-12-2021
   enddo                                                   !!>> HC 20-12-2021
   call random_seed(put=idum)                              !!>> HC 20-12-2021 Put the read random seed
   call system("cp rseed.dat rseed_original.dat")          !!>> HC 20-12-2021 Copy it just in case we need to check it
   rval=0.0d0                                              !!>> HC 20-12-2021
   do ich=1,wrsize                                         !!>> HC 20-12-2021 CREATE THE RANDOM SEED FOR THE NEXT FILE
      call random_number(rval)                             !!>> HC 20-12-2021
      wridum(ich)=int(rval*1000)                           !!>> HC 20-12-2021 make a new seed
   enddo                                                   !!>> HC 20-12-2021
   open(837,file="rseed.dat")                              !!>> HC 20-12-2021 store it
      write(837,*) wrsize                                 !!>> HC 20-12-2021 Now the children of this individual will have different mutations!
      write(837,*) wridum                                 !!>> HC 20-12-2021
   close(837)                                              !!>> HC 20-12-2021

   call random_number(a)  
   if (a>0.80d0)then
      mind=2
   else
      mind=1
   endif
   !first mutate this individual                   !!>> HC 8-9-2020 Mutate this individual 
   limit1=0                                        !!>> HC 8-9-2020
   remutate=0                                      !!>> HC 8-9-2020
   do while (limit1==0 .and. remutate<10)          !!>> HC 8-9-2020
      runornot = 0                                 !!>>AL 11-4-24
      call suremuta_no_kadh_functional_net(limit1,inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue,rangfile,runornot) !!>>AL 11-4-24
      remutate=remutate+1                          !!>> HC 8-9-2020
   end do                                          !!>> HC 8-9-2020
   if(remutate==10)then
           print*,"No mutation worked... :("
   endif
   call iniio
   call writesnap                                  !!>> HC 8-9-2020

   !!!!!!!!!!!
   if(runornot == 1)then 
      open(123, file="rundevelopment.txt")          !!>>AL 11-4-24
      write(123,*) "Let's develope, pal"
      close(123)  
      !the existance of rundevelopment.txt indicates to run_all.sh that EMaker should run this individual
   end if
   !!!!!!!!!!

   ! The mutation led to an inviable individual             !!>> HC 8-9-2020  
   !!>> AL 6-9-24 Inviable means you deleted the last gene interaction in the network or have no genes kindof(3<). This is highly unprobable
   if (inviable==1) then                                    !!>> HC 8-9-2020
      print*, "INVIABLE INDIVIDUAL: ABORTION"
      open(666, file="inviable.e")
      write(666,*) "Blood and tears"                        !!>>HC 24-11-2020
      close(666)
      open(777, file="mutacode.dat")                        !!>>HC 16-9-2021
      write(777,*) "	         666         666         666   0.0000000000000000       0.0000000000000000" !!>>HC 16-9-2021
      close(777)                                            !!>>HC 16-9-2021        
   else
      open(777, file="mutacode.dat")                        !!>>HC 16-9-2021
      write(777,*) mutacode, mutageni, mutagenj, prevalue, newvalue,runornot    !!>>HC 29-11-2021
      close(777)                                            !!>>HC 16-9-2021
   endif
   
   !!>> AL 6-9-24: not sure if individual.dat0.dat and if it does, where and when... 
   moving="mv individual.dat0.dat "//trim(embryo)!//" "//trim(embryoi)   !!>>HC 20-11-2020 This substitutes the original by the mutated individual####
   call system(moving)                                     !!>>HC 20-11-2020

end program mutar
