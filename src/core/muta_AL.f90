!gfortran -w -g -fexceptions -fno-underscoring -fcheck=all  general.mod.f90 genetic.mod.f90 io.mod.f90 test-mutation_AL.mod.f90  test-muta_AL.f90 -o t-muta_AL.e 
program mutar

   use io
   use mutation
   use general
   use genetic
   !use neighboring
   !use inicial

   implicit none

   character*140 embryo, evaparams, whomut, embryoi
   character*400 moving, rangfile,popdatfitness,paste                         !!>> HC 6-10-2021
   character*80 cu
   character*80 ind
   integer :: remutate, inviable, limit1, individual
   integer :: ich,jch,mind,wrsize,howmanymut
   integer :: runornot                                    !!>> AL 10-4-24: if = 1 the mutation affects functional subnetwork so development of individual is runned. 
                                                          !!>> AL 10-4-24: if = 0 no development is needed, this individual inherits parent fitness
   integer, allocatable, dimension(:) :: wridum, svidum
   integer, allocatable, dimension(:) :: whomutates, effu      !!>>HC 20-11-2020 Here we store what individuals are going to mutate
   integer :: mutacode, mutageni, mutagenj, nmuta    
   real*8  :: newvalue, prevalue, rval, pcpval,p_m            !!>> HC 29-11-2021 Output information of mutations
   logical :: exist

   call getarg(1,embryo)
   call getarg(2,evaparams)
   call getarg(3,rangfile)   !!>> HC 6-10-2021 file with the ranges of parameters

   aut    = 2                                              !!>> AL 26-11-24 no idea what this does 
   pcpval = 0.0d0
   carg   = embryo                                         !!>> HC 2-9-2020 We read the embryo !!>> AL 26-11-24 no idea what this does 
   p_m = 0.1                                               !!>> AL 12-2-25: probability of gene mutation 

   call random_seed(size = nseed)                          !!>> HC 2-9-2020 !!>> AL 26-11-24 nseed is initialized in general 
   
   !call getarg(2,cu)                                      !!>> HC 2-9-2020

   if(allocated(idum))deallocate(idum)                     !!>> HC 2-9-2020 random seed length
   allocate(idum(nseed))                                   !!>> HC 2-9-2020

   if(allocated(idumoriginal))deallocate(idumoriginal)     !!>> HC 2-9-2020
   allocate(idumoriginal(nseed))                           !!>> HC 2-9-2020

   call read_elli_evaparam(evaparams, effu)                !!>> HC 2-9-2020 Read parameterfile
   call iniread
   call readsnap(embryo)

   if (allocated(svidum)) deallocate(svidum)               !!>> HC 20-12-2021 allocate
   allocate(svidum(size(idum)))                            !!>> HC 20-12-2021 
   svidum=0                                                !!>> HC 20-12-2021
   svidum=idum                                             !!>> HC 20-12-2021

   idum=0; wrsize=0                                        !!>> HC 20-12-2021 READ THE RANDOM SEED TO MAKE MUTATIONS
   open(837,file="rseed.dat")                              !!>> HC 20-12-2021 the seed should be stored in the individual directory
      read(837,*) wrsize                                   !!>> HC 20-12-2021 read the sixe
      if (allocated(wridum)) deallocate(wridum)            !!>> HC 20-12-2021 allocate
      allocate(wridum(wrsize))                             !!>> HC 20-12-2021 
      read(837,*) wridum                                   !!>> HC 20-12-2021 read the seed
   close(837)                                              !!>> HC 20-12-2021
   do ich=1,size(idum)                                     !!>> HC 20-12-2021 Change the seed
      if (ich>wrsize)cycle                                 !!>> HC 20-12-2021
      idum(ich)=wridum(ich)                                !!>> HC 20-12-2021 !!>> AL 26-11-24 i dont see the point in modifying public variable idum...
   enddo                                                   !!>> HC 20-12-2021
   call random_seed(put=idum)                              !!>> HC 20-12-2021 Put the read random seed
   call system("cp rseed.dat rseed_original.dat")          !!>> HC 20-12-2021 Copy it just in case we need to check it !!>> AL 26-11-24: this might have been for a test.
   rval=0.0d0                                              !!>> HC 20-12-2021
   do ich=1,wrsize                                         !!>> HC 20-12-2021 CREATE THE RANDOM SEED FOR THE NEXT FILE
      call random_number(rval)                             !!>> HC 20-12-2021
      wridum(ich)=int(rval*1000)                           !!>> HC 20-12-2021 make a new seed
   enddo                                                   !!>> HC 20-12-2021
   open(837,file="rseed.dat")                              !!>> HC 20-12-2021 store it
      write(837,*) wrsize                                  !!>> HC 20-12-2021 Now the children of this individual will have different mutations! !!>> AL 27-11-24 in substitute we rewrite rseed so this is duplicated
      write(837,*) wridum                                  !!>> HC 20-12-2021
   close(837)                                              !!>> HC 20-12-2021
   
   howmanymut = 0
   runornot = 0                                            !!>>AL 25-2-25

   do ich = 1,ng                                           !!>>AL 14-2-25 each gene has its own probability of mutation
      
      if (int(gen(ich)%kindof) == 9) cycle 
      
      nmuta = 0
      call random_number(a)
      if (a < p_m) then  

         nmuta = 1
         call random_number(a)
         if (a < p_m**2) then
            nmuta = 2
         end if 

         howmanymut = howmanymut + nmuta

         do jch = 1,nmuta

            call random_number(a)                          !!>>AL 30-5-25 we force IS muta for fase II
            if (a > 0.80d0)then                             !!>>AL 14-2-25: choose type of mutation (T or IS)
               mind=2
            else
               mind=1
            end if
            
            mutacode = 0
            mutageni = 0
            mutagenj = 0 
            prevalue = 0
            newvalue = 0
            pcpval   = -1

            limit1   = 0                                         !!>> HC 8-9-2020
            remutate = 0                                         !!>> HC 8-9-2020
            do while (limit1 == 0 .and. remutate < 100)          !!>> HC 8-9-2020
               call suremuta_no_kadh_functional_net_rec_3(ich,limit1,inviable,mind,mutacode,mutageni,mutagenj,prevalue,newvalue &
               ,rangfile,runornot,pcpval) !!>>AL 11-4-24
               remutate=remutate+1                           !!>> HC 8-9-2020
            end do                                           !!>> HC 8-9-2020
            
            if(remutate == 100)then
               open(123, file="nomutacode.txt")              !!>>AL 11-4-24 if this file exists we run development 
                  write(123,*) "No mutation worked... :("
               close(123)  
               inviable = 1
            end if

            inquire(file = "mutacodetmp.dat", exist = exist)              !!>> AL 27-1-25            
            if (inviable == 1) then                                    
               open(666, file="inviable.e")
                  write(666,*) "Blood and tears"                        
               close(666)
               open(777, file="mutacodetmp.dat")                        
                  write(777,*) " 666 666 666 0 0 0 0 0" 
               close(777)                                                 
            else
               if (exist) then
                  !print*,ich
                  open(777, file="mutacodetmp.dat", position='append')                        !!>>HC 16-9-2021
                     write(777,*) ich,mutacode, mutageni, mutagenj, prevalue, newvalue,runornot,pcpval     !>>HC 29-11-2021
                  close(777)  
               else
                  open(777, file="mutacodetmp.dat")                        !!>>HC 16-9-2021
                     write(777,*) ich,mutacode, mutageni, mutagenj, prevalue, newvalue,runornot,pcpval     !>>HC 29-11-2021
                  close(777)  
               end if  
            end if 

         end do
      end if
   end do 

   if(howmanymut == 0)then 
      open(777, file="mutacodetmp.dat")                        !!>>AL 27-2-25
         write(777,*) "0 0 0 0 0 0 0 0" 
      close(777)                               
   end if 
   
   open(777, file="mutacode.dat", position='append')            !!>> AL 25-2-25       
      write(777,*) howmanymut    
   close(777)  

   paste = "paste "//"mutacodetmp.dat"//" >> "//"mutacode.dat"  !!>> AL 25-2-25  
   call system(trim(paste))    

   paste = "rm "//"mutacodetmp.dat"                             !!>> AL 25-2-25  
   call system(trim(paste))

   call iniio
   call writesnap                                               !!>> HC 8-9-2020

   inquire(file="rundevelopment.txt", exist=exist)              !!>> AL 27-1-25     !25-2-25: not sure when it is possible that this file already exists... check!       
   if(exist)then
      runornot = 1
   end if

   if(runornot == 1)then 
      open(123, file="rundevelopment.txt")                      !!>>AL 11-4-24 if this file exists we run development 
         write(123,*) "Let's develope, pal"
      close(123)  
   end if 

   !!>> AL 6-9-24: not sure if individual.dat0.dat and if it does, where and when... 
   moving="mv individual.dat0.dat "//trim(embryo)!//" "//trim(embryoi)   !!>>HC 20-11-2020 This substitutes the original by the mutated individual####
   !moving="mv icoepi_nuclei_24G_rec.dat0.dat "//trim(embryo)!//" "//trim(embryoi)   !! this option only for individual testing
   call system(moving)                                     !!>>HC 20-11-2020

end program mutar
