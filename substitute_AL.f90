program substitute
  implicit none
  integer :: i, j, nsubs, nmax, npop, theone, todie, thechosen, nseed,nminfits,largo
  character*700 :: bindir, population, num, finish, finished, counter, running, fitnesses, popsize,assnfit,fitnass
  character*700 :: directory, indfile, datfitness, allfin, newid, besty, best, datfin, parent, time, popids
  character*700 :: output, outind, robfile, mutacode, idfile, start, rseed, oldseed,outfile,voll,rseedsubs
  character*1400 :: cp, rm, cat
  character*10000 :: paste 
  character*5 :: id
  real*8, allocatable, dimension(:) :: fits
  real*8 :: indfit, minfit, maxfit, a, fitn, assn
  integer, allocatable, dimension(:) :: idum, newidum, popid
  logical :: exist
  
  call getarg(1,bindir)                                    !!>> HC 21-11-2020
  call getarg(2,finished)                                  !!>> HC 21-11-2020 This is the file that has been produced in the directory finished indicating
  open(1,file=trim(finished))                              !!>> HC 21-11-2020 that one CPU has finished running EMaker
      read(1,*,END=666) theone                             !!>> HC 21-11-2020 the one that has finished
  close(1)                                                 !!>> HC 21-11-2020
  write(id,'(I5.5)') theone                                !!>> HC 21-11-2020

  finish=trim(bindir)//"/finished"                         !!>> HC 21-11-2020 path to the directory where we save the CPUs that have finished
  population=trim(bindir)//"/population"                   !!>> HC 21-11-2020 path to the directory where we have the population files
  running=trim(bindir)//"/running"                         !!>> HC 21-11-2020 path to the directory where we have the CPUs running
  output=trim(bindir)//"/output_model.dat"                 !!>> HC 21-11-2020 path to the file where we save the general output of the evolutionary model
  rseedsubs=trim(bindir)//"/rseedsubs.dat"                 !!>> AL 5-4-24
  fitnesses=trim(population)//"/population.datfitness"     !!>> HC 21-11-2020 path to the file where we save the population fitnesses 
  outfile=trim(population)//"transout.dat"                 !!>> HC 21-11-2020 transient file
  popsize=trim(population)//"/npop.dat"                    !!>> HC 21-11-2020 path to the file where we save the number of individuals in the population
  popids=trim(population)//"/popids.dat"                   !!>> HC 21-11-2020 path to the file where we save the ids of the individuals in the population

  fitnass=trim(population)//"/fitnassy.txt"                !!>> AL 15-4-2024

  counter=trim(finish)//"/counter.dat"                     !!>> HC 21-11-2020 path to the file where we save the number of substitutions we have run
  
  directory=trim(running)//"/"//trim(id)                   !!>> HC 21-11-2020 path to the individual that has finished
  indfile=trim(directory)//"/individual.dat"               !!>> HC 21-11-2020 path to the individual's ic file
  datfitness=trim(directory)//"/individual.datfitness"     !!>> HC 21-11-2020 path to the individual's fitness
  datfin=trim(directory)//"/individual.dat*dat"            !!>> HC 21-11-2020 path to the output file (morphology)
  parent=trim(directory)//"/parent.dat"                    !!>> HC 21-11-2020 path to the individual's parent index
  robfile=trim(directory)//"/rob.val"                      !!>> HC 21-11-2020 path to the individual's robustness
  mutacode=trim(directory)//"/mutacode.dat"                !!>> HC 21-11-2020 path to the individual's code of mutation
  outind=trim(directory)//"/outind.dat"                    !!>> HC 21-11-2020 transient file
  idfile=trim(directory)//"/idfile.dat"                    !!>> HC 21-11-2020 path to the new ID of this individual
  start=trim(directory)//"/start.dat"                      !!>> HC 21-11-2020 path to the file that will allow EMaker to run again
  time=trim(directory)//"/development_time.dat"                        !!>> HC 21-11-2020 path to the individual's runtime
  rseed=trim(directory)//"/rseed.dat"                      !!>> HC 21-11-2020 path to the random seed
  oldseed=trim(directory)//"/rseed_original.dat"           !!>> HC 21-11-2020 path to the original random seed !!AL>> this actually is pointless

  voll=trim(directory)//"/individual.volume.txt"           !!>>AL 15-4-2024
  assnfit=trim(directory)//"/assynfit.txt"                 !!>>AL 15-4-2024

   !***************************************************************************!
   !!>>AL 1-4-2024 The following allows to repeat the sequence of substitutions 
   open(386,file=trim(rseedsubs))                         !!>>AL 1-4-2024
      read(386,*) nseed                                   !!>>AL 1-4-2024
      if (allocated(idum)) deallocate(idum)               !!>>AL 1-4-2024 
      allocate(idum(nseed))                               !!>>AL 1-4-2024 
      idum=0                                              !!>>AL 1-4-2024 
      read(386,*) idum                                    !!>>AL 1-4-2024 
   close(386)                                             !!>>AL 1-4-2024 
   call random_seed(put=idum)

   call random_seed(size=largo)
   if (allocated(idum)) deallocate(idum)                  !!>> HC 20-12-2021 Allocate the vector that will store the new
   allocate(idum(largo))                                  !!>> HC 20-12-2021 random seed for next sub
   idum=0     

   do i=1,size(idum)                                      !!>> AL 5-4-24  
      call random_number(a)                               !!>> AL 5-4-24  
      idum(i)=int(a*1000)                                 !!>> AL 5-4-24 
   end do          

   open(836,file=rseedsubs)                               !!>> AL 5-4-24  write the random seed
      write(836,*) size(idum)                             !!>> AL 5-4-24 
      write(836,*) idum                                   !!>> AL 5-4-24 
   close(836)                      
   !***************************************************************************! 

   inquire(file=trim(datfitness), exist=exist)       !!>> AL 6-9-24            
   if(exist)then
     open(6021, file=trim(datfitness))
       read(6021, *) indfit
     close(6021)
     print *,"[substitute_AL.sh]: this is fitness inside individual.datfitness ", indfit
   else
      print *,"[substitute_AL.f90]: error, no fitness file"
      open(1,file=trim(datfitness))
              write(1,*) 0
      close(1)
   end if

  !ios = 0                                                  !!>> AL 6-9-24 
  !open(1,file=trim(datfitness))                            
  !    read(1,*,iostat=ios) indfit                          
  !    if(ios .ne. 0)then
  !       indfit = 0
  !    endif
  !close(1)                                                 
  !if(ios .ne. 0)then
  !      print *,"[substitute_AL.f90]: error, no fitness file"
  !      open(1,file=trim(datfitness))
  !              write(1,*) 0
  !      close(1)
  !endif
 
  call getarg(3,num)                                       !!>> HC 21-11-2020
  read (num,*) nmax                                        !!>> HC 21-11-2020 maximum number of substitutions
  
  open(1,file=trim(popsize))                               !!>> HC 21-11-2020
        read(1,*) npop                                     !!>> HC 21-11-2020 number of individuals
  close(1)                                                 !!>> HC 21-11-2020

  open(1,file=trim(counter))                               !!>> HC 21-11-2020
      read(1,*) nsubs                                      !!>> HC 21-11-2020 number of substitutions we have already done
  close(1)                                                 !!>> HC 21-11-2020

  if(allocated(fits)) deallocate(fits)                     !!>> HC 21-11-2020
  allocate(fits(npop)); fits=0.0d0                         !!>> HC 21-11-2020
  open(1,file=trim(fitnesses))                             !!>> HC 21-11-2020
     do i=1,npop                                           !!>> HC 21-11-2020
        read(1,*) fits(i)                                  !!>> HC 21-11-2020 fitnesses in the present population
     enddo                                                 !!>> HC 21-11-2020
  close(1)                                                 !!>> HC 21-11-2020

  if(allocated(popid)) deallocate(popid)                   !!>> HC 21-11-2020
  allocate(popid(npop)); popid=0                           !!>> HC 21-11-2020
  open(1,file=trim(popids))                                !!>> HC 21-11-2020
     do i=1,npop                                           !!>> HC 21-11-2020
        read(1,*) popid(i)                                 !!>> HC 21-11-2020 population ids in the present population (aka which was your substitution number)
     enddo                                                 !!>> HC 21-11-2020
  close(1)                                                 !!>> HC 21-11-2020

  nsubs=nsubs+1                                            !!>> HC 21-11-2020 Number of this substitution

  if (nsubs.le.nmax)then                                   !!>> HC 21-11-2020 IF there are still substitutions to do      
      open(1,file=trim(counter))                            !!>> HC 21-11-2020 write the new number of done substitutions
         write(1,*) nsubs                                  !!>> HC 21-11-2020
      close(1)                                              !!>> HC 21-11-2020

      minfit=minval(fits)                                   !!>> HC 21-11-2020 minimum fitness value in the population
      maxfit=maxval(fits)                                   !!>> HC 21-11-2020
      if (indfit.ge.minfit)then                             !!>> HC 21-11-2020 This individual has a fitness higher or equal than the minimum
         nminfits=0                                         !!>> HC 21-11-2020 number of individuals with the minimum fitness
         do i=1,npop                                        !!>> HC 21-11-2020
           if(fits(i)>minfit)cycle                         !!>> HC 21-11-2020
           nminfits=nminfits+1                             !!>> HC 21-11-2020
        enddo                                              !!>> HC 21-11-2020
        
        call random_number(a)                              !!>> HC 21-11-2020 This needs to be instantiated
        nminfits=nint(nminfits*a)                          !!>> HC 21-11-2020 The chosen individual with the smallest fitness
        if(nminfits==0) nminfits=1                         !!>> HC 21-11-2020
        j=0                                                !!>> HC 21-11-2020
        do i=1,npop
           if(fits(i)>minfit)cycle                         !!>> HC 4-3-2024 Go through the individuals that have the smallest fitness value
           j=j+1                                           !!>> HC 4-3-2024 the number of those
           if(nminfits==j)then                             !!>> HC 4-3-2024 this is the randomly chosen one
             todie=i                                       !!>> HC 4-3-2024 id of the individual with the minimum fitness value 
             fits(i)=indfit                                !!>> HC 4-3-2024 substitute in the fitness list
             popid(i)=nsubs                                !!>> HC 4-3-2024 The id of this new individual will be the number of this substitution
             exit                                          !!>> HC 4-3-2024
           endif                                           !!>> HC 4-3-2024
        enddo                                              !!>> HC 4-3-2024
        write(id,'(I5.5)') todie                                  !!>> HC 21-11-2020 THE INDIVIDUAL TO BE SUSTITUTED BY
        newid=trim(population)//"/individual"//trim(id)//".dat"   !!>> HC 21-11-2020
        cp="cp "//trim(indfile)//" "//trim(newid)                 !!>> HC 21-11-2020
        call system(trim(cp))                                     !!>> HC 21-11-2020
        cp='sed -i "'//trim(id)//'s/.*/$(cat '//trim(assnfit)//')/" '//trim(fitnass) !!>>AL 15-4-2024: replace assymetry and fitness in population file (fitnass) for the new population integrant 
        call system(trim(cp))
        open(1,file=trim(idfile))                           !!>> HC 21-11-2020
           write(1,*) todie                                 !!>> HC 21-11-2020
        close(1)                                            !!>> HC 21-11-2020
        
        open(1,file=trim(fitnesses))                        !!>> HC 21-11-2020
            do i=1,npop                                     !!>> HC 21-11-2020
               write(1,*) fits(i)                           !!>> HC 21-11-2020 new fitnesses list
            enddo                                           !!>> HC 21-11-2020
        close(1)                                            !!>> HC 21-11-2020
        
        open(1,file=trim(popids))                           !!>> HC 21-11-2020
            do i=1,npop                                     !!>> HC 21-11-2020
               write(1,*) popid(i)                          !!>> HC 21-11-2020 new individual ids
            enddo                                           !!>> HC 21-11-2020
        close(1)                                            !!>> HC 21-11-2020
       
        if (indfit>maxfit)then                                            !!>> HC 21-11-2020 This individual has a maximum fitness value
           best=trim(bindir)//"/best"                                     !!>> HC 21-11-2020
           write(id,'(I5.5)') nsubs                                       !!>> HC 21-11-2020
           besty=trim(bindir)//"/best/individual"//trim(id)//".dat"       !!>> HC 21-11-2020
           cp="cp "//trim(datfin)//" "//trim(besty)                       !!>> HC 21-11-2020
           call system(trim(cp))                                          !!>> HC 21-11-2020 save it in the "best" directory
        endif                                                             !!>> HC 21-11-2020
     else                                                                 !!>> HC 21-11-2020 THIS INDIVIDUAL IS NOT BETTER THAN THE WORST IN THE POPULATION
        open(1,file=trim(idfile))                                         !!>> HC 21-11-2020
           write(1,*) "NA"                                                !!>> HC 21-11-2020 NA means this individual has been discarded
        close(1)                                                          !!>> HC 21-11-2020
     endif                                                                !!>> HC 21-11-2020
     
     paste="paste "//trim(counter)//" "//trim(idfile)//" "//trim(parent)//" "//trim(datfitness)// &     !!>> HC 21-11-2020
     " "//trim(robfile)//" "//trim(voll)//" "//trim(mutacode)//" "//trim(time)//" > "//trim(outind)     !!>> HC 21-11-2020
     call system(trim(paste))                                                                           !!>> HC 21-11-2020 output variables of this individual
     cat="cat "//trim(outind)//" >> "//trim(output)                                                     !!>> HC 21-11-2020
     call system(trim(cat))                                                                             !!>> HC 21-11-2020 save to the general output file

     if (nsubs.ne.nmax)then                             !!>> HC 21-11-2020 IF THERE ARE STILL SUBSTITUTIONS TO DO
        if (allocated(newidum)) deallocate(newidum)     !!>> HC 21-11-2020
        allocate(newidum(largo))                        !!>> HC 21-11-2020
        newidum=0                                       !!>> HC 21-11-2020
        do i=1,largo                                    !!>> HC 21-11-2020 RANDOM SEED FOR THE NEXT MUTATION
           call random_number(a)                        !!>> HC 21-11-2020
           newidum(i)=int(a*1000)                       !!>> HC 21-11-2020
        enddo                                           !!>> HC 21-11-2020
        cp="cp "//trim(rseed)//" "//trim(oldseed)       !!>> AL: kinda pointless really (?)
        call system(cp)                                 !!>> HC 21-11-2020
        open(386,file="rseed.dat")                      !!>> HC 21-11-2020 
            write(386,*) largo                          !!>> HC 21-11-2020
            write(386,*) newidum                        !!>> HC 21-11-2020
        close(386)                                      !!>> HC 21-11-2020
        
        call random_number(a)                                   !!>> HC 21-11-2020 CHOOSE THE NEW INDIVIDUAL TO MUTATE AND RUN 
        thechosen=ceiling(a*npop)                               !!>> HC 21-11-2020 
        if(thechosen==0) thechosen=1                            !!>> HC 21-11-2020 
        write(id,'(I5.5)') thechosen                            !!>> HC 21-11-2020 
        print *,"This individual will start EMaker:", thechosen
        !assnfit=trim(directory)//"/assynfit"//trim(id)//".txt"                                  !!>>AL 15-4-2024
        !cp='sed -i "1s/.*/$(awk NR=='//trim(id)//' '//trim(fitnass)//')/" '//trim(assnfit)      !!>>AL 16-4-2024 For an esoteric reason (to me), this gives problems at MareNostrum5 supercomputer
        !cp="awk '{if(NR=="//trim(id)//") print $1, $2}' "//trim(fitnass)//" > "//trim(assnfit)  !!>>AL 10-9-2024
        i = 0                                                    !!>>AL 10-9-2024
        open(666,file=trim(fitnass))                             !!>>AL 10-9-2024
            do while(i .ne. thechosen)                           !!>>AL 10-9-2024
               read(666,*) fitn, assn                           !!>>AL 10-9-2024
               i = i+1                                          !!>>AL 10-9-2024
            enddo                                                !!>>AL 10-9-2024
        close(666)                                               !!>>AL 10-9-2024       
        open(666,file=trim(assnfit))                             !!>>AL 10-9-2024
            write(666,*) fitn, assn                              !!>>AL 10-9-2024
        close(666)                                               !!>>AL 10-9-2024
        !call system(trim(cp))                                    !!>>AL 10-9-2024
        newid=trim(population)//"/individual"//trim(id)//".dat" !!>> HC 21-11-2020 

        cp="cp "//trim(newid)//" "//trim(indfile)               !!>> HC 21-11-2020 
        call system(trim(cp))                                   !!>> HC 21-11-2020  copy the new individual
        open(1,file=trim(parent))                               !!>> HC 21-11-2020 
            write(1,*) popid(thechosen), thechosen              !!>> HC 21-11-2020 
        close(1)                                                !!>> HC 21-11-2020 
        rm="rm "//trim(finished)                                !!>> HC 21-11-2020  remove the file indicating that this individual has finished
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        rm="rm "//trim(datfitness)                              !!>> HC 21-11-2020  remove the fitness of the substituted file
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        rm="rm "//trim(mutacode)                                !!>> HC 21-11-2020  remove the mutacode of the substituted file
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        rm="rm "//trim(robfile)                                 !!>> HC 21-11-2020  remove the rob file of the substituted file
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        rm="rm "//trim(idfile)                                  !!>> HC 21-11-2020  remove the id file of the substituted file
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        rm="rm "//trim(datfin)                                  !!>> HC 21-11-2020  remove the final file of the substituted file
        call system(trim(rm))                                   !!>> HC 21-11-2020 
        open(1, file=trim(start))                               !!>> HC 21-11-2020 
            write(1,*) "This individual can start running again"!!>> HC 21-11-2020 
        close(1)                                                !!>> HC 21-11-2020 
     else                                                       !!>> HC 21-11-2020  We reached the maximum number of substiturions: THE SIMULATION IS OVER
        allfin=trim(finish)//"/allfin.dat"                      !!>> HC 21-11-2020 
        open(1,file=trim(allfin))                               !!>> HC 21-11-2020 
            write(1,*) "This is the end, all the individuals have run"     !!>> HC 21-11-2020   
        close(1)                                                !!>> HC 21-11-2020 
     endif                                                      !!>> HC 21-11-2020 
  
  else                                                          !!>> HC 21-11-2020  We reached the maximum number of substitutions: THE SIMULATION IS OVER
     allfin=trim(finish)//"/allfin.dat"                         !!>> HC 21-11-2020 
     open(1,file=trim(allfin))                                  !!>> HC 21-11-2020 
         write(1,*) "This is the end, all the individuals have run"  !!>> HC 21-11-2020   
     close(1)                                                   !!>> HC 21-11-2020 
  endif                                                         !!>> HC 21-11-2020 
  
  if (nsubs.ge.nmax)then                                        !!>> HC 21-11-2020 
     allfin=trim(finish)//"/allfin.dat"                         !!>> HC 21-11-2020 
     open(1,file=trim(allfin))                                  !!>> HC 21-11-2020 
         write(1,*) "This is the end, all the individuals have run"  !!>> HC 21-11-2020   
     close(1)                                                   !!>> HC 21-11-2020 
  endif                                                         !!>> HC 21-11-2020 
  
  print*, "substitution", nsubs                                 !!>> HC 21-11-2020 
  
stop
  666 print*, "[substitute_AL.sh]: ERROR! 666 ERROR! individual.dat file. Something went wrong in EMaker", trim(finished)   !!>> HC 21-11-2020 
    cat="cat "//trim(finished)                                  !!>> HC 21-11-2020 
    call system(cat)                                            !!>> HC 21-11-2020
   if (nsubs.ne.nmax)then                             !!>> HC 21-11-2020 IF THERE ARE STILL SUBSTITUTIONS TO DO
         if (allocated(newidum)) deallocate(newidum)     !!>> HC 21-11-2020
         allocate(newidum(largo))                        !!>> HC 21-11-2020
         newidum=0                                       !!>> HC 21-11-2020
         do i=1,largo                                    !!>> HC 21-11-2020 RANDOM SEED FOR THE NEXT MUTATION
            call random_number(a)                        !!>> HC 21-11-2020
            newidum(i)=int(a*1000)                       !!>> HC 21-11-2020
         enddo                                           !!>> HC 21-11-2020
         cp="cp "//trim(rseed)//" "//trim(oldseed)        !!kinda pointless really
         call system(cp)                                 !!>> HC 21-11-2020
         open(386,file="rseed.dat")                      !!>> HC 21-11-2020
            write(386,*) largo                          !!>> HC 21-11-2020
            write(386,*) newidum                        !!>> HC 21-11-2020
         close(386)                                      !!>> HC 21-11-2020

         call random_number(a)                                   !!>> HC 21-11-2020 CHOOSE THE NEW INDIVIDUAL TO MUTATE AND RUN
         thechosen=ceiling(a*npop)                               !!>> HC 21-11-2020
         if(thechosen==0) thechosen=1                            !!>> HC 21-11-2020
         write(id,'(I5.5)') thechosen                            !!>> HC 21-11-2020

         !assnfit=trim(directory)//"/assynfit"//trim(id)//".txt"   !!>>AL 15-4-2024
         !cp='sed -i "1s/.*/$(awk NR=='//trim(id)//' '//trim(fitnass)//')/" '//trim(assnfit) !!>>AL 16-4-2024
         i = 0                                                    !!>>AL 10-9-2024
         open(666,file=trim(fitnass))                             !!>>AL 10-9-2024
             do while(i .ne. thechosen)                           !!>>AL 10-9-2024
                read(666,*) fitn, assn                           !!>>AL 10-9-2024
                i = i+1                                          !!>>AL 10-9-2024
             enddo                                                !!>>AL 10-9-2024
         close(666)                                               !!>>AL 10-9-2024       
         open(666,file=trim(assnfit))                             !!>>AL 10-9-2024
             write(666,*) fitn, assn                              !!>>AL 10-9-2024
         close(666)                                               !!>>AL 10-9-2024
         !call system(trim(cp))
         newid=trim(population)//"/individual"//trim(id)//".dat" !!>> HC 21-11-2020

         cp="cp "//trim(newid)//" "//trim(indfile)               !!>> HC 21-11-2020
         call system(trim(cp))                                   !!>> HC 21-11-2020  copy the new individual
         open(1,file=trim(parent))                               !!>> HC 21-11-2020
            write(1,*) popid(thechosen), thechosen              !!>> HC 21-11-2020
         close(1)                                                !!>> HC 21-11-2020
         rm="rm "//trim(finished)                                !!>> HC 21-11-2020  remove the file indicating that this individual has finished
         call system(trim(rm))                                   !!>> HC 21-11-2020
         rm="rm "//trim(datfitness)                              !!>> HC 21-11-2020  remove the fitness of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020
         rm="rm "//trim(mutacode)                                !!>> HC 21-11-2020  remove the mutacode of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020
         rm="rm "//trim(robfile)                                 !!>> HC 21-11-2020  remove the rob file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020
         rm="rm "//trim(idfile)                                  !!>> HC 21-11-2020  remove the id file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020
         rm="rm "//trim(datfin)                                  !!>> HC 21-11-2020  remove the final file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020
         open(1, file=trim(start))                               !!>> HC 21-11-2020
            write(1,*) "This individual can start running again"!!>> HC 21-11-2020
         close(1)
   end if
end program substitute
