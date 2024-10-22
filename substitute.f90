program substitute
  implicit none
  integer :: i, j, nsubs, nmax, npop, theone, todie, thechosen, nseed, toreproduce, nminfits
  character*700 :: bindir, childir, population, num, finish, finished, counter, running, fitnesses, popsize
  character*700 :: directory, indfile, datfitness, allfin, newid, besty, best, datfin, parent, time, popids, cfile
  character*700 :: output, outind, robfile, mutacode, idfile, start, rseed, oldseed,outfile,svgen
  character*1400 :: cp, rm, cat, mkdir
  character*10000 :: paste 
  character*5 :: id
  real*8, allocatable, dimension(:) :: fits
  real*8 :: indfit, minfit, maxfit, a
  integer, allocatable, dimension(:) :: idum, newidum, popid
  
  call getarg(1,bindir)                                    !!>> HC 21-11-2020
  finish=trim(bindir)//"/finished"                         !!>> HC 21-11-2020 path to the directory where we save the CPUs that have finished
  counter=trim(finish)//"/counter.dat"                     !!>> HC 21-11-2020 path to the file where we save the number of substitutions we have run
  population=trim(bindir)//"/population"                   !!>> HC 21-11-2020 path to the directory where we have the population files
  fitnesses=trim(population)//"/population.datfitness"     !!>> HC 21-11-2020 path to the file where we save the population fitnesses 
  running=trim(bindir)//"/running"                         !!>> HC 21-11-2020 path to the directory where we have the CPUs running
  output=trim(bindir)//"/output_model.dat"                 !!>> HC 21-11-2020 path to the file where we save the general output of the evolutionary model
  outfile=trim(population)//"transout.dat"                 !!>> HC 21-11-2020 transient file
  popsize=trim(population)//"/npop.dat"                    !!>> HC 21-11-2020 path to the file where we save the number of individuals in the population
  popids=trim(population)//"/popids.dat"                   !!>> HC 21-11-2020 path to the file where we save the ids of the individuals in the population
  
  call getarg(2,finished)                                  !!>> HC 21-11-2020 This is the file that has been produced in the directory finished indicating
  open(1,file=trim(finished))                              !!>> HC 21-11-2020 that one CPU has finished running EMaker
      read(1,*,END=666) theone                             !!>> HC 21-11-2020 the one that has finished
  close(1)                                                 !!>> HC 21-11-2020
  write(id,'(I5.5)') theone                                !!>> HC 21-11-2020
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
  time=trim(directory)//"/time.dat"                        !!>> HC 21-11-2020 path to the individual's runtime
  rseed=trim(directory)//"/rseed.dat"                      !!>> HC 21-11-2020 path to the random seed
  oldseed=trim(directory)//"/rseed_original.dat"           !!>> HC 21-11-2020 path to the original random seed
!  cfile=trim(directory)//"/complexity_vals.dat"
  
  open(1,file=trim(datfitness))                            !!>> HC 21-11-2020
      read(1,*) indfit                                     !!>> HC 21-11-2020 individual's fitness value
  close(1)                                                 !!>> HC 21-11-2020
  
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
        call random_number(a)                              !!>> HC 21-11-2020
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
     " "//trim(robfile)//" "//trim(mutacode)//" "//trim(time)//" > "//trim(outind)                      !!>> HC 21-11-2020
     call system(trim(paste))                                                                           !!>> HC 21-11-2020 output variables of this individual
     cat="cat "//trim(outind)//" >> "//trim(output)                                                     !!>> HC 21-11-2020
     call system(trim(cat))                                                                             !!>> HC 21-11-2020 save to the general output file

     if (nsubs.ne.nmax)then                             !!>> HC 21-11-2020 IF THERE ARE STILL SUBSTITUTIONS TO DO
     
        open(386,file=trim(rseed))                      !!>> HC 21-11-2020  it was saved here by seeding.f90
            read(386,*)nseed                            !!>> HC 21-11-2020  and it is calculated from the master random seed
            if (allocated(idum)) deallocate(idum)       !!>> HC 21-11-2020  (the one in the starting i/o Emaker file)
            allocate(idum(nseed))                       !!>> HC 21-11-2020
            idum=0                                      !!>> HC 21-11-2020
            read(386,*) idum                            !!>> HC 21-11-2020
        close(386)                                      !!>> HC 21-11-2020

        if (allocated(newidum)) deallocate(newidum)     !!>> HC 21-11-2020
        allocate(newidum(nseed))                        !!>> HC 21-11-2020
        newidum=0                                       !!>> HC 21-11-2020
        call random_seed(put=idum)                      !!>> HC 21-11-2020 
        do i=1,nseed                                    !!>> HC 21-11-2020 RANDOM SEED FOR THE NEXT MUTATION
           call random_number(a)                        !!>> HC 21-11-2020
           newidum(i)=int(a*1000)                       !!>> HC 21-11-2020
        enddo                                           !!>> HC 21-11-2020
        cp="cp "//trim(rseed)//" "//trim(oldseed)
        call system(cp)                                 !!>> HC 21-11-2020
        open(386,file="rseed.dat")                      !!>> HC 21-11-2020 
            write(386,*) nseed                          !!>> HC 21-11-2020
            write(386,*) newidum                        !!>> HC 21-11-2020
        close(386)                                      !!>> HC 21-11-2020
                            
        call random_number(a)                                   !!>> HC 21-11-2020 CHOOSE THE NEW INDIVIDUAL TO MUTATE AND RUN !!>>AL modificar para que la prob. 
        !! de mutar sea proporcional al fitness
        thechosen=ceiling(a*npop)                               !!>> HC 21-11-2020 
        if(thechosen==0) thechosen=1                            !!>> HC 21-11-2020 
        write(id,'(I5.5)') thechosen                            !!>> HC 21-11-2020 
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
  
666 print*, "THERE HAS BEEN AN ERROR READING", trim(finished)   !!>> HC 21-11-2020 
    cat="cat "//trim(finished)                                  !!>> HC 21-11-2020 
    call system(cat)                                            !!>> HC 21-11-2020 
  
  
 

end program substitute
