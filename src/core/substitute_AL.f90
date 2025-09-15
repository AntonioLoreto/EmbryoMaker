program substitute
   
   implicit none
   logical         :: exist
   integer         :: i, j, nsubs, nmax, npop, theone, todie, thechosen,thechosen_recombinator,nseed,nminfits,largo,nrecs,recnum
   integer         :: donorsubnum, errorr,inter,waittime,iniwaittime
   real*8          :: indfit, minfit, maxfit, a, fitn, assn,recombination_p,fitrechosen,assrechosen,tiempo
   character*6     :: id
   character*6     :: mostfit
   character*700   :: bindir, population, num, finish, finished, counter, running, fitnesses, popsize,assnfit,fitnass
   character*700   :: directory, indfile, datfitness, allfin, newid, besty, best, datfin, parent, time, popids
   character*700   :: output, outind, robfile, mutacode, idfile, start, rseed, oldseed,outfile,voll,rseedsubs,onlytest
   character*700   :: recombination_donor,donorid,recnumtxt,counter_recs,counter_subs,recinfo,recs_history,donorecinfo
   character*700   :: fitparentnrecomb,mutacodehis,traitshis, traitlcl
   character*1400  :: cp,rm,cat,all_morfs
   character*1400  :: paste 
   integer, allocatable, dimension(:) :: idum, newidum, popid
   real*8, allocatable, dimension(:)  :: fits

   character(len=500) :: filename
   character(len=500) :: latest_file
   integer :: ioss, numm,pos,file_count
   character(len=500) :: command               !!>> AL 9-4-21 

   real*8 :: last_most_fit                     !!>> AL 16-5-25
   integer :: dash_pos, dot_dat, len_f
   character(len=:), allocatable :: rawnum
   character(len=10) :: len_str
   character(len=20) :: read_format
   real*8 :: u = 0.00001d0                       !!>> AL 28-5-25

!important: !AL 22-11-24
   
   recombination_p = 0.1

   iniwaittime = 2000000
  
   call getarg(1,bindir)                                    !!>> HC 21-11-2020
   call getarg(2,finished)                                  !!>> HC 21-11-2020 This is the file that has been produced in the directory finished indicating
   call getarg(3,num)                                      

   open(1,file=trim(finished))                              !!>> HC 21-11-2020 that one CPU has finished running EMaker
      read(1,*,END=666) theone                              !!>> HC 21-11-2020 the one that has finished
   close(1)                                                 !!>> HC 21-11-2020

   write(id,'(I6.6)') theone                                !!>> HC 21-11-2020

   read (num,*) nmax 

   running     =trim(bindir)//"/running"                       !!>> HC 21-11-2020 path to the directory where we have the CPUs running
   output      =trim(bindir)//"/output_model.dat"              !>> HC 21-11-2020 path to the file where we save the general output of the evolutionary model
   rseedsubs   =trim(bindir)//"/rseedsubs.dat"                 !!>> AL 5-4-24
   mutacodehis =trim(bindir)//"/mutacodehis.dat"               !!>> AL 21-2-25
   traitshis   =trim(bindir)//"/traitshis.dat"                   !!>> AL 12-3-25
   
   population  =trim(bindir)//"/population"                   !!>> HC 21-11-2020 path to the directory where we have the population files
   fitnesses   =trim(population)//"/population.datfitness"    !!>> HC 21-11-2020 path to the file where we save the population fitnesses 
   outfile     =trim(population)//"transout.dat"              !!>> HC 21-11-2020 transient file
   popsize     =trim(population)//"/npop.dat"                 !!>> HC 21-11-2020 path to the file where we save the number of individuals in the population
   popids      =trim(population)//"/popids.dat"               !!>> HC 21-11-2020 path to the file where we save the ids of the individuals in the population
   fitnass     =trim(population)//"/fitnassy.txt"             !!>> AL 15-4-2024
   onlytest    =trim(population)//"/onlytest.txt"             !!>> AL 13-1-2025

   finish      =trim(bindir)//"/finished"                     !!>> HC 21-11-2020 path to the directory where we save the CPUs that have finished
   counter_subs=trim(finish)//"/counter_subs.dat"             !>> HC 21-11-2020 path to the file where we save the number of substitutions we have run  !!>> AL 26-11-24:changed name
   counter_recs=trim(finish)//"/counter_recombinations.dat"   !!>> AL 26-11-24 count how many recombinations. This is info is needed to reconstruct lineage
   recs_history=trim(finish)//"/recombination_history.txt"

   directory   =trim(running)//"/"//trim(id)                  !!>> HC 21-11-2020 path to the individual that has finished
   indfile     =trim(directory)//"/individual.dat"            !!>> HC 21-11-2020 path to the individual's ic file
   datfitness  =trim(directory)//"/individual.datfitness"     !!>> HC 21-11-2020 path to the individual's fitness
   datfin      =trim(directory)//"/individual.dat*dat"        !!>> HC 21-11-2020 path to the output file (morphology)
   parent      =trim(directory)//"/parent.dat"                !!>> HC 21-11-2020 path to the individual's parent index
   robfile     =trim(directory)//"/rob.val"                   !!>> HC 21-11-2020 path to the individual's robustness
   mutacode    =trim(directory)//"/mutacode.dat"              !!>> HC 21-11-2020 path to the individual's code of mutation
   outind      =trim(directory)//"/outind.dat"                !!>> HC 21-11-2020 transient file !!>> AL 9-1-2025: This is generated in this script and is the info for output_model.dat
   idfile      =trim(directory)//"/idfile.dat"                !!>> HC 21-11-2020 path to the new ID of this individual
   start       =trim(directory)//"/start.dat"                 !!>> HC 21-11-2020 path to the file that will allow EMaker to run again
   time        =trim(directory)//"/development_time.dat"      !!>> HC 21-11-2020 path to the individual's runtime
   rseed       =trim(directory)//"/rseed.dat"                 !!>> HC 21-11-2020 path to the random seed
   oldseed     =trim(directory)//"/rseed_original.dat"        !!>> HC 21-11-2020 path to the original random seed !!AL>> this actually is pointless
   voll        =trim(directory)//"/individual.volume.txt"     !!>>AL 15-4-2024
   assnfit     =trim(directory)//"/assynfit.txt"              !!>>AL 15-4-2024
   recinfo     =trim(directory)//"/recombination_info.txt"    !!>> AL 26-11-24
   recnumtxt   =trim(directory)//"/ifrec_numrec.txt"          !!>> AL 27-11-24
   donorecinfo =trim(directory)//"/donorecinfo.txt"           !!>> AL 27-11-24
   traitlcl    =trim(directory)//"/traits_local.dat"          !!>> AL 12-3-25
   fitparentnrecomb     =trim(directory)//"/fitparentnrecomb.txt"          !!>> AL 27-11-24
   recombination_donor  =trim(directory)//"/recombination_donor.dat"       !!>> AL 20-11-24

   !!>> AL 29-1-25 important: only for testing reconstruction 
   !all_morfs=trim(finish)//"/all_morfs/"           !!>> AL 29-1-25 

   !***************************************************************************! !>>AL 13-11-24: this ideally should be a subroutine to improve readability 
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

   !> AL 20-1-25: Error handling. If some of the next files doesnt exist for whatever mysterious reason, we totally ignore this subs.  
   errorr = 0
   call check_file(counter_subs,errorr)
   call check_file(parent,errorr)
   call check_file(datfitness,errorr)
   call check_file(robfile,errorr)
   call check_file(voll,errorr)
   call check_file(mutacode,errorr)
   call check_file(time,errorr)

   if(errorr == 1)then 
      print *,"Some file was not found. We ignore this subsitution as it never happened.", &
      "We will finish the simulation, look were the fockin problem is"
      nsubs = nmax
      go to 6666
   end if

   open(6021, file=trim(datfitness))
      read(6021, *) indfit
   close(6021)
   !print *,"[substitute_AL.sh]: this is fitness inside individual.datfitness ", indfit

   open(1,file=trim(popsize))                               !!>> HC 21-11-2020
         read(1,*) npop                                     !!>> HC 21-11-2020 number of individuals
   close(1)                                                 !!>> HC 21-11-2020

   open(1,file=trim(counter_subs))                          !!>> HC 21-11-2020
      read(1,*) nsubs                                       !!>> HC 21-11-2020 number of substitutions we have already done
   close(1)                                                 !!>> HC 21-11-2020

   if(recombination_p == 0)then 
      nrecs = 0
   else 
      open(1,file=trim(counter_recs))                          !!>> AL 26-11-2024 number of recombinations so far   
         read(1,*) nrecs                                       
      close(1)        
   end if                                         

   if(allocated(fits)) deallocate(fits)                     !!>> HC 21-11-2020
   allocate(fits(npop)); fits=666                           !!>> HC 21-11-2020
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

   !-- Check if 'best' folder is empty                      !!>> AL 2-6-25
   filename = trim(bindir)//"/best/"
   command = "ls -1 "//trim(filename)//" | wc -l > file_count.txt"
   call system(trim(command))

   !-- Read the count
   open(unit=11, file="file_count.txt", action="read")
      read(11, *) file_count
   close(11)
   call system("rm file_count.txt")

   if(file_count /= 0)then
      !— find the latest “individual…” file                    !!>> AL 2-6-25
      filename = trim(bindir)//"/best/"
      command = "ls -t "//trim(filename)//" | head -n 1 > latest_file.txt"
      call system(trim(command))

      !— read its name
      open(unit=10, file="latest_file.txt", action="read", iostat=ioss)
         if (ioss /= 0) then
            print *, "Error opening latest_file.txt"
            go to 19
         end if
         read(10,'(A)') latest_file
      close(10)

      !— parse out the integer “numm”
      pos = index(latest_file, "individual") + 10
      read(latest_file(pos:pos+5), '(I6)') numm

      !— now parse the fitness after the “-” and before “.dat”
      dash_pos = index(latest_file, "-")
      dot_dat  = index(latest_file, ".dat")
      
      if (dash_pos == 0 .or. dot_dat == 0) then   
         print *, "Filename format unexpected: ", latest_file
         go to 19
      else

         len_f = dot_dat - dash_pos - 1
         allocate(character(len=len_f) :: rawnum)
         rawnum = latest_file(dash_pos+1 : dash_pos+len_f)
         do i = 1, len(rawnum)
            if (rawnum(i:i) == '_') rawnum(i:i) = '.'
         end do
         ! Read the number from rawnum
         read(rawnum, '(F6.4)') last_most_fit
      end if
      call system("rm latest_file.txt")

   end if 

   !filtro de saturacion 16-5-25
   if(nsubs>iniwaittime)then                                       !!>> AL 24-3-25 we end simulation if no increase in fitness during arbitrary number of subs (we reached saturation)
      !— now you can use last_most_fit
      waittime = numm*2
      if((nsubs-numm)>waittime.and.(last_most_fit-last_most_fit*u) <= indfit)nsubs=nmax !the last condition is redundant since no morfs with less distance than this treshold are saved
19     continue
   end if

   if (nsubs .le. nmax) then                                !!>> HC 21-11-2020 IF there are still substitutions to do  

      open(1,file=trim(counter_subs))                       !!>> HC 21-11-2020 write the new number of done substitutions
         write(1,*) nsubs                                   !!>> HC 21-11-2020
      close(1)                                              !!>> HC 21-11-2020

      print*,"[substitute_AL]: starting substitution #",nsubs
      
      !!>> AL 26-11-24 check if there was a recombination in this subss 
      recnum = 0
      
      inquire(file=trim(recinfo), exist=exist)                     
      if(exist)then
         nrecs=nrecs+1
         open(1,file=trim(counter_recs))                    !!>> AL 26-11-2024 update number of recombinations  
            write(1,*) nrecs                                       
         close(1)       

         recnum = nrecs                                     !!>> AL 27-11-2024 

         open(1,file=trim(recs_history),status="old",position="append")                                
            write(1,*)"#",nrecs                               
         close(1)

         cat = "cat "//trim(recinfo)//" >> "//trim(recs_history)
         call system(trim(cat))                              !!>> AL 26-11-24 we save recombination info 

         cp = "rm "//trim(recinfo)                             !!>> AL 26-11-24 remove evidence of recombination
         call system(trim(cp))
         cp = "rm "//trim(recombination_donor)              
         call system(trim(cp))   
      else
         cp = "echo '-1' > "//trim(donorecinfo)                 !!>> AL 9-12-2024 this is for paste
         call system(trim(cp))         
      end if

      open(1,file=trim(recnumtxt))                          !!>> AL 27-11-2024 this is for paste
         write(1,*) recnum                                       
      close(1)  

      minfit=minval(fits)                                   !!>> HC 21-11-2020 minimum fitness value in the population
      maxfit=maxval(fits)                                   !!>> HC 21-11-2020

      open(1,file=trim(time))                               !!>> AL 4-2-2025
         read(1,*) tiempo                                   
      close(1)                                              

      if(tiempo /= 0 .and. tiempo < 1)then                  !!>>AL 20-11-25  important: if dev. way too fast there was an error
         go to 90909
      end if 

     ! !>>AL 4-2-2025: important, just for testing:we save info of the mutated morf
     ! write(id,'(I6.6)') nsubs                                  
     ! newid=trim(all_morfs)//"indatsubs-"//trim(id)//".dat"     
     ! cp="cp "//trim(indfile)//" "//trim(newid)                 
     ! call system(trim(cp))                                     
      
      !IF SUBSTITUTION IS ACCEPTED                           !!>>AL 20-11-24
      if (indfit .le. maxfit)then                            !!>> HC 21-11-2020 This individual has a fitness higher or equal than the minimum
      
         nminfits=0                                         !!>> HC 21-11-2020 number of individuals with the minimum fitness
         do i=1,npop                                        !!>> HC 21-11-2020
            if(fits(i) < maxfit)cycle                       !!>> HC 21-11-2020
            nminfits = nminfits + 1                         !!>> HC 21-11-2020
         enddo                                              !!>> HC 21-11-2020
         
         call random_number(a)                               !!>> HC 21-11-2020 This needs to be instantiated
         nminfits = nint(nminfits*a)                           !!>> HC 21-11-2020 The chosen individual with the smallest fitness
         if(nminfits==0) nminfits=1                          !!>> HC 21-11-2020
      
         j = 0                                                 !!>> HC 21-11-2020
         do i = 1,npop
            if(fits(i) < maxfit)cycle                         !!>> HC 4-3-2024 Go through the individuals that have the smallest fitness value
            j = j + 1                                           !!>> HC 4-3-2024 the number of those
            if(nminfits == j)then                             !!>> HC 4-3-2024 this is the randomly chosen one
               todie    = i                                       !!>> HC 4-3-2024 id of the individual with the minimum fitness value 
               fits(i)  = indfit                                !!>> HC 4-3-2024 substitute in the fitness list
               popid(i) = nsubs                                !!>> HC 4-3-2024 The id of this new individual will be the number of this substitution
               exit                                          !!>> HC 4-3-2024
            endif                                           !!>> HC 4-3-2024
         enddo                                              !!>> HC 4-3-2024

         write(id,'(I6.6)') todie                                  !!>> HC 21-11-2020 THE INDIVIDUAL TO BE SUSTITUTED BY
         newid=trim(population)//"/individual"//trim(id)//".dat"   !!>> HC 21-11-2020
         cp="cp "//trim(indfile)//" "//trim(newid)                 !!>> HC 21-11-2020
         call system(trim(cp))                                     !!>> HC 21-11-2020 !!>>AL 20-11-2024: we make substitution

         !!>>AL 4-2-25: the follwing is needlessly esoteric (but works). Change!
         cp='sed -i "'//trim(id)//'s/.*/$(cat '//trim(assnfit)//')/" '//trim(fitnass) !!>>AL 15-4-2024: replace assymetry and fitness in population file (fitnass) for the new population integrant 
         call system(trim(cp))                                                         
         
         open(1,file=trim(idfile))                           !!>> HC 21-11-2020
            write(1,*) todie                                 !!>> HC 21-11-2020
         close(1)                                            !!>> HC 21-11-2020
         
         open(1,file=trim(fitnesses))                        !!>> HC 21-11-2020 !!>>AL 4-2-25:you actually dont need to rewrite the whole array. Change!
            do i=1,npop                                      !!>> HC 21-11-2020
               write(1,*) fits(i)                            !!>> HC 21-11-2020 new fitnesses list
            enddo                                            !!>> HC 21-11-2020
         close(1)                                            !!>> HC 21-11-2020
         
         open(1,file=trim(popids))                           !!>> HC 21-11-2020 !!>>AL 4-2-25:you actually dont need to rewrite the whole array. Change!
            do i=1,npop                                      !!>> HC 21-11-2020
               write(1,*) popid(i)                           !!>> HC 21-11-2020 new individual ids
            enddo                                            !!>> HC 21-11-2020
         close(1)                                            !!>> HC 21-11-2020
         
         if(file_count /= 0)then 
            if (indfit < (last_most_fit - last_most_fit*u))then                                            !!>> HC 21-11-2020 This individual has a maximum fitness value
            
               best = trim(bindir)//"/best"                                     !!>> HC 21-11-2020
               write(id,'(I6.6)') nsubs                                         !!>> HC 21-11-2020
               write(mostfit,'(F6.4)') indfit                                   !!>> HC 21-11-2020
               
               do i = 1, len(mostfit)
                  if (mostfit(i:i) == '.') mostfit(i:i) = '_'
               end do
               
               besty = trim(bindir)//"/best/individual"//trim(id)//"-"//trim(mostfit)//".dat"       !!>> HC 21-11-2020
            
               cp = "cp "//trim(datfin)//" "//trim(besty)                       !!>> HC 21-11-2020
               call system(trim(cp))                                          !!>> HC 21-11-2020 save it in the "best" directory
            endif                                                             !!>> HC 21-11-2020
         else 
            if (indfit < minfit)then                                            !!>> HC 21-11-2020 This individual has a maximum fitness value
            
               best = trim(bindir)//"/best"                                     !!>> HC 21-11-2020
               write(id,'(I6.6)') nsubs                                         !!>> HC 21-11-2020
               write(mostfit,'(F6.4)') indfit                                   !!>> HC 21-11-2020
               
               do i = 1, len(mostfit)
                  if (mostfit(i:i) == '.') mostfit(i:i) = '_'
               end do
               
               besty = trim(bindir)//"/best/individual"//trim(id)//"-"//trim(mostfit)//".dat"       !!>> HC 21-11-2020
            
               cp = "cp "//trim(datfin)//" "//trim(besty)                       !!>> HC 21-11-2020
               call system(trim(cp))                                          !!>> HC 21-11-2020 save it in the "best" directory
            endif                                                             !!>> HC 21-11-2020
         end if 
      
      else                                                                 !!>> HC 21-11-2020 THIS INDIVIDUAL IS NOT BETTER THAN THE WORST IN THE POPULATION

90909    continue
         open(1,file=trim(idfile))                                         !!>> HC 21-11-2020
            write(1,*) "          NA"                                      !!>> HC 21-11-2020 NA means this individual has been discarded !!AL 13-11-24: because fitness was not better thatn the worst one
         close(1)                                                          !!>> HC 21-11-2020

      endif                                                                !!>> HC 21-11-2020
      
      !paste = "paste "//trim(counter_subs)//" "//trim(idfile)//" "//trim(parent)//" "//trim(datfitness)//" "// &
      !      trim(robfile)//" "//trim(voll)//" "//trim(mutacode)//" "//trim(time)//" "//trim(recnumtxt)//" "// &
      !      trim(donorecinfo)//" > "//trim(outind)     

      paste = "paste "//trim(counter_subs)//" "//trim(idfile)//" "//trim(parent)//" "//trim(datfitness)//" "// &
            trim(robfile)//" "//trim(voll)//" "//trim(time)//" "//trim(recnumtxt)//" "// &
            trim(donorecinfo)//" > "//trim(outind)   
      call system(trim(paste))                                                                           !!>> HC 21-11-2020 output variables of this individual
      
      cat = "cat "//trim(outind)//" >> "//trim(output)                                                   !!>> HC 21-11-2020
      call system(trim(cat))                                                                             !!>> HC 21-11-2020 save to the general output file
      
      open(1,file=trim(mutacodehis),position="append")                                                   !!>> AL 21-2-25       
         write(1,*)"#",nsubs                               
      close(1)

      cat = "cat "//trim(mutacode)//" >> "//trim(mutacodehis)                                            !!>> AL 21-2-25
      call system(trim(cat))    

      cp = "rm "//trim(donorecinfo)                                                                      !!>> AL 10-1-25 remove evidence of recombination
      call system(trim(cp))

      cp = "rm "//trim(outind)                                                                           !!>> AL 4-2-25 clean after yourself
      call system(trim(cp))

      cp = "rm "//trim(mutacode)                                                                         !!>> AL 4-2-25 clean after yourself
      call system(trim(cp))

      inquire(file=trim(traitlcl), exist=exist)                                                          !!>> AL 12-3-25            
      if(exist)then
         cat = "paste "//trim(counter_subs)//" "//trim(traitlcl)// " >> "//trim(traitshis)               !!>> AL 12-3-25
         call system(trim(cat))    
         
         cp = "rm "//trim(traitlcl)                                                                      !!>> AL 12-3-25 clean after yourself
         call system(trim(cp)) 
      end if

      !!>> HC 21-11-2020 IF THERE ARE STILL SUBSTITUTIONS TO DO
6666  continue 
      if (nsubs /= nmax)then                             

         if (allocated(newidum)) deallocate(newidum)     !!>> HC 21-11-2020
         allocate(newidum(largo))                        !!>> HC 21-11-2020
      
         newidum=0                                       !!>> HC 21-11-2020
      
         do i = 1, largo                                 !!>> HC 21-11-2020 RANDOM SEED FOR THE NEXT MUTATION
            call random_number(a)                        !!>> HC 21-11-2020
            newidum(i) = int(a*1000)                      !!>> HC 21-11-2020
         enddo                                           !!>> HC 21-11-2020

         cp="cp "//trim(rseed)//" "//trim(oldseed)       
         call system(cp)                                 !!>> HC 21-11-2020
         open(386,file="rseed.dat")                      !!>> HC 21-11-2020 
            write(386,*) largo                           !!>> HC 21-11-2020
            write(386,*) newidum                         !!>> HC 21-11-2020
         close(386)                                      !!>> HC 21-11-2020
 
         inquire(file=trim(recinfo), exist=exist)                 !!>> AL 20-1-25  !!>>AL 4-2-25: This should be necessary          
         if(exist)then
            cp = "rm "//trim(recinfo)                             !!>> AL AL 20-1-25  remove evidence of prior recombination
            call system(trim(cp))
         end if

         inquire(file=trim(recombination_donor), exist=exist)     !!>> AL 20-1-25      !!>>AL 4-2-25: This should be necessary        
         if(exist)then
            cp = "rm "//trim(recombination_donor)                  !!>> AL 20-1-25  remove evidence of prior recombination
            call system(trim(cp))
         end if

         inquire(file=trim(donorecinfo), exist=exist)              !!>> AL 20-1-25      !!>>AL 4-2-25: This should be necessary        
         if(exist)then
            cp = "rm "//trim(donorecinfo)                          !!>> AL 20-1-25  remove evidence of prior recombination
            call system(trim(cp))
         end if

         !choose if there will be recombination
         thechosen_recombinator = 0 
         call random_number(a)                                   
         if(a .lt. recombination_p)then

            call random_number(a)                                   !!>> AL 3-2-2025 choose the recombination receiver  
            thechosen = ceiling(a*npop)                                
            if(thechosen == 0) thechosen = 1                             

            open(666, file = trim(fitnass))                             
               do i = 1, npop
                  read(666,*) assn, fitn                                                                        
                  if(i == thechosen) exit
               end do 
            close(666) 
               
            call random_number(a)                                    !!>> AL 3-2-2025 choose the recombination donor
            thechosen_recombinator = ceiling(a*npop)
            if(thechosen_recombinator==0) thechosen_recombinator = 1  

            do while(thechosen_recombinator == thechosen)            !!>> AL 21-1-25 donor cannot be the same as the receiver
               call random_number(a)
               thechosen_recombinator = ceiling(a*npop)
               if(thechosen_recombinator == 0) thechosen_recombinator = 1 
            end do 

            open(666,file=trim(fitnass))                             
               do i = 1, npop
                  read(666,*) assrechosen, fitrechosen                                                                        
                  if(i == thechosen_recombinator) exit
               end do 
            close(666) 

            !!>>AL 7-3-25: Donor must always be of more or equal distance to optimum than reciever. 
            if(fitn > fitrechosen)then 
               inter = thechosen
               thechosen = thechosen_recombinator
               thechosen_recombinator = inter
            end if 
            
         end if 

         if(thechosen_recombinator == 0) then                       !!>>AL 3-2-2025 therefore no recombination   
            call random_number(a)                                   !!>> HC 21-11-2020 CHOOSE THE NEW INDIVIDUAL TO MUTATE AND RUN 
            thechosen=ceiling(a*npop)                               !!>> HC 21-11-2020 
            if(thechosen==0) thechosen=1                            !!>> HC 21-11-2020 
            write(id,'(I6.6)') thechosen                            !!>> HC 21-11-2020 
            
            print *,"This individual will start EMaker:", thechosen
            newid=trim(population)//"/individual"//trim(id)//".dat" !!>> HC 21-11-2020 
            cp="cp "//trim(newid)//" "//trim(indfile)               !!>> HC 21-11-2020 
            call system(trim(cp))                                   !!>> HC 21-11-2020  copy the new individual

            open(666,file=trim(fitnass))                             
               do i = 1, npop
                  read(666,*) assn, fitn                                                                        
                  if(i == thechosen) exit
               end do 
            close(666) 

            open(666,file=trim(assnfit))                             !!>>AL 10-9-2024
               write(666,*) assn, fitn                               !
            close(666)                                               !!>>AL 10-9-2024

            open(1,file=trim(parent))                               !!>> HC 21-11-2020 
               write(1,*) popid(thechosen), thechosen                
            close(1)                                                !!>> HC 21-11-2020 
         else                                                        !!>>AL 3-2-2025 intimacy   
            !print*,'[substitute_AL]: we will do a recombination at subs: ',nsubs

            write(id,'(I6.6)') thechosen_recombinator         
            donorid=trim(population)//"/individual"//trim(id)//".dat"  
            cp="cp "//trim(donorid)//" "//trim(recombination_donor)             
            call system(trim(cp))                                   !! AL 21-11-24 We copy donor ico into $directory 

            open(666,file=trim(fitnass))                             
               do i = 1, npop
                  read(666,*) assrechosen, fitrechosen                                                                        
                  if(i == thechosen_recombinator) exit
               end do 
            close(666) 

            open(1,file=trim(donorecinfo))                                        
               write(1,*) popid(thechosen_recombinator)        !AL 9-12-24 important, recdonor sub number (HC called it ind. id which i hate)
            close(1) 

            write(id,'(I6.6)') thechosen   
            newid=trim(population)//"/individual"//trim(id)//".dat" !!>> HC 21-11-2020 
            cp="cp "//trim(newid)//" "//trim(indfile)               !!>> HC 21-11-2020 
            call system(trim(cp))                                   !!>> HC 21-11-2020  copy the new individual

            open(666,file=trim(fitnass))                             
               do i = 1, npop
                  read(666,*) assn, fitn                                                                        
                  if(i == thechosen) exit
               end do 
            close(666) 

            open(666,file=trim(assnfit))                             !!>>AL 10-9-2024
               write(666,*) assn, fitn                               !
            close(666)                                               !!>>AL 10-9-2024

            open(1,file=trim(parent))                               !!>> HC 21-11-2020 
               write(1,*) popid(thechosen), thechosen                
            close(1)                                                !!>> HC 21-11-2020 
         end if

         rm="rm "//trim(finished)                                !!>> HC 21-11-2020  remove the file indicating that this individual has finished
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(datfitness)                              !!>> HC 21-11-2020  remove the fitness of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         !rm="rm "//trim(mutacode)                                !!>> HC 21-11-2020  remove the mutacode of the substituted file
         !call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(robfile)                                 !!>> HC 21-11-2020  remove the rob file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(idfile)                                  !!>> HC 21-11-2020  remove the id file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(datfin)                                  !!>> HC 21-11-2020  remove the final file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(voll)                                    !!>> HC 21-11-2020  remove the final file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 
         rm="rm "//trim(time)                                    !!>> HC 21-11-2020  remove the final file of the substituted file
         call system(trim(rm))                                   !!>> HC 21-11-2020 

         open(1, file=trim(start))                               !!>> HC 21-11-2020 
            write(1,*) "This individual can start running again" !!>> HC 21-11-2020 
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

      !!>> AL 16-4-25: we erease all /running & /population folders because of file limited space in Picassso HPC
      !command = "rm -r "//trim(running)
      !call system(trim(command))
      !command = "rm -r "//trim(population)
      !call system(trim(command))
   endif                                                         !!>> HC 21-11-2020 

   if (nsubs.ge.nmax)then                                        !!>> HC 21-11-2020 
   
      allfin=trim(finish)//"/allfin.dat"                         !!>> HC 21-11-2020 
      open(1,file=trim(allfin))                                  !!>> HC 21-11-2020 
         write(1,*) "This is the end, all the individuals have run"  !!>> HC 21-11-2020   
      close(1)                                                   !!>> HC 21-11-2020 
      !!>> AL 16-4-25: we erease all /running & /population folders because of file limited space in Picassso HPC
      command = "rm -r "//trim(running)
      call system(trim(command))
      command = "rm -r "//trim(population)
      call system(trim(command))
   endif                                                         !!>> HC 21-11-2020 

   stop
   
   666 print*, "[substitute_AL.sh]: ERROR! 666 ERROR! individual.dat file. Something went wrong in EMaker", trim(finished)   !!>> HC 21-11-2020 
      cat="cat "//trim(finished)                                  !!>> HC 21-11-2020 
      call system(cat)                                            !!>> HC 21-11-2020

   contains

   subroutine check_file(filename,errorr)
      character(len=*), intent(in) :: filename
      logical :: file_exists
      integer :: errorr
      inquire(file=filename, exist=file_exists)
      if (.not. file_exists) then
         print*, "Error: File does not exist: ", trim(filename)
         errorr = 1
      end if
   end subroutine check_file

end program substitute
