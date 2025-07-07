!gfortran -O2 -Wall -fbounds-check  general.mod.f90  io.mod.f90 genetic.mod.f90 recombination.f90 -o rec.e
module recombination

use io
use general
use genetic

contains
subroutine chose_genes_to_recombinate(individual2)

    implicit none
   
    integer             :: ihc,count,jhc 
    real*8              :: a
    character*700       :: individual2
    integer, allocatable, dimension(:)  :: lequel 
    real*8, allocatable, dimension(:,:) :: genes4recomb,es,params   
    CHARACTER(len=100) :: format_string1,format_string2
    
    call iniread
    call readsnap(trim(individual2))

    if(allocated(genes4recomb)) deallocate(genes4recomb)                     
    allocate(genes4recomb(ng,ng))
    genes4recomb=0   
   
    if(allocated(es)) deallocate(es)                     
    allocate(es(ng,nga))
    es=0   
    
    if(allocated(params)) deallocate(params)                     
    allocate(params(ng,4))
    params=0   
    
    if(allocated(lequel)) deallocate(lequel)                     
    allocate(lequel(ng))
    lequel=0

    count = 0

    do ihc = 1,ng
        call random_number(a)       
        if(a < 0.5)then         
            count               = count + 1
            lequel(ihc)         = 1
            genes4recomb(ihc,:) = gen(ihc)%t(:)
            es(ihc,:)           = gen(ihc)%e(:)
            params(ihc,1)       = gen(ihc)%diffu
            params(ihc,2)       = gen(ihc)%mu
            params(ihc,3)       = gen(ihc)%mich
            params(ihc,4)       = gen(ihc)%kindof
        end if
    end do 

    ! Construct the format string dynamically
    write(format_string1, '(A,I0,A)') '(I2, ', ng, 'E15.7)'
    write(format_string2, '(A,I0,A)') '(I2, ', nga, 'E15.7)'

    open(1,file="recombination_info.txt")      
        write(1,*) count                        
        do ihc = 1,ng    
            if(lequel(ihc)==0)cycle
            write(1,format_string1) ihc,genes4recomb(ihc, :)
            write(1,format_string2) ihc,es(ihc,:)
            write(1,'(I2,4E15.7)') ihc,params(ihc,1),params(ihc,2),params(ihc,3),params(ihc,4)                
        end do 
    close(1)                                                 

end subroutine chose_genes_to_recombinate

subroutine do_recombination_4_evo(individual1,individual2)

    implicit none
   
    integer         :: count,ihc,jhc,wrsize
    real*8          :: rval
    character*700   :: individual1, individual2, rename
    integer, allocatable, dimension(:)  :: lequel 
    real*8, allocatable, dimension(:,:) :: genes4recomb,es,params
    integer, allocatable, dimension(:)  :: wridum

    print*,"we are in do_recombination_4_evo"
    !RANDOM NUMBER THINGs                                   !!>> AL 27-11-24
    !we will use the file rseed to instantiate random number generation. 
    !This file is the same used in muta_AL.f90, which is generated in seeding and modify in substitute
    call random_seed(size = nseed)                          !!>> AL 27-11-24 this is defined as public variable in general.modf90 
    !call getarg(2,cu)                                       

    if(allocated(idum))deallocate(idum)                     !!>> AL 27-11-24 idum is public variable defined in general and writen to EMaker file in io.f90               
    allocate(idum(nseed))    

    idum=0; wrsize=0                                        
    open(837,file="rseed.dat")                              !!>> HC 20-12-2021 the seed should be stored in the individual directory
        read(837,*) wrsize                                   
        if (allocated(wridum)) deallocate(wridum)            
        allocate(wridum(wrsize))                            
        read(837,*) wridum                                   
    close(837)              

    do ihc=1,size(idum)                                     !!>> HC 20-12-2021 Change the seed
        if (ihc>wrsize)cycle                                !!>> AL 27-11-24 no sure why cycle instead of do till wrsize...
        idum(ihc)=wridum(ihc)                                
    enddo                                                   

    call random_seed(put=idum)                              !!>> AL 27-11-24 seed random numbers used in this module
    rval=0.0d0                                              
    do ihc=1,wrsize                                         !!>> HC 20-12-2021 CREATE THE RANDOM SEED FOR THE NEXT FILE
        call random_number(rval)                            
        wridum(ihc)=int(rval*1000)                          
    enddo                                                   

    open(837,file="rseed.dat")                              
        write(837,*) wrsize                                  !!>> HC 20-12-2021 Now the children of this individual will have different mutations! !!>> AL 27-11-24 in substitute we rewrite rseed so this is duplicated
        write(837,*) wridum                                  
    close(837)                                              

    print*,"about to choose genes"

    call chose_genes_to_recombinate(individual2)
    print*,"genes chosen"


    call iniread
    
    call readsnap(individual1)

    open(1,file="recombination_info.txt") 
        read(1,*) count                             
    close(1)  

    if(allocated(genes4recomb)) deallocate(genes4recomb)                     
    allocate(genes4recomb(count,ng))
    genes4recomb=0   
    
    if(allocated(es)) deallocate(es)                     
    allocate(es(count,nga))
    es=0   
    
    if(allocated(params)) deallocate(params)                     
    allocate(params(count,4))
    params=0   
    
    if(allocated(lequel)) deallocate(lequel)                     
    allocate(lequel(count))
    lequel=0   
     print*,"here still "
    open(1,file="recombination_info.txt") 
        read(1,*) count                             
        do ihc = 1,count    
            read(1,*) jhc,genes4recomb(ihc,:)
            read(1,*) jhc,es(ihc,:)    
            read(1,*) jhc,params(ihc,1),params(ihc,2),params(ihc,3),params(ihc,4)             
            lequel(ihc) = jhc     
        end do 
    close(1) 
    print*,"now here "
    !do ihc = 1,count  
    !    print*, lequel(ihc)   
    !    print*, genes4recomb(ihc,:)
    !    print*, es(ihc,:)    
    !    print*, params(ihc,1),params(ihc,2),params(ihc,3)          
    !end do 
    
    do ihc = 1,count                            !this is the actual recombination process
        jhc             = lequel(ihc)           !which gene
        gen(jhc)%t(:)   = genes4recomb(ihc,:)
        gen(jhc)%e(:)   = es(ihc,:)  
        gen(jhc)%diffu  = params(ihc,1)
        gen(jhc)%mu     = params(ihc,2)
        gen(jhc)%mich   = params(ihc,3)
        gen(jhc)%kindof = params(ihc,4)   
    end do 
    call iniio
    call writesnap
    
    rename=trim("mv fort.1 "//trim(individual1)) 
    call system(rename)

end subroutine do_recombination_4_evo

subroutine do_recombination_4_evopt1(individual2)

    implicit none
   
    integer         :: count,ihc,jhc,wrsize
    real*8          :: rval
    character*700   :: individual2
    integer, allocatable, dimension(:)  :: lequel 
    real*8, allocatable, dimension(:,:) :: genes4recomb,es,params
    integer, allocatable, dimension(:)  :: wridum

    !RANDOM NUMBER THINGs                                   !!>> AL 27-11-24
    !we will use the file rseed to instantiate random number generation. 
    !This file is the same used in muta_AL.f90, which is generated in seeding and modify in substitute
    call random_seed(size = nseed)                          !!>> AL 27-11-24 this is defined as public variable in general.modf90 
    !call getarg(2,cu)                                       

    if(allocated(idum))deallocate(idum)                     !!>> AL 27-11-24 idum is public variable defined in general and writen to EMaker file in io.f90               
    allocate(idum(nseed))    

    idum=0; wrsize=0                                        
    open(837,file="rseed.dat")                              !!>> HC 20-12-2021 the seed should be stored in the individual directory
        read(837,*) wrsize                                   
        if (allocated(wridum)) deallocate(wridum)            
        allocate(wridum(wrsize))                            
        read(837,*) wridum                                   
    close(837)              

    do ihc=1,size(idum)                                     !!>> HC 20-12-2021 Change the seed
        if (ihc>wrsize)cycle                                !!>> AL 27-11-24 no sure why cycle instead of do till wrsize...
        idum(ihc)=wridum(ihc)                                
    enddo                                                   

    call random_seed(put=idum)                              !!>> AL 27-11-24 seed random numbers used in this module
    rval=0.0d0                                              
    do ihc=1,wrsize                                         !!>> HC 20-12-2021 CREATE THE RANDOM SEED FOR THE NEXT FILE
        call random_number(rval)                            
        wridum(ihc)=int(rval*1000)                          
    enddo                                                   

    open(837,file="rseed.dat")                              
        write(837,*) wrsize                                  !!>> HC 20-12-2021 Now the children of this individual will have different mutations! !!>> AL 27-11-24 in substitute we rewrite rseed so this is duplicated
        write(837,*) wridum                                  
    close(837)                                              

    call chose_genes_to_recombinate(trim(individual2))

end subroutine do_recombination_4_evopt1

subroutine do_recombination_4_evopt2(individual1)

    implicit none
   
    integer         :: count,ihc,jhc,wrsize
    real*8          :: rval
    character*700   :: individual1, rename
    integer, allocatable, dimension(:)  :: lequel 
    real*8, allocatable, dimension(:,:) :: genes4recomb,es,params
    integer, allocatable, dimension(:)  :: wridum

    
    call iniread
    call readsnap(trim(individual1))

    open(1,file="recombination_info.txt") 
        read(1,*) count                             
    close(1)  

    if(allocated(genes4recomb)) deallocate(genes4recomb)                     
    allocate(genes4recomb(count,ng))
    genes4recomb=0   
    
    if(allocated(es)) deallocate(es)                     
    allocate(es(count,nga))
    es=0   
    
    if(allocated(params)) deallocate(params)                     
    allocate(params(count,4))
    params=0   
    
    if(allocated(lequel)) deallocate(lequel)                     
    allocate(lequel(count))
    lequel=0   
    open(1,file="recombination_info.txt") 
        read(1,*) count                             
        do ihc = 1,count    
            read(1,*) jhc,genes4recomb(ihc,:)
            read(1,*) jhc,es(ihc,:)    
            read(1,*) jhc,params(ihc,1),params(ihc,2),params(ihc,3),params(ihc,4)             
            lequel(ihc) = jhc     
        end do 
    close(1) 
    !do ihc = 1,count  
    !    print*, lequel(ihc)   
    !    print*, genes4recomb(ihc,:)
    !    print*, es(ihc,:)    
    !    print*, params(ihc,1),params(ihc,2),params(ihc,3)          
    !end do 
    
    do ihc = 1,count                            !this is the actual recombination process
        jhc             = lequel(ihc)           !which gene
        gen(jhc)%t(:)   = genes4recomb(ihc,:)
        gen(jhc)%e(:)   = es(ihc,:)  
        gen(jhc)%diffu  = params(ihc,1)
        gen(jhc)%mu     = params(ihc,2)
        gen(jhc)%mich   = params(ihc,3)
        gen(jhc)%kindof = params(ihc,4)   
    end do 
    call iniio
    call writesnap
    
    rename=trim("mv fort.1 "//trim(individual1)) 
    call system(rename)

end subroutine do_recombination_4_evopt2

subroutine do_recombination(individual1,individual2)    !AL 10-12-24 important check what happens when different number of real genes

    implicit none
    integer         :: count,ihc,jhc
    character*700   :: individual1, individual2, rename
    integer, allocatable, dimension(:)  :: lequel 
    real*8, allocatable, dimension(:,:) :: genes4recomb,es,params
    
    call chose_genes_to_recombinate(individual2)
    
    call iniread
    call readsnap(individual1)

    open(1,file="recombination_info.txt") 
        read(1,*) count                             
    close(1)  

    if(allocated(genes4recomb)) deallocate(genes4recomb)                     
    allocate(genes4recomb(count,ng))
    genes4recomb=0   
    if(allocated(es)) deallocate(es)                     
    allocate(es(count,nga))
    es=0   
    if(allocated(params)) deallocate(params)                     
    allocate(params(count,4))
    params=0   
    if(allocated(lequel)) deallocate(lequel)                     
    allocate(lequel(count))
    lequel=0   
    
    open(1,file="recombination_info.txt") 
        read(1,*) count                             
        do ihc = 1,count    
            read(1,*) jhc,genes4recomb(ihc,:)
            read(1,*) jhc,es(ihc,:)    
            read(1,*) jhc,params(ihc,1),params(ihc,2),params(ihc,3),params(ihc,4)             
            lequel(ihc) = jhc     
        end do 
    close(1) 
    
    !do ihc = 1,count  
    !    print*, lequel(ihc)   
    !    print*, genes4recomb(ihc,:)
    !    print*, es(ihc,:)    
    !    print*, params(ihc,1),params(ihc,2),params(ihc,3)          
    !end do 
    
    do ihc = 1,count                        !this is the actual recombination process
        jhc = lequel(ihc)                   !which gene
        gen(jhc)%t(:)   = genes4recomb(ihc,:)
        gen(jhc)%e(:)   = es(ihc,:)  
        gen(jhc)%diffu  = params(ihc,1)
        gen(jhc)%mu     = params(ihc,2)
        gen(jhc)%mich   = params(ihc,3)
        gen(jhc)%kindof = params(ihc,4)   
    end do 

    call iniio
    call writesnap
    
    rename=trim(individual1)//"_recombined.dat"
    rename=trim("cp fort.1 "//trim(rename)) 
    call system(rename)
    call system("rm fort.1")
    !call system("rm recombination_info.txt")


end subroutine do_recombination

subroutine create_dummy_icos(individual,num)

    character*700 :: individual
    character*200 :: write
    integer :: num

    call iniread
    call readsnap(individual)

    print*, 'creating dumy: ',individual        
    do ihc = 1,ng
        gen(ihc)%t(:)   = num
        gen(ihc)%e(:)   = num
        gen(ihc)%diffu  = num
        gen(ihc)%mu     = num
        gen(ihc)%mich   = num
        gen(ihc)%kindof   = num
    end do 

    call iniio
    call writesnap
    write=trim("cp fort.1 "//adjustl(individual)) 
    call system(write)
    
end subroutine create_dummy_icos

end module recombination