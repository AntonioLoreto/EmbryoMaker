! Reads individual fitnesses from stdin and prints out the selected
! parent indices (start from 0) to stdout.
!
! Build: gfortran -o select select.f90
!
! Use: ./select < fitnesses.dat > parents.dat
!
! Probability of choosing an individual is linearly proportional to it's
! fitness. It is easy to change in read_cumulative_sum(), for example.

program select
  implicit none

  integer, parameter :: m = 10000 !2048  ! Maximum array size !!>> HC 3-2-2020

  integer :: n,i         ! Number of individuals
  real*8    :: s(m)      ! Cumulative sum of fitnesses
  real*8    :: fit(m)      ! Cumulative sum of fitnesses
  
!  call initialise_randnr
!  call read_cumulative_sum(n, s)
!  call write_sample_indices(n, s)
   call best_first
!   call best_half
   
contains

  subroutine initialise_randnr
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    
    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
  
  end subroutine initialise_randnr


  subroutine read_cumulative_sum(n, s)
    integer, intent(out) :: n
    real*8, intent(out)    :: s(:)
    real*8                 :: f, sf
    integer              :: i
    n = 0
    sf = 0.0
    do
       read(UNIT = *, FMT = *, END = 11) f
       n = n + 1
       if (n > m) then
          error stop 22
       end if
       sf = sf + f
       s(n) = sf
    end do
    
11  continue
  end subroutine read_cumulative_sum

  subroutine write_sample_indices(n, s)
    integer, intent(in) :: n
    real*8, intent(in)    :: s(:)
    real*8, allocatable   :: r(:)
    integer             :: ir,g

    allocate(r(n))
    call random_number(r)
    do ir = 1, n
       write(*,'(I0)') linear_search(n, r(ir)*s(n), s)
    end do
  
  end subroutine write_sample_indices

  function linear_search(n, v, s) result(il)
    integer, intent(in) :: n
    real*8, intent(in)    :: v
    real*8, intent(in)    :: s(:)
    integer             :: il
    il = 0
    do
       if (v .le. s(il+1)) exit
       il = il + 1
       if (il .gt. n) error stop 33
    end do
  end function linear_search
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine best_first                     !!>> HC 5-10-2021 New selection algorithm written by Hugo Cano
    real*8                 :: f             !!>> HC 5-10-2021
    real*8, allocatable   :: r(:), rf(:)    !!>> HC 5-10-2021
    integer, allocatable  :: par(:)         !!>> HC 5-10-2021
    integer              :: i, j, ir        !!>> HC 5-10-2021
    real*8                :: a, b, tf       !!>> HC 5-10-2021
    n = 0                                   !!>> HC 5-10-2021
    do                                      !!>> HC 5-10-2021
       read(UNIT = *, FMT = *, END = 11) f  !!>> HC 5-10-2021  READ INPUT FITNESSES FROM TERMINAL
       n = n + 1                            !!>> HC 5-10-2021 counter of the number of individuals
       if (n > m) then                      !!>> HC 5-10-2021 execeded max num of individuals
          error stop 22                     !!>> HC 5-10-2021
       end if                               !!>> HC 5-10-2021
       s(n) = f                             !!>> HC 5-10-2021 save fitnesses
    end do                                  !!>> HC 5-10-2021
    
    11 continue                             !!>> HC 5-10-2021 finished reading
        
    allocate(r(1:n))                        !!>> HC 5-10-2021 these arrays are to store fitnesses and
    allocate(rf(1:n))                       !!>> HC 15-11-2021
    allocate(par(1:n))                      !!>> HC 5-10-2021 order of parents
    r=s(1:n)                                !!>> HC 5-10-2021
    par=0                                   !!>> HC 5-10-2021
    b=1.0d0/real(n)                         !!>> HC 15-11-2021 Proportion represented by 1 individual
    tf=0.0d0; rf=0.0d0                      !!>> HC 15-11-2021
    do i=1,n                                !!>> HC 15-11-2021
       tf=tf+r(i)                           !!>> HC 15-11-2021 Sum of absoulute fitness
    enddo                                   !!>> HC 15-11-2021
    do i=1,n                                !!>> HC 15-11-2021
       rf(i)=r(i)/tf                        !!>> HC 15-11-2021 Relative fitness
    enddo                                   !!>> HC 15-11-2021
    
    do i=1,n                                !!>> HC 5-10-2021 SOMEONE IS SELECTING THE CHILDREN HERE
       a=maxval(rf(1:n))                    !!>> HC 5-10-2021 maximum fitness
       do j=1,n                             !!>> HC 5-10-2021
          if (rf(j)==a)then                 !!>> HC 5-10-2021 look for the individual with maximum fitness
             par(i)=j-1                     !!>> HC 5-10-2021 Give a child to this parental
             rf(j)=rf(j)-b                  !!>> HC 5-10-2021 remove 1 to the fitness of the parental
          endif                             !!>> HC 5-10-2021
       enddo                                !!>> HC 5-10-2021
    enddo                                   !!>> HC 5-10-2021
    
    do ir = 1, n                            !!>> HC 5-10-2021 WRITING OUTPUT TO TERMINAL
       write(*,*) par(ir)                   !!>> HC 5-10-2021 for evolve.sh to read it
    end do                                  !!>> HC 5-10-2021
    
  end subroutine best_first                 !!>> HC 5-10-2021

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine best_half 
    real*8                 :: f             !!>> HC 5-10-2021
    real*8, allocatable   :: r(:)           !!>> HC 5-10-2021
    integer, allocatable  :: par(:),old(:)  !!>> HC 5-10-2021
    integer              :: i, j, ir, k     !!>> HC 5-10-2021
    real*8                :: a, b           !!>> HC 5-10-2021
    n = 0                                   !!>> HC 5-10-2021
    do                                      !!>> HC 5-10-2021
       read(UNIT = *, FMT = *, END = 11) f  !!>> HC 5-10-2021  READ INPUT FITNESSES FROM TERMINAL
       n = n + 1                            !!>> HC 5-10-2021 counter of the number of individuals
       if (n > m) then                      !!>> HC 5-10-2021 execeded max num of individuals
          error stop 22                     !!>> HC 5-10-2021
       end if                               !!>> HC 5-10-2021
       s(n) = f                             !!>> HC 5-10-2021 save fitnesses
    end do                                  !!>> HC 5-10-2021
    
    11 continue                             !!>> HC 5-10-2021 finished reading  

    allocate(r(1:n))                        !!>> HC 5-10-2021 these arrays are to store fitnesses and
    allocate(par(1:n), old(1:n))            !!>> HC 5-10-2021 order of parents
    r=s(1:n)                                !!>> HC 5-10-2021
    par=0; old=0                            !!>> HC 5-10-2021  
    
    do i=1,n                                !!>> HC 6-10-2021  ORDER FITNESS DATA
       a=0.0d0                              !!>> HC 6-10-2021
       do j=1,n                             !!>> HC 6-10-2021
          if (r(j)>a)then                   !!>> HC 6-10-2021 
             a=r(j)                         !!>> HC 6-10-2021 a stores the maximum value
             old(i)=j                       !!>> HC 6-10-2021 old(i) stores the old index in vector r
          endif                             !!>> HC 6-10-2021
       enddo                                !!>> HC 6-10-2021
       r(old(i))=0.0d0                      !!>> HC 6-10-2021 remove fitness in vector r 
    enddo                                   !!>> HC 6-10-2021 so we do not find the same individual twice
    
    ir=0;k=1                                !!>> HC 6-10-2021 SOMEONE IS GIVING CHILDREN HERE
    do i=1,n                                !!>> HC 6-10-2021 go through all the new individuals and
       ir=ir+1                              !!>> HC 6-10-2021 give 2 children to the best half of parentals
       if (ir>2)then                        !!>> HC 6-10-2021 
          ir=1                              !!>> HC 6-10-2021
          k=k+1                             !!>> HC 6-10-2021 this is the parental
       endif                                !!>> HC 6-10-2021
       par(i)=old(k)                        !!>> HC 6-10-2021
    enddo                                   !!>> HC 6-10-2021
    
    do ir = 1, n                            !!>> HC 5-10-2021 WRITING OUTPUT TO TERMINAL
       write(*,*) par(ir)                   !!>> HC 5-10-2021 for evolve.sh to read it
    end do                                  !!>> HC 5-10-2021
  
  end subroutine best_half

end program select
