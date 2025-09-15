program trait_measurement
   !gfortran -w -fexceptions -fno-underscoring -O2 -fbounds-check general.mod.f90 neighboring.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 get_fitness_traitsCS.f90 -o gfit.e
   use io
   use general
   use genetic
   use neighboring

   implicit none 

   logical :: signo
   integer :: ihc,jhc,khc,lhc,contador,contadora
   real*8  :: size_box,maxx,maxy,maxz,totmax,linex,liney,linez,comp,t,dot_PA,dot_vv,mag,dlta,sumdist
   real*8  :: x0(3),x2(3),mini(3),dire(3),projection(3),this(3),closest_projection(3),distances(26)
   real*8  :: traits_morpho(26,3),centroidx,centroidy,centroidz,suma,CS,dis
   character*152 :: path, input, output,input2
   real, allocatable, dimension(:,:)  :: ncoords,cual_caja,orden_final_x,orden_final_y,orden_final_z
   real, allocatable, dimension(:,:)  :: nodos_ep
   integer, allocatable, dimension(:) :: indx
   real*8, dimension(14,3) :: optimum_traits, target_traits

   call getarg(1,input)
   print*,'Analazing file: ', trim(input)
   call getarg(2,input2)
   print*,'Optimum trait proportions: ', trim(input2)
  
   call iniread
   call readsnap(input)
   call iniboxes
   call neighbor_build

   !>>AL 12-9-2024 Centroid
   centroidx = 0;centroidy = 0;centroidz = 0
   do ihc=1,nd
      centroidx = centroidx + node(ihc)%x
      centroidy = centroidy + node(ihc)%y
      centroidz = centroidz + node(ihc)%z
   enddo
   
   centroidx = centroidx/nd
   centroidy = centroidy/nd
   centroidz = centroidz/nd
   
   !centering morpho around cero
   do ihc=1,nd
      node(ihc)%x = node(ihc)%x - centroidx
      node(ihc)%y = node(ihc)%y - centroidy
      node(ihc)%z = node(ihc)%z - centroidz
   end do 

   dlta=nodeo(1)%add*2                             !AL 19-3-25 nodeo%add is node%add at dev. time = 0, so it is always the same
   traits_morpho = 0
   contadora = 0
   do ihc = 1,-1,-1
      do jhc = 1,-1,-1
         do khc = 1,-1,-1
         if(ihc==0 .and. jhc==0 .and. khc==0)cycle
            contadora = contadora + 1
            !with the point in x2 we generate the trait direction vector, which is (x2 - 0)

            x2(1) = ihc 
            x2(2) = jhc
            x2(3) = khc

            contador = 0

            dire = x2/sqrt(sum(x2**2)) !direction of trait vector 
            do lhc=1,nd                                  !Count how many apical epithelial nodes there are
               if(node(lhc)%tipus.ne.1)cycle             !only consider apical epithelial nodes

               x0(1) = node(lhc)%x; x0(2) = node(lhc)%y; x0(3) = node(lhc)%z

               if(ihc .ne. 0 ) then                       !only consider nodes in the same octant as direction vector
                  signo = (sign(x0(1), x2(1)) .eq. x0(1)) !this is true if x0 and x1 have same sign
                  if(.not. signo)cycle
               end if
               if(jhc .ne. 0 ) then
                  signo = (sign(x0(2), x2(2)) .eq. x0(2))
                  if(.not. signo)cycle
               end if
               if(khc .ne. 0 ) then
                  signo = (sign(x0(3), x2(3)) .eq. x0(3))
                  if(.not. signo)cycle
               end if
                
               t = dot_product(dire,x0)                   !distance between 0 and projection of x0 on the line x2 (scalar projection)
               projection = t*dire                        !projection vector of point x on trait direction, the magnitud will be trait value in case of closest and furthest point to trait line
               this = projection - x0
               
               if(sqrt(sum(this**2)) > dlta) cycle       !magnitude of vector from projection vector to node coordinate vector (how far away from the direction is the node in consideration)
               contador = contador + 1
               if(contador == 1) then
                  mini(1) = node(lhc)%x
                  mini(2) = node(lhc)%y
                  mini(3) = node(lhc)%z
                  comp = sqrt(sum(projection**2))
                  closest_projection = projection
               else 
                  mag = sqrt(sum(projection**2))
                  if(mag > comp)then
                     comp = mag
                     mini(1) = node(lhc)%x
                     mini(2) = node(lhc)%y
                     mini(3) = node(lhc)%z
                     closest_projection = projection
                  end if
               end if

            end do
            traits_morpho(contadora,1) = closest_projection(1)
            traits_morpho(contadora,2) = closest_projection(2)
            traits_morpho(contadora,3) = closest_projection(3)
         end do 
      end do 
   end do 


   optimum_traits(1,:)  = traits_morpho(13,:)     !(0,0,+) D1 
   optimum_traits(2,:)  = traits_morpho(11,:)     !(0,+,0) D2
   optimum_traits(3,:)  = traits_morpho(14,:)     !(0,0,-) D3 
   optimum_traits(4,:)  = traits_morpho(16,:)     !(0,-,0) D4
   optimum_traits(5,:)  = traits_morpho(1,:)      !(+,+,+) D5
   optimum_traits(6,:)  = traits_morpho(3,:)      !(+,+,-) D6
   optimum_traits(7,:)  = traits_morpho(5,:)      !(+,0,0) D7
   optimum_traits(8,:)  = traits_morpho(9,:)      !(+,-,-) D8
   optimum_traits(9,:)  = traits_morpho(7,:)      !(+,-,+) D9
   optimum_traits(10,:) = traits_morpho(18,:)     !(-,+,+) D10
   optimum_traits(11,:) = traits_morpho(20,:)     !(-,+,-) D11
   optimum_traits(12,:) = traits_morpho(22,:)     !(-,0,0) D12
   optimum_traits(13,:) = traits_morpho(26,:)     !(-,-,-) D13
   optimum_traits(14,:) = traits_morpho(24,:)     !(-,-,+) D14

   !calculate centroid size 
   suma = 0
   do ihc=1,14
      suma = suma + sum(optimum_traits(ihc,:)**2)
   end do

   CS = sqrt(suma)

   !normalize
   optimum_traits = optimum_traits/CS
   open(10, file=trim(input2)) !this distances should be already normalized
      do ihc = 1, 14
        read(10, *) target_traits(ihc, 1), target_traits(ihc, 2), target_traits(ihc, 3)  ! Read each row
      end do
   close(10)
   
   dis = 0
   do ihc = 1, 14
      dis = dis + sqrt(sum((target_traits(ihc,:) - optimum_traits(ihc,:))**2))
   end do

   open (123, file="fitnes.dat", status="replace", action="write")
      write(123,*) dis
   close(123)

end program trait_measurement