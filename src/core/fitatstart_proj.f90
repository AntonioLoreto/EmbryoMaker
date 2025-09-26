program proj_fit
! gfortran general.mod.f90 genetic.mod.f90  io.mod.f90  geompack3.f90 neighboring.mod.f90 fitatstart_proj.f90 -o fitatstartproj.e
   
   use io
   use general
   use genetic
   use neighboring

   implicit none 
   logical :: signo
   integer :: ihc,jhc,khc,lhc,contador,contadora
   real*8  :: size_box,maxx,maxy,maxz,totmax,linex,liney,linez,comp,t,dot_PA,dot_vv,mag,dlta,sumdist
   real*8  :: x0(3), x2(3),mini(3),dire(3),projection(3),this(3),closest_projection(3),distances(14)
   real*8  :: traits_morpho(26,3), fitelli,emd,centroidx,centroidy,centroidz
   character*1000 :: path, input, output, input2
   real, allocatable, dimension(:,:)  :: ncoords,cual_caja,orden_final_x,orden_final_y,orden_final_z
   real, allocatable, dimension(:,:)  :: nodos_ep
   integer, allocatable, dimension(:) :: indx
   real*8, dimension(14) :: target_distances,diff_distances

   call getarg(1,input)
   call getarg(2,input2)

   print*,"file to analize fit",trim(input)
   call iniread
   call readsnap(trim(input))
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
   
   !Generate the lines wich represent trait directions 
   maxx=0.0d0; maxy=0.0d0; maxz=0.0d0
   do ihc=1,nd
      if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
      if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
      if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
   enddo

   totmax=0.0d0
   totmax=maxx
   if (maxy>totmax) totmax=maxy
   if (maxz>totmax) totmax=maxz

   totmax = totmax + 0.1*totmax !with this value we'll generate the trait lines. They will be always longer than the morphology
   dlta=nodeo(1)%add*2
   traits_morpho = 0
   contadora = 0
   do ihc = 1,-1,-1
      do jhc = 1,-1,-1
         do khc = 1,-1,-1
         if(ihc==0 .and. jhc==0 .and. khc==0)cycle
            contadora = contadora + 1
            !with the point in x2 we generate the trait direction vector, which is (x2 - 0)
            x2(1)= totmax*ihc 
            x2(2) = totmax*jhc
            x2(3) = totmax*khc

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
                
               t = dot_product(dire,x0)                   !distance between 0 and projection of x0 on the line x2 (dot product)
               projection = t*dire                        !projection vector of point x on trait direction, the magnitud will be trait value in case of closest point to trait line
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

   distances = 0
   distances(1) = sqrt(sum((traits_morpho(13,:)**2)))
   distances(2) = sqrt(sum((traits_morpho(11,:)**2)))
   distances(3) = sqrt(sum((traits_morpho(14,:)**2)))
   distances(4) = sqrt(sum((traits_morpho(16,:)**2)))
   distances(5) = sqrt(sum((traits_morpho(1,:)**2)))
   distances(6) = sqrt(sum((traits_morpho(3,:)**2)))
   distances(7) = sqrt(sum((traits_morpho(5,:)**2)))
   distances(8) = sqrt(sum((traits_morpho(9,:)**2)))
   distances(9) = sqrt(sum((traits_morpho(7,:)**2)))
   distances(10) = sqrt(sum((traits_morpho(18,:)**2)))
   distances(11) = sqrt(sum((traits_morpho(20,:)**2)))
   distances(12) = sqrt(sum((traits_morpho(22,:)**2)))
   distances(13) = sqrt(sum((traits_morpho(26,:)**2)))
   distances(14) = sqrt(sum((traits_morpho(24,:)**2)))
   
   sumdist = 0
   do ihc=1,14
      sumdist = sumdist + distances(ihc)
   end do 

   distances = distances/sumdist

   open(10, file=trim(input2)) !this distances should be already normalized
      read(10,*) target_distances
   close(10)

   EMD = (sum(abs(target_distances(:) - distances(:))))

   fitelli = 1.0d0-EMD

   open (123, file="fitnessatstart.dat")
      write(123,*) fitelli
   close(123)

   !open (124, file="individual.volume.txt")
   !   write(124,*) 0.0
   !close(124)

   !Checar problema de necesitar un fichero llamado individual.volume.txt. Escribirlo al calcular la asimetria bilateral.
   !open (123, file="individual.volume.txt")
   !   write(123,*) volume                 !!>>AL 4-4-24: If no displacement, volume(1)=volume(2)=volume(3)
   !close(123)

   !open(616,file="rob.val") !checar por que necesitas esto. nuevamente, podria solo escribirlo en la parta de asimetria bilateral
   !   write(616,*) 0.00111   
   !close(616)
    
end program proj_fit