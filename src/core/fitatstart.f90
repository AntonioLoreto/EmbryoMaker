!TODO: automaticamente elegir cual funcion de fitness dependiendo de paraopc.par !AL 25-9-25

program proj_fit
! gfortran general.mod.f90 genetic.mod.f90  io.mod.f90  geompack3.f90 neighboring.mod.f90 fitatstart.f90 -o fitatstar.e
   
   use io
   use general
   use genetic
   use neighboring
   use analysis

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

   !print*,"file to analize fit",trim(input)
   call iniread
   call readsnap(trim(input))
   call iniboxes
   call neighbor_build
   
   call traits_intersection(fitelli)

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