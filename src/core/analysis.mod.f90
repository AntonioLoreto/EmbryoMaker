!    EmbryoMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Created by Renske 24-04-2018
! This module contains functions for analysing embryonic development. They could be called by auto in automaticon, 
! which would however require one to divide an elli run into subruns (otherwise there is not much to analyse)

module analysis

  use general
  use genetic
  use io
  use neighboring

  !parameters
  implicit none 
  real*8, allocatable :: storprev(:,:),stornow(:,:) !to generate a centered tissue
  real*8, allocatable ::maxgex(:,:) !store max expression per gene (and the time of maximum
  integer ::epi_prev, epi_now
  integer :: nrnods !the previous nr of nodes
  integer :: grofile=942
  integer :: OPCdistfile=337
  integer :: OPCfile=338
  integer :: gefile1=863  !file handle
  integer :: gefile2=864  !file handle
  integer :: gefile3=865 !file handle
  
  !************************************
  contains  
  ! gathers the several subroutines in one handy bundle
  subroutine do_analysis(which)
    integer :: which
    
    call fill_node_arrays
    
    if (which==1 .or. which==3) then 
      call analyse_growth
      call OPCval
    end if
    if (which==2 .or. which==3) then 
      call analyse_gene_expression(1)
      call analyse_gene_expression(2)
      call maxexpr
    else if (which/=1) then
      print *, "Error in analysis.mod: invalid option. Please use 1, 2 or 3. Exiting..."
      stop
    end if
     
    !copy the situation now for later 
    !if (which==1 .or. which==3) call copy_nodes
    
  end subroutine

  !to save the initial conditions
  subroutine initialise_analysis(filestart,which)
    integer which
    character*300 filestart
    
    print *, "starting analysis"  
    call fill_node_arrays  
      
    if (which==1 .or. which==3) then !analyse tissue growth and shape
      open(grofile, file=trim(filestart)//"growthprofile.dat")
      write (grofile, *) 0, nd, 0.0, 0, 0, 0, 0, 0, 0 
      open(OPCdistfile, file=trim(filestart)//"patchsizedist.dat")
      open(OPCfile, file=trim(filestart)//"patchcount.dat")
      call OPCval
      call copy_nodes
    end if
    
    if (which==2 .or. which==3) then !analyse gene expression
      !call copy_gex
      open(gefile1, file=trim(filestart)//"genexpression1.dat")
      call analyse_gene_expression(1)
      open(gefile2, file=trim(filestart)//"genexpression2.dat")
      call analyse_gene_expression(2)
      
      !store max expression
      open(gefile3,file=trim(filestart)//"maxgex.dat")
      if (allocated(maxgex)) deallocate(maxgex)
      allocate(maxgex(ng,2))
      maxgex=0.0d0
      call maxexpr
      
    else if (which/=1) then
      print *, "Error in analysis.mod: invalid option. Please use 1, 2 or 3. Exiting..."
      stop
    end if
  end subroutine 
  
  subroutine copy_nodes
    if(allocated(storprev)) deallocate(storprev)
    allocate(storprev(epi_now,6))
    storprev=stornow
  end subroutine
  
  !close the files that were written into, delete the used memory
  subroutine close_analysis(which)
    integer which
    print *, "ending analysis"
    if (which==1 .or. which==3) then
      close(grofile)
      close(OPCfile)
      close(OPCdistfile)
      !if (allocated(prevnods)) deallocate(prevnods)
    end if
    if (which==2 .or. which==3) then 
      close(gefile1)
      close(gefile2)
      !if (allocated(prevgex)) deallocate(prevgex)
      call writemax
      close(gefile3)
    end if
    
    if (allocated(storprev)) deallocate(storprev)
    if (allocated(stornow)) deallocate(stornow)
    if (allocated(maxgex)) deallocate(maxgex)
    
  end subroutine
  

  subroutine maxexpr
    
    do i=1,nd
      do j=1, ng
        if (gex(i,j)>maxgex(j,1)) then
          maxgex(j,1)=gex(i,j)
          maxgex(j,2)=rtime
        endif
      end do
    end do
  
  end subroutine
  
  subroutine writemax
    print *, "here"
    write(gefile3, *) "gene"//char(9)//"maxex"//char(9)//"timept"
    do j=1, ng
      write(gefile3, '(A,es10.3,es10.1)')gen(j)%label, maxgex(j,1), maxgex(j,2)
    end do
    
  end subroutine
  
  !new version: gotta deal with the fact that nodes are in different layers
  subroutine analyse_gene_expression(layer)
    integer::layer !are we looking at apical or basal layer?
    !real*8, allocatable linegex(:,:,:) !to store the expression on a line in a nr of directions
    integer :: i,j, gefile
    
    !get the right file handle
    if (layer==1) gefile=gefile1
    if (layer==2) gefile=gefile2
    
    
    !x is 0, y and z vary
    do i=1, epi_now
    !!.and. stornow(i,2)>(miny_x+cell-stornow(i,4)*0.5) & .and. stornow(i,2)<(miny_x+cell + stornow(i,4)*0.5)
      if (stornow(i,5)==layer .and. stornow(i,1)<stornow(i,4)*1.5 .and. stornow(i,1)>-stornow(i,4)*1.5  ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    
    !y is 0, x and z vary
    do i=1, epi_now
    !.and. stornow(i,1)>minx_y+cell-stornow(i,4)*0.5 & .and. stornow(i,1)<minx_y+cell+stornow(i,4)*0.5
      if(stornow(i,5)==layer .and. stornow(i,2)<stornow(i,4)*1.5 .and. stornow(i,2)>-stornow(i,4)*1.5 ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    
    !z is 0, x and y vary
    do i=1, epi_now
    !.and. stornow(i,1)>minx_z+cell-stornow(i,4)*0.5 & .and. stornow(i,1)<minx_z+cell+stornow(i,4)*0.5
      if(stornow(i,5)==layer .and.  stornow(i,3)<stornow(i,4)*1.5 .and. stornow(i,3)>-stornow(i,4)*1.5 ) then 
        write(gefile, *) rtime, stornow(i,1), stornow(i,2), stornow(i,3), gex(int(stornow(i,6)),:)
        !cell=cell+stornow(i,4)
      end if
    end do
    write(gefile, *)
    write(gefile, *)
    
  end subroutine
  
  
  subroutine analyse_growth
    real*8 :: d, dd, aa, a, b, c, req1, req2 !for emd
    integer :: rostral1, caudal1, dorsal1, ventral1, left1, right1 !to assess growth or displacement 
    integer :: rostral2, caudal2, dorsal2, ventral2, left2, right2 !in different areas
    integer :: i,j
    integer :: patchcount
      !now EMD (this does not seem to be working too well: always prints 1000000)
      
      d=0.0d0
      do i=1,epi_prev
        ! if(node(i)%tipus>2)cycle
        a=storprev(i,1) ; b=storprev(i,2) ; c=storprev(i,3) ; req1=storprev(i,4)
        aa=1000000
        do j=1,epi_now
          dd=((sqrt((stornow(j,1)-a)**2+(stornow(j,2)-b)**2+(stornow(j,3)-c)**2)))
          if (dd<aa) aa=dd
        end do          
        d=d+aa
      end do
      do i=1,epi_now
        a=stornow(i,1) ; b=stornow(i,2) ; c=stornow(i,3) ; req2=stornow(i,4)
        aa=1000000
        do j=1,epi_prev
          dd=((sqrt((storprev(j,1)-a)**2+(storprev(j,2)-b)**2+(storprev(j,3)-c)**2)))
          if (dd<aa) aa=dd
        end do          
        d=d+aa
      enddo
      d=d/real(epi_prev+epi_now)
    
      !use the centered tissues to identify where more growth has occurred: use the fact that new nodes are added at the end of the array!
      rostral1=0
      caudal1=0
      dorsal1=0
      ventral1=0
      left1=0
      right1=0
      do i=epi_prev+1, epi_now
        if (stornow(i,1)>0.0) then ; right1=right1+1; else; left1=left1+1; end if
        if (stornow(i,2)>0.0) then ; dorsal1=dorsal1+1; else; ventral1=ventral1+1; end if
        if (stornow(i,3)>0.0) then ; rostral1=rostral1+1; else; caudal1=caudal1+1; end if  
      end do
      
!       !use the centered tissues to identify where more growth has occurred
!       rostral2=0
!       caudal2=0
!       dorsal2=0
!       ventral2=0
!       left2=0
!       right2=0
!       do i=1, epi_now
!         if (stornow(i,1)>0.0) then ; right2=right2+1; else; left2=left2+1; end if
!         if (stornow(i,2)>0.0) then ; dorsal2=dorsal2+1; else; ventral2=ventral2+1; end if
!         if (stornow(i,3)>0.0) then ; rostral2=rostral2+1; else; caudal2=caudal2+1; end if  
!       end do
    
      
      write (grofile, *) rtime, nd, d, rostral1,caudal1,dorsal1,ventral1,left1,right1
    
     
  end subroutine

  
  !The arrays with the orientation of each node, and their neighbours, have already been filled in fill_node_arrays
  subroutine OPCval
    integer :: patch_count
    integer, allocatable :: checks(:), store1(:),store2(:)
    integer, allocatable :: orientation(:), epi_node_list(:),translate(:) !to store the epithelium orientation
    integer::patch_size,nstore1,nstore2,min_patch_size,max_patch, maxmaxpatch 
    integer,allocatable::epi_neigh_list(:,:), patchsize_dist(:)
    integer :: nepi, maxneigh !nr of epithelial cells (by counting only apical nodes), and nr of neighbouring cells
    real*8 :: a, b, c, d
    integer :: i, j, ii, jj, jjj, k, kkk
    
    maxmaxpatch=500000000 !!>>HC 3-9-2020 with large polyps we can have a large OPC sizes
    
    max_patch=0
    min_patch_size=25
    !print *, "starting OPC"
    
    !allocating memory
    if (allocated(orientation))deallocate(orientation)
    allocate(orientation(ncels))

    if (allocated(epi_node_list))deallocate(epi_node_list)
    allocate(epi_node_list(ncels))
    
    if (allocated(store1))deallocate(store1)
    allocate(store1(ncels))

    if (allocated(store2))deallocate(store2)
    allocate(store2(ncels))
    
    if (allocated(checks))deallocate(checks)
    allocate(checks(ncels))    

    if (allocated(epi_neigh_list)) deallocate(epi_neigh_list)
    maxneigh=maxval(nneigh) !; print*,"maxneigh",maxneigh
    allocate(epi_neigh_list(ncels,maxneigh))
    
    if (allocated(patchsize_dist))deallocate(patchsize_dist)
    allocate(patchsize_dist(ncels))    
    
    if (allocated(translate))deallocate(translate)
    allocate(translate(nda))
      
    nepi=0
    do i=1,nd
      if (node(i)%tipus==2) then
        nepi=nepi+1
        epi_node_list(nepi)=i
        translate(i)=nepi
        ii=0
        ii=node(i)%altre
        a=node(i)%x-node(ii)%x  !subtraction gives the vector spanning the epithelial cylinder
        b=node(i)%y-node(ii)%y
        c=node(i)%z-node(ii)%z
        d=1/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
           
       !we determine the orientation of the vector  
       if(c>=-0.95d0.and.c<=-epsilod)then  
         if(b<=-epsilod)then
           if(a<=-epsilod)then
             orientation(nepi)=1
           else
             orientation(nepi)=2
           end if
         else
           if(a<=-epsilod)then
             orientation(nepi)=3
           else
             orientation(nepi)=4
           end if
         end if
       elseif(c<=0.95d0.and.c>=epsilod)then
         if(b<=-epsilod)then
           if(a<=-epsilod)then
             orientation(nepi)=5
           else
             orientation(nepi)=6
           end if
         else
           if(a<=-epsilod)then
             orientation(nepi)=7
           else
             orientation(nepi)=8
           end if
         end if
       else
         orientation(nepi)=9
       end if !end of orientation determination
         
     end if !end of apical nodes conditional
        
    enddo
    
    epi_neigh_list=0
    !fill list of neighbouring nodes
    do ii=1,nepi
      jj=0
      i=epi_node_list(ii)
      do j=1,nneigh(i)
        if(node(neigh(i,j))%tipus==2)then
          jj=jj+1
          epi_neigh_list(ii,jj)=translate(neigh(i,j))
        end if 
      end do
    end do
   
    !count patches and record distribution
    checks=0
    patch_count=0
    patchsize_dist=0
    do i=1,nepi
    ii=epi_node_list(i)
    if(checks(i)==1) cycle
    checks(i)=1
    if(orientation(i)==0) cycle !it has not proper orientation
    !we start a patch
    store1=0 ; store2=0 ; nstore1=0 ; nstore2=0
    patch_size=1
    do j=1,maxneigh
      k=epi_neigh_list(i,j) 
      if(k==0) exit
    !  print*,"i",i,"j",j,"k",k
      if(checks(k)==1) cycle
      if(orientation(i)/=orientation(k))cycle
      patch_size=patch_size+1
      nstore2=nstore2+1
      store2(nstore2)=k
      checks(k)=1
    end do
    do while (nstore2/=0)
      store1=store2 ; store2=0
      nstore1=nstore2 ; nstore2=0
      do j=1,nstore1
        jjj=store1(j)
        do k=1,maxneigh
          kkk=epi_neigh_list(jjj,k)
          if(kkk==0) exit
          if(checks(kkk)==1) cycle
          if(orientation(jjj)/=orientation(kkk))cycle
          nstore2=nstore2+1
          patch_size=patch_size+1
          store2(nstore2)=kkk
          checks(kkk)=1
        end do
      end do
    end do
 
    if(patch_size>=min_patch_size)then 
      patch_count=patch_count+1
      patchsize_dist(patch_size)=patchsize_dist(patch_size)+1
      if(patch_size>max_patch) max_patch=patch_size
    endif
    
  end do
  
  if (max_patch>maxmaxpatch) print *, "OPC warning: exceeding maximum defined patch size"
  
  
  !write the size distribution to a file
  do i=min_patch_size,max_patch
    if (patchsize_dist(i)>0) then
      do j=1,patchsize_dist(i)
        write(OPCdistfile, fmt="(es10.2,I4)") rtime, i
      end do
    end if
  end do
!     do i=max_patch+1,maxmaxpatch
!       write(OPCdistfile, fmt="(I8,I4,I4)") timepoint, i, 0 
!     end do
!     write(OPCdistfile, *) ""

!    print*, "OPC", patch_count  !!>> HC 30-11-2020

      open(202,file="individual.datfitness")                              !!>> HC 11-11-2021 WRITING OUTPUT
      write(202,*) real(patch_count)                                      !!>> HC 11-11-2021
      close(202)                                                          !!>> HC 11-11-2021
      open(203,file="OPC.val")                                            !!>> HC 11-11-2021
      write(203,*) real(patch_count)                                      !!>> HC 11-11-2021

    
    write(OPCfile, fmt="(es10.2,I4)") rtime, patch_count 
    
    !deallocating memory, this important, otherwise mem leaks may occur!
    if (allocated(orientation))deallocate(orientation)

    if (allocated(epi_node_list))deallocate(epi_node_list)
     
    if (allocated(store1))deallocate(store1)
   
    if (allocated(store2))deallocate(store2)
      
    if (allocated(checks))deallocate(checks)
   
    if (allocated(epi_neigh_list)) deallocate(epi_neigh_list)
       
    if (allocated(patchsize_dist))deallocate(patchsize_dist)
        
    if (allocated(translate))deallocate(translate)
   
    
  end subroutine
  
    
  subroutine fill_node_arrays
    real*8 :: centx, centy, centz
    integer ::epnd
    integer :: i
     
    
    !store the node positions of the current time point
    epi_now=0
    do i=1,nd
      if(node(i)%tipus>2)cycle  
      epi_now=epi_now+1
    enddo
    
    if (epi_now>0) then 
    !now we transfer the coordinates of the individual to matrix storprev
      if (allocated(stornow)) deallocate(stornow)
      allocate(stornow(epi_now,6))
      epnd=0
      
      do i=1,nd
        if(node(i)%tipus>2)cycle  
        epnd=epnd+1
        stornow(epnd,1)=node(i)%x
        stornow(epnd,2)=node(i)%y
        stornow(epnd,3)=node(i)%z
        stornow(epnd,4)=node(i)%eqd  !!>>HC 30-6-2020
        stornow(epnd,5)=node(i)%tipus
        stornow(epnd,6)=i
      
      end do
      
    endif
  
    centx=sum(stornow(:,1))/epnd     !centroid
    centy=sum(stornow(:,2))/epnd
    centz=sum(stornow(:,3))/epnd
   
    do i=1,epnd
      stornow(i,1)=stornow(i,1)-centx   !centering 
      stornow(i,2)=stornow(i,2)-centy
      stornow(i,3)=stornow(i,3)-centz
    end do
  
 
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine joint_entropy_fourth_original_neighbors(fitelli)  !!>> HC 28-9-2021 This subroutine was mady by Hugo Cano and it calculates complexity
                                                      !!>> HC 28-9-2021 as the joint entropy of the curvature classes 
                                                      !!>> HC 28-9-2021 of the original nodes of the embryo with their 4th original neighbors
                                                      !!>> HC 28-9-2021 We assume that there is not apoptosis
                                                      !!>> HC 28-9-2021 We need to know the original number of apical nodes and
                                                      !!>> HC 28-9-2021 the original number of  relationships node-neighbor
                                                      !!>> HC 28-9-2021 These relationships must also be written in two files in the working directory
   implicit none                                      !!>> HC 28-9-2021 fourth_original_neighs.dat for the tipe 1 nodes
                                                      !!>> HC 28-9-2021 fourth_original_neighs2.dat for the tipe 2 nodes
   integer, parameter :: finitnn = 8280               !!>> HC 28-9-2021 Initial number of association node-4th neighbour
   integer, parameter :: napical = 372                !!>> HC 28-9-2021 Number of original nodes
   real*8 :: angle, upv, modu, modv, sumangle, numn, medianglee, entropy, entropy2, entropy1 !!>> HC 28-9-2021
   real*8 ::  pij,  classtemsum, rnapical, totalneightshc                                    !!>> HC 28-9-2021
   integer ::  ihc, jhc, khc, tiponodhc, order, order2, node2_1                            !!>> HC 28-9-2021
   integer:: tempclass, lowlimit, highlimit, steps, scaling, nclasseshc                    !!>> HC 28-9-2021
   integer, dimension(1:napical) :: classangles, classtem, oldorder                        !!>> HC 28-9-2021
   real*8, dimension(1:3) :: u, v                                                            !!>> HC 28-9-2021
   real*8,  dimension(1:napical) :: avangles, probclass                                      !!>> HC 28-9-2021
   real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs, classprod                 !!>> HC 28-9-2021
   integer, dimension(1:finitnn) :: foldnodes, foldnodesneight, angles1, angles2           !!>> HC 28-9-2021
   real*8 :: fitelli                                                                       !!>> HC 11-11-2021
   
   call neighbor_build                                                                     !!>> HC 28-9-2021



   lowlimit = -10                                            !!>> HC 28-9-2021 Lower limit of the interval in classes
   highlimit= 9                                              !!>> HC 28-9-2021 Higher limit of the interval in classes
   steps = 1                                                 !!>> HC 28-9-2021 width of the interval in classes
   nclasseshc = 20                                           !!>> HC 28-9-2021 Number of classes
   scaling = 10                                              !!>> HC 28-9-2021 Scaling of the avangle vector to fit into the classes

   if(allocated(cooccurrences))deallocate(cooccurrences)     !!>> HC 28-9-2021
   allocate( cooccurrences( nclasseshc, nclasseshc ) )       !!>> HC 28-9-2021
   if(allocated(cooprobs))deallocate(cooprobs)               !!>> HC 28-9-2021
   allocate( cooprobs( nclasseshc, nclasseshc ) )            !!>> HC 28-9-2021
   if(allocated(classprod))deallocate(classprod)             !!>> HC 28-9-2021
   allocate( classprod( nclasseshc, nclasseshc ) )           !!>> HC 28-9-2021
   rnapical = napical !!>> HC 28-9-2021 In case we have to operate with napical as a real number 

   entropy1=0; entropy2=0                                    !!>> HC 28-9-2021

   do tiponodhc=1, 2                                         !!>> HC 28-9-2021
      if (tiponodhc==1) then                                 !!>> HC 28-9-2021
      !!>> HC 28-9-2021 This reads the file were the initial 4th neightbour associations are written down
	 open (777,file="fourth_original_neighs.dat")             !!>> HC 28-9-2021 
	      do ihc=1, finitnn                                   !!>> HC 28-9-2021
	         read(777,*) foldnodes(ihc), foldnodesneight(ihc) !!>> HC 28-9-2021
	      enddo                                               !!>> HC 28-9-2021
	 close(777)                                               !!>> HC 28-9-2021
      else                                                        !!>> HC 28-9-2021
        !!>> HC 28-9-2021 This reads the file were the initial 4th neightbour associations are written down
	 open ( 777, file = "fourth_original_neighs2.dat" )       !!>> HC 28-9-2021
	      do ihc = 1, finitnn                                 !!>> HC 28-9-2021
		 read(777,*) foldnodes(ihc), foldnodesneight(ihc) !!>> HC 28-9-2021
	      enddo                                               !!>> HC 28-9-2021
	 close(777)                                               !!>> HC 28-9-2021
      endif                                                       !!>> HC 28-9-2021

 	!!>> HC 28-9-2021 this sets the dimension of the vector of average angles to the number of apical nodes

	!initializing 
      avangles=0; classangles=0; classtem=0; probclass=0           !!>> HC 28-9-2021
      cooccurrences=0                                              !!>> HC 28-9-2021
      oldorder=0; cooccurrences=0; cooprobs=0; tempclass=0         !!>> HC 28-9-2021
      totalneightshc=0; medianglee=0;  sumangle=0 ; numn=0         !!>> HC 28-9-2021
      pij=0; u=0; v=0; order=0; order2=0 ; angles1=0; angles2=0    !!>> HC 28-9-2021
	             

	!!!!!>> HC 28-9-2021 NEIGHBOURS !!!
	!!!!!>> HC 28-9-2021 This is the neightbour finding/angle calculating loop!!!
 
      do ihc=1, finitnn                                                   !!>> HC 28-9-2021 Look in the list of initial 2nd neightbourd associations
         u(1)=(node(node(foldnodes(ihc))%altre)%x-node(foldnodes(ihc))%x) !!>> HC 28-9-2021 Apical-basal vector i cell  
         u(2)=(node(node(foldnodes(ihc))%altre)%y-node(foldnodes(ihc))%y) !!>> HC 28-9-2021
         u(3)=(node(node(foldnodes(ihc))%altre)%z-node(foldnodes(ihc))%z) !!>> HC 28-9-2021
         v(1)=(node(foldnodesneight(ihc))%x-node(foldnodes(ihc))%x)       !!>> HC 28-9-2021 Vector between apical i cell node and the apical j cell node	
         v(2)=(node(foldnodesneight(ihc))%y-node(foldnodes(ihc))%y)       !!>> HC 28-9-2021
         v(3)=(node(foldnodesneight(ihc))%z-node(foldnodes(ihc))%z)       !!>> HC 28-9-2021

         upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)                                !!>> HC 28-9-2021 vectorial product u * v
         modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)                         !!>> HC 28-9-2021 modulus u
         modv=sqrt((v(1))**2+(v(2))**2+(v(3))**2)                         !!>> HC 28-9-2021 modulus v
         angle=upv/(modu*modv)                                            !!>> HC 28-9-2021 getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
		
         if (ihc==finitnn)then                     !!>> HC 28-9-2021
            sumangle=sumangle+angle                !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i (last node of the class must be taken into account)
            numn=numn +1                           !!>> HC 28-9-2021 Counter of the real number of neights (last node of the class must be taken into account)		
            order=order +1                         !!>> HC 28-9-2021 This is actually the order of the nodes, fom 1st to 271st
            medianglee=sumangle/numn               !!>> HC 28-9-2021 This obtains the average dot product of the cell i with its neightbours
            avangles(order)=medianglee             !!>> HC 28-9-2021 This stores this average in the avangles vector
            numn=0                                 !!>> HC 28-9-2021 This resets the number of neights
            sumangle=0                             !!>> HC 28-9-2021 This resets the summation of dot products
            oldorder(order)=foldnodes(ihc)         !!>> HC 28-9-2021 Saves the old node number to use it later when calculating the cooccurrence                          
         else
            if (foldnodes(ihc)==foldnodes(ihc+1)) then    !!>> HC 28-9-2021 If the next node of the list is still the same
               sumangle=sumangle+angle                    !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i 
               numn=numn+1                                !!>> HC 28-9-2021 Counter of the real number of neights 
            else                                          !!>> HC 28-9-2021 If the next node is different
               sumangle=sumangle+angle                    !!>> HC 28-9-2021 sum of angles of the neighbouts in a cell i (last node of the class must be taken into account)
               numn=numn+1                                !!>> HC 28-9-2021 Counter of the real number of neights (last node of the class must be taken into account)		
               order=order+1                              !!>> HC 28-9-2021 This is actually the order of the nodes, fom 1st to 271st
               medianglee=sumangle/numn                   !!>> HC 28-9-2021 This obtains the average dot product of the cell i with its neightbours
               avangles(order)=medianglee                 !!>> HC 28-9-2021 This stores this average in the avangles vector
               numn=0                                     !!>> HC 28-9-2021 This resets the number of neights
               sumangle=0                                 !!>> HC 28-9-2021 This resets the summation of dot products
               oldorder(order)=foldnodes(ihc)             !!>> HC 28-9-2021 Saves the old node number to use it later when calculating the cooccurrence
            endif                                         !!>> HC 28-9-2021
         endif                                            !!>> HC 28-9-2021
      enddo                                               !!>> HC 28-9-2021


	!!!!>> HC 28-9-2021 ANGLE CLASSES!!
	!!!>> HC 28-9-2021 This divides the averange angles in classes each 0.1
      do ihc=1,napical                                      !!>> HC 28-9-2021
         do jhc=lowlimit, highlimit, steps                  !!>> HC 28-9-2021 The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
            if (avangles(ihc)*scaling .ge. jhc ) then       !!>> HC 28-9-2021
               if (avangles(ihc)*scaling < jhc +1 ) then    !!>> HC 28-9-2021
                  classangles(ihc)=jhc                      !!>> HC 28-9-2021
               endif                                        !!>> HC 28-9-2021
            endif                                           !!>> HC 28-9-2021
         enddo                                              !!>> HC 28-9-2021
      enddo                                                 !!>> HC 28-9-2021


	!!! COOCURRENCE !!!!!>> HC 28-9-2021
	!!! This is the neightbour finding/ cooccurrence calculating loop!!!!!>> HC 28-9-2021

      order=0                                               !!>> HC 28-9-2021
      !!!>> HC 28-9-2021 This loop replaces the number of node by its curvature class
      !!!>> HC 28-9-2021 In the original lists of nodes and neigbours
      do ihc=1, finitnn                                     !!>> HC 28-9-2021
         do jhc=1, napical                                  !!>> HC 28-9-2021
	     if (foldnodes(ihc)==oldorder(jhc))then         !!>> HC 28-9-2021
                order=order+1                               !!>> HC 28-9-2021
                angles1(order)=classangles(jhc)             !!>> HC 28-9-2021
             endif                                          !!>> HC 28-9-2021
             if (foldnodesneight(ihc)==oldorder(jhc))then   !!>> HC 28-9-2021
                order2=order2+1                             !!>> HC 28-9-2021
                angles2(order2)=classangles(jhc)            !!>> HC 28-9-2021
             endif                                          !!>> HC 28-9-2021
         enddo                                              !!>> HC 28-9-2021
      enddo                                                 !!>> HC 28-9-2021

      ! We have to add 11 to transform the curvature class into the   !!>> HC 28-9-2021
      ! number of row/column in the matrix cooccurrence matrix        !!>> HC 28-9-2021
      angles1=angles1+11                                              !!>> HC 28-9-2021
      angles2=angles2+11                                              !!>> HC 28-9-2021

      !This loop assign each combination of curvatures to a cell in the coocurrence matrix   !!>> HC 28-9-2021
      do ihc=1, finitnn                                                                      !!>> HC 28-9-2021
         cooccurrences(angles1(ihc),angles2(ihc))=cooccurrences(angles1(ihc),angles2(ihc))+1 !!>> HC 28-9-2021
      enddo                                                                                  !!>> HC 28-9-2021

      !This loop gets the total number of node-neights in the embryo !!>> HC 28-9-2021
      do ihc=1, nclasseshc                                           !!>> HC 28-9-2021
         do jhc=1, nclasseshc                                        !!>> HC 28-9-2021
            totalneightshc=totalneightshc+cooccurrences(ihc,jhc)     !!>> HC 28-9-2021
         enddo                                                       !!>> HC 28-9-2021
      enddo                                                          !!>> HC 28-9-2021

      !The probabilities are the number of neights with a given combination of curvatures between the total num of neights !!>> HC 28-9-2021
      cooprobs=cooccurrences/totalneightshc                                                                                !!>> HC 28-9-2021


	
      if (tiponodhc==1)then                                                       !!>> HC 28-9-2021 CALCULATING JOINT ENTROPY
         !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)    !!>> HC 28-9-2021
         do ihc=1, nclasseshc                                                     !!>> HC 28-9-2021
            do jhc=1, nclasseshc                                                  !!>> HC 28-9-2021
	       if (jhc>ihc)cycle                                                  !!>> HC 11-11-2021
	       if (cooprobs(ihc,jhc)==0.0d0) cycle                                    !!>> HC 11-11-2021
	       if (jhc==ihc)then                                                  !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)                                           !!>> HC 11-11-2021
	       else                                                               !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)+cooprobs(jhc,ihc)                         !!>> HC 11-11-2021
               endif                                                              !!>> HC 28-9-2021
               entropy1=entropy1+(pij*log(pij)/log(2d0))                          !!>> HC 28-9-2021
            enddo                                                                 !!>> HC 28-9-2021
         enddo                                                                    !!>> HC 28-9-2021
         entropy1=-entropy1                                                       !!>> HC 28-9-2021
      else                                                                        !!>> HC 28-9-2021
         !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)    !!>> HC 28-9-2021
         do ihc=1, nclasseshc                                                     !!>> HC 28-9-2021
            do jhc=1, nclasseshc                                                  !!>> HC 28-9-2021
	       if (jhc>ihc)cycle                                                  !!>> HC 11-11-2021
	       if (cooprobs(ihc,jhc)==0.0d0) cycle                                !!>> HC 11-11-2021
	       if (jhc==ihc)then                                                  !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)                                           !!>> HC 11-11-2021
	       else                                                               !!>> HC 11-11-2021
	          pij=cooprobs(ihc,jhc)+cooprobs(jhc,ihc)                         !!>> HC 11-11-2021
               endif                                                              !!>> HC 28-9-2021                                             !!>> HC 28-9-2021
               entropy2=entropy2+(pij*log(pij)/log(2d0))                          !!>> HC 28-9-2021
            enddo                                                                 !!>> HC 28-9-2021
         enddo                                                                    !!>> HC 28-9-2021
         entropy2=-entropy2                                                       !!>> HC 28-9-2021
      endif                                                                       !!>> HC 28-9-2021
  enddo                                                                           !!>> HC 28-9-2021

    entropy=(entropy2+entropy1)/2          !!>> HC 20-10-2021
    fitelli=(entropy)**2                   !!>> HC 21-12-2021 
    
  end subroutine                           !!>> HC 28-9-2021
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
 
  subroutine flux_resistance_decrease(fitelli)                         !!>> HC 11-11-2021 This subroutine was written by HUgo Cano and it calculates the flux resistance of a morphology
    implicit none                                                      !!>> HC 11-11-2021 calculated as the average dot product between an arbritary flux vector and the spring vectors.
    integer ::  ihc, jhc, khc, nnodes, naltech                         !!>> HC 11-11-2021
    real*8, dimension(1:3) :: flux, u                                  !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu, fitelli        !!>> HC 11-11-2021
 
    call neighbor_build                                                !!>> HC 11-11-2021

    flux(1)=0.0d0; flux(2)=0.0d0; flux(3)=1.0d0                        !!>> HC 11-11-2021 This is the FLUX VECTOR in the environment (arbitrary)
    nnodes=0; sumd=0.0d0; average_flux_resistance=0.0d0; modu=0.0d0    !!>> HC 11-11-2021 it has to be a unit vector

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.1)cycle                                   !!>> HC 11-11-2021
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021

    average_flux_resistance=(1-(sumd/nnodes))                          !!>> HC 11-11-2021 this is the average fux resistance
    fitelli=average_flux_resistance
      
  end subroutine                                                       !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
 
  subroutine flux_resistance_increase(fitelli)                         !!>> HC 11-11-2021 This subroutine was written by HUgo Cano and it calculates the flux resistance of a morphology
    implicit none                                                      !!>> HC 11-11-2021 calculated as the average dot product between an arbritary flux vector and the spring vectors.
    integer ::  ihc, jhc, khc, nnodes, naltech                         !!>> HC 11-11-2021
    real*8, dimension(1:3) :: flux, u                                  !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu, fitelli        !!>> HC 11-11-2021
 
    call neighbor_build                                                !!>> HC 11-11-2021

    flux(1)=0.0d0; flux(2)=0.0d0; flux(3)=1.0d0                        !!>> HC 11-11-2021 This is the FLUX VECTOR in the environment (arbitrary)
    nnodes=0; sumd=0.0d0; average_flux_resistance=0.0d0; modu=0.0d0    !!>> HC 11-11-2021 it has to be a unit vector

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.1)cycle                                   !!>> HC 11-11-2021
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021

    average_flux_resistance=(sumd/nnodes)                              !!>> HC 11-11-2021 this is the average fux resistance
    fitelli=average_flux_resistance                                    !!>> HC 11-11-2021
  
  end subroutine                                                       !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021

  subroutine flux_half_increase_half_decrease(fitelli)                        !!>> HC 12-5-2022
  Implicit none                                                               !!>> HC 12-5-2022
  integer ::  ich, jch, kch, lch, epind                                       !!>> HC 12-5-2022
  real*8, allocatable, dimension(:,:) ::compind,  tarind, orinodes            !!>> HC 12-5-2022
  real*8 :: EMDval, centroidx, centroidy, centroidz, req1, storenorm, tarreq  !!>> HC 12-5-2022
  integer ::  ihc, jhc, khc, nnodes, naltech                                  !!>> HC 11-11-2021
  real*8, dimension(1:3) :: flux, u                                           !!>> HC 11-11-2021
  real*8 :: upv, sumd, average_flux_resistance, modu, fitelli,average_flux_resistance2        !!>> HC 11-11-2021
  integer, allocatable, dimension(:) :: nodid
  
  if (allocated(orinodes)) deallocate(orinodes)                 !!>> HC 24-9-2021
  allocate(orinodes(nd,3))                                      !!>> HC 24-9-2021
  
  if (allocated(nodid)) deallocate(nodid)                       !!>> HC 24-9-2021
  allocate(nodid(nd))                                           !!>> HC 24-9-2021
  nodid=0
  
  if (allocated(compind)) deallocate(compind)                   !!>> HC 24-9-2021
  allocate(compind(nd,3))                                       !!>> HC 24-9-2021

  if (allocated(tarind)) deallocate(tarind)                     !!>> HC 24-9-2021
  allocate(tarind(nd,3))                                        !!>> HC 24-9-2021


  compind=0.0d0; tarind=0.0d0                                   !!>> HC 12-5-2022
  lch=0; kch=0; epind=0; EMDval=0.0d0; orinodes=0.0d0           !!>> HC 12-5-2022

  do ich=1,nd                                                   !!>> HC 12-5-2022
     if(node(ich)%tipus.ne.2)cycle                              !!>> HC 12-5-2022
     epind=epind+1                                              !!>> HC 12-5-2022
     orinodes(ich,1)=node(ich)%x                                !!>> HC 24-9-2021
     orinodes(ich,2)=node(ich)%y                                !!>> HC 24-9-2021
     orinodes(ich,3)=node(ich)%z                                !!>> HC 24-9-2021
  enddo
       
   centroidx=sum(orinodes(:,1))/epind                           !!>> HC 24-9-2021 !centroid
   centroidy=sum(orinodes(:,2))/epind                           !!>> HC 24-9-2021
   centroidz=sum(orinodes(:,3))/epind                           !!>> HC 24-9-2021
   epind=0 
   do ich=1,nd                                                  !!>> HC 24-9-2021
      if(node(ich)%tipus.ne.2)cycle
      orinodes(ich,1)=orinodes(ich,1)-centroidx                        !!>> HC 24-9-2021 !centering 
      orinodes(ich,2)=orinodes(ich,2)-centroidy                        !!>> HC 24-9-2021
      orinodes(ich,3)=orinodes(ich,3)-centroidz                        !!>> HC 24-9-2021
   end do                                                              !!>> HC 24-9-2021    

  do ich=1,nd                                                   !!>> HC 12-5-2022
     if(node(ich)%tipus.ne.2)cycle                              !!>> HC 12-5-2022
     if (orinodes(ich,3)>centroidz)then                         !!>> HC 12-5-2022
        nodid(ich)=1                                            !!>> HC 12-5-2022
     else                                                       !!>> HC 12-5-2022
        nodid(ich)=2                                            !!>> HC 12-5-2022
     endif                                                      !!>> HC 12-5-2022
  enddo                                                         !!>> HC 12-5-2022

    flux(1)=0.0d0; flux(2)=0.0d0; flux(3)=1.0d0                        !!>> HC 11-11-2021 This is the FLUX VECTOR in the environment (arbitrary)
    nnodes=0; sumd=0.0d0; average_flux_resistance=0.0d0; modu=0.0d0    !!>> HC 11-11-2021 it has to be a unit vector

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.2)cycle                                   !!>> HC 11-11-2021
       if(nodid(ihc).ne.1)cycle                                        !!>> HC 12-5-2022
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021

    average_flux_resistance=(1-(sumd/nnodes))                          !!>> HC 11-11-2021 this is the average fux resistance
    
    nnodes=0; sumd=0.0d0; modu=0.0d0;average_flux_resistance2=0.0d0    !!>> HC 11-11-2021

    do ihc=1, nd                                                       !!>> HC 11-11-2021
       if(node(ihc)%tipus.ne.2)cycle                                   !!>> HC 11-11-2021
       if(nodid(ihc).ne.2)cycle                                        !!>> HC 12-5-2022
       nnodes=nnodes+1                                                 !!>> HC 11-11-2021
       naltech=node(ihc)%altre                                         !!>> HC 11-11-2021
       u(1)=node(ihc)%x-node(naltech)%x                                !!>> HC 11-11-2021
       u(2)=node(ihc)%y-node(naltech)%y                                !!>> HC 11-11-2021
       u(3)=node(ihc)%z-node(naltech)%z                                !!>> HC 11-11-2021
       modu=sqrt(u(1)**2+u(2)**2+u(3)**2)                              !!>> HC 11-11-2021
       u(1)=u(1)/modu                                                  !!>> HC 11-11-2021 u is the unit spring vector of all the epithelial cells
       u(2)=u(2)/modu                                                  !!>> HC 11-11-2021
       u(3)=u(3)/modu                                                  !!>> HC 11-11-2021
       upv=flux(1)*u(1)+flux(2)*u(2)+flux(3)*u(3)                      !!>> HC 11-11-2021 this is the dot producti between the spring vector and the flux vector
       sumd=sumd+abs(upv)                                              !!>> HC 11-11-2021
    enddo                                                              !!>> HC 11-11-2021
    average_flux_resistance2=(sumd/nnodes)
         open(206,file="complexity_vals.dat")                                   !!>> HC 11-11-2021
             write(206,*) average_flux_resistance, average_flux_resistance2     !!>> HC 11-11-2021
         close(206)                                                             !!>> HC 11-11-2021
    fitelli=average_flux_resistance+average_flux_resistance2                    !!>> HC 11-11-2021    
   
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine local_joint_entropy_second(fitelli)                                                   !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                           !!>> HC 11-11-2021 and it calculates the local joint entropy

    real*8 :: angle, upv, modu, modv, sumangle, medianglee, entropy, entropy1, disthc
    real*8 :: entropy2, pij, pij_pipj, entropy3,  rnapical, totalneightshc, fitelli
    integer ::  napical, napicali, ihc, jhc, khc, lhc, mhc, order, order2, node2_1, numn, ndouble, neichi
    integer:: tempclass, lowlimit, highlimit, steps, tiponodhc, scaling, nclasseshc, neich, classtemsum
    integer, allocatable, dimension(:) :: classangles, oldorder, dummycol
    integer, allocatable, dimension(:,:) :: double_neighs, trans_double_neighs
    integer, allocatable, dimension(:) :: double_nneighs, trans_double_nneighs
    real*8, dimension(1:3) :: u, v
    real*8, allocatable, dimension(:) :: avangles, probclass, xcoord, ycoord, zcoord
    real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs, classprod


    disthc=5     ! Distance of the neightbours
    tiponodhc=2  ! Type of node we are aiming 2= apical 1 = basal
    lowlimit=-10 ! Lower limit of the interval in classes
    highlimit=9  ! Higher limit of the interval in classes
    steps=1      ! width of the interval in classes
    nclasseshc=20!Number of classes
    scaling=10   ! Scaling of the avangle vector to fit into the classes
    ndouble=maxval(nneigh(1:nd))*2

    ! This calculates the number of apical cells
    napical=0
    do ihc=1,nd                            !This loop calculates the number of apical nodes that are going to be inspected
       if (node(ihc)%tipus==tiponodhc)then !Select only the apical nodes
          napical=napical+1                !Counter of the number of apical nodes inspected
       endif
    enddo

     ! this sets the dimension of the vector of average angles to the number of apical nodes
 
    if(allocated(avangles))deallocate(avangles)
    allocate(avangles(nd)) 

    if(allocated(classangles))deallocate(classangles)
    allocate(classangles(nd)) 

    if(allocated(probclass))deallocate(probclass)
    allocate(probclass(nclasseshc)) 

    if(allocated(cooccurrences))deallocate(cooccurrences)
    allocate(cooccurrences(nclasseshc,nclasseshc)) 
  
    if(allocated(cooprobs))deallocate(cooprobs)
    allocate(cooprobs(nclasseshc,nclasseshc)) 

    if(allocated(classprod))deallocate(classprod)
    allocate(classprod(nclasseshc,nclasseshc)) 

    if(allocated(double_neighs))deallocate(double_neighs)
    allocate(double_neighs(nd,ndouble)) 

    if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
    allocate(trans_double_neighs(nd,ndouble)) 

    if(allocated(double_nneighs))deallocate(double_nneighs)
    allocate(double_nneighs(nd)) 
  
    if(allocated(trans_double_nneighs))deallocate(trans_double_nneighs)
    allocate(trans_double_nneighs(napical)) 
  
    entropy=0.0d0; entropy1=0.0d0; entropy2=0.0d0; entropy3=0.0d0

    do tiponodhc=1,2
       avangles=666.0d0; classangles=666; probclass=0; napicali=0; rnapical=napical
       cooccurrences=0; cooprobs=0.0d0; tempclass=0
       totalneightshc=0.0d0; double_nneighs=0; double_neighs=0
       pij_pipj=0; pij=0; order=0  

       !!Finding second order neighbors (can be optimized...)
       order2=0; order=0
       do ihc=1,nd
          if(node(ihc)%tipus.ne.tiponodhc)cycle
          order=0
          do jhc=1, nneigh(ihc)
             neich=neigh(ihc,jhc)
             if(neich==ihc)cycle
             if(node(neich)%tipus.ne.tiponodhc)cycle
             do khc=1,nneigh(neich)
                neichi=neigh(neich,khc)
                if(node(neichi)%tipus.ne.tiponodhc)cycle
                if(neichi==neich)cycle
                if(neichi==ihc)cycle
                if(any(neigh(ihc,1:nneigh(ihc))==neichi))cycle
                if (order>0)then
                   if(any(double_neighs(ihc,1:double_nneighs(ihc))==neichi))cycle
                endif
                order=order+1
                if (order>ndouble)then
                   if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
                   allocate(trans_double_neighs(nd,ndouble))
                   trans_double_neighs=0
                   trans_double_neighs(1:nd,1:ndouble)=double_neighs(1:nd,1:ndouble)  
                   ndouble=ndouble+1    
                   if(allocated(double_neighs))deallocate(double_neighs)
                   allocate(double_neighs(nd,ndouble))
                   double_neighs=0
                   double_neighs(1:nd,1:ndouble-1)=trans_double_neighs(1:nd,1:ndouble-1)    
                endif
                double_neighs(ihc,order)=neichi
                double_nneighs(ihc)=double_nneighs(ihc)+1
             enddo
          enddo
       enddo

   !!! This is th neightbour finding/angle calculating loop!!!

       do ihc=1,nd                                       ! Inspecting all the nodes i of the embryo
          sumangle=0 ; numn=0                            ! initializing      
          if (node(ihc)%tipus.ne.tiponodhc )cycle        ! We are just looking for apical cells  
          do jhc=1,double_nneighs(ihc)                           ! Comaring them with its neighbors
             neich=double_neighs(ihc,jhc)
             if(ihc==neich)cycle                         ! Not themselves
             if(node(neich)%tipus.ne.tiponodhc)cycle     ! We are just looking for apical cells 
             v(1)=(node(neich)%x-node(ihc)%x)            ! Vector between apical i cell node and the apical j cell node
             v(2)=(node(neich)%y-node(ihc)%y)
             v(3)=(node(neich)%z-node(ihc)%z)	
             modv=sqrt( (v(1))**2+(v(2))**2+(v(3))**2)   ! modulus v
             u(1)=(node(node(ihc)%altre)%x-node(ihc)%x)  ! Apical-basal vector i cell  
             u(2)=(node(node(ihc)%altre)%y-node(ihc)%y)
             u(3)=(node(node(ihc)%altre)%z-node(ihc)%z)
             modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)    ! modulus u
             upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)           ! vectorial product u * v
             angle=upv/(modu*modv)                       ! getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
             sumangle=sumangle+angle                     ! sum of angles of the neighbouts in a cell i 
             numn=numn+1                                 ! Counter of the real number of neights (just apical nodes) of the cell i
          enddo
          if (numn.ne.0) then               ! Just in case a cell has no neights
             napicali=napicali+1            ! This is a counter of the number of cells i with neights
                                            ! it will be used after to get the probability of classes
             medianglee=sumangle/real(numn) !This obtains the average angle of the cell i with its neightbours
             avangles(ihc)=medianglee
          endif  
       enddo

       order=0

       ! This divides the average angles in classes each 0.1
       do ihc=1,nd
          if (avangles(ihc)==666.0d0)cycle
          do jhc=lowlimit,highlimit, steps   ! The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
             if (avangles(ihc)*scaling.ge.jhc)then 
                if (avangles(ihc)*scaling<jhc+1)then
	            classangles(ihc)=jhc
                endif
             endif
          enddo
       enddo

       !!! COOCURRENCE !!!

       !!! This is the neightbour finding/ cooccurrence calculating loop!!!
       v=0.0d0
       order=0; order2=0  

       do khc=lowlimit, highlimit, steps                         !! Go over all the classes
          order=order+1                                          !! This is the class i counter 
          do lhc=lowlimit, highlimit, steps                      !! Go over all the combinations of clases
             order2=lhc+11                                        !! This is the class j counter          
             if (khc>lhc)cycle                                      !! but we consider i,j the same category as j,i
             do ihc=1,nd                                         !! go over all nodes
                if (classangles(ihc).ne.khc)cycle                !! find those whose class is equal to khc
                do jhc=1, double_nneighs(ihc)                            !! go over all their neighbors
                   neich=double_neighs(ihc,jhc)                          !! this is the neighbor
                   if (neich>ihc)cycle                           !! We do not want to consider neighbor interaction twice (although the result should be the same)
                   if (classangles(neich).ne.lhc)cycle             !! check whether the angle class of neich is lhc, if yes we've got a match!!
                   cooccurrences(order,order2)=cooccurrences(order,order2)+1.0d0  !! write down the match in the coocurrences matrix
                   totalneightshc=totalneightshc+1.0d0           !! counter of the total number of coocurrences of all angle classes
                enddo
             enddo
          enddo
       enddo



       !The probabilities are the number of neights with a given combination of curvatures between the total num of neights
       cooprobs=cooccurrences/totalneightshc

       if (tiponodhc==1)then
          !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)
          do ihc=1, nclasseshc
             do jhc=1, nclasseshc
                if (cooprobs(ihc,jhc)==0.0d0)cycle
                pij=cooprobs(ihc,jhc)
                entropy1=entropy1+(pij*log(pij)/log(2.0000000000))
	     enddo
          enddo
          pij=0
          entropy1=-entropy1
       else
          !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)
          do ihc=1, nclasseshc
             do jhc=1, nclasseshc
                if (cooprobs(ihc,jhc)==0.0d0)cycle
                pij=cooprobs(ihc,jhc)
                entropy2=entropy2+(pij*log(pij)/log(2.0000000000))
             enddo
          enddo
          pij=0
          entropy2=-entropy2 
       endif
    enddo

    entropy=(entropy2+entropy1)/2.0d0                                 !!>> HC 11-11-2021
    fitelli=entropy                                                   !!>> HC 11-11-2021

  end subroutine                                                        !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  
  subroutine local_joint_entropy_fourth(fitelli)                                                   !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                           !!>> HC 11-11-2021 and it calculates the local joint entropy
    real*8 :: angle, upv, modu, modv, sumangle, medianglee, entropy, entropy1, disthc
    real*8 :: entropy2, pij, pij_pipj, entropy3,  rnapical, totalneightshc, fitelli
    integer ::  napical, napicali, ihc, jhc, khc, lhc, mhc, order, order2, node2_1, numn, neichi
    integer :: ndouble, ntriple, nfourth
    integer:: tempclass, lowlimit, highlimit, steps, tiponodhc, scaling, nclasseshc, neich, classtemsum
    integer, allocatable, dimension(:) :: classangles, oldorder, dummycol
    integer, allocatable, dimension(:,:) :: double_neighs, trans_double_neighs
    integer, allocatable, dimension(:) :: double_nneighs, trans_double_nneighs
    integer, allocatable, dimension(:,:) :: triple_neighs, trans_triple_neighs
    integer, allocatable, dimension(:) :: triple_nneighs
    integer, allocatable, dimension(:,:) :: fourth_neighs, trans_fourth_neighs
    integer, allocatable, dimension(:) :: fourth_nneighs
    real*8, dimension(1:3) :: u, v
    real*8, allocatable, dimension(:) :: avangles, probclass, xcoord, ycoord, zcoord
    real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs, classprod


    disthc=5     ! Distance of the neightbours
    tiponodhc=2  ! Type of node we are aiming 2= apical 1 = basal
    lowlimit=-10 ! Lower limit of the interval in classes
    highlimit=9  ! Higher limit of the interval in classes
    steps=1      ! width of the interval in classes
    nclasseshc=20!Number of classes
    scaling=10   ! Scaling of the avangle vector to fit into the classes
    ndouble=maxval(nneigh(1:nd))*2
    ntriple=maxval(nneigh(1:nd))*2
    nfourth=maxval(nneigh(1:nd))*2

    ! This calculates the number of apical cells
    napical=0
    do ihc=1,nd                            !This loop calculates the number of apical nodes that are going to be inspected
       if (node(ihc)%tipus==tiponodhc)then !Select only the apical nodes
          napical=napical+1                !Counter of the number of apical nodes inspected
       endif
    enddo

    if(allocated(avangles))deallocate(avangles)
    allocate(avangles(nd)) 

    if(allocated(classangles))deallocate(classangles)
    allocate(classangles(nd)) 

    if(allocated(probclass))deallocate(probclass)
    allocate(probclass(nclasseshc)) 

    if(allocated(cooccurrences))deallocate(cooccurrences)
    allocate(cooccurrences(nclasseshc,nclasseshc)) 
  
    if(allocated(cooprobs))deallocate(cooprobs)
    allocate(cooprobs(nclasseshc,nclasseshc)) 

    if(allocated(classprod))deallocate(classprod)
    allocate(classprod(nclasseshc,nclasseshc)) 

    if(allocated(double_neighs))deallocate(double_neighs)
    allocate(double_neighs(nd,ndouble)) 

    if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
    allocate(trans_double_neighs(nd,ndouble)) 

    if(allocated(double_nneighs))deallocate(double_nneighs)
    allocate(double_nneighs(nd)) 
  
    if(allocated(triple_neighs))deallocate(triple_neighs)
    allocate(triple_neighs(nd,ntriple)) 

    if(allocated(trans_triple_neighs))deallocate(trans_triple_neighs)
    allocate(trans_triple_neighs(nd,ntriple)) 

    if(allocated(triple_nneighs))deallocate(triple_nneighs)
    allocate(triple_nneighs(nd)) 
  
    if(allocated(fourth_neighs))deallocate(fourth_neighs)
    allocate(fourth_neighs(nd,nfourth)) 

    if(allocated(trans_fourth_neighs))deallocate(trans_fourth_neighs)
    allocate(trans_fourth_neighs(nd,nfourth)) 

    if(allocated(fourth_nneighs))deallocate(fourth_nneighs)
    allocate(fourth_nneighs(nd))
  
  
    entropy=0.0d0; entropy1=0.0d0; entropy2=0.0d0; entropy3=0.0d0

    do tiponodhc=1,2
       avangles=666.0d0; classangles=666; probclass=0; napicali=0; rnapical=napical
       cooccurrences=0; cooprobs=0.0d0; tempclass=0
       totalneightshc=0.0d0; double_nneighs=0; double_neighs=0
       triple_nneighs=0; triple_neighs=0; fourth_nneighs=0; fourth_neighs=0
       pij_pipj=0; pij=0; order=0

       !!Finding second order neighbors (can be optimized...)
       order2=0; order=0
       do ihc=1,nd
          if(node(ihc)%tipus.ne.tiponodhc)cycle
          order=0
          do jhc=1, nneigh(ihc)
             neich=neigh(ihc,jhc)
             if(neich==ihc)cycle
             if(node(neich)%tipus.ne.tiponodhc)cycle
             do khc=1,nneigh(neich)
                neichi=neigh(neich,khc)
                if(node(neichi)%tipus.ne.tiponodhc)cycle
                if(neichi==neich)cycle
                if(neichi==ihc)cycle
                if(any(neigh(ihc,1:nneigh(ihc))==neichi))cycle
                if (order>0)then
                   if(any(double_neighs(ihc,1:double_nneighs(ihc))==neichi))cycle
                endif
                order=order+1
                if (order>ndouble)then
                   if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
                   allocate(trans_double_neighs(nd,ndouble))
                   trans_double_neighs=0
                   trans_double_neighs(1:nd,1:ndouble)=double_neighs(1:nd,1:ndouble)  
                   ndouble=ndouble+1    
                   if(allocated(double_neighs))deallocate(double_neighs)
                   allocate(double_neighs(nd,ndouble))
                   double_neighs=0
                   double_neighs(1:nd,1:ndouble-1)=trans_double_neighs(1:nd,1:ndouble-1)    
                endif
                double_neighs(ihc,order)=neichi
                double_nneighs(ihc)=double_nneighs(ihc)+1
             enddo
          enddo
       enddo
 
     order=0  
       do ihc=1,nd
          order=0
          do jhc=1, double_nneighs(ihc)
             neich=double_neighs(ihc,jhc)
             do lhc=1,nneigh(neich)
                neichi=neigh(neich,lhc)
                if(node(neichi)%tipus.ne.tiponodhc)cycle
                if(neichi==ihc)cycle
                if(neichi==neich)cycle
                if(any(neigh(ihc,1:nneigh(ihc))==neichi))cycle
                if(any(double_neighs(ihc,1:double_nneighs(ihc))==neichi))cycle
                if (order>0)then
                   if(any(triple_neighs(ihc,1:triple_nneighs(ihc))==neichi))cycle
                endif
                order=order+1
                if (order>ntriple)then
                   if(allocated(trans_triple_neighs))deallocate(trans_triple_neighs)
                   allocate(trans_triple_neighs(nd,ntriple))
                   trans_triple_neighs=0
                   trans_triple_neighs(1:nd,1:ntriple)=triple_neighs(1:nd,1:ntriple)  
                   ntriple=ntriple+1    
                   if(allocated(triple_neighs))deallocate(triple_neighs)
                   allocate(triple_neighs(nd,ntriple))
                   triple_neighs=0
                   triple_neighs(1:nd,1:ntriple-1)=trans_triple_neighs(1:nd,1:ntriple-1)    
                endif
                triple_neighs(ihc,order)=neichi
                triple_nneighs(ihc)=triple_nneighs(ihc)+1
             enddo
          enddo
       enddo
   
       order=0 
       do ihc=1,nd
          order=0
          do jhc=1, triple_nneighs(ihc)
             neich=triple_neighs(ihc,jhc)
             do lhc=1,nneigh(neich)
                neichi=neigh(neich,lhc)
                if(node(neichi)%tipus.ne.tiponodhc)cycle
                if(neichi==ihc)cycle
                if(neichi==neich)cycle
                if(any(neigh(ihc,1:nneigh(ihc))==neichi))cycle
                if(any(double_neighs(ihc,1:double_nneighs(ihc))==neichi))cycle
                if(any(triple_neighs(ihc,1:triple_nneighs(ihc))==neichi))cycle
                if (order>0)then
                   if(any(fourth_neighs(ihc,1:fourth_nneighs(ihc))==neichi))cycle
                endif
                order=order+1
                if (order>nfourth)then
                   if(allocated(trans_fourth_neighs))deallocate(trans_fourth_neighs)
                   allocate(trans_fourth_neighs(nd,nfourth))
                   trans_fourth_neighs=0
                   trans_fourth_neighs(1:nd,1:nfourth)=fourth_neighs(1:nd,1:nfourth)  
                   nfourth=nfourth+1    
                   if(allocated(fourth_neighs))deallocate(fourth_neighs)
                   allocate(fourth_neighs(nd,nfourth))
                   fourth_neighs=0
                   fourth_neighs(1:nd,1:nfourth-1)=trans_fourth_neighs(1:nd,1:nfourth-1)    
                endif
                fourth_neighs(ihc,order)=neichi
                fourth_nneighs(ihc)=fourth_nneighs(ihc)+1
             enddo
          enddo
       enddo
  
       !!! This is th neightbour finding/angle calculating loop!!!

       do ihc=1,nd                                       ! Inspecting all the nodes i of the embryo
          sumangle=0 ; numn=0                            ! initializing      
          if (node(ihc)%tipus.ne.tiponodhc )cycle        ! We are just looking for apical cells  
          do jhc=1,fourth_nneighs(ihc)                           ! Comaring them with its neighbors
             neich=fourth_neighs(ihc,jhc)
             if(ihc==neich)cycle                         ! Not themselves
             if(node(neich)%tipus.ne.tiponodhc)cycle     ! We are just looking for apical cells 
             v(1)=(node(neich)%x-node(ihc)%x)            ! Vector between apical i cell node and the apical j cell node
             v(2)=(node(neich)%y-node(ihc)%y)
             v(3)=(node(neich)%z-node(ihc)%z)	
             modv=sqrt( (v(1))**2+(v(2))**2+(v(3))**2)   ! modulus v
             u(1)=(node(node(ihc)%altre)%x-node(ihc)%x)  ! Apical-basal vector i cell  
             u(2)=(node(node(ihc)%altre)%y-node(ihc)%y)
             u(3)=(node(node(ihc)%altre)%z-node(ihc)%z)
             modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)    ! modulus u
             upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)           ! vectorial product u * v
             angle=upv/(modu*modv)                       ! getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
             sumangle=sumangle+angle                     ! sum of angles of the neighbouts in a cell i 
             numn=numn+1                                 ! Counter of the real number of neights (just apical nodes) of the cell i
          enddo
          if (numn.ne.0) then               ! Just in case a cell has no neights
             napicali=napicali+1            ! This is a counter of the number of cells i with neights
                                            ! it will be used after to get the probability of classes
             medianglee=sumangle/real(numn) !This obtains the average angle of the cell i with its neightbours
             avangles(ihc)=medianglee
          endif  
       enddo

       order=0


       ! This divides the average angles in classes each 0.1
       do ihc=1,nd
          if (avangles(ihc)==666.0d0)cycle
          do jhc=lowlimit,highlimit, steps   ! The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
             if (avangles(ihc)*scaling.ge.jhc)then 
                if (avangles(ihc)*scaling<jhc+1)then
	            classangles(ihc)=jhc
	         endif
             endif
          enddo
       enddo


       !!! COOCURRENCE !!!

       !!! This is the neightbour finding/ cooccurrence calculating loop!!!
       v=0.0d0
       order=0; order2=0  

       do khc=lowlimit, highlimit, steps                         !! Go over all the classes
          order=order+1                                          !! This is the class i counter 
          do lhc=lowlimit, highlimit, steps                      !! Go over all the combinations of clases
             order2=lhc+11                                        !! This is the class j counter          
             if (khc>lhc)cycle                                      !! but we consider i,j the same category as j,i
             do ihc=1,nd                                         !! go over all nodes
                if (classangles(ihc).ne.khc)cycle                !! find those whose class is equal to khc
                do jhc=1, fourth_nneighs(ihc)                            !! go over all their neighbors
                   neich=fourth_neighs(ihc,jhc)                          !! this is the neighbor
                   if (neich>ihc)cycle                           !! We do not want to consider neighbor interaction twice (although the result should be the same)
                   if (classangles(neich).ne.lhc)cycle             !! check whether the angle class of neich is lhc, if yes we've got a match!!
                   cooccurrences(order,order2)=cooccurrences(order,order2)+1.0d0  !! write down the match in the coocurrences matrix
                   totalneightshc=totalneightshc+1.0d0           !! counter of the total number of coocurrences of all angle classes
                enddo
             enddo
          enddo
       enddo



       !The probabilities are the number of neights with a given combination of curvatures between the total num of neights
       cooprobs=cooccurrences/totalneightshc



       if (tiponodhc==1)then
          !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)
          do ihc=1, nclasseshc
             do jhc=1, nclasseshc
                if (cooprobs(ihc,jhc)==0.0d0)cycle
                pij=cooprobs(ihc,jhc)
                entropy1=entropy1+(pij*log(pij)/log(2.0000000000))
	     enddo
          enddo
          pij=0
          entropy1=-entropy1      
       else
          !This calculates the joint entropy as -sum( pij*log(pij)) (Sole 2004)
          do ihc=1, nclasseshc
             do jhc=1, nclasseshc
                if (cooprobs(ihc,jhc)==0.0d0)cycle
                pij=cooprobs(ihc,jhc)
                entropy2=entropy2+(pij*log(pij)/log(2.0000000000))
	     enddo
          enddo
          pij=0
          entropy2=-entropy2 
       endif
    enddo

    entropy=(entropy1+entropy2)/2.0d0
    fitelli=entropy                                                   !!>> HC 11-11-2021

  end subroutine                                                        !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  
  subroutine surface_volume_ratio(fitelli)                                                       !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                                !!>> HC 11-11-2021 and if calculates the surface/volume ratio 
    integer ::  ihc, jhc, khc, nnodes, naltech, surface, volume, outside, total, inside          !!>> HC 11-11-2021 of a morphology using a volume integration method
    integer ::  iihc,jjhc,kkhc, ord, lhc                                                         !!>> HC 11-11-2021
    real*8  :: suvol, fitelli                                                                    !!>> HC 11-11-2021
    real*8 :: upv, sumd, average_flux_resistance, modu                                           !!>> HC 11-11-2021
    integer, allocatable, dimension(:,:,:) :: bfillz, bfillz2, bfillx, bfillx2, bfilly, bfilly2  !!>> HC 11-11-2021
    call extrem                                                                                  !!>> HC 11-11-2021

    a=2*maxval(node(:nd)%add)                                           !!>> HC 11-11-2021 maximal interaction distance between mesenchymal nodes
    rv=a+1d-3                                                           !!>> HC 11-11-2021
    urv=1.0d0/a                                                         !!>> HC 11-11-2021

    nboxes=nint(extre*urv)+int(maxval(node(:nd)%dmo)+dmax)+1            !!>> HC 11-11-2021
    if (allocated(boxes)) deallocate(boxes)                             !!>> HC 11-11-2021
    allocate(boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))       !!>> HC 11-11-2021
    boxes=0                                                             !!>> HC 11-11-2021
    do ihc=1,nd                                                         !!>> HC 11-11-2021
      if(node(ihc)%tipus.ne.1)cycle                                     !!>> HC 11-11-2021
      iihc=nint(node(ihc)%x*urv);jjhc=nint(node(ihc)%y*urv);kkhc=nint(node(ihc)%z*urv)  !!>> HC 11-11-2021
      boxes(iihc,jjhc,kkhc)=ihc                                         !!>> HC 11-11-2021
    end do                                                              !!>> HC 11-11-2021

    if(allocated(bfillz))deallocate(bfillz)                             !!>> HC 11-11-2021 This saves whether the box is
    allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))      !!>> HC 11-11-2021 1=surface; 2=inside (volume); 0=outside
    bfillz=666                                                          !!>> HC 11-11-2021

    if(allocated(bfillz2))deallocate(bfillz2)                           !!>> HC 11-11-2021
    allocate(bfillz2(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))     !!>> HC 11-11-2021
    bfillz2=666                                                         !!>> HC 11-11-2021
    surface=0; volume=0; total=0; outside=0                             !!>> HC 11-11-2021

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!  Z AXIS  !!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      !!>> HC 11-11-2021
    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021 This goes over all the boxes in
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021 one of the directions of the z axis
          inside=0                                                      !!>> HC 11-11-2021
          do khc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             if (boxes(ihc,jhc,khc)>0)then                              !!>> HC 11-11-2021 if the box is full of nodes, then it is surface
                bfillz(ihc,jhc,khc)=1                                   !!>> HC 11-11-2021
                if (inside==0)then; inside=1; else; inside=0; endif     !!>> HC 11-11-2021 and if it is the firs one we find, then we start counting
             else                                                       !!>> HC 11-11-2021 the volume (inner boxes)
                if (inside==1)then                                      !!>> HC 11-11-2021 If the box is empty, it can be
                   bfillz(ihc,jhc,khc)=2                                !!>> HC 11-11-2021 inside the embryo = volume
                else                                                    !!>> HC 11-11-2021
                   bfillz(ihc,jhc,khc)=0                                !!>> HC 11-11-2021 outside the embryo = outer space
                endif                                                   !!>> HC 11-11-2021
             endif                                                      !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021


    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021  This goes over all the boxes in
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021  the other direction of the z axis
          inside=0; ord=0                                               !!>> HC 11-11-2021  an it does the same as above
          do lhc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             khc=nboxes-ord                                             !!>> HC 11-11-2021
             ord=ord+1                                                  !!>> HC 11-11-2021
             if (boxes(ihc,jhc,khc)>0)then                              !!>> HC 11-11-2021
                bfillz2(ihc,jhc,khc)=1                                  !!>> HC 11-11-2021
                if (inside==0)then; inside=1; else; inside=0; endif     !!>> HC 11-11-2021
             else                                                       !!>> HC 11-11-2021
                if (inside==1)then                                      !!>> HC 11-11-2021
                   bfillz2(ihc,jhc,khc)=2                               !!>> HC 11-11-2021
                else                                                    !!>> HC 11-11-2021
                   bfillz2(ihc,jhc,khc)=0                               !!>> HC 11-11-2021
                endif                                                   !!>> HC 11-11-2021
             endif                                                      !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!  CALCULATIONS  !!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                !!>> HC 11-11-2021

    do ihc=-nboxes,nboxes                                               !!>> HC 11-11-2021
       do jhc=-nboxes,nboxes                                            !!>> HC 11-11-2021
          do khc=-nboxes,nboxes                                         !!>> HC 11-11-2021
             total=total+1                                              !!>> HC 11-11-2021
             if(bfillz2(ihc,jhc,khc)==0) bfillz(ihc,jhc,khc)=0          !!>> HC 11-11-2021 the outer space has to be the same in both directions
             if(bfillz(ihc,jhc,khc)==0) outside=outside+1               !!>> HC 11-11-2021 we count the boxes in each category
             if(bfillz(ihc,jhc,khc)==2) volume=volume+1                 !!>> HC 11-11-2021
             if(bfillz(ihc,jhc,khc)==1) surface=surface+1               !!>> HC 11-11-2021
          enddo                                                         !!>> HC 11-11-2021
       enddo                                                            !!>> HC 11-11-2021
    enddo                                                               !!>> HC 11-11-2021
    
    if (volume>100)then                                                 !!>> HC 11-11-2021
       suvol=real(surface)/real(volume)                                 !!>> HC 11-11-2021 SURFACE/VOLUME
       fitelli=suvol                                                    !!>> HC 11-11-2021
    else                                                                !!>> HC 11-11-2021
       fitelli=0.0d0                                                    !!>> HC 11-11-2021 This is completely broken
    endif                                                               !!>> HC 11-11-2021

  end subroutine                                                        !!>> HC 11-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine curvature_patch_count(fitelli)                                                   !!>> HC 11-11-2021 This subroutine was created by Hugo Cano
    Implicit none                                                                           !!>> HC 11-11-2021 and it calculates the local joint entropy

    real*8 :: angle, upv, modu, modv, sumangle, medianglee, entropy, entropy1, disthc
    real*8 :: entropy2, pij, pij_pipj, entropy3,  rnapical, totalneightshc, fitelli
    integer ::  napical, napicali, ihc, jhc, khc, lhc, mhc, order, order2, node2_1, numn, ndouble, neichi
    integer:: tempclass, lowlimit, highlimit, steps, tiponodhc, scaling, nclasseshc, neich, classtemsum
    integer, allocatable, dimension(:) :: classangles, oldorder, dummycol, patched
    integer, allocatable, dimension(:,:) :: double_neighs, trans_double_neighs
    integer, allocatable, dimension(:) :: double_nneighs, trans_double_nneighs
    real*8, dimension(1:3) :: u, v
    real*8, allocatable, dimension(:) :: avangles, probclass, xcoord, ycoord, zcoord
    real*8, allocatable, dimension(:,:) :: cooccurrences, cooprobs, classprod
    integer :: CPC1, CPC2, minpatch, nneich, patchsize
    real :: CPC

 

    disthc=5     ! Distance of the neightbours
    tiponodhc=2  ! Type of node we are aiming 2= apical 1 = basal
    lowlimit=-10 ! Lower limit of the interval in classes
    highlimit=9  ! Higher limit of the interval in classes
    steps=1      ! width of the interval in classes
    nclasseshc=20!Number of classes
    scaling=10   ! Scaling of the avangle vector to fit into the classes
    ndouble=maxval(nneigh(1:nd))*2
    minpatch=3

    ! This calculates the number of apical cells
    napical=0
    do ihc=1,nd                            !This loop calculates the number of apical nodes that are going to be inspected
       if (node(ihc)%tipus==tiponodhc)then !Select only the apical nodes
          napical=napical+1                !Counter of the number of apical nodes inspected
       endif
    enddo

     ! this sets the dimension of the vector of average angles to the number of apical nodes

    if(allocated(patched))deallocate(patched)
    allocate(patched(nd)) 
     
    if(allocated(avangles))deallocate(avangles)
    allocate(avangles(nd)) 

    if(allocated(classangles))deallocate(classangles)
    allocate(classangles(nd)) 

    if(allocated(probclass))deallocate(probclass)
    allocate(probclass(nclasseshc)) 

    if(allocated(cooccurrences))deallocate(cooccurrences)
    allocate(cooccurrences(nclasseshc,nclasseshc)) 
  
    if(allocated(cooprobs))deallocate(cooprobs)
    allocate(cooprobs(nclasseshc,nclasseshc)) 

    if(allocated(classprod))deallocate(classprod)
    allocate(classprod(nclasseshc,nclasseshc)) 

    if(allocated(double_neighs))deallocate(double_neighs)
    allocate(double_neighs(nd,ndouble)) 

    if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
    allocate(trans_double_neighs(nd,ndouble)) 

    if(allocated(double_nneighs))deallocate(double_nneighs)
    allocate(double_nneighs(nd)) 
  
    if(allocated(trans_double_nneighs))deallocate(trans_double_nneighs)
    allocate(trans_double_nneighs(napical)) 
  
    entropy=0.0d0; entropy1=0.0d0; entropy2=0.0d0; entropy3=0.0d0
    CPC=0.0d0;    CPC1=0;    CPC2=0

    do tiponodhc=1,2
       avangles=666.0d0; classangles=666; probclass=0; napicali=0; rnapical=napical
       cooccurrences=0; cooprobs=0.0d0; tempclass=0; patched=0
       totalneightshc=0.0d0; double_nneighs=0; double_neighs=0
       pij_pipj=0; pij=0; order=0  

       !!Finding second order neighbors (can be optimized...)
       order2=0; order=0
       do ihc=1,nd
          if(node(ihc)%tipus.ne.tiponodhc)cycle
          order=0
          do jhc=1, nneigh(ihc)
             neich=neigh(ihc,jhc)
             if(neich==ihc)cycle
             if(node(neich)%tipus.ne.tiponodhc)cycle
             do khc=1,nneigh(neich)
                neichi=neigh(neich,khc)
                if(node(neichi)%tipus.ne.tiponodhc)cycle
                if(neichi==neich)cycle
                if(neichi==ihc)cycle
                if(any(neigh(ihc,1:nneigh(ihc))==neichi))cycle
                if (order>0)then
                   if(any(double_neighs(ihc,1:double_nneighs(ihc))==neichi))cycle
                endif
                order=order+1
                if (order>ndouble)then
                   if(allocated(trans_double_neighs))deallocate(trans_double_neighs)
                   allocate(trans_double_neighs(nd,ndouble))
                   trans_double_neighs=0
                   trans_double_neighs(1:nd,1:ndouble)=double_neighs(1:nd,1:ndouble)  
                   ndouble=ndouble+1    
                   if(allocated(double_neighs))deallocate(double_neighs)
                   allocate(double_neighs(nd,ndouble))
                   double_neighs=0
                   double_neighs(1:nd,1:ndouble-1)=trans_double_neighs(1:nd,1:ndouble-1)    
                endif
                double_neighs(ihc,order)=neichi
                double_nneighs(ihc)=double_nneighs(ihc)+1
             enddo
          enddo
       enddo

   !!! This is th neightbour finding/angle calculating loop!!!

       do ihc=1,nd                                       ! Inspecting all the nodes i of the embryo
          sumangle=0 ; numn=0                            ! initializing      
          if (node(ihc)%tipus.ne.tiponodhc )cycle        ! We are just looking for apical cells  
          do jhc=1,double_nneighs(ihc)                           ! Comaring them with its neighbors
             neich=double_neighs(ihc,jhc)
             if(ihc==neich)cycle                         ! Not themselves
             if(node(neich)%tipus.ne.tiponodhc)cycle     ! We are just looking for apical cells 
             v(1)=(node(neich)%x-node(ihc)%x)            ! Vector between apical i cell node and the apical j cell node
             v(2)=(node(neich)%y-node(ihc)%y)
             v(3)=(node(neich)%z-node(ihc)%z)	
             modv=sqrt( (v(1))**2+(v(2))**2+(v(3))**2)   ! modulus v
             u(1)=(node(node(ihc)%altre)%x-node(ihc)%x)  ! Apical-basal vector i cell  
             u(2)=(node(node(ihc)%altre)%y-node(ihc)%y)
             u(3)=(node(node(ihc)%altre)%z-node(ihc)%z)
             modu=sqrt((u(1))**2+(u(2))**2+(u(3))**2)    ! modulus u
             upv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)           ! vectorial product u * v
             angle=upv/(modu*modv)                       ! getting the dot product between u and v ranges from -1 to 1 and is directly proportional to the angle
             sumangle=sumangle+angle                     ! sum of angles of the neighbouts in a cell i 
             numn=numn+1                                 ! Counter of the real number of neights (just apical nodes) of the cell i
          enddo
          if (numn.ne.0) then               ! Just in case a cell has no neights
             napicali=napicali+1            ! This is a counter of the number of cells i with neights
                                            ! it will be used after to get the probability of classes
             medianglee=sumangle/real(numn) !This obtains the average angle of the cell i with its neightbours
             avangles(ihc)=medianglee
          endif  
       enddo

       order=0

       ! This divides the average angles in classes each 0.1
       do ihc=1,nd
          if (avangles(ihc)==666.0d0)cycle
          do jhc=lowlimit,highlimit, steps   ! The class -10 is form -1 to -0.9 and so on... The class 9 is form 0.9 to 1
             if (avangles(ihc)*scaling.ge.jhc)then 
                if (avangles(ihc)*scaling<jhc+1)then
	            classangles(ihc)=jhc
                endif
             endif
          enddo
       enddo
       
       do jhc=lowlimit,highlimit, steps
          if (any(classangles==jhc))then
             do ihc=1,nd-1
                if (classangles(ihc).ne.jhc)cycle
                if (patched(ihc)==1)cycle
                patched(ihc)=1
                nneich=nneigh(ihc)
                patchsize=1
                do khc=ihc+1,nd
                   if (classangles(khc).ne.jhc)cycle
                   if (patched(khc)==1)cycle
                   if (any(neigh(ihc,1:nneich)==khc))then
                      patched(khc)=1
                      patchsize=patchsize+1
                   endif
                enddo
                if(patchsize>minpatch.and.tiponodhc==1) CPC1=CPC1+1
                if(patchsize>minpatch.and.tiponodhc==2) CPC2=CPC2+1
             enddo
          endif
       enddo
       
    enddo

    CPC=(real(CPC1)+real(CPC2))/2.0d0                                 !!>> HC 11-11-2021
    open(206,file="complexity_vals.dat")                                   !!>> HC 11-11-2021
         write(206,*) CPC1, CPC2     !!>> HC 11-11-2021
    close(206)                                                             !!>> HC 11-11-2021
    
    fitelli=CPC                                                   !!>> HC 11-11-2021

  end subroutine                                                        !!>> HC 11-11-2021
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 11-11-2021  
  
subroutine directional_selection(fitelli)                                                   !!>> HC 19-2-2024 This subroutine was created by Hugo Cano
    Implicit none                                                                           !!>> HC 19-2-2024 and 
   
integer ::  ihc, jhc, khc,  naltech, surface, volume, outside, total, inside, differ, centrow
integer ::  iihc,jjhc,kkhc, ord1,  lhc, prev, changes, mhc, nhc, neich, neichi, newdots, phc
real*8 :: maxcoord, fitelli
real*8 :: upv, sumd, average_flux_resistance, modu, ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, pershared
real*8 :: maxx, minx, maxy, miny, maxz, minz, maxadd, minadd, anchx, anchy, anchz, anchadd, totmax, totmin, totanch
real, allocatable, dimension(:,:) :: ncoords
integer, allocatable, dimension(:,:,:) :: bfillz, bfillz2, bfillx, bfillx2, bfilly, bfilly2 
real*8, dimension(26) :: corners
real*8, dimension(26,3) :: cornerscoord


maxx=0.0d0; minx=10000; maxy=0.0d0; miny=10000; maxz=0.0d0; minz=10000
anchx=0.0d0; anchy=0.0d0; anchz=0.0d0; totmax=0.0d0; totmin=0.0d0; totanch=0.0d0
maxadd=0.0d0


newdots=0
do ihc=1,nd
   if(node(ihc)%tipus.ne.1)cycle
   newdots=newdots+1
   do jhc=1,nneigh(ihc)
      neich=neigh(ihc,jhc)
      if(node(neich)%tipus.ne.1)cycle
      do khc=1,nneigh(neich)
         neichi=neigh(neich,khc)
         if(node(neichi)%tipus.ne.1)cycle
         if(neichi==ihc)cycle
            newdots=newdots+10
      enddo
   enddo
enddo

if(allocated(ncoords))deallocate(ncoords)
allocate(ncoords(1:newdots,1:3))
ncoords=0.0d0

ord1=0
do ihc=1,nd
   if(node(ihc)%tipus.ne.1)cycle
   ord1=ord1+1
   ncoords(ord1,1)=node(ihc)%x
   ncoords(ord1,2)=node(ihc)%y
   ncoords(ord1,3)=node(ihc)%z
enddo

do ihc=1,nd
   if(node(ihc)%tipus.ne.1)cycle
   do jhc=1,nneigh(ihc)
      neich=neigh(ihc,jhc)
      if(node(neich)%tipus.ne.1)cycle
      do khc=1,nneigh(neich)
         neichi=neigh(neich,khc)
         if(node(neichi)%tipus.ne.1)cycle
         if(neichi==ihc)cycle
            do phc=1,10
                ahc = node(neich)%x - node(ihc)%x
                bhc = node(neich)%y - node(ihc)%y
                chc = node(neich)%z - node(ihc)%z
                
                ahc2 = node(neichi)%x - node(ihc)%x
                bhc2 = node(neichi)%y - node(ihc)%y
                chc2 = node(neichi)%z - node(ihc)%z
                
                call random_number(u1); call random_number(u2);
                if ( (u1+u2) > 1.0d0)then
                   u1=1-u1
                   u2=1-u2
                endif
                
                ahc3 = u1*ahc + u2*ahc2
                bhc3 = u1*bhc + u2*bhc2
                chc3 = u1*chc + u2*chc2
                ord1=ord1+1
                ncoords(ord1,1)= node(ihc)%x + ahc3
                ncoords(ord1,2)= node(ihc)%y + bhc3
                ncoords(ord1,3)= node(ihc)%z + chc3
            enddo      
      enddo
   enddo 
enddo

do ihc=1,nd
   if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
   if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
   if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
enddo

totmax=maxx
if (maxy>totmax) totmax=maxy
if (maxz>totmax) totmax=maxz

maxadd=nodeo(1)%add*2
urv=1/maxadd
nboxes=(urv*totmax) +2

if(allocated(bfillz))deallocate(bfillz)
allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
bfillz=0
do ihc=1,ord1
   jhc=nint(ncoords(ihc,1)*urv)
   khc=nint(ncoords(ihc,2)*urv)
   nhc=nint(ncoords(ihc,3)*urv)
   bfillz(jhc,khc,nhc)=1
enddo
bfillz(nboxes,nboxes,nboxes)=2
changes=11; ord1=0
do while(changes>0)
   changes=0
   do ihc=-nboxes,nboxes
      do jhc=-nboxes,nboxes
         do khc=-nboxes,nboxes
            if( bfillz(ihc,jhc,khc).ne.2)cycle
            do nhc=-1,1
               ord1=ihc+nhc
               if(ord1>nboxes.or.ord1<(-nboxes))cycle
               if(bfillz(ord1,jhc,khc).ne.0)cycle
               bfillz(ord1,jhc,khc)=2
               changes=changes+1
            enddo
            do nhc=-1,1
               ord1=jhc+nhc
               if(ord1>nboxes.or.ord1<(-nboxes))cycle
               if(bfillz(ihc,ord1,khc).ne.0)cycle
               bfillz(ihc,ord1,khc)=2
               changes=changes+1
            enddo
            do nhc=-1,1
               ord1=khc+nhc
               if(ord1>nboxes.or.ord1<(-nboxes))cycle
               if(bfillz(ihc,jhc,ord1).ne.0)cycle
               bfillz(ihc,jhc,ord1)=2
               changes=changes+1
            enddo
            bfillz(ihc,jhc,khc)=3
         enddo
      enddo
   enddo
enddo

surface=0; volume=0; outside=0; total=0
do ihc=-nboxes,nboxes
   do jhc=-nboxes,nboxes
      do khc=-nboxes,nboxes
         if (bfillz(ihc,jhc,khc)==1)then
            surface=surface+1
         elseif(bfillz(ihc,jhc,khc)==0)then
            volume=volume+1
         else
            outside=outside+1
         endif
      enddo
   enddo
enddo

total=surface+volume

do ihc=-nboxes,nboxes
   do jhc=-nboxes,nboxes
      do khc=-nboxes,nboxes
         if ( bfillz(ihc,jhc,khc).le.1)then
            bfillz(ihc,jhc,khc)=1
         else
            bfillz(ihc,jhc,khc)=0
         endif
      enddo
   enddo
enddo

 centrow=0
do ihc=-nboxes,nboxes
   do jhc=-nboxes,nboxes
      if ( bfillz(ihc,0,jhc)==0)cycle
      centrow=centrow+1
   enddo
enddo

differ=0
do ihc=1,nboxes
   jhc=-ihc
   do khc=-nboxes,nboxes
      do nhc=-nboxes,nboxes
         if( bfillz(ihc,khc,nhc) == bfillz(jhc,khc,nhc))cycle
         differ=differ+1
      enddo
   enddo
enddo
pershared=real(differ)/real(total-centrow)


corners=0.0d0
nhc=0; mhc=0
do khc=-1,1
   do lhc=-1,1
      do jhc=-1,1
         if(khc==0.and.lhc==0.and.jhc==0)cycle
         do ihc=-nboxes,0
            if (bfillz(ihc*khc,ihc*lhc,ihc*jhc).ne.1)cycle
               nhc=nhc+1
               cornerscoord(nhc,1)=real(ihc*khc)*maxadd
               cornerscoord(nhc,2)=real(ihc*lhc)*maxadd
               cornerscoord(nhc,3)=real(ihc*jhc)*maxadd
               exit
         enddo      
      enddo
   enddo
enddo

maxcoord=0.0d0
do ihc=1,nd
   if(node(ihc)%tipus.ne.1)cycle
   if(abs(node(ihc)%x)>maxcoord) maxcoord=abs(node(ihc)%x)
   if(abs(node(ihc)%y)>maxcoord) maxcoord=abs(node(ihc)%y)
   if(abs(node(ihc)%z)>maxcoord) maxcoord=abs(node(ihc)%z)
enddo

maxcoord=1/maxcoord
do ihc=1,26
   cornerscoord(ihc,1)=cornerscoord(ihc,1)*maxcoord
   cornerscoord(ihc,2)=cornerscoord(ihc,2)*maxcoord
   cornerscoord(ihc,3)=cornerscoord(ihc,3)*maxcoord
enddo

corners=-666.0d0
do ihc=1,26
   corners(ihc)=sqrt( cornerscoord(ihc,1)**2 + cornerscoord(ihc,2)**2 + cornerscoord(ihc,3)**2  )
enddo

if (volume==0)then
   fitelli=0.0010d0
elseif(pershared>0.1500d0)then
   fitelli=0.0010d0
else
   ahc=0.0d0; bhc=0.0d0
   do ihc=1,26
      if(ihc==2)cycle
      if(ihc==10)cycle
      if(ihc==11)cycle
      if(ihc==12)cycle
      if(ihc==19)cycle
      ahc=ahc+corners(ihc)
   enddo
   bhc=corners(2)+corners(10)+corners(11)+corners(12)+corners(19)
   fitelli=bhc/(ahc+bhc)
endif

 open (123, file="individual.datfitness")
      write(123,*) fitelli
 close(123)
 
 open (123, file="rob.val")
      write(123,*) pershared
 close(123)
    
end subroutine  
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 19-2-2024
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 19-2-2024
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!>> HC 19-2-2024 

subroutine complexity(fitelli, volumenes) !!>> AL 12-3-2024
   integer                 :: ihc,jhc,khc,surface,volume,outside,total,differ,shared,ord1,changes,nhc,ijk,puntos_fill,coco,iter,cb
   integer                 :: neich,neichi,newdots,cuantos,w_sect,compa_tot,secciones,phc,comparaciones,size_box,numrowspersect,ca
   integer                 :: guardar,voll
   real*8                  :: ahc,bhc,chc,ahc2,bhc2,chc2,ahc3,bhc3,chc3,maxx,minx,maxy,miny,maxz,maxadd,totmax,crompio,conta_n_nodos
   real*8                  :: u1,u2,pershared,pershared_z,pershared_x,pershared_y,per,alfa,fitelli,cociente,crompiopokito
   character*300           :: input,output
   character*152           :: pathandfile, crear
   character(len=3)        :: str1
   character(len=80)       :: fichero
   character(len=80)       :: NAME
   logical                 :: exists
   real, allocatable, dimension(:,:)      :: ncoords,cual_caja,orden_final_x,orden_final_y,orden_final_z
   integer, allocatable, dimension(:,:,:) :: bfillz
   integer, allocatable, dimension(:)     :: volumenes
   !******************************************************************************************************!
   size_box        = 1
   numrowspersect  = 3     !Number of rows per section 
   puntos_fill     = 100
   iter            = 1
   cociente        = 0
!   call iniread
!   call readsnap(pathandfile)
!   call iniboxes
   call neighbor_build
   if(allocated(volumenes))deallocate(volumenes)
   allocate(volumenes(3)) 
   volumenes=0

   !! 1. Calcular el n de puntos que vamos a usar
   crompio = 1
   crompiopokito = 1
   newdots = 0
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      newdots=newdots+1
      do jhc=1,nneigh(ihc) !AL> this tells you how many neighbours ihc has.
         neich=neigh(ihc,jhc) !AL> this tells you the node number of neighbours
         if(node(neich)%tipus.ne.1)cycle
         do khc=1,nneigh(neich)
               neichi=neigh(neich,khc)
               if(node(neichi)%tipus.ne.1)cycle
               if(neichi==ihc)cycle
                  newdots=newdots+puntos_fill
         enddo
      enddo
   enddo

   !!!!!!!!!calculate complexity z axis with displacement!!!!!!!!!
   alfa = 0
   pershared = 0
   do ijk = 1,iter
      !! ncoords guarda los puntos que vamos a usar para describir la morph
      if(allocated(ncoords))deallocate(ncoords)
      allocate(ncoords(1:newdots,1:3)) 
      ncoords=0.0d0

      !! 2. coordenadas de los nodos de la morphología
      ord1=0
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         ord1=ord1+1
         ncoords(ord1,1) = (node(ihc)%x ) 
         ncoords(ord1,2) = (node(ihc)%y)
         ncoords(ord1,3) = (node(ihc)%z + alfa) !Desplazamiento en eje z
      enddo
      cuantos = ord1
      !! 3. Puntos extra entre las células de la morphología
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         do jhc=1,nneigh(ihc)
               neich=neigh(ihc,jhc)
               if(node(neich)%tipus.ne.1)cycle
               do khc=1,nneigh(neich)
                  neichi=neigh(neich,khc)
                  if(node(neichi)%tipus.ne.1)cycle
                  if(neichi==ihc)cycle
                  do phc=1,puntos_fill
                     ahc = node(neich)%x - node(ihc)%x
                     bhc = node(neich)%y - node(ihc)%y
                     chc = (node(neich)%z + alfa) - (node(ihc)%z + alfa)
                     
                     ahc2 = node(neichi)%x - node(ihc)%x
                     bhc2 = node(neichi)%y - node(ihc)%y
                     chc2 = (node(neichi)%z + alfa) - (node(ihc)%z + alfa)
                     
                     call random_number(u1); call random_number(u2); !random real between 0 <= 1
                     if ( (u1+u2) > 1.0d0)then
                           u1=1-u1
                           u2=1-u2
                     endif
                     
                     ahc3 = u1*ahc + u2*ahc2
                     bhc3 = u1*bhc + u2*bhc2
                     chc3 = u1*chc + u2*chc2
                     ord1=ord1+1
                     ncoords(ord1,1)= node(ihc)%x + ahc3
                     ncoords(ord1,2)= node(ihc)%y + bhc3
                     ncoords(ord1,3)= (node(ihc)%z + alfa) + chc3 
                  enddo      
               enddo
         enddo 
      enddo
      !! 4. Calcular el tamaño de las cajas
      maxx=0.0d0 
      maxy=0.0d0 
      maxz=0.0d0
      do ihc=1,nd
         if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
         if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
         if(abs(node(ihc)%z + alfa)>maxz) maxz=abs(node(ihc)%z + alfa)
      enddo

      totmax=0.0d0
      totmax=maxx
      if (maxy>totmax) totmax=maxy
      if (maxz>totmax) totmax=maxz

      maxadd=0.0d0
      maxadd=nodeo(1)%add*size_box    !! TAMAÑO DE LAS CAJAS !!ESTO ES INDEPENDIENTE DE LA MORFO YA QUE nodeo%add NO CAMBIA
      urv=1/maxadd                    !! COEFICIENTE DE NORMALIZACIÓN PARA PASAR DE COORDENADAS DE NODO A ÍNDICE DE CAJA
      nboxes=(urv*totmax) + 2         !! NUMERO DE CAJAS !!AL>> Notar el +2. los extremos de bfillz (i.e.bfillz(+/-nboxes,+/-nboxes,+/-nboxes)) no tienen morfo adentro 

      !! bfillz es la matriz con todas las cajas
      if(allocated(bfillz))deallocate(bfillz)
      allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes)) 
      bfillz=0
   
      !! Llenamos bfillz con las coordenadas de los puntos 
      !0 = dentro morfo !1 = epitelio !3 = fuera epitelio !NOTA:Inicialmente toda caja que no es epitelio = 0
      do ihc=1,ord1
         jhc=nint(ncoords(ihc,1)*urv)
         khc=nint(ncoords(ihc,2)*urv)
         nhc=nint(ncoords(ihc,3)*urv)
         bfillz(jhc,khc,nhc) = 1
      end do

      changes  = 11
      ord1     = 0
      !!5. RELLENAMOS EL EXTERIOR DE LA MORPHOLOGÍA FOREST FIRE
      bfillz(nboxes,nboxes,nboxes) = 2
      do while(changes>0) 
         changes=0
         do ihc=-nboxes,nboxes
               do jhc=-nboxes,nboxes
                  do khc=-nboxes,nboxes
                     if(bfillz(ihc,jhc,khc).ne.2)cycle
                     do nhc=-1,1
                           ord1=ihc+nhc
                           if(ord1>nboxes .or. ord1<(-nboxes))cycle
                           if(bfillz(ord1,jhc,khc) .ne. 0)cycle
                           bfillz(ord1,jhc,khc)=2
                           changes=changes+1
                     enddo
                     do nhc=-1,1
                           ord1=jhc+nhc
                           if(ord1>nboxes.or.ord1<(-nboxes))cycle
                           if(bfillz(ihc,ord1,khc).ne.0)cycle
                           bfillz(ihc,ord1,khc)=2
                           changes=changes+1
                     enddo
                     do nhc=-1,1
                           ord1=khc+nhc
                           if(ord1>nboxes.or.ord1<(-nboxes))cycle
                           if(bfillz(ihc,jhc,ord1).ne.0)cycle
                           bfillz(ihc,jhc,ord1)=2
                           changes=changes+1
                     enddo
                     bfillz(ihc,jhc,khc)=3
                  enddo
               enddo
         enddo
      enddo

      !! CONTAMOS EL NÚMERO DE CAJAS
      !! surface= cubiertas por puntos (células)
      !! outside= fuera de la morphología
      !! volume = dentro de la morphología
      surface=0; volume=0; outside=0; total=0
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if (bfillz(ihc,jhc,khc)==1)then
                  surface=surface+1
                  elseif(bfillz(ihc,jhc,khc)==0)then
                  volume=volume+1
                  else
                  outside=outside+1
                  endif
               enddo
         enddo
      enddo
      total=surface+volume
      volumenes(1) = volume
      if(volume==0)then
         crompio = 0
      end if
      cociente = surface/volume
      if(cociente>1)then
         crompiopokito = 0
      endif

      open (123, file="individual.volume.txt")
         write(123,*) volume                 !!>>AL 4-4-24: If no displacement, volume(1)=volume(2)=volume(3)
      close(123)
      !! Lo pasamos a 0=fuera; 1=dentro
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if (bfillz(ihc,jhc,khc).le.1)then
                     bfillz(ihc,jhc,khc)=1
                  else
                     bfillz(ihc,jhc,khc)=0
                  endif
               enddo
         enddo
      enddo

      !!Calculamos el volumen no compartido entre secciones en eje z !>>AL 17/01/24
      per            = 0
      compa_tot      = 0
      differ         = 0
      shared         = 0
      pershared_z    = 0
      comparaciones  = 0
      secciones      = 0
      do khc = (nboxes - 1), (-nboxes + 1 + numrowspersect), -numrowspersect !>>AL 17/01/24: top to bottom comparison
         secciones = secciones + 1
         do nhc = khc - numrowspersect, (-nboxes + 1), -numrowspersect
               do phc = 0, (numrowspersect - 1) !AL 17/01/24:this moves you through rows in sections being compared
                  if((nhc - phc) < (-nboxes+1))then; cycle; end if !AL 17/01/24:in case you are trying to compare a section with i rows against a section with j rows where i>j
                  do ihc = (nboxes - 1), (-nboxes + 1), -1 
                     do jhc = (nboxes - 1), (-nboxes + 1), -1
                           if(bfillz(ihc,jhc,khc - phc) == 0 .and. bfillz(ihc,jhc,nhc - phc) == 0 ) cycle; !AL 17/01/24>> boxes of both sections are outside morphology
                           if(bfillz(ihc,jhc,khc - phc) == 1 .and. bfillz(ihc,jhc,nhc - phc) == 1 ) then;
                              shared = shared + 1
                           else
                              differ = differ + 1
                           end if 
                     end do
                  end do
                  if(shared == 0)then; differ = 0; end if; !AL 17/01/24:this is in case you compare an empty section with a nonempty one.
               end do   
               if(shared == 0 .and. differ == 0)cycle              !>>AL 17/01/24 both sections are empty
                                                                  !!this can happen because the limits of the grid of boxes are determined
                                                                  !!by variable totmax which only considers max coordinate in one component
                                                                  !!so maybe in other component max coordinate is smaller so anything bigger is empty. also only one could be empty
               per = (real(differ) / real(shared + differ)) + per  !>>AL 17/01/24 nonsharedvolume
               comparaciones = comparaciones + 1                   
               shared = 0
               differ = 0
         end do
         if(comparaciones == 0)then;cycle; end if                !>>AL 17/01/24 this for the case in which all sections compared were empty (boundary sects)
         pershared_z = per/real(comparaciones) + pershared_z
         compa_tot = comparaciones + compa_tot
         comparaciones = 0
         per = 0
      end do

      pershared = pershared + pershared_z 
      alfa = alfa + 0.064
   end do
   pershared = pershared/iter
   pershared_z = pershared

   alfa = 0
   pershared = 0
   do ijk = 1,iter
      !! ncoords guarda los puntos que vamos a usar para describir la morph
      if(allocated(ncoords))deallocate(ncoords)
      allocate(ncoords(1:newdots,1:3)) 
      ncoords=0.0d0
      !! 2. coordenadas de los nodos de la morphología
      ord1=0
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         ord1=ord1+1
         ncoords(ord1,1) = (node(ihc)%x + alfa) !Desplazamiento en eje x
         ncoords(ord1,2) = (node(ihc)%y)
         ncoords(ord1,3) = (node(ihc)%z) 
      enddo
      cuantos = ord1
      !! 3. Puntos extra entre las células de la morphología
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         do jhc=1,nneigh(ihc)
               neich=neigh(ihc,jhc)
               if(node(neich)%tipus.ne.1)cycle
               do khc=1,nneigh(neich)
                  neichi=neigh(neich,khc)
                  if(node(neichi)%tipus.ne.1)cycle
                  if(neichi==ihc)cycle
                  do phc=1,puntos_fill
                     ahc = (node(neich)%x + alfa) - (node(ihc)%x + alfa)
                     bhc = node(neich)%y - node(ihc)%y
                     chc = node(neich)%z - node(ihc)%z
                     
                     ahc2 = (node(neichi)%x + alfa) - (node(ihc)%x + alfa)
                     bhc2 = node(neichi)%y - node(ihc)%y
                     chc2 = node(neichi)%z - node(ihc)%z 
                     
                     call random_number(u1); call random_number(u2); !random real between 0 <= 1
                     if ( (u1+u2) > 1.0d0)then
                           u1=1-u1
                           u2=1-u2
                     endif
                     
                     ahc3 = u1*ahc + u2*ahc2
                     bhc3 = u1*bhc + u2*bhc2
                     chc3 = u1*chc + u2*chc2
                     ord1=ord1+1
                     ncoords(ord1,1)= (node(ihc)%x + alfa) + ahc3
                     ncoords(ord1,2)= node(ihc)%y + bhc3
                     ncoords(ord1,3)= node(ihc)%z + chc3 
                  enddo      
               enddo
         enddo 
      enddo

      !! 4. Calcular el tamaño de las cajas
      maxx=0.0d0 
      maxy=0.0d0 
      maxz=0.0d0
      do ihc=1,nd
         if(abs(node(ihc)%x + alfa)>maxx) maxx=abs(node(ihc)%x + alfa)
         if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
         if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
      enddo

      totmax=0.0d0
      totmax=maxx
      if (maxy>totmax) totmax=maxy
      if (maxz>totmax) totmax=maxz

      maxadd=0.0d0
      maxadd=nodeo(1)%add*size_box  !! TAMAÑO DE LAS CAJAS !!ESTO ES INDEPENDIENTE DE LA MORFO YA QUE nodeo%add NO CAMBIA
      urv=1/maxadd           !! COEFICIENTE DE NORMALIZACIÓN PARA PASAR DE COORDENADAS DE NODO A ÍNDICE DE CAJA
      nboxes=(urv*totmax) +2 !! NUMERO DE CAJAS !!AL>> Notar el +2. los extremos de bfillz (i.e.bfillz(+/-nboxes,+/-nboxes,+/-nboxes)) no tienen morfo adentro 

      !! bfillz es la matriz con todas las cajas
      if(allocated(bfillz))deallocate(bfillz)
      allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes)) 
      bfillz=0
      !! Llenamos bfillz con las coordenadas de los puntos
      !0 = dentro morfo !1 = epitelio !3 = fuera epitelio
      do ihc=1,ord1
         jhc=nint(ncoords(ihc,1)*urv)
         khc=nint(ncoords(ihc,2)*urv)
         nhc=nint(ncoords(ihc,3)*urv)
         bfillz(jhc,khc,nhc) = 1
      end do

      changes  = 11
      ord1     = 0
      !!5. RELLENAMOS EL EXTERIOR DE LA MORPHOLOGÍA COMO EN EL PAINT
      bfillz(nboxes,nboxes,nboxes) = 2
      do while(changes>0)     !esto sobra
         changes=0
         do ihc=-nboxes,nboxes
               do jhc=-nboxes,nboxes
                  do khc=-nboxes,nboxes
                  if(bfillz(ihc,jhc,khc).ne.2)cycle
                  do nhc=-1,1
                     ord1=ihc+nhc
                     if(ord1>nboxes .or. ord1<(-nboxes))cycle
                     if(bfillz(ord1,jhc,khc) .ne. 0)cycle
                     bfillz(ord1,jhc,khc)=2
                     changes=changes+1
                  enddo
                  do nhc=-1,1
                     ord1=jhc+nhc
                     if(ord1>nboxes.or.ord1<(-nboxes))cycle
                     if(bfillz(ihc,ord1,khc).ne.0)cycle
                     bfillz(ihc,ord1,khc)=2
                     changes=changes+1
                  enddo
                  do nhc=-1,1
                     ord1=khc+nhc
                     if(ord1>nboxes.or.ord1<(-nboxes))cycle
                     if(bfillz(ihc,jhc,ord1).ne.0)cycle
                     bfillz(ihc,jhc,ord1)=2
                     changes=changes+1
                  enddo
                  bfillz(ihc,jhc,khc)=3
                  enddo
               enddo
         enddo
      enddo
      !! CONTAMOS EL NÚMERO DE CAJAS
      !! surface= cubiertas por puntos (células)!! outside= fuera de la morphología !! volume = dentro de la morphología
      surface=0; volume=0; outside=0; total=0
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if (bfillz(ihc,jhc,khc)==1)then
                  surface=surface+1
                  elseif(bfillz(ihc,jhc,khc)==0)then
                  volume=volume+1
                  else
                  outside=outside+1
                  endif
               enddo
         enddo
      enddo
      total=surface+volume
      volumenes(2) = volume
      if(volume==0)then
         crompio = 0
      end if
      cociente = surface/volume
      if(cociente>1)then
         crompiopokito = 0
      endif

      !! Lo pasamos a 0=fuera; 1=dentro
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if ( bfillz(ihc,jhc,khc).le.1)then
                  bfillz(ihc,jhc,khc)=1
                  else
                  bfillz(ihc,jhc,khc)=0
                  endif
               enddo
         enddo
      enddo

      !!Calculamos el volumen no compartido entre secciones en eje z !>>AL 17/01/24
      per            = 0
      compa_tot      = 0
      differ         = 0
      shared         = 0
      pershared_x    = 0
      comparaciones  = 0
      ca             = 0
      cb             = 0
      coco           = 0

      do ihc = (nboxes - 1), (-nboxes + 1 + numrowspersect), -numrowspersect !>>AL 17/01/24 top to bottom comparison
         secciones = secciones + 1
         do nhc = ihc - numrowspersect , (-nboxes + 1), -numrowspersect
               do phc = 0, (numrowspersect - 1)
                  if((nhc - phc) < (-nboxes+1))then;cycle; end if !in case you are trying to compare a section with i rows against a section with j rows where i>j
                  do khc = (nboxes - 1), (-nboxes + 1), -1
                     do jhc = (nboxes - 1), (-nboxes + 1), -1
                           if(bfillz(ihc - phc,jhc,khc) == 0 .and. bfillz(nhc - phc,jhc,khc) == 0 ) cycle; !AL 17/01/24>> boxes of both sections are outside morphology
                           if(bfillz(ihc - phc,jhc,khc) == 1 .and. bfillz(nhc - phc,jhc,khc) == 1 ) then;
                              shared = shared + 1
                           else
                              differ = differ + 1
                           end if
                     end do
                  end do
                  if(shared == 0)then;differ = 0;end if; !when one of the sections is empty differ is max 
               end do
               if(shared == 0 .and. differ == 0)cycle                !>>AL 17/01/24 both sections are empty
                                                                     !this can happen because the limits of the grid of boxes are determined
                                                                     !by variable totmax which only considers max coordinate in one component
                                                                     !so maybe in other component max coordinate is smaller so anything bigger is empty
               per = (real(differ) / real(shared + differ)) + per    !>>AL 17/01/24 nonsharedvolume
               shared = 0
               differ = 0
               comparaciones = comparaciones + 1
         end do
         if(comparaciones == 0)then;cycle; end if  !>>AL 17/01/24 this for the case in which all sections compared were empty 
         pershared_x = per/real(comparaciones) + pershared_x
         compa_tot = comparaciones + compa_tot
         comparaciones = 0
         per = 0
      end do

      pershared = pershared + pershared_x 
      alfa = alfa + 0.064
   end do 

   pershared = pershared/iter
   pershared_x = pershared

   alfa = 0
   pershared = 0
   do ijk = 1,iter
      !! ncoords guarda los puntos que vamos a usar para describir la morph
      if(allocated(ncoords))deallocate(ncoords)
      allocate(ncoords(1:newdots,1:3)) 
      ncoords=0.0d0
      !! 2. coordenadas de los nodos de la morphología
      ord1=0
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         ord1=ord1+1
         ncoords(ord1,1) = (node(ihc)%x ) 
         ncoords(ord1,2) = (node(ihc)%y + alfa)
         ncoords(ord1,3) = (node(ihc)%z) !Desplazamiento en eje z
      enddo
      cuantos = ord1
      !! 3. Puntos extra entre las células de la morphología
      do ihc=1,nd
         if(node(ihc)%tipus.ne.1)cycle
         do jhc=1,nneigh(ihc)
               neich=neigh(ihc,jhc)
               if(node(neich)%tipus.ne.1)cycle
               do khc=1,nneigh(neich)
                  neichi=neigh(neich,khc)
                  if(node(neichi)%tipus.ne.1)cycle
                  if(neichi==ihc)cycle
                  do phc=1,puntos_fill
                     ahc = node(neich)%x - node(ihc)%x
                     bhc = (node(neich)%y + alfa) - (node(ihc)%y + alfa)
                     chc = node(neich)%z - node(ihc)%z
                     
                     ahc2 = node(neichi)%x - node(ihc)%x
                     bhc2 = (node(neichi)%y + alfa) - (node(ihc)%y + alfa)
                     chc2 = node(neichi)%z - node(ihc)%z 
                     
                     call random_number(u1); call random_number(u2); !random real between 0 <= 1
                     if ( (u1+u2) > 1.0d0)then
                           u1=1-u1
                           u2=1-u2
                     endif
                     
                     ahc3 = u1*ahc + u2*ahc2
                     bhc3 = u1*bhc + u2*bhc2
                     chc3 = u1*chc + u2*chc2
                     ord1=ord1+1
                     ncoords(ord1,1)= node(ihc)%x + ahc3
                     ncoords(ord1,2)= (node(ihc)%y + alfa) + bhc3
                     ncoords(ord1,3)= node(ihc)%z + chc3 
                  enddo      
               enddo
         enddo 
      enddo

      !! 4. Calcular el tamaño de las cajas
      maxx=0.0d0 
      maxy=0.0d0 
      maxz=0.0d0
      do ihc=1,nd
         if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
         if(abs(node(ihc)%y + alfa)>maxy) maxy=abs(node(ihc)%y + alfa)
         if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
      enddo

      totmax=0.0d0
      totmax=maxx
      if (maxy>totmax) totmax=maxy
      if (maxz>totmax) totmax=maxz

      maxadd=0.0d0
      maxadd=nodeo(1)%add*size_box  !! TAMAÑO DE LAS CAJAS !!ESTO ES INDEPENDIENTE DE LA MORFO YA QUE nodeo%add NO CAMBIA
      urv=1/maxadd           !! COEFICIENTE DE NORMALIZACIÓN PARA PASAR DE COORDENADAS DE NODO A ÍNDICE DE CAJA
      nboxes=(urv*totmax) +2 !! NUMERO DE CAJAS !!AL>> Notar el +2. los extremos de bfillz (i.e.bfillz(+/-nboxes,+/-nboxes,+/-nboxes)) no tienen morfo adentro 

      !! bfillz es la matriz con todas las cajas
      if(allocated(bfillz))deallocate(bfillz)
      allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
      bfillz=0

      !! Llenamos bfillz con las coordenadas de los puntos
      !0 = dentro morfo !1 = epitelio !3 = fuera epitelio
      do ihc=1,ord1
         jhc=nint(ncoords(ihc,1)*urv)
         khc=nint(ncoords(ihc,2)*urv)
         nhc=nint(ncoords(ihc,3)*urv)
         bfillz(jhc,khc,nhc) = 1
      end do

      changes  = 11
      ord1     = 0
      !!5. RELLENAMOS EL EXTERIOR DE LA MORPHOLOGÍA COMO EN EL PAINT
      bfillz(nboxes,nboxes,nboxes) = 2
      do while(changes>0)     !esto sobra
         changes=0
         do ihc=-nboxes,nboxes
               do jhc=-nboxes,nboxes
                  do khc=-nboxes,nboxes
                     if(bfillz(ihc,jhc,khc).ne.2)cycle
                     do nhc=-1,1
                           ord1=ihc+nhc
                           if(ord1>nboxes .or. ord1<(-nboxes))cycle
                           if(bfillz(ord1,jhc,khc) .ne. 0)cycle
                           bfillz(ord1,jhc,khc)=2
                           changes=changes+1
                     enddo
                     do nhc=-1,1
                           ord1=jhc+nhc
                           if(ord1>nboxes.or.ord1<(-nboxes))cycle
                           if(bfillz(ihc,ord1,khc).ne.0)cycle
                           bfillz(ihc,ord1,khc)=2
                           changes=changes+1
                     enddo
                     do nhc=-1,1
                           ord1=khc+nhc
                           if(ord1>nboxes.or.ord1<(-nboxes))cycle
                           if(bfillz(ihc,jhc,ord1).ne.0)cycle
                           bfillz(ihc,jhc,ord1)=2
                           changes=changes+1
                     enddo
                     bfillz(ihc,jhc,khc)=3
                  enddo
               enddo
         enddo
      enddo
      !! CONTAMOS EL NÚMERO DE CAJAS
      !! surface= cubiertas por puntos (células) !! outside= fuera de la morphología !! volume = dentro de la morphología
      surface=0; volume=0; outside=0; total=0
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if (bfillz(ihc,jhc,khc)==1)then
                  surface=surface+1
                  elseif(bfillz(ihc,jhc,khc)==0)then
                  volume=volume+1
                  else
                  outside=outside+1
                  endif
               enddo
         enddo
      enddo
      total=surface+volume
      volumenes(3) = volume
      if(volume==0)then
         crompio = 0
      end if
      cociente = surface/volume
      if(cociente>1)then
         crompiopokito = 0
      endif
      !! Lo pasamos a 0=fuera; 1=dentro
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
               do khc=-nboxes,nboxes
                  if ( bfillz(ihc,jhc,khc).le.1)then
                  bfillz(ihc,jhc,khc)=1
                  else
                  bfillz(ihc,jhc,khc)=0
                  endif
               enddo
         enddo
      enddo

      !!Calculamos el volumen no compartido entre secciones en eje y !>>AL 17/01/24
      per            = 0
      compa_tot      = 0
      differ         = 0
      shared         = 0
      pershared_y    = 0
      comparaciones  = 0    
      do jhc = (nboxes - 1), (-nboxes + 1 + numrowspersect), -numrowspersect !>>AL 17/01/24 top to bottom comparison
         secciones = secciones + 1
         do nhc = jhc - numrowspersect , (-nboxes + 1), -numrowspersect
               do phc = 0, (numrowspersect - 1)
                  if((nhc - phc) < (-nboxes+1))then; cycle; end if !in case you are trying to compare a section with i rows against a section with j rows where i>j
                  do ihc = (nboxes - 1), (-nboxes + 1), -1
                  do khc = (nboxes - 1), (-nboxes + 1), -1
                     if(bfillz(ihc,jhc - phc,khc) == 0 .and. bfillz(ihc,nhc - phc,khc) == 0 ) cycle; !AL 17/01/24>> boxes of both sections are outside morphology
                     if(bfillz(ihc,jhc - phc,khc) == 1 .and. bfillz(ihc,nhc - phc,khc) == 1 ) then;
                           shared = shared + 1
                     else
                           differ = differ + 1
                     end if
                  end do
                  end do
                  if(shared == 0)then;differ = 0;end if; !when one of the sections is empty differ is max 
               end do
               if(shared == 0 .and. differ == 0)cycle                !>>AL 17/01/24 both sections are empty
                                                                  !this can happen because the limits of the grid of boxes are determined
                                                                  !by variable totmax which only considers max coordinate in one component
                                                                  !so maybe in other component max coordinate is smaller so anything bigger is empty
               per = (real(differ) / real(shared + differ)) + per    !>>AL 17/01/24 nonsharedvolume
               shared = 0
               differ = 0
               comparaciones = comparaciones + 1
         end do
         if(comparaciones == 0)then;cycle; end if  !>>AL 17/01/24 this for the case in which all sections compared were empty 
         pershared_y = per/real(comparaciones) + pershared_y
         compa_tot = comparaciones + compa_tot
         comparaciones = 0
         per = 0
      end do

      pershared = pershared + pershared_y 
      alfa = alfa + 0.064
   end do 

   pershared = pershared/iter
   pershared_y = pershared

   fitelli = 0
   fitelli = (pershared_x +pershared_y+pershared_z)/3
   if(crompio==0)fitelli=0
   if(crompiopokito==0)then 
      fitelli=0                                        !!>>AL 4-4-24 this is in case volumen is not 0 but there is a whole in morphology's surface
      !volumenes(1) = 1
      !volumenes(2) = 1
      !volumenes(3) = 1
   end if
   print *, "We are here!!"
   print *, "fitelli:", fitelli 
   open (123, file="individual.datfitness")
         write(123,*) fitelli
   close(123)
   
end subroutine 

subroutine trait_distance(fitelli)
   Implicit none

   integer ::  ihc, jhc, khc, surface, volume, outside, total, landmark_number
   integer ::  ord1,  lhc, changes, mhc, nhc, neich, neichi, newdots, phc
   real*8  :: maxcoord, centroidx, centroidy, centroidz, centroid_size, sumx, sumy, sumz
   real*8  ::  ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, EMD,totmax2
   real*8  :: maxx, maxy,  maxz,  maxadd, fitelli
   real*8, dimension(26)   :: corners, landmarks
   real*8, dimension(26,3) :: cornerscoord
   real*8, dimension(26)   :: target_corners, distances
   real*8, dimension(26,3) :: target_landmarks
   real, allocatable, dimension(:,:)      :: ncoords
   integer, allocatable, dimension(:,:,:) :: bfillz 
   !******************************************************************************************************!
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
   
   !centering landmarks around cero
   do ihc=1,nd
      node(ihc)%x = node(ihc)%x - centroidx
      node(ihc)%y = node(ihc)%y - centroidy
      node(ihc)%z = node(ihc)%z - centroidz
   end do 
   
   newdots=0
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      newdots=newdots+1
      do jhc=1,nneigh(ihc)
         neich=neigh(ihc,jhc)
         if(node(neich)%tipus.ne.1)cycle
         do khc=1,nneigh(neich)
            neichi=neigh(neich,khc)
            if(node(neichi)%tipus.ne.1)cycle
            if(neichi==ihc)cycle
               newdots=newdots+10
         enddo
      enddo
   enddo
   
   if(allocated(ncoords))deallocate(ncoords)
   allocate(ncoords(1:newdots,1:3))
   ncoords=0.0d0
   
   ord1=0
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      ord1=ord1+1
      ncoords(ord1,1)=node(ihc)%x
      ncoords(ord1,2)=node(ihc)%y
      ncoords(ord1,3)=node(ihc)%z
   enddo
   
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      do jhc=1,nneigh(ihc)
         neich=neigh(ihc,jhc)
         if(node(neich)%tipus.ne.1)cycle
         do khc=1,nneigh(neich)
            neichi=neigh(neich,khc)
            if(node(neichi)%tipus.ne.1)cycle
            if(neichi==ihc)cycle
               do phc=1,10
                   ahc = node(neich)%x - node(ihc)%x
                   bhc = node(neich)%y - node(ihc)%y
                   chc = node(neich)%z - node(ihc)%z
                   
                   ahc2 = node(neichi)%x - node(ihc)%x
                   bhc2 = node(neichi)%y - node(ihc)%y
                   chc2 = node(neichi)%z - node(ihc)%z
                   
                   call random_number(u1); call random_number(u2);
                   if ( (u1+u2) > 1.0d0)then
                      u1=1-u1
                      u2=1-u2
                   endif
                   
                   ahc3 = u1*ahc + u2*ahc2
                   bhc3 = u1*bhc + u2*bhc2
                   chc3 = u1*chc + u2*chc2
                   ord1=ord1+1
                   ncoords(ord1,1)= node(ihc)%x + ahc3
                   ncoords(ord1,2)= node(ihc)%y + bhc3
                   ncoords(ord1,3)= node(ihc)%z + chc3
               enddo      
         enddo
      enddo 
   enddo
   
   maxx=0;maxy=0;maxz=0;totmax2=0
   do ihc=1,nd
      if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
      if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
      if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
   enddo
   
   totmax2=maxx
   if (maxy>totmax2) totmax2=maxy
   if (maxz>totmax2) totmax2=maxz
   
   maxadd=nodeo(1)%add*2
   urv=1/maxadd
   nboxes=(urv*totmax2) +2
   
   if(allocated(bfillz))deallocate(bfillz)
   allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
   bfillz=0
   
   do ihc=1,ord1
      jhc=nint(ncoords(ihc,1)*urv)
      khc=nint(ncoords(ihc,2)*urv)
      nhc=nint(ncoords(ihc,3)*urv)
      bfillz(jhc,khc,nhc)=1
   enddo
   
   bfillz(nboxes,nboxes,nboxes)=2
   changes=11; ord1=0
   do while(changes>0)
      changes=0
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
            do khc=-nboxes,nboxes
               if( bfillz(ihc,jhc,khc).ne.2)cycle
               do nhc=-1,1
                  ord1=ihc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ord1,jhc,khc).ne.0)cycle
                  bfillz(ord1,jhc,khc)=2
                  changes=changes+1
               enddo
               do nhc=-1,1
                  ord1=jhc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ihc,ord1,khc).ne.0)cycle
                  bfillz(ihc,ord1,khc)=2
                  changes=changes+1
               enddo
               do nhc=-1,1
                  ord1=khc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ihc,jhc,ord1).ne.0)cycle
                  bfillz(ihc,jhc,ord1)=2
                  changes=changes+1
               enddo
               bfillz(ihc,jhc,khc)=3
            enddo
         enddo
      enddo
   enddo
   
   surface=0; volume=0; outside=0; total=0
   do ihc=-nboxes,nboxes
      do jhc=-nboxes,nboxes
         do khc=-nboxes,nboxes
            if (bfillz(ihc,jhc,khc)==1)then
               surface=surface+1
            elseif(bfillz(ihc,jhc,khc)==0)then
               volume=volume+1
            else
               outside=outside+1
            endif
         enddo
      enddo
   enddo
   
   total=surface+volume

   open (123, file="individual.volume.txt")
      write(123,*) volume                 !!>>AL 4-4-24: If no displacement, volume(1)=volume(2)=volume(3)
   close(123)

   !>>AL 13-9-2024: get landmark coordinates in cartesian coordinate system
   corners=0.0d0
   nhc=0; mhc=0
   do khc=-1,1
      do lhc=-1,1
         do jhc=-1,1
            if(khc==0.and.lhc==0.and.jhc==0)cycle
            do ihc=-nboxes,0
               if (bfillz(ihc*khc,ihc*lhc,ihc*jhc).gt.1)cycle
                  nhc=nhc+1
                  cornerscoord(nhc,1)=real(ihc*khc)*maxadd
                  cornerscoord(nhc,2)=real(ihc*lhc)*maxadd
                  cornerscoord(nhc,3)=real(ihc*jhc)*maxadd
                  exit
            enddo      
         enddo
      enddo
   enddo
   
   open(10, file="target_1.dat") !this are already normalized
   do ihc=1,26
      read(10,*) target_landmarks(ihc,1), target_landmarks(ihc,2), target_landmarks(ihc,3)
   end do 
   close(10)

   landmarks = 1
   landmark_number = 0

   do ihc=1,26
      if(target_landmarks(ihc,1) == 0 .and. target_landmarks(ihc,2) == 0 &
      .and. target_landmarks(ihc,3) == 0)then 
         landmarks(ihc) = 0; 
      else 
         landmark_number = landmark_number + 1; 
      endif
   end do 

   !>>AL 12-9-2024 Centroid
   !centroidx = 0;centroidy = 0;centroidz = 0
   !do ihc=1,26
   !   if(landmarks(ihc) == 0)then; cycle; endif
   !   centroidx = cornerscoord(ihc,1) + centroidx
   !   centroidy = cornerscoord(ihc,2) + centroidy
   !   centroidz = cornerscoord(ihc,3) + centroidz
   !end do 
   
   !centroidx = centroidx/real(landmark_number)
   !centroidy = centroidy/real(landmark_number)
   !centroidz = centroidz/real(landmark_number)

   !print *, "Centroids before centering:"
   !print *, "x:",centroidx,"y:",centroidy,"z:",centroidz
   
   !centering landmarks around cero
   !do ihc=1,26
   !   if(landmarks(ihc) == 0)then; cycle; endif
   !   cornerscoord(ihc,1) = cornerscoord(ihc,1) - centroidx
   !   cornerscoord(ihc,2) = cornerscoord(ihc,2) - centroidy
   !   cornerscoord(ihc,3) = cornerscoord(ihc,3) - centroidz
   !end do 
   
   sumx = 0;sumy = 0;sumz = 0
   do ihc=1, 26
      if(landmarks(ihc) == 0)then; cycle; endif
      sumx = sumx + (cornerscoord(ihc,1))**2
      sumy = sumy + (cornerscoord(ihc,2))**2
      sumz = sumz + (cornerscoord(ihc,3))**2
   end do 
   
   !>>AL 12-9-2024 Centroid size 
   centroid_size = sqrt(sumx + sumy + sumz)
   
   cornerscoord(:,1) = cornerscoord(:,1)/centroid_size
   cornerscoord(:,2) = cornerscoord(:,2)/centroid_size
   cornerscoord(:,3) = cornerscoord(:,3)/centroid_size
   
   maxcoord=0.0d0
   corners=6666.0d0
   !normalizing distances of landmarks with respect to centroid size
   !corners = corners/centroid_size

   distances = 0
   do ihc=1,26
      if(landmarks(ihc) == 0)then; cycle; endif
      distances(ihc)=sqrt((cornerscoord(ihc,1)-target_landmarks(ihc,1))**2 + (cornerscoord(ihc,2) -&
      target_landmarks(ihc,2))**2 &
      +(cornerscoord(ihc,3) -target_landmarks(ihc,3))**2)
   enddo

   EMD = sum(distances)
   if(EMD==0.0d0) EMD=0.0010d0
   fitelli = 1.0d0/EMD

   open(616,file="rob.val")
   if (volume==0)then
      fitelli=0.0d0
      write(616,*) 0.0d0
   else
      write(616,*) 1.0d0
   endif
   close(616)

   open (123, file="individual.datfitness")
         write(123,*) fitelli
   close(123)
end subroutine trait_distance

subroutine trait_distancenormalization(fitelli) !AL 11-11-24: this version is deprecated but part of legacy now
   Implicit none

   integer ::  ihc, jhc, khc, surface, volume, outside, total
   integer ::  ord1,  lhc, changes, mhc, nhc, neich, neichi, newdots, phc
   real*8  :: maxcoord, centroidx, centroidy, centroidz, centroid_size, sumx, sumy, sumz
   real*8  ::  ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, EMD,totmax2
   real*8  :: maxx, maxy,  maxz,  maxadd, fitelli,sumdist
   real*8, dimension(26)   :: corners
   real*8, dimension(26,3) :: cornerscoord
   real*8, dimension(26)   :: distances
   real*8, dimension(14)   :: target_distances,diff_distances
   real*8, dimension(26)   :: target_corners
   real, allocatable, dimension(:,:)      :: ncoords
   integer, allocatable, dimension(:,:,:) :: bfillz 
   !******************************************************************************************************!
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
   
   !centering landmarks around cero
   do ihc=1,nd
      node(ihc)%x = node(ihc)%x - centroidx
      node(ihc)%y = node(ihc)%y - centroidy
      node(ihc)%z = node(ihc)%z - centroidz
   end do 
   
   newdots=0
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      newdots=newdots+1
      do jhc=1,nneigh(ihc)
         neich=neigh(ihc,jhc)
         if(node(neich)%tipus.ne.1)cycle
         do khc=1,nneigh(neich)
            neichi=neigh(neich,khc)
            if(node(neichi)%tipus.ne.1)cycle
            if(neichi==ihc)cycle
               newdots=newdots+10
         enddo
      enddo
   enddo
   
   if(allocated(ncoords))deallocate(ncoords)
   allocate(ncoords(1:newdots,1:3))
   ncoords=0.0d0
   
   ord1=0
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      ord1=ord1+1
      ncoords(ord1,1)=node(ihc)%x
      ncoords(ord1,2)=node(ihc)%y
      ncoords(ord1,3)=node(ihc)%z
   enddo
   
   do ihc=1,nd
      if(node(ihc)%tipus.ne.1)cycle
      do jhc=1,nneigh(ihc)
         neich=neigh(ihc,jhc)
         if(node(neich)%tipus.ne.1)cycle
         do khc=1,nneigh(neich)
            neichi=neigh(neich,khc)
            if(node(neichi)%tipus.ne.1)cycle
            if(neichi==ihc)cycle
               do phc=1,10
                   ahc = node(neich)%x - node(ihc)%x
                   bhc = node(neich)%y - node(ihc)%y
                   chc = node(neich)%z - node(ihc)%z
                   
                   ahc2 = node(neichi)%x - node(ihc)%x
                   bhc2 = node(neichi)%y - node(ihc)%y
                   chc2 = node(neichi)%z - node(ihc)%z
                   
                   call random_number(u1); call random_number(u2);
                   if ( (u1+u2) > 1.0d0)then
                      u1=1-u1
                      u2=1-u2
                   endif
                   
                   ahc3 = u1*ahc + u2*ahc2
                   bhc3 = u1*bhc + u2*bhc2
                   chc3 = u1*chc + u2*chc2
                   ord1=ord1+1
                   ncoords(ord1,1)= node(ihc)%x + ahc3
                   ncoords(ord1,2)= node(ihc)%y + bhc3
                   ncoords(ord1,3)= node(ihc)%z + chc3
               enddo      
         enddo
      enddo 
   enddo
   
   maxx=0;maxy=0;maxz=0;totmax2=0
   do ihc=1,nd
      if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)
      if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)
      if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)
   enddo
   
   totmax2=maxx
   if (maxy>totmax2) totmax2=maxy
   if (maxz>totmax2) totmax2=maxz
   
   maxadd=nodeo(1)%add*2
   urv=1/maxadd
   nboxes=(urv*totmax2) +2
   
   if(allocated(bfillz))deallocate(bfillz)
   allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
   bfillz=0
   
   do ihc=1,ord1
      jhc=nint(ncoords(ihc,1)*urv)
      khc=nint(ncoords(ihc,2)*urv)
      nhc=nint(ncoords(ihc,3)*urv)
      bfillz(jhc,khc,nhc)=1
   enddo
   
   bfillz(nboxes,nboxes,nboxes)=2
   changes=11; ord1=0
   do while(changes>0)
      changes=0
      do ihc=-nboxes,nboxes
         do jhc=-nboxes,nboxes
            do khc=-nboxes,nboxes
               if( bfillz(ihc,jhc,khc).ne.2)cycle
               do nhc=-1,1
                  ord1=ihc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ord1,jhc,khc).ne.0)cycle
                  bfillz(ord1,jhc,khc)=2
                  changes=changes+1
               enddo
               do nhc=-1,1
                  ord1=jhc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ihc,ord1,khc).ne.0)cycle
                  bfillz(ihc,ord1,khc)=2
                  changes=changes+1
               enddo
               do nhc=-1,1
                  ord1=khc+nhc
                  if(ord1>nboxes.or.ord1<(-nboxes))cycle
                  if(bfillz(ihc,jhc,ord1).ne.0)cycle
                  bfillz(ihc,jhc,ord1)=2
                  changes=changes+1
               enddo
               bfillz(ihc,jhc,khc)=3
            enddo
         enddo
      enddo
   enddo
   
   surface=0; volume=0; outside=0; total=0
   do ihc=-nboxes,nboxes
      do jhc=-nboxes,nboxes
         do khc=-nboxes,nboxes
            if (bfillz(ihc,jhc,khc)==1)then
               surface=surface+1
            elseif(bfillz(ihc,jhc,khc)==0)then
               volume=volume+1
            else
               outside=outside+1
            endif
         enddo
      enddo
   enddo
   
   total=surface+volume

   open (123, file="individual.volume.txt")
      write(123,*) volume                 !!>>AL 4-4-24: If no displacement, volume(1)=volume(2)=volume(3)
   close(123)

   !>>AL 13-9-2024: get landmark coordinates in cartesian coordinate system
   corners=0.0d0
   nhc=0; mhc=0
   do khc=-1,1
      do lhc=-1,1
         do jhc=-1,1
            if(khc==0.and.lhc==0.and.jhc==0)cycle
            do ihc=-nboxes,0
               if (bfillz(ihc*khc,ihc*lhc,ihc*jhc).gt.1)cycle
                  nhc=nhc+1
                  cornerscoord(nhc,1)=real(ihc*khc)*maxadd
                  cornerscoord(nhc,2)=real(ihc*lhc)*maxadd
                  cornerscoord(nhc,3)=real(ihc*jhc)*maxadd
                  exit
            enddo      
         enddo
      enddo
   enddo
   
   open(10, file="target_1.dat") !this distances should be already normalized
      read(10,*) target_distances
   close(10)

!the following is not correct. the actual reference system is different. check old template
   sumdist        = 0
   distances      = 0                                    
   distances(5)   = sqrt(sum((cornerscoord(5,:)**2)))
   distances(22)  = sqrt(sum((cornerscoord(22,:)**2)))
   distances(13)  = sqrt(sum((cornerscoord(13,:)**2)))
   distances(14)  = sqrt(sum((cornerscoord(14,:)**2)))
   distances(11)  = sqrt(sum((cornerscoord(11,:)**2)))
   distances(18)  = sqrt(sum((cornerscoord(18,:)**2)))
   distances(1)   = sqrt(sum((cornerscoord(1,:)**2)))
   distances(3)   = sqrt(sum((cornerscoord(3,:)**2)))
   distances(20)  = sqrt(sum((cornerscoord(20,:)**2)))
   distances(16)  = sqrt(sum((cornerscoord(16,:)**2)))
   distances(24)  = sqrt(sum((cornerscoord(24,:)**2)))
   distances(7)   = sqrt(sum((cornerscoord(7,:)**2)))
   distances(9)   = sqrt(sum((cornerscoord(9,:)**2)))
   distances(26)  = sqrt(sum((cornerscoord(26,:)**2)))
   
   do ihc=1,14
      sumdist = sumdist + distances(ihc)
   end do 

   distances = distances/sumdist

   diff_distances = 0d0
   diff_distances(1) = abs(target_distances(1)- distances(5))
   diff_distances(2) = abs(target_distances(2)- distances(22))
   diff_distances(3) = abs(target_distances(3)- distances(13))
   diff_distances(4) = abs(target_distances(4)- distances(14))
   diff_distances(5) = abs(target_distances(5)- distances(11))
   diff_distances(6) = abs(target_distances(6)- distances(18))
   diff_distances(7) = abs(target_distances(7)- distances(1))
   diff_distances(8) = abs(target_distances(8)- distances(3))
   diff_distances(9) = abs(target_distances(9)- distances(20))
   diff_distances(10) = abs(target_distances(10)- distances(16))
   diff_distances(11) = abs(target_distances(11)- distances(24))
   diff_distances(12) = abs(target_distances(12)- distances(7))
   diff_distances(13) = abs(target_distances(13)- distances(9))
   diff_distances(14) = abs(target_distances(14)- distances(26))

   EMD = sum(diff_distances)
   fitelli = 1.0d0-EMD

   open(616,file="rob.val")
   if (volume==0)then
      fitelli=0.0d0
      write(616,*) 0.0d0
   else
      write(616,*) 1.0d0
   endif
   close(616)

   open (123, file="individual.datfitness")
         write(123,*) fitelli
   close(123)

end subroutine trait_distancenormalization

subroutine traits_with_proyections(fitelli) !AL 11-11-24
   
   implicit none 
   logical :: signo
   integer :: ihc,jhc,khc,lhc,contador,contadora
   real*8  :: size_box,maxx,maxy,maxz,totmax,linex,liney,linez,comp,t,dot_PA,dot_vv,mag,dlta,sumdist
   real*8  :: x0(3), x2(3),mini(3),dire(3),projection(3),this(3),closest_projection(3),distances(14)
   real*8  :: traits_morpho(26,3), fitelli,emd,centroidx,centroidy,centroidz
   character*152 :: path, input, output
   real, allocatable, dimension(:,:)  :: ncoords,cual_caja,orden_final_x,orden_final_y,orden_final_z
   real, allocatable, dimension(:,:)  :: nodos_ep
   integer, allocatable, dimension(:) :: indx
   real*8, dimension(14) :: target_distances,diff_distances

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
            !print*,"ihc:",ihc,"jhc:",jhc,"khc:",khc
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

   open(10, file="target_1.dat") !this distances should be already normalized
      read(10,*) target_distances
   close(10)

   EMD = (sum(abs(target_distances(:) - distances(:))))

   fitelli = 1.0d0-EMD

   open (123, file="individual.datfitness")
      write(123,*) fitelli
   close(123)

   open (124, file="individual.volume.txt")
      write(124,*) 0.0
   close(124)

   !Checar problema de necesitar un fichero llamado individual.volume.txt. Escribirlo al calcular la asimetria bilateral.
   !open (123, file="individual.volume.txt")
   !   write(123,*) volume                 !!>>AL 4-4-24: If no displacement, volume(1)=volume(2)=volume(3)
   !close(123)

   open(616,file="rob.val") !checar por que necesitas esto. nuevamente, podria solo escribirlo en la parta de asimetria bilateral
      write(616,*) 0.00111   
   close(616)

end subroutine traits_with_proyections

subroutine traits_with_proyections_CS(fitelli) !AL 11-11-24
   
   implicit none 
   logical :: signo
   integer :: ihc,jhc,khc,lhc,contador,contadora
   real*8  :: size_box,maxx,maxy,maxz,totmax,linex,liney,linez,comp,t,dot_PA,dot_vv,mag,dlta,sumdist
   real*8  :: x0(3), x2(3),mini(3),dire(3),projection(3),this(3),closest_projection(3),distances(14)
   real*8  :: traits_morpho(26,3), fitelli,emd,centroidx,centroidy,centroidz,suma,CS,dis
   character*152 :: path, input, output
   real, allocatable, dimension(:,:)  :: ncoords,cual_caja,orden_final_x,orden_final_y,orden_final_z
   real, allocatable, dimension(:,:)  :: nodos_ep
   integer, allocatable, dimension(:) :: indx
   !real*8, dimension(14) :: target_distances,diff_distances
   real*8, dimension(14,3) :: morpho_traits, target_traits
   real*8, dimension(14)   :: trait_magnitudes

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

   morpho_traits(1,:)  = traits_morpho(13,:)     !(0,0,+) D1 
   morpho_traits(2,:)  = traits_morpho(11,:)     !(0,+,0) D2
   morpho_traits(3,:)  = traits_morpho(14,:)     !(0,0,-) D3 
   morpho_traits(4,:)  = traits_morpho(16,:)     !(0,-,0) D4
   morpho_traits(5,:)  = traits_morpho(1,:)      !(+,+,+) D5
   morpho_traits(6,:)  = traits_morpho(3,:)      !(+,+,-) D6
   morpho_traits(7,:)  = traits_morpho(5,:)      !(+,0,0) D7
   morpho_traits(8,:)  = traits_morpho(9,:)      !(+,-,-) D8
   morpho_traits(9,:)  = traits_morpho(7,:)      !(+,-,+) D9
   morpho_traits(10,:) = traits_morpho(18,:)     !(-,+,+) D10
   morpho_traits(11,:) = traits_morpho(20,:)     !(-,+,-) D11
   morpho_traits(12,:) = traits_morpho(22,:)     !(-,0,0) D12
   morpho_traits(13,:) = traits_morpho(26,:)     !(-,-,-) D13
   morpho_traits(14,:) = traits_morpho(24,:)     !(-,-,+) D14

   !calculate centroid size 
   suma = 0
   do ihc=1,14
      suma = suma + sum(morpho_traits(ihc,:)**2)
   end do

   CS = sqrt(suma)

  !normalize
  morpho_traits = morpho_traits/CS

   !open(10, file=target) important
   open(10, file="target_1.dat") !this distances should be already normalized
      do ihc = 1, 14
        read(10, *) target_traits(ihc, 1), target_traits(ihc, 2), target_traits(ihc, 3)  ! Read each row
      end do
   close(10)

   dis = 0
   do ihc = 1, 14
      dis = dis + sqrt(sum((target_traits(ihc,:) - morpho_traits(ihc,:))**2))
   end do

   fitelli = dis

   open (123, file="individual.datfitness")
      write(123,*) fitelli
   close(123)

   open (124, file="individual.volume.txt")
      write(124,*) 0.0
   close(124)

   open(616,file="rob.val") !checar por que necesitas esto. nuevamente, podria solo escribirlo en la parta de asimetria bilateral
      write(616,*) 0.00111   
   close(616)

   !save trait magnitudes
   do ihc = 1, 14
      trait_magnitudes(ihc) = sqrt(sum(morpho_traits(ihc,:)**2))
   end do

   open(616,file="traits_local.dat")
      write(616,*) trait_magnitudes   
   close(616)

end subroutine traits_with_proyections_CS

subroutine intersection(neigh0,neigh1,neigh2,Tr,res,intrs)
   ! Inputs
   integer, intent(in) :: neigh0, neigh1, neigh2
   real(8), intent(in) :: Tr(3)
   ! Outputs
   integer, intent(out) :: res
   real(8), intent(out) :: intrs(3)

   ! Locals
   real(8) :: p0(3),p1(3),p2(3),lab(3),p01(3),p02(3),sc1,sc2
   real(8) :: t,u,v,cr1(3),cr(3)

   ! Möller–Trumbore ray–triangle intersection algorithm
   p0(1) = node(neigh0)%x; p0(2) = node(neigh0)%y; p0(3) = node(neigh0)%z
   p1(1) = node(neigh1)%x; p1(2) = node(neigh1)%y; p1(3) = node(neigh1)%z
   p2(1) = node(neigh2)%x; p2(2) = node(neigh2)%y; p2(3) = node(neigh2)%z

   p01(1)=p1(1) - p0(1); p01(2)=p1(2) - p0(2); p01(3)=p1(3) - p0(3)
   p02(1)=p2(1) - p0(1); p02(2)=p2(2) - p0(2); p02(3)=p2(3) - p0(3)
   
   ! intersection between Tr and plane formed by p0,p1,p2 is given by lab*t

   !p01 x p02 = n (plane normal)
   cr(1)=p01(2)*p02(3) - p01(3)*p02(2)
   cr(2)=p01(3)*p02(1) - p01(1)*p02(3)
   cr(3)=p01(1)*p02(2) - p01(2)*p02(1)

   sc1 = dot_product(-p0,cr)
   sc2 = dot_product(-Tr,cr)

   !if sc2 = 0 then the plane is parallel or (maybe) collinear to trait direction
   if(sc2 .eq. 0)then 
      print*,"WARNING: no interesection"
      res = 0
      open(616,file="nointeresection.dat") !AL: 25-9-25 
         write(616,*) "nointeresection.dat"   
      close(616)
   end if 

   t = sc1/sc2
   intrs = t*Tr

   !getting u 
   !p02 x -Tr
   cr1(1)=p02(2)*-Tr(3) - p02(3)*-Tr(2)
   cr1(2)=p02(3)*-Tr(1) - p02(1)*-Tr(3)
   cr1(3)=p02(1)*-Tr(2) - p02(2)*-Tr(1)

   sc1 = dot_product(-p0,cr1)
   sc2 = dot_product(-Tr,cr)
   
   u = sc1/sc2

   !getting v
   !-Tr x p01
   cr1(1)=Tr(2)*-p01(3) - Tr(3)*-p01(2)
   cr1(2)=Tr(3)*-p01(1) - Tr(1)*-p01(3)
   cr1(3)=Tr(1)*-p01(2) - Tr(2)*-p01(1)

   sc1 = dot_product(-p0,cr1)
   sc2 = dot_product(-Tr,cr)
   
   v = sc1/sc2

   if(((u + v) .le. 1) .and. (u >= 0.d0 .and. u <= 1.d0 .and. v >= 0.d0 .and. v <= 1.d0))then
      !print*,"TRIANGULO!"
      res = 1
   else if (u >= 0.d0 .and. u <= 1.d0 .and. v >= 0.d0 .and. v <= 1.d0) then
      !print*,"PARALELOGRAMO!"
      res = 2
   else
      res = 0
   end if

   return 

end subroutine intersection

subroutine traits_intersection(fitelli) !AL 25-9-25

   implicit none 

   logical :: signo,triangulo
   integer :: ihc,jhc,khc,lhc,mhc,nhc,whc,contador,contadora,vecino,idd,neigh0,neigh1,neigh2
   integer :: mini_nd(26),res,cou
   real*8  :: size_box,maxx,maxy,maxz,totmax,linex,liney,linez,comp,t,dot_PA,dot_vv,dlta,sumdist
   real*8  :: x0(3),x2(3),maxii(3),dire(3),projection(3),this(3),closest_projection(3),distances(26)
   real*8  :: traits_morpho(26,3),centroidx,centroidy,centroidz,suma,CS,dis,t_prom,t0
   real*8  :: p0(3),p1(3),p2(3),Tr(3),intrs(3),mag,fitelli
   character*152 :: path
   real, allocatable, dimension(:,:) :: nudos,nodos_ep,nudos1
   real, dimension(26,8) :: nudos_gud
   real*8, dimension(14,3) :: optimum_traits, target_traits
   real*8, dimension(14)   :: trait_magnitudes

   !>>AL 12-9-2024 Centroid
   centroidx = 0;centroidy = 0;centroidz = 0
   do ihc=1,nd
      centroidx = centroidx + node(ihc)%x
      centroidy = centroidy + node(ihc)%y
      centroidz = centroidz + node(ihc)%z
   end do
   
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
   end do

   totmax=0.0d0
   totmax=maxx
   if (maxy>totmax) totmax=maxy
   if (maxz>totmax) totmax=maxz

   totmax = totmax + 0.1*totmax !with this value we'll generate the trait lines. They will be always longer than the morphology
   
   dlta=nodeo(1)%add*2                             !AL 19-3-25 nodeo%add is node%add at dev. time = 0, so it is always the same
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
            cou = 0
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

               neigh0 = lhc
               Tr = x2

               do mhc = 1,nneigh(neigh0)                       !AL> this tells you how many neighbours lhc has.
                  neigh1 = neigh(neigh0,mhc)                    !AL> this tells you the node number of neighbours
                  if(node(neigh1)%tipus.ne.1)cycle
                  do nhc = mhc+1,nneigh(neigh0)
                     neigh2 = neigh(neigh0,nhc) 
                     if(node(neigh2)%tipus.ne.1)cycle
                     
                     res = 0
                     call intersection(neigh0,neigh1,neigh2,Tr,res,intrs)

                     if(res == 1 .or. res == 2) cou = cou + 1

                  end do 
               end do
            end do

            allocate(nudos(cou,5))
            allocate(nudos1(cou,3))

            cou = 0
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

               neigh0 = lhc
               Tr = x2

               do mhc = 1,nneigh(neigh0)                       !AL> this tells you how many neighbours lhc has.
                  neigh1 = neigh(neigh0,mhc)                    !AL> this tells you the node number of neighbours
                  if(node(neigh1)%tipus.ne.1)cycle
                  do nhc = mhc+1,nneigh(neigh0)
                     neigh2 = neigh(neigh0,nhc) 
                     if(node(neigh2)%tipus.ne.1)cycle

                     res = 0
                     call intersection(neigh0,neigh1,neigh2,Tr,res,intrs)

                     if(res == 1 .or. res == 2) then 
                        cou = cou + 1
                        mag = sqrt(dot_product(intrs, intrs)) !guardar interesccion
                        nudos(cou,1) = neigh0
                        nudos(cou,2) = neigh1
                        nudos(cou,3) = neigh2
                        nudos(cou,4) = res
                        nudos(cou,5) = mag     !test
                        nudos1(cou,1) = intrs(1)
                        nudos1(cou,2) = intrs(2)
                        nudos1(cou,3) = intrs(3)
                     end if 
                  end do 
               end do
            end do 

            comp = 0
            triangulo = .false.

            if(cou == 0)then
               print*,"PELIGRO PELIGRO PELIGRO: NO HAY INTERESECCION TRIAN PARALEL"
               open(616,file="nointeresection.dat",position="append") !AL: 25-9-25 
                  write(616,*) "nointeresection res = 0"   
               close(616)
            end if 
            do whc = 1, cou
               if(whc == 1)then 
                  do nhc = 1, cou                           !we prioritize triangle intersection 
                     if(int(nudos(nhc,4)) .ne. 1)cycle 
                     triangulo = .true.
                     if(nudos(nhc,5) > comp) then 
                        comp = nudos(nhc,5) 
                        nudos_gud(contadora,1) = nudos(nhc,1)  
                        nudos_gud(contadora,2) = nudos(nhc,2) 
                        nudos_gud(contadora,3) = nudos(nhc,3) 
                        nudos_gud(contadora,4) = nudos(nhc,4) 
                        nudos_gud(contadora,5) = nudos(nhc,5)
                        nudos_gud(contadora,6) = nudos1(nhc,1)
                        nudos_gud(contadora,7) = nudos1(nhc,2)
                        nudos_gud(contadora,8) = nudos1(nhc,3)

                     end if
                  end do 
               end if 

               if(triangulo .eqv. .true.) exit

               if(nudos(whc,5) > comp) then 
                  comp = nudos(whc,5) 
                  nudos_gud(contadora,1) = nudos(whc,1)  
                  nudos_gud(contadora,2) = nudos(whc,2) 
                  nudos_gud(contadora,3) = nudos(whc,3) 
                  nudos_gud(contadora,4) = nudos(whc,4) 
                  nudos_gud(contadora,5) = nudos(whc,5)
                  nudos_gud(contadora,6) = nudos1(nhc,1)
                  nudos_gud(contadora,7) = nudos1(nhc,2)
                  nudos_gud(contadora,8) = nudos1(nhc,3)
               end if 

            end do
      
            deallocate(nudos)
            deallocate(nudos1)

            maxii(1) = node(int(nudos_gud(contadora,1)))%x 
            maxii(2) = node(int(nudos_gud(contadora,1)))%y
            maxii(3) = node(int(nudos_gud(contadora,1)))%z

            traits_morpho(contadora,1) = dire(1)*nudos_gud(contadora,5)
            traits_morpho(contadora,2) = dire(2)*nudos_gud(contadora,5)
            traits_morpho(contadora,3) = dire(3)*nudos_gud(contadora,5)

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

   open(10, file="target_1.dat") !this distances should be already normalized
      do ihc = 1, 14
        read(10, *) target_traits(ihc, 1), target_traits(ihc, 2), target_traits(ihc, 3)  ! Read each row
      end do
   close(10)
   
   dis = 0
   do ihc = 1, 14
      dis = dis + sqrt(sum((target_traits(ihc,:) - optimum_traits(ihc,:))**2))
   end do

   fitelli = dis

   open (123, file="individual.datfitness", status="replace", action="write")
      write(123,*) dis
   close(123)

   !save trait magnitudes
   do ihc = 1, 14
      trait_magnitudes(ihc) = sqrt(sum(optimum_traits(ihc,:)**2))
   end do

   open(616,file="traits_local.dat") !AL: 25-9-25 not normalized
      write(616,*) trait_magnitudes   
   close(616)

   open (124, file="individual.volume.txt") !AL: 29-9-25 this shouldn't be necessary but haven't had time to change
      write(124,*) 0.0
   close(124)


end subroutine traits_intersection

end module  