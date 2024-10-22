!    GNOMO software (General Node Model)
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




! Is 20-5-13
module ic	!the mesenchimal and epithelial i.c. subroutines are now unified. 
use general
use genetic
use shell ! miguel4-11-13
use death
!use nexus

integer  :: nodecela
integer  :: l, ciclo, cont
real*8   :: x, xx, y, yy, z, zz, de, di, hip              
real*8   :: dx1, dx2, dx3, dy1, dy2, dy3
integer, public  :: radi,radicel,layer,mradi,mradicel,xradi,xlayer
real*8, public   ::  zepi,zmes,radius,zx
integer, public  :: packed    !>>>Miquel10-1-14
integer          :: iccag
integer,public,allocatable :: mesradicel(:)

contains

!*************************************************************************************

subroutine default_values  !IS 2-1-14

    ndmax=100000
    ecmmax=0.25d0
    realtime=0
    ttalone=10
    mnn=500
    dif_req=1.0d0
    min_comp=-1.0d-1
    screen_radius=1.0d0
    reqmin=0.05d0 
    df_reqmax=0.50d0
    deltamin=1.0d-3
    prec=0.1  ! Is 26-8-14
    angletor=0.00 !Miquel15-9-14
    ldi=epsilod
    maxad=15.0d0 !!>> HC 16-6-2020

end subroutine


!*************************************************************************************

subroutine blastula_fibonacci                                            !!>>HC 18-3-2021 This subroutine creates an epithelial sphere full of mesenchymal cells
implicit none                                                            !!>>HC 18-3-2021 using the Fibonacci algorithm (made by Hugo Cano 18-3-2021)
integer ::  ich, jch, ord, ord2, newcell                                 !!>>HC 18-3-2021
real*8 :: radch, reqch, reqch2, dich, addch, addch2                      !!>>HC 18-3-2021
real*8 :: ratio_add_eqd, ratio_req2_req1                                 !!>>HC 18-3-2021
real*8 :: goldenmean, goldenangle                                        !!>>HC 18-3-2021
real*8, allocatable, dimension(:) :: lonch, latch, xch, ych, zch         !!>>HC 18-3-2021
real*8, allocatable, dimension(:) :: mlonch, mlatch, mxch, mych, mzch    !!>>HC 18-3-2021
integer :: nmesch, ncelepch, nepch, naltrech                             !!>>HC 18-3-2021
real*8 :: xmax, xmax1, ymax, zmax, xmin, ymin, zmin, xmin1               !!>>HC 18-3-2021
integer :: npilech, pilech                                               !!>>HC 19-3-2021
real*8 :: parpil, dismech, dismch                                        !!>>HC 19-3-2021

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!!!!!!CONSTANTS!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
goldenmean=(sqrt(5.0d0)+1.0)/2.0                                         !!>>HC 18-3-2021 Golden mean, used to calculate the golden angle
goldenangle=(2.0-goldenmean)*(2.0*pi)                                    !!>>HC 18-3-2021 Golden angle, used to calculate the spheric coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021 of the blastula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!PARAMETERS OF THE BLASTULA!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
radch=1.20d0                                                             !!>>HC 18-3-2021 blastula's radius
reqch=0.12                                                               !!>>HC 18-3-2021 EQD for external epitelial nodes and mesenchyme (type 1 and 3)
ratio_add_eqd=1.50d0                                                     !!>>HC 18-3-2021 ADD/EQD of nodes
ratio_req2_req1=0.90d0                                                   !!>>HC 19-3-2021 Ration EQD external nodes/EQD internal nodes
dich=0.90d0                                                              !!>>HC 18-3-2021 distance between upper and lower epitelial nodes
nmesch=0!300                                                               !!>>HC 18-3-2021 number of mesenchymal cells
npilech=3                                                                !!>>HC 19-3-2021 Number of mesenchymal layers
dismch=0.20d0                                                            !!>>HC 19-3-2021 distance between mesenchimal layers
ncelepch=1000                                                             !!>>HC 18-3-2021 number of epithelial cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!!!!DERIVED PARAMETERS!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021

addch=ratio_add_eqd*reqch                                                !!>>HC 18-3-2021 ADD for external epitelial nodes and mesenchyme (type 1 and 3)
reqch2=ratio_req2_req1*reqch                                             !!>>HC 19-3-2021 EQD for internal epitelial nodes (type 2)
addch2=ratio_add_eqd*reqch2                                              !!>>HC 19-3-2021 ADD for internal epitelial nodes (type 2)
nepch=ncelepch*2                                                         !!>>HC 18-3-2021 Number of epitelial nodes
ncels=ncelepch+nmesch                                                    !!>>HC 18-3-2021 Total number of cells
ncals=ncels+10                                                           !!>>HC 18-3-2021 number of cells +10
nd=nepch+nmesch                                                          !!>>HC 18-3-2021 total number of nodes
nda=nd+10                                                                !!>>HC 18-3-2021 number of nodes+10
pilech=nmesch/npilech                                                    !!>>HC 19-3-2021 Number of mesenchymal cells per layer
parpil=dich/(real(npilech)+2.0d0)                                        !!>>HC 19-3-2021 Distance between layers (percentage from the center)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!ALLOCATING MEMORY!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021

if(allocated(xch)) deallocate(xch); allocate(xch(1:ncelepch))            !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(ych)) deallocate(ych); allocate(ych(1:ncelepch))            !!>>HC 18-3-2021 cartesian coordinates
if(allocated(zch)) deallocate(zch); allocate(zch(1:ncelepch))            !!>>HC 18-3-2021
if(allocated(lonch)) deallocate(lonch); allocate(lonch(1:ncelepch))      !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(latch)) deallocate(latch); allocate(latch(1:ncelepch))      !!>>HC 18-3-2021 spheric coordinates (latch=latitude lonch=longitude)
lonch=0.0d0; latch=0.0d0                                                 !!>>HC 18-3-2021
xch=0.0d0; ych=0.0d0; zch=0.0d0                                          !!>>HC 18-3-2021

if(allocated(mxch)) deallocate(mxch); allocate(mxch(1:pilech))           !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(mych)) deallocate(mych); allocate(mych(1:pilech))           !!>>HC 18-3-2021 cartesian coordinates
if(allocated(mzch)) deallocate(mzch); allocate(mzch(1:pilech))           !!>>HC 18-3-2021
if(allocated(mlonch)) deallocate(mlonch); allocate(mlonch(1:pilech))     !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(mlatch)) deallocate(mlatch); allocate(mlatch(1:pilech))     !!>>HC 18-3-2021 spheric coordinates (latch=latitude lonch=longitude)
mlonch=0.0d0; mlatch=0.0d0                                               !!>>HC 18-3-2021
mxch=0.0d0; mych=0.0d0; mzch=0.0d0                                       !!>>HC 18-3-2021


if (allocated(node))deallocate(node); allocate(node(1:nda))              !!>>HC 18-3-2021 
if (allocated(cels))deallocate(cels); allocate(cels(1:ncals))            !!>>HC 18-3-2021
if(allocated(nodeo)) deallocate(nodeo); allocate(nodeo(nda))             !!>>HC 18-3-2021
call iniarrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
   !******* #2 DEFINING MODEL PARAMETERS *******                         !!>>HC 18-3-2021 Here we set the global parameters so that the program can run

    !implementation and initializations                                  !!>>HC 18-3-2021
    getot=0                                                              !!>>HC 18-3-2021
    itacc=0                                                              !!>>HC 18-3-2021
    nparti=1000                                                          !!>>HC 18-3-2021
    idum=-11111                                                          !!>>HC 18-3-2021
    idumoriginal=idum                                                    !!>>HC 18-3-2021
    nvarglobal_out=5                                                     !!>>HC 18-3-2021

    !physical                                                            !!>>HC 18-3-2021
    temp=0.1d1 !low value is low temperature	                           !!>>HC 18-3-2021
    desmax=0.001                                                         !!>>HC 18-3-2021
    resmax=1d-3                                                          !!>>HC 18-3-2021
    prop_noise=0.10d0                                                    !!>>HC 18-3-2021
    deltamax=1d-2 ! miguel 14-10-13                                      !!>>HC 18-3-2021
    dmax=1                                                               !!>>HC 18-3-2021
    screen_radius=1.0d0                                                  !!>>HC 18-3-2021
    
    maxbox=10.0d0                                                        !!>>HC 19-3-2021
    maxad=5.0d0                                                          !!>>HC 18-3-2021
    maxcycl=1.50d0                                                       !!>>HC 18-3-2021

    ffu=0                                                                !!>>HC 18-3-2021 
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                        !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                        !!>>HC 18-3-2021


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                           !!>>HC 18-3-2021
!!!!!!!!!NODE!!!!!!!!!!!!!!!!!                                           !!>>HC 18-3-2021 HERE WE FILL THE NODE MATRIX
!!!!!NODE TYPES AND ALTRE/CELL ASIGNATION!!!!!!!!!!                      !!>>HC 18-3-2021

newcell=0                                                                !!>>HC 18-3-2021 We decide node type, altres and cells
do ich=1,nd                                                              !!>>HC 18-3-2021 
   if (ich.le.nepch)then                                                 !!>>HC 18-3-2021 This is true for all EPITELIAL nodes
      if (modulo(ich,2).ne.0)then                                        !!>>HC 18-3-2021 This is true for odd nodes (type 2 nodes go always first)
         node(ich)%tipus=2                                               !!>>HC 18-3-2021 So this is an odd epitelial node (1,3,5...) and it must be tipus 2
         node(ich)%altre=ich+1                                           !!>>HC 18-3-2021 Its altre node will be the next node (tipus 1)
         newcell=newcell+1                                               !!>>HC 18-3-2021 Tipus 2 always goes first, so this node starts a new cell
         node(ich)%icel=newcell                                          !!>>HC 18-3-2021 Assign the cell to de node
      else                                                               !!>>HC 18-3-2021
         node(ich)%tipus=1                                               !!>>HC 18-3-2021 This node is an even node (2,4,6..) and it must be tipus 1
         node(ich)%altre=ich-1                                           !!>>HC 18-3-2021 Its altre must be the previous node (tipus 1)
         node(ich)%icel=newcell                                          !!>>HC 18-3-2021 and it is in the same cell than his mate (ie the previous node)
      endif                                                              !!>>HC 18-3-2021
   else                                                                  !!>>HC 18-3-2021 This is for MESENCHIMAL nodes
      newcell=newcell+1                                                  !!>>HC 18-3-2021 Every node belongs to a single cell
      node(ich)%tipus=3                                                  !!>>HC 18-3-2021 The typus of mesenchymal cells
      node(ich)%icel=newcell                                             !!>>HC 18-3-2021 Assing the cell
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

!!!!!NODE PROPERTIES!!!!!!!!!!                                           !!>>HC 18-3-2021
do ich=1,nd                                                              !!>>HC 18-3-2021 Here we give the nodes PHYSICAL PROPERTIES
   if (node(ich)%tipus==1.or.node(ich)%tipus==3)then                     !!>>HC 18-3-2021 External epithelial nodes and mesenchymal nodes have the same size
      node(ich)%eqd=reqch                                                !!>>HC 18-3-2021 
      node(ich)%add=addch                                                !!>>HC 18-3-2021 
   else                                                                  !!>>HC 18-3-2021
      node(ich)%eqd=reqch2                                               !!>>HC 18-3-2021 
      node(ich)%add=addch2                                               !!>>HC 18-3-2021 
   endif                                                                 !!>>HC 18-3-2021
   node(ich)%grd=node(ich)%eqd                                           !!>>HC 18-3-2021
   node(ich)%adh=5.0d0; node(ich)%rec=30.0d0                             !!>>HC 18-3-2021
   node(ich)%erp=80.0d0; node(ich)%est=100.0d0                               !!>>HC 18-3-2021
   node(ich)%eqs=0.150d0; node(ich)%hoo=10.0d0                           !!>>HC 18-3-2021
   node(ich)%mov=0.010d0; node(ich)%dmo=0.0010d0                         !!>>HC 18-3-2021
   node(ich)%pla=0.0d0                                                   !!>>HC 18-3-2021
   node(ich)%cod=0.0d0;  node(ich)%kvol=0.0d0; node(ich)%pld=0.0d0       !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021
  
!!!!!NODE POSITIONS!!!!!!!!!!                                            !!>>HC 18-3-2021 HERE WE DECIDE NODE POSITIONS (FIBONACCI SPHERE)
do ich=1,ncelepch                                                        !!>>HC 18-3-2021 EPITHELUM
   latch(ich)=asin(-1.0+2.0*real(ich)/(real(ncelepch)+1))                !!>>HC 18-3-2021 This is the latitude calculated with the golden angle (fibonacci sphere) 
   lonch(ich)=goldenangle*real(ich)                                      !!>>HC 18-3-2021 This is the longitude calculated with the golden angle (fibonacci sphere) 
   xch(ich)=cos(lonch(ich))*cos(latch(ich))                              !!>>HC 18-3-2021 This transforms spherical coordinates to cartesian coordinates
   ych(ich)=sin(lonch(ich))*cos(latch(ich))                              !!>>HC 18-3-2021
   zch(ich)=sin(latch(ich))                                              !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

ord=0                                                                    !!>>HC 18-3-2021 We give the cartesian coordinates to the external nodes
do ich=1,nd                                                              !!>>HC 18-3-2021
   if (node(ich)%tipus==1)then                                           !!>>HC 18-3-2021
      ord=ord+1                                                          !!>>HC 18-3-2021
      node(ich)%x=xch(ord)*radch                                         !!>>HC 18-3-2021 Multiply the coordinate by the radius, otherwise, radius will be =1
      node(ich)%y=ych(ord)*radch                                         !!>>HC 18-3-2021
      node(ich)%z=zch(ord)*radch                                         !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

do ich=1,nd                                                              !!>>HC 18-3-2021 We give the cartesian coordinates to the internal nodes
   if (node(ich)%tipus==2)then                                           !!>>HC 18-3-2021 Basically we multiply the external radius coordinates by the distance between spheres (dich)
      naltrech=node(ich)%altre                                           !!>>HC 18-3-2021
      node(ich)%x=node(naltrech)%x*dich                                  !!>>HC 18-3-2021
      node(ich)%y=node(naltrech)%y*dich                                  !!>>HC 18-3-2021
      node(ich)%z=node(naltrech)%z*dich                                  !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

do ich=1,pilech                                                          !!>>HC 18-3-2021 MESENCHYMAL NODES
   mlatch(ich)=asin(-1.0+2.0*real(ich)/(real(pilech)+1))                 !!>>HC 18-3-2021 We create the Fibonacci sphere 
   mlonch(ich)=goldenangle*real(ich)                                     !!>>HC 18-3-2021 (same as epithelium, but with different node number)
   mxch(ich)=cos(mlonch(ich))*cos(mlatch(ich))                           !!>>HC 18-3-2021
   mych(ich)=sin(mlonch(ich))*cos(mlatch(ich))                           !!>>HC 18-3-2021
   mzch(ich)=sin(mlatch(ich))                                            !!>>HC 18-3-2021
enddo

ord=nepch                                                                !!>>HC 19-3-2021
do ich=1,npilech                                                         !!>>HC 19-3-2021 We assign the coordinates to each spherical layer of mesenchymal cells
   dismech=(radch*dich)-(dismch*ich)                                     !!>>HC 19-3-2021 This is the radius of each mesenchymal spherical layer
   ord2=0                                                                !!>>HC 19-3-2021
   do jch=1,pilech                                                       !!>>HC 19-3-2021
      ord=ord+1                                                          !!>>HC 19-3-2021
      node(ord)%x=mxch(jch)*dismech                                      !!>>HC 19-3-2021
      node(ord)%y=mych(jch)*dismech                                      !!>>HC 19-3-2021
      node(ord)%z=mzch(jch)*dismech                                      !!>>HC 19-3-2021
   enddo                                                                 !!>>HC 19-3-2021
enddo                                                                    !!>>HC 19-3-2021

nodeo=node                                                               !!>>HC 19-3-2021 set the original node properties


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  !!>>HC 18-3-2021
!!!!!!!!!!!!CELLS!!!!!!!!!!!!!!!!!!!!!!                                  !!>>HC 18-3-2021 FILL THE CELS MATRIX

newcell=0                                                                !!>>HC 18-3-2021
do ich=1,nd                                                              !!>>HC 18-3-2021 EPITHELIUM
   if (node(ich)%tipus==2)then                                           !!>>HC 18-3-2021 epithelial internal nodes
      newcell=newcell+1                                                  !!>>HC 18-3-2021 each one belongs to one cell
      cels(newcell)%ctipus=1                                             !!>>HC 18-3-2021 this cell is ctipus 1=epithelial
      cels(newcell)%nunodes=2                                            !!>>HC 18-3-2021 it has 2 nodes
      cels(newcell)%nodela=3                                             !!>>HC 18-3-2021 but the node array has 3 positions
      call random_number(a)                                              !!>>HC 18-3-2021
      cels(newcell)%fase=a                                               !!>>HC 18-3-2021 Give a random number to the cells fase
      if(allocated(cels(newcell)%node)) deallocate(cels(newcell)%node)   !!>>HC 18-3-2021 Allocate the node array inside cells
      allocate (cels(newcell)%node(3))                                   !!>>HC 18-3-2021
      cels(newcell)%node(1)=ich                                          !!>>HC 18-3-2021 The node ich belongs to this cell
      naltrech=node(ich)%altre                                           !!>>HC 18-3-2021 
      cels(newcell)%cex=(node(ich)%x+node(naltrech)%x)/2                 !!>>HC 18-3-2021 Calculate the cell's centroid as the 
      cels(newcell)%cey=(node(ich)%y+node(naltrech)%y)/2                 !!>>HC 18-3-2021 midpoint between nodes
      cels(newcell)%cez=(node(ich)%z+node(naltrech)%z)/2                 !!>>HC 18-3-2021
   elseif (node(ich)%tipus==1)then                                       !!>>HC 18-3-2021 If this is an external epithelial node
          cels(newcell)%node(2)=ich                                      !!>>HC 18-3-2021 it belongs to the previously described cell
          cels(newcell)%node(3)=0                                        !!>>HC 18-3-2021 empty
   else                                                                  !!>>HC 18-3-2021 MESENCHYME
          newcell=newcell+1                                              !!>>HC 18-3-2021 each node has its own cell
          cels(newcell)%ctipus=3                                         !!>>HC 18-3-2021 ctipus for mesenchymal cells
          cels(newcell)%nunodes=1                                        !!>>HC 18-3-2021 mesenchymal cells only have one node
          cels(newcell)%nodela=1                                         !!>>HC 18-3-2021 and the actual size of the node array is 1
          if(allocated(cels(newcell)%node)) deallocate(cels(newcell)%node)  !!>>HC 18-3-2021 Allocate the node array inside cells
          allocate (cels(newcell)%node(1))                               !!>>HC 18-3-2021
          cels(newcell)%node(1)=ich                                      !!>>HC 18-3-2021 Assing cell
          call random_number(a)                                          !!>>HC 18-3-2021
          cels(newcell)%fase=a                                           !!>>HC 18-3-2021 random cell fase
          cels(newcell)%cex=node(ich)%x                                  !!>>HC 18-3-2021 centroid coordinates are the coordinates of the node
          cels(newcell)%cey=node(ich)%y                                  !!>>HC 18-3-2021
          cels(newcell)%cez=node(ich)%z                                  !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!GENES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     !!>>HC 18-3-2021  CREATING THE GENE NETWORK

ng=3                                                                     !!>>HC 18-3-2021
call initiate_gene                                                       !!>>HC 18-3-2021

print *,nparam_per_node,"nparam si"
gen(3)%e=1.0d0
gen(3)%e(nparam_per_node+2)=1.d10

do ich=1,ng                                                              !!>>HC 18-3-2021
   if(allocated(gen(ich)%t)) deallocate(gen(ich)%t)                      !!>>HC 18-3-2021
   allocate(gen(ich)%t(1:ng))                                            !!>>HC 18-3-2021
   gen(ich)%t(:)=0.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%e(:)=0.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%diffu=0.0010d0                                               !!>>HC 18-3-2021
   gen(ich)%mu=0.20d0                                                    !!>>HC 18-3-2021
   gen(ich)%mich=1.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%kindof=4.0d0                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

gen(3)%kindof=4                                                          !!>>HC 18-3-2021
gen(3)%diffu=0.0d0                                                       !!>>HC 18-3-2021
gen(3)%mu=0.0d0                                                          !!>>HC 18-3-2021

if (allocated(gex)) deallocate(gex)                                      !!>>HC 18-3-2021
allocate(gex(1:nda,1:ng))                                                 !!>>HC 18-3-2021
gex(:,:)=0                                                               !!>>HC 18-3-2021
 
if (allocated(kadh)) deallocate(kadh)                                    !!>>HC 18-3-2021
ntipusadh=0                                                              !!>>HC 18-3-2021


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    !!>>HC 18-3-2021
  ! WE CREATE THE INITIAL GRADIENTS !                                    !!>>HC 18-3-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    !!>>HC 18-3-2021
  zmax = maxval(node(:)%z); zmin = minval(node(:)%z)                     !!>>HC 18-3-2021
  xmax = maxval(node(:)%x); xmin = minval(node(:)%x)                     !!>>HC 18-3-2021
  ymax = maxval(node(:)%x); ymin = minval(node(:)%x)                     !!>>HC 18-3-2021

  ! Gene 1 is expressed at the top only.                                 !!>>HC 18-3-2021
  do ich=1, nd                                                           !!>>HC 18-3-2021
    if(node(ich)%z-zmin .ge. 0.9*(zmax-zmin))then                        !!>>HC 18-3-2021
      gex(ich,1)=1                                                       !!>>HC 18-3-2021
      gex(node(ich)%altre,1)=1                                           !!>>HC 18-3-2021
    end if                                                               !!>>HC 18-3-2021
  end do                                                                 !!>>HC 18-3-2021

! Gene 2 is expressed in a triangle on side of the embryo.               !!>>HC 18-3-2021
  xmin1=0.35*(xmax-xmin)+xmin; xmax1=0.65*(xmax-xmin)+xmin               !!>>HC 18-3-2021
  do ich=1, nd                                                           !!>>HC 18-3-2021
    if(node(ich)%tipus > 1)cycle                                         !!>>HC 18-3-2021
    if(node(ich)%x .gt. xmax1 .or. node(ich)%x .lt. xmin1)cycle          !!>>HC 18-3-2021
    if(node(ich)%y-ymin > 0.4*(ymax-ymin))cycle                          !!>>HC 18-3-2021
    gex(ich,2) = 10*(node(ich)%x-xmin1)*(xmax1-node(ich)%x)              !!>>HC 18-3-2021
    gex(node(ich)%altre,2)=gex(ich,2)                                    !!>>HC 18-3-2021
  end do                                                                 !!>>HC 18-3-2021
  gex(:,:)=gex(:,:)/maxval(gex(:,:))                                     !!>>HC 18-3-2021
  
  do ich=1,nd                                                            !!>>HC 18-3-2021
     if (node(ich)%tipus.ge.3)then                                       !!>>HC 18-3-2021
        gex(ich,3)=0.250d0                                               !!>>HC 18-3-2021
     else                                                                !!>>HC 18-3-2021
       gex(ich,3)=1.0d0                                                  !!>>HC 18-3-2021
       if (gex(ich,1)>0.0d0 .and. gex(ich,2)>0.0d0 ) gex(ich,2)=0.0d0    !!>>HC 18-3-2021
       if (gex(ich,1)>0.0d0) gex(ich,1)=1.0d0                            !!>>HC 18-3-2021
     endif                                                               !!>>HC 18-3-2021
  enddo                                                                  !!>>HC 18-3-2021
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine
!*************************************************************************************


!*************************************************************************************

subroutine blastula_icosahedron                                          !!>>HC 18-3-2021 This subroutine creates an epithelial sphere full of mesenchymal cells
implicit none                                                            !!>>HC 18-3-2021 using coordinates already calculated for geodesic icosahedra 
integer ::  ich, jch, ord, ord2, newcell                                 !!>>HC 18-3-2021 stored in a directory called "ico"
real*8 :: radch, reqch, reqch2, dich, addch, addch2                      !!>>HC 18-3-2021 (made by Hugo Cano 30-3-2021)
real*8 :: ratio_add_eqd, ratio_req2_req1                                 !!>>HC 18-3-2021
real*8 :: goldenmean, goldenangle                                        !!>>HC 18-3-2021
real*8, allocatable, dimension(:) :: lonch, latch, xch, ych, zch         !!>>HC 18-3-2021
real*8, allocatable, dimension(:) :: mlonch, mlatch, mxch, mych, mzch    !!>>HC 18-3-2021
integer :: nmesch, ncelepch, nepch, naltrech                             !!>>HC 18-3-2021
real*8 :: xmax, xmax1, ymax, zmax, xmin, ymin, zmin, xmin1               !!>>HC 18-3-2021
integer :: npilech, pilech                                               !!>>HC 19-3-2021
real*8 :: parpil, dismech, dismch                                        !!>>HC 19-3-2021
character*100 :: icfile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!PARAMETERS OF THE BLASTULA!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
radch=1.0d0                                                              !!>>HC 18-3-2021 blastula's radius
reqch=0.12                                                               !!>>HC 18-3-2021 EQD for external epitelial nodes and mesenchyme (type 1 and 3)
ratio_add_eqd=1.20d0                                                     !!>>HC 18-3-2021 ADD/EQD of nodes
ratio_req2_req1=0.90d0                                                   !!>>HC 19-3-2021 Ration EQD external nodes/EQD internal nodes
dich=0.90d0                                                              !!>>HC 18-3-2021 distance between upper and lower epitelial nodes
nmesch=80*int((4.0/3.0)*pi*(radch+reqch)**3)!300 !miguel 2021 automatic  !!>>HC 18-3-2021 number of mesenchymal cells
npilech=3                                                                !!>>HC 19-3-2021 Number of mesenchymal layers
dismch=0.20d0                                                            !!>>HC 19-3-2021 distance between mesenchimal layers
!ncelepch=212;  icfile="../ico/geo_ic_coordinates212.dat"                 !!>>HC 31-3-2021 number of epithelial cells
ncelepch=372;  icfile="../ico/geo_ic_coordinates372.dat"                 !!>>HC 31-3-2021 and path to the file where the icosahedrum coordinates are stored
!ncelepch=482;  icfile="../ico/geo_ic_coordinates482.dat"                 !!>>HC 31-3-2021 use one of these lines, comment the rest
!ncelepch=522;  icfile="../ico/geo_ic_coordinates522.dat"                 !!>>HC 31-3-2021
!ncelepch=632:  icfile="../ico/geo_ic_coordinates632.dat"                 !!>>HC 31-3-2021
!ncelepch=792;  icfile="../ico/geo_ic_coordinates792.dat"                 !!>>HC 31-3-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 31-3-2021
!!!!!!!!!!!DERIVED PARAMETERS!!!!!!!!!!!!!!!!!                           !!>>HC 31-3-2021

addch=ratio_add_eqd*reqch                                                !!>>HC 18-3-2021 ADD for external epitelial nodes and mesenchyme (type 1 and 3)
reqch2=ratio_req2_req1*reqch                                             !!>>HC 19-3-2021 EQD for internal epitelial nodes (type 2)
addch2=ratio_add_eqd*reqch2                                              !!>>HC 19-3-2021 ADD for internal epitelial nodes (type 2)
nepch=ncelepch*2                                                         !!>>HC 18-3-2021 Number of epitelial nodes
ncels=ncelepch+nmesch                                                    !!>>HC 18-3-2021 Total number of cells
ncals=ncels+10                                                           !!>>HC 18-3-2021 number of cells +10
nd=nepch+nmesch                                                          !!>>HC 18-3-2021 total number of nodes
nda=nd+10                                                                !!>>HC 18-3-2021 number of nodes+10
pilech=nmesch/npilech                                                    !!>>HC 19-3-2021 Number of mesenchymal cells per layer
parpil=dich/(real(npilech)+2.0d0)                                        !!>>HC 19-3-2021 Distance between layers (percentage from the center)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!ALLOCATING MEMORY!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021

if(allocated(xch)) deallocate(xch); allocate(xch(1:ncelepch))            !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(ych)) deallocate(ych); allocate(ych(1:ncelepch))            !!>>HC 18-3-2021 cartesian coordinates
if(allocated(zch)) deallocate(zch); allocate(zch(1:ncelepch))            !!>>HC 18-3-2021
if(allocated(lonch)) deallocate(lonch); allocate(lonch(1:ncelepch))      !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(latch)) deallocate(latch); allocate(latch(1:ncelepch))      !!>>HC 18-3-2021 spheric coordinates (latch=latitude lonch=longitude)
lonch=0.0d0; latch=0.0d0                                                 !!>>HC 18-3-2021
xch=0.0d0; ych=0.0d0; zch=0.0d0                                          !!>>HC 18-3-2021

if(allocated(mxch)) deallocate(mxch); allocate(mxch(1:pilech))           !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(mych)) deallocate(mych); allocate(mych(1:pilech))           !!>>HC 18-3-2021 cartesian coordinates
if(allocated(mzch)) deallocate(mzch); allocate(mzch(1:pilech))           !!>>HC 18-3-2021
if(allocated(mlonch)) deallocate(mlonch); allocate(mlonch(1:pilech))     !!>>HC 18-3-2021 These vectors will store the Fibonacci's sphere 
if(allocated(mlatch)) deallocate(mlatch); allocate(mlatch(1:pilech))     !!>>HC 18-3-2021 spheric coordinates (latch=latitude lonch=longitude)
mlonch=0.0d0; mlatch=0.0d0                                               !!>>HC 18-3-2021
mxch=0.0d0; mych=0.0d0; mzch=0.0d0                                       !!>>HC 18-3-2021


if (allocated(node))deallocate(node); allocate(node(1:nda))              !!>>HC 18-3-2021 
if (allocated(cels))deallocate(cels); allocate(cels(1:ncals))            !!>>HC 18-3-2021
if(allocated(nodeo)) deallocate(nodeo); allocate(nodeo(nda))             !!>>HC 18-3-2021
call iniarrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           !!>>HC 18-3-2021
   !******* #2 DEFINING MODEL PARAMETERS *******                         !!>>HC 18-3-2021 Here we set the global parameters so that the program can run

    !implementation and initializations                                  !!>>HC 18-3-2021
    getot=0                                                              !!>>HC 18-3-2021
    itacc=0                                                              !!>>HC 18-3-2021
    nparti=1000                                                          !!>>HC 18-3-2021
    idum=-11111                                                          !!>>HC 18-3-2021
    idumoriginal=idum                                                    !!>>HC 18-3-2021
    nvarglobal_out=5                                                     !!>>HC 18-3-2021

    !physical                                                            !!>>HC 18-3-2021
    temp=0.1d1 !low value is low temperature	                           !!>>HC 18-3-2021
    desmax=0.001                                                         !!>>HC 18-3-2021
    resmax=1d-3                                                          !!>>HC 18-3-2021
    prop_noise=0.10d0                                                    !!>>HC 18-3-2021
    deltamax=1d-2 ! miguel 14-10-13                                      !!>>HC 18-3-2021
    dmax=1                                                               !!>>HC 18-3-2021
    screen_radius=1.0d0                                                  !!>>HC 18-3-2021
    
    maxbox=10.0d0                                                        !!>>HC 19-3-2021
    maxad=5.0d0                                                          !!>>HC 18-3-2021
    maxcycl=1.50d0                                                       !!>>HC 18-3-2021

    ffu=0                                                                !!>>HC 18-3-2021 
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta
    
!!!!!!!!!NODE!!!!!!!!!!!!!!!!!                                           !!>>HC 18-3-2021 HERE WE FILL THE NODE MATRIX
!!!!!NODE TYPES AND ALTRE/CELL ASIGNATION!!!!!!!!!!                      !!>>HC 18-3-2021

newcell=0                                                                !!>>HC 18-3-2021 We decide node type, altres and cells
do ich=1,nd                                                              !!>>HC 18-3-2021 
   if (ich.le.nepch)then                                                 !!>>HC 18-3-2021 This is true for all EPITELIAL nodes
      if (modulo(ich,2).ne.0)then                                        !!>>HC 18-3-2021 This is true for odd nodes (type 2 nodes go always first)
         node(ich)%tipus=2                                               !!>>HC 18-3-2021 So this is an odd epitelial node (1,3,5...) and it must be tipus 2
         node(ich)%altre=ich+1                                           !!>>HC 18-3-2021 Its altre node will be the next node (tipus 1)
         newcell=newcell+1                                               !!>>HC 18-3-2021 Tipus 2 always goes first, so this node starts a new cell
         node(ich)%icel=newcell                                          !!>>HC 18-3-2021 Assign the cell to de node
      else                                                               !!>>HC 18-3-2021
         node(ich)%tipus=1                                               !!>>HC 18-3-2021 This node is an even node (2,4,6..) and it must be tipus 1
         node(ich)%altre=ich-1                                           !!>>HC 18-3-2021 Its altre must be the previous node (tipus 1)
         node(ich)%icel=newcell                                          !!>>HC 18-3-2021 and it is in the same cell than his mate (ie the previous node)
      endif                                                              !!>>HC 18-3-2021
   else                                                                  !!>>HC 18-3-2021 This is for MESENCHIMAL nodes
      newcell=newcell+1                                                  !!>>HC 18-3-2021 Every node belongs to a single cell
      node(ich)%tipus=3                                                  !!>>HC 18-3-2021 The typus of mesenchymal cells
      node(ich)%icel=newcell                                             !!>>HC 18-3-2021 Assing the cell
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

!!!!!NODE PROPERTIES!!!!!!!!!!                                           !!>>HC 18-3-2021
do ich=1,nd                                                              !!>>HC 18-3-2021 Here we give the nodes PHYSICAL PROPERTIES
   if (node(ich)%tipus==1.or.node(ich)%tipus==3)then                     !!>>HC 18-3-2021 External epithelial nodes and mesenchymal nodes have the same size
      node(ich)%eqd=reqch                                                !!>>HC 18-3-2021 
      node(ich)%add=addch                                                !!>>HC 18-3-2021 
   else                                                                  !!>>HC 18-3-2021
      node(ich)%eqd=reqch2                                               !!>>HC 18-3-2021 
      node(ich)%add=addch2                                               !!>>HC 18-3-2021 
   endif                                                                 !!>>HC 18-3-2021
   node(ich)%grd=node(ich)%eqd                                           !!>>HC 18-3-2021
   node(ich)%adh=5.0d0; node(ich)%rec=50.0d0                             !!>>HC 18-3-2021
   node(ich)%erp=20.0d0; node(ich)%est=70.0d0                            !!>>HC 30-3-2021
   node(ich)%eqs=0.150d0; node(ich)%hoo=10.0d0                           !!>>HC 18-3-2021
   node(ich)%mov=0.010d0; node(ich)%dmo=0.0010d0                         !!>>HC 18-3-2021
   node(ich)%pla=0.050d0                                                   !!>>HC 18-3-2021
   node(ich)%cod=0.0d0;  node(ich)%kvol=0.50d0; node(ich)%pld=0.0d0      !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021
  
!!!!!NODE POSITIONS!!!!!!!!!!                                            !!>>HC 30-3-2021 ICOSAHEDRUM COORDINATES  
open(391, file=trim(icfile))                                             !!>>HC 30-3-2021 already calculated and written in the file                   
 do ich=1,ncelepch                                                       !!>>HC 30-3-2021
    read(391,*) xch(ich), ych(ich), zch(ich)                             !!>>HC 30-3-2021
 enddo                                                                   !!>>HC 30-3-2021
 close(391)                                                              !!>>HC 30-3-2021

ord=0                                                                    !!>>HC 18-3-2021 We give the cartesian coordinates to the external nodes
do ich=1,nd                                                              !!>>HC 18-3-2021
   if (node(ich)%tipus==1)then                                           !!>>HC 18-3-2021
      ord=ord+1                                                          !!>>HC 18-3-2021
      node(ich)%x=xch(ord)*radch                                         !!>>HC 18-3-2021 Multiply the coordinate by the radius, otherwise, radius will be =1
      node(ich)%y=ych(ord)*radch                                         !!>>HC 18-3-2021
      node(ich)%z=zch(ord)*radch                                         !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

do ich=1,nd                                                              !!>>HC 18-3-2021 We give the cartesian coordinates to the internal nodes
   if (node(ich)%tipus==2)then                                           !!>>HC 18-3-2021 Basically we multiply the external radius coordinates by the distance between spheres (dich)
      naltrech=node(ich)%altre                                           !!>>HC 18-3-2021
      node(ich)%x=node(naltrech)%x*dich                                  !!>>HC 18-3-2021
      node(ich)%y=node(naltrech)%y*dich                                  !!>>HC 18-3-2021
      node(ich)%z=node(naltrech)%z*dich                                  !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

!do ich=1,pilech                                                          !!>>HC 18-3-2021 MESENCHYMAL NODES
!   mlatch(ich)=asin(-1.0+2.0*real(ich)/(real(pilech)+1))                 !!>>HC 18-3-2021 We create the Fibonacci sphere 
!   mlonch(ich)=goldenangle*real(ich)                                     !!>>HC 18-3-2021 (same as epithelium, but with different node number)
!   mxch(ich)=cos(mlonch(ich))*cos(mlatch(ich))                           !!>>HC 18-3-2021
!   mych(ich)=sin(mlonch(ich))*cos(mlatch(ich))                           !!>>HC 18-3-2021
!   mzch(ich)=sin(mlatch(ich))                                            !!>>HC 18-3-2021
      
!enddo

!ord=nepch                                                                !!>>HC 19-3-2021
!do ich=1,npilech                                                         !!>>HC 19-3-2021 We assign the coordinates to each spherical layer of mesenchymal cells
!   dismech=(radch*dich)-(dismch*ich)                                     !!>>HC 19-3-2021 This is the radius of each mesenchymal spherical layer
!   ord2=0                                                                !!>>HC 19-3-2021
!   do jch=1,pilech                                                       !!>>HC 19-3-2021
!      ord=ord+1                                                          !!>>HC 19-3-2021
!      node(ord)%x=mxch(jch)*dismech                                      !!>>HC 19-3-2021
!      node(ord)%y=mych(jch)*dismech                                      !!>>HC 19-3-2021
!      node(ord)%z=mzch(jch)*dismech                                      !!>>HC 19-3-2021
!   enddo                                                                 !!>>HC 19-3-2021
!enddo                                                                    !!>>HC 19-3-2021

do ich=nepch+1,nd                                 ! Miguel 2021 randomized mesenchyme
   7757 call random_number(x) ; call random_number(y) ; call random_number(z)
   xx=radch-reqch ; yy=(0.9*xx)**2
   x=x*2.0*xx-xx   ; y=y*2.0*xx-xx   ; z=z*2.0*xx-xx
   if((x**2+y**2+z**2).lt.yy)then 
     node(ich)%x=x ;  node(ich)%y=y ;  node(ich)%z=z
   else ; goto 7757 ; end if
end do

nodeo=node                                                               !!>>HC 19-3-2021 set the original node properties


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  !!>>HC 18-3-2021
!!!!!!!!!!!!CELLS!!!!!!!!!!!!!!!!!!!!!!                                  !!>>HC 18-3-2021 FILL THE CELS MATRIX

newcell=0                                                                !!>>HC 18-3-2021
do ich=1,nd                                                              !!>>HC 18-3-2021 EPITHELIUM
   if (node(ich)%tipus==2)then                                           !!>>HC 18-3-2021 epithelial internal nodes
      newcell=newcell+1                                                  !!>>HC 18-3-2021 each one belongs to one cell
      cels(newcell)%ctipus=1                                             !!>>HC 18-3-2021 this cell is ctipus 1=epithelial
      cels(newcell)%nunodes=2                                            !!>>HC 18-3-2021 it has 2 nodes
      cels(newcell)%nodela=3                                             !!>>HC 18-3-2021 but the node array has 3 positions
      call random_number(a)                                              !!>>HC 18-3-2021
      cels(newcell)%fase=a                                               !!>>HC 18-3-2021 Give a random number to the cells fase
      if(allocated(cels(newcell)%node)) deallocate(cels(newcell)%node)   !!>>HC 18-3-2021 Allocate the node array inside cells
      allocate (cels(newcell)%node(3))                                   !!>>HC 18-3-2021
      cels(newcell)%node(1)=ich                                          !!>>HC 18-3-2021 The node ich belongs to this cell
      naltrech=node(ich)%altre                                           !!>>HC 18-3-2021 
      cels(newcell)%cex=(node(ich)%x+node(naltrech)%x)/2                 !!>>HC 18-3-2021 Calculate the cell's centroid as the 
      cels(newcell)%cey=(node(ich)%y+node(naltrech)%y)/2                 !!>>HC 18-3-2021 midpoint between nodes
      cels(newcell)%cez=(node(ich)%z+node(naltrech)%z)/2                 !!>>HC 18-3-2021
   elseif (node(ich)%tipus==1)then                                       !!>>HC 18-3-2021 If this is an external epithelial node
          cels(newcell)%node(2)=ich                                      !!>>HC 18-3-2021 it belongs to the previously described cell
          cels(newcell)%node(3)=0                                        !!>>HC 18-3-2021 empty
   else                                                                  !!>>HC 18-3-2021 MESENCHYME
          newcell=newcell+1                                              !!>>HC 18-3-2021 each node has its own cell
          cels(newcell)%ctipus=3                                         !!>>HC 18-3-2021 ctipus for mesenchymal cells
          cels(newcell)%nunodes=1                                        !!>>HC 18-3-2021 mesenchymal cells only have one node
          cels(newcell)%nodela=1                                         !!>>HC 18-3-2021 and the actual size of the node array is 1
          if(allocated(cels(newcell)%node)) deallocate(cels(newcell)%node)  !!>>HC 18-3-2021 Allocate the node array inside cells
          allocate (cels(newcell)%node(1))                               !!>>HC 18-3-2021
          cels(newcell)%node(1)=ich                                      !!>>HC 18-3-2021 Assing cell
          call random_number(a)                                          !!>>HC 18-3-2021
          cels(newcell)%fase=a                                           !!>>HC 18-3-2021 random cell fase
          cels(newcell)%cex=node(ich)%x                                  !!>>HC 18-3-2021 centroid coordinates are the coordinates of the node
          cels(newcell)%cey=node(ich)%y                                  !!>>HC 18-3-2021
          cels(newcell)%cez=node(ich)%z                                  !!>>HC 18-3-2021
   endif                                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     !!>>HC 18-3-2021
!!!!!!!!!!!!!!!!GENES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     !!>>HC 18-3-2021  CREATING THE GENE NETWORK

ng=3                                                                     !!>>HC 18-3-2021
call initiate_gene                                                       !!>>HC 18-3-2021

do ich=1,ng                                                              !!>>HC 18-3-2021 This sets up gene properties
   if(allocated(gen(ich)%t)) deallocate(gen(ich)%t)                      !!>>HC 18-3-2021
   allocate(gen(ich)%t(1:ng))                                            !!>>HC 18-3-2021
   gen(ich)%t(:)=0.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%e(:)=0.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%diffu=0.000d0                                                !!>>HC 18-3-2021
   gen(ich)%mu=0.0d0                                                     !!>>HC 18-3-2021
   gen(ich)%mich=1.0d0                                                   !!>>HC 18-3-2021
   gen(ich)%kindof=4.0d0                                                 !!>>HC 18-3-2021
enddo                                                                    !!>>HC 18-3-2021

gen(3)%kindof=4                                                          !!>>HC 18-3-2021
gen(3)%diffu=0.0d0                                                       !!>>HC 18-3-2021
gen(3)%mu=0.0d0                                                          !!>>HC 18-3-2021

if (allocated(gex)) deallocate(gex)                                      !!>>HC 18-3-2021 allocating matrix of gene expression
allocate(gex(1:nda,1:ng))                                                 !!>>HC 18-3-2021
gex(:,:)=0                                                               !!>>HC 18-3-2021
 
if (allocated(kadh)) deallocate(kadh)                                    !!>>HC 18-3-2021 allocating matrix of specific adhesion
ntipusadh=0                                                              !!>>HC 18-3-2021


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    !!>>HC 18-3-2021
  ! WE CREATE THE INITIAL GRADIENTS !                                    !!>>HC 18-3-2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    !!>>HC 18-3-2021 We create two gradients
  zmax = maxval(node(:)%z); zmin = minval(node(:)%z)                     !!>>HC 18-3-2021 
  xmax = maxval(node(:)%x); xmin = minval(node(:)%x)                     !!>>HC 18-3-2021 
  ymax = maxval(node(:)%x); ymin = minval(node(:)%x)                     !!>>HC 18-3-2021 

  ! Gene 1 is expressed at the top only.                                 !!>>HC 18-3-2021 Gene 1 is expressed in the animal pole
  do ich=1, nd                                                           !!>>HC 18-3-2021
    if(node(ich)%z-zmin .ge. 0.9*(zmax-zmin))then                        !!>>HC 18-3-2021
      gex(ich,1)=1.0d0                                                   !!>>HC 18-3-2021
      gex(node(ich)%altre,1)=1                                           !!>>HC 18-3-2021
    end if                                                               !!>>HC 18-3-2021
  end do                                                                 !!>>HC 18-3-2021

! Gene 2 is expressed in a triangle on side of the embryo.               !!>>HC 18-3-2021 Gene 2 is expressed in the dorsal area
  xmin1=0.35*(xmax-xmin)+xmin; xmax1=0.65*(xmax-xmin)+xmin               !!>>HC 18-3-2021
  do ich=1, nd                                                           !!>>HC 18-3-2021
    if(node(ich)%tipus > 1)cycle                                         !!>>HC 18-3-2021
    if(node(ich)%x .gt. xmax1 .or. node(ich)%x .lt. xmin1)cycle          !!>>HC 18-3-2021
    if(node(ich)%y-ymin > 0.4*(ymax-ymin))cycle                          !!>>HC 18-3-2021
    if(gex(ich,1)>0.0d0)cycle
    gex(ich,2) = 1.0d0 !10*(node(ich)%x-xmin1)*(xmax1-node(ich)%x)              !!>>HC 18-3-2021
    gex(node(ich)%altre,2)=gex(ich,2)                                    !!>>HC 18-3-2021
  end do                                                                 !!>>HC 18-3-2021
  gex(:,:)=gex(:,:)/maxval(gex(:,:))                                     !!>>HC 18-3-2021
  
  do ich=1,nd                                                            !!>>HC 18-3-2021 Gene 3 is expressed in all the embryo 
     if (node(ich)%tipus.ge.3)then                                       !!>>HC 18-3-2021 0.5 in mesenchymal nodes
        gex(ich,3)=0.50d0                                                !!>>HC 23-3-2021 1 in epithelial nodes
     else                                                                !!>>HC 18-3-2021 (it is the one that will activate division)
       gex(ich,3)=1.0d0                                                  !!>>HC 18-3-2021
       if (gex(ich,1)>0.0d0 .and. gex(ich,2)>0.0d0 ) gex(ich,2)=0.0d0    !!>>HC 18-3-2021
       if (gex(ich,1)>0.0d0) gex(ich,1)=1.0d0                            !!>>HC 18-3-2021
     endif                                                               !!>>HC 18-3-2021
  enddo                                                                  !!>>HC 18-3-2021
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end subroutine

!*************************************************************************************

subroutine epi_mes

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=10     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-1.0d0  !z-position of uppermost layer
        prec=0.01    !accuracy for the optional adaptive stepsize

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

    !print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(7)=0
    ffu(8)=0
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
    
    
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=1d1    	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0;	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.50; 
        node(i)%hoo=5d1
        node(i)%erp=5d1
        !!node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%pla=1d0
        node(i)%kvol=1d1
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.50 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !!node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    !node(4)%altre=13 ; node(4)%icel=7 ; node(13)%altre=4
    !do i=5,nd
    !  if(node(i)%tipus==1)then
    !    node(i)%altre=i-3 ; node(i-3)%altre=i
    !    node(i)%icel=node(i-3)%icel
    !  end if
    !end do
      
    
    
    
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

       !****a rellenar****


    !Gene-behavior interactions

       !****a rellenar****

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

         !****a rellenar****

    call update_npag
   
    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3
 
    node(:)%dif=0.0d0

end subroutine

!**************************************************************************************************************
subroutine epi_growth_and_division(re,rm)

integer re,rm

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=re       !number of radial layers of nodes per cell
	radicel=rm    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	packed=1
    mradi=0     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-1   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters


!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    
    ffu(7)=0
    ffu(8)=0 !neighboring by triangulation
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.30  
        node(i)%hoo=1d1
        node(i)%erp=3d1
        !!node(i)%est=3d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d0;node(i)%adh=1d0 	!>>Miquel 26-10-12
        node(i)%rep=1d0;node(i)%rec=1d0	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ;node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.3
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=3d1  !only for epithelium
        !!node(i)%est=3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=cels(i)%nunodes*4 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=cels(i)%nunodes*4 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0.5d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.1d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+1)=1d0
      gen(1)%e(nparam_per_node+2)=1d-1
      cels(:)%maxsize_for_div=28

    !Gene-gene interactions
      gen(1)%t(1)=1.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      do i=1,nd
        gex(i,1)=1.0d0
      end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************************************

subroutine epi_mes_growth_and_division

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=14     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(5)=0
    ffu(12)=1
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.30  
        node(i)%hoo=1d1
        node(i)%erp=1d0
        !!node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.2  
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !!node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=20 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=20 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0.5d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+1)=1d-2
      gen(1)%e(nparam_per_node+2)=1d-2

    !Gene-gene interactions
      gen(1)%t(1)=1.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      do i=1,nd
        gex(i,1)=1d0
      end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine


!**************************************************************************************************************

subroutine mes_polar_growth

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
        !packed=1
	mradi=8     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 ;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters



    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-4 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    deltamin=1d-3
    ecmmax=0.25d0
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=0
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25 
        node(i)%hoo=1d1
        node(i)%erp=5d-1
        !!node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    ! Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.40; 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=mradi*2 !>>> Is 5-2-14
      end do
    end if


    do i=1,ncels
      call random_number(a)
      cels(i)%fase=a
    end do



   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters


    !Gene-behavior interactions
    gen(2)%e(nparam_per_node+1)=5.0d-1 
    gen(2)%e(nparam_per_node+2)=2.0d-2  
    gen(1)%e(nparam_per_node+8)=1.0d0  
    gen(1)%e(nparam_per_node+11)=0.0d0
    gen(2)%e(nparam_per_node+9)=1.0d0

    gen(:)%diffu=0.1d0

    !min_comp=-1d0

    !Adhesion molecules

    ntipusadh=5
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecule interactions
!        kadh(1,1)=1d1
!        kadh(1,2)=1d1
!        kadh(1,3)=1d1
!        kadh(2,1)=1d1
        kadh(2,2)=1d0
        kadh(2,3)=1d0
!        kadh(3,1)=1d1
        kadh(3,2)=1d0
        kadh(3,3)=1d0

    end if

    gen(1)%e(1)=1
    gen(2)%e(1)=2
    gen(3)%e(1)=3
    gen(4)%e(1)=4
    gen(5)%e(1)=5

    !Gene expression on nodes
      do i=1,nd
        if(node(i)%tipus==3)then
          gex(i,1)=1.0d-2
          gex(i,2)=1d0
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine 

!*****************************************************************************************

subroutine mes_polar_growth_old



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=10     !number of nodes per cell
	mradicel=1   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(nodecel<mradi)then;nodecel=mradi;nodecela=2*nodecel+1;end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters



    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3  
    prop_noise=0.1d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d4;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=1d4;node(i)%rec=1d6	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.3 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2    
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

  !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13 diffu0.5d3 en epi
      gen(2)%diffu=1d1 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions

    gen(2)%e(nparam_per_node+1)=2.0d-5  
    gen(2)%e(nparam_per_node+2)=2.0d-4  
    gen(1)%e(nparam_per_node+8)=1.0d0  
!    gen(1)%e(nparam_per_node+11)=0.0d0
    gen(2)%e(nparam_per_node+9)=0.5d0


    !Gene-gene interactions
      gen(2)%t(2)=1.0d0

    !Adhesion molecules

	ntipusadh=0
      if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)
      gex(:,2)=1d0

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!**************************************************************************************************************

subroutine epi_growth




!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=2       !number of radial layers of nodes per cell
    radicel=1    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=0 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.70; 
        node(i)%hoo=1d1
        node(i)%erp=5d1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
    !  j=0
    !  do i=1,radicel-2
    !    j=j+i
    !  end do
    !  j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    !  do i=j+1,ndepi
    !    node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !!node(i)%est=1d1
    !  end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0.5d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.1d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+1)=1d-4

    !Gene-gene interactions
      !gen(1)%t(2)=1.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      do i=1,nd
        gex(i,1)=1.0d0
      end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine 
!*************************************************************************************************************************************

subroutine invagination  !Changed some parameters, now there is only contraction on the apical side (no expansion on the basal).
                         !Instead volume conservation and plasticity do the trick. But in order to work properly the apical contracting
                         !nodes have to have %pla=0 and %kvol=0 (only those)     !>>Miquel28-7-14


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=6    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !!nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta
    maxad=5.0d0 !!>> HC 13-7-2021 
     !>>> TT
    
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=5.0d0  ;node(i)%adh=5.0d0         	!>>Miquel 26-10-12
        node(i)%rep=5d1    ;node(i)%rec=5d1	        !>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0     	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.60  
        node(i)%hoo=1d1
        node(i)%erp=20d0
        node(i)%est=60d0
        node(i)%pla=2d0
        node(i)%kvol=5d-1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0

   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0.0d0 ; gen(1)%kindof=1 ; gen(1)%mu=0d0 ;! gen(1)%label="basal expressed gene"
      gen(2)%diffu=0.0d0 ; gen(2)%kindof=1 ; gen(2)%mu=0d0 ;! gen(2)%label="apically expressed gene, contraction"

    !Gene-behavior interactions
      a=node(1)%add/node(1)%eqd
!      gen(1)%e(5)=0.1d0  !that is the maximal req possible
!      gen(2)%e(5)=-0.05d0

      !gen(1)%e(5)=9.0d0   ! >>> Is 28-6-14 !that is the maximal req possible

      gen(2)%e(21)=-0.10d0 ! >>> Is 28-6-14
      gen(2)%e(6)=0.10d0
      !gen(2)%e(27)=-1d-1  !this affects plasticity and
      !gen(2)%e(28)=-1d2  !volume conservation of the contracting nodes
      gen(1)%e(21)=0.10d0 ! >>> Is 28-6-14
      gen(1)%e(6)=0.15d0


!      gen(1)%e(7)=gen(1)%e(5)*1.1  !that is the maximal da possible
!      gen(2)%e(7)=gen(2)%e(5)*0.1


    !Gene-gene interactions

      gen(1)%t=0d0
      gen(2)%t=0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

      !we get the cell in the middle to express different genes and these are different in the upper and lower part of the cell
      gex=0.0d0
      do i=1,nd
        j=node(i)%icel
        if (j<8) then
          node(i)%pla=0d0
          node(i)%kvol=0d0
          if (node(i)%tipus==1) then
            gex(i,1)=1.0d0
          else
            gex(i,2)=1.0d0
          end if
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine
!*************************************************************************************************************************************

subroutine polarized



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=3    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d1 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.01d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2 


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1  ;node(i)%adh=1.0d1         	!>>Miquel 26-10-12
        node(i)%rep=1d1    ;node(i)%rec=1d1	        !>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0     	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.20  
        node(i)%hoo=1d1
        node(i)%erp=1d-1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=0 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0

    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+8)=1d0

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!*************************************************************************

subroutine mesenchyme

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=8     !number of nodes per cell
	mradicel=3   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3   
    prop_noise=0.1d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=1d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

!******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters


    !Gene-behavior interactions


    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

         !****a rellenar****

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine mesenchyme

!**********************************************************************************************************

subroutine mes_ecm_degradation

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=8      !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=0d0     !z-position of uppermost layer

    !ECM dimension parameters
	xradi=8  !number of radial node layers
    xlayer=5   !number of ECM node layers
	zx=-0.5d0   !z-position of uppermost layer


!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if

    if(xradi>0.and.xlayer>0)then
      j=0
      do i=1,xradi-1
        j=j+i
      end do
      ndx=(6*j+1)*xlayer	!number of ECM nodes
    else
      ndx=0
    end if
    !End initializing dimension parameters

    nd=ndepi+ndmes+ndx
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10

!  print*,"nd",nd,"ndx",ndx,"ndmes",ndmes
!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3   
    prop_noise=0.0d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=0 !
    ffu(5)=0 !>>> Is 5-2-14
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.40 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!ECM
    if(xradi>0.and.xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=1d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.20 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
    if(xradi>0.and.xlayer>0)               call matrix(xradi,xlayer,zx)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

!******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=3
    call initiate_gene
    !Gene parameters
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))

      gen(1)%diffu=1.0d1 ; gen(1)%kindof=2 ; gen(1)%mu=0d0 !primary form ; adhesion protein for mesenchyme
      gen(1)%post(1)=2 
      gen(2)%diffu=5.0d-1 ; gen(2)%kindof=4 ; gen(2)%mu=1d0 !extracellular protease
      gen(2)%pre(1)=1
      gen(3)%diffu=0.0d0 ; gen(3)%kindof=1 ; gen(3)%mu=0d0 !ECM adhesion protein


    !Gene-behavior interactions
      gen(1)%e(1)=1
      gen(3)%e(1)=2
      gen(2)%e(nparam_per_node+14)=5d-3

    !Gene-gene interactions
      gen(1)%t(1)=1.0d1

      gen(1)%nww=1
      gen(1)%r(1,1)=1
      gen(1)%r(1,2)=2
      gen(1)%r(1,3)=1.0d0

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         kadh(1,1)=1d0
         kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
         kadh(2,2)=1d0

    end if

    !Gene expression on nodes
      gex(1:ndmes,1)=1d0
      gex(ndmes+1:nd,3)=1d0

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************************************************************

subroutine mes_cell_sorting
!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=8     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=3      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameteers

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=5d0 !low value is low temperature	
    desmax=5d-2
    resmax=1d-3  
    prop_noise=0.5
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=0.7d0
    deltamin=1d-3

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(7)=1
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise

    ffu(5)=0 !integration method 0=euler , 1=runge-kutta forces , 2=runge-kutta forces+genes
    ffu(12)=0 !adaptive time step 0=no , 1=yes for forces , 2=yes for forces+genes
    ffu(20)=0      !!>> HC 16-12-2020
    maxad=50.0d0  !!>> HC 16-12-2020

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 ;node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=0d0 
        node(i)%hoo=0d0
        node(i)%erp=1d0
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=6d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=6d1;node(i)%rec=8d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ;node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*2.0d0 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=5d0
        node(i)%dmo=5d-2
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)!mesenq_cell_sorting(mradi,ncels)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

  !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0.0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0d0 ! miguel 14-10-13
      gen(2)%diffu=0.0d0 ; gen(2)%kindof=0 ; gen(2)%mu=0d0 ! miguel 14-10-13


    !Gene-behavior interactions
       gen(1)%e(1)=1d0
       gen(2)%e(1)=2d0


    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

        kadh(1,1)=1d2            ! "type 2" cells surround "type 1" cells
        kadh(1,2)=1d1            !
        kadh(2,1)=kadh(1,2)      !
        kadh(2,2)=1d0		     !

    end if

    !Gene expression on nodes

      do i=1,ncels!ndepi+1,nd      ! here we put the differential adhesion (cells to type 1 or type 2 at random) 
            call random_number(a)        
           if(a>0.7d0)then
            do j=1,cels(i)%nunodes 
              k=cels(i)%node(j)
              gex(k,1)=1d0            ! type 1 cells express gene 1   
              !gex(i,2)=0.1d0  
            end do                
           else            
            do j=1,cels(i)%nunodes
              k=cels(i)%node(j)
              !gex(i,1)=0.1d0
              gex(k,2)=1d0            ! type 2 cells express gene 2 
            end do
           endif                   
       end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine mes_cell_sorting


!***********************************************************************************************************************

subroutine teloblast ! miguel4-11-13 (subtitution of the previous mes_shell)

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=100     !number of nodes per cell
	mradicel=1   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=1d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-19989
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-1 !low value is low temperature	
    desmax=0.001
    resmax=0.001   
    prop_noise=0.0d0
    reqmin=0.05d0 
    deltamax=2d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    ffu(4)=0
    ffu(5)=0 !>>> Is 5-2-14
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=1d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.60d0 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      !write(*,*)'node',i,'xyx',node(i)%x,node(i)%y,node(i)%z
      !node(i)%z=abs(node(i)%z)
    !write(*,*)'ndz',node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* eggshell initialization 
   ! call shellinicial
   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=10 !cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

!******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters

    
    do i=ndepi+1,nd   
      if(node(i)%z.le.0.6d0)then !parametrizar
        gex(i,2)=1d0
      endif
    end do

    !Gene-behavior interactions
    gen(1)%diffu=0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13
    gen(2)%diffu=0d0 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0

    gen(2)%e(nparam_per_node+1)=1d-3 !  miguel 14-10-13 (growth)
    gen(2)%e(nparam_per_node+2)=1d-2 !  miguel 14-10-13 (cell cycle)    

    gen(1)%e(nparam_per_node+12)=1d0 ! div asimétrica
    gen(1)%e(nparam_per_node+11)=1d0 ! dependencia polarizacion quimica
    gen(1)%e(nparam_per_node+8)=1d0  ! polarizacion quimica (grad Z)
    !gen(3)%e(nparam_per_node+8)=1d0  ! polarizacion quimica (cell-contact)
    
    !Gene-behavior interactions
    !gen(2)%e(1)=1

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

         !****a rellenar****

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine teloblast

!************************************************************

subroutine mesenchyme_apop

call mes_ecm

    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=1

goto 24

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=15     !number of nodes per cell
	mradicel=3   !number of radial cell layers
	layer=3      !number of planar cell layers
	zmes=0.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.001
    resmax=1d-3   
    prop_noise=0.0d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(7)=1
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=1d0 	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.30d0 !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.30 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

  !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
       gen(1)%diffu=0.0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0d0 ! miguel 14-10-13

    !Gene-behavior interactions
24     do i=1,ng
         gen(i)%e=0.0d0
       end do
       gen(1)%e(nparam_per_node+3)=5d-1   !apoptosis

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

      do i=1,cels(1)%nunodes ; gex(cels(1)%node(i),1)=1.0d0 ; end do
      do i=1,cels(2)%nunodes ; gex(cels(2)%node(i),1)=1.0d0  ; end do
      do i=1,cels(4)%nunodes ; gex(cels(4)%node(i),1)=1.0d0  ; end do
      do i=1,cels(6)%nunodes ; gex(cels(6)%node(i),1)=1.0d0  ; end do

      do i=1,cels(19+1)%nunodes ; gex(cels(19+1)%node(i),1)=1.0d0 ; end do
      do i=1,cels(19+2)%nunodes ; gex(cels(19+2)%node(i),1)=1.0d0  ; end do
      do i=1,cels(19+4)%nunodes ; gex(cels(19+4)%node(i),1)=1.0d0  ; end do
      !do i=1,cels(19+6)%nunodes ; gex(cels(19+6)%node(i),1)=1.0d0  ; end do

      !do i=1,cels(19+19+1)%nunodes ; gex(cels(19+19+1)%node(i),1)=1.0d0 ; end do
      !do i=1,cels(19+19+2)%nunodes ; gex(cels(19+19+2)%node(i),1)=1.0d0  ; end do
      !do i=1,cels(19+19+4)%nunodes ; gex(cels(19+19+4)%node(i),1)=1.0d0  ; end do
      !do i=1,cels(19+19+6)%nunodes ; gex(cels(19+19+6)%node(i),1)=1.0d0  ; end do





!    do i=1,cels(57)%nunodes ; gex(cels(57)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(56)%nunodes ; gex(cels(56)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(54)%nunodes ; gex(cels(54)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(45)%nunodes ; gex(cels(45)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(44)%nunodes ; gex(cels(44)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(39)%nunodes ; gex(cels(39)%node(i),1)=1.0d0  ; end do

!    do i=1,cels(38)%nunodes ; gex(cels(38)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(37)%nunodes ; gex(cels(37)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(35)%nunodes ; gex(cels(35)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(25)%nunodes ; gex(cels(25)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(26)%nunodes ; gex(cels(26)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(20)%nunodes ; gex(cels(20)%node(i),1)=1.0d0  ; end do

!    do i=1,cels(94)%nunodes ; gex(cels(94)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(95)%nunodes ; gex(cels(95)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(92)%nunodes ; gex(cels(92)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(83)%nunodes ; gex(cels(83)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(82)%nunodes ; gex(cels(82)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(77)%nunodes ; gex(cels(77)%node(i),1)=1.0d0  ; end do

!    do i=1,cels(113)%nunodes ; gex(cels(113)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(114)%nunodes ; gex(cels(114)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(111)%nunodes ; gex(cels(111)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(102)%nunodes ; gex(cels(102)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(101)%nunodes ; gex(cels(101)%node(i),1)=1.0d0  ; end do
!    do i=1,cels(96)%nunodes  ; gex(cels(96)%node(i),1)=1.0d0  ; end do


    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************************************************************

subroutine epiteli(radi,radicel,zepi,dreq,dreqs)		!zepi is the z-position of the bottom layer of nodes
integer            ::radi,radicel,valcel
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi,dreq,dreqs

    di=dreqs*2
    de=dreq*2  !miguel 14-10-13
    vector=0.0d0
    
    node(:)%marge=1
 

	beta=2d0*pi/6d0


	do i=1,ncelsepi
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


	ii=0
	do i=1,2
		ii=ii+1
!		alt=real(i-1)+zepi
        if(i==1)then
          alt=zepi
        else
          alt=zepi+di
        end if
!		print*,"alt",alt
		node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
        if(i==1) node(ii)%marge=0
		if(i==1)then
			node(ii)%tipus=2
		else if(i==2)then
			node(ii)%tipus=1
		end if
		cels(1)%nunodes=cels(1)%nunodes+1
		cels(1)%node(cels(1)%nunodes)=ii
!                node(ii)%marge=1
		node(ii)%icel=1
		do j=2,radi
			angle=beta
!			print*,"angle",angle
			d=de*(j-1d0)
			ii=ii+1
			node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
			if(i==1)then
				node(ii)%tipus=2
			else
				node(ii)%tipus=1
			end if
			cels(1)%nunodes=cels(1)%nunodes+1
			cels(1)%node(cels(1)%nunodes)=ii
			node(ii)%icel=1
			p1(1)=d;p1(2)=0;p1(3)=alt
			p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
!				print*,"jj",jj
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					if(i==1)then
						node(ii)%tipus=2
					else
						node(ii)%tipus=1
					end if
					cels(1)%nunodes=cels(1)%nunodes+1
					cels(1)%node(cels(1)%nunodes)=ii
					node(ii)%icel=1
				end do
			end if
			do k=2,5
				angle=angle+beta
				p1=p2
				ii=ii+1
				node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
				if(i==1)then
					node(ii)%tipus=2
				else
					node(ii)%tipus=1
				end if
				cels(1)%nunodes=cels(1)%nunodes+1
				cels(1)%node(cels(1)%nunodes)=ii
				node(ii)%icel=1
				p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
				jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
				if(jj>0)then
					vector=p2-p1
					modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
					iii=jj+1
					vector=vector/iii
					do kk=1,jj
						ii=ii+1
						node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
						if(i==1)then
							node(ii)%tipus=2
						else
							node(ii)%tipus=1
						end if
						cels(1)%nunodes=cels(1)%nunodes+1
						cels(1)%node(cels(1)%nunodes)=ii
						node(ii)%icel=1
					end do
				end if
			end do
			angle=angle+beta
!			print*,"angle",angle
			p1=p2
			ii=ii+1
			node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
			if(i==1)then
				node(ii)%tipus=2
			else
				node(ii)%tipus=1
			end if
			cels(1)%nunodes=cels(1)%nunodes+1
			cels(1)%node(cels(1)%nunodes)=ii
			node(ii)%icel=1
			p2(1)=d;p2(2)=0;p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					if(i==1)then
						node(ii)%tipus=2
					else
						node(ii)%tipus=1
					end if
					cels(1)%nunodes=cels(1)%nunodes+1
					cels(1)%node(cels(1)%nunodes)=ii
					node(ii)%icel=1	
				end do
			end if
		end do
	end do

!hem fet la 1a cel al centre 

!ara farem la resta copiant i pegant



if(radicel>1)then
	a=2*sqrt(((radi-1)*de)**2-(0.5*de*(radi-1))**2)+de	!distancia de centre de cel·lula a centre de cel·lula
	valcel=1
	do i=2,radicel
!		angle=0.5*pi		!"ok" hexagonal configuration
		angle=0.5*pi+pi/12	!perfect hexagonal configuration
!		print*,"i",i,"angle",angle
!		vector(1)=0;vector(2)=a*(i-1);vector(3)=0			!"ok" hexagonal configuration
		vector(1)=-a*(i-1)*sin(pi/12);vector(2)=a*(i-1)*cos(pi/12);vector(3)=0	!perfect hexagonal configuration
		valcel=valcel+1
		do j=1,cels(1)%nunodes
			ii=ii+1
            if (j==1) node(ii)%marge=0
			node(ii)%x=node(j)%x+vector(1);node(ii)%y=node(j)%y+vector(2);node(ii)%z=node(j)%z
			node(ii)%tipus=node(j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do

		p1(1)=a*(i-1)*dcos(angle-beta);p1(2)=a*(i-1)*dsin(angle-beta);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
			iii=jj+1
			vector=vector/iii
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (j==1) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
		p1=p2
!if(1==2)then
		do k=2,5
			angle=angle+beta
!			print*,"angle",angle
			p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
			valcel=valcel+1
			do j=1,cels(1)%nunodes
				ii=ii+1
                if (j==1) node(ii)%marge=0
				node(ii)%x=node(j)%x+p2(1);node(ii)%y=node(j)%y+p2(2);node(ii)%z=node(j)%z
				node(ii)%tipus=node(j)%tipus
				cels(valcel)%nunodes=cels(valcel)%nunodes+1
				cels(valcel)%node(cels(valcel)%nunodes)=ii
				node(ii)%icel=valcel
			end do
			jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!				iii=jj+1
				vector=vector/(jj+1)
				do kk=1,jj
					valcel=valcel+1
					do j=1,cels(1)%nunodes
						kkk=cels(valcel-1)%node(j)
						ii=ii+1
                                                if (j==1) node(ii)%marge=0
						node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
						node(ii)%tipus=node(j)%tipus
						cels(valcel)%nunodes=cels(valcel)%nunodes+1
						cels(valcel)%node(cels(valcel)%nunodes)=ii
						node(ii)%icel=valcel
					end do
				end do
			end if
			p1=p2
		end do
		angle=angle+beta
!		p2(1)=0;p2(2)=a*(i-1);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		valcel=valcel+1
		do j=1,cels(1)%nunodes
			ii=ii+1
            if (j==1) node(ii)%marge=0
			node(ii)%x=node(j)%x+p2(1);node(ii)%y=node(j)%y+p2(2);node(ii)%z=node(j)%z
			node(ii)%tipus=node(j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!			iii=jj+1
			vector=vector/(jj+1)
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                                        if (j==1) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
!end if
	end do
end if




	!define the cell's centroid

	do i=1,ncelsepi
		cels(i)%ctipus=1
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
		do j=1,cels(i)%nunodes
			k=cels(i)%node(j)
            if(node(k)%tipus==1)then
	 			cels(i)%cex=cels(i)%cex+node(k)%x
 				cels(i)%cey=cels(i)%cey+node(k)%y
				cels(i)%cez=cels(i)%cez+node(k)%z
			end if
		end do
		cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
		cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
		cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
	end do


	do i=1,ndepi		!veins i parametres mecanics
		ii=node(i)%icel
		if(node(i)%tipus==1)then
   		  node(i)%altre=i-cels(1)%nunodes/2
		else
		  node(i)%altre=i+cels(1)%nunodes/2
		end if
	end do

!        node(1)%marge=0.0d0

        !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster >>> Is 14-9-13 
	do i=1,ncelsepi
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (node(k)%marge==0) then
              ii=k
              exit
            end if
	  end do										
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (k==ii) then
              jj=j
              exit
            end if
	  end do									
          cels(i)%node(jj)=jjj
	end do											


    !some funny rotation to make a less biased diffusion lost
!    do i=1,nd
!      c=0.75
!      a=node(i)%x*cos(c)-node(i)%y*sin(c)
!      b=node(i)%x*sin(c)+node(i)%y*cos(c)
!      node(i)%x=a
!      node(i)%y=b
!    end do

    node(:)%talone=0.0d0

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine mesenq(radi,radicel,layer,zmes,dreq)									!!! miguel 20.11
integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,ent,ont       ! number of layers and concentric hexagons 
real*8   ::rad,der,zmes,de,di,dreq
real*8   :: xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
de=dreq*2                    ! call radius
!di=2.0d0*de                 ! distance between cells and layers (it has to be >2 to avoid cell contacts)
di=2*de+2*de*cos(60d0*rad)

	do i=ncelsepi+1,ncels               
!		cels(i)%nunodes=radi+1
		cels(i)%nunodes=radi				!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
!		allocate(cels(i)%node(radi+1))
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%node=0
    end do

node(ndepi+1:)%marge=1

!	print*,"node",size(node),"cels",size(cels),"celsnode",size(cels(1)%node)


	kkk=ndmes/layer
! radi=radi+1
	cont=ncelsepi+1

	ii=ndepi		!node counter
	jj=ncelsepi		!cel counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! central nodes (=nucleus) !!!!

   do l=1,layer  
!	jj=(kkk*(l-1))+radi+ndepi+1
	ii=ii+1
	jj=jj+1
!	node(ii)%x=(sqrt(di/2d0))*mod(l+1,2) ; node(ii)%y=(sqrt(di/2d0))*mod(l+1,2) ! Euler packing algorhytm
	node(ii)%x=0d0 ; node(ii)%y=0d0 ! Euler packing algorhytm
	node(ii)%z=zmes-di*(l-1)	!zmes marks the z-position of the uppermost layer
!	print*,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
    node(ii)%tipus=3 ; node(ii)%icel=jj                        ! origin	
	cels(jj)%node(1)=ii
!	print*,"cel central, node central. ii:",ii
	ent=ii
	ont=ii
	!fill with the external nodes
	if(radi>1)then
        signo=1                                                              ! miguel 4-6-13
		do k=2,radi
			ii=ii+1
            !print*,"cel perif. ent:",ent,ii,k-1, 
             call random_number(a)
			der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
121         call random_number(a)
            xx=der*(1d0-2*a) ; 
            call random_number(a)
            yy=der*(1d0-2*a)            ! miguel 4-6-13 
            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
            if(zz.lt.0)then;goto 121;endif                                   ! miguel 4-6-13
            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz                ! miguel 4-6-13
            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz                ! miguel 4-6-13
  	        else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if       ! miguel 4-6-13       
			!a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de   ! miguel 4-6-13     
			node(ii)%x=node(ent)%x+a
			node(ii)%y=node(ent)%y+b
			node(ii)%z=node(ent)%z+c
		    node(ii)%tipus=3 ; node(ii)%icel=jj
			cels(jj)%node(k)=ii
		end do
	end if

	do j=2,radicel                                                        ! "cicle"
		dx1=0.0 ; dx2=0.0 ; dx3=0.0 ; dy1=0.0 ; dy2=0.0 ; dy3=0.0 
		do i=1,6                                                        ! "sector" of the hexagon
			if(i.eq.1)then   

				dx2=0!*cos(real(i)*60d0*rad) 
				dy2=di*(real(j)-1d0)!*sin(real(i)*60d0*rad);print*,"dx2",dx2,"dy2",dy2
				dx1=di*(real(j)-1d0)*sin(-60d0*rad) ; dy1=di*(real(j)-1d0)*cos(-60d0*rad)!;print*,"dx1",dx1

			else
				hip=di*(real(j)-1d0)                                    ! hipotenusa
				dx1=dx2 ; dy1=dy2
				dx2=hip*sin(real(i-1)*60d0*rad) ; dy2=hip*cos(real(i-1)*60d0*rad)          
			end if            
			ii=ii+1
			jj=jj+1
!			print*,"cel perif, node central. ent:",ent,ii
			node(ii)%x=node(ont)%x+dx2 ; node(ii)%y=node(ont)%y+dy2 ; node(ii)%z=node(ont)%z
!			print*,"vertex i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
!			print*,i,"ii",node(ii)%x,"ent",node(ent)%x
			node(ii)%icel=jj ; node(ii)%tipus=3
			cels(jj)%node(1)=ii
			ent=ii
			if(radi>1)then
                signo=1                                                                  ! miguel 4-6-13
				do k=2,radi
					ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
122                         call random_number(a)
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 122;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13

!					a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13     
					node(ii)%x=node(ent)%x+a
					node(ii)%y=node(ent)%y+b
					node(ii)%z=node(ent)%z+c
				    node(ii)%tipus=3 ; node(ii)%icel=jj
					cels(jj)%node(k)=ii
				end do
			end if

			if(j.gt.2)then                                              ! intermediate points
				dx3=dx2-dx1       ; dy3=dy2-dy1                         ! vectors which link "extreme" points
				dx3=dx3/real(j-1) ; dy3=dy3/real(j-1)                   ! sub-vectors                    
				do k=1,j-2                               
					ii=ii+1
					jj=jj+1
					node(ii)%x=node(ont)%x+(dx1+(real(k)*dx3)) 
					node(ii)%y=node(ont)%y+(dy1+(real(k)*dy3))    
					node(ii)%z=node(ont)%z                                                
!					print*,"intra i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
					node(ii)%icel=jj ; node(ii)%tipus=3 ; cels(jj)%node(1)=ii
					ent=ii
					if(radi>1)then
                        signo=1											            ! miguel 4-6-13
						do kk=2,radi
							ii=ii+1
!							print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
							der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a) 
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
                        
							!a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13     
							node(ii)%x=node(ent)%x+a
							node(ii)%y=node(ent)%y+b
							node(ii)%z=node(ent)%z+c
						    node(ii)%tipus=3 ; node(ii)%icel=jj
							cels(jj)%node(kk)=ii
						end do
					end if
				end do
			end if

		end do
	end do
  end do




	!define the cell's centroid and nucleus (margin)			!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels											!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%ctipus=3									!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0			!>>>>>>>>>>>>>> Miquel 14-4-13
!		print*,"i",i,"node",cels(i)%node(:),"ncels",ncels
		do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Miquel 14-4-13
			k=cels(i)%node(j)								!>>>>>>>>>>>>>> Miquel 14-4-13
!			print*,"k",k,"nunodes",cels(i)%nunodes
			cels(i)%cex=cels(i)%cex+node(k)%x				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cey=cels(i)%cey+node(k)%y				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cez=cels(i)%cez+node(k)%z				!>>>>>>>>>>>>>> Miquel 14-4-13
		end do												!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
	end do													!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels
          b=1.0d8											!>>>>>>>>>>>>>> Is 14-9-13
   	  do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Is 14-9-13
            k=cels(i)%node(j)								!>>>>>>>>>>>>>> Is 14-9-13
   	    a=sqrt((cels(i)%cex-node(k)%x)**2+(cels(i)%cey-node(k)%y)**2+(cels(i)%cez-node(k)%z)**2)
            if (b>a) then ; b=a ; ii=k ; jj=j ; end if
	  end do												!>>>>>>>>>>>>>> Is 14-9-13
          node(ii)%marge=0

          !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster  
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
          cels(i)%node(jj)=jjj
	end do													!>>>>>>>>>>>>>> Is 14-9-13

    node(:)%talone=0.0d0

end subroutine mesenq

!**********************************************************************************************************

subroutine mesenq_packed(radi,radicel,layer,zepi,dreq)
integer            ::radi,radicel,valcel,layer,valcelo
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi,dreq

    di=2d0*dreq
    de=2d0*dreq  !miguel 14-10-13
    vector=0.0d0
    
    node(ndepi+1:nd)%marge=1
 

	beta=2d0*pi/6d0


	do i=ncelsepi+1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


	ii=ndepi
	do i=1,radi
		ii=ii+1
		alt=real(i-1)*di+zepi
!        if(i==1)then
!          alt=zepi
!        else
!          alt=zepi+di
!        end if
!		print*,"alt",alt
		node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
        if(i==1) node(ii)%marge=0
        node(ii)%tipus=3
		cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
		cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
		node(ii)%icel=ncelsepi+1
		do j=2,radi
			angle=beta
!			print*,"angle",angle
			d=de*(j-1d0)
			ii=ii+1
			node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
			node(ii)%tipus=3
			cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
			cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
			node(ii)%icel=ncelsepi+1
			p1(1)=d;p1(2)=0;p1(3)=alt
			p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
!				print*,"jj",jj
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					node(ii)%tipus=3
					cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
					cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
					node(ii)%icel=ncelsepi+1
				end do
			end if
			do k=2,5
				angle=angle+beta
				p1=p2
				ii=ii+1
				node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
				node(ii)%tipus=3
				cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
				cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
				node(ii)%icel=ncelsepi+1
				p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
				jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
				if(jj>0)then
					vector=p2-p1
					modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
					iii=jj+1
					vector=vector/iii
					do kk=1,jj
						ii=ii+1
						node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
						node(ii)%tipus=3
						cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
						cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
						node(ii)%icel=1
					end do
				end if
			end do
			angle=angle+beta
!			print*,"angle",angle
			p1=p2
			ii=ii+1
			node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
			node(ii)%tipus=3
			cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
			cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
			node(ii)%icel=ncelsepi+1
			p2(1)=d;p2(2)=0;p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					node(ii)%tipus=3
					cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
					cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
					node(ii)%icel=ncelsepi+1	
				end do
			end if
		end do
	end do

!hem fet la 1a cel al centre 

!ara farem la resta copiant i pegant


	valcel=ncelsepi+1  !>>>Miquel28-1-14
if(radicel>1)then
	a=2*sqrt(((radi-1)*de)**2-(0.5*de*(radi-1))**2)+de	!distancia de centre de cel·lula a centre de cel·lula
	do i=2,radicel
!		angle=0.5*pi		!"ok" hexagonal configuration
		angle=0.5*pi+pi/12	!perfect hexagonal configuration
!		print*,"i",i,"angle",angle
!		vector(1)=0;vector(2)=a*(i-1);vector(3)=0			!"ok" hexagonal configuration
		vector(1)=-a*(i-1)*sin(pi/12);vector(2)=a*(i-1)*cos(pi/12);vector(3)=0	!perfect hexagonal configuration
		valcel=valcel+1
		do j=1,cels(ncelsepi+1)%nunodes
			ii=ii+1
            if (node(ndepi+j)%marge==0) node(ii)%marge=0
			node(ii)%x=node(ndepi+j)%x+vector(1);node(ii)%y=node(ndepi+j)%y+vector(2);node(ii)%z=node(ndepi+j)%z
			node(ii)%tipus=node(ndepi+j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do

		p1(1)=a*(i-1)*dcos(angle-beta);p1(2)=a*(i-1)*dsin(angle-beta);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
			iii=jj+1
			vector=vector/iii
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(ncelsepi+1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (node(ndepi+j)%marge==0) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(ndepi+j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
		p1=p2
		do k=2,5
			angle=angle+beta
!			print*,"angle",angle
			p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
			valcel=valcel+1
			do j=1,cels(ncelsepi+1)%nunodes
				ii=ii+1
                if (node(ndepi+j)%marge==0) node(ii)%marge=0
				node(ii)%x=node(ndepi+j)%x+p2(1);node(ii)%y=node(ndepi+j)%y+p2(2);node(ii)%z=node(ndepi+j)%z
				node(ii)%tipus=node(ndepi+j)%tipus
				cels(valcel)%nunodes=cels(valcel)%nunodes+1
				cels(valcel)%node(cels(valcel)%nunodes)=ii
				node(ii)%icel=valcel
			end do
			jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!				iii=jj+1
				vector=vector/(jj+1)
				do kk=1,jj
					valcel=valcel+1
					do j=1,cels(ncelsepi+1)%nunodes
						kkk=cels(valcel-1)%node(j)
						ii=ii+1
                        if (node(ndepi+j)%marge==0) node(ii)%marge=0
						node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
						node(ii)%tipus=node(ndepi+j)%tipus
						cels(valcel)%nunodes=cels(valcel)%nunodes+1
						cels(valcel)%node(cels(valcel)%nunodes)=ii
						node(ii)%icel=valcel
					end do
				end do
			end if
			p1=p2
		end do
		angle=angle+beta
!		p2(1)=0;p2(2)=a*(i-1);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		valcel=valcel+1
		do j=1,cels(ncelsepi+1)%nunodes
			ii=ii+1
            if (node(ndepi+j)%marge==0) node(ii)%marge=0
			node(ii)%x=node(ndepi+j)%x+p2(1);node(ii)%y=node(ndepi+j)%y+p2(2);node(ii)%z=node(ndepi+j)%z
			node(ii)%tipus=node(ndepi+j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!			iii=jj+1
			vector=vector/(jj+1)
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(ncelsepi+1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (node(ndepi+j)%marge==0) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(ndepi+j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
	end do
end if

valcelo=valcel
if(layer>1)then
  do i=2,layer
    kk=0
    do jj=1,mesradicel(i)-1
      kk=kk+jj
    end do
    kk=(6*kk+1)  !number of mesenchymal cells
    do j=ncelsepi+1,ncelsepi+kk
      valcel=valcel+1
      do k=1,cels(j)%nunodes
        kk=cels(j)%node(k) !;print*,"cel",valcel,"node",kk
        ii=ii+1
        if (node(kk)%marge==0) node(ii)%marge=0
        node(ii)%x=node(kk)%x;node(ii)%y=node(kk)%y;node(ii)%z=node(kk)%z-radi*(i-1)*di
        node(ii)%tipus=node(kk)%tipus
        node(ii)%icel=valcel
        cels(valcel)%nunodes=cels(valcel)%nunodes+1
        cels(valcel)%node(cels(valcel)%nunodes)=ii
      end do
    end do
  end do
end if


	!define the cell's centroid

	do i=ncelsepi+1,ncels
		cels(i)%ctipus=3
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
		do j=1,cels(i)%nunodes
			k=cels(i)%node(j)
 			cels(i)%cex=cels(i)%cex+node(k)%x
			cels(i)%cey=cels(i)%cey+node(k)%y
			cels(i)%cez=cels(i)%cez+node(k)%z
		end do
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)
	end do


!	do i=1,ndepi		!veins i parametres mecanics
!		ii=node(i)%icel
!		if(node(i)%tipus==1)then
!   		  node(i)%altre=i-cels(1)%nunodes/2
!		else
!		  node(i)%altre=i+cels(1)%nunodes/2
!		end if
!	end do

!        node(1)%marge=0.0d0

        !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster >>> Is 14-9-13 
	do i=ncelsepi+1,ncels
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (node(k)%marge==0) then
              ii=k
              exit
            end if
	  end do										
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (k==ii) then
              jj=j
              exit
            end if
	  end do									
          cels(i)%node(jj)=jjj
	end do											


    !some funny rotation to make a less biased diffusion lost
!    do i=1,nd
!      c=0.75
!      a=node(i)%x*cos(c)-node(i)%y*sin(c)
!      b=node(i)%x*sin(c)+node(i)%y*cos(c)
!      node(i)%x=a
!      node(i)%y=b
!    end do

    node(:)%talone=0.0d0


!do i=ndepi+1,nd
! print*,"tipus",node(i)%tipus
!end do

end subroutine




subroutine matrix(radi,layer,zepi)		!to place a block of ECM
integer            ::radi,radicel,valcel,layer
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi

    di=1.0d0
    de=2d0*node(1)%eqd  !miguel 14-10-13
    vector=0.0d0
    beta=2d0*pi/6d0



    ii=ndepi+ndmes
    node(ii+1:nd)%marge=1  !so there are no nuclei among ECM

    do i=1,layer
      ii=ii+1
      alt=zepi-real(i-1)*2*node(1)%eqd
      node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
      node(ii)%tipus=4
      node(ii)%icel=-ii
      do j=2,radi
        angle=beta
        d=de*(j-1d0)
        ii=ii+1
        node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
        node(ii)%tipus=4
        node(ii)%icel=-ii
        p1(1)=d;p1(2)=0;p1(3)=alt
        p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
        jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
        if(jj>0)then
          vector=p2-p1
          modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
          iii=jj+1
          vector=vector/iii
          do kk=1,jj
            ii=ii+1
            node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
            node(ii)%tipus=4
            node(ii)%icel=-ii
          end do
        end if
        do k=2,5
          angle=angle+beta
          p1=p2
          ii=ii+1
          node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
          node(ii)%tipus=4
          node(ii)%icel=-ii
          p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
          jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
          if(jj>0)then
            vector=p2-p1
            modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
            iii=jj+1
            vector=vector/iii
            do kk=1,jj
              ii=ii+1
              node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
              node(ii)%tipus=4
              node(ii)%icel=-ii
            end do
          end if
        end do
        angle=angle+beta
        p1=p2
        ii=ii+1
        node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
        node(ii)%tipus=4
        node(ii)%icel=-ii
        p2(1)=d;p2(2)=0;p2(3)=alt
        jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
        if(jj>0)then
          vector=p2-p1
          modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
          iii=jj+1
          vector=vector/iii
          do kk=1,jj
            ii=ii+1
            node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
            node(ii)%tipus=4
            node(ii)%icel=-ii
          end do
        end if
      end do
    end do

    node(:)%talone=0.0d0

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine mesenq_cell_sorting(radi,radicel)									!!! miguel 8-7-2013 

integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,signo1,ent,ont       ! number of layers and concentric hexagons 
real*8   ::rad,der,dar,zmes,de,di
real*8   :: ka,kb,kc,xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
de=0.4d0                     ! cell radius
di=5.5d0*de                    ! radius of the bulk of cells

	do i=ncelsepi+1,ncels               !	
		cels(i)%nunodes=radi				!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela     !    
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%node=0
    !cels(i)%tipcel=1       
  end do
         signo1=1 ; ii=1
         do i=1,radicel						
                            call random_number(a)
						dar=(sqrt(di)*sqrt(di*a))                            ! miguel 4-6-13
122         call random_number(a)         
            xx=dar*(1d0-2*a) ; 
            call random_number(a)
            yy=dar*(1d0-2*a)         ! miguel 4-6-13 
            zz=(dar**2)-(xx**2)-(yy**2) 			               						  ! miguel 4-6-13
            if(zz.lt.0)then;goto 122;endif                                ! miguel 4-6-13
            zz=signo1*sqrt(zz) ; signo1=-1*signo1                            ! miguel 4-6-13  
            if(mod(i,3).eq.0)then      ; kb=xx ; kc=yy ; ka=zz              ! miguel 4-6-13
            else if(mod(i,3).eq.1)then ; ka=xx ; kc=yy ; kb=zz              ! miguel 4-6-13
  	        else if(mod(i,3).eq.2)then ; ka=xx ; kb=yy ; kc=zz ; end if     ! miguel 4-6-13
            signo=1
            call random_number(a) 
						do kk=1,radi			
            call random_number(a)			
							der=(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a)           
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
							node(ii)%x=ka+a
							node(ii)%y=kb+b
							node(ii)%z=kc+c
					    node(ii)%tipus=3 ; node(ii)%icel=i !;write(*,*)'cel',i,'node',ii,'en orden',kk,'r,rcel',radi,radicel
							cels(i)%node(kk)=ii ; 	ii=ii+1
						end do
         end do

	!define the cell's centroid								!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels 							!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%ctipus=3									!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0			!>>>>>>>>>>>>>> Miquel 14-4-13
		do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Miquel 14-4-13
			k=cels(i)%node(j)								!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cex=cels(i)%cex+node(k)%x				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cey=cels(i)%cey+node(k)%y				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cez=cels(i)%cez+node(k)%z				!>>>>>>>>>>>>>> Miquel 14-4-13
		end do												!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
	end do													!>>>>>>>>>>>>>> Miquel 14-4-13

    prop_noise=0.1d0    
    realtime=0
    reqmin=0.05d0 

    mmae=node(1)%eqd

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine mesenq_cell_sorting
!*************************************************************************************************************

!subroutine epiteli_marge(radi,radicel)
!integer            ::radi,radicel,valcel
!real*8,dimension(:)::p1(3),p2(3),vector(3)
!real*8             ::beta,angle,alt,de,di
!
!    call epiteli(radi,radicel,1d0)
!
!
!
!	a=sqrt(0.05d0)	!quantitat convenient de molecules d'adhesio intercel pk la força sigui la meitat que la dels nodes de la mateixa cell 
!
!	do i=1,nd		!veins i parametres mecanics
!		node(i)%you=1d1;node(i)%adh=1d1
!		node(i)%rep=1d2;node(i)%rec=1d2
!		node(i)%eqd=0.5d0 !;node(i)%codel=0.5d0 !de
!                node(i)%eqs=1d0
!		node(i)%add=node(i)%eqd*1.25!+da 
!                node(i)%hoo=1d1
!		ii=node(i)%icel
!		if(node(i)%tipus==1)then
!   		  node(i)%altre=i-nodecel/2
!		else
!		  node(i)%altre=i+nodecel/2
!		end if
!                node(i)%erp=1d0
!	end do
!        do i=1,nd
!          node(i)%marge=1
!          if (node(i)%tipus==1) then
!            d=sqrt(node(i)%x**2+node(i)%y**2+(node(i)%z-1)**2)
!            if (d<4d0) then 
!              node(i)%marge=0
!              ii=node(i)%icel
!              do j=1,nd
!                if (node(j)%icel==ii) node(j)%marge=0
!              end do
!            end if
!          else
!            d=sqrt(node(i)%x**2+node(i)%y**2+node(i)%z**2)
!            if (d<4d0) then 
!              node(i)%marge=0
!              ii=node(i)%icel
!              do j=1,nd
!                if (node(j)%icel==ii) node(j)%marge=0
!              end do
!            end if
!          end if
!        end do
!        do i=1,nd
!          ii=node(i)%icel
!          if (node(i)%marge==0) then
!            do j=1,nd
!              if (node(j)%icel==ii) node(j)%marge=0
!            end do
!          end if
!        end do
!
!  !physical
!  getot=0
!  temp=1.1  !high value is low temperature
!  !mathematical
!  itacc=0
!  nparti=1000
!  desmax=0.0001
!        resmax=1d-3
!  idum=-11111
!  idumoriginal=idum
!  rv=0.75d0	!el maxim rang d'abast que tindra un node
!  interfmax=1d0
!
!
!  !biological
!    !functions used
!    
!    
!    ffu=0
!     !spring of the ellipse
!    ffu(2)=0 !screening
!    ffu(3)=0 !torsion
!     !miguel4-11-13
!    
!    !number of params
!    !nparam=32
!    nvarglobal_out=5
!      angleto=0      !no es el parametre en si si no una forma de veure'l
!    angltormax=dcos(angleto*2*pi/360) !the param itself, a dcos not an angle   
!
!    !
!    ntipusadh=1
!
!
!    if (allocated(kadh)) deallocate(kadh)
!    allocate(kadh(ntipusadh,ntipusadh))
!    kadh=0d0
!    gen(1)%e(1)=1
!    prop_noise=0.01d0    
!    realtime=0
!    reqmin=0.05d0 
!
!    mmae=node(1)%eqd
!
!    cels(:)%fase=0.0d0
!    cels(:)%minsize_for_div=cels(:)%nunodes*2
!    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
!    dmax=2
!    node(:)%talone=0.0d0
!    ramax=maxval(node(:)%add)*3
!    node(:)%dif=0.0d0
!end subroutine

!***************************************************************************************************************

subroutine epitelitest					!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
integer            ::radi,radicel,valcel
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di


	radi=2
	radicel=1


	nodecel=4
 	ncels=1
	ncals=ncels+10
    nd=4
    nda=nd+2
    if (allocated(node)) deallocate(node)
    allocate(node(nda))
    if (allocated(cels)) deallocate(cels)
    allocate(cels(ncals))

	node(:)%x=0d0;node(:)%y=0d0;node(:)%z=0d0

	do i=1,ncals
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecel))
	end do


   tacre=1000

	node(1)%x=-0.125d0;node(1)%y=0d0;node(1)%z=0d0;	node(1)%tipus=2;	node(1)%altre=3
	node(2)%x=0.625d0;node(2)%y=0d0;node(2)%z=0d0;	node(2)%tipus=2;	node(2)%altre=4
	node(3)%x=0.0d0;node(3)%y=0d0;node(3)%z=1.0d0;	node(3)%tipus=1;	node(3)%altre=1
	node(4)%x=0.5d0;node(4)%y=0d0;node(4)%z=1.0d0;	node(4)%tipus=1;	node(4)%altre=2
	node(:)%icel=1


	do i=1,nd		!veins i parametres mecanics
		node(i)%you=1d4;node(i)%adh=1d2
		node(i)%rep=1d4;node(i)%rec=1d3
		node(i)%eqd=0.5d0 !;node(i)%codel=0.5d0
		node(i)%add=node(i)%eqd*1.5d0
		node(i)%eqs=1.0d0
        node(i)%hoo=1d5
        node(i)%erp=0d0
        !node(i)%est=1d3
		node(i)%eqd=0.5d0!;node(i)%codel=0.5d0
	end do




	!!!!!!!posició experimental dels nodes
	node(1)%eqd=0.75d0;node(2)%eqd=0.75d0


!		node(1)%x=node(1)%x-0.2d0;node(2)%x=node(2)%x+0.2d0

!		node(2)%x=0.1d0;node(2)%y=0d0;node(2)%z=-0.5d0
!		node(4)%x=0.1d0;node(4)%y=0d0;node(4)%z=0.5d0

!		node(3)%eqd=0.3d0;node(4)%eqd=0.3d0;node(3)%codel=0.3d0;node(4)%codel=0.3d0
!		node(1)%eqd=1d0;node(2)%eqd=1d0;node(1)%codel=1d0;node(2)%codel=1d0
!		node(1)%add=node(1)%eqd*1.5;node(2)%add=node(2)%eqd*1.5
!		node(:)%add=node(1)%eqd*1.5



!		node(3)%eqd=0.5d0*node(3)%eqd;node(3)%codel=0.5d0*node(3)%codel
!		node(1)%eqd=1.2d0*node(1)%eqd;node(1)%codel=1.2d0*node(1)%codel


	!define the cell's centroid


		cels(1)%ctipus=1
		cels(1)%cex=0;cels(1)%cey=0;cels(1)%cez=0
		do j=1,cels(1)%nunodes
			k=cels(1)%node(j)
			cels(1)%cex=cels(1)%cex+node(k)%x
			cels(1)%cey=cels(1)%cey+node(k)%y
			cels(1)%cez=cels(1)%cez+node(k)%z
		end do
		cels(1)%cex=cels(1)%cex/4d0
		cels(1)%cey=cels(1)%cey/4d0
		cels(1)%cez=cels(1)%cez/4d0





  !physical
  getot=0
  temp=0.1d-1  !high value is low temperature
  !mathematical
  itacc=0
  nparti=1000
  desmax=0.01
        resmax=1d-3
  idum=-11111
  idumoriginal=idum
  rv=1.5d0	!el maxim rang d'abast que tindra un node
  interfmax=1d0


  !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
     ! miguel4-11-13
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise

    !number of params
    !nparam=32
    nvarglobal_out=5

    ntipusadh=1
    gen(1)%e(1)=1
    if (allocated(kadh)) deallocate(kadh)
    allocate(kadh(ntipusadh,ntipusadh))
realtime=0
    reqmin=0.05d0 

    mmae=node(1)%eqd

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!*********************************************************************************************************************

subroutine basal(side)
integer   :: side, step, contt, ndq 
real*8    :: limx, limy                            
               
limx=(side*sqrt(7.5d-1*(de**2)))/2d0 ; limy=((de/2d0)*mod(side,2)+side*de)/2d0 ! centre of basal nodes square
ndq=side**2                                         ! total number of nodes
 
      cels(ncels+1)%nunodes=ndq
      allocate(cels(ncels+1)%node(ndq))

   cont=ndmes+1
   do i=1,side                                      ! rows
   xx=(de/2d0)*mod(i,2)
       do j=1,side                                  ! columns
       node(cont)%x=(i*sqrt(7.5d-1*(de**2)))-limx 
       node(cont)%y=(xx+(j*de))-limy 
       node(cont)%z=di-de-0.01                      ! -differential (in order to avoid mesenchimal and basal nodes in the same plane)
       if(j.gt.1)then;node(cont)%altre=cont-1 
       else ;         node(cont)%altre=cont+1 ;endif 
       node(cont)%icel=ncels+1   
       contt=cont-ndmes
       cels(ncels+1)%node(contt)=cont         
           !!!!!!!!!!!!!!!!!!!!!!!!!              
  	   step=contt-1              		
           step=contt+1		
           step=contt+side		
           step=contt-side  
             if(mod(i,2).eq.0)then
             step=(contt+side)-1            
	     step=(contt-side)-1         	
             else 
             step=(contt+side)+1        
             step=(contt-side)+1
	     end if                       
           node(cont)%tipus=4 
	   !!!!!!!!!!!!!!!!!!!!!!!!!          
           cont=cont+1
       end do
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=ndmes+1,nd 		               !veins i parametres mecanics
		node(i)%you=1d0;node(i)%adh=1d0
		node(i)%rep=1d0;node(i)%rec=1d0
		node(i)%eqd=1d0 !;node(i)%codel=0.1d0 !de
                node(i)%eqs=1d0
		node(i)%add=2*de
                node(i)%hoo=1d1		
	    end do

    prop_noise=0.01d0    
realtime=0
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine basal

!***********************************************************************

subroutine escribe
 open(666,file="nodoscels.dat",action='write')
    do i=1,size(node)
    write(666,*)'node',i,'tipus',node(i)%tipus,'cel',node(i)%icel,'xyz',node(i)%x,node(i)%y,node(i)%z
    end do
    do i=1,size(cels)
    write(666,*)'cel',i,'nodes',cels(i)%node(:)
    end do
  close(666)
end subroutine escribe

!***********************************************************************

subroutine migration

print *,""
print *,"WARNING THIS IS GOING TO TAKE A WHILE"
print *,""


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=3    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=14     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=1.4d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.0 !1d-3 !low value is low temperature	
    desmax=0.0d0
    resmax=1d-3
    prop_noise=0.05d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.9d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 1= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta
    maxad=5.00d0  !!>> HC 16-12-2020

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=5d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1 !>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.28d0	!>>Miquel 26-10-12!de
        node(i)%add=node(i)%eqd*1.40 
        node(i)%eqs=0.25d0
        node(i)%hoo=1d1
        node(i)%erp=1d1
        node(i)%est=1d1
        node(i)%mov=0d0
        node(i)%dmo=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=6d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=2d1;node(i)%rec=3.0d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ;node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*2.0d0 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=1d0
        node(i)%dmo=1d-1
      end do
    end if




    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
      !Special modifications
      do i=ndepi+1,nd
        node(i)%y=node(i)%y+1.50d0
      end do


    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=3
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0
       gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0   !mesenchyme
       gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0   !epithelium


    !Gene-behavior interactions
      gen(1)%e(1)=1  !adhesion
      gen(2)%e(1)=2  !adhesion
      gen(3)%e(1)=3  !adhesion

      !gen(2)%e(6)=0.1d0           !this is da (keep expression levels at one so this will be the actual value
      !gen(2)%e(16)=0.1d0   !this is dmo (keep expression levels at one so this will be the actual value

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

       !kadh(1,2)=3d1 !adhesion table
       kadh(2,1)=10d0 ; kadh(1,2)=kadh(2,1)!adhesion table THIS HAS TO BE SYMMETRIC
       kadh(3,3)=10.0d1
       kadh(2,2)=10d0  ! HI HA ALGO QUE NO VA
       kadh(3,2)=0d0  ! HI HA ALGO QUE NO VA
       kadh(2,3)=0d0  ! HI HA ALGO QUE NO VA


    end if

    !Gene expression on nodes
       do i=1,ndepi
         gex(i,3)=1.0d0
       end do
       do i=ndepi+1,nd
         gex(i,2)=1d0
         !gex(i,3)=0d0
       end do

       !adhesive concentration gradient in epithelium
       a=minval(node(:)%x)
       do i=1,ndepi
         if (node(i)%tipus==1) then
           gex(i,1)=1d0*((node(i)%y+a))**2
         end if
       end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
 
    do i=1,ndepi
      node(i)%fix=2
    end do

!    node(ndepi+1:nd)%z=node(ndepi+1:nd)%z-0.3

end subroutine

!************************************************************************************************************************

subroutine epi_apoptosis   !>>>>>>>>>>>> Miquel 18-6-13


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters




    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d1 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.3d0
    deltamax=1d-2 ! miguel 14-10-13
    dmax=1
    screen_radius=1.0d0


    !biological
    !functions used
    
    
    ffu=0
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta
    maxad=5.0d0
     !>>> TT
    
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=5.0d0;node(i)%adh=5d0	!>>Miquel 26-10-12  !!>>HC 13-7-2021
        node(i)%rep=5d1;node(i)%rec=5d1	!>>Miquel 8-10-12   !!>>HC 13-7-2021
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
		node(i)%add=node(i)%eqd*1.50 
        node(i)%hoo=1d1
        node(i)%erp=20d0  !!>>HC 13-7-2021
        node(i)%est=60d0  !!>>HC 13-7-2021
        node(i)%mov=1d0
        node(i)%dmo=5d-3
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=0d0 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0


    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+3)=5d-1   !apoptosis

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      do i=1,cels(1)%nunodes ; gex(cels(1)%node(i),1)=1.0d0 ; end do
      do i=1,cels(8)%nunodes ; gex(cels(8)%node(i),1)=1.0d0  ; end do
      do i=1,cels(10)%nunodes ; gex(cels(10)%node(i),1)=1.0d0  ; end do
      do i=1,cels(12)%nunodes ; gex(cels(12)%node(i),1)=1.0d0  ; end do
      do i=1,cels(14)%nunodes ; gex(cels(14)%node(i),1)=1.0d0  ; end do
      do i=1,cels(16)%nunodes ; gex(cels(16)%node(i),1)=1.0d0  ; end do
      do i=1,cels(18)%nunodes ; gex(cels(18)%node(i),1)=1.0d0  ; end do

      gex(:14,1)=1.0d0
      do i=1,cels(1)%nunodes
        j=cels(1)%node(i)
        gex(j,1)=1.0d0
      end do



    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
    

end subroutine epi_apoptosis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine epiteli_sphere(radi,radicel,radius)		!radi is the number of cell rows from one pole to the equator
integer            ::radi,radicel,valcel

integer::trobat,cont,rep,n,i,j,k,ii,jj,kk,lev,numdos,ncels2,val,val2,ll,val3,val4,val5,val6,iiii,jjjj,kkkk,suma,cocels
real*8::modul,pescab,pescac,pescbc,modula,modulb,modulc,aax,aay,aaz,bbx,bby,bbz,ccx,ccy,ccz,costat,rhex,lad,ax,ay,az
real*8::alf,bet,gam,a,b,c,l,aa,bb,cc,aaa,bbb,ccc,angle,rpent,modmin,pesc,radius,radius2,bx,by,bz,cost,ucost,sint,thet
real*8,dimension(:)::minu(3),vec1(3),vec2(3),vec3(3),vec4(3),vec5(3)
real*8,allocatable::cmalla(:,:),cveci(:,:),primers(:)




    di=1.0d0
    de=0.5d0

    radicel=(radicel-1)*2

    !SPHERE CODE, FIRST HEMISPHERE
	cont=0 ; cocels=0
	ii=radicel/2
!	radius=5d0  !radius of the sphere
    radius2=radius+di

	alf=pi/radicel

!!!pole cell
	cont=cont+1 ; cocels=cocels+1
    node(1)%x=0d0 ; node(1)%y=0d0 ; node(1)%z=radius
    node(1)%altre=2 ; node(1)%tipus=2 ; node(1)%icel=cocels
    cont=cont+1
    node(2)%x=0d0 ; node(2)%y=0d0 ; node(2)%z=radius2
    node(2)%altre=1 ; node(2)%tipus=1 ; node(2)%icel=cocels

    !el radi 2 de la cèlula

    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius
    node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius2
    node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

    gam=2*pi/6d0
    do i=1,5
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius
      node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius2
      node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
    end do

!!!!! end pole cell


	do i=1,ii

        !!!!! spine cell
!  print*,"*********espina******",i
		cont=cont+1 ; cocels=cocels+1
		suma=0
		jj=cont
		angle=alf*i
!		print*,"cont",cont
		vec1(1)=radius*sin(angle)
		vec1(2)=0
		vec1(3)=radius*cos(angle)
		node(cont)%x=vec1(1) ; node(cont)%y=vec1(2) ; node(cont)%z=vec1(3)
        node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!        print*,"altre",node(cont)%altre

        cont=cont+1
		vec2(1)=radius2*sin(angle)
		vec2(2)=0
		vec2(3)=radius2*cos(angle)
		node(cont)%x=vec2(1) ; node(cont)%y=vec2(2) ; node(cont)%z=vec2(3)
        node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

        ux=vec2(1)-vec1(1) ; uy=vec2(2)-vec1(2) ; uz=vec2(3)-vec1(3)
        a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a

        ax=-de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"vora"," a",ax,ay,az


        do j=1,6
          cont=cont+1
          thet=j*gam
!   print*,"thet",thet
          cost=cos(thet); ucost=1-cost ; sint=sin(thet)
!  print*,"cos",cost,"ucost",ucost,"sint",sint
          bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
          by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
          bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
          node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
          node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!   print*,"b",bx,by,bz
!  print*,"vec1",vec1
!  print*,"nodecont",node(cont)%x,node(cont)%y,node(cont)%z
          cont=cont+1

          node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
          node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
        end do

        !!!!!!!!! end spine cell



!        print*,"altre",node(cont)%altre

		if(i+1<=ii+1)then
			kk=6+6*(i-1)
		else
			kk=6+6*(radicel-i-1)
		end if
		bet=2*pi/kk

!		print*,"i",i,"nombre de cels del paralel",kk

		do j=1,kk-1

            !!!!!!!!!!! rib cell

			angle=alf

			cont=cont+1 ; cocels=cocels+1
			a=node(1)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
			b=node(1)%y
			c=node(1)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
			a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
			b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
			c=node(cont)%z
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
            node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
            vec1(1)=node(cont)%x ; vec1(2)=node(cont)%y ; vec1(3)=node(cont)%z

			cont=cont+1
			a=node(2)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
			b=node(2)%y
			c=node(2)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
			a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
			b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
			c=node(cont)%z
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
            node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

            ux=node(cont)%x-node(cont-1)%x
            uy=node(cont)%y-node(cont-1)%y
            uz=node(cont)%z-node(cont-1)%z
            a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a !;print*,"spring",a
            vec2(1)=node(cont)%x ; vec2(2)=node(cont)%y ; vec2(3)=node(cont)%z

            ax=-de*sin(j*bet) ; ay=de*cos(j*bet) ; az=0d0

!        ax=de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"bet",j*bet,"a",ax,ay,"modul",sqrt(ax**2+ay**2)
            do k=1,6
              cont=cont+1
              thet=k*gam
              cost=cos(thet); ucost=1-cost ; sint=sin(thet)
              bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
              by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
              bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
              node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
              node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels

!              print*,"b rib",bx,by,bz,"modul",sqrt(bx**2+by**2+bz**2)

              cont=cont+1

              node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
              node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
           end do


            !!!!!!!!!!! end rib cell



		end do
	end do

!print*,"cont",cont,"cocels",cocels


!the 2ond hemisphere

    j=0
    do i=1,(radicel/2+1)-2
      j=j+i
    end do
    ii=(6*j+1)*2*7   !number of nodes on the 2ond hemisphere

    print*,"ii",ii
    do i=1,ii
      cont=cont+1
      node(cont)%x=node(i)%x
      node(cont)%y=node(i)%y
      node(cont)%z=-node(i)%z
      node(cont)%tipus=node(i)%tipus
      if(node(i)%tipus==2)then
        node(cont)%altre=cont+1
      else
        node(cont)%altre=cont-1
      end if
      node(cont)%icel=node(i)%icel+cocels
    end do

!print*,"cont complet",cont,"cocels",cocels


!	cont=cont+1
!	malla(cont,1)=0;malla(cont,2)=0;malla(cont,3)=-radius
!	primers(radi+1)=cont
!	print*,i,"malla",malla(cont,:)
!	print*,"nombre de cels creades",cont
!	print*,"nombre de paral·lels",lev





	!define the cell's centroid

    do i=1,ncelsepi
      kk=0
      cels(i)%ctipus=1
      cels(i)%nunodes=nodecel
      cels(i)%nodela=nodecela
      cels(i)%cex=0d0 ; cels(i)%cey=0d0 ; cels(i)%cez=0d0
      cels(i)%polx=0d0 ; cels(i)%poly=0d0 ; cels(i)%polz=0d0
      allocate(cels(i)%node(nodecela))
      do j=1,ndepi
        k=node(j)%icel
        if(k==i)then
          kk=kk+1
          cels(i)%node(kk)=j
          if(node(j)%tipus==1)then
            cels(i)%cex=cels(i)%cex+node(j)%x
            cels(i)%cey=cels(i)%cey+node(j)%y
            cels(i)%cez=cels(i)%cez+node(j)%z
          end if
        end if
      end do
      cels(i)%cex=2*cels(i)%cex/cels(i)%nunodes
      cels(i)%cey=2*cels(i)%cey/cels(i)%nunodes
      cels(i)%cez=2*cels(i)%cez/cels(i)%nunodes
    end do


!	do i=1,ncelsepi
!		cels(i)%ctipus=1
!		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
!		do j=1,cels(i)%nunodes
!			k=cels(i)%node(j)
!            if(node(k)%tipus==1)then
!	 			cels(i)%cex=cels(i)%cex+node(k)%x
! 				cels(i)%cey=cels(i)%cey+node(k)%y
!				cels(i)%cez=cels(i)%cez+node(k)%z
!			end if
!		end do
!		cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
!		cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
!		cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
!	end do


	do i=1,ndepi		!veins i parametres mecanics
		node(i)%you=1d1;node(i)%adh=1d1
		node(i)%rep=1d2;node(i)%rec=1d2
		node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0 !de
            node(i)%eqs=0.5d0
		node(i)%add=node(i)%eqd*1.25!+da 
        node(i)%hoo=1d1
!        node(i)%tipus=3
!        node(i)%icel=i
	end do


    realtime=0
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine epiteli_sphere

!**********************************************************************************************

subroutine epi_sphere


	!epithelium's dimension parameters
	radi=1       !number of radial layers of nodes per cell
	radicel=5    !number of radial layers of cells
    radius=3d0   !radius of the sphere
!	zepi=0d0     !z-position of basal layer


	j=0
	do i=1,radi-1
	  j=j+i
	end do
	nodecel=(6*j+1)*2	!number of nodes per cell

	nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

!    nodecel=1
!    nodecela=1


	j=0					!cell count
	do i=1,radicel-1
  	  j=j+i
	end do

    k=0
	do i=1,radicel-2
  	  k=k+i
	end do



	ncelsepi=(6*j+1)+(6*k+1) ; print*,"ncelsepi",ncelsepi,"j",j
    ndepi=nodecel*ncelsepi

	ncels=ncelsepi+ncelsmes

	ncelsmes=0

    nx=0  !number of ECM nodes

    !ndepi=ncelsepi*2

	ndmes=0

	nd=ndepi+ndmes+ndx ; print*,"nd",nd

    nda=nd+10
	ncals=ncels+10



	print*,"nd",nd,"ncels",ncels,"nodecela",nodecela

    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 

    call iniarrays


    call epiteli_sphere(radi,radicel,radius)


  ng=0       !>>>>>> Is 29-4-13

  !physical
  getot=0
  temp=0.1d-5 !low value is low temperature	
  !mathematical
  itacc=0
  nparti=1000
  desmax=0.01
        resmax=1d-3
  idum=-11111
  idumoriginal=idum
  rv=0.5d0*1.5d0	!el maxim rang d'abast que tindra un node   !it is redefined below as function of da
	print*,"rv",rv
  interfmax=1d0
    dmax=1
    screen_radius=1.0d0

  !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
     !miguel4-11-13
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


    !number of params
    !nparam=32
    nvarglobal_out=5

	geu=0;rtime=0
    tacre=5	!growth rate and differential	>>Miquel 8-10-12

                    !no es el parametre en si si no una forma de veure'l

    ntipusadh=1


    if (allocated(kadh)) deallocate(kadh)
    allocate(kadh(ntipusadh,ntipusadh))
    kadh=0d0

    do i=1,nd		!veins i parametres mecanics
      node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
      node(i)%rep=1.3d1;node(i)%rec=1.8d1	!>>Miquel 8-10-12
      node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
      node(i)%eqs=0.5d0
      node(i)%add=node(i)%eqd*1.35  
      node(i)%hoo=1d2
      node(i)%erp=1d0
      !node(i)%est=1d1
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      ii=node(i)%icel

      if(node(i)%tipus==1)then
        node(i)%eqd=0.35d0
!        node(i)%codel=0.35d0
        node(i)%add=node(i)%eqd*1.35d0
      end if
      node(i)%dmo=desmax
      node(i)%mov=temp

    end do

    rv=2*node(1)%add;  print*,"rv",rv

    ng=1

    call initiate_gene
    gen(1)%e(1)=1

    gen(1)%e(15)=1d-4
    do i=1,nd
      gex(i,1)=0.0d0	!no growth
    end do

    call update_npag

    prop_noise=0.1d0    


    realtime=0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine epi_sphere

!***************************************************************************************************************

subroutine epi_mes_ecm


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=15     !number of nodes per cell
	mradicel=4   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=mradi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-4 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    ecmmax=0.25d0
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3

    !biological
    !functions used
    
    
    ffu=0
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta

    ffu(5)=0 !integration method 0=euler , 1=runge-kutta forces , 2=runge-kutta forces+genes
    ffu(12)=0 !adaptive time step 0=no , 1=yes for forces , 2=yes for forces+genes
     !>>> TT
    maxad=0.50d0  !!>> HC 16-12-2020

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=2d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.50
        node(i)%hoo=1d1
        node(i)%erp=5d0
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%pla=1d0
        node(i)%kvol=1d0
        node(i)%ecm=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=2d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
		node(i)%add=node(i)%eqd*1.90; 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=0d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters

      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=0d0 ; gen(1)%npost=0 ; gen(1)%npre=0
!       gen(1)%npre=1   ; allocate(gen(1)%pre(1)) ; gen(1)%pre(1)=5
       gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 ; gen(2)%npost=0 ; gen(2)%npre=0
       gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 ; gen(3)%npost=0 ; gen(3)%npre=0
       gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 ; gen(4)%npost=0 ; gen(4)%npre=0
!       gen(5)%kindof=1 ; gen(5)%diffu=5d0 ; gen(5)%mu=0d0 ; 
!       gen(5)%npost=1   ; allocate(gen(5)%post(1)) ; gen(5)%post(1)=1

    !Gene-behavior interactions
      gen(1)%e(1)=1
      gen(2)%e(1)=2
      !gen(3)%e(1)=3
      gen(4)%e(nparam_per_node+4)=2d-1 
      gen(1)%e(nparam_per_node+5)=1d0
      gen(4)%e(nparam_per_node+6)=2d1
      gen(4)%e(nparam_per_node+7)=0.10d0

    !Gene-gene interactions
      gen(1)%t(4)=1.0d4
      !gene 4 induces the synthesis of the secreted gene

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecule interactions

        kadh(1,1)=1d2
        kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
        !kadh(1,3)=1d0;kadh(3,1)=kadh(1,3)
        kadh(2,2)=1d1
        !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
        !kadh(3,3)=1d1


    end if

    !Gene expression on nodes
      do i=1,nd
        if(node(i)%tipus==3)then
          gex(i,2)=1d0
          !gex(i,1)=1.0d0
          !gex(i,4)=1.0d0;!print*,"i",i
        elseif(node(i)%icel<8.and.node(i)%tipus==2) then
          !gex(i,2)=1d0
          !gex(i,1)=1.0d0
          gex(i,4)=1.0d0;!print*,"i",i
        end if
        if(node(i)%tipus<3)then
          gex(i,2)=1d0
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine epi_mes_ecm

!***************************************************************************************************

subroutine mes_ecm

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=14     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=2      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=mradi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-4 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    ecmmax=0.25d0
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=1
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.25 
        node(i)%hoo=1d1
        node(i)%erp=5d-1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    ! Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
	node(i)%add=node(i)%eqd*1.40; 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters


    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+4)=1d1 
      gen(1)%e(nparam_per_node+5)=1d1
      gen(1)%e(nparam_per_node+6)=1d2
      gen(1)%e(nparam_per_node+7)=0.15d0


    !Adhesion molecules

    ntipusadh=5
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecule interactions
!        kadh(1,1)=1d1
!        kadh(1,2)=1d1
!        kadh(1,3)=1d1
!        kadh(2,1)=1d1
        kadh(2,2)=1d1
        kadh(2,3)=1d1
!        kadh(3,1)=1d1
        kadh(3,2)=1d1
        kadh(3,3)=1d1

    end if

    gen(1)%e(1)=1
    gen(2)%e(1)=2
    gen(3)%e(1)=3
    gen(4)%e(1)=4
    gen(5)%e(1)=5

    !Gene expression on nodes
      do i=1,nd
        if(node(i)%tipus==3)then
          gex(i,1)=1.0d-2
          gex(i,2)=1d-2
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine mes_ecm

!*************************************************************************************************

subroutine neg_adh

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=10     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=mradi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-4 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.3d0
    deltamax=1d-2 ! miguel 14-10-13   
    ecmmax=0.25d0
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d2;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d2;node(i)%rec=1d2	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25 
        node(i)%hoo=1d2
        node(i)%erp=5d-1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d2;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d2;node(i)%rec=1d2	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
		node(i)%add=node(i)%eqd*1.40; 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%ecm=0d0
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters


    !Gene-behavior interactions
    gen(1)%e(1)=1
    gen(2)%e(1)=2
    gen(3)%e(1)=3
    gen(4)%e(1)=4
    gen(5)%e(1)=5

    !Adhesion molecules

      ntipusadh=5
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0

      !Adhesion molecule interactions
      kadh(1,1)=1d1
      kadh(2,2)=1d1
      kadh(1,2)=-1d1
      kadh(2,1)=-1d1
    end if

    node(:)%x=node(:)%x*0.7
    node(:)%y=node(:)%y*0.7

    !Gene expression on nodes
      do i=1,nd
        if(node(i)%x>0.5)then
          gex(i,2)=1d0
        else
          gex(i,1)=1d0
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine


!*******************************************************************************************************

subroutine mes_ecmold
integer::i,j,k,ii,jj,kk,sel !,radi,radicel,layer,mradi,mradicel

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=7     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=3      !number of planar cell layers
	zmes=-1.0d0   !z-position of uppermost layer



    call invagination

    ng=3

	ntipusadh=3
	if (allocated(kadh)) deallocate(kadh)
	allocate(kadh(ntipusadh,ntipusadh))
	kadh=0d0

    kadh(1,:)=(/5d0,1d1,1d1/)
    kadh(2,:)=(/1d1,1d1,5d0/)
    kadh(3,:)=(/1d1,5d0,5d1/)
    gen(1)%e(1)=1           ! Is 29-10-13
    gen(2)%e(1)=2           ! Is 29-10-13
    gen(3)%e(1)=3           ! Is 29-10-13

    call initiate_gene
    gen(1)%e(nparam_per_node+4)=5d-4
    gen(1)%e(nparam_per_node+5)=1d0
    gen(1)%e(nparam_per_node+6)=1d2
    gen(1)%e(nparam_per_node+7)=2*1.25d0

    node(:)%ecm=0.0d0

    do i=1,nd
      if(node(i)%tipus==3)then
        gex(i,1)=1d0
        gex(i,3)=1d0
      elseif(node(i)%icel<8.and.node(i)%tipus==2)then
        gex(i,2)=1.0d0
      end if
    end do

    call update_npag

    prop_noise=0.3d0
    ecmmax=0.25d0
    realtime=0
    reqmin=0.05d0 

    mmae=node(1)%eqd

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************

subroutine epi_polar_growth_conc !here the polarization comes from a signal emmited from cells 2 and 3

  call invagination
  !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
!    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


    do i=1,nd		!veins i parametres mecanics
      node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
      node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
      node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
      node(i)%eqs=0.5d0
      node(i)%add=node(i)%eqd*1.30  
      node(i)%hoo=1d2
      node(i)%erp=1d0
      !node(i)%est=1d1
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      ii=node(i)%icel
      node(i)%dmo=desmax
      node(i)%mov=temp
    end do

    rv=node(1)%add*2.0d0 !miguel 22-5-13
    urv=1d0/rv
    ng=6                !miguel 22-5-13

    call initiate_gene  !miguel 22-5-13

    mmae=node(1)%eqd

    gen(1)%e=0.0d0
    gen(2)%e=0.0d0
    gen(3)%e=0.0d0
    gen(4)%e=0.0d0
    gen(5)%e=0.0d0
    gen(6)%e(nparam_per_node+1)=5.0d3  
    gen(6)%e(nparam_per_node+2)=5.0d-2  
    gen(3)%e(nparam_per_node+8)=1.0d0  
    gen(6)%e(nparam_per_node+9)=1.0d0

    gen(1)%t(1)=1.0d0
    gen(1)%mu=0.1d0
    gen(1)%idiffu=1.0d0
    gen(1)%kindof=0
    gen(1)%nww=2
    gen(1)%r(1,1)=2
    gen(1)%r(1,2)=3
    gen(1)%r(1,3)=0.1d0
    gen(1)%r(2,1)=4
    gen(1)%r(2,2)=5
    gen(1)%r(2,3)=0.1d0


    gen(2)%t(1)=1.0d0
    gen(2)%mu=0.1d0
    gen(2)%idiffu=1.0d0
    gen(2)%kindof=1
    gen(2)%npost=1
    allocate(gen(2)%post(1))
    gen(2)%post(1)=3

    gen(3)%t(2)=1.0d0
    gen(3)%mu=0.1d0
    gen(3)%idiffu=0.0d0
    gen(3)%diffu=0.5d0
    gen(3)%kindof=2
    gen(3)%npost=0
    gen(3)%npre=1
    allocate(gen(3)%pre(1))
    gen(3)%pre(1)=2

    gex(:,4)=1.0d0
    gen(4)%t(4)=1.0d0
    gen(4)%mu=0.1d0
    gen(4)%idiffu=0.0d0
    gen(4)%diffu=0.0d0
    gen(4)%kindof=1
    gen(4)%npost=1
    allocate(gen(4)%post(1))
    gen(4)%post(1)=5
    gen(4)%npre=0

    gen(5)%t(3)=1.0d0
    gen(5)%mu=0.0d0
    gen(5)%idiffu=0.0d0
    gen(5)%diffu=0.0d0
    gen(5)%kindof=2
    gen(5)%npost=0
    gen(5)%npre=1
    allocate(gen(5)%pre(1))
    gen(5)%pre(1)=4
    
    gex(:,6)=1.0d0

    do i=1,cels(2)%nunodes
      ii=cels(2)%node(i)
      gex(ii,1)=1.0d0
      gex(ii,4)=0.0d0
      gex(ii,6)=0.0d0
    end do
    do i=1,cels(3)%nunodes
      ii=cels(3)%node(i)
      gex(ii,1)=1.0d0
      gex(ii,4)=0.0d0
      gex(ii,6)=0.0d0
    end do

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14

    call update_npag    !miguel 22-5-13

    prop_noise=0.1d0    

    realtime=0

    reqmin=0.05d0 
    mmae=node(1)%eqd

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14

    resmax=1d-1  
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!*******************************************************************************************

subroutine epi_polar_growth			!>>Miquel 14-10-12

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.9d0
    min_comp=-1d-1
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1   !!>> HC 15-7-2021 1= multi-node cells
    ffu(5)=0   !!>> HC 13-7-2021 0= Runge-Kutta integration 
    ffu(6)=0   !!>> HC 13-7-2021 0= plasticity
    ffu(9)=0   !!>> HC 13-7-2021 0= NEW VERSION 2021
    ffu(10)=0  !!>> HC 13-7-2021 0= volume conservation
    ffu(15)=0  !!>> HC 13-7-2021 0= unbiased random noise (faster)
    ffu(19)=0  !!>> HC 13-7-2021 0= polarization and polarized division controled by same wa
    ffu(20)=0  !!>> HC 13-7-2021 0= Maximum value for adhesion (maxad)
    ffu(23)=0  !!>> HC 13-7-2021 0= only neighbors are the neighbors in ADD range
    ffu(24)=0  !!>> HC 13-7-2021 0= We recover the lost ADD neighbors
    ffu(25)=0  !!>> HC 13-7-2021 0= Nutrients limited, maximum global cell cycle increase (maxcycl)
    ffu(26)=0  !!>> HC 13-7-2021 0= Gradual implementation of cell properties and behaviors
    ffu(27)=0  !!>> HC 13-7-2021 0= We call to neighbor_build only once in Runge-Kutta
    maxad=0.5

!ffu(8)=1
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=6d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.60  
        node(i)%hoo=1d1
        node(i)%erp=1d0
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1.3d1;node(i)%adh=6d0 	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.40
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=0d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13
      gen(2)%diffu=1d1 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(2)%e(nparam_per_node+1)=5.0d-1
      gen(2)%e(nparam_per_node+2)=3.0d-1  
      gen(1)%e(nparam_per_node+8)=1.0d0  
      gen(2)%e(nparam_per_node+9)=0.5d0   
      gen(2)%e(nparam_per_node+11)=1.0d0   

    !Gene-gene interactions
      gen(2)%t(2)=0.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)
      gex(:,2)=1d0

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine



subroutine epi_polar_growth_single			!>>Miquel 14-10-12

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=1       !number of radial layers of nodes per cell
	radicel=1    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations
    if(radi==1.or.mradi==1) ffu(1)=0


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.5d0
    min_comp=0d0!-1d-1
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(7)=1
    ffu(5)=0
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
    ffu(17)=1 

    ffu(18)=0 


!ffu(8)=1
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=6d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.60  
        node(i)%hoo=1d1
        node(i)%erp=1d0
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1.3d1;node(i)%adh=6d0 	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.40
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=0d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13
      gen(2)%diffu=1d1 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(2)%e(nparam_per_node+1)=5.0d-1
      gen(2)%e(nparam_per_node+2)=1.0d-1  
      !gen(1)%e(nparam_per_node+8)=1.0d0  
      !gen(2)%e(nparam_per_node+9)=0.5d0   
      !gen(2)%e(nparam_per_node+11)=1.0d0   

    !Gene-gene interactions
      gen(2)%t(2)=0.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)
      gex(:,2)=1d0

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine




!***********************************************************************

subroutine differential_growth_mes


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=10     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=4      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters


   print*,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ncels",ncels,"ncals",ncals

    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

  print*,"cels",size(cels)

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3   
    prop_noise=0.0d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(13)=1
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
		node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
		node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
		node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0 !only for epithelium
		node(i)%add=node(i)%eqd*1.3 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=mradi*2 !>>> Is 5-2-14
      end do
    end if

!******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters
    gen(1)%t(1)=0.1d0
    gen(1)%kindof=0 ; gen(1)%diffu=0.5d1 ; gen(1)%mu=1d-1
    gen(2)%kindof=0 ; gen(2)%diffu=0.5d1 ; gen(1)%mu=0d0


    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+1)=1d-1 !0!d2
      gen(1)%e(nparam_per_node+2)=1d-4
      gen(2)%e(nparam_per_node+1)=0
      gen(2)%e(nparam_per_node+2)=1d-2
      gen(3)%e(1)=1
      gen(4)%e(1)=2

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

        kadh(1,1)=1d0
        kadh(1,2)=0.5d0 ; kadh(2,1)=kadh(1,2)
        kadh(2,2)=1d0


    end if

    !Gene expression on nodes

      gex(:nd/2,1)=1.0d0
      gex(:nd/2,3)=1.0d0

      gex(nd/2+1:nd,2)=1.0d0
      gex(nd/2+1:nd,4)=1.0d0


    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0

    min_comp=-1d0

end subroutine

!***********************************************************************

subroutine differential_growth_epi_mes

  !epithelium's dimension parameters
  radi=2       !number of radial layers of nodes per cell
  radicel=2    !number of radial layers of cells
  zepi=0d0     !z-position of basal layer

  !mesenchyme's dimension parameters
  mradi=25     !number of nodes per cell
  mradicel=2   !number of radial cell layers
  layer=1      !number of planar cell layers
  zmes=-1.0d0   !z-position of uppermost layer

  call epi_mes

  do i=1,nd
    if (node(i)%tipus<3) then
      gex(i,1)=1.0d0
    else
      gex(i,2)=1.0d0
    end if
  end do

  kadh(1,2)=1.0d8
  kadh(2,1)=1.0d8
  gen(1)%e=0.0d0
  gen(2)%e=0.0d0
  gen(1)%e(1)=1.d0
  gen(2)%e(1)=1.d0
  gen(1)%e(nparam_per_node+1)=1.0d-3  
  gen(1)%e(nparam_per_node+2)=1.0d-3  
  gen(2)%e(nparam_per_node+1)=1.0d-7  
  gen(2)%e(nparam_per_node+2)=1.0d-7  
  gen(1)%diffu=0.5d0 ; gen(1)%kindof=0  ! miguel 14-10-13   
  gen(2)%diffu=0.5d0 ; gen(2)%kindof=0  ! miguel 14-10-13    

  gen(1)%e(1)=1           ! Is 29-10-13
  gen(2)%e(1)=2           ! Is 29-10-13
  call update_npag

  node(:)%dmo=desmax*4
  node(:)%mov=temp*3

  do i=1,nd		!veins i parametres mecanics
    if (node(i)%tipus<3) cycle
    node(i)%erp=1d2
    !node(i)%est=1d2
    node(i)%dmo=desmax*4
    node(i)%mov=temp*2
  end do

!  do i=15,15!26
!    do j=1,cels(i)%nunodes
!      k=cels(i)%node(j)
!      node(k)%eqd=0.1d-8
!      node(k)%add=0.1d-8
!      gex(k,:)=0.0d0
!    end do
!  end do

!  node(:)%z=node(:)%z*0.5d0

  resmax=3d-1
  deltamax=0.1d0*rv
  prop_noise=0.4
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************

subroutine differential_growth_epi




!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=3    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used

    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%rec=1.5d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.30  
        node(i)%hoo=1d1
        node(i)%erp=1d1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d2  !only for epithelium
        !node(i)%est=1d2 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%kindof=0 ; gen(1)%diffu=0.5d1 ; gen(1)%mu=0d0 ! miguel 14-10-13
      gen(2)%kindof=0 ; gen(2)%diffu=0.5d1 ; gen(2)%mu=0d0 ! miguel 14-10-13
      gen(4)%kindof=0 ; gen(3)%diffu=0.5d1 ; gen(3)%mu=0d0 ! miguel 14-10-13
      gen(3)%kindof=0 ; gen(4)%diffu=0.5d1 ; gen(4)%mu=0d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(1)%e(nparam_per_node+1)=3.0d-1  
      gen(1)%e(nparam_per_node+2)=3.0d-4  
      gen(2)%e(nparam_per_node+1)=1.0d-4  
      gen(2)%e(nparam_per_node+2)=1.0d-4 
      gen(3)%e(1)=1
      gen(4)%e(1)=2

    !Gene-gene interactions

         !****a rellenar****

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions


          kadh(1,1)=1.0d0
          kadh(1,2)=0.5d0
          kadh(2,1)=0.5d0
          kadh(2,2)=1.0d0


    end if

    !Gene expression on nodes
      do i=1,ncels/2
        do j=1,cels(i)%nunodes
          k=cels(i)%node(j)
          gex(k,1)=1d0
          gex(k,3)=1d0
        end do
      end do
      do i=ncels/2+1,ncels
        do j=1,cels(i)%nunodes
          k=cels(i)%node(j)
          gex(k,2)=1d0
          gex(k,4)=1d0
        end do
      end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!****************************************************************************************

subroutine directed_mit_mes

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=8     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

call mes_polar_growth


  !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(2)%diffu=5d1 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0 ! miguel 14-10-13
      gen(1)%diffu=1d1 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions
    gen(2)%e(nparam_per_node+11)=1d0! miguel 14-10-13 (dependance of physical or chemical polarization)
    gen(2)%e(nparam_per_node+1)=1d-1! miguel 14-10-13 (growth)
    gen(2)%e(nparam_per_node+2)=1d-1! miguel 14-10-13 (cell cycle)    
    gen(1)%e(nparam_per_node+12)=1d4  ! miguel 14-10-13 ! assymetric division
    gen(1)%e(nparam_per_node+8)=1d0  ! miguel 14-10-13 ! cell polarization


    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)
      gex(:,2)=1d0                     


    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!******************************************************************************************

subroutine mes_shell ! miguel4-11-13 (subtitution of the previous mes_shell)

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=99     !number of nodes per cell
	mradicel=1   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-1 !low value is low temperature	
    desmax=0.05
    resmax=1d-2   
    prop_noise=0.2d0
    reqmin=0.05d0 
    deltamax=1d-2
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=0d0;node(i)%adh=0d0    	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d3;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=1d4;node(i)%rec=1d5	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.5d0 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
!******* eggshell initialization 
call shellinicial

!******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
    gex(:,1)=1d-2 
    
    !Gene-behavior interactions
    gen(1)%diffu=1d-5 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13
    gen(1)%e(nparam_per_node+1)=1d-6! miguel 14-10-13 (growth)
    gen(1)%e(nparam_per_node+2)=1d-6! miguel 14-10-13 (cell cycle)    
    
    !Gene-behavior interactions


    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

         !****a rellenar****

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine mes_shell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************

subroutine grn_test
  call epi_polar_growth
  gen(1)%e(:)=0
  gen(2)%e(:)=0
!  gen(1)%mu=0.1d0
  gen(2)%mu=0.1d0
  gex=0.0d0
  gex(40:60,1)=1.0d0
  gex(40:50,2)=1.0d0
  gen(1)%t(1)=1d0
  gen(1)%t(2)=-1d0
  call update_npag    !miguel 22-5-13
  dmax=2
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!****************************************************************************

subroutine reaction_test
!  call epi_polar_growth
  call invagination
  ng=3
  call initiate_gene  !Is 15-5-13 It has to be called here because we want to give values to g and w
  node(:)%you=1.2d4
  gen(1)%e=0.0d0
  gen(2)%e=0.0d0
!  gen(3)%e=0.0d0
  ffu(5)=0
!  gen(1)%t(1)=1.1d0
  gen(:)%mu=0.01d0
  gex=0.0d0
  gex(nd,1)=1.0d0
!  gex(:nd,2)=1.0d0
  gen(:)%diffu=0.0d0
  gen(1)%diffu=0.1d0
  gen(:)%idiffu=0.1d0
  gen(1)%kindof=0
  gen(1)%nww=1
  gen(1)%r(1,1)=2
  gen(1)%r(1,2)=3
  gen(1)%r(1,3)=0.1d0
  gen(2)%kindof=1
  gen(3)%kindof=2
  gen(:)%npost=0
  gen(:)%npre=0 
  gen(2)%npost=1
  allocate(gen(2)%post(1)) 
  gen(2)%post(1)=3
  gen(3)%npost=0
  gen(3)%npre=1
  allocate(gen(3)%pre(1))
  gen(3)%pre(1)=2
  gen(1)%t(1)=10.1d0
  gen(2)%t(1)=10.1d0
  gen(3)%t(2)=10.1d0
  deltamax=1d-2
  call update_npag    !miguel 22-5-13
  dmax=2
  screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!************************************************************************************************************


subroutine epi_reaction_diffusion



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=-1.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.01d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1;node(i)%adh=8d1    	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d5	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0;	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25; 
        node(i)%hoo=1d2
        node(i)%erp=1d1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=9
    call initiate_gene


    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=1d0 ; gen(1)%mu=1d-2 ; gen(1)%post=2  !the activator intracell form
       gen(2)%kindof=4 ; gen(2)%diffu=1d0 ; gen(2)%mu=1d-2 ; gen(2)%pre=1   !the activator extracell form
       gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d-2 ; gen(3)%post=4  !the activator receptor inactive form
       gen(4)%kindof=2 ; gen(4)%diffu=1d0 ; gen(4)%mu=1d-2 ; gen(4)%pre=3   !the activator receptor active form
       gen(5)%kindof=1 ; gen(5)%diffu=1d0 ; gen(5)%mu=1d-2 ; gen(5)%post=6  !the inhibitor intracell form
       gen(6)%kindof=4 ; gen(6)%diffu=5d0 ; gen(6)%mu=1d-2 ; gen(6)%pre=5   !the inhibitor extracell form
       gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=1d-2 ; gen(7)%post=8  !the inhibitor receptor inactive form
       gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=1d-2 ; gen(8)%pre=7   !the inhibitor receptor active form
       gen(9)%kindof=0 ; gen(9)%diffu=1d0 ; gen(9)%mu=0d0                   !house-keeping gene, induces expression of receptors

    !Gene-behavior interactions

       !****a rellenar****

    !Gene-gene interactions
      gen(1)%t(2)=1.0d0  !activator induces own secretion
      gen(2)%t(4)=1.0d0  !activator extracell activates activator receptor
      gen(4)%t(1)=5.0d0  !activator active receptor activates transcription of activator
      gen(4)%t(5)=1.0d0  !activator active receptor activates transcription of inhibitor
      gen(5)%t(6)=1.0d0  !inhibitor induces own secretion
      gen(6)%t(8)=1.0d0  !inhibitor extracell activates inhibitor receptor
      gen(8)%t(1)=-5.0d0 !inhibitor active receptor inhibits activator transcription
      gen(9)%t(3)=1.0d0  !house-keeping induces expression of the activator receptor
      gen(9)%t(7)=1.0d0  !house-keeping induces expression of the inhibitor receptor


    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

      !one cell expressing activator in the centre
      gex(1,1)=1.0d0
      gex(:,3)=1.0d0
      gex(:,7)=1.0d0
      gex(:,9)=1.0d0

    
    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!***********************************************************************************

subroutine epi_active_transport


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=-1.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.01d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1;node(i)%adh=8d1    	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d5	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0;	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25; 
        node(i)%hoo=1d2
        node(i)%erp=1d1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene


    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=5 ; gen(1)%diffu=1d0 ; gen(1)%mu=0d0 !; gen(1)%post=2  !the activator intracell form


    !Gene-behavior interactions

       !****a rellenar****

    !Gene-gene interactions
 !     gen(1)%t(1)=1.0d0  !activator induces own secretion


    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

      !one cell expressing activator in the centre
      where(node(:nd)%marge==0) gex(:nd,1)=10d0
    
    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine


!************************************************************************************************************

subroutine diffusion_test

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=-1.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-5
    prop_noise=0.01d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    ffu(7)=0 !constant delta
    ffu(8)=0 !neighboring by triangulation
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1;node(i)%adh=8d1    	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%rec=1d1	!>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0;	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.25; 
        node(i)%hoo=1d1
        node(i)%erp=1d1
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%rec=0d0	!>>Miquel 8-10-12
        node(i)%eqd=0d0 !;node(i)%codel=0d0	!>>Miquel 26-10-12!de
        node(i)%eqs=0d0
        node(i)%add=0d0
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
       gen(:)%kindof=1 ; gen(:)%diffu=1d-1 ; gen(:)%mu=0.1d0  

      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=2 ; gen(1)%diffu=1d1 ; gen(1)%mu=0d0 ; gen(1)%npost=1; 
       allocate(gen(1)%post(1))
       gen(1)%post(1)=2
       
       gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=0d0 ; gen(2)%npre=1 ; 
       allocate(gen(2)%pre(1))
       gen(2)%pre(1)=1


    !Gene-behavior interactions

       !****a rellenar****

    !Gene-gene interactions
      gen(1)%t(1)=1.0d1
      gen(1)%t(2)=0.0d0
      gen(2)%t(1)=1.0d1
!      do i=3,ng
!        gen(i)%t(i)=1.0d1
!      end do
      gen(1)%nww=1
      !gen(1)%r(1,1)=1
      gen(1)%r(1,1)=1
      gen(1)%r(1,2)=2
      gen(1)%r(1,3)=1d1
      

    !Adhesion molecules

    ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(nd-1,1)=1.0d0
    !  gex(nd-1,2)=1.0d0
    
    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!*************************************************************************************

subroutine epi_mes_primordium



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
	mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
	mradicel=4   !number of radial cell layers
	layer=2      !number of planar cell layers
	zmes=-1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-2
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d6

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.30; 
        node(i)%hoo=1d1
        node(i)%erp=5d0
        !node(i)%est=5d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1 ; node(i)%adh=1d1
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.15 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%fix=1
        node(i)%rep=1d5;node(i)%rec=1d5
        node(i)%hoo=1d5
        node(i)%erp=1d5
        !node(i)%est=1d5
      end do

!      j=0
!      do i=1,mradicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+ndmes
        node(i)%fix=1
        node(i)%rep=1d5;node(i)%rec=1d5
!        node(i)%eqd=0.30 ; node(i)%add=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do
    !end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=10
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=2 ; gen(1)%diffu=1d2 ; gen(1)%mu=5d-1 !activator transcript
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
      gen(1)%post(1)=2

      gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=5d-1 !activator signal
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=1

      gen(3)%kindof=2 ; gen(3)%diffu=1d2 ; gen(3)%mu=0d0 !free receptor
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=3 ; gen(4)%diffu=1d2 ; gen(4)%mu=5d-1 !activated receptor
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=3 ; gen(4)%pre(2)=2

      gen(5)%kindof=1 ; gen(5)%diffu=1d2 ; gen(5)%mu=0d0 !TF for activator (inactive)
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=3 ; gen(6)%diffu=1d2 ; gen(6)%mu=5d-1 !activated TF
      gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=5

!      gen(7)%kindof=1 ; gen(7)%diffu=1d1 ; gen(7)%mu=5d-1 !cell proliferation factor

      gen(7)%kindof=1 ; gen(7)%diffu=1d2 ; gen(7)%mu=0d0 !inactive cell proliferation factor (mesenchyme)
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=8

      gen(8)%kindof=3 ; gen(8)%diffu=1d2 ; gen(8)%mu=5d-1 !active cell proliferation factor (mesenchyme)
      gen(8)%npre=1 ; allocate(gen(8)%pre(gen(8)%npre))
      gen(8)%pre(1)=7

      gen(9)%kindof=1 ; gen(9)%diffu=1d2 ; gen(9)%mu=0d0 !inactive cell proliferation factor (epithelium)
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=10

      gen(10)%kindof=3 ; gen(10)%diffu=1d2 ; gen(10)%mu=5d-1 !active cell proliferation factor (epithelium)
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=9


    !Gene-behavior interactions

       gen(8)%e(nparam_per_node+1)=1d-3  !growth
       gen(8)%e(nparam_per_node+2)=1d-3  !cell division

       gen(10)%e(nparam_per_node+1)=1d-5  !growth
       gen(10)%e(nparam_per_node+2)=1d-5  !cell division


    !Gene-gene interactions
       gen(1)%nww=1
       gen(1)%r(1,2)=1d1 !1 autocatalyzes transcription
       gen(2)%nww=1
       gen(2)%r(3,4)=1d1 !signal activates receptor
       gen(4)%nww=3
       gen(4)%r(5,6)=1d1 !active receptor activates TF
       gen(4)%r(7,8)=1d1 !active receptor activates cell prolif. factor
       gen(4)%r(9,10)=1d1 !active receptor activates cell prolif. factor
       gen(6)%t(1)=1d1 !active TF promotes transcirption


    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

    do i=1,nd
      gex(i,3)=1d3
      if(node(i)%tipus==3)then
        gex(i,5)=1d3
        gex(i,7)=1d3
        if(node(i)%marge==0)then
          gex(i,1)=1d0
        end if
      else
        gex(i,9)=1d3
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine


!*************************************************************************************

subroutine mes_primordium
integer:: val


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
        packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
	mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
	mradicel=4   !number of radial cell layers
	layer=3      !number of planar cell layers
	zmes=0.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi ; print*,"nodecel",nodecel	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells      
      val=(6*j+1)*nodecel ;print*,"val",val    !number of nodes per layer
      ndmes=nodecel*ncelsmes  !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        if(packed==0)then
          nodecel=mradi
        else
    !      nodecel=(6*j+1)*mradi
        end if
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10
    !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes,"nd",nd,"ncels",ncels


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-2
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d2

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=5d0 ; node(i)%rec=5d0
        else                      !apical
          node(i)%you=5d0 ; node(i)%adh=5d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25d0
        node(i)%add=node(i)%eqd*1.30; 
        node(i)%hoo=1d2
        node(i)%erp=1d0
        !node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0 	
        node(i)%rep=1d2;node(i)%rec=1d2	
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.15 
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
      do i=1,val !this is the "epithelial like part"
        node(i)%add=node(i)%eqd*1.35
      end do
    end if
    
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !let's do it for the most external layer of cells
      !j=0
      !do i=1,radicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      !do i=j+1,ndepi
      !  node(i)%fix=1
      !end do

      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of mesenchymal nodes wich are not external

      do i=ndepi+j+1,ndepi+val
        node(i)%fix=1
        node(i)%rep=1d6;node(i)%rec=1d6
!        node(i)%eqd=0.30 ; node(i)%add=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do
      do i=ndepi+val+j+1,ndepi+ndmes
        node(i)%fix=1
        node(i)%rep=1d6;node(i)%rec=1d6
!        node(i)%eqd=0.30 ; node(i)%add=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do


    !end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=12
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=2 ; gen(1)%diffu=1d1 ; gen(1)%mu=5d-1 !activator transcript
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
      gen(1)%post(1)=2

      gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=5d-1 !activator signal
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=1

      gen(3)%kindof=2 ; gen(3)%diffu=1d1 ; gen(3)%mu=0d0 !free receptor
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=3 ; gen(4)%diffu=1d1 ; gen(4)%mu=5d-1 !activated receptor
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=3 ; gen(4)%pre(2)=2

      gen(5)%kindof=1 ; gen(5)%diffu=1d1 ; gen(5)%mu=0d0 !TF for activator (inactive)
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=3 ; gen(6)%diffu=1d1 ; gen(6)%mu=5d-1 !activated TF
      gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d1 ; gen(7)%mu=0d0 !inactive cell proliferation factor (mesenchyme)
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=8

      gen(8)%kindof=3 ; gen(8)%diffu=1d1 ; gen(8)%mu=5d-1 !active cell proliferation factor (mesenchyme)
      gen(8)%npre=1 ; allocate(gen(8)%pre(gen(8)%npre))
      gen(8)%pre(1)=7

      gen(9)%kindof=1 ; gen(9)%diffu=1d1 ; gen(9)%mu=0d0 !inactive cell proliferation factor (epithelium)
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=10

      gen(10)%kindof=3 ; gen(10)%diffu=1d1 ; gen(10)%mu=5d-1 !active cell proliferation factor (epithelium)
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=9

      gen(11)%kindof=1 ; gen(11)%diffu=1d1 ; gen(11)%mu=0d0 !epithelial adhesion molecule
      
      gen(12)%kindof=1 ; gen(12)%diffu=1d1 ; gen(12)%mu=0d0 !mesenchymal adhesion molecule

    !Gene-behavior interactions

       gen(8)%e(nparam_per_node+1)=1d-2  !growth
       gen(8)%e(nparam_per_node+2)=5d-2  !cell division

       gen(10)%e(nparam_per_node+1)=5d-4  !growth
       gen(10)%e(nparam_per_node+2)=1d-4  !cell division
       
       gen(11)%e(1)=1
       gen(12)%e(1)=2


    !Gene-gene interactions

       gen(1)%nww=1
       gen(1)%r(1,2)=1d1 !1 autocatalyzes transcription
       gen(2)%nww=1
       gen(2)%r(3,4)=1d1 !signal activates receptor
       gen(4)%nww=3
       gen(4)%r(5,6)=1d1 !active receptor activates TF
       gen(4)%r(7,8)=1d1 !active receptor activates cell prolif. factor
       gen(4)%r(9,10)=1d1 !active receptor activates cell prolif. factor
       gen(6)%t(1)=1d1 !active TF promotes transcirption



    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=1d1
      kadh(1,2)=-1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d1


    end if

    !Gene expression on nodes

    do i=1,nd
      gex(i,3)=1d3
    !  if(node(i)%tipus==3)then
    !    gex(i,5)=1d3
    !    gex(i,7)=1d3
    !    if(node(i)%marge==0)then
    !      gex(i,1)=1d0
    !    end if
    !  else
    !    gex(i,9)=1d3
    !  end if
    end do
    
    do i=1,val
      gex(i,9)=1d3
      gex(i,11)=1d0
    end do
    do i=val+1,nd
      gex(i,5)=1d3
      gex(i,7)=1d3
      gex(i,12)=1d0
      if(node(i)%marge==0) gex(i,1)=1d0
    end do

    call update_npag

    node(:)%talone=0.0d0

    !print*,"w de 1",gen(1)%t(:)
    !print*,"w de 2",gen(2)%t(:)
    !print*,"w de 3",gen(3)%t(:)
    !print*,"w de 4",gen(4)%t(:)
    !print*,"w de 5",gen(5)%t(:)
    !print*,"w de 6",gen(6)%t(:)
    !print*,"w de 7",gen(7)%t(:)
    !print*,"w de 8",gen(8)%t(:)
    !print*,"w de 9",gen(9)%t(:)
    !print*,"w de 10",gen(10)%t(:)
    !print*,"w de 11",gen(11)%t(:)
    !print*,"w de 12",gen(12)%t(:)
    ramax=maxval(node(:)%add)*3    
    node(:)%dif=0.0d0
end subroutine

!!!!!!!!!***********************************SUBROUTINE********************************

subroutine epi_mes_bud



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=2       !number of radial layers of nodes per cell
    radicel=4    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=-1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d0


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=5d0 ; node(i)%rec=5d0
        else                      !apical
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=5d0 ; node(i)%rec=5d0
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.40; 
        node(i)%hoo=1d2
        node(i)%erp=1d0
        !node(i)%est=1d0
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%pla=1d-2 ; node(i)%kvol=3d-5
        node(i)%kfi=khold
        node(i)%dif=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d0 ; node(i)%adh=0d0
        node(i)%rep=1d0 ; node(i)%rec=5d0
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if


!    do i=2,7
!      a=node(1)%x-node(i)%x ; b=node(1)%y-node(i)%y ; c=node(1)%z-node(i)%z
!      d=sqrt(a**2+b**2+c**2) ; a=a*0.15/d ; b=b*0.15/d ; c=c*0.15/d
!      node(i)%x=node(i)%x+a ; node(i)%y=node(i)%y+b ; node(i)%z=node(i)%z+c
!    end do

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
!    node(2:7)%fix=1 !; node(2:7)%rep=1d3 ; node(2:7)%you=1d3
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%fix=1 ;node(i)%kfi=khold
        node(i)%oriz=node(i)%oriz-20
        if(node(i)%tipus==2)then;node(i)%orix=0d0 ; node(i)%oriy=0d0;end if
        !node(i)%rep=1d1;node(i)%rec=1d1
        !node(i)%hoo=1d1
        !node(i)%erp=1d1
        !!node(i)%est=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      !j=(6*j+1)*mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndmes
      !  node(i)%fix=1
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
    
    
!    do i=ndepi+1,nd ; node(i)%fix=1;enddo
    
    
    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=0d0 !activator transcript
      
      gen(2)%kindof=1 ; gen(2)%diffu=1d1 ; gen(2)%mu=0d0 !epithelial basal adhesion molecule
      gen(3)%kindof=1 ; gen(3)%diffu=1d1 ; gen(3)%mu=0d0 !epithelial apical adhesion molecule
      gen(4)%kindof=1 ; gen(4)%diffu=1d1 ; gen(4)%mu=0d0 !mesenchymal adhesion molecule

      gen(5)%kindof=1 ; gen(5)%diffu=1d1 ; gen(5)%mu=0d0 !activator transcript



    !Gene-behavior interactions
 
       gen(2)%e(1)=1  !epithelial adhesion molecule
       gen(3)%e(1)=2  !epithelial-mesenchymal adhesion molecule
       gen(4)%e(1)=3 !mesenchymal adhesion molecule

    !Gene-gene interactions



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d0
      kadh(1,2)=5d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0
      kadh(1,3)=5d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,2)=5d0
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d1


    end if

    !Gene expression on nodes

    do i=1,nd
      !gex(i,3)=1d3 !;gex(i,10)=1d-5
      if(node(i)%tipus==3)then
        !gex(i,5)=1d3
        !gex(i,7)=1d3
        gex(i,4)=1d0 ;gex(i,5)=1d0
        !if(node(i)%marge==0)then
        !  gex(i,1)=1d0
        !end if
         gex(i,1)=1d0
      else
        !gex(i,9)=1d3
        gex(i,1)=1d0
        if(node(i)%tipus==2)then
          gex(i,2)=1d0
        else
          gex(i,3)=1d0
        end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3

    
    
end subroutine


!!!!!!!!!***********************************SUBROUTINE********************************

subroutine hair_placode



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=0    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=9   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.0d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 ;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.5d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.7d0
    khold=1d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    

    ffu(6)=0 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(15)=0 !0 = unbiased random noise / 1 = noise biased by energies

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.70; 
        node(i)%hoo=1d2
        node(i)%erp=5d1
        !node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%rec=1d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if



    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
    !  j=0
    !  do i=1,radicel-2
    !    j=j+i
    !  end do
    !  j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    !  do i=j+1,ndepi
    !    node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !!node(i)%est=1d1
    !  end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
        node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=6
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !E-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 !P-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !basal lamina

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 !migration gene

      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0 !morphogen transcript
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=4 ; gen(6)%diffu=1d0 ; gen(6)%mu=5d-1 !morphogen diffusible form
      gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=5

      
    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !epithelial adhesion molecule
       gen(2)%e(1)=2  !epithelial-mesenchymal adhesion molecule
       gen(3)%e(1)=3  !basal lamina
       gen(4)%e(5)=-0.05 !this is req (emulating a migratory behavior)
       gen(4)%e(6)=0.10 !this is da (emulating a migratory behavior)
       gen(4)%e(16)=0.01 !this is dmo (migratory cells)
       gen(6)%e(nparam_per_node+8)=1d0
       gen(6)%e(nparam_per_node+16)=1d0 !this makes random noise biased towards the gradient of the gene



    !Gene-gene interactions

      !gen(5)%t(5)=1.0d3
      !gen(5)%nww=1
      !gen(5)%r(1,1)=5
      !gen(5)%r(1,2)=6
      !gen(5)%r(1,3)=1d2

    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=1d1
      kadh(1,2)=2d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d1
      kadh(1,3)=-1d1 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=-1d1 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d1
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-3
      j=j+ii
    end do
    j=(6*j+1) ;print*,"jota",j

    do i=1,nd
      !a=sqrt(node(i)%x**2+node(i)%y**2)
      if(i<=ndepi)then
        if(node(i)%icel<=j)then
          gex(i,2)=1d0!(3.5-a)/3.5 !; print*,"gex(i,2)",gex(i,2)
          !gex(i,4)=1d0
          gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        else
          gex(i,1)=1d0
        end if
        !if(node(i)%tipus==2)then
        !  gex(i,3)=1d0
        !  gex(i,1)=0d0 ; gex(i,2)=0d0
        !end if
      else
        j=0
        do ii=1,mradicel-3
          j=j+ii
        end do
        j=(6*j+1) !territory 1 (inner ring)
        k=0
        do ii=1,mradicel-1
          k=k+ii
        end do
        k=(6*k+1) !whole layer
        if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
          gex(i,2)=1d0!(3.5-a)/3.5
          gex(i,4)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
          gex(i,1)=1d0
        elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
          gex(i,2)=1d0
          gex(i,4)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
          gex(i,1)=1d0
        !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
        end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3
   
end subroutine

!************************************************************************************


subroutine feather_placode



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=0
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=1d0
    angletor=0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation (epithelial)
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60; 
        node(i)%hoo=5d1
        node(i)%erp=5d0
        node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=5d-2
        node(i)%kvol=5d-2
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=3d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%fix=1 ;node(i)%rec=1d2
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%fix=1;node(i)%rec=1d2
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      do i=ndepi+k+1,nd
        node(i)%fix=1;node(i)%rec=1d2
      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=14
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=1d0 !Epithelial-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=1d0 !Mesenchymal-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d0 !housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 !epi. morphogen transcript (like SHH)
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=5

      gen(5)%kindof=4 ; gen(5)%diffu=1d0 ; gen(5)%mu=1d0 !epi. morphogen diffusible form (like SHH)
      gen(5)%npre=2 ; allocate(gen(5)%pre(gen(5)%npre))
      gen(5)%pre(1)=4 ; gen(5)%pre(2)=9
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=9

      gen(6)%kindof=2 ; gen(6)%diffu=0d0 ; gen(6)%mu=0d0 !mes. morphogen transcript (like BMP)
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=7

      gen(7)%kindof=4 ; gen(7)%diffu=1d0 ; gen(7)%mu=5d-1 !mes. morphogen diffusible form (like BMP)
      gen(7)%npre=2 ; allocate(gen(7)%pre(gen(7)%npre))
      gen(7)%pre(1)=6 ; gen(7)%pre(2)=11
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=11
      
      gen(8)%kindof=2 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d0 !SHH receptor
      gen(8)%npost=1 ; allocate(gen(8)%post(gen(8)%npost))
      gen(8)%post(1)=9
      
      gen(9)%kindof=8 ; gen(9)%diffu=0d0 ; gen(9)%mu=0d0 !SHH receptor activated
      gen(9)%npre=2 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=8 ; gen(9)%pre(2)=5
      gen(9)%npost=2 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=8 ; gen(9)%post(2)=5
      
      gen(10)%kindof=2 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d0 !BMP receptor
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=0d0 !BMP receptor activated
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=10 ; gen(11)%pre(2)=5
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=10 ; gen(11)%post(2)=5

      gen(12)%kindof=1 ; gen(12)%diffu=0d0 ; gen(12)%mu=1d0 !housekeeping gene mesenchymal
      gen(13)%kindof=1 ; gen(13)%diffu=1d0 ; gen(13)%mu=1d0 !housekeeping gene signalling center
      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=0d0 !adhesion molecule hold
      
      
    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !epithelial adhesion molecule
       gen(2)%e(1)=2  !epithelial-mesenchymal adhesion molecule
       gen(14)%e(1)=3  !basal lamina
       !gen(4)%e(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%e(6)=0.10 !this is da (emulating a migratory behavior)
       !gen(4)%e(16)=0.03 !this is dmo (migratory cells)
       !gen(6)%e(nparam_per_node+8)=1d0
       !gen(6)%e(nparam_per_node+16)=1d-2 !this makes random noise biased towards the gradient of the gene
       gen(9)%e(nparam_per_node+2)=1d-2 !SHH effect on epithelial growth
       gen(11)%e(nparam_per_node+2)=5d-3 !BMP effect on mesenchymal growth



    !Gene-gene interactions

      !gen(11)%nww=2      !BMP activated receptor induces production of BMP
      !gen(11)%r(1,1)=6  
      !gen(11)%r(1,2)=7
      !gen(11)%r(1,3)=1d1
      !gen(11)%r(2,1)=4   !BMP activated receptor induces production of SHH
      !gen(11)%r(2,2)=5
      !gen(11)%r(2,3)=1d1
      
      !gen(7)%nww=1      !BMP morphogen activates BMP receptor
      !gen(7)%r(1,1)=10  
      !gen(7)%r(1,2)=11
      !gen(7)%r(1,3)=1d0

      gen(3)%t(3)=1d0 !housekeeping activates itself
      gen(1)%t(3)=1d0  !activates epi cadherin
      gen(8)%t(3)=1d0  !activates ssh receptor

      gen(12)%t(12)=1d0 !housekeeping activates itself
      gen(2)%t(12)=1d0   !activates mes cadherin
      gen(10)%t(12)=1d0  !activates BMP receptor

      gen(13)%t(13)=1d0 !housekeeping activates itself
      gen(1)%t(13)=1d0  !activates epi cadherin
      gen(4)%t(13)=1d1  !activates ssh transcript



      gen(13)%nww=1    !housekeeping epi mediates secretion of shh
      gen(13)%r(1,1)=4
      gen(13)%r(1,2)=5
      gen(13)%r(1,3)=1d1

      
      gen(9)%nww=4      !SHH morphogen activates SHH receptor
      gen(9)%r(1,1)=5  
      gen(9)%r(1,2)=9
      gen(9)%r(1,3)=1d0
      gen(9)%r(2,1)=9  
      gen(9)%r(2,2)=5
      gen(9)%r(2,3)=1d0

      gen(9)%r(3,1)=8  
      gen(9)%r(3,2)=9
      gen(9)%r(3,3)=1d0
      gen(9)%r(4,1)=9  
      gen(9)%r(4,2)=8
      gen(9)%r(4,3)=1d0


      gen(11)%nww=4      !SHH morphogen activates SHH receptor
      gen(11)%r(1,1)=5  
      gen(11)%r(1,2)=11
      gen(11)%r(1,3)=1d0
      gen(11)%r(2,1)=11  
      gen(11)%r(2,2)=5
      gen(11)%r(2,3)=1d0

      gen(11)%r(3,1)=10
      gen(11)%r(3,2)=11
      gen(11)%r(3,3)=1d0
      gen(11)%r(4,1)=11  
      gen(11)%r(4,2)=10
      gen(11)%r(4,3)=1d0



                        !SHH morphogen activates BMP receptor
      !gen(5)%r(3,1)=5  
      !gen(5)%r(3,2)=11
      !gen(5)%r(3,3)=1d0
      !gen(5)%r(4,1)=11  
      !gen(5)%r(4,2)=5
      !gen(5)%r(4,3)=1d0

      !gen(8)%nww=2      !ssh receptor activates SHH receptor
      !gen(8)%r(1,1)=8  
      !gen(8)%r(1,2)=9
      !gen(8)%r(1,3)=1d0
      !gen(8)%r(2,1)=9  
      !gen(8)%r(2,2)=8
      !gen(8)%r(2,3)=1d0

      !gen(10)%nww=2               !bmp receptor activates BMP receptor
      !gen(10)%r(1,1)=10  
      !gen(10)%r(1,2)=11
      !gen(10)%r(1,3)=1d0
      !gen(10)%r(2,1)=11  
      !gen(10)%r(2,2)=10
      !gen(10)%r(2,3)=1d0



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=2d1
      kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d2

    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,nd
      if(node(i)%fix==1) gex(i,14)=1d0
      !a=sqrt(node(i)%x**2+node(i)%y**2)
      if(i<=2)then; gex(i,4)=1d0; gex(i,13)=1;end if
      if(i<=ndepi)then
        gex(i,1)=1d0
        if(node(i)%icel>1)then
          gex(i,3)=1d0
          gex(i,8)=1d0
        !else
        !  gex(i,1)=1d0
        end if
      else
        gex(i,12)=1d0
        gex(i,10)=1d0
        j=0
        do ii=1,mradicel-2
          j=j+ii
        end do
        j=(6*j+1) !territory 1 (inner ring)
        k=0
        do ii=1,mradicel-1
          k=k+ii
        end do
        k=(6*k+1) !whole layer
        if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
          gex(i,2)=1d0!(3.5-a)/3.5
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
          !gex(i,1)=1d0
        elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
          gex(i,2)=1d0
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
          !gex(i,1)=1d0
        !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
        end if
      end if
    end do

    


    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    
end subroutine

!************************************************************************************



subroutine epi_mes_bud_ingrowth



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+30
    
    if(radi==1.or.mradi==1) ffu(1)=0
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-4
    prop_noise=0.3d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=5d0
    mnn=700

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=0 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=2d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
          node(i)%erp=3d1
          !node(i)%est=1d1
        else                      !apical
          node(i)%you=2d1 ; node(i)%adh=0d0
          node(i)%rep=2d1 ; node(i)%rec=2d1
          node(i)%erp=1d1
          !node(i)%est=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.40; 
        node(i)%hoo=1d1
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=5d-2
        node(i)%kvol=1d-2
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=2d1 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.60
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if



    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%oriz=node(i)%oriz-100
        !node(i)%rep=1d2;node(i)%rec=1d2
        !node(i)%hoo=1d1
        !node(i)%erp=1d1
        !!node(i)%est=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+ndmes
        node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d2;node(i)%rec=1d2
      end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !E-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 !P-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !basal lamina

      gen(4)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !motility gene

      !gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0 !growth gene epithelium

    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !interplacodal adhesion molecule
       gen(2)%e(1)=2  !placodal adhesion molecule
       gen(3)%e(1)=3  !basal lamina
       !gen(4)%e(nparam_per_node+1)=4d-4 !growth mesenchyme
       !gen(4)%e(nparam_per_node+2)=1d-3 !cell cycle
       !gen(5)%e(nparam_per_node+1)=1d-5 !growth epithelium
       !gen(5)%e(nparam_per_node+2)=1d-3 !cell cycle

       gen(4)%e(5)=0.00 !this is req (emulating a migratory behavior)
       gen(4)%e(6)=0.15 !this is da (emulating a migratory behavior)
       !gen(4)%e(9)=1d1 ; gen(4)%e(10)=1d1 !this is rep (so the nodes don't collapse)
       gen(4)%e(16)=0.01 !this is dmo (migratory cells)       
       

    !Gene-gene interactions



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=8d1
      kadh(1,2)=5d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=8d1
      kadh(1,3)=0d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=2d2


    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) ;print*,"jota",j

    do i=1,nd
      !gex(i,1)=1d0
      if(node(i)%tipus==2) gex(i,3)=1d0
      !if(i<=ndepi)then
        if(node(i)%icel<=j.or.(i>ndepi.and.node(i)%icel<=ncelsepi+j))then
          gex(i,4)=1d0
          gex(i,2)=1d0 !P-cadherin
        !if(node(i)%tipus==2)then
          !gex(i,3)=1d0
        else
          gex(i,1)=1d0 !E-cadherin
        end if
      !else
        !gex(i,1)=1d0
        !if(node(i)%icel<=ncelsepi+j) gex(i,4)=1d0
      !end if

      !if(i<=ndepi)then
      !  if(node(i)%icel<=61)then
      !    gex(i,2)=1d0!(3.5-a)/3.5 !; print*,"gex(i,2)",gex(i,2)
      !    gex(i,4)=1d0
      !  else
      !    gex(i,1)=1d0
      !  end if
      !  !if(node(i)%tipus==2)then
      !  !  gex(i,3)=1d0
      !  !  gex(i,1)=0d0 ; gex(i,2)=0d0
      !  !end if
      !else
      !  j=0
      !  do ii=1,mradicel-2
      !    j=j+ii
      !  end do
      !  j=(6*j+1) !;print*,"jota",j
      !  if(node(i)%icel<=ncelsepi+j)then
      !    gex(i,2)=1d0!(3.5-a)/3.5
      !    gex(i,4)=1d0
      !  elseif(node(i)%icel>ncelsepi+ncelsmes/2 .and.node(i)%icel<=ncelsepi+ncelsmes/2+j)then
      !    gex(i,2)=1d0!(3.5-a)/3.5
      !    gex(i,4)=1d0
      !  else
      !    gex(i,1)=1d0
      !  end if
      !end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
end subroutine


!*************************************************************************************

subroutine ic_emt

  call epi_mes

  ng=1
  call initiate_gene

  do j=1,cels(1)%nunodes
    i=cels(1)%node(j)
    gex(i,1)=1.0d0
    !node(cels(1)%node(j))%tipus=3
  end do
  gen(1)%e(nparam_per_node+13)=1.0d0

  call update_npag

  !cels(1)%ctipus=3
  dmax=2
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!************************************************************************************

subroutine twoeps !is to test what happens with face-to-face epithelial interactions

call invagination

do i=1,ncels
  if (node(cels(i)%node(1))%x<0) then 
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      node(ii)%z=node(ii)%z+1.2
      node(ii)%x=abs(node(ii)%x)
      if (node(ii)%tipus==1) then ; node(ii)%tipus=2 ; else ; node(ii)%tipus=1 ; end if
    end do
  end if
end do

do i=1,ng
  gen(i)%e=0.0
end do

!node(:)%you=1.2d2  ;node(:)%adh=1.0d2         
prop_noise=0.1d0
temp=0.001

!node(:)%erp=0.0d1
!node(:)%est=0.0d3

call iniboxesll
  dmax=2
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!************************************************************************************

subroutine twoeps_rev !is to test what happens with face-to-face epithelial interactions
integer ind

call invagination

do i=1,ncels
  if (node(cels(i)%node(1))%x<0) then 
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      node(ii)%z=node(ii)%z+1.2
      node(ii)%x=abs(node(ii)%x)
!      if (node(ii)%tipus==1) then ; node(ii)%tipus=2 ; else ; node(ii)%tipus=1 ; end if
    end do
  end if
end do

do i=1,ng
  gen(i)%e=0.0
end do

!node(:)%you=1.2d3  ;node(:)%adh=1.0d3         
prop_noise=0.1d0

!node(:)%erp=0.0d1
!node(:)%est=0.0d3
call iniboxesll
  dmax=2
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine

!****************************************************************************************************

subroutine founding_father

  character*200 carga

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=1    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    !End initializing dimension parameters

    !print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-2
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    ndmax=5000

    !biological
    !functions used    
    ffu=0
     !spring of the ellipse
    ffu(2)=1 !to quite if there are too many cells
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=0 !lonely cells and nodes die
    ffu(5)=0
    ffu(8)=0
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1  ;node(i)%adh=1.0d1         	!>>Miquel 26-10-12
        node(i)%rep=1d1    ;node(i)%rec=1d2	        !>>Miquel 8-10-12
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd !;node(i)%codel=0.25d0     	!>>Miquel 26-10-12!de
        node(i)%eqs=0.5d0
        node(i)%add=node(i)%eqd*1.30  
        node(i)%hoo=1d1
        node(i)%erp=1d-1
        !node(i)%est=1d3
        node(i)%mov=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
!    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0

   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=30 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=11
    call initiate_gene
    do i=1,8
      gen(i)%kindof=i
      gen(i)%diffu=1.0d-2
      gen(i)%mu=1.0d-2
    end do
    gen(1)%nww=4
    gen(1)%r(1,1)=2
    gen(1)%r(1,2)=3
    gen(1)%r(1,3)=1.0d0
    gen(1)%r(2,1)=3
    gen(1)%r(2,2)=4
    gen(1)%r(2,3)=1.0d0
    gen(1)%r(3,1)=9
    gen(1)%r(3,2)=5
    gen(1)%r(3,3)=1.0d0
    gen(1)%r(4,1)=9
    gen(1)%r(4,2)=6
    gen(1)%r(4,3)=1.0d0

    gen(2)%npost=1
    allocate(gen(2)%post(1))
    gen(2)%post(1)=3   !1 t-> 2 -> 3 -> 4
    gen(2)%e(1)=1.0d0

    gen(3)%npre=1
    allocate(gen(3)%pre(1))
    gen(3)%pre(1)=2
    gen(3)%npost=1
    allocate(gen(3)%post(1))
    gen(3)%post(1)=4
    gen(4)%nww=1
    gen(4)%r(1,1)=11 !inactive form of the receptor
    gen(4)%r(1,2)=8  !active (bound) form of the receptor
    gen(4)%r(1,3)=1.0d0

    gen(9)%kindof=2
    gen(9)%npost=2
    allocate(gen(9)%post(2))
    gen(9)%post(1)=5
    gen(9)%post(2)=6   !1 t-> 9 -> 5
                       !      9 -> 6
    gen(5)%npre=1
    allocate(gen(5)%pre(1))
    gen(5)%pre(1)=9
    gen(6)%npre=1
    allocate(gen(6)%pre(1))
    gen(6)%pre(1)=9

    gen(10)%kindof=2
    gen(10)%npost=1
    allocate(gen(10)%post(1))
    gen(10)%post(1)=7  !1 t-> 10 -> 7
    gen(7)%npre=1
    allocate(gen(7)%pre(1))
    gen(7)%pre(1)=10
    gen(7)%nww=1
    gen(7)%r(1,1)=10  !this is like homotypic binding
    gen(7)%r(1,2)=7
    gen(7)%r(1,3)=1.0d0 !notch needs to catalyze its own synthesis

    gen(11)%kindof=2
    gen(11)%npost=1
    allocate(gen(11)%post(1))
    gen(11)%post(1)=8  !1 t-> 11 <-> 8
    gen(11)%npre=1            
    allocate(gen(11)%pre(1))
    gen(11)%pre(1)=8  
    gen(8)%npre=1
    allocate(gen(8)%pre(1))
    gen(8)%pre(1)=11
    gen(8)%npost=1
    allocate(gen(8)%post(1))
    gen(8)%post(1)=11        !so that it can revert

    ! the active receptor catalyzes its own inactivation
    gen(8)%nww=1
    gen(8)%r(1,1)=8
    gen(8)%r(1,2)=11
    gen(8)%r(1,3)=0.1d0

    do i=1,ng 
      do j=1,nga
        gen(i)%e(j)=0.0d0
      end do
      do j=1,ng
        gen(i)%t(j)=0.0d0
      end do
    end do

    gen(1)%t(1)=1.0d0  ! 1 activates itself and all the primary forms
    gen(2)%t(1)=1.0d0
    gen(9)%t(1)=1.0d0
    gen(10)%t(1)=1.0d0
    gen(11)%t(1)=1.0d0

    gen(8)%t(8)=-0.1d0

    !Gene expression on nodes
    gex=0.0d0
    gex(:,1)=1.0d0

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
    gen(2)%e(1)=1.0d0   !adhesion molecule
    gen(1)%e(23)=9.0d-4
    gen(2)%e(23)=9.0d-4
    gen(1)%e(nparam_per_node+1)=5d-5

    ntipusadh=1
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
     
         !****a rellenar****
         kadh(1,1)=1d-3
    end if


    call update_npag
    getot=1d4

    call getarg(1,carga)
    if (len_trim(carga)==0.or.carga=="0") then
      if (len(ccag)>0.or.ccag=="0") then
        read(ccag,*,ERR=333) iccag                   !>>> Is 2-3-14
        print *,-iccag,"notice for this ic the random seed idum changes with the date"
        idumoriginal(1)=-iccag               !>>> Is 2-3-14
        call random_seed(put=idumoriginal)   !>>> Is 2-3-14
      end if
    end if
333 continue
end subroutine

!*********************************************************************************************************

subroutine invaginacio_diff

call epi_polar_growth
gen(1)%e(23)=1.0d-4
gen(2)%e(23)=1.0d-4
call update_npag

end subroutine

!**************************************************************************************************************************

subroutine blastuloid
  type(nod) :: cnode
  real*8 vertices(12,3)

  call founding_father

  call icosahedra(vertices)

  cnode=node(1)

  nd=24
  nda=nd+10
  ncels=12
  ncals=ncels+10 

  deallocate(node)
  allocate(node(nda))
  deallocate(nneigh)
  allocate(nneigh(nda))
  if (allocated(neigh)) deallocate(neigh)
  allocate(neigh(nda,mnn))
  if (allocated(dneigh)) deallocate(dneigh)
  allocate(dneigh(nda,mnn))

  node(:nd)=cnode

  do i=1,12
    node(i)%marge=1
    node(i)%tipus=1
    node(i)%altre=i+12
    node(i)%x=vertices(i,1)*0.75d0
    node(i)%y=vertices(i,2)*0.75d0
    node(i)%z=vertices(i,3)*0.75d0
    node(i)%add=node(i)%add*2.8
    node(i)%eqd=node(i)%eqd*2.8
    node(i)%cod=node(i)%cod*2.8
    node(i)%grd=node(i)%grd*2.8
    node(i)%pld=node(i)%pld*2.8
  end do
  do i=13,24
    node(i)%marge=0
    node(i)%tipus=2
    node(i)%altre=i-12
    node(i)%x=vertices(i-12,1)*0.25d0
    node(i)%y=vertices(i-12,2)*0.25d0
    node(i)%z=vertices(i-12,3)*0.25d0
  end do

  do i=1,nd
    ii=node(i)%altre
    node(i)%hoo=sqrt((node(i)%x-node(ii)%x)**2+(node(i)%y-node(ii)%y)**2+(node(i)%z-node(ii)%z)**2)
  end do

  deallocate(gex)
  allocate(gex(nda,ng))
  deallocate(agex)
  allocate(agex(nda,ng))

  if (allocated(px)) deallocate(px)
  if (allocated(py)) deallocate(py)
  if (allocated(pz)) deallocate(pz)
  if (allocated(dex)) deallocate(dex)
  allocate(px(nda),py(nda),pz(nda),dex(nda))
  px=0 ; py=0 ; pz=0 ; dex=0

   
  if (allocated(vcilx)) deallocate(vcilx)
  if (allocated(vcily)) deallocate(vcily)
  if (allocated(vcilz)) deallocate(vcilz)
  allocate(vcilx(nda),vcily(nda),vcilz(nda))
  vcilx=0 ; vcily=0 ; vcilz=0

  if (allocated(vtorx)) deallocate(vtorx)
  if (allocated(vtory)) deallocate(vtory)
  if (allocated(vtorz)) deallocate(vtorz)
  allocate(vtorx(nda),vtory(nda),vtorz(nda))
  vtorx=0 ; vtory=0 ; vtorz=0
  if (allocated(vstorx)) deallocate(vstorx)
  if (allocated(vstory)) deallocate(vstory)
  if (allocated(vstorz)) deallocate(vstorz)
  allocate(vstorx(nda),vstory(nda),vstorz(nda))
  vstorx=0 ; vstory=0 ; vstorz=0
  if (allocated(vsprx)) deallocate(vsprx)
  if (allocated(vspry)) deallocate(vspry)
  if (allocated(vsprz)) deallocate(vsprz)
  allocate(vsprx(nda),vspry(nda),vsprz(nda))
  vsprx=0 ; vspry=0 ; vsprz=0

  if (allocated(erep)) deallocate(erep)
  if (allocated(erepcel)) deallocate(erepcel)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(eyou)) deallocate(eyou)
  if (allocated(espring)) deallocate(espring)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(etor)) deallocate(etor)
  allocate(erep(nda),erepcel(nda),eadh(nda),eyou(nda),espring(nda),etor(nda))

  if (allocated(fmeanl)) deallocate(fmeanl)  !>>Miquel23-1-14
  allocate(fmeanl(nda))                     !>>Miquel23-1-14
  fmeanl=0                                !>>Miquel23-1-14
   
  if (allocated(fmeanv)) deallocate(fmeanv)  !>>Miquel23-1-14
  allocate(fmeanv(nda))                     !>>Miquel23-1-14
  fmeanv=0                                !>>Miquel23-1-14

  do i=1,12
    node(i)%icel=i
  end do

  do i=13,24
    node(i)%icel=node(i)%altre
  end do

  if (allocated(cels)) deallocate(cels)
  allocate(cels(ncals))
  ncels=nd/2
  do i=1,ncels
    cels(i)%fase=0d0
    cels(i)%minsize_for_div=cels(i)%nunodes*2
    cels(i)%maxsize_for_div=30 !>>> Is 5-2-14
    cels(i)%nunodes=2
    cels(i)%nodela=4
    if (allocated(cels(i)%node)) deallocate(cels(i)%node)
    allocate(cels(i)%node(4))
    cels(i)%node(1)=i
    cels(i)%node(2)=i+12
    mmae=node(cels(i)%node(1))%eqd
  end do

  gex(:nd,:)=0.0d0
  gex(:nd,1)=1.0d0


!  ncals=ncels+10 
!  ffu(8)=1
  gen(1)%e(nparam_per_node+2)=1d-3
  call update_npag
  prop_noise=1.0d-1
  
  getot=1d4
  
end subroutine blastuloid

!***************************************************************************************************************

subroutine blastula
  integer vv(12,24),vvv(12),eli(12)
  real*8 vertices(12,3)
  integer :: i,j,k,ii,jj,kk,iii,jjj,kkk,nr(6)
  integer :: margs(12)
  real*8  :: d,a,e  
  real*8  :: npx(12,6),npy(12,6),npz(12,6)
  real*8  :: snpx(12,6),snpy(12,6),snpz(12,6)


  call epi_growth_and_division(2,3)
  deallocate(gex)
  allocate(gex(nda,ng))
  deallocate(agex)
  allocate(agex(nda,ng))

  if (allocated(px)) deallocate(px)
  if (allocated(py)) deallocate(py)
  if (allocated(pz)) deallocate(pz)
  if (allocated(dex)) deallocate(dex)
  allocate(px(nda),py(nda),pz(nda),dex(nda))
  px=0 ; py=0 ; pz=0 ; dex=0

   
  if (allocated(vcilx)) deallocate(vcilx)
  if (allocated(vcily)) deallocate(vcily)
  if (allocated(vcilz)) deallocate(vcilz)
  allocate(vcilx(nda),vcily(nda),vcilz(nda))
  vcilx=0 ; vcily=0 ; vcilz=0

  if (allocated(vtorx)) deallocate(vtorx)
  if (allocated(vtory)) deallocate(vtory)
  if (allocated(vtorz)) deallocate(vtorz)
  allocate(vtorx(nda),vtory(nda),vtorz(nda))
  vtorx=0 ; vtory=0 ; vtorz=0
  if (allocated(vstorx)) deallocate(vstorx)
  if (allocated(vstory)) deallocate(vstory)
  if (allocated(vstorz)) deallocate(vstorz)
  allocate(vstorx(nda),vstory(nda),vstorz(nda))
  vstorx=0 ; vstory=0 ; vstorz=0
  if (allocated(vsprx)) deallocate(vsprx)
  if (allocated(vspry)) deallocate(vspry)
  if (allocated(vsprz)) deallocate(vsprz)
  allocate(vsprx(nda),vspry(nda),vsprz(nda))
  vsprx=0 ; vspry=0 ; vsprz=0

  if (allocated(erep)) deallocate(erep)
  if (allocated(erepcel)) deallocate(erepcel)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(eyou)) deallocate(eyou)
  if (allocated(espring)) deallocate(espring)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(etor)) deallocate(etor)
  allocate(erep(nda),erepcel(nda),eadh(nda),eyou(nda),espring(nda),etor(nda))

  if (allocated(fmeanl)) deallocate(fmeanl)  !>>Miquel23-1-14
  allocate(fmeanl(nda))                     !>>Miquel23-1-14
  fmeanl=0                                !>>Miquel23-1-14
   
  if (allocated(fmeanv)) deallocate(fmeanv)  !>>Miquel23-1-14
  allocate(fmeanv(nda))                     !>>Miquel23-1-14
  fmeanv=0                                !>>Miquel23-1-14

  ng=10

  call initiate_gene   

!call founding_father
    ffu=0
     !spring of the ellipse
    ffu(2)=1 !to quite if there are too many cells
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=0 !lonely cells and nodes die
    ffu(5)=0
    ffu(6)=1
    ffu(7)=0
    ffu(8)=1
    ffu(10)=0
    ffu(11)=1 ! diffusion of req
    ffu(12)=1 ! cell growth outside the cell

  nd=14*12
  do i=13,ncels
    deallocate(cels(i)%node)
  end do
  ncels=12
  
do i=1,ncels
  !cels(i)%polx=1.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0
  cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0
  cels(i)%maxsize_for_div=cels(i)%nunodes+1
end do

do i=1,ng
  gex(:,i)=0.0d0
  gen(i)%t=0.0
  gen(i)%r=0.0
  gen(i)%e=0.0
  gen(i)%npre=0
  gen(i)%npost=0
  gen(i)%diffu=10
  gen(i)%e(nparam_per_node+8)=0.0d0   
  gen(i)%e(nparam_per_node+1)=1d-4 !0
end do

! FOR SIMPLICITY we make that all genes are pres and posts of others so that it is simpler and that only which ww are not zero matter
do i=1,ng
  if (allocated(gen(i)%pre)) deallocate(gen(i)%pre)
  allocate(gen(i)%pre(ng))
  do j=1,ng
    gen(i)%pre(j)=j
  end do
  if (allocated(gen(i)%post)) deallocate(gen(i)%post)
  allocate(gen(i)%post(ng))
  do j=1,ng
    gen(i)%post(j)=j
  end do
end do

! we also make that all genes are either growth factors or transcriptional factors that can be affected by non-transcriptional catalysis too
do i=1,5
  gen(i)%kindof=6
end do
do i=6,10
  gen(i)%kindof=2
end do

!call update_npag

! now we build the ico

call icosahedra(vertices)

node(:)%marge=1
  do ii=1,12
    i=cels(ii)%node(1)
    node(i)%marge=1
    node(i)%x=vertices(ii,1)*0.6d0
    node(i)%y=vertices(ii,2)*0.6d0
    node(i)%z=vertices(ii,3)*0.6d0
    node(i)%add=node(i)%add*1.4
!    node(i)%eqd=node(i)%eqd*2.8
!    node(i)%cod=node(i)%cod*2.8
!    node(i)%grd=node(i)%grd*2.8
!    node(i)%pld=node(i)%pld*2.8
  end do

  do ii=1,12
    i=node(cels(ii)%node(1))%altre
    node(i)%marge=0
    margs(ii)=i
!    node(i)%tipus=2
    iii=node(i)%altre
    node(i)%x=node(iii)%x*0.25d0
    node(i)%y=node(iii)%y*0.25d0
    node(i)%z=node(iii)%z*0.25d0
    node(i)%add=node(i)%add*1.4
!    node(i)%x=vertices(ii,1)*0.25d0
!    node(i)%y=vertices(ii,2)*0.25d0
!    node(i)%z=vertices(ii,3)*0.25d0
  end do

desmax=1d-5

vv=0
vvv=0
do i=1,ncels
  do j=1,ncels
    ii=cels(i)%node(1)
    iii=cels(j)%node(1)
    a=sqrt((node(ii)%x-node(iii)%x)**2+(node(ii)%y-node(iii)%y)**2+(node(ii)%z-node(iii)%z)**2)
    if (a<=1.51d0) then
      vvv(i)=vvv(i)+1
      vv(i,vvv(i))=j
    end if
  end do
end do

do i=1,ncels
  k=1
  do j=1,vvv(i)
    if (vv(i,j)>0) then
      k=k+1
      ii=cels(i)%node(k)
      node(ii)%x=(2*node(cels(i)%node(1))%x+node(cels(vv(i,j))%node(1))%x)/3.0d0
      node(ii)%y=(2*node(cels(i)%node(1))%y+node(cels(vv(i,j))%node(1))%y)/3.0d0
      node(ii)%z=(2*node(cels(i)%node(1))%z+node(cels(vv(i,j))%node(1))%z)/3.0d0
    end if
  end do
end do

do i=1,ncels
  do j=1,cels(i)%nunodes/2 ! >>> Is 7-6-14
    iii=cels(i)%node(j)
    ii=node(iii)%altre
    a=0.75
    b=0.5
    node(ii)%x=node(iii)%x*b
    node(ii)%y=node(iii)%y*b
    node(ii)%z=node(iii)%z*b
    node(ii)%add=node(ii)%add*a
    node(ii)%eqd=node(ii)%eqd*a
    node(ii)%cod=node(ii)%eqd*a
    node(ii)%grd=node(ii)%grd*a
    node(ii)%pld=node(ii)%pld*a
    !node(ii)%erp=12d0 !node(ii)%erp*1
    !node(ii)%est=node(ii)%est*10
    node(ii)%eqs=0.5d0*sqrt((node(ii)%x-node(iii)%x)**2+(node(ii)%y-node(iii)%y)**2+(node(ii)%z-node(iii)%z)**2)
  end do
  do j=1,cels(i)%nunodes/2
    ii=cels(i)%node(j)
    iii=node(ii)%altre
a=1.2
    node(ii)%add=node(ii)%add*a
    node(ii)%eqd=node(ii)%eqd*a
    node(ii)%cod=node(ii)%eqd*a
    node(ii)%grd=node(ii)%grd*a
    node(ii)%pld=node(ii)%pld*a
    !node(ii)%erp=12d0
    !node(ii)%est=node(ii)%est*10
    node(ii)%eqs=0.5d0*sqrt((node(ii)%x-node(iii)%x)**2+(node(ii)%y-node(iii)%y)**2+(node(ii)%z-node(iii)%z)**2)
  end do
end do

!node(:)%pla=1.d-3

call iniboxes

do i=1,ncels
  call apoptosis(cels(i)%node(1))
end do

do i=1,cels(1)%nunodes
  ii=cels(1)%node(i)
  gex(ii,1)=1.0d0   ! The simplest initial conditions so that there is some asymmetry
end do

gex(:,1)=1.0d0

call update_npag

! now some rotation of some of the cells to make the perfect match

goto 67

do i=1,2
do j=2,6
  ii=cels(i)%node(j)
  e=1d6
  do jj=2,6
    if (j==jj.or.nr(jj)==j) cycle
    iii=cels(i)%node(jj)
    d=sqrt((node(ii)%x-node(iii)%x)**2+(node(ii)%y-node(iii)%y)**2+(node(ii)%z-node(iii)%z)**2)
    if (d<e) then
      e=d ; kk=iii ; kkk=jj
    end if
  end do
  nr(j)=kkk

  npx(i,j)=(node(ii)%x+node(kk)%x)*0.5d0
  npy(i,j)=(node(ii)%y+node(kk)%y)*0.5d0
  npz(i,j)=(node(ii)%z+node(kk)%z)*0.5d0
  snpx(i,j)=(node(node(ii)%altre)%x+node(node(kk)%altre)%x)*0.5d0
  snpy(i,j)=(node(node(ii)%altre)%y+node(node(kk)%altre)%y)*0.5d0
  snpz(i,j)=(node(node(ii)%altre)%z+node(node(kk)%altre)%z)*0.5d0
end do
end do

do i=1,2
do j=2,6
  ii=cels(i)%node(j)
  node(ii)%x=npx(i,j)
  node(ii)%y=npy(i,j)
  node(ii)%z=npz(i,j)
  iii=node(ii)%altre
  node(iii)%x=snpx(i,j)
  node(iii)%y=snpy(i,j)
  node(iii)%z=snpz(i,j)
end do
end do

67 continue

node(:)%marge=1
node(1)%marge=0
node(13)%marge=0
node(26)%marge=0
node(38)%marge=0
node(51)%marge=0
node(63)%marge=0
node(76)%marge=0
node(88)%marge=0
node(101)%marge=0
node(114)%marge=0
node(125)%marge=0
node(138)%marge=0

dif_req=1.d-1    ! Is 11-6-14
node(:)%pla=1d-3
node(:)%kvol=1d-1
prop_noise=1d-2

reqmin=1d-3

do i=1,nd
  node(i)%dmo=desmax
  node(i)%mov=1d2
end do


getot=1d4

end subroutine

!***************************************************************************************************************

subroutine blastula_ensemble

  character*200 carga

  call blastula

  ng=8

!  ffu(11)=0

!  goto 15

     !spring of the ellipse
    ffu=0
    ffu(2)=1 !to quite if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=1 !lonely cells and nodes die
    ffu(5)=0
    ffu(7)=1
    ffu(8)=0
    ffu(15)=1

    !Number of genes
    ng=11
    call initiate_gene
    do i=1,8
      gen(i)%kindof=i
      gen(i)%diffu=1.0d-1
      gen(i)%mu=1.0d-1
    end do


    do i=1,ng 
      do j=1,nga
        gen(i)%e(j)=0.0d0
      end do
      do j=1,ng
        gen(i)%t(j)=0.0d0
      end do
    end do

    gen(3)%kindof=1
    gen(1)%t(1)=1.0d0
    gen(4)%t(1)=1.0d0
    gen(3)%t(4)=-1.0d0
    gen(3)%t(3)=1.0d-1
    gen(2)%t(4)=1.0d0
    gen(4)%diffu=1.0d0

    !Gene expression on nodes
    gex=0.0d0
!    gex(:,1)=1.0d0
!    gex(1:2,2)=1.0d0
    gex(1,1)=1.0d0
    gex(:,3)=1.0d0
    gex(:,5)=1.0d0

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
    gen(2)%e(1)=1.0d0   !adhesion molecule
    gen(1)%e(23)=1.0d-1
    gen(2)%e(23)=1.0d-1
    gen(1)%e(nparam_per_node+1)=1d1 !0
    gen(1)%e(nparam_per_node+1)=1d1 !5d-5
    gen(5)%mu=0.0d0

    ntipusadh=1
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
     
         !****a rellenar****
         kadh(1,1)=1d-3
    end if
    cels(:)%maxsize_for_div=cels(:)%nunodes*2
    
    do i=1,nd
      if (node(i)%tipus==2) then
        if (node(i)%marge==1) then
!          node(i)%marge=0
!          node(node(i)%altre)%marge=1
        end if
      end if
    end do

!    node(:nd)%adh=node(:nd)%adh*5
!    node(:nd)%rec=node(:nd)%rec*2
!    node(:nd)%rep=node(:nd)%rep*0.5

    node(:nd)%erp=node(:nd)%erp*0.5d0
!    node(:nd)%est=node(:nd)%est*10
!    node(:nd)%hoo=node(:nd)%hoo*10

    call update_npag
15  getot=1.0d1

    call getarg(1,carga)
    if (len_trim(carga)==0.or.carga=="0") then
      if (len(ccag)>0.or.ccag=="0") then
        read(ccag,*,ERR=3333) iccag                   !>>> Is 2-3-14
        print *,-iccag,"notice for this ic the random seed idum changes with the date"
        idumoriginal(1)=-iccag               !>>> Is 2-3-14
        call random_seed(put=idumoriginal)   !>>> Is 2-3-14
      end if
    end if
3333 continue

!print *,cels(:ncels)%eqdmax,mmae,"req"

end subroutine blastula_ensemble

!***************************************************************************************************************

subroutine blastula_ensembleold  !makes a gradient from the animal pole and that inhibits an homogeneous gene

  character*200 carga

  call blastula

  ng=11

!  ffu(11)=0

!  goto 15

     !spring of the ellipse
    ffu=0
    ffu(2)=1 !to quite if there are too many cells
    ffu(2)=0 !screening
    ffu(3)=0 !torsion
    
    
    ffu(4)=1 !lonely cells and nodes die
    ffu(5)=0 !ACHTUNG, trial, put that at least to 1
    ffu(7)=1
    ffu(8)=1
    ffu(15)=1

    !Number of genes
    ng=9
    call initiate_gene
    do i=1,ng-1
      gen(i)%kindof=i
      gen(i)%diffu=1.0d-1
      gen(i)%mu=1.0d-1
    end do

    ! catalyzation of the growth factor gene 4
    gen(1)%nww=4
    gen(1)%r(1,1)=2
    gen(1)%r(1,2)=4
    gen(1)%r(1,3)=1.0d-2
!    gen(1)%r(2,1)=3
!    gen(1)%r(2,2)=4
!    gen(1)%r(2,3)=1.0d0
!    gen(1)%r(3,1)=3
!    gen(1)%r(3,2)=5
!    gen(1)%r(3,3)=1.0d0
!    gen(1)%r(4,1)=3
!    gen(1)%r(4,2)=6
!    gen(1)%r(4,3)=1.0d0

    gen(2)%npost=1
    allocate(gen(2)%post(1))
    gen(2)%post(1)=4   !1 t-> 2 -> 4
    gen(2)%e(1)=1.0d0
    gen(2)%t(2)=1.d0
    gen(2)%t(4)=-9.d2 ! direct inhibition by the growth factor 
    gen(2)%diffu=1.0d0
    gen(4)%diffu=1.0d0

!    gen(4)%nww=1
!    gen(4)%r(1,1)=11 !inactive form of the receptor
!    gen(4)%r(1,2)=8  !active (bound) form of the receptor
!    gen(4)%r(1,3)=1.0d0

    gen(9)%kindof=2
    gen(9)%t(1)=1.d0

!    gen(5)%npre=1
!    allocate(gen(5)%pre(1))
!    gen(5)%pre(1)=2
!    gen(6)%npre=1
!    allocate(gen(6)%pre(1))
!    gen(6)%pre(1)=2

!    gen(7)%npre=1
!    allocate(gen(7)%pre(1))
!    gen(7)%pre(1)=2
!    gen(7)%nww=1
!    gen(7)%r(1,1)=10  !this is like homotypic binding
!    gen(7)%r(1,2)=7
!    gen(7)%r(1,3)=1.0d0 !notch needs to catalyze its own synthesis

!    gen(8)%npre=1
!    allocate(gen(8)%pre(1))
!    gen(8)%pre(1)=2
!    gen(8)%npost=1
!    allocate(gen(8)%post(1))
!    gen(8)%post(1)=2        !so that it can revert

    ! the active receptor catalyzes its own inactivation
!    gen(8)%nww=1
!    gen(8)%r(1,1)=8
!    gen(8)%r(1,2)=2
!    gen(8)%r(1,3)=0.1d0

    do i=1,ng 
      do j=1,nga
        gen(i)%e(j)=0.0d0
      end do
    end do

    !Gene expression on nodes
    gex=0.0d0
    do i=1,cels(1)%nunodes
      j=cels(1)%node(i)
      gex(j,1)=1.0d0
    end do
    gex(:,9)=1.0d0
    gex(:,2)=1.0d-2

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
    gen(9)%e(1)=1.0d0   !adhesion molecule
    gen(1)%e(23)=2.d0
    gen(9)%e(23)=2.0d0
    gen(9)%e(nparam_per_node+1)=1d1 !0
    gen(9)%e(nparam_per_node+1)=1d1 !5d-5

    ntipusadh=1
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
     
         !****a rellenar****
         kadh(1,1)=1d-3
    end if
    cels(:)%maxsize_for_div=cels(:)%nunodes*2
    
!    node(:nd)%adh=node(:nd)%adh*5
!    node(:nd)%rec=node(:nd)%rec*2
!    node(:nd)%rep=node(:nd)%rep*0.5

    node(:nd)%erp=node(:nd)%erp*0.5d0
!    node(:nd)%est=node(:nd)%est*10
!    node(:nd)%hoo=node(:nd)%hoo*10

    call update_npag
15  getot=1d4

    call getarg(1,carga)
    if (len_trim(carga)==0.or.carga=="0") then
      if (len(ccag)>0.or.ccag=="0") then
        read(ccag,*,ERR=3333) iccag                   !>>> Is 2-3-14
        print *,-iccag,"notice for this ic the random seed idum changes with the date"
        idumoriginal(1)=-iccag               !>>> Is 2-3-14
        call random_seed(put=idumoriginal)   !>>> Is 2-3-14
      end if
    end if
3333 continue

!print *,cels(:ncels)%eqdmax,mmae,"req"

end subroutine

!***************************************************************************************************************

subroutine icosahedra(vertices)

real*8 vertices(12,3)
real*8 a,b,c,aa,bb,cc
real*8 gold_ratio

gold_ratio=0.5d0*(1+sqrt(5.0d0))

vertices(1,1)=0.0d0 ; vertices(1,2)=1.0d0  ; vertices(1,3)=gold_ratio
vertices(2,1)=0.0d0 ; vertices(2,2)=1.0d0  ; vertices(2,3)=-gold_ratio
vertices(3,1)=0.0d0 ; vertices(3,2)=-1.0d0 ; vertices(3,3)=gold_ratio
vertices(4,1)=0.0d0 ; vertices(4,2)=-1.0d0 ; vertices(4,3)=-gold_ratio

vertices(5,1)=1.0d0  ; vertices(5,2)=gold_ratio  ; vertices(5,3)=0.0d0
vertices(6,1)=1.0d0  ; vertices(6,2)=-gold_ratio ; vertices(6,3)=0.0d0
vertices(7,1)=-1.0d0 ; vertices(7,2)=gold_ratio  ; vertices(7,3)=0.0d0
vertices(8,1)=-1.0d0 ; vertices(8,2)=-gold_ratio ; vertices(8,3)=0.0d0

vertices(9,1)=gold_ratio   ; vertices(9,2)=0.0d0  ; vertices(9,3)=1.0d0
vertices(10,1)=gold_ratio  ; vertices(10,2)=0.0d0 ; vertices(10,3)=-1.0d0
vertices(11,1)=-gold_ratio ; vertices(11,2)=0.0d0 ; vertices(11,3)=1.0d0
vertices(12,1)=-gold_ratio ; vertices(12,2)=0.0d0 ; vertices(12,3)=-1.0d0


end subroutine

!**************************************************************************************************

subroutine tooth_bud



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=1d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-1 !low value is low temperature	
    !desmax=0.001
    resmax=1d-3
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=5d0

    !biological
    !functions used
    
    
    ffu=0
    if(radi==1.or.mradi==1) ffu(1)=0

     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=1 !neighboring by triangulation
    
    ffu(9)=1 !fixed gene 1 on borders
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60; 
        node(i)%hoo=5d1
        node(i)%erp=5d0
        !node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=1d-2
        node(i)%kvol=1d-1
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=3d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        !node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%fix=1 ;node(i)%rec=5d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%fix=1;node(i)%rec=5d1
        node(i)%add=node(i)%eqd*1.50
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      !do i=ndepi+k+1,nd !lower mesenchymal layer all
      !  node(i)%fix=1;node(i)%rec=5d1
      !  node(i)%add=node(i)%eqd*1.50
      !end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=11
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=4 ; gen(1)%diffu=1d-1 ; gen(1)%mu=1d0 !Wnt7b ligand, comes from the epithelial borders
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre)) ; gen(1)%pre(1)=2

      gen(2)%kindof=2 ; gen(2)%diffu=1d0 ; gen(2)%mu=5d-1 !wnt7b transcript
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost)) ; gen(2)%post(1)=2
      
      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=5d-1 !Housekeeping gene (it keeps itself at stable levels, may activate other genes)

      gen(4)%kindof=1 ; gen(4)%diffu=1d0 ; gen(4)%mu=5d-1 !Adhesion molecule epithelial (activated by housekeeping)

      gen(5)%kindof=4 ; gen(5)%diffu=1d-1 ; gen(5)%mu=1d0 !Wnt10b ligand
      gen(5)%npre=1 ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=6
      
      gen(6)%kindof=2 ; gen(6)%diffu=1d0 ; gen(6)%mu=5d-1 !Wnt10b transcript
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost)) ; gen(6)%post(1)=5
      
      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=0d0 !Adhesion molecule external

      gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=5d-1 !Wnt receptor inactive
      gen(8)%npost=2 ; allocate(gen(8)%post(gen(8)%npost)) ; gen(8)%post(1)=9 ; gen(8)%post(2)=10
      
      gen(9)%kindof=3 ; gen(9)%diffu=1d0 ; gen(9)%mu=5d-1 !Wnt receptor activated by Wnt7b
      gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre)) ; gen(9)%pre(1)=8

      gen(10)%kindof=3 ; gen(10)%diffu=1d0 ; gen(10)%mu=5d-1 !Wnt receptor activated by Wnt10b
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre)) ; gen(10)%pre(1)=8
      
      gen(11)%kindof=1 ; gen(11)%diffu=1d0 ; gen(11)%mu=1d-1 !intermediate mediator for Wnt7b (slow but resilient)
     

      
    !Gene-behavior interactions
    

       gen(4)%e(1)=1  !Adhesion molecule epithelial
       gen(7)%e(1)=2  !Adhesion molecule external

       gen(10)%e(nparam_per_node+2)=4d-4 !Wnt10b promotes division



    !Gene-gene interactions

      gen(3)%t(3)=1d0  !housekeeping activates himself
    
      gen(4)%t(3)=1d0  !housekeeping produces adhesion molecule
      gen(8)%t(3)=1d0  !housekeeping produces Wnt receptor inactive

      gen(2)%t(10)=-4.0d0  !receptor activated by Wnt10b inhibits Wnt7b production
      gen(6)%t(10)=4.0d0  !receptor activated by Wnt10b promotes Wnt10b production

      gen(11)%t(9)=1d-1  !wnt7b activated receptor produces mediator at a slow rate

      gen(2)%t(11)=1d0  !mediator produced by Wnt7b promotes Wnt7b production
      gen(6)%t(11)=-1d0  !mediatior produced by Wnt7b inhibits Wnt10b production

      gen(3)%nww=2
      gen(3)%r(1,1)=2  !housekeeping processes Wnt7b secretion
      gen(3)%r(1,2)=1
      gen(3)%r(1,3)=1d0
      gen(3)%r(2,1)=6  !housekeeping processes Wnt10b secretion
      gen(3)%r(2,2)=5
      gen(3)%r(2,3)=1d0

      gen(1)%nww=1
      gen(1)%r(1,1)=8  !Wnt7b activates Wnt receptor into specific activated form
      gen(1)%r(1,2)=9
      gen(1)%r(1,3)=1d1

      gen(5)%nww=1
      gen(5)%r(1,1)=8  !Wnt10b activates Wnt receptor into specific activated form
      gen(5)%r(1,2)=10
      gen(5)%r(1,3)=1d1


    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d1  !external adhesion between hold nodes

    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,nd
      if(node(i)%fix==1)then
        gex(i,1)=1d0
        gex(i,7)=1d0
      elseif(i<=ndepi)then
        gex(i,3)=1d0 !housekeeping
        gex(i,4)=1d0 !adhesion molecule
        gex(i,5)=5d-2 !Wnt10b ligand
        gex(i,6)=5d-2 !Wnt10b transcript
        gex(i,8)=1d0 !receptor
        
      else
        gex(i,3)=1d0
        gex(i,4)=1d0
        gex(i,5)=5d-2 !Wnt10b ligand
        gex(i,6)=5d-2 !Wnt10b transcript
        gex(i,8)=1d0
        !j=0
        !do ii=1,mradicel-2
        !  j=j+ii
        !end do
        !j=(6*j+1) !territory 1 (inner ring)
        !k=0
        !do ii=1,mradicel-1
        !  k=k+ii
        !end do
        !k=(6*k+1) !whole layer
        !if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
        !elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
        !elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
        !elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
        !end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    
end subroutine

subroutine delta_notch



!******* #1 DEFINING SPATIAL DIMENSIONS *******

    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=2    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=3      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
   
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    mesradicel(3)=mradicel

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10    
    
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature   
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=1
    screen_radius=1.0d0
    khold=1d0

    !biological
    !functions used
   
   
    ffu=0 ; if (radi==1.or.mradi==1) ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=1 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=0 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation (epithelial)
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60;
        node(i)%hoo=5d1
        node(i)%erp=5d0
        node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=5d-2
        node(i)%kvol=5d-2
      end do
    end if

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=3d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
   
    !do i=1,nd  !ACTIVATE THIS IF YOU WANT THAT THE INITIAL POSITIONS OF NODES ARE A BIT NOISY
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
!      j=0
!      do i=1,radicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!      do i=j+1,ndepi
!        node(i)%fix=1 ;node(i)%rec=1d2
!      end do
!      j=0
!      do i=1,mradicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      k=0
!      do i=1,mradicel-1
!        k=k+i
!      end do
!      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
!        node(i)%fix=1;node(i)%rec=1d2
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
!      end do
!      do i=ndepi+k+1,nd
!        node(i)%fix=1;node(i)%rec=1d2
!      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
   
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=7
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=7 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0!5d-1 !activated notch
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre)); gen(1)%pre(1)=4
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost)); gen(1)%post(1)=4

     
      gen(2)%kindof=7 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0!5d-1 !activated delta
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre)); gen(2)%pre(1)=5
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost)); gen(2)%post(1)=5


      gen(3)%kindof=1 ; gen(3)%diffu=5d0 ; gen(3)%mu=3d-1 !Housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0!5d-1 !inactive notch
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost)); gen(4)%post(1)=1
      gen(4)%npre=1 ; allocate(gen(4)%pre(gen(4)%npre)); gen(4)%pre(1)=1

     
      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0!5d-1 !inactive delta
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost)); gen(5)%post(1)=2
      gen(5)%npre=1 ; allocate(gen(5)%pre(gen(5)%npre)); gen(5)%pre(1)=2

     
      gen(6)%kindof=2 ; gen(6)%diffu=5d0 ; gen(6)%mu=5d-1 !effector gene inactive
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost)); gen(6)%post(1)=7
     
      gen(7)%kindof=3 ; gen(7)%diffu=5d0 ; gen(7)%mu=5d-1 !effector gene active
      gen(7)%npre=1 ; allocate(gen(7)%pre(gen(7)%npre)); gen(7)%pre(1)=6
     
     
    !Gene-behavior interactions
   

       !gen(1)%e(1)=1  !epithelial adhesion molecule
       !gen(2)%e(1)=2  !mesenchymal adhesion molecule
       !gen(3)%e(1)=3  !basal lamina


    !Gene-gene interactions

      !gen(3)%t(3)=1d0 !housekeeping epi. maintains its levels of expressions
      !gen(4)%t(3)=1d0
      !gen(5)%t(3)=1d0
      !gen(6)%t(3)=1d0

      gen(4)%nww=2    !inactive notch mediates binding of delta
      gen(4)%r(1,1)=5
      gen(4)%r(1,2)=2
      gen(4)%r(1,3)=1d0
      gen(4)%r(2,1)=2
      gen(4)%r(2,2)=5
      gen(4)%r(2,3)=1d0

      gen(5)%nww=2    !inactive delta mediates binding of notch
      gen(5)%r(1,1)=4
      gen(5)%r(1,2)=1
      gen(5)%r(1,3)=1d0
      gen(5)%r(2,1)=1
      gen(5)%r(2,2)=4
      gen(5)%r(2,3)=1d0

     
      !gen(1)%nww=2
      !gen(1)%r(1,1)=6     !activated notch activates effector
      !gen(1)%r(1,2)=7
      !gen(1)%r(1,3)=1d-1
      !gen(1)%r(2,1)=2     !activated notch mediates unbinding of delta
      !gen(1)%r(2,2)=5
      !gen(1)%r(2,3)=1d0

      !gen(2)%nww=1     !active delta mediates unbinding of notch
      !gen(2)%r(1,1)=1
      !gen(2)%r(1,2)=4
      !gen(2)%r(1,3)=1d0    
     
     
      !gen(5)%t(7)=-1d1


    !Adhesion molecules

    ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !kadh(1,1)=2d1
      !kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      !kadh(2,2)=2d1
      !kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d2
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=3,ndepi
      !if(node(i)%tipus==1)  gex(i,1)=1d0
      !if(node(i)%tipus==2)  gex(i,3)=1d0

      gex(i,4)=1d0
      !gex(i,5)=1d0
      !gex(i,3)=1d0
      !gex(i,6)=1d0
    end do

    gex(1:2,5)=1d0 ; !gex(1,5)=1d0
    !gex(3,5)=1d0 ; !gex(3,4)=1d0

    
    !do i=ndepi+1,nd
    !  gex(i,2)=1d0
    !  gex(i,5)=1d0
    !end do



    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3
mmae=node(1)%eqd
   
end subroutine

subroutine receptor_ligand_test



!******* #1 DEFINING SPATIAL DIMENSIONS *******

    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=4    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=3      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
   
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    mesradicel(3)=mradicel

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10    
    
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature   
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=1
    screen_radius=1.0d0
    khold=1d0

    !biological
    !functions used
   
   
    ffu=0 ; if (radi==1.or.mradi==1) ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=1 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation (epithelial)
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60;
        node(i)%hoo=5d1
        node(i)%erp=5d0
        node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=5d-2
        node(i)%kvol=5d-2
      end do
    end if

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=3d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
   
    !do i=1,nd  !ACTIVATE THIS IF YOU WANT THAT THE INITIAL POSITIONS OF NODES ARE A BIT NOISY
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
!      j=0
!      do i=1,radicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!      do i=j+1,ndepi
!        node(i)%fix=1 ;node(i)%rec=1d2
!      end do
!      j=0
!      do i=1,mradicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      k=0
!      do i=1,mradicel-1
!        k=k+i
!      end do
!      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
!        node(i)%fix=1;node(i)%rec=1d2
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
!      end do
!      do i=ndepi+k+1,nd
!        node(i)%fix=1;node(i)%rec=1d2
!      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
   
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=8 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 !bound receptor-ligand
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre)); gen(1)%pre(1)=4 !; gen(1)%pre(2)=2
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost)); gen(1)%post(1)=4 !; gen(1)%post(2)=2

      !gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost)); gen(1)%post(1)=4

     
      gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=5d-1 !ligand
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre)); gen(2)%pre(1)=5
      !gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost)); gen(2)%post(1)=1


      gen(3)%kindof=1 ; gen(3)%diffu=5d0 ; gen(3)%mu=3d-1 !Housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=5d0 ; gen(4)%mu=5d-1 !inactive receptor
      !gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost)); gen(4)%post(1)=1
      !gen(4)%npre=1 ; allocate(gen(4)%pre(gen(4)%npre)); gen(4)%pre(1)=1

     
      gen(5)%kindof=2 ; gen(5)%diffu=5d0 ; gen(5)%mu=5d-1 !ligand transcript
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost)); gen(5)%post(1)=2
      !gen(5)%npre=1 ; allocate(gen(5)%pre(gen(5)%npre)); gen(5)%pre(1)=2

     
      !gen(6)%kindof=2 ; gen(6)%diffu=5d0 ; gen(6)%mu=5d-1 !effector gene inactive
      !gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost)); gen(6)%post(1)=7
     
      !gen(7)%kindof=3 ; gen(7)%diffu=5d0 ; gen(7)%mu=5d-1 !effector gene active
      !gen(7)%npre=1 ; allocate(gen(7)%pre(gen(7)%npre)); gen(7)%pre(1)=6
     
     
    !Gene-behavior interactions
   

       !gen(1)%e(1)=1  !epithelial adhesion molecule
       !gen(2)%e(1)=2  !mesenchymal adhesion molecule
       !gen(3)%e(1)=3  !basal lamina


    !Gene-gene interactions

      gen(3)%t(3)=1d0 !housekeeping epi. maintains its levels of expressions
      gen(4)%t(3)=1d0
      gen(5)%t(3)=1d0

      gen(3)%nww=1    !housekeeping epi. processes ligand secretion
      gen(3)%r(1,1)=5
      gen(3)%r(1,2)=2
      gen(3)%r(1,3)=1d1


      gen(1)%nww=2
      gen(1)%r(1,1)=2  !assotiation
      gen(1)%r(1,2)=4
      gen(1)%r(1,3)=1d1
      gen(1)%r(2,1)=4  !dissociation
      gen(1)%r(2,2)=2
      gen(1)%r(2,3)=1d1



    !Adhesion molecules

    ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !kadh(1,1)=2d1
      !kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      !kadh(2,2)=2d1
      !kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d2
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,ndepi
      !if(node(i)%tipus==1)  gex(i,1)=1d0
      !if(node(i)%tipus==2)  gex(i,3)=1d0

      gex(i,4)=1d0
      gex(i,5)=1d0
      gex(i,3)=1d0
      gex(i,6)=1d0

    end do

    !do i=ndepi+1,nd
    !  gex(i,2)=1d0
    !  gex(i,5)=1d0
    !end do



    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3
mmae=node(1)%eqd
   
end subroutine

!*****************************************************************
!**********************************************************************************************************************************************
subroutine tooth40a
real*8::pericel,radio2,req1,req2,adjustreq,di
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=6    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=1.0d0
    ecmmax=req
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=0 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(10)=0 !conservacion del volumen
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if
    print*,"node(ndepi+1)%eqd",node(ndepi+1)%eqd
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
          node(i)%eqd=req*adjustreq ; node(i)%grd=node(i)%eqd
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
      node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%eqd1",node(i)%eqd
!    else       
!   
!    print*,"node(i)%eqd2",node(i)%eqd       
!    end if
   
    node(i)%eqs=req
        node(i)%add=node(i)%eqd*1.80;
        node(i)%hoo=1d1
        node(i)%erp=43d0
        node(i)%est=53d0
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=1d1
    node(i)%kvol=1d0
    node(i)%pla=1.9d-1
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)
    maxdidare=node(1)%eqd*0.1d0


    !Distribution of nodes in space
!    call epiteli_sphere(radi,radicel,radius)
!    do i=1,nd
!      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
!    end do
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    do i=j+1,ndepi
        node(i)%fix=1 ;node(i)%rec=1d1!
      end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%fix=1 ;node(i)%rec=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !node(i)%est=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=8
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=4d0 ; gen(1)%kindof=1 ; gen(1)%mu=1d-1     !adhesion epi., house keeping
    gen(2)%diffu=2d0 ; gen(2)%kindof=1 ; gen(2)%mu=1d-1    !adhesion mesen., house keeping
    gen(8)%diffu=1d0 ; gen(8)%kindof=1 ; gen(8)%mu=1d0
    gen(6)%diffu=1d0 ; gen(6)%kindof=1 ; gen(6)%mu=1d0    !celular div. mesen.
   
    gen(4)%diffu=1d0 ; gen(4)%kindof=2 ; gen(4)%mu=1d0    !pre extracelular signal
    gen(4)%npost=1  ; allocate(gen(4)%post(gen(4)%npost)) ; gen(4)%post(1)=5
   
    gen(5)%diffu=5d1 ; gen(5)%kindof=4 ; gen(5)%mu=1d1
    gen(5)%npre=1  ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=4    !extracel signal, promotes division from epi

    gen(3)%diffu=1d1 ; gen(3)%kindof=2 ; gen(3)%mu=1d0    !pre division epi
    gen(3)%npost=1  ; allocate(gen(3)%post(gen(3)%npost)) ; gen(3)%post(1)=7

    gen(7)%diffu=1d1 ; gen(7)%kindof=3 ; gen(7)%mu=5d1    !promotes epi division (activated by 5)
    gen(7)%npre=1  ; allocate(gen(7)%pre(gen(7)%npre)) ; gen(7)%pre(1)=3   

   
    !Gene-behavior interactions
    gen(1)%e(1)=1   
    gen(2)%e(1)=2
    gen(7)%e(nparam_per_node+2)=2.8999d-1
       
    gen(6)%e(nparam_per_node+2)=5.9d-4!40********************************************2.4d-4
    gen(8)%e(nparam_per_node+2)=4d-5!40*********************************************5d-5

   

    !Gene-gene interactions
        gen(1)%t(1)=5.0d1
    gen(3)%t(1)=1.0d0
    gen(8)%t(1)=1.0d2
        gen(2)%t(2)=5.0d1
    gen(6)%t(2)=1.0d2
    gen(4)%t(2)=1d1
           
    gen(2)%nww=1
    gen(2)%r(1,1)=4
    gen(2)%r(1,2)=5
    gen(2)%r(1,3)=1d2 

    gen(5)%nww=1
    gen(5)%r(1,1)=3
    gen(5)%r(1,2)=7
    gen(5)%r(1,3)=1d2 

    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=1d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=10d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d0
   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,ndepi
    gex(i,1)=5d0
        !gex(i,3)=5d0
    !if (node(i)%tipus==2) gex(i,7)=2d0
    end do

    do i=ndepi+1,nd
        gex(i,2)=5d0
   
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
   
 

end subroutine tooth40a
!*************************************************************************************************************************************
!*****************************************************************
!*****************************************************************
!**********************************************************************************************************************************************
subroutine single_node_sphere
real*8::pericel,radio2,req1,req2,adjustreq,di,ax,ay,az,a,b,c,aa,bb,cc
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=15    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      k=0
      do i=1,radicel-2
        k=k+i
      end do

      ncelsepi=(6*j+1)+(6*k+1) ; print*,"ncelsepi",ncelsepi,"j",j

      !ncelsepi=(6*j+1)*2
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if




    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if

    ndx=ncelsepi



    nd=ndepi+ndmes+ndx
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    !angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=1d1
    ecmmax=0.10d0
    screen_radius=0.7d0
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=0 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(10)=0 !conservacion del volumen
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if
    print*,"node(ndepi+1)%eqd",node(ndepi+1)%eqd
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=5d0
          node(i)%eqd=req*adjustreq ; node(i)%grd=node(i)%eqd
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=5d0
      node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%eqd1",node(i)%eqd
!    else       
!   
!    print*,"node(i)%eqd2",node(i)%eqd       
!    end if
   
    node(i)%eqs=req
        node(i)%add=node(i)%eqd*1.80;
        node(i)%hoo=1d1
        node(i)%erp=1d1
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=1d1
    node(i)%kvol=5d0
    node(i)%pla=5.0d0
        node(i)%ecm=0d0
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%ecm=0d0

      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)
    maxdidare=node(1)%eqd*0.1d0


    !Distribution of nodes in space


      !putting some ECM nodes representing the hyaline layer
      do i=1,ncels
        ii=ndepi+i
        k=cels(i)%node(2)
        a=node(k)%x ; b=node(k)%y ; c=node(k)%z
        kk=cels(i)%node(1)
        aa=node(kk)%x ; bb=node(kk)%y ; cc=node(kk)%z
        ax=a-aa ; ay=b-bb ; az=c-cc
        d=sqrt(ax**2+ay**2+az**2)
        ax=2*node(k)%eqd*ax/d ; ay=2*node(k)%eqd*ay/d ; az=2*node(k)%eqd*az/d

        node(ii)%x=a+ax ; node(ii)%y=b+ay ; node(ii)%z=c+az
        !nodeo(ii)%x=a+ax ; nodeo(ii)%y=b+ay ; nodeo(ii)%z=c+az
        node(ii)%orix=a+ax ; node(ii)%oriy=b+ay ; node(ii)%oriz=c+az

        node(ii)%tipus=4 ; node(ii)%icel=-ii 
        node(ii)%eqd=0.15d0 ; node(ii)%grd=node(ii)%eqd
        node(ii)%rec=5d1 ; node(ii)%add=node(ii)%eqd*1.7d0
        node(ii)%fix=2
      end do


   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    !node(:)%fix=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    !do i=j+1,ndepi
        !node(i)%fix=1 ;node(i)%rec=1d1!
     ! end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%fix=1 ;node(i)%rec=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !node(i)%est=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=0d0 ; gen(1)%kindof=1 ; gen(1)%mu=0d0     !adhesion epi., house keeping
    gen(2)%diffu=0d0 ; gen(2)%kindof=1 ; gen(2)%mu=0d0    !adhesion mesen., house keeping

    gen(3)%diffu=0d0 ; gen(3)%kindof=1 ; gen(3)%mu=0d0    !pre division epi
   
    gen(4)%diffu=0d0 ; gen(4)%kindof=1 ; gen(4)%mu=0d0    !pre extracelular signal
   
    gen(5)%diffu=0d0 ; gen(5)%kindof=1 ; gen(5)%mu=0d0

    !Gene-behavior interactions
    gen(1)%e(1)=1   
    gen(2)%e(1)=2
    gen(5)%e(1)=3
    gen(4)%e(1)=4

    !gen(3)%e(21)=-0.05d0
    !gen(3)%e(27)=-1d2
    !gen(3)%e(28)=-1d2
    !gen(1)%e(27)=-1d2

    gen(3)%e(nparam_per_node+4)=5d-1
    gen(4)%e(nparam_per_node+5)=1d0
    gen(3)%e(nparam_per_node+6)=5d0
    gen(3)%e(nparam_per_node+7)=0.20d0



 
    !Gene-gene interactions

    gen(4)%t(3)=1d0
 
    !Adhesion molecules

    ntipusadh=4
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=5d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=5d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0

      kadh(1,3)=5d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=5d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=5d1

      kadh(1,4)=5d0 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d0 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=1d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=5d0 

   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,j
      k=cels(i)%node(2)
      gex(k,1)=1d0
      gex(node(k)%altre,1)=1d0
      gex(k,3)=1d0
    end do
    do i=j+1,ncels
      k=cels(i)%node(2)
      gex(k,2)=1d0
      gex(node(k)%altre,2)=1d0
    end do
    do i=ndepi+1,nd
      gex(i,5)=1d0
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
   
 

end subroutine single_node_sphere

subroutine tooth40_ligand
real*8::pericel,radio2,req1,req2,adjustreq,di
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=6    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    !angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=0d0
    ecmmax=req
    mnn=1000
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=0 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    ffu(10)=0 !conservacion del volumen
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if
    print*,"node(ndepi+1)%eqd",node(ndepi+1)%eqd
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
          node(i)%eqd=req*adjustreq ; node(i)%grd=node(i)%eqd
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
      node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%eqd1",node(i)%eqd
!    else       
!   
!    print*,"node(i)%eqd2",node(i)%eqd       
!    end if
   
    node(i)%eqs=req
        node(i)%add=node(i)%eqd*1.80;
        node(i)%hoo=1d1
        node(i)%erp=43d0
        node(i)%est=43d0
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=1d1
    node(i)%kvol=1d0
    node(i)%pla=1.9d-1
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%rec=1d1
        node(i)%eqd=req ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=desmax
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)
    maxdidare=node(1)%eqd*0.1d0


    !Distribution of nodes in space
!    call epiteli_sphere(radi,radicel,radius)
!    do i=1,nd
!      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
!    end do
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    do i=j+1,ndepi
        node(i)%fix=1 ;node(i)%rec=1d1!
      end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%fix=1 ;node(i)%rec=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !node(i)%est=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%eqdmax=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=10
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=4d0 ; gen(1)%kindof=1 ; gen(1)%mu=1d-1     !adhesion epi., house keeping
    gen(2)%diffu=2d0 ; gen(2)%kindof=1 ; gen(2)%mu=1d-1    !adhesion mesen., house keeping
    gen(8)%diffu=1d0 ; gen(8)%kindof=1 ; gen(8)%mu=1d0
    gen(6)%diffu=1d0 ; gen(6)%kindof=1 ; gen(6)%mu=1d0    !celular div. mesen.
   
    gen(4)%diffu=1d0 ; gen(4)%kindof=2 ; gen(4)%mu=1d0    !pre extracelular signal
    gen(4)%npost=1  ; allocate(gen(4)%post(gen(4)%npost)) ; gen(4)%post(1)=5
   
    gen(5)%diffu=5d1 ; gen(5)%kindof=4 ; gen(5)%mu=1d1
    gen(5)%npre=1  ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=4    !extracel signal, promotes division from epi

    gen(3)%diffu=1d1 ; gen(3)%kindof=2 ; gen(3)%mu=1d0    !pre division epi
    !gen(3)%npost=1  ; allocate(gen(3)%post(gen(3)%npost)) ; gen(3)%post(1)=7
    !gen(3)%npre=1  ; allocate(gen(3)%pre(gen(3)%npre)) ; gen(3)%pre(1)=7


    gen(7)%diffu=0d0 ; gen(7)%kindof=8 ; gen(7)%mu=0d0    !receptor activated by 5 promotes epi division (activated by 5)
    gen(7)%npre=2  ; allocate(gen(7)%pre(gen(7)%npre)) ; gen(7)%pre(1)=3 ; gen(7)%pre(2)=5
    gen(7)%npost=2  ; allocate(gen(7)%post(gen(7)%npost)) ; gen(7)%post(1)=3 ; gen(7)%post(2)=5

   

    gen(9)%diffu=5d-2 ; gen(9)%kindof=4 ; gen(9)%mu=5d-1  !signal from bead

    gen(10)%diffu=0d0 ; gen(10)%kindof=8 ; gen(10)%mu=0d0    !receptor activated by 5 promotes epi division (activated by 5)
    gen(10)%npre=2  ; allocate(gen(10)%pre(gen(10)%npre)) ; gen(10)%pre(1)=3 ; gen(10)%pre(2)=9
    gen(10)%npost=2  ; allocate(gen(10)%post(gen(10)%npost)) ; gen(10)%post(1)=3 ; gen(10)%post(2)=9

    !Gene-behavior interactions
    gen(1)%e(1)=1   
    gen(2)%e(1)=2
    gen(7)%e(nparam_per_node+2)=2.8999d-1
       
    gen(6)%e(nparam_per_node+2)=2.4d-4!40********************************************2.4d-4
    gen(8)%e(nparam_per_node+2)=5d-5!40*********************************************5d-5

   

    !Gene-gene interactions
        gen(1)%t(1)=5.0d1
    gen(3)%t(1)=5.0d1
    gen(8)%t(1)=1.0d2
        gen(2)%t(2)=5.0d1
    gen(6)%t(2)=1.0d2
    gen(4)%t(2)=1d1
           
    gen(2)%nww=1
    gen(2)%r(1,1)=4
    gen(2)%r(1,2)=5
    gen(2)%r(1,3)=1d2 

    !gen(5)%nww=1
    !gen(5)%r(1,1)=3
    !gen(5)%r(1,2)=7
    !gen(5)%r(1,3)=1d2 

    !gen(9)%nww=1
    !gen(9)%r(1,1)=3
    !gen(9)%r(1,2)=7
    !gen(9)%r(1,3)=1d2 


    gen(7)%nww=4
    gen(7)%r(1,1)=3
    gen(7)%r(1,2)=7
    gen(7)%r(1,3)=1d2 
    gen(7)%r(2,1)=7
    gen(7)%r(2,2)=3
    gen(7)%r(2,3)=1d2
    gen(7)%r(3,1)=5
    gen(7)%r(3,2)=7
    gen(7)%r(3,3)=1d2 
    gen(7)%r(4,1)=7
    gen(7)%r(4,2)=5
    gen(7)%r(4,3)=1d2 


    gen(10)%nww=4
    gen(10)%r(1,1)=3
    gen(10)%r(1,2)=10
    gen(10)%r(1,3)=1d2 
    gen(10)%r(2,1)=10
    gen(10)%r(2,2)=3
    gen(10)%r(2,3)=1d2
    gen(10)%r(3,1)=9
    gen(10)%r(3,2)=10
    gen(10)%r(3,3)=1d2 
    gen(10)%r(4,1)=10
    gen(10)%r(4,2)=9
    gen(10)%r(4,3)=1d2 


    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=3d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=10d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d0
   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,ndepi
    gex(i,1)=5d0
    gex(i,3)=5d-1

        !gex(i,3)=5d0
    !if (node(i)%tipus==2) gex(i,7)=2d0
    end do

    do i=ndepi+1,nd
        gex(i,2)=5d0
   
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
   
 

end subroutine tooth40_ligand




!*************************************************************************************************************************************
!*****************************************************************


subroutine epiteli_sphere1(radi,radicel,radius,di)        !radi is the number of cell rows from one pole to the equator
integer            ::radi,radicel,valcel

integer::trobat,cont,rep,n,i,j,k,ii,jj,kk,lev,numdos,ncels2,val,val2,ll,val3,val4,val5,val6,iiii,jjjj,kkkk,suma,cocels
real*8::modul,pescab,pescac,pescbc,modula,modulb,modulc,aax,aay,aaz,bbx,bby,bbz,ccx,ccy,ccz,costat,rhex,lad,ax,ay,az,di
real*8::alf,bet,gam,a,b,c,l,aa,bb,cc,aaa,bbb,ccc,angle,rpent,modmin,pesc,radius,radius2,bx,by,bz,cost,ucost,sint,thet,radicelepisphere
real*8,dimension(:)::minu(3),vec1(3),vec2(3),vec3(3),vec4(3),vec5(3)
real*8,allocatable::cmalla(:,:),cveci(:,:),primers(:)




    di=di*0.80
    de=0.5d0

    radicelepisphere=(radicel-1)*2

    !SPHERE CODE, FIRST HEMISPHERE
    cont=0 ; cocels=0
    ii=radicelepisphere/2
!    radius=5d0  !radius of the sphere
    radius2=radius+di

    alf=pi/radicelepisphere

!!!pole cell
    cont=cont+1 ; cocels=cocels+1
    node(1)%x=0d0 ; node(1)%y=0d0 ; node(1)%z=radius
    node(1)%altre=2 ; node(1)%tipus=2 ; node(1)%icel=cocels
    cont=cont+1
    node(2)%x=0d0 ; node(2)%y=0d0 ; node(2)%z=radius2
    node(2)%altre=1 ; node(2)%tipus=1 ; node(2)%icel=cocels

    !el radi 2 de la cèlula!
    if (radi>1) then !111
    cont=cont+1    !pone un nodo adyacente al inicial
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius
    node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius2
    node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

    gam=2*pi/6d0
    do i=1,5  !pone resto de nodos
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius
      node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius2
      node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
    end do
    end if !111
!!!!! end pole cell

!!!!! spine cell
    do i=1,ii
       
!  print*,"*********espina******",i
        cont=cont+1 ; cocels=cocels+1    !hace columna de celulas, 1. nodo
        suma=0
        jj=cont
        angle=alf*i
!        print*,"cont",cont
        vec1(1)=radius*sin(angle)
        vec1(2)=0
        vec1(3)=radius*cos(angle)
        node(cont)%x=vec1(1) ; node(cont)%y=vec1(2) ; node(cont)%z=vec1(3)
        node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!        print*,"altre",node(cont)%altre

        cont=cont+1
        vec2(1)=radius2*sin(angle)
        vec2(2)=0
        vec2(3)=radius2*cos(angle)
        node(cont)%x=vec2(1) ; node(cont)%y=vec2(2) ; node(cont)%z=vec2(3)
        node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

        ux=vec2(1)-vec1(1) ; uy=vec2(2)-vec1(2) ; uz=vec2(3)-vec1(3)
        a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a

        ax=-de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"vora"," a",ax,ay,az

    if (radi>1) then !111
        do j=1,6    !pone nodos restantes alrededor 1. nodo
          cont=cont+1
          thet=j*gam
   print*,"thet",thet
          cost=cos(thet); ucost=1-cost ; sint=sin(thet)
  print*,"cos",cost,"ucost",ucost,"sint",sint
          bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
          by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
          bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
          node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
          node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!   print*,"b",bx,by,bz
!  print*,"vec1",vec1
!  print*,"nodecont",node(cont)%x,node(cont)%y,node(cont)%z
          cont=cont+1

          node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
          node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
        end do
    end if !111
        !!!!!!!!! end spine cell



!        print*,"altre",node(cont)%altre
      if (radi>1) then  !111
        if(i+1<=ii+1)then
            kk=6+6*(i-1)
        else
            kk=6+6*(radicelepisphere-i-1)
    end if !111
      else    !111
    kk=6+6*(i-1)    !111
      end if
   
        bet=2*pi/kk

!        print*,"i",i,"nombre de cels del paralel",kk

        do j=1,kk-1

            !!!!!!!!!!! rib cell

            angle=alf

            cont=cont+1 ; cocels=cocels+1
            a=node(1)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
            b=node(1)%y
            c=node(1)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
            b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
            c=node(cont)%z
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
            vec1(1)=node(cont)%x ; vec1(2)=node(cont)%y ; vec1(3)=node(cont)%z

            cont=cont+1
            a=node(2)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
            b=node(2)%y
            c=node(2)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
            b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
            c=node(cont)%z
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

            ux=node(cont)%x-node(cont-1)%x
            uy=node(cont)%y-node(cont-1)%y
            uz=node(cont)%z-node(cont-1)%z
            a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a !;print*,"spring",a
            vec2(1)=node(cont)%x ; vec2(2)=node(cont)%y ; vec2(3)=node(cont)%z

            ax=-de*sin(j*bet) ; ay=de*cos(j*bet) ; az=0d0

!        ax=de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"bet",j*bet,"a",ax,ay,"modul",sqrt(ax**2+ay**2)
    if (radi>1) then           
    do k=1,6
              cont=cont+1
              thet=k*gam
              cost=cos(thet); ucost=1-cost ; sint=sin(thet)
              bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
              by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
              bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
              node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
              node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels

!              print*,"b rib",bx,by,bz,"modul",sqrt(bx**2+by**2+bz**2)

              cont=cont+1

              node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
              node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
           end do
    end if


            !!!!!!!!!!! end rib cell



        end do
    end do

!print*,"cont",cont,"cocels",cocels


!the 2ond hemisphere

    j=0
!    do i=1,(radicel/2+1)-2
    do i=1,radicel-2
      j=j+i
    end do
    ii=(6*j+1)*2   !number of nodes on the 2ond hemisphere, si se quita sale la mitad de la esfera

    print*,"cont",cont,"ii",ii,"(radicel/3+1)-2",(radicel/2+1)-2
    do i=1,ii
      cont=cont+1
      node(cont)%x=node(i)%x
      node(cont)%y=node(i)%y
      node(cont)%z=-node(i)%z
      node(cont)%tipus=node(i)%tipus
      if(node(i)%tipus==2)then
        node(cont)%altre=cont+1
      else
        node(cont)%altre=cont-1
      end if
      node(cont)%icel=node(i)%icel+cocels
    end do

!print*,"cont complet",cont,"cocels",cocels


!    cont=cont+1
!    malla(cont,1)=0;malla(cont,2)=0;malla(cont,3)=-radius
!    primers(radi+1)=cont
!    print*,i,"malla",malla(cont,:)
!    print*,"nombre de cels creades",cont
!    print*,"nombre de paral·lels",lev





    !define the cell's centroid

    do i=1,ncelsepi
      kk=0
      cels(i)%ctipus=1    !tipo celular 1, epitelio
      cels(i)%nunodes=nodecel !
      cels(i)%nodela=nodecela
      cels(i)%cex=0d0 ; cels(i)%cey=0d0 ; cels(i)%cez=0d0
      cels(i)%polx=0d0 ; cels(i)%poly=0d0 ; cels(i)%polz=0d0
      allocate(cels(i)%node(nodecela))
      do j=1,ndepi
        k=node(j)%icel
        if(k==i)then
          kk=kk+1
          cels(i)%node(kk)=j
          if(node(j)%tipus==1)then
            cels(i)%cex=cels(i)%cex+node(j)%x
            cels(i)%cey=cels(i)%cey+node(j)%y
            cels(i)%cez=cels(i)%cez+node(j)%z
          end if
        end if
      end do
      cels(i)%cex=2*cels(i)%cex/cels(i)%nunodes
      cels(i)%cey=2*cels(i)%cey/cels(i)%nunodes
      cels(i)%cez=2*cels(i)%cez/cels(i)%nunodes
    end do

    do i=1,ndepi
      if(node(i)%tipus==1) node(i)%marge=0
      if(node(i)%tipus==2) node(i)%marge=1
    end do


!    do i=1,ncelsepi
!        cels(i)%ctipus=1
!        cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
!        do j=1,cels(i)%nunodes
!            k=cels(i)%node(j)
!            if(node(k)%tipus==1)then
!                 cels(i)%cex=cels(i)%cex+node(k)%x
!                 cels(i)%cey=cels(i)%cey+node(k)%y
!                cels(i)%cez=cels(i)%cez+node(k)%z
!            end if
!        end do
!        cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
!        cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
!        cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
!    end do


 


    realtime=0
    maxdidare=node(1)%eqd*0.1d0
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%add)*3
    node(:)%dif=0.0d0
end subroutine epiteli_sphere1
!*************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesenq1(radi,radicel,layer,zmes,dreq)                                    !!! miguel 20.11
integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,ent,ont       ! number of layers and concentric hexagons
real*8   ::rad,der,zmes,de,di,dreq
real*8   :: xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
dreq=0.05
de=dreq*2                    ! call radius
!di=2.0d0*de                 ! distance between cells and layers (it has to be >2 to avoid cell contacts)
di=2*de+2*de*cos(60d0*rad)
print*,"dreq",dreq
    do i=ncelsepi+1,ncels              
!        cels(i)%nunodes=radi+1
        cels(i)%nunodes=radi                !>>>>>>>>>>>>>>Miquel 21-3-13
        cels(i)%nodela=nodecela
!        allocate(cels(i)%node(radi+1))
        allocate(cels(i)%node(nodecela))        !>>>>>>>>>>>>>>Miquel 21-3-13
        cels(i)%node=0
    end do

node(ndepi+1:)%marge=1

!    print*,"node",size(node),"cels",size(cels),"celsnode",size(cels(1)%node)


    kkk=ndmes/layer
! radi=radi+1
    cont=ncelsepi+1

    ii=ndepi        !node counter
    jj=ncelsepi        !cel counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! central nodes (=nucleus) !!!!

   do l=1,layer 
!    jj=(kkk*(l-1))+radi+ndepi+1
    ii=ii+1
    jj=jj+1
!    node(ii)%x=(sqrt(di/2d0))*mod(l+1,2) ; node(ii)%y=(sqrt(di/2d0))*mod(l+1,2) ! Euler packing algorhytm
    node(ii)%x=0d0 ; node(ii)%y=0d0 ! Euler packing algorhytm
    node(ii)%z=zmes-di*(l-1)    !zmes marks the z-position of the uppermost layer
!    print*,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
    node(ii)%tipus=3 ; node(ii)%icel=jj                        ! origin   
    cels(jj)%node(1)=ii
!    print*,"cel central, node central. ii:",ii
    ent=ii
    ont=ii
    !fill with the external nodes
    if(radi>1)then
        signo=1                                                              ! miguel 4-6-13
        do k=2,radi
            ii=ii+1
            !print*,"cel perif. ent:",ent,ii,k-1,
             call random_number(a)
            der=(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
121         call random_number(a)
            xx=der*(1d0-2*a) ;
            call random_number(a)
            yy=der*(1d0-2*a)            ! miguel 4-6-13
            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
            if(zz.lt.0)then;goto 121;endif                                   ! miguel 4-6-13
            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz                ! miguel 4-6-13
            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz                ! miguel 4-6-13
              else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if       ! miguel 4-6-13      
            !a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de   ! miguel 4-6-13    
            node(ii)%x=node(ent)%x+a
            node(ii)%y=node(ent)%y+b
            node(ii)%z=node(ent)%z+c
            node(ii)%tipus=3 ; node(ii)%icel=jj
            cels(jj)%node(k)=ii
        end do
    end if

    do j=2,radicel                                                        ! "cicle"
        dx1=0.0 ; dx2=0.0 ; dx3=0.0 ; dy1=0.0 ; dy2=0.0 ; dy3=0.0
        do i=1,6                                                        ! "sector" of the hexagon
            if(i.eq.1)then  

                dx2=0!*cos(real(i)*60d0*rad)
                dy2=di*(real(j)-1d0)!*sin(real(i)*60d0*rad);print*,"dx2",dx2,"dy2",dy2
                dx1=di*(real(j)-1d0)*sin(-60d0*rad) ; dy1=di*(real(j)-1d0)*cos(-60d0*rad)!;print*,"dx1",dx1

            else
                hip=di*(real(j)-1d0)                                    ! hipotenusa
                dx1=dx2 ; dy1=dy2
                dx2=hip*sin(real(i-1)*60d0*rad) ; dy2=hip*cos(real(i-1)*60d0*rad)         
            end if           
            ii=ii+1
            jj=jj+1
!            print*,"cel perif, node central. ent:",ent,ii
            node(ii)%x=node(ont)%x+dx2 ; node(ii)%y=node(ont)%y+dy2 ; node(ii)%z=node(ont)%z
!            print*,"vertex i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
!            print*,i,"ii",node(ii)%x,"ent",node(ent)%x
            node(ii)%icel=jj ; node(ii)%tipus=3
            cels(jj)%node(1)=ii
            ent=ii
            if(radi>1)then
                signo=1                                                                  ! miguel 4-6-13
                do k=2,radi
                    ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
122                         call random_number(a)
                            xx=der*(1d0-2*a) ;
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13
                            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
                            if(zz.lt.0)then;goto 122;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
                            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
                              else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13

!                    a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13    
                    node(ii)%x=node(ent)%x+a
                    node(ii)%y=node(ent)%y+b
                    node(ii)%z=node(ent)%z+c
                    node(ii)%tipus=3 ; node(ii)%icel=jj
                    cels(jj)%node(k)=ii
                end do
            end if

            if(j.gt.2)then                                              ! intermediate points
                dx3=dx2-dx1       ; dy3=dy2-dy1                         ! vectors which link "extreme" points
                dx3=dx3/real(j-1) ; dy3=dy3/real(j-1)                   ! sub-vectors                   
                do k=1,j-2                              
                    ii=ii+1
                    jj=jj+1
                    node(ii)%x=node(ont)%x+(dx1+(real(k)*dx3))
                    node(ii)%y=node(ont)%y+(dy1+(real(k)*dy3))   
                    node(ii)%z=node(ont)%z                                               
!                    print*,"intra i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
                    node(ii)%icel=jj ; node(ii)%tipus=3 ; cels(jj)%node(1)=ii
                    ent=ii
                    if(radi>1)then
                        signo=1                                                        ! miguel 4-6-13
                        do kk=2,radi
                            ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a)
                            xx=der*(1d0-2*a) ;
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13
                            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
                              else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
                       
                            !a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13    
                            node(ii)%x=node(ent)%x+a
                            node(ii)%y=node(ent)%y+b
                            node(ii)%z=node(ent)%z+c
                            node(ii)%tipus=3 ; node(ii)%icel=jj
                            cels(jj)%node(kk)=ii
                        end do
                    end if
                end do
            end if

        end do
    end do
  end do




    !define the cell's centroid and nucleus (margin)            !>>>>>>>>>>>>>> Miquel 14-4-13

    do i=ncelsepi+1,ncels                                            !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%ctipus=3                                    !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0            !>>>>>>>>>>>>>> Miquel 14-4-13
!        print*,"i",i,"node",cels(i)%node(:),"ncels",ncels
        do j=1,cels(i)%nunodes                                !>>>>>>>>>>>>>> Miquel 14-4-13
            k=cels(i)%node(j)                                !>>>>>>>>>>>>>> Miquel 14-4-13
!            print*,"k",k,"nunodes",cels(i)%nunodes
            cels(i)%cex=cels(i)%cex+node(k)%x                !>>>>>>>>>>>>>> Miquel 14-4-13
            cels(i)%cey=cels(i)%cey+node(k)%y                !>>>>>>>>>>>>>> Miquel 14-4-13
            cels(i)%cez=cels(i)%cez+node(k)%z                !>>>>>>>>>>>>>> Miquel 14-4-13
        end do                                                !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
    end do                                                    !>>>>>>>>>>>>>> Miquel 14-4-13

    do i=ncelsepi+1,ncels
          b=1.0d8                                            !>>>>>>>>>>>>>> Is 14-9-13
         do j=1,cels(i)%nunodes                                !>>>>>>>>>>>>>> Is 14-9-13
            k=cels(i)%node(j)                                !>>>>>>>>>>>>>> Is 14-9-13
           a=sqrt((cels(i)%cex-node(k)%x)**2+(cels(i)%cey-node(k)%y)**2+(cels(i)%cez-node(k)%z)**2)
            if (b>a) then ; b=a ; ii=k ; jj=j ; end if
      end do                                                !>>>>>>>>>>>>>> Is 14-9-13
          node(ii)%marge=0

          !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster 
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
          cels(i)%node(jj)=jjj
    end do                                                    !>>>>>>>>>>>>>> Is 14-9-13

    node(:)%talone=0.0d0

end subroutine mesenq1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine epi_hierarchic



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=15    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=0      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=0      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    !allocate(mesradicel(layer))
    !mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=0
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=1d0
    angletor=0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation (epithelial)
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60; 
        node(i)%hoo=5d1
        node(i)%erp=5d0
        node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=0d0
        node(i)%kvol=1d-1
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%rec=3d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.30
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=1d0  !only for epithelium
        node(i)%est=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        !node(i)%fix=1 ;node(i)%rec=1d2
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        !node(i)%fix=1;node(i)%rec=1d2
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      do i=ndepi+k+1,nd
        !node(i)%fix=1;node(i)%rec=1d2
      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=14
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=1d0 !Epithelial-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=1d0 !Effector

      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d-1 !housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=1d0 ; gen(4)%mu=0d0 !FC1 transcript 
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=5

      gen(5)%kindof=4 ; gen(5)%diffu=1.5d0 ; gen(5)%mu=1d0 !FC1
      gen(5)%npre=2 ; allocate(gen(5)%pre(gen(5)%npre))
      gen(5)%pre(1)=4 ; gen(5)%pre(2)=9
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=9

      gen(6)%kindof=2 ; gen(6)%diffu=1d0 ; gen(6)%mu=0d0 !FC2 transcript
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=7

      gen(7)%kindof=4 ; gen(7)%diffu=3d0 ; gen(7)%mu=5d-1 !FC2
      gen(7)%npre=2 ; allocate(gen(7)%pre(gen(7)%npre))
      gen(7)%pre(1)=6 ; gen(7)%pre(2)=11
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=11
      
      gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=1d-1 !R1 receptor
      gen(8)%npost=1 ; allocate(gen(8)%post(gen(8)%npost))
      gen(8)%post(1)=9
      
      gen(9)%kindof=8 ; gen(9)%diffu=0d0 ; gen(9)%mu=0d0 !R1* receptor activated
      gen(9)%npre=2 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=8 ; gen(9)%pre(2)=5
      gen(9)%npost=2 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=8 ; gen(9)%post(2)=5
      
      gen(10)%kindof=2 ; gen(10)%diffu=1d0 ; gen(10)%mu=1d-1 !R2 receptor
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=0d0 !R2* receptor activated
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=10 ; gen(11)%pre(2)=7
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=10 ; gen(11)%post(2)=7

      gen(12)%kindof=1 ; gen(12)%diffu=0d0 ; gen(12)%mu=1d-1 !FT2
      gen(13)%kindof=1 ; gen(13)%diffu=1d0 ; gen(13)%mu=1d-1 !housekeeping gene signalling center (also FT1)
      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=1d-1 !FT3
      
      
    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !epithelial adhesion molecule
       gen(14)%e(5)=-0.05
       gen(14)%e(28)=-1d2

       gen(12)%e(5)=0.05
       gen(12)%e(28)=-1d2


       !gen(14)%e(1)=3  !basal lamina
       !gen(4)%e(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%e(6)=0.10 !this is da (emulating a migratory behavior)
       !gen(4)%e(16)=0.03 !this is dmo (migratory cells)
       !gen(6)%e(nparam_per_node+8)=1d0
       !gen(6)%e(nparam_per_node+16)=1d-2 !this makes random noise biased towards the gradient of the gene
       !gen(9)%e(nparam_per_node+2)=1d-2 !SHH effect on epithelial growth
       !gen(11)%e(nparam_per_node+2)=5d-3 !BMP effect on mesenchymal growth



    !Gene-gene interactions

      !gen(11)%nww=2      !BMP activated receptor induces production of BMP
      !gen(11)%r(1,1)=6  
      !gen(11)%r(1,2)=7
      !gen(11)%r(1,3)=1d1
      !gen(11)%r(2,1)=4   !BMP activated receptor induces production of SHH
      !gen(11)%r(2,2)=5
      !gen(11)%r(2,3)=1d1
      
      !gen(7)%nww=1      !BMP morphogen activates BMP receptor
      !gen(7)%r(1,1)=10  
      !gen(7)%r(1,2)=11
      !gen(7)%r(1,3)=1d0

      gen(13)%t(13)=1d0 !housekeeping central activates itself
      gen(1)%t(13)=1d0  !activates epi cadherin
      gen(4)%t(13)=1d0  !activates FC1 transcript
      gen(8)%t(13)=1d0  !activates ssh receptor
      gen(10)%t(13)=1d0  !activates ssh receptor
      !gen(14)%t(13)=1.5d0  !activates FT3
      gen(12)%t(13)=-1.0d0  !activates FT3

      gen(3)%t(3)=1d0 !housekeeping activates itself
      gen(1)%t(3)=1d0  !activates epi cadherin
      gen(8)%t(3)=1d0  !activates ssh receptor
      gen(10)%t(3)=1d0  !activates ssh receptor

      gen(6)%t(12)=1d0   !FT2 activates FC2 transcript
      gen(14)%t(12)=-1.5d0  !FT2 inhibits FT3

      gen(12)%t(9)=1d0  !R1* activates FT2
      gen(14)%t(11)=1d0  !R2* activates FT3



      gen(13)%nww=1    !housekeeping epi mediates secretion of FC1
      gen(13)%r(1,1)=4
      gen(13)%r(1,2)=5
      gen(13)%r(1,3)=1d1

      gen(3)%nww=1    !housekeeping epi mediates secretion of FC2
      gen(3)%r(1,1)=6
      gen(3)%r(1,2)=7
      gen(3)%r(1,3)=1d1

      
      gen(9)%nww=4      !FC1 morphogen binds R1
      gen(9)%r(1,1)=5  
      gen(9)%r(1,2)=9
      gen(9)%r(1,3)=1d0
      gen(9)%r(2,1)=9  
      gen(9)%r(2,2)=5
      gen(9)%r(2,3)=1d0

      gen(9)%r(3,1)=8  
      gen(9)%r(3,2)=9
      gen(9)%r(3,3)=1d0
      gen(9)%r(4,1)=9  
      gen(9)%r(4,2)=8
      gen(9)%r(4,3)=1d0


      gen(11)%nww=4      !FC2 morphogen activates R2
      gen(11)%r(1,1)=7  
      gen(11)%r(1,2)=11
      gen(11)%r(1,3)=1d0
      gen(11)%r(2,1)=11  
      gen(11)%r(2,2)=7
      gen(11)%r(2,3)=1d0

      gen(11)%r(3,1)=10
      gen(11)%r(3,2)=11
      gen(11)%r(3,3)=1d0
      gen(11)%r(4,1)=11  
      gen(11)%r(4,2)=10
      gen(11)%r(4,3)=1d0



                        !SHH morphogen activates BMP receptor
      !gen(5)%r(3,1)=5  
      !gen(5)%r(3,2)=11
      !gen(5)%r(3,3)=1d0
      !gen(5)%r(4,1)=11  
      !gen(5)%r(4,2)=5
      !gen(5)%r(4,3)=1d0

      !gen(8)%nww=2      !ssh receptor activates SHH receptor
      !gen(8)%r(1,1)=8  
      !gen(8)%r(1,2)=9
      !gen(8)%r(1,3)=1d0
      !gen(8)%r(2,1)=9  
      !gen(8)%r(2,2)=8
      !gen(8)%r(2,3)=1d0

      !gen(10)%nww=2               !bmp receptor activates BMP receptor
      !gen(10)%r(1,1)=10  
      !gen(10)%r(1,2)=11
      !gen(10)%r(1,3)=1d0
      !gen(10)%r(2,1)=11  
      !gen(10)%r(2,2)=10
      !gen(10)%r(2,3)=1d0



    !Adhesion molecules

	ntipusadh=1
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      !kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      !kadh(2,2)=2d1
      !kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d2

    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,nd
      !if(node(i)%fix==1) gex(i,14)=1d0
      !a=sqrt(node(i)%x**2+node(i)%y**2)
      gex(i,1)=1d0
      if(node(i)%icel<=61)then;
        gex(i,4)=0d0; gex(i,13)=1d0;
        gex(i,8)=1d0; gex(i,10)=1d0
      elseif(i<=ndepi)then
          gex(i,3)=1d0
          gex(i,8)=1d0
          gex(i,10)=1d0
      else
        !gex(i,12)=1d0
        !gex(i,10)=1d0
        j=0
        do ii=1,mradicel-2
          j=j+ii
        end do
        j=(6*j+1) !territory 1 (inner ring)
        k=0
        do ii=1,mradicel-1
          k=k+ii
        end do
        k=(6*k+1) !whole layer
        if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
          gex(i,2)=1d0!(3.5-a)/3.5
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
          !gex(i,1)=1d0
        elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
          gex(i,2)=1d0
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
          !gex(i,1)=1d0
        !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
        end if
      end if
    end do

    


    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    
end subroutine

!*********************************************************************************************************

subroutine fig4

call blastula_ensemble

   !******* #4 DEFINING GENETIC PARAMETERS *******
    !Number of genes
    ng=14
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=1d0 !Epithelial-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=1d0 !Effector

      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d-1 !housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=1d0 ; gen(4)%mu=0d0 !FC1 transcript 
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=5

      gen(5)%kindof=4 ; gen(5)%diffu=1.0d-1 ; gen(5)%mu=1d0 !FC1
      gen(5)%npre=2 ; allocate(gen(5)%pre(gen(5)%npre))
      gen(5)%pre(1)=4 ; gen(5)%pre(2)=9
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=9

      gen(6)%kindof=2 ; gen(6)%diffu=1d0 ; gen(6)%mu=0d0 !FC2 transcript
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=7

      gen(7)%kindof=4 ; gen(7)%diffu=1.0d-1 ; gen(7)%mu=1d0 !FC2
      gen(7)%npre=2 ; allocate(gen(7)%pre(gen(7)%npre))
      gen(7)%pre(1)=6 ; gen(7)%pre(2)=11
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=11
      
      gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=1.0d0 !R1 receptor
      gen(8)%npost=1 ; allocate(gen(8)%post(gen(8)%npost))
      gen(8)%post(1)=9
      
      gen(9)%kindof=8 ; gen(9)%diffu=0d0 ; gen(9)%mu=0d0 !R1* receptor activated
      gen(9)%npre=2 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=8 ; gen(9)%pre(2)=5
      gen(9)%npost=2 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=8 ; gen(9)%post(2)=5
      
      gen(10)%kindof=2 ; gen(10)%diffu=1d0 ; gen(10)%mu=1.0d0 !R2 receptor
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=0.0 !R2* receptor activated
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=10 ; gen(11)%pre(2)=7
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=10 ; gen(11)%post(2)=7

      gen(12)%kindof=1 ; gen(12)%diffu=0.0 ; gen(12)%mu=3d-2 !1d-1 !FT2
      gen(13)%kindof=1 ; gen(13)%diffu=0.0 ; gen(13)%mu=1d-1 !housekeeping gene signalling center (also FT1)
      gen(14)%kindof=1 ; gen(14)%diffu=0.0 ; gen(14)%mu=3d-2 !1d-1 !FT3
      
      
    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !epithelial adhesion molecule
       gen(14)%e(21)=-0.02 !-0.05
       gen(14)%e(28)=-1d2

       gen(12)%e(21)=0.02 !0.05
       gen(12)%e(28)=-1d2


       !gen(14)%e(1)=3  !basal lamina
       !gen(4)%e(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%e(6)=0.10 !this is da (emulating a migratory behavior)
       !gen(4)%e(16)=0.03 !this is dmo (migratory cells)
       !gen(6)%e(nparam_per_node+8)=1d0
       !gen(6)%e(nparam_per_node+16)=1d-2 !this makes random noise biased towards the gradient of the gene
       !gen(9)%e(nparam_per_node+2)=1d-2 !SHH effect on epithelial growth
       !gen(11)%e(nparam_per_node+2)=5d-3 !BMP effect on mesenchymal growth

    !Gene-gene interactions

      !gen(11)%nww=2      !BMP activated receptor induces production of BMP
      !gen(11)%r(1,1)=6  
      !gen(11)%r(1,2)=7
      !gen(11)%r(1,3)=1d1
      !gen(11)%r(2,1)=4   !BMP activated receptor induces production of SHH
      !gen(11)%r(2,2)=5
      !gen(11)%r(2,3)=1d1
      
      gen(12)%nww=1      !FT2 activates FC2
      gen(12)%r(1,1)=6  
      gen(12)%r(1,2)=7
      gen(12)%r(1,3)=1d3

      gen(13)%t(13)=1d0 !housekeeping central activates itself
      gen(1)%t(13)=1d0  !activates epi cadherin
      gen(4)%t(13)=1d0  !activates FC1 transcript
      gen(8)%t(8)=1d0  !R1 activates its own transcription
      gen(12)%t(13)=-1d3 !FT1 inhibits FT2
      gen(10)%t(10)=1d0  ! R2 activates its own trans
      gen(8)%t(8)=1d0  !R1 activates its own transcription
      !gen(14)%t(13)=1.5d0  !activates FT3
      gen(12)%t(13)=-1.0d5  !activates FT3
!      gen(12)%t(12)=1.500d1  !FT3 activates itself

      gen(3)%t(3)=1d0 !housekeeping activates itself
      gen(1)%t(3)=1d0  !activates epi cadherin
      !gen(8)%t(3)=1d0  !activates FC1 receptor
      gen(10)%t(12)=1d0  !activates ssh receptor

      gen(6)%t(12)=1d4   !FT2 activates FC2 transcript
      gen(14)%t(12)=-1.5d0  !FT2 inhibits FT3
      gen(14)%t(14)=0.50d0  !FT3 activates itself

      gen(12)%t(9)=1.0d4  !R1* activates FT2
      gen(14)%t(11)=0.1d0  !R2* activates FT3

      gen(13)%nww=1    !housekeeping epi mediates secretion of FC1
      gen(13)%r(1,1)=4
      gen(13)%r(1,2)=5
      gen(13)%r(1,3)=1d1

      gen(3)%nww=1    !housekeeping epi mediates secretion of FC2
      gen(3)%r(1,1)=6
      gen(3)%r(1,2)=7
      gen(3)%r(1,3)=1d1

      gen(9)%nww=4      !FC1 morphogen binds R1
      gen(9)%r(1,1)=5  
      gen(9)%r(1,2)=9
      gen(9)%r(1,3)=1d0
      gen(9)%r(2,1)=9  
      gen(9)%r(2,2)=5
      gen(9)%r(2,3)=1d0

      gen(9)%r(3,1)=8  
      gen(9)%r(3,2)=9
      gen(9)%r(3,3)=1d0
      gen(9)%r(4,1)=9  
      gen(9)%r(4,2)=8
      gen(9)%r(4,3)=1d0


      gen(11)%nww=4      !FC2 morphogen activates R2
      gen(11)%r(1,1)=7  
      gen(11)%r(1,2)=11
      gen(11)%r(1,3)=1d3
      gen(11)%r(2,1)=11  
      gen(11)%r(2,2)=7
      gen(11)%r(2,3)=1d0

      gen(11)%r(3,1)=10
      gen(11)%r(3,2)=11
      gen(11)%r(3,3)=1d0
      gen(11)%r(4,1)=11  
      gen(11)%r(4,2)=10
      gen(11)%r(4,3)=1d0



                        !SHH morphogen activates BMP receptor
      !gen(5)%r(3,1)=5  
      !gen(5)%r(3,2)=11
      !gen(5)%r(3,3)=1d0
      !gen(5)%r(4,1)=11  
      !gen(5)%r(4,2)=5
      !gen(5)%r(4,3)=1d0

      !gen(8)%nww=2      !ssh receptor activates SHH receptor
      !gen(8)%r(1,1)=8  
      !gen(8)%r(1,2)=9
      !gen(8)%r(1,3)=1d0
      !gen(8)%r(2,1)=9  
      !gen(8)%r(2,2)=8
      !gen(8)%r(2,3)=1d0

      !gen(10)%nww=2               !bmp receptor activates BMP receptor
      !gen(10)%r(1,1)=10  
      !gen(10)%r(1,2)=11
      !gen(10)%r(1,3)=1d0
      !gen(10)%r(2,1)=11  
      !gen(10)%r(2,2)=10
      !gen(10)%r(2,3)=1d0



    !Adhesion molecules

	ntipusadh=1
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      !kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      !kadh(2,2)=2d1
      !kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d2

    end if

    gen(8)%e(nparam_per_node+2)=5d-2 !0
    gen(8)%e(nparam_per_node+2)=5d-2 !5d-5

    do i=1,cels(1)%nunodes
      ii=cels(1)%node(i)
      gex(ii,13)=1.0d0
    end do
    gex(:,8)=1.d0
    gex(:,10)=1.0d0

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    ffu(1)=0
    ffu(11)=1

    dif_req=1.20d2

end subroutine

!*****************************************************************************************************

subroutine fig4_pre ! it makes a blastula made of cells made of a single cylinder

    call blastula_ensemble

    ng=14
    call initiate_gene

    gen(8)%e(nparam_per_node+2)=4d-3 !0
    gex(:,8)=1.0d0
    gex(:,10)=1.0d0
    gen(8)%diffu=1.0d0
    gen(10)%diffu=1.0d0
    gen(8)%mu=1.0d0
    gen(10)%mu=1.0d0
    gen(8)%t(8)=1.1d0
    gen(10)%t(10)=1.1d0
    gen(8)%kindof=2
    gen(10)%kindof=2

    do i=1,ng
      gen(i)%e(23)=0.0
    end do

    node(:)%add=node(:)%add*1.2
    node(:)%erp=node(:)%erp!*10.0
    node(:)%est=node(:)%est*1000.0
    cels(:)%maxsize_for_div=24 !cels(:)%nunodes*2

    ncels=nd/2
    ncals=ncels+10
print *,nd/2,ncels,"nd/2 ncels"
    if (allocated(cels)) deallocate(cels)
    allocate(cels(ncals))
    k=0
    do i=1,nd
      if (node(i)%tipus==1) then
        k=k+1
        node(i)%icel=k
        if (allocated(cels(k)%node)) deallocate(cels(k)%node)
        allocate(cels(k)%node(2))
        cels(k)%nunodes=2
        cels(k)%node(1)=i
        cels(k)%node(2)=node(i)%altre
        node(i)%marge=0
      end if
    end do
    do i=1,nd
      if (node(i)%tipus==2) then
        node(i)%icel=node(node(i)%altre)%icel
        node(i)%marge=1
      end if
    end do
    cels(:)%ctipus=1
    cels(:)%nodela=2
!do i=1,nd
!  print *,i,"i",node(i)%icel
!end do
!do i=1,ncels
!  print *,i,cels(i)%node
!end do

    min_comp=100

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    ffu(1)=0
    ffu(11)=1

    dif_req=1.20d2

end subroutine

!***************************************************************************************************************
!***************************************************************************************************************

subroutine father ! it makes a blastula made of cells made of a single cylinder

    call blastula_ensemble

  if (allocated(px)) deallocate(px)
  if (allocated(py)) deallocate(py)
  if (allocated(pz)) deallocate(pz)
  if (allocated(dex)) deallocate(dex)
  allocate(px(nda),py(nda),pz(nda),dex(nda))
  px=0 ; py=0 ; pz=0 ; dex=0

   
  if (allocated(vcilx)) deallocate(vcilx)
  if (allocated(vcily)) deallocate(vcily)
  if (allocated(vcilz)) deallocate(vcilz)
  allocate(vcilx(nda),vcily(nda),vcilz(nda))
  vcilx=0 ; vcily=0 ; vcilz=0

  if (allocated(vtorx)) deallocate(vtorx)
  if (allocated(vtory)) deallocate(vtory)
  if (allocated(vtorz)) deallocate(vtorz)
  allocate(vtorx(nda),vtory(nda),vtorz(nda))
  vtorx=0 ; vtory=0 ; vtorz=0
  if (allocated(vstorx)) deallocate(vstorx)
  if (allocated(vstory)) deallocate(vstory)
  if (allocated(vstorz)) deallocate(vstorz)
  allocate(vstorx(nda),vstory(nda),vstorz(nda))
  vstorx=0 ; vstory=0 ; vstorz=0
  if (allocated(vsprx)) deallocate(vsprx)
  if (allocated(vspry)) deallocate(vspry)
  if (allocated(vsprz)) deallocate(vsprz)
  allocate(vsprx(nda),vspry(nda),vsprz(nda))
  vsprx=0 ; vspry=0 ; vsprz=0

  if (allocated(erep)) deallocate(erep)
  if (allocated(erepcel)) deallocate(erepcel)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(eyou)) deallocate(eyou)
  if (allocated(espring)) deallocate(espring)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(etor)) deallocate(etor)
print *,nda,"nda"
  allocate(erep(nda),erepcel(nda),eadh(nda),eyou(nda),espring(nda),etor(nda))

  if (allocated(fmeanl)) deallocate(fmeanl)  !>>Miquel23-1-14
  allocate(fmeanl(nda))                     !>>Miquel23-1-14
  fmeanl=0                                !>>Miquel23-1-14
   
  if (allocated(fmeanv)) deallocate(fmeanv)  !>>Miquel23-1-14
  allocate(fmeanv(nda))                     !>>Miquel23-1-14
  fmeanv=0                                !>>Miquel23-1-14

    gen(5)%e(25)=1.7d-1
!    gen(5)%e(nparam_per_node+2)=4.0d-2
    gen(5)%e(nparam_per_node+2)=5.0d-1
!    gen(1)%e(nparam_per_node+2)=1.0d-5
    gen(4)%e(nparam_per_node+8)=4.0d1
    gen(5)%e(nparam_per_node+11)=4.0d50
    !gen(2)%e(23)=1.0d-1

    node(:)%add=node(:)%add*1.2
!    node(:)%erp=node(:)%est*1.0d-1
    node(:)%est=node(:)%est*1000.0
    node(:)%est=node(:)%erp*1.0d1
    cels(:)%maxsize_for_div=24 !cels(:)%nunodes*2

    ncels=nd/2
    ncals=ncels+10
    if (allocated(cels)) deallocate(cels)
    allocate(cels(ncals))
    k=0
    do i=1,nd
      if (node(i)%tipus==1) then
        k=k+1
        node(i)%icel=k
        if (allocated(cels(k)%node)) deallocate(cels(k)%node)
        allocate(cels(k)%node(2))
        cels(k)%nunodes=2
        cels(k)%node(1)=i
        cels(k)%node(2)=node(i)%altre
        node(i)%marge=1
      end if
    end do
    do i=1,nd
      if (node(i)%tipus==2) then
        node(i)%icel=node(node(i)%altre)%icel
        node(i)%marge=0
      end if
    end do
    cels(:)%ctipus=1
    cels(:)%nodela=2

    min_comp=100

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

    ffu(1)=0
    ffu(11)=1
    ffu(7)=0
    ffu(14)=1
!    ffu(12)=0 !ACHTUNG

!    node(:)%pla=1.0d0

    deltamin=1d-4

!    dif_req=1.20d2

    ic_load=9

    prop_noise=0.05d-1
!    prop_noise=0.0d0

print *,ng,"ng"

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine father_inva ! it makes a blastula made of cells made of a single cylinder

    call blastula_ensemble

  if (allocated(px)) deallocate(px)
  if (allocated(py)) deallocate(py)
  if (allocated(pz)) deallocate(pz)
  if (allocated(dex)) deallocate(dex)
  allocate(px(nda),py(nda),pz(nda),dex(nda))
  px=0 ; py=0 ; pz=0 ; dex=0

   
  if (allocated(vcilx)) deallocate(vcilx)
  if (allocated(vcily)) deallocate(vcily)
  if (allocated(vcilz)) deallocate(vcilz)
  allocate(vcilx(nda),vcily(nda),vcilz(nda))
  vcilx=0 ; vcily=0 ; vcilz=0

  if (allocated(vtorx)) deallocate(vtorx)
  if (allocated(vtory)) deallocate(vtory)
  if (allocated(vtorz)) deallocate(vtorz)
  allocate(vtorx(nda),vtory(nda),vtorz(nda))
  vtorx=0 ; vtory=0 ; vtorz=0
  if (allocated(vstorx)) deallocate(vstorx)
  if (allocated(vstory)) deallocate(vstory)
  if (allocated(vstorz)) deallocate(vstorz)
  allocate(vstorx(nda),vstory(nda),vstorz(nda))
  vstorx=0 ; vstory=0 ; vstorz=0
  if (allocated(vsprx)) deallocate(vsprx)
  if (allocated(vspry)) deallocate(vspry)
  if (allocated(vsprz)) deallocate(vsprz)
  allocate(vsprx(nda),vspry(nda),vsprz(nda))
  vsprx=0 ; vspry=0 ; vsprz=0

  if (allocated(erep)) deallocate(erep)
  if (allocated(erepcel)) deallocate(erepcel)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(eyou)) deallocate(eyou)
  if (allocated(espring)) deallocate(espring)
  if (allocated(eadh)) deallocate(eadh)
  if (allocated(etor)) deallocate(etor)
print *,nda,"nda"
  allocate(erep(nda),erepcel(nda),eadh(nda),eyou(nda),espring(nda),etor(nda))

  if (allocated(fmeanl)) deallocate(fmeanl)  !>>Miquel23-1-14
  allocate(fmeanl(nda))                     !>>Miquel23-1-14
  fmeanl=0                                !>>Miquel23-1-14
   
  if (allocated(fmeanv)) deallocate(fmeanv)  !>>Miquel23-1-14
  allocate(fmeanv(nda))                     !>>Miquel23-1-14
  fmeanv=0                                !>>Miquel23-1-14

    gen(5)%e(25)=0.0d0 !1.7d-1
    !gen(5)%e(nparam_per_node+2)=2.0d-2
!    gen(5)%e(nparam_per_node+2)=1.4d0
    gen(5)%e(nparam_per_node+2)=2.5d-1
!    gen(1)%e(nparam_per_node+2)=1.0d-5
!    gen(4)%e(nparam_per_node+8)=1d10 !4.0d1
    gen(4)%e(nparam_per_node+8)=4.0d1
!    gen(5)%e(nparam_per_node+11)=4.0d50
    gen(6)%e(nparam_per_node+1)=5.0d1
    gen(5)%e(21)=-4.0d-4
!!    gen(4)%e(21)=-9.0d-1  !!!
    gen(4)%e(28)=1.0d1
    
    gen(7)%t(7)=1.0d0 !this is to make that kplast is not high from the beginning
    gen(7)%diffu=1.0d0
    gen(7)%mu=0.1d0
!    gen(7)%e(27)=1.0d-1 !bo
    gen(7)%e(27)=4.0d-1 !bo
!    gen(7)%e(27)=0.0d0
    gen(7)%kindof=1
    gex(:,7)=1.0d-2

    !gen(2)%e(23)=1.0d-1;

    node(:)%add=node(:)%add*1.2
!    node(:)%erp=node(:)%est*1.0d-1
!    node(:)%est=node(:)%est*1000.0
    node(:)%est=node(:)%erp*1.0d3
!    node(:)%erp=node(:)%erp*5.0d-1
    node(:)%erp=node(:)%erp*9.0d-1
!    node(:)%est=node(:)%erp*1.0d0
!    node(:)%est=node(:)%erp!*1.0d-3
!    cels(:)%maxsize_for_div=24 !cels(:)%nunodes*2

    ncels=nd/2
    ncals=ncels+10
    if (allocated(cels)) deallocate(cels)
    allocate(cels(ncals))
    k=0
    do i=1,nd
      if (node(i)%tipus==1) then
        k=k+1
        node(i)%icel=k
        if (allocated(cels(k)%node)) deallocate(cels(k)%node)
        allocate(cels(k)%node(2))
        cels(k)%nunodes=2
        cels(k)%node(1)=i
        cels(k)%node(2)=node(i)%altre
        node(i)%marge=1
        node(i)%eqd=node(node(i)%altre)%eqd 
        node(i)%add=node(node(i)%altre)%add
        node(i)%x=node(i)%x*1.5d0
        node(i)%y=node(i)%y*1.5d0
        node(i)%z=node(i)%z*1.5d0
        node(i)%grd=node(node(i)%altre)%grd*0.7d0
        node(i)%cod=node(node(i)%altre)%cod*0.7d0
      else
        node(i)%x=node(i)%x*1.5d0 
        node(i)%y=node(i)%y*1.5d0
        node(i)%z=node(i)%z*1.5d0
!        node(i)%grd=node(node(i)%altre)%grd 
!        node(i)%cod=node(node(i)%altre)%cod
      end if
    end do
    do i=1,nd
      if (node(i)%tipus==2) then
        node(i)%icel=node(node(i)%altre)%icel
        node(i)%marge=0
        node(i)%pla=0.0d0
      end if
    end do
    cels(:)%ctipus=1
    cels(:)%nodela=2

    min_comp=100

    gen(6)%mu=0
    gen(6)%diffu=0
    gex(:,6)=1.0d0

    do i=1,ncels
      cels(i)%maxsize_for_div=node(cels(i)%node(1))%eqd*0.8d0
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3

!    resmax=1d-3

ffu=0
    ffu(1)=0
    ffu(2)=0 !>>> TT no screening
    ffu(11)=1
!    ffu(10)=0 ! volume conservation
    ffu(14)=1
    ffu(6)=1
    ffu(17)=1
    ffu(18)=0
!    node(:)%pla=1.0d-1

    deltamin=1d-6

!    dif_req=1.20d2

    df_reqmax=1d10

    ic_load=10

    reqmin=1.0d-4

!    desmax=1d-5
!    do i=1,nd
!      node(i)%dmo=desmax
!    end do

    prop_noise=0.0!1d-5
!    prop_noise=0.0d0

print *,ng,"ng"

end subroutine

!!!!!!!!!***********************************SUBROUTINE********************************

subroutine test_RD

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=0       !number of radial layers of nodes per cell
    radicel=0    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=20   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.0d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 ;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

    if(radi==1.or.mradi==1) ffu(1)=0
   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d0 ! miguel 14-10-13   
    deltamin=1d-1
    dmax=1
    screen_radius=1.0d0
    khold=1d0
    angletor=0.05

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
   ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation
    ffu(15)=1 !0 = unbiased random noise / 1 = noise biased by energies


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%rec=1d1
        end if
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.70; 
        node(i)%hoo=1d1
        node(i)%erp=5d1
        node(i)%est=5d0
        node(i)%mov=temp
        node(i)%dmo=0.00
        node(i)%dif=0d0
        node(i)%pla=1d0
        node(i)%kvol=1d0
        node(i)%kfi=khold
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%rec=1d1
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.70
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=0d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.00
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%eqd,node(1)%eqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%eqd) ; end if ; end if



    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%fix=0

    do i=1,nd
     node(i)%fix=2
    end do

    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
    !    node(i)%fix=1
    !    !node(i)%add=node(i)%eqd*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%rec=1d1
    !    !node(i)%hoo=1d1
    !    !node(i)%erp=1d1
    !    !!node(i)%est=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        !node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
        !node(i)%fix=1
        !node(i)%add=node(i)%eqd*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%rec=1d1
      end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%fix=1
      !  !node(i)%add=node(i)%eqd*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%rec=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   do i=1,ncels
     call random_number(a)
     cels(i)%fase=a
   end do



   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=13
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d0 ; gen(1)%mu=0d0 !E-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=1d0 ; gen(2)%mu=5d-1 !P-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 !housekeeping

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=1d-3 !activator factor

      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-3 !activator transcript
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=7

      gen(6)%kindof=2 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-3 !inhibitor transcript
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=8


      gen(7)%kindof=4 ; gen(7)%diffu=4.0d-1 ; gen(7)%mu=5d-1 !activator diffusible form
      gen(7)%npre=1 ; allocate(gen(7)%pre(gen(7)%npre))
      gen(7)%pre(1)=5

      gen(8)%kindof=4 ; gen(8)%diffu=2d0 ; gen(8)%mu=5d-1!inhibitor diffusible form
      gen(8)%npre=1 ; allocate(gen(8)%pre(gen(8)%npre))
      gen(8)%pre(1)=6


      gen(9)%kindof=2 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 !activator receptor unbound
      !gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre))
      !gen(9)%pre(1)=6

      gen(10)%kindof=2 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 !inhibitor receptor unbound
      !gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      !gen(10)%pre(1)=6

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=1d-3 !activator receptor complex
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre)) ; gen(11)%pre(1)=7 ; gen(11)%pre(2)=9
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost)) ; gen(11)%post(1)=7 ; gen(11)%post(2)=9

      gen(12)%kindof=8 ; gen(12)%diffu=0d0 ; gen(12)%mu=1d-3 !inhibitor receptor complex
      gen(12)%npre=2 ; allocate(gen(12)%pre(gen(12)%npre)) ; gen(12)%pre(1)=8 ; gen(12)%pre(2)=10
      gen(12)%npost=2 ; allocate(gen(12)%post(gen(12)%npost)) ; gen(12)%post(1)=8 ; gen(12)%post(2)=10

      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=1d-3 !inhibitor factor
      
    !Gene-behavior interactions
    

       gen(1)%e(1)=1  !epithelial adhesion molecule
       gen(2)%e(1)=2  !epithelial-mesenchymal adhesion molecule
       !gen(3)%e(1)=3  !basal lamina
       !gen(4)%e(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%e(6)=0.10 !this is da (emulating a migratory behavior)
       
       !migration
       !gen(4)%e(16)=2d-1 !this is dmo (migratory cells)
       !gen(6)%e(nparam_per_node+8)=1d0
       !gen(4)%e(nparam_per_node+16)=1d0 !this makes random noise biased towards the gradient of the gene

       !proliferation
       !gen(4)%e(nparam_per_node+2)=1d-3 !no division for placodal cells
       !gen(13)%e(nparam_per_node+2)=1d-3 !no division for placodal cells



    !Gene-gene interactions


      gen(3)%t(3)=1d0
      !gen(1)%t(3)=1d0

!!!!!NO RD NETWORK
      !gen(4)%t(4)=1d1
      !gen(1)%t(4)=1d0
      !gen(13)%t(13)=1d-1
      !gen(2)%t(13)=1d-2

      !gen(4)%t(13)=-1d1

!!!!!RD NETWORK
      gen(9)%t(3)=1.0d0
      gen(10)%t(3)=1d0
      !gen(5)%t(3)=1d-3


      !gen(4)%t(4)=1d-2
      gen(5)%t(4)=1d-2
      gen(6)%t(4)=1d-2

      gen(13)%t(12)=1d4

      gen(4)%t(11)=1.3d0
      gen(4)%t(12)=-1.0d0

      !gen(4)%t(13)=-1d0

      gen(3)%nww=2
      gen(3)%r(1,1)=5
      gen(3)%r(1,2)=7
      gen(3)%r(1,3)=1d0
      gen(3)%r(2,1)=6
      gen(3)%r(2,2)=8
      gen(3)%r(2,3)=1d0

      gen(11)%nww=4
      gen(11)%r(1,1)=7
      gen(11)%r(1,2)=11
      gen(11)%r(1,3)=1d0
      gen(11)%r(2,1)=11
      gen(11)%r(2,2)=7
      gen(11)%r(2,3)=1d0
      gen(11)%r(3,1)=9
      gen(11)%r(3,2)=11
      gen(11)%r(3,3)=1d0
      gen(11)%r(4,1)=11
      gen(11)%r(4,2)=9
      gen(11)%r(4,3)=1d0

      gen(12)%nww=4
      gen(12)%r(1,1)=8
      gen(12)%r(1,2)=12
      gen(12)%r(1,3)=1d0
      gen(12)%r(2,1)=12
      gen(12)%r(2,2)=8
      gen(12)%r(2,3)=1d0
      gen(12)%r(3,1)=10
      gen(12)%r(3,2)=12
      gen(12)%r(3,3)=1d0
      gen(12)%r(4,1)=12
      gen(12)%r(4,2)=10
      gen(12)%r(4,3)=1d0
!!!!!!!!!!!!!!!!!!

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=3d1
      kadh(1,2)=3d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d1
      !kadh(1,3)=-1d1 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=-1d1 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d1
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-3
      j=j+ii
    end do
    j=(6*j+1) ;print*,"jota",j

    !gex(1,4)=1d-2
    !gex(1,13)=1d-1
    do i=1,nd
       !gex(i,1)=1d0

       !if(node(i)%icel<128.or.(i>ndepi.and.node(i)%icel<ncelsepi+128))then
       !  gex(i,1)=1d0
       !else
       !  gex(i,2)=1d0
       !end if
       gex(i,1)=1d0
       !if(node(i)%tipus==1) cycle
       gex(i,3)=1d0
       gex(i,9)=1d0
       gex(i,10)=1d0
       !gex(i,13)=5d-2

       call random_number(a)
       gex(i,4)=a*0.01d0

   !if(i==1) gex(i,4)=0.1d0

       !gex(i,9)=1d-1
       !gex(i,10)=1d-1

       !call random_number(a)
       !gex(i,4)=a*1d-2

      !a=sqrt(node(i)%x**2+node(i)%y**2)
      if(i<=ndepi)then
      !  if(node(i)%icel<=j)then
          !gex(i,3)=1d-1
          !gex(i,9)=1d-1
          !gex(i,10)=1d-1
        !else
          !gex(i,1)=1d0
      !  end if
        !if(node(i)%tipus==2)then
        !  gex(i,3)=1d0
        !  gex(i,1)=0d0 ; gex(i,2)=0d0
        !end if
        !if(node(i)%icel<128)then
        !  !gex(i,1)=1d0
        !  gex(i,4)=1d-1
        !else
        !  !gex(i,2)=1d0
        !  gex(i,13)=1d-1
        !end if

      else
        !if(node(i)%icel<ncelsepi+128)then
        !  !gex(i,1)=1d0
        !  gex(i,4)=1d-1
        !else
        !  !gex(i,2)=1d0
        !  gex(i,13)=1d-1
        !end if
      !  j=0
      !  do ii=1,mradicel-3
      !    j=j+ii
      !  end do
      !  j=(6*j+1) !territory 1 (inner ring)
      !  k=0
      !  do ii=1,mradicel-1
      !    k=k+ii
      !  end do
      !  k=(6*k+1) !whole layer
      !  if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
      !    gex(i,2)=1d0!(3.5-a)/3.5
      !    gex(i,4)=1d0
      !    !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
      !  elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
      !    gex(i,1)=1d0
      !  elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
      !    gex(i,2)=1d0
      !    gex(i,4)=1d0
      !    !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
      !  elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
      !    gex(i,1)=1d0
      !  !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
      !  end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%add)*3
   
end subroutine

subroutine epidermal_rectangle

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=0       !number of radial layers of nodes per cell
    radicel=0    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=0      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=0      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

    
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************

    !if(mod(ly,2)==0)then;
    !  ncelsepi=(ly/2)*(lx+lx-1)
    !else
    !  ncelsepi=((ly-1)/2)*(lx+lx-1)+lx-1
    !endif
    ncelsepi=lx*ly

    !ncelsepi=(lx-1)*(ly-1)
    !if(mod(ly,2)==0)then; ncelsepi=ncelsepi+lx-2 ;else;ncelsepi=ncelsepi+lx; endif
    ndepi=ncelsepi*2
    ncelsmes=ncelsepi*layer
    ndmes=ncelsmes
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nodecel=2 ; nodecela=5
    

    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=0
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=2d0
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=0
    ffu(2)=0 !to quit if there are too many cells
    ffu(2)=0 !>>> TT no screening !screening
    ffu(3)=0 !torsion
    

    ffu(6)=1 !epithelial node plastic deformation
    ffu(7)=1 !0 dynamic delta, 1 fixed delta
    ffu(8)=0 !neighboring by triangulation
    
    ffu(10)=0 !volume conservation (epithelial)
    ffu(15)=1 !1 = noise biased by energies , 0 = unbiased noise
    ffu(17)=0


 

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%rec=5d0
          !node(i)%pld=-0.05d0
          !node(i)%eqd=0.20d0
        !node(i)%pla=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%rec=5d0
        !node(i)%pla=1d1
        !  !node(i)%pld=0.05d0
        !  !node(i)%eqd=0.30d0
        end if
        node(i)%eqd=0.25d0
        node(i)%grd=0.25d0 !; node(i)%grd=node(i)%eqd
        node(i)%eqs=0.25
        node(i)%add=node(i)%eqd*1.60; 
        node(i)%hoo=5d0
        node(i)%erp=3d1
        node(i)%est=1d1
        node(i)%mov=temp
        node(i)%dmo=0
        node(i)%dif=0d0
        node(i)%kfi=khold
        node(i)%pla=1d1
        node(i)%kvol=1d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=5d0 ; node(i)%rec=5d0
        node(i)%eqd=0.25d0 ; node(i)%grd=node(i)%eqd
        node(i)%eqs=0d0
        node(i)%add=node(i)%eqd*1.60
        node(i)%hoo=0d0   !only for epithelium
        node(i)%erp=0d0  !only for epithelium
        node(i)%est=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mov=temp
        node(i)%dmo=0.0
        node(i)%dif=0d0
        node(i)%kfi=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%add)

    node(:)%fix=0

	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do

    !epithelial cells
    !origin of rectangle
    aa=node(1)%eqd*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%eqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%eqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%fix=1;node(ii-1)%fix=1 ;
          node(ii)%rec=1d2;node(ii-1)%rec=1d2;
        end if
      end do
    end do

    !mesenchymals
    if(layer>0)then
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%fix=2 ;node(ii)%rec=1d2 ; 
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%fix=1 ;node(ii)%rec=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%fix=1 ;node(ii)%rec=1d2 ;end if
        end do
      end do
    end do
    end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%fix=2
  !node(:)%fix=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%eqd
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=15
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !
      
      gen(2)%kindof=4 ; gen(2)%diffu=5d-1 ; gen(2)%mu=1.0d-1 !activator signal
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost))
      gen(2)%post(1)=4
      
      gen(3)%kindof=2 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 !activator receptor inactive epithelial
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !activator receptor active epithelial
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 !activator receptor inactive mesench
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !activator receptor active mesench
      gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 !

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 !

      gen(9)%kindof=4 ; gen(9)%diffu=2.5d1 ; gen(9)%mu=5.0d-1 !inhibitor signal, epithelial
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=11
      
      gen(10)%kindof=2 ; gen(10)%diffu=0d0 ; gen(10)%mu=1.0d0 !inhibitor receptor inactive epithelial
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !inhibitor receptor active epithelial
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 !epithelial adhesion molecule
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 !mesench adhesion molecule

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 !
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  !
      
      
      
    !Gene-behavior interactions
    

      !gen(7)%e(1)=1
      !gen(8)%e(1)=2
      gen(12)%e(1)=1
      gen(13)%e(1)=2

      gen(12)%e(nparam_per_node+2)=1d-1

    !Gene-gene interactions

     !wavefront setting
      !gen(2)%t(1)=1d2  !this will spontaneously generate activator
      !gen(15)%t(1)=-1d1  !this will spontaneously generate activator
      !gen(1)%t(14)=1d0  !this will produce a ware
      !!gen(1)%t(15)=-1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(14)%t(14)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(15)%t(15)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(2)%t(15)=-1d5  !this is everywhere, inhibiting production of activator (in a threshold manner)



      !gen(2)%t(1)=-1d5
      !gen(2)%t(1)=-1d5

      !gen(2)%t(2)=1d0  !the activated receptor directly transcripts signal
      !gen(9)%t(2)=1d0

      !gen(7)%t(2)=1d0
      !gen(8)%t(2)=1d0
      !gen(7)%t(3)=-6d-2
      !gen(8)%t(3)=-1d3
      !gen(7)%t(5)=-1d3
      !gen(8)%t(5)=-6d-2
      !gen(14)%t(7)=-1d2
      !gen(15)%t(8)=-1d2
      
      !gen(2)%t(9)=-1d0
      
      gen(3)%t(3)=1d0  !autoactivation of receptors
      gen(5)%t(5)=1d0
      gen(10)%t(10)=1d0
      
      gen(12)%t(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%t(13)=1d0
      !gen(14)%t(3)=1d0
      !gen(15)%t(5)=1d0

      
      !gen(4)%nww=4      !activator signal activates activator receptor epi
      !gen(4)%r(1,1)=2  
      !gen(4)%r(1,2)=4
      !gen(4)%r(1,3)=1d0
      !gen(4)%r(2,1)=4  
      !gen(4)%r(2,2)=2
      !gen(4)%r(2,3)=1d0
      !
      !gen(4)%r(3,1)=3  
      !gen(4)%r(3,2)=4
      !gen(4)%r(3,3)=1d0
      !gen(4)%r(4,1)=4  
      !gen(4)%r(4,2)=3
      !gen(4)%r(4,3)=1d0
      !
      !
      !gen(6)%nww=4      !activator morphogen activates activator receptor mesench
      !gen(6)%r(1,1)=2  
      !gen(6)%r(1,2)=6
      !gen(6)%r(1,3)=1d0
      !gen(6)%r(2,1)=6  
      !gen(6)%r(2,2)=2
      !gen(6)%r(2,3)=1d0
      !
      !gen(6)%r(3,1)=5
      !gen(6)%r(3,2)=6
      !gen(6)%r(3,3)=1d0
      !gen(6)%r(4,1)=6  
      !gen(6)%r(4,2)=5
      !gen(6)%r(4,3)=1d0
      !
      !gen(11)%nww=4      !inhibitor signal activates inhibitor receptor epi
      !gen(11)%r(1,1)=9  
      !gen(11)%r(1,2)=11
      !gen(11)%r(1,3)=1d0
      !gen(11)%r(2,1)=11  
      !gen(11)%r(2,2)=9
      !gen(11)%r(2,3)=1d0
      !
      !gen(11)%r(3,1)=10  
      !gen(11)%r(3,2)=11
      !gen(11)%r(3,3)=1d0
      !gen(11)%r(4,1)=11  
      !gen(11)%r(4,2)=10
      !gen(11)%r(4,3)=1d0


    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d0
      kadh(1,2)=3d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0
      !kadh(3,3)=5d0
      !kadh(3,4)=2d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,4)=5d0
    end if

    !Gene expression on nodes

    
    
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          !node(k)%fix=3
          !node(node(k)%altre)%fix=3
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            node(k)%fix=3
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!


   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
          gex(i,3)=1d0  !receptors
          gex(i,10)=1d0
          gex(i,14)=1d0
          call random_number(a)
          gex(i,2)=a!*0.1d0
          gex(i,9)=a!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0
      end if
      if(node(i)%tipus==3)then
        gex(i,5)=1d0
        gex(i,13)=1d0
        gex(i,15)=1d0
      end if
      !if(node(i)%fix>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%add)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module
