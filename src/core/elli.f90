! TO DO:

!-sONERA
!-Revertir lo del email de nordea
!-Backup
!-see if artefact in genetic
!-add ecm and mesen
!-explicar al text que posem un limit a la quantitat difosa
!-explicar al text que la difusio es fa amb neigh
!-explicar que la divisio no fa dilucio
!-add nodeo to input files online
!-chose node and cell by cursor
!-chose with mouse color ranges dynamically
!-3d plot on epithelium of properties
!-add ldi to model parameters

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



!*****************************************
!            PROGRAMA
!*****************************************
program veu

!use opengl
use opengl_gl
use opengl_glu
use opengl_glut
use opengl_kinds
use view_modifier
use function_plotter
use inicial
use automaticon

implicit none
character*400 caordre,ordre !!>>HC 18-11-2020 these should be a bit longer for evolution
character*200 precarg       !!>>HC 18-11-2020 (evolve calls with complete paths)
character*80 cu
character*10 pgm,xcad

integer :: winid, submenuid, mode=ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH))

integer*4 iu,id

integer*8 iv
type(c_ptr) point

!this is necessary to accomadate to different versions of linux
call random_seed(size = nseed)
print*,nseed
call getarg(2,cu)
if (len_trim(cu)==1) then
  read (cu,*) aut
else
  if(allocated(idum))deallocate(idum)
  allocate(idum(nseed))
  if(allocated(idumoriginal))deallocate(idumoriginal)
  allocate(idumoriginal(nseed))  
end if


!print *,"1st commandline input: name of the input file, 0 if no input file"
!print *,"2nd commandline input: 1 automatic, 2 automatic and saving, "
!print *,"3 the same but showing the images, 4 to just make a gif of the file, 5 if called from reva.f90" !>>> Is 22-1-14
!print *,"3th commandline input: iteration step (only for automatics)"
!print *,"4th commandline input: number of snapshots (only for automatics)"

call cpu_time(start_time) !;print*,"start",start_time !!>>HC 17-11-2020 This is necessary for filtering
whichend=0 !!>>HC 26-2-2021 means it has not been filtered (yet)

call getarg(2,cu)
if (len_trim(cu)==1) then
  read (cu,*) aut
end if
!if (aut==1) then ; print *,"automatic mode" ; end if
call getarg(1,carg)
precarg=carg
if (len(precarg)<2) then
  precarg(198:200)="kkk"
end if
do i=1,len(precarg)
  if (precarg(i:i)==" ") precarg(i:i)="_"
end do
ordre="rm -f *kk " !2> kok" 4-3-2020
!call system(ordre) !!>> HC 18-11-2020 This gives problems with evolution
call getarg(0,cazero)
caordre="cksum "//cazero//" > "//precarg//"kk"
call system(caordre) !assumes this file is called elli.e
open(1,file=precarg//"kk")
read(1,*) precaa !i
precaa=adjustl(precaa)
close(1)
ordre="date > "//precarg//"date"
call system(ordre)
open(1,file=precarg//"date")
read(1,*) cad,cae,caf,cag  
close(1)
ordre="rm -f "//precarg//"kk"
call system(ordre)
ordre="rm -f "//precarg//"date"
call system(ordre)
!print *,precaa,cab,cacc,cad,"cad",cae,"cae",caf,"caf",cag,"cag"
iu=1
id=2

ccag=cag
do i=1,len(cag)
  if (ccag(i:i)==":") then ; do j=i,len(ccag)-1 ; ccag(j:j)=ccag(j+1:j+1) ; end do ; end if
end do
ccag(7:8)=" "

winame=precaa//cazero(3:10)//cae//caf//cag//char(0)

caa=precaa(:10)//"_"//ccag(:6)
do i=1,len(caa) ; if (caa(i:i)==" ") caa(i:i)="_" ; end do 

!print *,ccag,"ccag elli"

call getarg(2,cu)    !>>> Is 22-1-14
!print *,2,cu,"that" !>>> Is 22-1-14
if (len_trim(cu)<1) then   !>>> Is 22-1-14
  aut=0               !>>> Is 22-1-14
else                  !>>> Is 22-1-14
  read (cu,*) aut     !>>> Is 22-1-14
end if                !>>> Is 22-1-14

if (aut==5) eva=1

if (aut/=5) then
 ordre="ls config_file.txt ; echo $? > notafile"  !>>Miquel1-10-14
 call system(ordre)
 open(1,file="notafile")
 read(1,*) ii
print*,"ii",ii
 if(ii==0)then           !>>Miquel1-10-14
   call read_config_file !>>Miquel8-9-14
 else
   call no_config_file
 end if
 close(1);
 ordre="rm notafile"
 call system(ordre)
else
   call no_config_file
end if

call initials
winame=adjustl(winame)
!print *,winame,"winame"

print *,nd,"nd"
if (aut/=1.and.aut/=5) then !>>> Is 22-1-14


  !Initializations
!  allocate(point)
  pgm="thePgm"
  point=c_null_ptr
!  call glutInit(1,iv)
  call glutinit(1,point) 
  call glutInitDisplayMode(mode) 
  call glutInitWindowPosition(1000,1000)        !>>>>>>Miquel21-2-14
  call glutInitWindowSize(windW,windH)        !>>>>>>Miquel26-11-13

  winid = glutCreateWindow(winame)
  !winid = glutCreateWindow("w")

  !initialize view_modifier, receiving the id for it's submenu
  submenuid = view_modifier_init()

  ! create the menu
  call make_menu(submenuid) 

  ! Set the display callback
  call glutDisplayFunc(display) 

  ! Set the window reshape callback

  call glutReshapeFunc(Reshape) 


  ! set the lighting conditions
  call glClearColor(0.0_glclampf, 0.0_glclampf, 0.0_glclampf, 1.0_glclampf)
  !call glClearColor(1.0_glclampf, 1.0_glclampf, 1.0_glclampf, 1.0_glclampf)
  call glLightfv(gl_light0, gl_diffuse, (/1.,1.,1.,1./))
  call glLightfv(gl_light0, gl_position, (/1.5,-.5,2.0,0.0/))
  call glEnable(gl_lighting)
  call glEnable(gl_light0)
  call glLightModelfv(gl_light_model_ambient, (/.5,.5,.5,1./))
  call glDepthFunc(gl_lequal)
  call glEnable(gl_depth_test)
  call inivisualitzacio 

end if

fgif=0 ; fgifmovie=0

call getarg(2,cu)

if (len_trim(cu)<1) then
  aut=0
else
  read (cu,*) aut
end if

!<<< Is 13-3-15 nodeo=node   !Is 25-5-13 this is just to save the initial conditions

!nodeoo=nodeo  !>>> Is 14-3-15 this is non-sense but without it the whole thing crushes, I think it has to do with a point pointer !!>> HC 4-12-2020

select case(aut)
case(1)
  call auto
case(5) !evolution
  eva=1 !>>> Is 22-1-14
  call auto
case default

  print*, "      " !!>> HC 16-11-2020 new versions of gfortran require some output here
  print*, "      " !!>> HC 16-11-2020
  print*,"EmbryoMaker software (General Node Model)"
  print*,"Computational model to simulate morphogenetic processes in living organs and tissues."
  print*,"Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad"
  print*,"This program comes with ABSOLUTELY NO WARRANTY"
  print*,"This is free software, and you are welcome to redistribute it under certain conditions;"
  print*,"read the LICENSE.txt file for details."
  print*, "      " !!>> HC 16-11-2020
  print*, "      " !!>> HC 16-11-2020

  call glutMainLoop
end select

call exit(1)

end program veu


