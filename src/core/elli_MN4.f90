
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
!!>> HC 30-11-2020 This is an alternative version of elli written by Hugo Cano to run in the supercomputer MareNostrum4
!!>> HC 30-11-2020 it does not include viewer and it does not write confirmation files

!use opengl             !!>> HC 30-11-2020 Viewer disabled
!use opengl_gl          !!>> HC 30-11-2020
!use opengl_glu         !!>> HC 30-11-2020
!use opengl_glut        !!>> HC 30-11-2020
!use opengl_kinds       !!>> HC 30-11-2020
!use view_modifier      !!>> HC 30-11-2020
!use function_plotter   !!>> HC 30-11-2020
use inicial
use automaticon

implicit none
character*400 caordre,ordre !!>>HC 18-11-2020 these should be a bit longer for evolution
character*200 precarg       !!>>HC 18-11-2020 (evolve calls with complete paths)
character*80 cu
character*10 pgm,xcad

integer :: winid, submenuid!, mode=ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)) !!>> HC 30-11-2020 no drawing

integer*4 iu,id

integer*8 iv
!type(c_ptr) point !!>> HC 30-11-2020

!this is necessary to accomadate to different versions of linux
call random_seed(size = nseed)
!print*,nseed  !!>> HC 30-11-2020 reduce the number of prints
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
!if (len(precarg)<2) then                                                       !!>> HC 30-11-2020 These checking files give problems in evolution
!  precarg(198:200)="kkk"                                                       !!>> HC 30-11-2020 when we are running multiple individuals in the same
!end if                                                                         !!>> HC 30-11-2020 directory and at the same time
!do i=1,len(precarg)                                                            !!>> HC 30-11-2020
!  if (precarg(i:i)==" ") precarg(i:i)="_"                                      !!>> HC 30-11-2020
!end do                                                                         !!>> HC 30-11-2020
!ordre="rm -f *kk " !2> kok" 4-3-2020                                           !!>> HC 30-11-2020
!call system(ordre) !!>> HC 18-11-2020 This gives problems with evolution       !!>> HC 30-11-2020
!call getarg(0,cazero)                                                          !!>> HC 30-11-2020
!caordre="cksum "//cazero//" > "//precarg//"kk"                                 !!>> HC 30-11-2020
!call system(caordre) !assumes this file is called elli.e                       !!>> HC 30-11-2020
!open(1,file=precarg//"kk")                                                     !!>> HC 30-11-2020
!read(1,*) precaa !i                                                            !!>> HC 30-11-2020
!precaa=adjustl(precaa)                                                         !!>> HC 30-11-2020
!close(1)                                                                       !!>> HC 30-11-2020
!ordre="date > "//precarg//"date"                                               !!>> HC 30-11-2020
!call system(ordre)                                                             !!>> HC 30-11-2020
!open(1,file=precarg//"date")                                                   !!>> HC 30-11-2020
!read(1,*) cad,cae,caf,cag                                                      !!>> HC 30-11-2020
!close(1)                                                                       !!>> HC 30-11-2020
!ordre="rm -f "//precarg//"kk"                                                  !!>> HC 30-11-2020
!call system(ordre)                                                             !!>> HC 30-11-2020
!ordre="rm -f "//precarg//"date"                                                !!>> HC 30-11-2020
!call system(ordre)                                                             !!>> HC 30-11-2020
!print *,precaa,cab,cacc,cad,"cad",cae,"cae",caf,"caf",cag,"cag"                !!>> HC 30-11-2020
iu=1
id=2

!ccag=cag                                                                                         !!>> HC 30-11-2020
!do i=1,len(cag)                                                                                  !!>> HC 30-11-2020
!  if (ccag(i:i)==":") then ; do j=i,len(ccag)-1 ; ccag(j:j)=ccag(j+1:j+1) ; end do ; end if      !!>> HC 30-11-2020
!end do                                                                                           !!>> HC 30-11-2020
!ccag(7:8)=" "                                                                                    !!>> HC 30-11-2020

!winame=precaa//cazero(3:10)//cae//caf//cag//char(0)                         !!>> HC 30-11-2020

!caa=precaa(:10)//"_"//ccag(:6)                                              !!>> HC 30-11-2020
!do i=1,len(caa) ; if (caa(i:i)==" ") caa(i:i)="_" ; end do                  !!>> HC 30-11-2020

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
! ordre="ls config_file.txt ; echo $? > notafile"  !>>Miquel1-10-14        !!>> HC 30-11-2020 this version does not have drawer so we do not need
! call system(ordre)                                                       !!>> HC 30-11-2020 config_file
! open(1,file="notafile")                                                  !!>> HC 30-11-2020
! read(1,*) ii                                                             !!>> HC 30-11-2020
!print*,"ii",ii                                                            !!>> HC 30-11-2020
! if(ii==0)then           !>>Miquel1-10-14                                 !!>> HC 30-11-2020
!   call read_config_file !>>Miquel8-9-14                                  !!>> HC 30-11-2020
! else                                                                     !!>> HC 30-11-2020
   call no_config_file
! end if                                                                   !!>> HC 30-11-2020
! close(1);                                                                !!>> HC 30-11-2020
! ordre="rm notafile"                                                      !!>> HC 30-11-2020
! call system(ordre)                                                       !!>> HC 30-11-2020
else
   call no_config_file
end if

call initials
!winame=adjustl(winame)    !!>> HC 30-11-2020
!print *,winame,"winame"

!print *,nd,"nd"  !!>> HC 30-11-2020 NO more prints
!if (aut/=1.and.aut/=5) then !>>> Is 22-1-14                                  !!>> HC 30-11-2020 NO aut=0 and no drawing


  !Initializations                                                            !!>> HC 30-11-2020
!  allocate(point)                                                            !!>> HC 30-11-2020
!  pgm="thePgm"                                                               !!>> HC 30-11-2020
!  point=c_null_ptr                                                           !!>> HC 30-11-2020
!  call glutInit(1,iv)                                                        !!>> HC 30-11-2020
!  call glutinit(1,point)                                                     !!>> HC 30-11-2020
!  call glutInitDisplayMode(mode)                                             !!>> HC 30-11-2020
!  call glutInitWindowPosition(1000,1000)        !>>>>>>Miquel21-2-14         !!>> HC 30-11-2020
!  call glutInitWindowSize(windW,windH)        !>>>>>>Miquel26-11-13          !!>> HC 30-11-2020

!  winid = glutCreateWindow(winame)                                           !!>> HC 30-11-2020
  !winid = glutCreateWindow("w")                                              !!>> HC 30-11-2020

  !initialize view_modifier, receiving the id for it's submenu                !!>> HC 30-11-2020
!  submenuid = view_modifier_init()                                           !!>> HC 30-11-2020

  ! create the menu                                                           !!>> HC 30-11-2020
!  call make_menu(submenuid)                                                  !!>> HC 30-11-2020

  ! Set the display callback                                                  !!>> HC 30-11-2020
!  call glutDisplayFunc(display)                                              !!>> HC 30-11-2020

  ! Set the window reshape callback                                           !!>> HC 30-11-2020

!  call glutReshapeFunc(Reshape)                                              !!>> HC 30-11-2020


  ! set the lighting conditions                                               !!>> HC 30-11-2020
!  call glClearColor(0.0_glclampf, 0.0_glclampf, 0.0_glclampf, 1.0_glclampf)  !!>> HC 30-11-2020
  !call glClearColor(1.0_glclampf, 1.0_glclampf, 1.0_glclampf, 1.0_glclampf)  !!>> HC 30-11-2020
!  call glLightfv(gl_light0, gl_diffuse, (/1.,1.,1.,1./))                     !!>> HC 30-11-2020
!  call glLightfv(gl_light0, gl_position, (/1.5,-.5,2.0,0.0/))                !!>> HC 30-11-2020
!  call glEnable(gl_lighting)                                                 !!>> HC 30-11-2020
!  call glEnable(gl_light0)                                                   !!>> HC 30-11-2020
!  call glLightModelfv(gl_light_model_ambient, (/.5,.5,.5,1./))               !!>> HC 30-11-2020
!  call glDepthFunc(gl_lequal)                                                !!>> HC 30-11-2020
!  call glEnable(gl_depth_test)                                               !!>> HC 30-11-2020
!  call inivisualitzacio                                                      !!>> HC 30-11-2020

!end if                                                                       !!>> HC 30-11-2020

!fgif=0 ; fgifmovie=0                                                         !!>> HC 30-11-2020

call getarg(2,cu)

if (len_trim(cu)<1) then
  aut=0
else
  read (cu,*) aut
end if

!<<< Is 13-3-15 nodeo=node   !Is 25-5-13 this is just to save the initial conditions

!nodeoo=nodeo  !>>> Is 14-3-15 this is non-sense but without it the whole thing crushes, I think it has to do with a point pointer !!>> HC 4-12-3030

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

!  call glutMainLoop   !!>> HC 30-11-2020 No drawing
end select

call exit(1)

end program veu


