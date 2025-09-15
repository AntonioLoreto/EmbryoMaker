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




!************************************************************************************************************************************************
!              MODUL VIEW_MODIFIER
!*************************************************************************************************************************************************


module view_modifier
use model
!use opengl
use opengl_gl
use opengl_glu
use opengl_glut
use opengl_kinds
use editor
use io
implicit none
private
      public :: view_modifier_init,reset_view,arrow_key_func,mouseswitch,pressed_left,switch,& !Tommi
                reset_to_init,colorswitch,view_from_front,left_button_func,CURPAN,CURPANZ,middle_button_func, & !Tommi 13.8.2013
                ROTATE, ZOOM, SELNO, SELCEL, middle_press, COLORMAX , COLORMIN ,& !>>Miquel27-10-14
                ARROWMAX , ARROWMIN , SPHEREMAX, SPHEREMIN, ARROWSCALE, SPHERESCALE,&  !>>Miquel27-10-14
                NOPAN, NOPANZ, CEPAN, CEPANZ !>>Miquel31-10-14

private :: PAN, SCALEX, SCALEY, SCALEZ, RESET, ABOVE, QUIT, PI, &
           init_lookat, init_lookfrom, &
           init_xscale_factor, init_yscale_factor, init_zscale_factor, & 
           anglep, shift, xscale_factor, yscale_factor, zscale_factor, &
           moving_left, moving_middle, begin_left, begin_middle, & !Tommi 10.9.2013
           cart2sphere, sphere2cart, cart3D_plus_cart3D, cart3D_minus_cart3D, &
           mouse, motion, arrows, &
           menu_handler, set_left_button, set_middle_button, set_arrow_keys
integer(kind=glint), parameter :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6, PAS=17, CURPAN=18, CURPANZ=19, SELNO=20, SELCEL=21, &
                      COLORMAX=22 , COLORMIN=23,ARROWMAX=24 , ARROWMIN=25 , SPHEREMAX=26, &
                      SPHEREMIN=27, ARROWSCALE=28, SPHERESCALE=29, & !>>Miquel27-10-14
                      NOPAN=40, NOPANZ=41, CEPAN=42, CEPANZ=43 !>>Miquel31-10-14
integer(kind=glint), parameter :: RESET = 10, ABOVE = 11, FRONT = 12, pers = 14, QUIT = 13
!real(kind=gldouble), parameter :: PI = 3.141592653589793_gldouble
real(kind=gldouble), public :: esca
integer(glint), public :: iqpas


type, private :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type, private :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type, private :: sphere3D ! 3D spherical coordinates
   real(kind=gldouble) :: theta, phi, rho
end type sphere3D

type(cart2D), save :: anglep
type(cart3D), save :: shift
real(kind=gldouble), save :: xscale_factor, yscale_factor, zscale_factor
logical, save :: moving_left, moving_middle, pressed_left,mouseswitch, mouseswitch2, & !Tommi
		 cellselectswitch, switch, altreswitch, colorswitch,cellswitch !Tommi 7.8.2013!

type(cart2D), save :: begin_left, begin_middle

interface operator(+)
   module procedure cart3D_plus_cart3D
end interface
interface operator(-)
   module procedure cart3D_minus_cart3D
end interface

! ------- Initial configuration -------

! Set the initial operation performed by each button and the arrow keys.
! The operations are ZOOM, PAN, ROTATE, SCALEX, SCALEY, and SCALEZ

integer, save ::   left_button_func = ROTATE, &
                 middle_button_func = ZOOM, &
                    arrow_key_func = PAS, &
                    middle_press= 0  !>>Miquel29-10-14


! Set the initial view as the point you are looking at, the point you are
! looking from, and the scale factors


type(cart3D), parameter :: &
   init_lookat = cart3D(0.0_gldouble, 0.0_gldouble, 0.0_gldouble), &
   init_lookfrom = cart3D(0.0_gldouble, 0.0_gldouble, 5.0_gldouble)	!display: distance of the camera from origin


!type(cart3D) :: init_lookat,init_lookfrom	!display: distance of the camera from origin


real(kind=gldouble), parameter :: &
   init_xscale_factor = 1.0_gldouble, &
   init_yscale_factor = 1.0_gldouble, &
   init_zscale_factor = 1.0_gldouble

!variables
real*8 ,public                ::mix,miy,miz,mx,my,mz,start,startxmax,startxmin,startymax,startymin,&
                                copy,tempo,factorx,ymaxtemp,zmaxtemp,converymax,converzmax, &!Tommi 16.9.2013
                                set_mix,set_miy,set_miz,set_mx,set_my,set_mz !>>Miquel30-9-14
real*8 ,public                ::maan,mapo,mesc,xcoord,ycoord,zcoord, xcoordmous, ycoordmous,&
                                yconver, zconver, ymax, zmax,converyconver,converzconver,change !Tommi 16.9.2013
integer,public                ::node_to_move,buttstate,butt,choice,plane,pasteselection!Tommi 27.9.2013
integer,public                ::make_flag_twe_five

integer(GLint),public :: windW = 800, windH = 800                    !>>>>>>Miquel26-11-13
real(GLdouble),public :: coolperspec=10d0    !perspective parameter     !>>>>>>Miquel26-11-13

real*8, public :: cursx,cursy,cursz

integer,public                ::ki,nki       !>>Miquel29-10-14
integer,public,allocatable    ::oopp(:)      !warning: arbitrary size

integer, public::dyn, fix_run !>>Miquel2-10-14
integer,public :: custom_colorselection,custom_arrowselection,custom_sphereselection
real*8,public :: custom_minval,custom_maxval,custom_aminval,custom_amaxval,custom_sminval,custom_smaxval
real*8,public ::amax_scale,smax_scale  !>>>Miquel30-10-14

! -------- end of Initial configuration ------

contains
!*******************************************SUBROUTINE********************************************************
!          -------------
subroutine reset_to_init
!          -------------

! This resets the view to the initial configuration

type(sphere3D) :: slookfrom
type(cart3D):: lookat, lookfrom

   !lookat = cart3D(conf_lookatx, conf_lookaty, conf_lookatz)
   !lookfrom = cart3D(conf_lookfromx, conf_lookfromy, conf_lookfromz)	!display: distance of the camera from origin

!slookfrom = cart2sphere(lookfrom-lookat)
!slookfrom = cart2sphere(init_lookfrom-init_lookat)
anglep%x = 0.0*slookfrom%theta/PI + conf_anglex !- 45.0_gldouble			!display, angle of camera respect the origin
anglep%y = 0.0*slookfrom%phi/PI + conf_angley !- 60.0_gldouble
shift%x = conf_shiftx  !>>Miquel5-11-14
shift%y = conf_shifty
shift%z = -conf_shiftz !-slookfrom%rho
xscale_factor = init_xscale_factor
yscale_factor = init_yscale_factor
zscale_factor = init_zscale_factor

cursx=maxval(node(:)%x)+0.1 ; cursy=maxval(node(:)%y)+0.1 ; cursz=minval(node(:)%z)+0.1

call glutPostRedisplay

node_to_move=0
make_flag_twe_five=0 

return
end subroutine reset_to_init
!*******************************************SUBROUTINE********************************************************
!          ---------------
subroutine view_from_above
!          ---------------

! This sets the view to be from straight above

type(sphere3D) :: slookfrom

slookfrom = cart2sphere(cart3D(0.0,0.0,1.0))
anglep%x = -200.0_gldouble*slookfrom%theta/PI
anglep%y = -200.0_gldouble*slookfrom%phi/PI

call glutPostRedisplay

return
end subroutine view_from_above
!*******************************************SUBROUTINE******************************************************** >>>>>!Tommi 17.9.2013
!          ---------------
subroutine view_from_front
!          ---------------

! This sets the view to be from straight ahead for the editor from an arbitrary distance. And in different planes depending on the selection.
! Also gives values to some variables for different planes for the editor. 

type(sphere3D) :: slookfrom

select case (plane)

case(1) !minimal y plane
 slookfrom = cart2sphere(cart3D(1.0,0.0,0.0))
 anglep%x = -180_gldouble*slookfrom%theta/PI
 anglep%y = -180.0_gldouble*slookfrom%phi/PI
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0094522189
 zconver=0.0092592593
 ymax=725 !voit joutua siirtämään tommi
 zmax=402
 make_flag_twe_five=0
case(2) !maximal x plane
 slookfrom = cart2sphere(cart3D(0.0,1.0,0.0))
 anglep%x = -180.0_gldouble*slookfrom%theta/PI
 anglep%y = -180.0_gldouble*slookfrom%phi/PI
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0094087937
 zconver=0.0093457944
 ymax=726 ! These are "handpicked", so might be a might a pixel off. But quite accurate.
 zmax=401
 !startmx=mx
 make_flag_twe_five=0
case(3) !minimal x plane
 slookfrom = cart2sphere(cart3D(0.0,1.0,0.0))
 anglep%x = 180.0_gldouble*slookfrom%theta/PI
 anglep%y = -180.0_gldouble*slookfrom%phi/PI
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0094087937
 zconver=0.0093457944
 ymax=726 
 zmax=401
 !startmix=mix
 make_flag_twe_five=0
case(4) !maximal y plane
 slookfrom = cart2sphere(cart3D(1.0,0.0,0.0))
 anglep%x = -180.0_gldouble
 anglep%y = -90.0_gldouble
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0094522189
 zconver=0.0092592593
 ymax=725
 zmax=402
 make_flag_twe_five=0
case(5)!maximal z plane
 slookfrom = cart2sphere(cart3D(0.0,0.0,1.0))
 anglep%x = -200.0_gldouble*slookfrom%theta/PI
anglep%y = -200.0_gldouble*slookfrom%phi/PI
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0098936268
 zconver=0.0099009901
 ymax=711
 zmax=402
 make_flag_twe_five=0
case(6)!minimal z plane
 slookfrom = cart2sphere(cart3D(0.0,0.0,1.0))
 anglep%x = 0
 anglep%y = -180.0_gldouble
 shift%x = 0.0_gldouble
 shift%y = 0.0_gldouble
 shift%z=-4.617117932711601
 startxmax=3.071971135656903
 startxmin=-3.071971135656904
 startymax=1.0
 startymin=0.0
 yconver=0.0100885752
 zconver=0.0102040816
 ymax=704
 zmax=501
 make_flag_twe_five=0
end select
xscale_factor = init_xscale_factor
yscale_factor = init_yscale_factor
zscale_factor = init_zscale_factor
colorswitch = .true. !turn on the selection plane coloring

call glutPostRedisplay

return
end subroutine view_from_front
!*******************************************SUBROUTINE********************************************************<<<<<<<<Tommi

!          ----------
subroutine reset_view
!          ----------
real(8),save :: x,y
! This routine resets the view to the current orientation and scale

call glMatrixMode(GL_MODELVIEW)
call glPopMatrix
call glPushMatrix
call glTranslated(shift%x, shift%y, shift%z)
call glRotated(anglep%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
 x=dcos(PI*anglep%x/180.0_gldouble)
 y=-dsin(PI*anglep%x/180.0_gldouble)
 call glRotated(anglep%y, x, y, 0.0_gldouble)
 call glTranslated(-init_lookat%x, -init_lookat%y, -init_lookat%z)
call glScaled(xscale_factor,yscale_factor,zscale_factor)

return
end subroutine reset_view
!*******************************************SUBROUTINE********************************************************
!          -----
subroutine mouse(button, state, x, y) bind(c) !>>>>>>Tommi 17.9.2013, heavily modified
!          -----
integer(kind=glint), value :: button, state, x, y

! This gets called when a mouse button changes
 
buttstate=state
butt=button
xcoordmous=x
ycoordmous=y
!print*,x,"xcoord",y,"ycoord"
pressed_left=.false.
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving_left = .true.
    begin_left = cart2D(x,y)
  endif
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving_left = .false.
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_DOWN) then
    moving_middle = .true.
    begin_middle = cart2D(x,y)
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_UP) then
    moving_middle = .false.

    if(middle_press==SELNO)then
      call selectnode(cursx,cursy,cursz)
      print*,"position",node(nodeindex)%x,node(nodeindex)%y,node(nodeindex)%z
      print*,"tipus",node(nodeindex)%tipus
      if(node(nodeindex)%tipus<3) print*,"down node",node(nodeindex)%altre
      !print*,"energy",node(nodeindex)%e
      print*,"cell",node(nodeindex)%icel
      print *,node(nodeindex)%altre,"altre"
      !node(nodeindex)%e=guardaene
      print*,"gene expression"
      print*,gex(nodeindex,:ng)
      nki=nki+1
      oopp(nki)=nodeindex
    elseif(middle_press==SELCEL)then
      call selectcell(cursx,cursy,cursz)
      if (cellid>0) then
        !do kk=1,cels(cellid)%nunodes
          nki=nki+1
          !ki=cels(cellid)%node(kk)
          !print *,kk,"th node in cell",cellid,"is",ki,"its tipus is",node(ki)%tipus
          oopp(nki)=cellid
        !end do
      end if
    end if

  endif

end subroutine mouse !<<<<<<<<<Tommi

!*******************************************SUBROUTINE********************************************************!

!          ------
subroutine motion(x, y) bind(c)
!          ------
integer(kind=glint), value :: x, y

! This gets called when the mouse moves

integer :: button_function
type(cart2D) :: begin
real(kind=gldouble) :: factor
real(kind=gldouble) :: norientx,norienty !>>Miquel28-10-14
type(cart3D) :: orient !>>Miquel28-10-14
type(sphere3D) :: sorient !>>Miquel28-10-14


! Determine and apply the button function


if (moving_left) then
   button_function = left_button_func
   begin = begin_left
else if(moving_middle) then
   button_function = middle_button_func
   begin = begin_middle
end if


select case(button_function)
case (ZOOM)
   if (y < begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(begin%y-y))
   else if (y > begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(y-begin%y)
   else
      factor = 1.0_gldouble
   end if
   copy=shift%z !<<<<<<Tommi 11.9.2013 make the selection/adding window to behave correctly when zooming
   shift%z = factor*shift%z
   tempo=(copy-shift%z)/copy
   start=startxmax
   startxmax = startxmax-(startxmax*tempo)
   change=(start-startxmax)/0.3326719033
   startxmin = startxmin-(startxmin*tempo)
   startymax = startymax-(startymax*tempo)
   startymin = startymin-(startymin*tempo)
   if (shift%z<-4.6171179327) then
     converyconver=change*(-0.0010917615)
     converzconver=change*(-0.0011125007)
     converymax=change*1.6
    else
   converyconver=change*(-0.0010931038)
   converzconver=change*(-0.0010686245)
   if (shift%z > -2.6171178818) then
   converymax=change*1
   else
   converymax=change*7
   endif
   endif
   ymax=ymax+converymax
   yconver=yconver+converyconver
   zconver=zconver+converzconver
   !print*,factor,"factor",yconver,"yconver"
   !print*,shift%z,"shiftz"!, yconver,"yconver"
   !print*,shift%z,"shiftz",startxmax,"startxmax",startxmin,"startxmin",startymax,"startymax",startymin,"startymin"
case (PAN) !<<<<<<Tommi, make the selection/adding window to stay in place while moving and also update the x and y coordinates according to each plane
   shift%x = shift%x + .01*(x - begin%x)
   shift%y = shift%y - .01*(y - begin%y)
   select case(plane)
   case(1)
   startxmax = startxmax - .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin - .1*(x - begin%x)
   startymax = startymax + .1*(y - begin%y)
   startymin = startymin + .1*(y - begin%y) 
   case(2)
   startxmax = startxmax - .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin - .1*(x - begin%x)
   startymax = startymax + .1*(y - begin%y)
   startymin = startymin + .1*(y - begin%y) 
   case(3)
   startxmax = startxmax + .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin + .1*(x - begin%x)
   startymax = startymax + .1*(y - begin%y)
   startymin = startymin + .1*(y - begin%y) 
   case(4)
   startxmax = startxmax + .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin + .1*(x - begin%x)
   startymax = startymax + .1*(y - begin%y)
   startymin = startymin + .1*(y - begin%y)
   case(5)
   startxmax = startxmax - .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin - .1*(x - begin%x)
   startymax = startymax + .1*(y - begin%y)
   startymin = startymin + .1*(y - begin%y) 
   case(6)
   startxmax = startxmax - .1*(x - begin%x) !We keep the selection window still, while the embryo moves!
   startxmin = startxmin - .1*(x - begin%x)
   startymax = startymax - .1*(y - begin%y)
   startymin = startymin - .1*(y - begin%y)
   end select !>>>>>Tommi
case (CURPAN)

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)

   if(orient%z>=0)then
     b=.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a + .01*(x - begin%x)*norienty*a
   else
     b=-.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a - .01*(x - begin%x)*norienty*a
   end if

   cursx = cursx + b
   cursy = cursy + c


!print*,"increments x",(x-begin%x),"y",(y-begin%y)
!print*,"orient",orient%x,orient%y
!print*,"norient",norientx,norienty
!print*,"b",b,"c",c
!print*,"cursx",cursx,"cursy",cursy
!print*,

case (CURPANZ)

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)
   if(orient%z>=0)then
     cursz = cursz - .01*(y - begin%y)
   else
     cursz = cursz + .01*(y - begin%y)
   endif

case (NOPAN)

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)

   if(orient%z>=0)then
     b=.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a + .01*(x - begin%x)*norienty*a
   else
     b=-.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a - .01*(x - begin%x)*norienty*a
   end if

   node(nodeindex)%x = node(nodeindex)%x + b
   node(nodeindex)%y = node(nodeindex)%y + c


!print*,"increments x",(x-begin%x),"y",(y-begin%y)
!print*,"orient",orient%x,orient%y
!print*,"norient",norientx,norienty
!print*,"b",b,"c",c
!print*,"cursx",cursx,"cursy",cursy
!print*,
case (NOPANZ)
   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)
   if(orient%z>=0)then
    node(nodeindex)%z = node(nodeindex)%z - .01*(y - begin%y)
   else
    node(nodeindex)%z = node(nodeindex)%z + .01*(y - begin%y)
   end if

case (CEPAN)

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)

   if(orient%z>=0)then
     b=.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a + .01*(x - begin%x)*norienty*a
   else
     b=-.01*(x - begin%x)*norientx*a - .01*(y - begin%y)*orient%x*a
     c=-.01*(y - begin%y)*orient%y*a - .01*(x - begin%x)*norienty*a
   end if

   do i=1,cels(cellid)%nunodes
     j=cels(cellid)%node(i)
     node(j)%x=node(j)%x+b
     node(j)%y=node(j)%y+c
   end do
   cels(cellid)%cex=cels(cellid)%cex+b
   cels(cellid)%cey=cels(cellid)%cey+c

!print*,"increments x",(x-begin%x),"y",(y-begin%y)
!print*,"orient",orient%x,orient%y
!print*,"norient",norientx,norienty
!print*,"b",b,"c",c
!print*,"cursx",cursx,"cursy",cursy
!print*,
case (CEPANZ) 

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)   
   orient%x=-orient%x ; orient%y=-orient%y
   
   norientx=orient%y
   norienty=-orient%x
   a=1/sqrt(norientx**2+norienty**2)
   if(orient%z>=0)then
     c= - .01*(y - begin%y)
   else
     c= + .01*(y - begin%y)
   end if

   do i=1,cels(cellid)%nunodes
     j=cels(cellid)%node(i)
     node(j)%z=node(j)%z+c
   end do
   cels(cellid)%cez=cels(cellid)%cez+c


case (ROTATE)
   anglep%x = anglep%x + (x - begin%x)
   anglep%y = anglep%y + (y - begin%y)

   sorient%theta=anglep%y*(2*PI)/360
   sorient%phi=-(anglep%x-90)*(2*PI)/360
   sorient%rho=-shift%z
   orient=sphere2cart(sorient)


!print*,"anglex",-anglep%x-90
!print*,"angley",anglep%y
!print*,"sorient",sorient%theta,sorient%phi,sorient%rho
!print*,"orient",orient%x,orient%y,orient%z
!print*,

case (COLORMIN)
  if(custom_minval<epsilod) custom_minval=1d-4
  custom_minval=custom_minval- 0.1*custom_minval*(y - begin%y)
  print*,"minimum color value",custom_minval
  print*, "     " !!>> HC 16-11-2020

case (COLORMAX)
  custom_maxval=custom_maxval- 0.1*custom_maxval*(y - begin%y)
  print*,"maximum color value",custom_maxval
  print*, "     " !!>> HC 16-11-2020
case (ARROWSCALE)
  amax_scale=amax_scale - 0.1*amax_scale*(y - begin%y)
  print*,"maximum arrow scale",amax_scale
  print*, "     " !!>> HC 16-11-2020
case (SPHERESCALE)
  smax_scale=smax_scale - 0.1*smax_scale*(y - begin%y)
  print*,"maximum sphere scale",smax_scale
  print*, "     " !!>> HC 16-11-2020
case (ARROWMIN)
  if(custom_aminval<epsilod) custom_aminval=1d-4
  custom_aminval=custom_aminval- 0.1*custom_aminval*(y - begin%y)
  print*,"minimum arrow value",custom_aminval
  print*, "     " !!>> HC 16-11-2020
case (ARROWMAX)
  custom_amaxval=custom_amaxval- 0.1*custom_amaxval*(y - begin%y)
  print*,"maximum arrow value",custom_amaxval
  print*, "     " !!>> HC 16-11-2020
case (SPHEREMIN)
  if(custom_sminval<epsilod) custom_sminval=1d-4
  custom_sminval=custom_sminval- 0.1*custom_sminval*(y - begin%y)
  print*,"minimum sphere value",custom_sminval
  print*, "     " !!>> HC 16-11-2020
case (SPHEREMAX)
  custom_smaxval=custom_smaxval- 0.1*custom_smaxval*(y - begin%y)
  print*,"maximum sphere value",custom_smaxval
  print*, "     " !!>> HC 16-11-2020
case (SCALEX)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   xscale_factor = xscale_factor * factor
case (SCALEY)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   yscale_factor = yscale_factor * factor
case (SCALEZ)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   zscale_factor = zscale_factor * factor
case (8)
  ! copy=mx
   dyn=0
   make_flag_twe_five=2
    set_mix=set_mix+extre*mesc*(x-begin%x)
   !mx=9.1389584338
   !print*,mix,"mix"
   tempo=(extre*mesc*(x-begin%x))/1.0
   if (set_mix <-3.1389584338) then
   converymax=tempo*8.8333333333
   converyconver=tempo*(-0.0002190013)
   if (set_mix <-4.1389584541) then
   converzconver=tempo*(-0.0002258356)
   else
   converzconver=tempo*(-8.65351407407424E-005)
   endif
   else
   converymax=tempo*6.8
   converyconver=tempo*(-0.0002223311)
   converzconver=tempo*(-0.0002361043)
   endif
   ymax=ymax-converymax
   yconver=yconver-converyconver
   zconver=zconver-converzconver
   !print*,ymax,"ymax",yconver,"yconver"
   !yconver=yconver*tempo
   ! print *,"mx",mx, "tempo",tempo,"yconver",yconver
case (9)
  ! copy=mx 
   dyn=0
   make_flag_twe_five=1
   set_mx=set_mx+extre*mesc*(x-begin%x)
   !mx=9.1389584338
   !print*,mx,"mx"
   tempo=(extre*mesc*(x-begin%x))/1.0
   if (set_mx >3.1389584338) then
   converymax=tempo*8.8333333333
   converyconver=tempo*(-0.0002190013)
   if (set_mx > 4.1389584541) then
   converzconver=tempo*(-0.0002258356)
   else
   converzconver=tempo*(-8.65351407407424E-005)
   endif
   else
   converymax=tempo*6.8
   converyconver=tempo*(-0.0002223311)
   converzconver=tempo*(-0.0002361043)
   endif
   ymax=ymax+converymax
   yconver=yconver+converyconver
   zconver=zconver+converzconver
   !print*,ymax,"ymax",xcoordmous,"xcoordmous",yconver,"yconver",zconver,"zconver"
   !yconver=yconver*tempo
   ! print *,"mx",mx, "tempo",tempo,"yconver",yconver
case (10)
   dyn=0
   make_flag_twe_five=4
   set_miy=set_miy+extre*mesc*(y-begin%y)
   !miy=-8.071971135656904
   tempo=(extre*mesc*(y-begin%y))/1.0
   if (set_miy < -3.0719711357) then
   converymax=tempo*8.6
   converyconver=tempo*(-0.0002208942)
   converzconver=tempo*(-0.0002162663)
   else
   converymax=tempo*6.75
   converyconver=tempo*(-0.0002141023)
   converzconver=tempo*(-0.0002362056)
   endif
   ymax=ymax-converymax
   yconver=yconver-converyconver
   zconver=zconver-converzconver
   !print*,ymax,"ymax",xcoordmous,"xcoordmous",yconver,"yconver",zconver,"zconver"
   !print *,"miy",miy
case (11)
   dyn=0
   make_flag_twe_five=3
   set_my=set_my+extre*mesc*(y-begin%y)
   !my=-0.9280288643
   tempo=(extre*mesc*(y-begin%y))/1.0
   if (set_my >3.0719711357) then
   converymax=tempo*8.6
   converyconver=tempo*(-0.0002208942)
   converzconver=tempo*(-0.0002216539)
   else
   converymax=tempo*6.75
   converyconver=tempo*(-0.0002141023)
   converzconver=tempo*(-0.0002362056)
   endif
   ymax=ymax+converymax
   yconver=yconver+converyconver
   zconver=zconver+converzconver
   !print*,ymax,"ymax",yconver,"yconver",zconver,"zconver"
   !print *,"my",my
case (12)
   dyn=0
   make_flag_twe_five=6
   set_miz=set_miz+extre*mesc*(y-begin%y)
   tempo=(extre*mesc*(y-begin%y))/1.0
   if (set_miz < 0.0) then
   converymax=tempo*7.25
   converzmax=tempo*2.25
   converyconver=tempo*(-0.0002200074)
   converzconver=tempo*(-0.0002104377)
   else
   converymax=tempo*5.8
   converzmax=tempo*2.25
   converyconver=tempo*(-0.0002131272)
   converzconver=tempo*(-0.0002525253)
   endif
   ymax=ymax-converymax
   zmax=zmax-converzmax
   yconver=yconver-converyconver
   zconver=zconver-converzconver
   !print*,ymax,"ymax",yconver,"yconver",zconver,"zconver"
   !miz=-4
   !print *,"miz",miz
case (13)
   dyn=0
   make_flag_twe_five=5
   set_mz=set_mz+extre*mesc*(y-begin%y)
   !mz=3.0
   tempo=(extre*mesc*(y-begin%y))/1.0
   if (set_mz > 1.0) then
   converymax=tempo*8
   converyconver=tempo*(-0.0002265048)
   converzconver=tempo*(-0.0002063983)
   else
   converymax=tempo*6
   converyconver=tempo*(-0.0002123482)
   converzconver=tempo*(-0.0002131287)
   endif
   ymax=ymax+converymax
   yconver=yconver+converyconver
   zconver=zconver+converzconver
   !print*,ymax,"ymax",yconver,"yconver",zconver,"zconver"
   !mz=-4.0
   !print *,"mz",mz
case(14)
   maan=maan+mesc*(y-begin%y)
   if (maan<0) maan=0
   call glMatrixMode(GL_PROJECTION)
   call glloadidentity()
   call gluPerspective(10.0_gldouble, 1.0_gldouble, maan, mapo)
   call glMatrixMode(GL_MODELVIEW)
   print *,"maan",maan
case(15)
   mapo=mapo+mesc*(y-begin%y)
   if (maan<0) mapo=0
   call glMatrixMode(GL_PROJECTION)
   call glloadidentity()
   call gluPerspective(10.0_gldouble, 1.0_gldouble, maan, mapo)
   call glMatrixMode(GL_MODELVIEW)
   print *,"mapo",mapo
case(16)
   mesc=mesc+0.001d0*(y-begin%y)
   print *,"mesc",mesc
case(30) !set the main display box back to dynamic
    dyn=1
case(31) !set the main box as fixed even after running iterations
    fix_run=1
end select

! update private variables and redisplay

if (moving_left) then
   begin_left = cart2D(x,y)
else if(moving_middle) then
   begin_middle = cart2D(x,y)
endif

if (moving_left .or. moving_middle) then
   call glutPostRedisplay
endif

return
end subroutine motion
!*******************************************SUBROUTINE******************************************************** !>>>>>>Tommi
!subroutine motionadd(x, y) bind(c) !Tommi 7.8.2013
!          ------
!integer(kind=glint), intent(in), value :: x, y

!print*,x,"x",y,"y"
!end subroutine motionadd
!*******************************************SUBROUTINE********************************************************
!          ------
subroutine arrows(key, x, y) bind(c)
!          ------
integer(glint), value :: key, x, y

! This routine handles the arrow key operations

real(kind=gldouble) :: factor

select case(arrow_key_func)
case(ZOOM)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble + .02_gldouble
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case default
      factor = 1.0_gldouble
   end select
   shift%z = factor*shift%z
case(PAN)
   select case(key)
   case(GLUT_KEY_LEFT)
      shift%x = shift%x - .02
   case(GLUT_KEY_RIGHT)
      shift%x = shift%x + .02
   case(GLUT_KEY_DOWN)
      shift%y = shift%y - .02
   case(GLUT_KEY_UP)
      shift%y = shift%y + .02
   end select
case(PAS)
   factor=1.0_gldouble
   select case(key)
   case(GLUT_KEY_DOWN)
     if (left_button_func == ZOOM) then
       left_button_func = ROTATE
     elseif(left_button_func == ROTATE) then
      left_button_func = ZOOM
     end if
     if(left_button_func == CURPAN) then
       left_button_func = CURPANZ
     elseif(left_button_func == CURPANZ)then
       left_button_func = CURPAN
     end if

     if(left_button_func == NOPAN) then
       left_button_func = NOPANZ
     elseif(left_button_func == NOPANZ)then
       left_button_func = NOPAN
     end if
     if(left_button_func == CEPAN) then
       left_button_func = CEPANZ
     elseif(left_button_func == CEPANZ)then
       left_button_func = CEPAN
     end if

   case(GLUT_KEY_UP)
     iqpas=iqpas*10

   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if

     call iteracio(iqpas)
   case(GLUT_KEY_RIGHT)

   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if

     call iteracio(iqpas)
   case(GLUT_KEY_LEFT)
     call go_iteration_back(1)
   end select
case(ROTATE)
   select case(key)
   case(GLUT_KEY_LEFT)
      anglep%x = anglep%x - 1.0_gldouble
   case(GLUT_KEY_RIGHT)
      anglep%x = anglep%x + 1.0_gldouble
   case(GLUT_KEY_DOWN)
      anglep%y = anglep%y + 1.0_gldouble
   case(GLUT_KEY_UP)
      anglep%y = anglep%y - 1.0_gldouble
   end select
case(SCALEX)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   xscale_factor = xscale_factor * factor
case(SCALEY)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   yscale_factor = yscale_factor * factor
case(SCALEZ)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   zscale_factor = zscale_factor * factor
end select
   
call glutPostRedisplay

return
end subroutine arrows
!*******************************************SUBROUTINE********************************************************
!          ------------
subroutine menu_handler(value) bind(c)
!          ------------
integer(kind=glint), value :: value

! This routine handles the first level entries in the menu

select case(value)
case(ROTATE)
   middle_press=0
case(ZOOM)
   middle_press=0
case(PAN)
   middle_press=0
case(RESET)
   call reset_to_init
case(ABOVE)
   call view_from_above
case(FRONT)
   call view_from_front
case(pers)
   print*,"****PRECAUTION****"
   print*,"WHEN CHANGING THE PERSPECTIVE VALUE, THE DISPLAY MAY GET FUCKED UP,"
   print*,"DON'T PANIC, RESHAPE THE WINDOW AND YOU'LL SEE YOUR NODES AGAIN"
   print*,"actual perspective value = ",coolperspec
   print*,"enter new perspective value"
   read(*,*)coolperspec
   call Reshape(windW,windH)
   call reset_to_init
case(QUIT)
   call exit(1)
case(30) !set dynamic display box boundaries again
  dyn=1
case(31) !set fixed display box boundaries while running iterations
  if(fix_run==0)then
    fix_run=1
    dyn=0
    print*,"Fixed dysplay box boundaries while running iterations,"
    print*,"sections will be kept over simulating time"
  else
    fix_run=0
    !dyn=1
    print*,"Now the sections will be undone when next iteration is run"
  end if

end select

return
end subroutine menu_handler
!*******************************************SUBROUTINE********************************************************
!          ---------------
subroutine set_left_button(value) bind(c)
!          ---------------
integer(kind=glint), value :: value

! This routine sets the function of the left button as given by menu selection

left_button_func = value

return
end subroutine set_left_button
!*******************************************SUBROUTINE********************************************************
!          -----------------
subroutine set_middle_button(value) bind(c)
!          -----------------
integer(kind=glint), value :: value

! This routine sets the function of the middle button as given by menu selection

middle_button_func = value

return
end subroutine set_middle_button
!*******************************************SUBROUTINE********************************************************
!          --------------
subroutine set_arrow_keys(value) bind(c)
!          --------------
integer(kind=glint), value :: value

! This routine sets the function of the arrow keys as given by menu selection

arrow_key_func = value

return
end subroutine set_arrow_keys
!*******************************************SUBROUTINE********************************************************!
!        ------------------
function view_modifier_init() result(menuid)
!        ------------------
integer(kind=glint) :: menuid

! This initializes the view modifier variables and sets initial view.
! It should be called immediately after glutCreateWindow

integer(kind=glint) :: button_left, button_middle, arrow_keys

! set the callback functions

call glutMouseFunc(mouse)
call glutMotionFunc(motion)
call glutSpecialFunc(arrows)


! create the menu

button_left = glutCreateMenu(set_left_button)
call glutAddMenuEntry("rotate"//char(0),ROTATE)
call glutAddMenuEntry("zoom"//char(0),ZOOM)
call glutAddMenuEntry("pan"//char(0),PAN)
!call glutAddMenuEntry("scale x"//char(0),SCALEX)
!call glutAddMenuEntry("scale y"//char(0),SCALEY)
!call glutAddMenuEntry("scale z"//char(0), SCALEZ)
call glutAddMenuEntry("Section from minimal x plane"//char(0),8)
call glutAddMenuEntry("Section from maximal x plane"//char(0),9)
call glutAddMenuEntry("Section from minimal y plane"//char(0),10)
call glutAddMenuEntry("Section from maximal y plane"//char(0),11)
call glutAddMenuEntry("Section from minimal z plane"//char(0),12)
call glutAddMenuEntry("Section from maximal z plane"//char(0),13)

!call glutAddMenuEntry("proximal clipping plane"//char(0),14)
!call glutAddMenuEntry("distal clipping  plane"//char(0),15)
!call glutAddMenuEntry("change the scale of mouse for clipping"//char(0),16)

!button_middle = glutCreateMenu(set_middle_button)
!call glutAddMenuEntry("rotate"//char(0),ROTATE)
!call glutAddMenuEntry("zoom"//char(0),ZOOM)
!call glutAddMenuEntry("pan"//char(0),PAN)
!call glutAddMenuEntry("scale x"//char(0),SCALEX)
!call glutAddMenuEntry("scale y"//char(0),SCALEY)
!call glutAddMenuEntry("scale z"//char(0), SCALEZ)
arrow_keys = glutCreateMenu(set_arrow_keys)
!call glutAddMenuEntry("rotate"//char(0),ROTATE)
!call glutAddMenuEntry("zoom"//char(0),ZOOM)
!call glutAddMenuEntry("pan"//char(0),PAN)
!call glutAddMenuEntry("next step"//char(0),PAS)
!call glutAddMenuEntry("scale x"//char(0),SCALEX)
!call glutAddMenuEntry("scale y"//char(0),SCALEY)
!call glutAddMenuEntry("scale z"//char(0), SCALEZ)
menuid = glutCreateMenu(menu_handler)
call glutAddSubMenu("left mouse button"//char(0),button_left)
!call glutAddSubMenu("middle mouse button"//char(0),button_middle)
!call glutAddSubMenu("arrow keys"//char(0),arrow_keys)
call glutAddMenuEntry("reset to initial view"//char(0),RESET)
call glutAddMenuEntry("view from above"//char(0),ABOVE)
call glutAddMenuEntry("view from front"//char(0),FRONT)

call glutAddMenuEntry("Undo sections"//char(0),30)
call glutAddMenuEntry("Toggle Fixed sections when running iterations"//char(0),31)


!call glutAddMenuEntry("change perspective param"//char(0),pers)
!call glutAddMenuEntry("quit"//char(0),QUIT)

! set the perspective

call glMatrixMode(GL_PROJECTION)
call gluPerspective(10.0_gldouble, 1.0_gldouble, 0.1_gldouble, 200.0_gldouble)

! set the initial view

call glMatrixMode(GL_MODELVIEW)
call glPushMatrix
call reset_to_init

return
end function view_modifier_init
!*******************************************SUBROUTINE********************************************************
!        -----------
function sphere2cart(spoint) result(cpoint)
!        -----------
type(sphere3D), intent(in) :: spoint
type(cart3D) :: cpoint

! This converts a 3D point from spherical to cartesean coordinates

real(kind=gldouble) :: t,p,r

t=spoint%theta
p=spoint%phi
r=spoint%rho

cpoint%x = r*dsin(t)*dcos(p)
cpoint%y = r*dsin(t)*dsin(p)
cpoint%z = r*dcos(t) !>>Miquel28-10-14

return
end function sphere2cart
!*******************************************SUBROUTINE********************************************************
!        -----------
function cart2sphere(cpoint) result(spoint)
!        -----------
type(cart3D), intent(in) :: cpoint
type(sphere3D) :: spoint

! This converts a 3D point from cartesean to spherical coordinates

real(kind=gldouble) :: x,y,z

x=cpoint%x
y=cpoint%y
z=cpoint%z

spoint%rho = sqrt(x*x+y*y+z*z)
if (x==0.0_gldouble .and. y==0.0_gldouble) then
   spoint%theta = 0.0_gldouble
else
   spoint%theta = atan2(y,x)
end if
if (spoint%rho == 0.0_gldouble) then
   spoint%phi = 0.0_gldouble
else
   spoint%phi = acos(z/spoint%rho)
endif

return
end function cart2sphere

!        ------------------
function cart3D_plus_cart3D(cart1,cart2) result(cart3)
!        ------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the sum of two 3D cartesean points

cart3%x = cart1%x + cart2%x
cart3%y = cart1%y + cart2%y
cart3%z = cart1%z + cart2%z

return
end function cart3D_plus_cart3D

!        -------------------
function cart3D_minus_cart3D(cart1,cart2) result(cart3)
!        -------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the difference of two 3D cartesean points

cart3%x = cart1%x - cart2%x
cart3%y = cart1%y - cart2%y
cart3%z = cart1%z - cart2%z

return
end function cart3D_minus_cart3D

end module view_modifier

!---------------------------------------------------------------------------
!**********************************************************************************************************************************************************************
!              MODUL FUNCTION PLOTTER
!**********************************************************************************************************************************************************************
module function_plotter
!use opengl
use opengl_gl
use opengl_glu
use opengl_glut
use opengl_kinds
use general
use neighboring
use model
use view_modifier
use mitosis
use gnuplotter  !>>>Miguel8-10-14
!use del
use io
use editor
use iso_c_binding
!use gifer
!use fitmo

implicit none
private
public :: display,Reshape,draw_func,menu_handler,make_menu,inivisualitzacio,menuu_handler,menuu2_handler,menud_handler, &
          menuq_handler,menus_handler,menuse_handler, menueditor_handler, planemenu_handler,nodemenu_handler, &
          cellmenu_handler,nodecopymenu_handler,nodepastemenu_handler,cellcopymenu_handler,sel_menu_handler, &
          cursor_menu_handler,color_menu_handler,menuplot_handler,&
          cmenuu_handler,cmenuu2_handler,amenuu_handler,amenuu2_handler,smenuu_handler,smenuu2_handler,palette_menu_handler,&
          arrow_menu_handler,sphere_menu_handler, move_menu_handler


! symbolic constants
! symbolic constants
real(GLFLOAT) :: color(4), normal(3), &
                 red(4) = (/1.0,0.0,0.0,1.0/), &
                 gris(4) = (/0.5,0.5,0.5,1.0/), &
                 black(4) = (/0.0,0.0,0.0,1.0/), &
                 verd(4) = (/0.0,1.0,0.0,1.0/), &
                 verdd(4) = (/0.0,0.7,0.0,1.0/), &
                 blau(4) = (/0.0,0.0,1.0,1.0/), &
                 white(4) = (/1.0,1.0,1.0,1.0/)  

integer, parameter :: surfgrid_toggle = 1, &
                      surfsolid_toggle = 2, &
                      contour_toggle = 3, &
                      quit_selected = 4, &
                      knots_toggle=5, &
                      gu_toggle=6,gd_toggle=7,gt_toggle=8,gq_toggle=9,dif_toggle=10,vu_toggle=11, &
                      ab_toggle=12, miquel_toggle=13,un_toggle=14,deu_toggle=15,cent_toggle=16,mil_toggle=17,&
					  cubs_toggle=18,x_toggle=19,e_toggle=20,order_toggle=21,cell_toggle=22,mitosi_toggle=23,&
					  node_toggle=24,vb_toggle=25,vd_toggle=26,move_toggle=27,pol_toggle=28,no_toggle=29,&
                                          polside_toggle=30

integer, parameter :: set_nx = 1, &
                      set_ny = 2, &
                      set_ncontour = 3, &
                      set_contour_val = 4, &
                      set_xrange = 5, &
                      set_yrange = 6, &
                      reset_params = 7, &
                      p_size =8
                

integer, parameter :: black_contour = 1, &
                      rainbow_contour = 2

integer, parameter :: white_surface = 1, &
                      red_surface = 2, &
                      rainbow_surface = 3

! Default initial settings

integer, parameter :: init_ngridx = 40, &
                      init_ngridy = 40, &
                      init_num_contour = 20, &
                      init_contour_color = black_contour, &
                      init_surface_color = rainbow_surface

real(GLDOUBLE), parameter :: init_minx = 0.0_GLDOUBLE, &
                             init_maxx = 1.0_GLDOUBLE, &
                             init_miny = 0.0_GLDOUBLE, &
                             init_maxy = 1.0_GLDOUBLE

logical, parameter :: init_draw_surface_grid = .false., &
                      init_draw_surface_solid = .true., &
                      init_draw_contour = .true.
		     
		      	
integer :: ngridx = init_ngridx, &
           ngridy = init_ngridy, &
           num_contour = init_num_contour, &
           contour_color = init_contour_color, &
           surface_color = init_surface_color

real(GLDOUBLE) :: minx = init_minx, &
                  maxx = init_maxx, &
                  miny = init_miny, &
                  maxy = init_maxy, &
                  minz = 0.0_GLDOUBLE, &
                  maxz = 0.0_GLDOUBLE
real(GLFLOAT)  :: pzz   = 20.0_GLFLOAT

logical :: draw_surface_grid = init_draw_surface_grid, &
           draw_surface_solid = init_draw_surface_solid, &
           draw_contour = init_draw_contour, &
           contour_values_given = .false.
          

real(GLDOUBLE), allocatable :: actual_contours(:)

!variables de control de visualitzacio
!flags 
integer,private                ::flag(41)
                               !dim 1=elipse lines; 2=boxes ; 3=cell polygons ; 4=polygon sides
                               !5=balls req ; 6=balls da
integer,private                ::kii   !>>>Miguel 8-10-14 
integer                        ::geness(10),nodoss(10) !>>>Miguel 8-10-14 
real*8 ,private                ::maxe,mine
real*8 ,private                ::amaxe,amine,smaxe,smine !>>Miquel30-10-14


integer,private ::colorselection,property !Tommi 23.9.2013
integer,private ::arrowselection,sphereselection !>>Miquel30-10-14
integer,private ::chogen
integer,private ::achogen,schogen !>>Miquel31-10-14
real*8 ,private, allocatable   ::theslice(:,:)
integer,private ::nq,decoy
real*8, private ::mixi,miyi,mizi

!type(nod) :: tempnod !Tommi 23.9.2013
!type(cel) :: tempcel !Tommi 23.9.2013

integer, public :: fgif,fgifmovie
integer, public :: it_step,nu_it


integer,private:: fiite !>>> Is 24-10-14 

contains

!*******************************************SUBROUTINE********************************************************

subroutine inivisualitzacio

  arrowselection=conf_arrowselection
  sphereselection=conf_sphereselection

  colorselection=conf_colorselection
  iqpas=1
  flag(1)=conf_flag(1) !elipse lines ATTENTION 0!!! Tommi 
  flag(2)=conf_flag(2) !boxes
  !flag(3)=0 !selecting for node
  !flag(4)=0 !selecting for cells
  flag(5)=conf_flag(5) !balls req					!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
  flag(6)=conf_flag(6) !balls da					!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
  !flag(7)=0 !marked selected nodes
  flag(8)=conf_flag(8) !no balls
  flag(9)=conf_flag(9) !upper balls
  flag(10)=conf_flag(10)!lower balls 
  flag(11)=conf_flag(11)!small balls
  flag(12)=conf_flag(12)!connexions between cells
  flag(13)=conf_flag(13)!connexions within cells
  flag(18)=conf_flag(18)!main box
  flag(19)=conf_flag(19)!polarization vectors 
  flag(20)=0! 
  flag(21)=conf_flag(21)!centroids 
  flag(22)=conf_flag(22)!cylinders 
  flag(30)=conf_flag(30)!epithelium
  flag(31)=conf_flag(31)!mesenchyme
  flag(32)=conf_flag(32)!ecm
!  flag(33)=conf_flag(33)!epigrid
  !flag(34)=0
  flag(35)=conf_flag(35)!eggshell !miguel4-11-13
  flag(36)=conf_flag(36)!plot the "hold" nodes (node()%fix=1) !>>>>Miquel15-1-14
  flag(37)=conf_flag(37)!plot cell contour  !>>Miquel29-8-14
  flag(38)=conf_flag(38)!plot intercellular contour  !>>Miquel29-8-14
  flag(39)=conf_flag(39)!plot displacement of nodes respect i.c.  >>>>Miquel26-3-14
  !flag(41)=conf_flag(41)!dynamic/fixed display box !>>Miquel30-9-14
  !plotting window
 

  dyn=1 ; fix_run=0  !dynamic display box flag ; fixed dysplay box after running iterations flag !>>Miquel2-10-14

 
  set_mix=minval(node(:nd)%x) ; set_mx=maxval(node(:nd)%x)  !>>Miquel30-9-14
  !set_mix=0.15d0 ; set_mx=maxval(node(:nd)%x)                   !uncomment this to make a slice on the center

  set_miy=minval(node(:nd)%y) ; set_my=maxval(node(:nd)%y)
  !set_miy=-0.2d0 ; set_my=maxval(node(:nd)%y)                   !uncomment this to make a slice on the center
  set_miz=minval(node(:nd)%z) ; set_mz=maxval(node(:nd)%z)
  mesc=1d-3
  maan=0
  mapo=2500

  custom_colorselection=conf_custom_colorselection
  custom_minval=conf_custom_min ;custom_maxval=conf_custom_max
  chogen=colorselection-43

  amax_scale=conf_custom_arrowscale
  custom_arrowselection=conf_custom_arrowselection
  custom_aminval=conf_custom_aminval ; custom_amaxval=conf_custom_amaxval
  achogen=arrowselection-42

  smax_scale=conf_custom_spherescale
  custom_sphereselection=conf_custom_sphereselection
  custom_sminval=conf_custom_sminval ; custom_smaxval=conf_custom_smaxval
  schogen=sphereselection-42

  select case (conf_select_what) !>>Miquel5-11-14
  case(0)
    flag(7)=0
    flag(3)=0 ; flag(4)=0
    nki=0
  case(1)
    flag(7)=1
    flag(3)=1 ; flag(4)=0
    nki=conf_nki
  case(2)
    flag(7)=1
    flag(3)=0 ; flag(4)=1
    nki=conf_nki
  end select

  if(nki>0)then
    allocate(oopp(nki))
    oopp(1:nki)=conf_oopp(1:nki)
    if(flag(3)==1) nodeindex=oopp(nki)
    if(flag(4)==1) cellid=oopp(nki)
  end if


  nq=450
  if (allocated(theslice)) deallocate(theslice)
  allocate(theslice(0:nq,0:nq))
  theslice=0
  decoy=nd  !this should not be larger than nd since it is a node to visualize the field

end subroutine

!*******************************************SUBROUTINE********************************************************

subroutine Reshape(width, height) bind(c)    !>>>>>Miquel26-11-13
integer(glcint), value :: width, height
real(gldouble):: ratio
  windW = width
  windH = height


  ratio =  width * 1.0 / height;

  !Use the Projection Matrix
  call glMatrixMode(GL_PROJECTION);

  !Reset Matrix
  call glLoadIdentity();

  !Set the viewport to be the entire window
  call glViewport(0, 0, windW, windH);

  !Set the correct perspective.
  call gluPerspective(coolperspec, ratio, 0.1_gldouble, 1000.0_gldouble);

  !Get Back to the Modelview
  call glMatrixMode(GL_MODELVIEW);

!alternate (doesn't work)
!  call glViewport(0_GLint, 0_GLint, windW, windH)
!  call glMatrixMode(GL_PROJECTION)
!  call glLoadIdentity()
!  call gluOrtho2D(-0.5_gldouble, windW + 0.5_gldouble, -0.5_gldouble, windH + 0.5_gldouble)
!  call glMatrixMode(GL_MODELVIEW)
return
end subroutine Reshape



!*******************************************SUBROUTINE********************************************************


subroutine display
  integer irf
  character*100 org
!  call fit(1,a)
!print *,carg,"file"
!print *,nd,"nd"

  !if (fgifmovie>0) then
  !  call glutHideWindow()
  !  it_step=100
  !  nu_it=10
  !  call get_its(it_step,nu_it)
  !  do irf=1,nu_it
  !    call iteracio(it_step)
  !    call reset_view ; call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT)) ;  
  !    call glCallList(1)
  !    call draw_func;
  !    if (fgifmovie==2) then !we keeping images in the screen meanwhile
  !      call glutSwapBuffers ;  
  !      call glutPostRedisplay
  !    end if
  !    call savegif
  !  end do
!    print *,"this is the end my friend,  or do you want to more? (input 1 if so)"
!    read(*,*) i
!    if (i/=1) stop
  !   stop
  !end if
  call the_plotting
!>>> Is 24-10-14
  if (passed>0) then
!    passed=0
    if (fmovie==1) then
      if (passed==1) then
        call glutPostRedisplay
        call takeimage
        call iteracio(freqsnap)
        if (getot>fiite) then
          do i=1,5 ; print *,"" ; end do         
          print *,"WARNING WARNING WARNING WARNING WARNING WARNING"
          print *,""
          print *,"you need to have the command convert from imagemagick: if not make a sudo apt-get install imagemagick"
          print *,""
          print *,""
          print *,"the individual gif files are stored in folder image/"//trim(caa)
          print *,""
          print *,""
          print *,"the movie is stored in the current directory under the name "//caa
          print *,""
          print *,"if you do not have eog this will give and error but the movie is sitll there"
          print *,""
          print *,"WARNING WARNING WARNING WARNING WARNING WARNING"
          do i=1,5 ; print *,"" ; end do         
          org="convert -delay 50 images/"//trim(caa)//"/*.gif -loop 0 "//caa//".gif"
          print *,org
          call system(org)
          print *,"if you do not have eog this will give and error but the movie is sitll there"
          org="eog "//trim(caa)//"*.gif &"
          print *,org
          call system(org)
          call exit
        end if
      else
!        call savegif
      end if
    end if
  else
    if (fmovie==1) call iteracio(freqsnap)
  end if
!>>> Is 24-10-14
end subroutine

subroutine the_plotting
  call reset_view ; call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT)) ;  
  call glCallList(1)
  call draw_func;
  !if (fgif==1) then
!print *,"ki"
  !  call savegif
  !  stop
  !end if
  call glutSwapBuffers ;  
  call glutPostRedisplay
end subroutine

!*******************************************SUBROUTINE********************************************************

!>>> Is 25-10-14
subroutine takeimage
  character*150 pimp,mop
  character*30 nomfitx,nofi
    call glutPopWindow
    call glutShowWindow
    nomfitx=caa
    do i=1,30
      if (nomfitx(i:i)==" ") nomfitx(i:i)="_"
    end do
    write (nofi,*) getot
    nofi=adjustl(nofi)
    nofi(len_trim(nofi)+1:len_trim(nofi)+5)=".gif"
    print *,winame,"winame"
    pimp='import -window "'//winame(1:33)//'" kk.gif' 
!    pimp='import -window ""w_wre"" kk.gif' 
    !pimp(16:16)='"'
    !pimp(55:55)='"'
    pimp=trim(pimp)
    print *,pimp,"   ff"
    call system(pimp)
    mop="cp kk.gif images/"//trim(caa)//"/"//nomfitx//nofi
    print *,mop
    call system(mop) 
!>>> Is 25-10-14
end subroutine

!*******************************************SUBROUTINE********************************************************
subroutine draw_func
integer primer
real(GLDOUBLE) :: a,b,c,vamax,vamin,maxh,r,coo,cooc,aa,bb,cc,aaa,bbb,ccc
real(GLDOUBLE) :: ux,uy,uz,dx,dy,dz,tx,ty,tz,cu,cd,ct
real(GLFLOAT)  :: normal(3)
real*8         :: disx(nd),disy(nd),disz(nd),distot(nd),ecomp(nd,6),discentch(nd)
integer i,j,k,ii,jj,kk,iii,jjj,kkk,iiii,ic,dii,djj,dkk,diii,djjj,dkkk,icc,nnv,ik,iccc,p
integer sampling,zlevel
integer ll,mm,nn,lll,mmm,nnn    !>>>> Miquel 25-9-13
real*8  oldx,oldy,oldz

real(GLFLOAT)  :: memcolor(nd,4)  !>>>Miquel10-12-13
real(GLFLOAT)  :: memarrow(nd)  !>>>Miquel30-10-14
real(GLFLOAT)  :: memsphere(nd)  !>>>Miquel30-10-14

!!!triangulation variables !>>Miquel4-3-14
integer:: npt !number of points
integer:: sizht !size of hash table
integer:: maxbf,maxfc
integer::nbf,nfc,nface,ntetra,ierr   !size of arrays
real*8,dimension(:,:) :: vcl(3,nd) !point coordinates
real(GLDOUBLE)::al,bl,cl,as,bs,cs,ax,bx,cx
real(GLDOUBLE):: cpolx,cpoly,cpolz,myfactor,umyfactor
real(GLFLOAT)  :: anglech,  rvecx,rvecy,rvecz
real*8 :: mpol
integer :: eloch, ich, jch

esca=0.1d0
call reset_view
call glDeleteLists(1_gluint,1_glsizei)
call glNewList(1_gluint,gl_compile_and_execute)


  !pintem els eixos	que a la vegada faran d'escala
  call glBegin(gl_lines)

   !axis
   !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,blau)
   !call glvertex3d(0.0_gldouble,0.0_gldouble,2*esca)
   !call glvertex3d(esca,0.0_gldouble,2*esca)
   !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
   !call glvertex3d(0.0_gldouble,0.0_gldouble,2*esca)
   !call glvertex3d(0.0_gldouble,esca,2*esca)
   !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,verd)
   !call glvertex3d(0.0_gldouble,0.0_gldouble,2*esca)
   !call glvertex3d(0.0_gldouble,0.0_gldouble,3*esca)

  !display box

  if(dyn==1)then !dynamic box spanning the whole system
    set_mix=minval(node(:nd)%x) ; set_mx=maxval(node(:nd)%x)   !>>>>Miquel30-9-14
    set_miy=minval(node(:nd)%y) ; set_my=maxval(node(:nd)%y)
    set_miz=minval(node(:nd)%z) ; set_mz=maxval(node(:nd)%z)
  end if

  mix=set_mix ; mx=set_mx   !>>>>Miquel30-9-14
  miy=set_miy ; my=set_my
  miz=set_miz ; mz=set_mz

  if (flag(18)==1) then
!    select case(make_flag_twe_five)
!    case(0) !no selection plane
!      mix=set_mix ; mx=set_mx   !>>>>Miquel30-9-14
!      miy=set_miy ; my=set_my
!      miz=set_miz ; mz=set_mz
!    case(1)!selection plane is mx
!      mix=set_mix ;
!      miy=set_miy ; my=set_my
!      miz=set_miz ; mz=set_mz
!    case(2)!selection plane is mix
!      mx=set_mx
!      miy=set_miy ; my=set_my
!      miz=set_miz ; mz=set_mz
!    case(3)!selection plane is my
!      mix=set_mix ; mx=set_mx
!      miy=set_miy ;
!      miz=set_miz ; mz=set_mz
!    case(4)!selection plane is miy
!      mix=set_mix ; mx=set_mx
!      my=set_my
!      miz=set_miz ; mz=set_mz
!    case(5)!selection plane is mz
!      mix=set_mix ; mx=set_mx
!      miy=set_miy ; my=set_my
!      my=set_my ;
!    case(6)!selection plane miz
!      mix=set_mix ; mx=set_mx
!      miy=set_miy ; my=set_my
!      mz=set_mz
!    end select
         

    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,verd)
    call glvertex3d(mix*esca,miy*esca,miz*esca) ; call glvertex3d(mx*esca,miy*esca,miz*esca)
    call glvertex3d(mix*esca,miy*esca,miz*esca) ; call glvertex3d(mix*esca,my*esca,miz*esca)
    call glvertex3d(mix*esca,my*esca,mz*esca)   ; call glvertex3d(mix*esca,my*esca,miz*esca)

    call glvertex3d(mix*esca,my*esca,mz*esca)   ; call glvertex3d(mx*esca,my*esca,mz*esca)
    call glvertex3d(mx*esca,miy*esca,mz*esca)   ; call glvertex3d(mx*esca,my*esca,mz*esca)
    call glvertex3d(mx*esca,my*esca,mz*esca)    ; call glvertex3d(mx*esca,my*esca,miz*esca)
	
    call glvertex3d(mix*esca,my*esca,miz*esca)  ; call glvertex3d(mx*esca,my*esca,miz*esca)
    call glvertex3d(mix*esca,miy*esca,mz*esca)  ; call glvertex3d(mx*esca,miy*esca,mz*esca)
    call glvertex3d(mix*esca,my*esca,mz*esca)   ; call glvertex3d(mix*esca,miy*esca,mz*esca)

    call glvertex3d(mx*esca,miy*esca,miz*esca)  ; call glvertex3d(mx*esca,my*esca,miz*esca)
    call glvertex3d(mx*esca,miy*esca,mz*esca)   ; call glvertex3d(mx*esca,miy*esca,miz*esca)
    call glvertex3d(mix*esca,my*esca,mz*esca)   ; call glvertex3d(mix*esca,my*esca,miz*esca)

    call glvertex3d(mix*esca,miy*esca,mz*esca)   ; call glvertex3d(mix*esca,miy*esca,miz*esca)
  end if !<<<<<<<Tommi

  !polarization vectors
  if (flag(19)==1) then			!>>>>>>>>>>>>>>>>>>>>>>>>>> Miquel 3-6-13
    do i=1,ncels
      k=cels(i)%node(1)
      if(node(k)%fix==1 .and. flag(36)==0) cycle !>>Miquel22-1-14
      if(node(k)%tipus==1.and.flag(9)==0) cycle
      if(node(k)%tipus==2.and.flag(10)==0) cycle
      if(node(k)%tipus<3.and.flag(30)==0) cycle
      if(node(k)%tipus==3.and.flag(31)==0) cycle
      if(node(k)%tipus>3.and.flag(32)==0) cycle
      a=cels(i)%cex ; b=cels(i)%cey ; c=cels(i)%cez
      !if(single==1)then; a=node(i)%x ; b=node(i)%y ; c=node(i)%z ;end if
      aa=a+cels(i)%polx ; bb=b+cels(i)%poly ; cc=c+cels(i)%polz
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
        call glvertex3d(a*esca,b*esca,c*esca)
        call glvertex3d(aa*esca,bb*esca,cc*esca)
      end if
    end do
  end if

  !movement vectors
  if (flag(23)==1.or.flag(27)==1) then			!>>>>>>>>>>>>>>>>>>>>>>>>>> Miquel 18-6-13
    a=0.0d0
    do i=1,nd
      if (dex(i)>a) then ; a=dex(i) ; ii=i ; end if
    end do
    aaa=20d1*esca/a !this way the vectors escalate so the greatest force is plotted with a constant length
    do i=1,nd
      if (node(i)%fix==1 .and. flag(36)==0) cycle !>>Miquel22-1-14
      if (flag(30)==0.and.node(i)%tipus<3) cycle
      if (flag(31)==0.and.node(i)%tipus==3) cycle
      if (flag(32)==0.and.node(i)%tipus>3) cycle
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      if (flag(27)==1) then 
        aa=a+px(i)*aaa ; bb=b+py(i)*aaa ; cc=c+pz(i)*aaa
      else
        bbb=1d0/sqrt(px(i)**2+py(i)**2+pz(i)**2)
        aa=a+px(i)*bbb ; bb=b+py(i)*bbb ; cc=c+pz(i)*bbb
      end if
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
        call glvertex3d(a*esca,b*esca,c*esca)
        call glvertex3d(aa*esca,bb*esca,cc*esca)
      end if
    end do
  end if

  !force components: repulsion-adhesion
  if (flag(24)==1) then			!>>>>>>>>>>>>>>>>>>>>>>>>>> Miquel 20-6-13
    a=0.0d0
    do i=1,nd
      aa=vcilx(i) ; bb=vcily(i) ; cc=vcilz(i)
      d=sqrt(aa**2+bb**2+cc**2)
      if (d>a) then ; a=d ; ii=i ; end if
    end do
    aaa=20d0*esca/a !this way the vectors escalate so the greatest force is plotted with a constant length
    color=0 ; color(1)=1 ; color(2)=1
    do i=1,nd
      if(node(i)%tipus==1.and.flag(9)==0) cycle
      if(node(i)%tipus==2.and.flag(10)==0) cycle
      if(node(i)%tipus<3.and.flag(30)==0) cycle
      if(node(i)%tipus==3.and.flag(31)==0) cycle
      if(node(i)%tipus>3.and.flag(32)==0) cycle
      if(node(i)%fix==1 .and. flag(36)==0) cycle !>>Miquel22-1-14
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      aa=a+vcilx(i)*aaa ; bb=b+vcily(i)*aaa ; cc=c+vcilz(i)*aaa
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
        call glvertex3d(a*esca,b*esca,c*esca)
        call glvertex3d(aa*esca,bb*esca,cc*esca)
      end if
    end do
  end if

  !force components: epithelial surface tension lateral
  if (flag(25)==1) then			!>>>>>>>>>>>>>>>>>>>>>>>>>> Miquel 20-6-13
    a=0.0d0
    do i=1,nd
      aa=vtorx(i) ; bb=vtory(i) ; cc=vtorz(i)
      d=sqrt(aa**2+bb**2+cc**2)
      if (d>a) then ; a=d ; ii=i ; end if
    end do
    aaa=20d0*esca/a !this way the vectors escalate so the greatest force is plotted with a constant length
    color=0 ; color(1)=0.5 ; color(2)=1
    do i=1,nd
      if(node(i)%tipus==1.and.flag(9)==0) cycle
      if(node(i)%tipus==2.and.flag(10)==0) cycle
      if(node(i)%tipus<3.and.flag(30)==0) cycle
      if(node(i)%tipus==3.and.flag(31)==0) cycle
      if(node(i)%tipus>3.and.flag(32)==0) cycle
      if(node(i)%fix==1 .and. flag(36)==0) cycle !>>Miquel22-1-14
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      aa=a+vtorx(i)*aaa ; bb=b+vtory(i)*aaa ; cc=c+vtorz(i)*aaa
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
        call glvertex3d(a*esca,b*esca,c*esca)
        call glvertex3d(aa*esca,bb*esca,cc*esca)
      end if
    end do
  end if

  !force components: epithelial surface tension apical-basal
  if (flag(26)==1) then			!>>>>>>>>>>>>>>>>>>>>>>>>>> Miquel 20-6-13
    a=0.0d0
    do i=1,nd
      aa=vstorx(i) ; bb=vstory(i) ; cc=vstorz(i)
      d=sqrt(aa**2+bb**2+cc**2)
      if (d>a) then ; a=d ; ii=i ; end if
    end do
    aaa=20d0*esca/a !this way the vectors escalate so the greatest force is plotted with a constant length
    color=0 ; color(1)=1 ; color(2)=0
    do i=1,nd
      if(node(i)%tipus==1.and.flag(9)==0) cycle
      if(node(i)%tipus==2.and.flag(10)==0) cycle
      if(node(i)%tipus<3.and.flag(30)==0) cycle
      if(node(i)%tipus==3.and.flag(31)==0) cycle
      if(node(i)%tipus>3.and.flag(32)==0) cycle
      if(node(i)%fix==1 .and. flag(36)==0) cycle !>>Miquel22-1-14
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      aa=a+vstorx(i)*aaa ; bb=b+vstory(i)*aaa ; cc=c+vstorz(i)*aaa
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
        call glvertex3d(a*esca,b*esca,c*esca)
        call glvertex3d(aa*esca,bb*esca,cc*esca)
      end if
    end do
  end if
  
  !%fix nodes' spring
  color=0 ; color(1)=1
  do i=1,nd
    if(node(i)%tipus==1.and.flag(9)==0) cycle
    if(node(i)%tipus==2.and.flag(10)==0) cycle
    if(node(i)%tipus<3.and.flag(30)==0) cycle
    if(node(i)%tipus==3.and.flag(31)==0) cycle
    if(node(i)%tipus>3.and.flag(32)==0) cycle
    if((flag(36)==1.and.node(i)%fix==1).or.(flag(39)==1.and.node(i)%fix==0))then
      a=node(i)%x ; b=node(i)%y ; c=node(i)%z
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        aa=node(i)%orix ; bb=node(i)%oriy ; cc=node(i)%oriz
        !if (aa>=mix.and.aa<=mx.and.bb>=miy.and.bb<=my.and.cc>=miz.and.cc<=mz) then
          call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
          call glvertex3d(a*esca,b*esca,c*esca)
          call glvertex3d(aa*esca,bb*esca,cc*esca)
        !end if
      end if
    end if
  end do
  
  call glEnd

   !!!!!!!!!!!!!!!!!!!!!! miguel4-11-13 eggshell   
  if (flag(35)==1) then                                                         
    ! call shellplot
    do i=1,int(pts*pts*2) 
      ii=trs(i,1) ; jj=trs(i,2) ; kk=trs(i,3)         
      a=shellp(ii,1)  ; b=shellp(ii,2)  ;c=shellp(ii,3)
      aa=shellp(jj,1) ; bb=shellp(jj,2) ;cc=shellp(jj,3)
      aaa=shellp(kk,1); bbb=shellp(kk,2);ccc=shellp(kk,3)
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
        if (aa>=mix.and.aa<=mx.and.bb>=miy.and.bb<=my.and.cc>=miz.and.cc<=mz) then
          if (aaa>=mix.and.aaa<=mx.and.bbb>=miy.and.bbb<=my.and.ccc>=miz.and.ccc<mz) then
            call glBegin(gl_triangles)		!MAKING THE SURFACE OF THE TRIANGLES
             cooc=ncels
             normal=normcrossprod((/aaa,a,aa/),(/bbb,b,bb/),(/ccc,c,cc/)) 
             normal=abs(normal)
             call glnormal3fv(normal)

             coo=i
             call get_rainbow(coo,1d0,cooc,color)
             call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
             call glvertex3d(a*esca,b*esca,c*esca)

             coo=i
             call get_rainbow(coo,1d0,cooc,color)
             call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
             call glvertex3d(aa*esca,bb*esca,cc*esca)

             coo=i
             call get_rainbow(coo,1d0,cooc,color)
             call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
             call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
            call glEnd
          end if
        end if
      end if
    end do						  
  end if									
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (flag(1)==1) then  !drawing epithelial springs
    do i=1,nd
      if(node(i)%tipus<3.and.flag(30)==0) cycle
      if(node(i)%fix==1 .and. flag(36)==0) cycle !>>>Miquel17-1-14
      call glBegin(gl_lines)
      if(node(i)%tipus<3)then						!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
!!print *,i,node(i)%altre,node(i)%tipus
        a=node(i)%x;b=node(i)%y;c=node(i)%z
        ii=node(i)%altre
        aa=node(ii)%x;bb=node(ii)%y;cc=node(ii)%z
        if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
          if (aa>=mix.and.aa<=mx.and.bb>=miy.and.bb<=my.and.cc>=miz.and.cc<=mz) then
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,verd)
            call glvertex3d(a*esca,b*esca,c*esca)
            call glvertex3d(aa*esca,bb*esca,cc*esca)
          endif
        end if
      end if
      call glend
    end do
  endif

  !GRID BOXES
  if(flag(2)==1)then
    !pintarem els boxes
    color=1
    d=rv
    call glBegin(gl_lines)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    do j=-nboxes-1,nboxes															!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
      do k=-nboxes-1,nboxes														!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        a=-d*(nboxes+0.5d0);aa=d*(nboxes+0.5d0)										!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        b=d*(j+0.5d0)
        c=d*(k+0.5d0)
        call glvertex3d(esca*a,esca*b,esca*c)
        call glvertex3d(esca*aa,esca*b,esca*c)
      end do
    end do
    do i=-nboxes-1,nboxes															!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
      do k=-nboxes-1,nboxes														!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        a=d*(i+0.5d0)
        b=-d*(nboxes+0.5d0);bb=d*(nboxes+0.5d0)										!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        c=d*(k+0.5d0)
        call glvertex3d(esca*a,esca*b,c*esca)
        call glvertex3d(esca*a,esca*bb,c*esca)
      end do
    end do
    do i=-nboxes-1,nboxes															!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
      do j=-nboxes-1,nboxes														!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        a=d*(i+0.5d0)
        b=d*(j+0.5d0)
        c=-d*(nboxes+0.5d0);cc=d*(nboxes+0.5d0)										!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
        call glvertex3d(esca*a,esca*b,esca*c)
        call glvertex3d(esca*a,esca*b,esca*cc)
      end do
    end do
    call glend
  end if

  !BALLS
  r=0.05*esca
 
 
  do i=1,nd                                                                      !!>> HC 12-2-2024
     discentch(i)=(sqrt( (node(i)%x)**2 + (node(i)%y)**2 + (node(i)%z)**2))**2   !!>> HC 12-2-2024
  enddo                                                                          !!>> HC 12-2-2024
  !BALL COLORS MAX AND MINS

  if(custom_colorselection==colorselection)then
    mine=custom_minval ; maxe=custom_maxval
  else
    select case(colorselection)	
    case(1) ; mine=minval(node(:nd)%x) ; maxe=maxval(node(:nd)%x)
    case(2) ; mine=minval(node(:nd)%y) ; maxe=maxval(node(:nd)%y)
    case(3) ; mine=minval(node(:nd)%z) ; maxe=maxval(node(:nd)%z)
    case(4) ; mine=minval(node(:nd)%e) ; maxe=maxval(node(:nd)%e)
    case(5) ; mine=minval(node(:nd)%eqd) ; maxe=maxval(node(:nd)%eqd) 
    case(6) ; mine=minval(node(:nd)%add)  ; maxe=maxval(node(:nd)%add)
    case(7) ; mine=minval(node(:nd)%you) ; maxe=maxval(node(:nd)%you)
    case(8) ; mine=minval(node(:nd)%adh) ; maxe=maxval(node(:nd)%adh)
    case(9) ; mine=minval(node(:nd)%rep) ; maxe=maxval(node(:nd)%rep)	
    case(10) ; mine=minval(node(:nd)%rec) ; maxe=maxval(node(:nd)%rec)
    case(11) ; mine=minval(node(:nd)%erp) ; maxe=maxval(node(:nd)%erp)
    case(12) ; mine=minval(node(:nd)%est) ; maxe=maxval(node(:nd)%est)	
    case(13) ; mine=minval(node(:nd)%eqs) ; maxe=maxval(node(:nd)%eqs)
    case(14) ; mine=minval(node(:nd)%hoo) ; maxe=maxval(node(:nd)%hoo) 
    case(15) ; mine=minval(node(:nd)%mov) ; maxe=maxval(node(:nd)%mov)
    case(16) ; mine=minval(node(:nd)%dmo) ; maxe=maxval(node(:nd)%dmo)
    case(17) ; mine=minval(node(:nd)%orix) ; maxe=maxval(node(:nd)%orix) 
    case(18) ; mine=minval(node(:nd)%oriy) ; maxe=maxval(node(:nd)%oriy)
    case(19) ; mine=minval(node(:nd)%oriz) ; maxe=maxval(node(:nd)%oriz)
    case(20) ; mine=minval(node(:nd)%ecm) ; maxe=maxval(node(:nd)%ecm)
    case(21) ; mine=minval(node(:nd)%cod) ; maxe=maxval(node(:nd)%cod)
    case(22) ; mine=minval(node(:nd)%grd) ; maxe=maxval(node(:nd)%grd)
    case(23) ; mine=minval(node(:nd)%pld) ; maxe=maxval(node(:nd)%pld)
    case(24) ; mine=minval(node(:nd)%vod) ; maxe=maxval(node(:nd)%vod)
    case(25) ; mine=minval(node(:nd)%dif) ; maxe=maxval(node(:nd)%dif)
    case(26) ; mine=minval(node(:nd)%kfi) ; maxe=maxval(node(:nd)%kfi)
    case(27) ; mine=minval(node(:nd)%pla) ; maxe=maxval(node(:nd)%pla)
    case(28) ; mine=minval(node(:nd)%kvol) ; maxe=maxval(node(:nd)%kvol)
    case(29) ; mine=1 ;maxe=7 !mine=minval(node(:nd)%tipus) ; maxe=maxval(node(:nd)%tipus) !this is tipus
    case(30) ; mine=minval(node(:nd)%icel) ; maxe=maxval(node(:nd)%icel)
    case(31) ; mine=minval(node(:nd)%altre) ; maxe=maxval(node(:nd)%altre)
    case(32) ; mine=minval(node(:nd)%marge) ; maxe=maxval(node(:nd)%marge) 
    case(33) ; mine=minval(node(:nd)%talone) ; maxe=maxval(node(:nd)%talone) 
    case(34) ; mine=0;maxe=1!minval(node(:nd)%fix) ; maxe=maxval(node(:nd)%fix) 
    case(35) ; mine=minval(node(:nd)%ndiv) ; maxe=maxval(node(:nd)%ndiv)        !!>> HC 16-3-2022
    case(36) ; mine=minval(discentch(:nd)) ; maxe=maxval(discentch(:nd))        !!>> HC 12-2-2024
    !case(36) ; if (chogen==0) then
    !           chogen=1
    !           mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !           print*,"maxim",maxe,"minim",mine; end if
    !           mine=minval(gex(:nd,1)) ; maxe=maxval(gex(:nd,1))
    !case(37) ; if (chogen==0) then
    !             if (ng>1) then
    !               chogen=2
    !                mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !               print*,"maxim",maxe,"minim",mine; 
    !             else
    !               print *,"WARNING, not such a gene"
    !             end if
    !           else
    !             mine=minval(gex(:nd,2)) ; maxe=maxval(gex(:nd,2))
    !           end if
    !case(38) ; if (chogen==0) then
    !           print *,"which gene" ; read(*,*) chogen ; print *,chogen,"chosen regulatory molecule"
    !           if (chogen>ng) then 
    !             print *,"WARNING, not such a gene"
    !           else
    !             mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !             print*,"maxim",maxe,"minim",mine; 
    !           end if
    !             mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !           end if
    case(37) ; mine=0d0 ; maxe=1d0
    case(38) ; do i=1,nd ; disx(i)=node(i)%x-node(i)%orix ; end do 
               mine=minval(disx) ; maxe=maxval(disx)
    case(39) ; do i=1,nd ; disy(i)=node(i)%y-node(i)%oriy ; end do 
               mine=minval(disy) ; maxe=maxval(disy)
    case(40) ; do i=1,nd ; disz(i)=node(i)%z-node(i)%oriz ; end do 
               mine=minval(disz) ; maxe=maxval(disz)
!    case(40) ; !do i=1,nd 
               !distot(i)=sqrt((node(i)%x-node(i)%orix)**2+(node(i)%y-node(i)%oriy)**2+(node(i)%z-node(i)%oriz)**2) 
               !end do 
               !mine=minval(distot) ; maxe=maxval(distot)
    case(41) ; do i=1,nd 
               distot(i)=sqrt(px(i)**2+py(i)**2+pz(i)**2) 
               end do 
               mine=minval(distot) ; maxe=maxval(distot)


    case(42) ; mine=1 ; maxe=nd
    case(43) ; mine=1 ; maxe=nd
    case(44) ; mine=1; maxe=maxval(nodeo(:nd)%icel)  !!>> HC 15-12-2021
    end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !>>>Miguel30-10-14 
    if(colorselection.gt.44)then
      chogen=colorselection-44 
      mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
      !print*,"chosen gene",chogen,"maxim",maxe,"minim",mine
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if


  !ARROW MAX AND MINS
  if(custom_arrowselection==arrowselection)then
    amine=custom_aminval ; amaxe=custom_amaxval
  else
    select case(arrowselection)	
    case(1) ; amine=minval(node(:nd)%x) ; amaxe=maxval(node(:nd)%x)
    case(2) ; amine=minval(node(:nd)%y) ; amaxe=maxval(node(:nd)%y)
    case(3) ; amine=minval(node(:nd)%z) ; amaxe=maxval(node(:nd)%z)
    case(4) ; amine=minval(node(:nd)%e) ; amaxe=maxval(node(:nd)%e)
    case(5) ; amine=minval(node(:nd)%eqd) ; amaxe=maxval(node(:nd)%eqd) 
    case(6) ; amine=minval(node(:nd)%add)  ; amaxe=maxval(node(:nd)%add)
    case(7) ; amine=minval(node(:nd)%you) ; amaxe=maxval(node(:nd)%you)
    case(8) ; amine=minval(node(:nd)%adh) ; amaxe=maxval(node(:nd)%adh)
    case(9) ; amine=minval(node(:nd)%rep) ; amaxe=maxval(node(:nd)%rep)	
    case(10) ; amine=minval(node(:nd)%rec) ; amaxe=maxval(node(:nd)%rec)
    case(11) ; amine=minval(node(:nd)%erp) ; amaxe=maxval(node(:nd)%erp)
    case(12) ; amine=minval(node(:nd)%est) ; amaxe=maxval(node(:nd)%est)	
    case(13) ; amine=minval(node(:nd)%eqs) ; amaxe=maxval(node(:nd)%eqs)
    case(14) ; amine=minval(node(:nd)%hoo) ; amaxe=maxval(node(:nd)%hoo) 
    case(15) ; amine=minval(node(:nd)%mov) ; amaxe=maxval(node(:nd)%mov)
    case(16) ; amine=minval(node(:nd)%dmo) ; amaxe=maxval(node(:nd)%dmo)
    case(17) ; amine=minval(node(:nd)%orix) ; amaxe=maxval(node(:nd)%orix) 
    case(18) ; amine=minval(node(:nd)%oriy) ; amaxe=maxval(node(:nd)%oriy)
    case(19) ; amine=minval(node(:nd)%oriz) ; amaxe=maxval(node(:nd)%oriz)
    case(20) ; amine=minval(node(:nd)%ecm) ; amaxe=maxval(node(:nd)%ecm)
    case(21) ; amine=minval(node(:nd)%cod) ; amaxe=maxval(node(:nd)%cod)
    case(22) ; amine=minval(node(:nd)%grd) ; amaxe=maxval(node(:nd)%grd)
    case(23) ; amine=minval(node(:nd)%pld) ; amaxe=maxval(node(:nd)%pld)
    case(24) ; amine=minval(node(:nd)%vod) ; amaxe=maxval(node(:nd)%vod)
    case(25) ; amine=minval(node(:nd)%dif) ; amaxe=maxval(node(:nd)%dif)
    case(26) ; amine=minval(node(:nd)%kfi) ; amaxe=maxval(node(:nd)%kfi)
    case(27) ; amine=minval(node(:nd)%pla) ; amaxe=maxval(node(:nd)%pla)
    case(28) ; amine=minval(node(:nd)%kvol) ; amaxe=maxval(node(:nd)%kvol)
    case(29) ; amine=1 ;amaxe=7 !amine=minval(node(:nd)%tipus) ; amaxe=maxval(node(:nd)%tipus) !this is tipus
    case(30) ; amine=minval(node(:nd)%icel) ; amaxe=maxval(node(:nd)%icel)
    case(31) ; amine=minval(node(:nd)%altre) ; amaxe=maxval(node(:nd)%altre)
    case(32) ; amine=minval(node(:nd)%marge) ; amaxe=maxval(node(:nd)%marge) 
    case(33) ; amine=minval(node(:nd)%talone) ; amaxe=maxval(node(:nd)%talone) 
    case(34) ; amine=0;amaxe=1!minval(node(:nd)%fix) ; amaxe=maxval(node(:nd)%fix) 
    case(35) ; amine=minval(node(:nd)%ndiv) ; amaxe=maxval(node(:nd)%ndiv)   !!>> HC 16-3-2022
    case(36) ; mine=minval(discentch(:nd)) ; maxe=maxval(discentch(:nd))        !!>> HC 12-2-2024
    !case(36) ; if (chogen==0) then
    !           chogen=1
    !           amine=minval(gex(:nd,chogen)) ; amaxe=maxval(gex(:nd,chogen))
    !           print*,"maxim",amaxe,"minim",amine; end if
    !           amine=minval(gex(:nd,1)) ; amaxe=maxval(gex(:nd,1))
    !case(37) ; if (chogen==0) then
    !             if (ng>1) then
    !               chogen=2
    !               amine=minval(gex(:nd,chogen)) ; amaxe=maxval(gex(:nd,chogen))
    !               print*,"maxim",amaxe,"minim",amine; 
    !             else
    !               print *,"WARNING, not such a gene"
    !             end if
    !           else
    !             amine=minval(gex(:nd,2)) ; amaxe=maxval(gex(:nd,2))
    !           end if
    !case(38) ; if (chogen==0) then
    !           print *,"which gene" ; read(*,*) chogen ; print *,chogen,"chosen regulatory molecule"
    !           if (chogen>ng) then 
    !             print *,"WARNING, not such a gene"
    !           else
    !             amine=minval(gex(:nd,chogen)) ; amaxe=maxval(gex(:nd,chogen))
    !             print*,"maxim",amaxe,"minim",amine; 
    !           end if
    !             amine=minval(gex(:nd,chogen)) ; amaxe=maxval(gex(:nd,chogen))
    !           end if
    case(37) ; amine=0d0 ; amaxe=1d0
    case(38) ; do i=1,nd ; disx(i)=node(i)%x-node(i)%orix ; end do 

               amine=minval(disx) ; amaxe=maxval(disx)
    case(39) ; do i=1,nd ; disy(i)=node(i)%y-node(i)%oriy ; end do 
               amine=minval(disy) ; amaxe=maxval(disy)
    case(40) ; do i=1,nd ; disz(i)=node(i)%z-node(i)%oriz ; end do 
               amine=minval(disz) ; amaxe=maxval(disz)
    case(41) ; do i=1,nd 
               distot(i)=sqrt((node(i)%x-node(i)%orix)**2+(node(i)%y-node(i)%oriy)**2+(node(i)%z-node(i)%oriz)**2) 
               end do 
               amine=minval(distot) ; amaxe=maxval(distot)
    case(42) ; amine=1 ; amaxe=nd
    case(43) ; amine=1 ; amaxe=nd
    end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !>>>Miguel30-10-14 
    if(arrowselection.gt.43)then
      achogen=arrowselection-43 
      amine=minval(gex(:nd,achogen)) ; amaxe=maxval(gex(:nd,achogen))
      !print*,"chosen gene",chogen,"maxim",maxe,"minim",mine
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endif

  !TRANSPARENT SPHERES MAX AND MINS
  if(custom_sphereselection==sphereselection)then
    smine=custom_sminval ; smaxe=custom_smaxval
  else
    select case(sphereselection)	
    case(1) ; smine=minval(node(:nd)%x) ; smaxe=maxval(node(:nd)%x)
    case(2) ; smine=minval(node(:nd)%y) ; smaxe=maxval(node(:nd)%y)
    case(3) ; smine=minval(node(:nd)%z) ; smaxe=maxval(node(:nd)%z)
    case(4) ; smine=minval(node(:nd)%e) ; smaxe=maxval(node(:nd)%e)
    case(5) ; smine=minval(node(:nd)%eqd) ; smaxe=maxval(node(:nd)%eqd) 
    case(6) ; smine=minval(node(:nd)%add)  ; smaxe=maxval(node(:nd)%add)
    case(7) ; smine=minval(node(:nd)%you) ; smaxe=maxval(node(:nd)%you)
    case(8) ; smine=minval(node(:nd)%adh) ; smaxe=maxval(node(:nd)%adh)
    case(9) ; smine=minval(node(:nd)%rep) ; smaxe=maxval(node(:nd)%rep)	
    case(10) ; smine=minval(node(:nd)%rec) ; smaxe=maxval(node(:nd)%rec)
    case(11) ; smine=minval(node(:nd)%erp) ; smaxe=maxval(node(:nd)%erp)
    case(12) ; smine=minval(node(:nd)%est) ; smaxe=maxval(node(:nd)%est)	
    case(13) ; smine=minval(node(:nd)%eqs) ; smaxe=maxval(node(:nd)%eqs)
    case(14) ; smine=minval(node(:nd)%hoo) ; smaxe=maxval(node(:nd)%hoo) 
    case(15) ; smine=minval(node(:nd)%mov) ; smaxe=maxval(node(:nd)%mov)
    case(16) ; smine=minval(node(:nd)%dmo) ; smaxe=maxval(node(:nd)%dmo)
    case(17) ; smine=minval(node(:nd)%orix) ; smaxe=maxval(node(:nd)%orix) 
    case(18) ; smine=minval(node(:nd)%oriy) ; smaxe=maxval(node(:nd)%oriy)
    case(19) ; smine=minval(node(:nd)%oriz) ; smaxe=maxval(node(:nd)%oriz)
    case(20) ; smine=minval(node(:nd)%ecm) ; smaxe=maxval(node(:nd)%ecm)
    case(21) ; smine=minval(node(:nd)%cod) ; smaxe=maxval(node(:nd)%cod)
    case(22) ; smine=minval(node(:nd)%grd) ; smaxe=maxval(node(:nd)%grd)
    case(23) ; smine=minval(node(:nd)%pld) ; smaxe=maxval(node(:nd)%pld)
    case(24) ; smine=minval(node(:nd)%vod) ; smaxe=maxval(node(:nd)%vod)
    case(25) ; smine=minval(node(:nd)%dif) ; smaxe=maxval(node(:nd)%dif)
    case(26) ; smine=minval(node(:nd)%kfi) ; smaxe=maxval(node(:nd)%kfi)
    case(27) ; smine=minval(node(:nd)%pla) ; smaxe=maxval(node(:nd)%pla)
    case(28) ; smine=minval(node(:nd)%kvol) ; smaxe=maxval(node(:nd)%kvol)
    case(29) ; smine=1 ;smaxe=7 !smine=minval(node(:nd)%tipus) ; smaxe=maxval(node(:nd)%tipus) !this is tipus
    case(30) ; smine=minval(node(:nd)%icel) ; smaxe=maxval(node(:nd)%icel)
    case(31) ; smine=minval(node(:nd)%altre) ; smaxe=maxval(node(:nd)%altre)
    case(32) ; smine=minval(node(:nd)%marge) ; smaxe=maxval(node(:nd)%marge) 
    case(33) ; smine=minval(node(:nd)%talone) ; smaxe=maxval(node(:nd)%talone) 
    case(34) ; smine=0;smaxe=1!minval(node(:nd)%fix) ; smaxe=maxval(node(:nd)%fix) 
    case(35) ; smine=minval(node(:nd)%ndiv) ; smaxe=maxval(node(:nd)%ndiv)   !!>> HC 16-3-2022
    case(36) ; mine=minval(discentch(:nd)) ; maxe=maxval(discentch(:nd))        !!>> HC 12-2-2024
    !case(36) ; if (chogen==0) then
    !           chogen=1
    !           smine=minval(gex(:nd,chogen)) ; smaxe=maxval(gex(:nd,chogen))
    !           print*,"maxim",smaxe,"minim",smine; end if
    !           smine=minval(gex(:nd,1)) ; smaxe=maxval(gex(:nd,1))
    !case(37) ; if (chogen==0) then
    !             if (ng>1) then
    !               chogen=2
    !               smine=minval(gex(:nd,chogen)) ; smaxe=maxval(gex(:nd,chogen))
    !               print*,"maxim",smaxe,"minim",smine; 
    !             else
    !               print *,"WARNING, not such a gene"
    !             end if
    !           else
    !             smine=minval(gex(:nd,2)) ; smaxe=maxval(gex(:nd,2))
    !           end if
    !case(38) ; if (chogen==0) then
    !           print *,"which gene" ; read(*,*) chogen ; print *,chogen,"chosen regulatory molecule"
    !           if (chogen>ng) then 
    !             print *,"WARNING, not such a gene"
    !           else
    !             smine=minval(gex(:nd,chogen)) ; smaxe=maxval(gex(:nd,chogen))
    !             print*,"maxim",smaxe,"minim",smine; 
    !           end if
    !             smine=minval(gex(:nd,chogen)) ; smaxe=maxval(gex(:nd,chogen))
    !           end if
    case(37) ; smine=0d0 ; smaxe=1d0
    case(38) ; do i=1,nd ; disx(i)=node(i)%x-node(i)%orix ; end do 
               smine=minval(disx) ; smaxe=maxval(disx)
    case(39) ; do i=1,nd ; disy(i)=node(i)%y-node(i)%oriy ; end do 
               smine=minval(disy) ; smaxe=maxval(disy)
    case(40) ; do i=1,nd ; disz(i)=node(i)%z-node(i)%oriz ; end do 
               smine=minval(disz) ; smaxe=maxval(disz)
    case(41) ; do i=1,nd 
               distot(i)=sqrt((node(i)%x-node(i)%orix)**2+(node(i)%y-node(i)%oriy)**2+(node(i)%z-node(i)%oriz)**2) 
               end do 
               smine=minval(distot) ; smaxe=maxval(distot)
    case(42) ; smine=1 ; smaxe=nd
    case(43) ; smine=1 ; smaxe=nd
    end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !>>>Miguel30-10-14 
    if(sphereselection.gt.43)then
      schogen=sphereselection-43
      smine=minval(gex(:nd,schogen)) ; smaxe=maxval(gex(:nd,schogen))
      !print*,"chosen gene",chogen,"maxim",maxe,"minim",mine
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if


  do i=1,nd
    !BALL COLOR VALUE
    color=0
    select case (colorselection)
    case(1) ; a=node(i)%x
    case(2) ; a=node(i)%y
    case(3) ; a=node(i)%z
    case(4) ; a=node(i)%e
    case(5) ; a=node(i)%eqd
    case(6) ; a=node(i)%add
    case(7) ; a=node(i)%you
    case(8) ; a=node(i)%adh
    case(9) ; a=node(i)%rep
    case(10) ; a=node(i)%rec
    case(11) ; a=node(i)%erp
    case(12) ; a=node(i)%est
    case(13) ; a=node(i)%eqs
    case(14) ; a=node(i)%hoo
    case(15) ; a=node(i)%mov
    case(16) ; a=node(i)%dmo
    case(17) ; a=node(i)%orix
    case(18) ; a=node(i)%oriy
    case(19) ; a=node(i)%oriz
    case(20) ; a=node(i)%ecm
    case(21) ; a=node(i)%cod
    case(22) ; a=node(i)%grd
    case(23) ; a=node(i)%pld
    case(24) ; a=node(i)%vod
    case(25) ; a=node(i)%dif
    case(26) ; a=node(i)%kfi
    case(27) ; a=node(i)%pla
    case(28) ; a=node(i)%kvol
    case(29) ; a=node(i)%tipus 
    case(30) ; a=node(i)%icel
    case(31) ; a=node(i)%altre
    case(32) ; a=node(i)%marge
    case(33) ; a=node(i)%talone
    case(34) ; a=node(i)%fix
    case(35) ; a=node(i)%ndiv
    case(36) ; a=discentch(i)  !!>> HC 12-2-2024
    !case(36) ; a=gex(i,1)
    !case(37) ; a=gex(i,2)
    !case(38) ; a=gex(i,chogen)
    case(37) ; if(node(i)%tipus<4)then; a=cels(node(i)%icel)%fase ; else; a=0d0 ; end if !>>Miquel18-9-14
    case(38) ; a=disx(i)
    case(39) ; a=disy(i)
    case(40) ; a=disz(i)
    case(41) ; a=distot(i)
    case(42) ; a=real(boxes(nint(node(i)%x*urv),nint(node(i)%y*urv),nint(node(i)%z*urv)))	     
    case(43) ; a=i 
    case(44) ; a=nodeo(i)%icel                          !!>> HC 15-12-2021
    end select
    if(colorselection.gt.44)then  !>>>Miguel30-10-14    !!>> HC 15-12-2021
       a=gex(i,chogen)            !>>>Miguel30-10-14
    end if                        !>>>Miguel30-10-14
    !MOD!!!!!!!!!!!!!!!!!!!!!!11
    !if(node(i)%tipus==3)then; a=node(i)%icel-19 ;end if
    !!!!!!!!!!!!!!!!!!!!!!!
    !hhh=a    !>>>Miguel20-10-14
    select case(conf_rainbow)
    case(1); call get_rainbow(a,mine,maxe,color)
    case(2); call get_rainbow3(a,mine,maxe,color)
    case(3); call get_rainbow4(a,mine,maxe,color)
    case(4); call get_rainbow4blue(a,mine,maxe,color)
    case default; call get_rainbow(a,mine,maxe,color)
    end select
    !MOD!!!!!!!!!!!!!!!!!!!!!!11
    !if(node(i)%tipus==3) call get_rainbow(a,1.0_gldouble,7.0_gldouble,color)

    !if(node(i)%tipus<3)then
    !  if(gex(i,3)>epsilod)then
    !    color=0 ; color(1)=1d0 ; color(2)=1d0
    !  else
    !    color=0 ; color(1)=0.25d0 ; color(3)=0.75d0
    !  end if
    !elseif(node(i)%tipus==4)then
    !  if(gex(i,4)>epsilod)then
    !    color=0 ; color(2)=0.5d0 ; color(3)=1.0d0
    !  else
    !    color=0 ; color(1)=1d0 ; color(2)=0.5d0
    !  end if
    !end if
  

    !!!!!!!!!!!!!!!!!!!!!!!

    if(colorselection==32.or.colorselection==34) call get_rainbow(a,mine,maxe,color)

    memcolor(i,:)=color(:)  !>>>Miquel10-12-13


  !ARROW LENGTH VALUE
!print*,"arrowselection",arrowselection
  select case (arrowselection)  !>>Miquel30-10-14
    case(1) ; a=node(i)%x
    case(2) ; a=node(i)%y
    case(3) ; a=node(i)%z
    case(4) ; a=node(i)%e
    case(5) ; a=node(i)%eqd
    case(6) ; a=node(i)%add
    case(7) ; a=node(i)%you
    case(8) ; a=node(i)%adh
    case(9) ; a=node(i)%rep
    case(10) ; a=node(i)%rec
    case(11) ; a=node(i)%erp
    case(12) ; a=node(i)%est
    case(13) ; a=node(i)%eqs
    case(14) ; a=node(i)%hoo
    case(15) ; a=node(i)%mov
    case(16) ; a=node(i)%dmo
    case(17) ; a=node(i)%orix
    case(18) ; a=node(i)%oriy
    case(19) ; a=node(i)%oriz
    case(20) ; a=node(i)%ecm
    case(21) ; a=node(i)%cod
    case(22) ; a=node(i)%grd
    case(23) ; a=node(i)%pld
    case(24) ; a=node(i)%vod
    case(25) ; a=node(i)%dif
    case(26) ; a=node(i)%kfi
    case(27) ; a=node(i)%pla
    case(28) ; a=node(i)%kvol
    case(29) ; a=node(i)%tipus 
    case(30) ; a=node(i)%icel
    case(31) ; a=node(i)%altre
    case(32) ; a=node(i)%marge
    case(33) ; a=node(i)%talone
    case(34) ; a=node(i)%fix
    case(35) ; a=node(i)%ndiv  !!>> HC 16-3-2022
    case(36) ; a=discentch(i)  !!>> HC 12-2-2024
    !case(36) ; a=gex(i,1)
    !case(37) ; a=gex(i,2)
    !case(38) ; a=gex(i,chogen)
    case(37) ; if(node(i)%tipus<4)then; a=cels(node(i)%icel)%fase ; else; a=0d0 ; end if !>>Miquel18-9-14
    case(38) ; a=disx(i)
    case(39) ; a=disy(i)
    case(40) ; a=disz(i)
    case(41) ; a=distot(i)
    case(42) ; a=real(boxes(nint(node(i)%x*urv),nint(node(i)%y*urv),nint(node(i)%z*urv)))	     
    case(43) ; a=i 
  end select
  if(arrowselection.gt.43)then  !>>>Miguel30-10-14
     a=gex(i,achogen)            !>>>Miguel30-10-14
  end if                        !>>>Miguel30-10-14

  memarrow(i)=a



  !TRANSPARENT SPHERE RADIUS VALUE
  select case (sphereselection)  !>>Miquel30-10-14
    case(1) ; a=node(i)%x
    case(2) ; a=node(i)%y
    case(3) ; a=node(i)%z
    case(4) ; a=node(i)%e
    case(5) ; a=node(i)%eqd
    case(6) ; a=node(i)%add
    case(7) ; a=node(i)%you
    case(8) ; a=node(i)%adh
    case(9) ; a=node(i)%rep
    case(10) ; a=node(i)%rec
    case(11) ; a=node(i)%erp
    case(12) ; a=node(i)%est
    case(13) ; a=node(i)%eqs
    case(14) ; a=node(i)%hoo
    case(15) ; a=node(i)%mov
    case(16) ; a=node(i)%dmo
    case(17) ; a=node(i)%orix
    case(18) ; a=node(i)%oriy
    case(19) ; a=node(i)%oriz
    case(20) ; a=node(i)%ecm
    case(21) ; a=node(i)%cod
    case(22) ; a=node(i)%grd
    case(23) ; a=node(i)%pld
    case(24) ; a=node(i)%vod
    case(25) ; a=node(i)%dif
    case(26) ; a=node(i)%kfi
    case(27) ; a=node(i)%pla
    case(28) ; a=node(i)%kvol
    case(29) ; a=node(i)%tipus 
    case(30) ; a=node(i)%icel
    case(31) ; a=node(i)%altre
    case(32) ; a=node(i)%marge
    case(33) ; a=node(i)%talone
    case(34) ; a=node(i)%fix
    case(35) ; a=node(i)%ndiv  !!>> HC 16-3-2022
    case(36) ; a=discentch(i)  !!>> HC 12-2-2024
    !case(36) ; a=gex(i,1)
    !case(37) ; a=gex(i,2)
    !case(38) ; a=gex(i,chogen)
    case(37) ; if(node(i)%tipus<4)then; a=cels(node(i)%icel)%fase ; else; a=0d0 ; end if !>>Miquel18-9-14
    case(38) ; a=disx(i)
    case(39) ; a=disy(i)
    case(40) ; a=disz(i)
    case(41) ; a=distot(i)
    case(42) ; a=real(boxes(nint(node(i)%x*urv),nint(node(i)%y*urv),nint(node(i)%z*urv)))	     
    case(43) ; a=i 
  end select
  if(sphereselection.gt.43)then  !>>>Miguel30-10-14
     a=gex(i,schogen)            !>>>Miguel30-10-14
  end if                        !>>>Miguel30-10-14

  memsphere(i)=a

    !MARK SELECTED NODES
    if(flag(7)==1)then
      r=0.05*esca
      if(flag(3)==1)then !mark the last selected node with a different color !>>Miquel31-10-14
        do j=1,nki
          ki=oopp(j)
          if(i==ki)then
            color=1
            if(i==nodeindex)then
              color=1 ; color(3)=0
            end if
            r=node(i)%add*0.5d0*esca
          end if
        enddo
      end if
      if(flag(4)==1)then !mark the last selected cell with a different color 
        do j=1,nki
          ki=oopp(j)
          if(node(i)%icel==ki)then
            color=1
            if(node(i)%icel==cellid)then
              color=1 ; color(3)=0
            end if
            r=node(i)%add*0.5d0*esca
          end if
        enddo
      end if
    end if

    !DETERMINING SIZE OF NODE
    if(flag(5)==1)then
      r=node(i)%eqd*esca				!>>>>>>>>>>>>>>>>>>>MIQUEL 4-3-13
    else
      if (flag(6)==1)then
        r=node(i)%add*esca
      end if
    end if
    if (flag(8)==1) r=0
    if (flag(11)==1) r=node(i)%add*esca*1d-1

  ! cursor for the editor
  if (flag(40)==1) then
    a=cursx ; b=cursy ; c=cursz
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
    call glTranslated(a*esca,b*esca,c*esca)
    d=0.3d-1
    call glutsolidsphere(d, 4, 4)
    call glTranslated(-a*esca,-b*esca,-c*esca)
  end if


    !DRAWING THE NODES
    if(flag(36)==0.and.node(i)%fix==1) cycle  !>>>>Miquel15-1-14
    if ((flag(30)==1.and.node(i)%tipus<3).or.(flag(31)==1.and.node(i)%tipus==3).or.(flag(32)==1.and.node(i)%tipus==4)) then
      if ((flag(9)==1.and.node(i)%tipus==1).or.(flag(10)==1.and.node(i)%tipus==2) .or.(node(i)%tipus>=3)) then      
        a=node(i)%x ; b=node(i)%y ;  c=node(i)%z      
        if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then 
          if (flag(22)==1.and.node(i)%tipus<3)then
            call cylinder(i,10,color)   !>>>>>>>>>> Miquel 11-6-13
          else
            if (flag(22)==1)then                 !!>> HC 15-12-2021 This precludes the spheres from turning black when we have cylinders
               normal=0.010d0                    !!>> HC 15-12-2021  I put the same normal to all spheres
               call glnormal3fv(normal)          !!>> HC 15-12-2021  The right normals should be calculated depending on the viewpoint
            endif                                !!>> HC 15-12-2021
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
            
            if (flag(33)==1)then                                                     !!>> HC 8-4-2022  ELONGATED CELLS AS ELLIPSES
               eloch=0                                                               !!>> HC 8-4-2022 
               if (npag(nparam_per_node+17) > 0) then                                !!>> HC 8-4-2022               
                  myfactor=1.0d0                                                     !!>> HC 8-4-2022  Compute the factor by which ADD and EQD are modified
                  do jch=1, npag(nparam_per_node+17)                                 !!>> HC 8-4-2022  Genes affecting elongation
                     kk=whonpag(nparam_per_node+17,jch)                              !!>> HC 8-4-2022 
                     if (gex(i,kk)>0.0d0) then                                       !!>> HC 8-4-2022 
                        myfactor=myfactor+gex(i,kk)*gen(kk)%e(nparam_per_node+17)    !!>> HC 8-4-2022  expresion of those in node i times the strength of activation
                     end if                                                          !!>> HC 8-4-2022 
                  enddo                                                              !!>> HC 8-4-2022 
                  if(myfactor > maxelong)myfactor=maxelong                           !!>> HC 8-4-2022  there is a maximum elongation value (free parameter)
                  umyfactor=1.0d0/myfactor                                           !!>> HC 8-4-2022 
                  cpolx=cels(node(i)%icel)%polx                                      !!>> HC 8-4-2022  Polarization vector of the cell
                  cpoly=cels(node(i)%icel)%poly                                      !!>> HC 8-4-2022 
                  cpolz=cels(node(i)%icel)%polz                                      !!>> HC 8-4-2022 
                  if (cpolx.ne.0.0d0.or.cpoly.ne.0.0d0.or.cpolz.ne.0.0d0)then        !!>> HC 8-4-2022  If there is a polarization vector (any of these .ne. 0)
                     mpol=sqrt(cpolx**2+cpoly**2+cpolz**2)                           !!>> HC 8-4-2022  module of the polarization vector
                     anglech=acos(cpolx/mpol)                                        !!>> HC 8-4-2022  angle between the polarization vector and vector (1,0,0)
                     anglech=anglech*180.0d0/pi                                      !!>> HC 8-4-2022  angle in degrees
                     rvecy=-cpolz                                                    !!>> HC 8-4-2022  Cross product between the polarization vector and vector (1,0,0)
                     rvecz=cpoly                                                     !!>> HC 8-4-2022  (x coordinate will always be =0.0d0)
                     eloch=1                                                         !!>> HC 8-4-2022  This node is elongated
                  endif                                                              !!>> HC 8-4-2022 
               endif                                                                 !!>> HC 8-4-2022 
               call glTranslated(a*esca,b*esca,c*esca)                               !!>> HC 8-4-2022  Move to the scaled position of the node
               if (eloch==1)then                                                     !!>> HC 8-4-2022 
                  call glRotatef(anglech,0.0000, rvecy, rvecz)                       !!>> HC 8-4-2022  rotate the node so that the x axis has the same direction of the polarization vector
                  call glScaled(myfactor, umyfactor, umyfactor)                      !!>> HC 8-4-2022  elongate in the new x axis (=direction of the polarization vector)
               endif                                                                 !!>> HC 8-4-2022 
               call glutsolidsphere(r, 10, 10)                                       !!>> HC 8-4-2022  make the elongated sphere
               if (eloch==1)then                                                     !!>> HC 8-4-2022 
                  anglech=360-anglech                                                !!>> HC 8-4-2022  new rotation angle to go back to the initial orientation
                  call glScaled(umyfactor, myfactor, myfactor)                       !!>> HC 8-4-2022  rescale to make round spheres again
                  call glRotatef(anglech,0.0000, rvecy, rvecz)                       !!>> HC 8-4-2022  recover old orientation
               endif                                                                 !!>> HC 8-4-2022 
                call glTranslated(-a*esca,-b*esca,-c*esca)                           !!>> HC 8-4-2022  move back to the origin
            else                                                                     !!>> HC 8-4-2022 
               call glTranslated(a*esca,b*esca,c*esca)                               !!>> HC 8-4-2022  NORMAL SPHERICAL NODES
               call glutsolidsphere(r, 10, 10)                                       !!>> HC 8-4-2022 
               call glTranslated(-a*esca,-b*esca,-c*esca)                            !!>> HC 8-4-2022 
            endif                                                                    !!>> HC 8-4-2022 
                

          end if
          !!!!!!!!!!!!!!!!!!!
          if((flag(28)==0).and.(flag(29)==1))then ; flag(29)=2 ; endif   ! initial switcher
          if((flag(29)==0).and.(flag(28)==1))then ; flag(28)=2 ; endif   ! initial switcher
          if((flag(28)==2).and.(flag(29)==1))then ; flag(28)=0 ; endif   ! second  switcher
          if((flag(28)==1).and.(flag(29)==2))then ; flag(29)=0 ; endif   ! second  switcher
               
          !DRAWING THE ARROWS
          if(arrowselection==0)then
            !amax_scale=1d0  !this resets the scale in case it has changed    !  ARROWS                 !!>>>Miguel20-10-14 made flag        
          else                                                                                                          
            if(node(i)%tipus.le.2)then        
              aa=node(node(i)%altre)%x ; bb=node(node(i)%altre)%y ; cc=node(node(i)%altre)%z                 
              aa=a-aa ; bb=b-bb ; cc=c-cc      
              d=sqrt(aa**2+bb**2+cc**2) ; aa=aa/d ; bb=bb/d ; cc=cc/d          
              !if(flag(33).eq.1)then         ! change scaling
              !  write(*,*)'Enter new value for the arrow length corresponding to the maximum value of the property being printed'
              !  write(*,*)'Actual maximum value=',amaxe,'actual length=',amax_scale,'-1 to reset default value (1)'
              !  read(*,*)amax_scale ; flag(33)=0
              !end if  
              !if(amax_scale.eq.-1)then;amax_scale=1d0;endif
              !dd=maxe
              r=amax_scale*(memarrow(i)-amine)/(amaxe-amine)
              color(1:4)=memcolor(i,1:4)
              call arrows(a,b,c,a+(aa*r),b+(bb*r),c+(cc*r),color)        
              call glEnd             
            end if        
          end if
          !DRAWING THE SEMITRANSPARENT BALLS
          if(sphereselection==0)then
            !smax_scale=nodeo(1)%add  !this resets the scale in case it has changed      !!>>>Miguel20-10-14 made flag        
          else                                       !  SEMITRANSPARENT BALLS   !!>>>Miguel20-10-14 made flag  
            if(node(i)%tipus.le.2)then                
              !if(flag(33).eq.1)then              ! change scaling           
              !  write(*,*)'Enter new value for the sphere radius corresponding to the maximum value of the property being printed'
              !  write(*,*)'Actual maximum value=,',smaxe,'Actual radius=',smax_scale,'-1 to reset default value (nodeo(1)%add)'
              !  read(*,*)smax_scale ; flag(33)=0
              !end if        
              !dd=maxe
              r=smax_scale*(memsphere(i)-smine)/(smaxe-smine)
              r=esca*r
              !if(fff.eq.-1)then;fff=nodeo(1)%add;endif       
              color(1:3)=memcolor(i,1:3)
              color(4)=0.25   
              call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)   
              call glTranslated(a*esca,b*esca,c*esca)                                      
              call glEnable(GL_BLEND)    
              call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)             
              call glutsolidsphere(r, 10, 30)          
              call gldisable(GL_BLEND)                    
              call glTranslated(-a*esca,-b*esca,-c*esca) 
              call glEnd                                  
              color(:)=memcolor(i,:)
            end if    							      
          end if
          !!!!!!!!!!!!!!!!!!!!!!          
        endif
       end if
    end if
  end do

  !we draw the intercellular connexions with lines (nodes further than node()%add won't be connected)
  call glBegin(gl_lines)
   !>>>>>>>>>>>>BEGIN Miquel 6-5-13 modified Is 15-5-13
  if(flag(12)==1.or.flag(13)==1)then   
    do i=1,nd
      if(node(i)%tipus==1.and.flag(9)==0) cycle
      if(node(i)%tipus==2.and.flag(10)==0) cycle
      if(node(i)%tipus<3.and.flag(30)==0) cycle
      if(node(i)%tipus==3.and.flag(31)==0) cycle
      if(node(i)%tipus>3.and.flag(32)==0) cycle
      if(flag(36)==0.and.node(i)%fix==1) cycle  !>>>>Miquel15-1-14
      if((flag(30)==1.and.node(i)%tipus<3).or.(flag(31)==1.and.node(i)%tipus==3).or.(flag(32)==1.and.node(i)%tipus==4))then
        a=node(i)%x ; b=node(i)%y ; c=node(i)%z
        if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then
          do j=1,nneigh(i)
            ic=neigh(i,j)
            if((flag(12)==1.and.node(i)%icel/=node(ic)%icel).or.(flag(13)==1 .and. node(i)%icel==node(ic)%icel))then
              aa=node(ic)%x ; bb=node(ic)%y ; cc=node(ic)%z
              if((flag(30)==1.and.node(ic)%tipus<3).or.(flag(31)==1.and.node(ic)%tipus==3).or.&
                                                         &(flag(32)==1.and.node(ic)%tipus==4))then
                if (aa>=mix.and.aa<=mx.and.bb>=miy.and.bb<=my.and.cc>=miz.and.cc<=mz) then
                  d=dneigh(i,j)
                  if(d<node(i)%add+node(ic)%add)then     !>>>>>> Is 9-5-13
                    if (flag(12)==1) then
                      if(node(i)%tipus==1)then;color=0;color(1)=1
                      else;color=0;color(2)=1;end if
                    else
                      if(node(i)%tipus==node(ic)%tipus)then;color=0;color(2)=0.75
                      else;color=0;color(1)=0.75;end if
                    end if
                    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
                    call glvertex3d(esca*a,esca*b,esca*c)
                    call glvertex3d(esca*aa,esca*bb,esca*cc)
                  end if
                end if
              end if
            end if
          end do
        end if
      end if
    end do
  end if
  !>>>>>>>>>>>>END Miquel 6-5-13
  call glEnd


  
  if(flag(37)==1)then !drawing the cell contours !>>Miquel29-8-14
    print*, "SORRY THIS FLAG IS DISABLED IN THIS VERSION DUE TO RAM MEMORY ISSUES" !!>> HC 18-8-21
  end if

  if (flag(14)==1) then  !2D plot x plane
    mine=minval(theslice) ; maxe=maxval(theslice)
    do i=0,nq
      do j=0,nq
        call glBegin(gl_points)
         d=1d0/nq
         a=mixi;b=miy+i*(my-miy)*d;c=miz+j*(mz-miz)*d
         coo=theslice(i,j)
         if (flag(17)==1) a=(mixi+coo)*esca*1d-6 !this is to make an actual plot (3D)
         call get_rainbow(coo,mine,maxe,color)
         call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         call glvertex3d(a*esca,b*esca,c*esca)
        call glEnd
      end do
    end do
  end if

  if (flag(15)==1) then  !2D plot y plane
    mine=minval(theslice) ; maxe=maxval(theslice)
    do i=0,nq
      do j=0,nq
        call glBegin(gl_points)
         d=1d0/nq
         a=mix+i*(mx-mix)*d;b=miyi;c=miz+j*(mz-miz)*d
         coo=theslice(i,j)
         if (flag(17)==1) b=(miyi+coo)*esca*1d-6 !this is to make an actual plot (3D) 
         call get_rainbow(coo,mine,maxe,color)
         call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         call glvertex3d(a*esca,b*esca,c*esca)
        call glEnd
      end do
    end do
  end if

  if (flag(16)==1) then  !2D plot z plane
    mine=minval(theslice) ; maxe=maxval(theslice)
    do i=0,nq
      do j=0,nq
        call glBegin(gl_points)
         d=1d0/nq
         a=mix+i*(mx-mix)*d;b=miy+j*(my-miy)*d;c=mizi
         coo=theslice(i,j)
         if (flag(17)==1) c=(mizi+coo)*esca*1d-4 !this is to make an actual plot (3D)
         call get_rainbow(coo,mine,maxe,color)
         call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         call glvertex3d(a*esca,b*esca,c*esca)
        call glEnd
      end do
    end do
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! miguel 30-5-13
  if (flag(21)==1) then                                                         ! miguel 30-5-13
    do i=1,ncels                                                              ! miguel 30-5-13 
      a=cels(i)%cex ; b=cels(i)%cey ; c=cels(i)%cez ; r=1d-2                ! miguel 30-5-13 
      if (a>=mix.and.a<=mx.and.b>=miy.and.b<=my.and.c>=miz.and.c<=mz) then  ! miguel 30-5-13   
        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,white)   ! miguel 30-5-13        
        call glTranslated(a*esca,b*esca,c*esca)                             ! miguel 30-5-13
        call glutsolidsphere(r, 10, 10)                                     ! miguel 30-5-13   
        call glTranslated(-a*esca,-b*esca,-c*esca)                          ! miguel 30-5-13
      end if								    ! miguel 30-5-13
    end do     								    ! miguel 30-5-13
  end if								    ! miguel 30-5-13
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
call glEndList
call glutPostRedisplay

end subroutine draw_func
!*******************************************
!*******************************************SUBROUTINE********************************************************
subroutine arrows(a,b,c,aa,bb,cc,color)                          !!>>>Miguel20-10-14 made subroutine
real*8 :: a,b,c,d,aa,bb,cc,dd,aaa,bbb,ccc,ddd,dddd,diametre,kx,ky,kz
real(GLFLOAT)::color(4) 
aaa=0d0 ; bbb=0d0 ; ccc=0d0 ; ddd=0d0 ; dddd=0d0
diametre=3d-2
 kx=aa ;  ky=bb ;  kz=cc     ! so as the arrowhead does not add extra lenght
call arrow(a,b,c,kx,ky,kz,diametre,diametre,color)   !cylinder (main body of the arrow) 
aaa=aa-a ; bbb=bb-b ; ccc=cc-c
ddd=sqrt(aaa**2+bbb**2+ccc**2)
aaa=aaa/ddd ; bbb=bbb/ddd ; ccc=ccc/ddd ! unit vector
dddd=0.3 ! "correction factor"
call arrow(kx,ky,kz,kx+(aaa*dddd),ky+(bbb*dddd),kz+(ccc*dddd),2d0*diametre,1d-3,color)   !arrowhead
end subroutine arrows

subroutine arrow(ix,iy,iz,vx,vy,vz,r,riv,color)  !!>>>Miguel20-10-14 made subroutine
integer::nod,iv,res,i,j,k,ii,jj,kk
real*8::r,riv,cx,cy,cz,ix,iy,iz,vx,vy,vz
real*8::rad,ax,ay,az,ux,uy,uz,thet
real*8::bx,by,bz,ox,oy,oz,cost,sint,ucost
real*8::a,b,c,d,aa,bb,cc,dd,aaa,bbb,ccc,ddd,aaaa,bbbb,cccc,dddd
real(GLFLOAT)::color(4),calor(4)
res=6 ! sides of the cylinder

  cx=vx-ix ; cy=vy-iy ; cz=vz-iz
  d=sqrt(cx**2+cy**2+cz**2)
  ux=cx/d ; uy=cy/d ; uz=cz/d

  !vector ortogonal al spring, fem el cross product amb el vector arbitrari (1,0,0)
  ax=0 ; ay=-cz ; az=cy !vector inicial, a partir d'aqui rotem
  d=sqrt(ay**2+az**2)
  ay=ay/d ; az=az/d

  rad=2*pi/real(res)

  ox=ax ; oy=ay ; oz=az
  do i=1,res-1
    thet=i*rad
    cost=cos(thet); ucost=1-cost ; sint=sin(thet)
    bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
    by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
    bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az

    call glBegin(gl_triangles) !the upper and lower face triangles
!      color=0 ; color(1)=1
      aa=ix+ox*r ; bb=iy+oy*r ; cc=iz+oz*r
      aaa=ix+bx*r ; bbb=iy+by*r ; ccc=iz+bz*r
      normal=normcrossprod((/ix,aa,aaa/),(/iy,bb,bbb/),(/iz,cc,ccc/)) 
      normal=abs(normal)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
      call glvertex3d(ix*esca,iy*esca,iz*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)

!      color=0 ; color(1)=1
      aa=vx+ox*riv ; bb=vy+oy*riv ; cc=vz+oz*riv
      aaa=vx+bx*riv ; bbb=vy+by*riv ; ccc=vz+bz*riv
      normal=normcrossprod((/vx,aa,aaa/),(/vy,bb,bbb/),(/vz,cc,ccc/)) 
      normal=abs(normal)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
      call glvertex3d(vx*esca,vy*esca,vz*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
    call glEnd

    call glBegin(gl_quads)
      a=ix+ox*r ; b=iy+oy*r ; c=iz+oz*r
      aa=ix+bx*r ; bb=iy+by*r ; cc=iz+bz*r
      aaa=vx+bx*riv ; bbb=vy+by*riv ; ccc=vz+bz*riv
      normal=normcrossprod((/a,aa,aaa/),(/b,bb,bbb/),(/c,cc,ccc/)) 
      normal=abs(normal)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)

      call glvertex3d(a*esca,b*esca,c*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
      call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
    call glEnd

    call glBegin(gl_line_loop)
      calor=0 !; color(1)=1
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!      call glvertex3d(ix*esca,iy*esca,iz*esca)
      call glvertex3d((ix+ox*r)*esca,(iy+oy*r)*esca,(iz+oz*r)*esca)
      call glvertex3d((ix+bx*r)*esca,(iy+by*r)*esca,(iz+bz*r)*esca)
    call glEnd
    call glBegin(gl_line_loop)
      calor=0 !; color(1)=1
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!      call glvertex3d(vx*esca,vy*esca,vz*esca)
      call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
      call glvertex3d((vx+bx*riv)*esca,(vy+by*riv)*esca,(vz+bz*riv)*esca)
    call glEnd
    ox=bx ; oy=by ; oz=bz
  end do

  call glBegin(gl_triangles) !the upper and lower face triangles
!    color=0 ; color(1)=1
    aa=ix+ox*r ; bb=iy+oy*r ; cc=iz+oz*r
    aaa=ix+ax*r ; bbb=iy+ay*r ; ccc=iz+az*r
    normal=normcrossprod((/ix,aa,aaa/),(/iy,bb,bbb/),(/iz,cc,ccc/)) 
    normal=abs(normal)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(ix*esca,iy*esca,iz*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)

!    color=0 ; color(1)=1
    aa=vx+ox*riv ; bb=vy+oy*riv ; cc=vz+oz*riv
    aaa=vx+ax*riv ; bbb=vy+ay*riv ; ccc=vz+az*riv
    normal=normcrossprod((/vx,aa,aaa/),(/vy,bb,bbb/),(/vz,cc,ccc/)) 
    normal=abs(normal)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(vx*esca,vy*esca,vz*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
  call glEnd

  call glBegin(gl_quads)
    a=ix+ox*r ; b=iy+oy*r ; c=iz+oz*r
    aa=ix+ax*r ; bb=iy+ay*r ; cc=iz+az*r
    aaa=vx+ax*riv ; bbb=vy+ay*riv ; ccc=vz+az*riv
    normal=normcrossprod((/a,aa,aaa/),(/b,bb,bbb/),(/c,cc,ccc/)) 
    normal=abs(normal)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(a*esca,b*esca,c*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
    call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
  call glEnd

  call glBegin(gl_line_loop)
    calor=0 !; color(1)=1
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!    call glvertex3d(ix*esca,iy*esca,iz*esca)
    call glvertex3d((ix+ox*r)*esca,(iy+oy*r)*esca,(iz+oz*r)*esca)
    call glvertex3d((ix+ax*r)*esca,(iy+ay*r)*esca,(iz+az*r)*esca)
  call glEnd
  call glBegin(gl_line_loop)
    calor=0 !; color(1)=1
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!    call glvertex3d(vx*esca,vy*esca,vz*esca)
    call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
    call glvertex3d((vx+ax*riv)*esca,(vy+ay*riv)*esca,(vz+az*riv)*esca)
  call glEnd

end subroutine arrow

!*******************************************SUBROUTINE********************************************************

subroutine get_rainbow(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   if (f/=0.0) then
     f=0.0_glfloat
   else
     f = 0.5_glfloat
   end if
endif


!if (f < .07) then
!   c(1) = 0.6_glfloat
!   c(2) = 0.6_glfloat 
!   c(3) = 0.6_glfloat
!   c(4) = 0.8_glfloat
if (f < 0.5) then
   c(1) = f*2.0
   c(2) = 0
   c(3) = 1-(f*2.0)
   c(4) = 0.5_glfloat
elseif (f < 1.0) then
   c(1) = 1.0_glfloat 
   c(2) = 2*f-0.5
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
else
   c(1) = 1.0_glfloat 
   c(2) = 1.0_glfloat 
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
endif

end subroutine 




!**************************************************************************************************
!THIS IS A YELLOW-RED RAINBOW
subroutine get_rainbow3(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

   c(1) = 1.0_glfloat 
   c(2) = f!1.0_glfloat
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat

end subroutine 

!*******************************************SUBROUTINE********************************************************

!THIS IS A GREEN RAINBOW
subroutine get_rainbow4(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

   c(1) = 0.0_glfloat !1.0_glfloat 
   c(2) = f!1.0_glfloat
   c(3) = 1.0_glfloat
   c(4) = 1.0_glfloat

end subroutine 


!THIS IS A BLUE RAINBOW
subroutine get_rainbow4blue(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

   c(1) = 0.0_glfloat !1.0_glfloat 
   c(2) = 1.0_glfloat
   c(3) = f!1.0_glfloat
   c(4) = 1.0_glfloat

end subroutine 


!*******************************************SUBROUTINE********************************************************

!THIS IS A TOTAL RAINBOW
subroutine get_rainbow5(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f,fu,fd,ft

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

fu=f
if (fu>0.333) f=0.333

fd=f-0.333
if (fd>0.666) f=0.666

ft=f-0.666
if (ft>1) f=1

c(1) = fu !1.0_glfloat 
c(2) = fd
c(3) = ft
 c(4) = 1.0_glfloat

end subroutine 

!*******************************************SUBROUTINE********************************************************

!THIS IS A FOR TIPUS RAINBOW
subroutine get_rainbow6(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f,fu,fd,ft

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

a=16581375

if (f<=0.34) then 
  fu=0.1  ; fd=(f+0.3)*1.1  ; ft=0!f*2
end if

if (f>0.44.and.f<0.75) then 
  fu=0.6; fd=0.1 ; ft=0.1
end if

if (f>0.75) then 
  fu=f ; fd=f ; ft=1-f
end if

c(1) = fu
c(2) = fd
c(3) = ft
c(4) = 1.0_glfloat

end subroutine 


!*******************************************SUBROUTINE********************************************************

!THIS IS A FOR ICEL RAINBOW
subroutine get_rainbow7(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f,fu,fd,ft

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

a=16581375

fu=0 ; fd=1-f*0.8 ; ft=0

c(1) = fu
c(2) = fd
c(3) = ft
c(4) = 1.0_glfloat

end subroutine 


!*******************************************SUBROUTINE********************************************************


function normcrossprod(x,y,z)
real(glfloat), dimension(3) :: normcrossprod
real(gldouble), dimension(3), intent(in) :: x,y,z
real(glfloat) :: t1(3),t2(3),norm
t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)
normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)
norm = sqrt(dot_product(normcrossprod,normcrossprod))
if (norm /= 0._glfloat) normcrossprod = normcrossprod/norm
end function normcrossprod
!*******************************************SUBROUTINE********************************************************

subroutine cylinder(nod,res,color)  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Miquel 15-5-13
integer::nod,iv,res,i,j,k,ii,jj,kk
real*8::r,riv,cx,cy,cz,ix,iy,iz,vx,vy,vz,d
real*8::rad,ax,ay,az,ux,uy,uz,thet
real*8::bx,by,bz,ox,oy,oz,cost,sint,ucost
real(GLFLOAT)::color(4),calor(4)

  iv=node(nod)%altre

  if(flag(5)==1)then
    r=node(nod)%eqd;riv=node(iv)%eqd
  else
    r=node(nod)%add;riv=node(iv)%add
  end if

  ix=node(nod)%x ; iy=node(nod)%y ; iz=node(nod)%z
  vx=node(iv)%x ; vy=node(iv)%y ; vz=node(iv)%z
  cx=vx-ix ; cy=vy-iy ; cz=vz-iz
  d=sqrt(cx**2+cy**2+cz**2)
  ux=cx/d ; uy=cy/d ; uz=cz/d

  !this will be half a cylinder, so we cut it in a half
  vx=vx-0.5d0*cx ; vy=vy-0.5d0*cy ; vz=vz-0.5d0*cz
  riv=riv+0.5d0*(r-riv)

  !vector ortogonal al spring, fem el cross product amb el vector arbitrari (1,0,0)
  ax=0 ; ay=-cz ; az=cy !vector inicial, a partir d'aqui rotem
  d=sqrt(ay**2+az**2)
  ay=ay/d ; az=az/d
!  ay=ay*r ; az=az*r
!  print*,"U",ux,uy,uz,"A",ax,ay,az


!  if(flag(6)==1)then  !the form of the cylinder encloses the volume above and below it corresponding to node()%add
    ix=ix-ux*r ; iy=iy-uy*r ; iz=iz-uz*r
!    vx=vx+ux*riv ; vy=vy+uy*riv ; vz=vz+uz*riv
!  end if

  rad=2*pi/real(res)

  ox=ax ; oy=ay ; oz=az
  do i=1,res-1
    thet=i*rad
    cost=cos(thet); ucost=1-cost ; sint=sin(thet)
    bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
    by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
    bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az

    call glBegin(gl_triangles) !the upper and lower face triangles
!      color=0 ; color(1)=1
      aa=ix+ox*r ; bb=iy+oy*r ; cc=iz+oz*r
      aaa=ix+bx*r ; bbb=iy+by*r ; ccc=iz+bz*r
      normal=normcrossprod((/ix,aa,aaa/),(/iy,bb,bbb/),(/iz,cc,ccc/)) 
      normal=abs(normal)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
      call glvertex3d(ix*esca,iy*esca,iz*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)

!      color=0 ; color(1)=1
      aa=vx+ox*riv ; bb=vy+oy*riv ; cc=vz+oz*riv
      aaa=vx+bx*riv ; bbb=vy+by*riv ; ccc=vz+bz*riv
      normal=normcrossprod((/vx,aa,aaa/),(/vy,bb,bbb/),(/vz,cc,ccc/)) 
      normal=abs(normal)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
      call glvertex3d(vx*esca,vy*esca,vz*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
    call glEnd

    call glBegin(gl_quads)
      a=ix+ox*r ; b=iy+oy*r ; c=iz+oz*r
      aa=ix+bx*r ; bb=iy+by*r ; cc=iz+bz*r
      aaa=vx+bx*riv ; bbb=vy+by*riv ; ccc=vz+bz*riv
      normal=normcrossprod((/a,aa,aaa/),(/b,bb,bbb/),(/c,cc,ccc/)) 
      !normal=abs(normal)     !>>>Miguel29-10-14 (comment)
      call glnormal3fv(normal)
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
      call glvertex3d(a*esca,b*esca,c*esca)
      call glvertex3d(aa*esca,bb*esca,cc*esca)
      call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
      call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
    call glEnd

    call glBegin(gl_line_loop)
      calor=0 !; color(1)=1
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!      call glvertex3d(ix*esca,iy*esca,iz*esca)
      call glvertex3d((ix+ox*r)*esca,(iy+oy*r)*esca,(iz+oz*r)*esca)
      call glvertex3d((ix+bx*r)*esca,(iy+by*r)*esca,(iz+bz*r)*esca)
    call glEnd
    call glBegin(gl_line_loop)
      calor=0 !; color(1)=1
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!      call glvertex3d(vx*esca,vy*esca,vz*esca)
      call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
      call glvertex3d((vx+bx*riv)*esca,(vy+by*riv)*esca,(vz+bz*riv)*esca)
    call glEnd
!    call glBegin(gl_lines)
!      color=0 !; color(1)=1
!      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
!      call glvertex3d((ix+bx)*esca,(iy+by)*esca,(iz+bz)*esca)
!      call glvertex3d((vx+bx)*esca,(vy+by)*esca,(vz+bz)*esca)
!    call glEnd
    ox=bx ; oy=by ; oz=bz
  end do

  call glBegin(gl_triangles) !the upper and lower face triangles
!    color=0 ; color(1)=1
    aa=ix+ox*r ; bb=iy+oy*r ; cc=iz+oz*r
    aaa=ix+ax*r ; bbb=iy+ay*r ; ccc=iz+az*r
    normal=normcrossprod((/ix,aa,aaa/),(/iy,bb,bbb/),(/iz,cc,ccc/)) 
    normal=abs(normal)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(ix*esca,iy*esca,iz*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)

!    color=0 ; color(1)=1
    aa=vx+ox*riv ; bb=vy+oy*riv ; cc=vz+oz*riv
    aaa=vx+ax*riv ; bbb=vy+ay*riv ; ccc=vz+az*riv
    normal=normcrossprod((/vx,aa,aaa/),(/vy,bb,bbb/),(/vz,cc,ccc/)) 
    normal=abs(normal)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(vx*esca,vy*esca,vz*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
  call glEnd

  call glBegin(gl_quads)
    a=ix+ox*r ; b=iy+oy*r ; c=iz+oz*r
    aa=ix+ax*r ; bb=iy+ay*r ; cc=iz+az*r
    aaa=vx+ax*riv ; bbb=vy+ay*riv ; ccc=vz+az*riv
    normal=normcrossprod((/a,aa,aaa/),(/b,bb,bbb/),(/c,cc,ccc/)) 
   ! normal=abs(normal)       !>>>Miguel29-10-14 (comment)
    call glnormal3fv(normal)
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
    call glvertex3d(a*esca,b*esca,c*esca)
    call glvertex3d(aa*esca,bb*esca,cc*esca)
    call glvertex3d(aaa*esca,bbb*esca,ccc*esca)
    call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
  call glEnd

  call glBegin(gl_line_loop)
    calor=0 !; color(1)=1
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!    call glvertex3d(ix*esca,iy*esca,iz*esca)
    call glvertex3d((ix+ox*r)*esca,(iy+oy*r)*esca,(iz+oz*r)*esca)
    call glvertex3d((ix+ax*r)*esca,(iy+ay*r)*esca,(iz+az*r)*esca)
  call glEnd
  call glBegin(gl_line_loop)
    calor=0 !; color(1)=1
    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,calor)
!    call glvertex3d(vx*esca,vy*esca,vz*esca)
    call glvertex3d((vx+ox*riv)*esca,(vy+oy*riv)*esca,(vz+oz*riv)*esca)
    call glvertex3d((vx+ax*riv)*esca,(vy+ay*riv)*esca,(vz+az*riv)*esca)
  call glEnd
!  call glBegin(gl_lines)
!    color=0 !; color(1)=1
!    call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
!    call glvertex3d((ix+ax)*esca,(iy+ay)*esca,(iz+az)*esca)
!    call glvertex3d((vx+ax)*esca,(vy+ay)*esca,(vz+az)*esca)
!  call glEnd

end subroutine cylinder

!*******************************************SUBROUTINE********************************************************

subroutine menu_handler(selection)
integer(kind=glint), value :: selection
real*8 :: guardaene

select case (selection)

case (3)
	iqpas=1
 	call iteracio(iqpas)
   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if
   call draw_func
case (4)
	iqpas=100
 	call iteracio(iqpas)
   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if
   call draw_func

case (5)
	iqpas=1000
 	call iteracio(iqpas)

   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if
   call draw_func

case (6)
	iqpas=50000
 	call iteracio(iqpas)
   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if
   call draw_func

case (7)
	print*,"how many iterations?"
	read(*,*)iqpas
 	call iteracio(iqpas)
   if(fix_run==0)then !>>Miquel2-10-14
     dyn=1
   else
     dyn=0
   end if
   call draw_func

!case (2)
!	print*,"which node?"
!	read(*,*)ki
!        print *,"to which x (from x=",node(ki)%x,") -1 keep the same value"
!        read(*,*) a
!        if (a/=-1) node(ki)%x=a
!        print *,"to which y (from y=",node(ki)%y,") -1  keep the same value"
!        read(*,*) a
!        if (a/=-1) node(ki)%y=a
!        print *,"to which z (from z=",node(ki)%z,") -1 keep the same value"
!        read(*,*) a
!        if (a/=-1) node(ki)%z=a
!        call iniboxesll
!	call draw_func
!case(18)                      ! >>>Miguel8-10-14  made "case" 
!        nodoss=0 ; kii=0 ; ki=0
!        do while(ki.ne.-1)
!          kii=kii+1
!          print*,"which node(s)? (-1 to skip select node)"        ! >>>Miguel8-10-14  
!          read(*,*)ki
!          if(kii.gt.10)then;write(*,*)'too many nodes chosen(max=10)';exit;endif  
!          if(ki.gt.nd)then;write(*,*)'wrong input:there are only',nd,'nodes';kii=kii-1;cycle;endif 
!          if((ki.eq.0).or.(ki.lt.-1))then;write(*,*)'wrong input: please choose a positive value';kii=kii-1;cycle;endif	                   
!          nodoss(kii)=ki          
!        end do
!        if(ki.eq.-1)then;nodoss(kii)=0;endif
!       geness=0 ;kii=0 ; ki=0
!        do while(ki.ne.-1)
!          kii=kii+1
!          print*,"which genes(s)? (-1 to skip select gene)"        ! >>>Miguel8-10-14  
!          read(*,*)ki
!          if(kii.gt.10)then;write(*,*)'too many genes chosen(max=10)';exit;endif  
!          if(ki.gt.ng)then;write(*,*)'wrong input: there are only',ng,'genes';kii=kii-1;cycle;endif 	                   
!          if((ki.eq.0).or.(ki.lt.-1))then;write(*,*)'wrong input: please choose a positive value';kii=kii-1;cycle;endif 
!          geness(kii)=ki          
!        end do
!        if(ki.eq.-1)then;geness(kii)=0;endif
!        
!        call pintagnufor(nodoss,geness)    ! >>>Miguel8-10-14
!case (1)
      !  if (flag(7)==1) then 
      !    flag(7)=0
      !  else 
      !    nki=0
 	  !if (allocated(oopp)) then ; deallocate(oopp) ; endif
      !    allocate(oopp(nd))
!56	  print*,"which node? (-1 to stop)"
	  !read(*,*) ki
      !    if (ki>0) then
      !      guardaene=node(ki)%e
      !      call energia(ki) 
  	  !  print*,"position",node(ki)%x,node(ki)%y,node(ki)%z
      !      print*,"down node",node(ki)%altre
	  !  print*,"energy",node(ki)%e
      !      print*,"altre",node(ki)%altre
      !      print*,"tipus",node(ki)%tipus
      !      print*,"cell",node(ki)%icel
      !      print *,node(ki)%altre,"altre"
      !      node(ki)%e=guardaene
      !      print*,"gene expression"
      !      print*,gex(ki,1:ng)
      !      print *,"node properties"
      !      print *,node(ki)
      !    end if
      !    nki=nki+1
      !    oopp(nki)=ki
      !    if (ki/=-1) goto 56
      !    nki=nki-1
      !    flag(7)=1
      !  end if
	!call draw_func
case (8)
        stop
!case (13)
      !  if (flag(7)==1) then 
      !    flag(7)=0
      !  else 
 	  !if (allocated(oopp)) then ; deallocate(oopp) ; endif
      !    allocate(oopp(nd))
      !    nki=0
!566	  print*,"which cell? (-1 to stop)"
	  !read(*,*) k
      !    if (k>0) then
      !      do kk=1,cels(k)%nunodes
      !        ki=cels(k)%node(kk)
      !        print *,kk,"th node in cell",k,"is",ki,"its tipus is",node(ki)%tipus
      !        oopp(kk)=ki
      !      end do
      !      nki=nki+cels(k)%nunodes
      !    end if
      !    if (k/=-1) goto 566
      !    flag(7)=1
      !  end if
	!call draw_func
case(14) 
  print *,"which node to move"
  read(*,*) node_to_move
  print *,"now do something about"
case(15)
  !call fit(1,a)

!  print*,"changing the node color code by setting other maximum and minimum values"
!  print*,"actual colorselection is",colorselection
!  print*,"enter new minimum and maximum values:"
!  read(*,*) custom_minval,custom_maxval
!  custom_colorselection=colorselection
!case(16)
!  print*,"changing the color palette used to paint nodes"
!  print*,"actual color palette is",conf_rainbow
!  print*,"choose new color palette"
!  print*,"1 = multicolor palette"
!  print*,"2 = yellow-red palette"
!  print*,"3 = blue palette"
!  print*,"4 = green palette"
!  read(*,*) conf_rainbow
end select

return
end subroutine menu_handler
!*******************************************SUBROUTINE********************************************************
subroutine param_handler(selection)
integer(kind=glint), value :: selection


end subroutine param_handler

!*******************************************SUBROUTINE********************************************************
subroutine menuu_handler(selection)
integer(kind=glint), value :: selection

!colorselection=selection

!if(colorselection>=36 .and. colorselection<=38) chogen=0

  select case(selection)
  case(1)
    arrowselection=0
  case(2)
    sphereselection=0
  case(3)
    write(*,*)'Enter new value for the arrow length corresponding to the maximum value of the property being printed'
    write(*,*)'Actual maximum value=',amaxe,'actual scale=',amax_scale,'-1 to reset default value (1)'
    read(*,*)amax_scale
    if(amax_scale==-1d0) amax_scale=1d0
  case(4)
    write(*,*)'Enter new value for the transparent sphere radius corresponding to the maximum value of the property being printed'
    write(*,*)'Actual maximum value=',smaxe,'actual scale=',smax_scale,'-1 to reset default value (1)'
    read(*,*)smax_scale
    if(smax_scale==-1d0) smax_scale=1d0
  end select

end subroutine

subroutine menuu2_handler(selection)
integer(kind=glint), value :: selection

!colorselection=selection

!if(colorselection>=nparam_per_node+1 .and. colorselection<=nparam_per_node+3)then
!  if(custom_colorselection/=colorselection)then
!    chogen=0
!  end if
!end if

  select case(selection)
  case(1)
    arrowselection=0
  case(2)
    sphereselection=0
  case(3)
    write(*,*)'Enter new value for the arrow length corresponding to the maximum value of the property being printed'
    write(*,*)'Actual maximum value=',amaxe,'actual scale=',amax_scale,'-1 to reset default value (1)'
    read(*,*)amax_scale
    if(amax_scale==-1d0) amax_scale=1d0
  case(4)
    write(*,*)'Enter new value for the transparent sphere radius corresponding to the maximum value of the property being printed'
    write(*,*)'Actual maximum value=',smaxe,'actual scale=',smax_scale,'-1 to reset default value (1)'
    read(*,*)smax_scale
    if(smax_scale==-1d0) smax_scale=1d0
  end select




end subroutine
!*******************************************SUBROUTINE********************************************************
subroutine cmenuu_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

colorselection=selection


call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe

!if(colorselection>=36 .and. colorselection<=38) chogen=0

end subroutine
!*******************************************SUBROUTINE********************************************************
subroutine cmenuu2_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

colorselection=selection

call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe

!if(colorselection>=nparam_per_node+1 .and. colorselection<=nparam_per_node+3)then
!  if(custom_colorselection/=colorselection)then
!    chogen=0
!  end if
!end if

end subroutine

!*******************************************SUBROUTINE********************************************************
subroutine amenuu_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

arrowselection=selection


call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe


end subroutine
!*******************************************SUBROUTINE********************************************************
subroutine amenuu2_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

arrowselection=selection


call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe


end subroutine
!*******************************************SUBROUTINE********************************************************
subroutine smenuu_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

sphereselection=selection


call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe

end subroutine
!*******************************************SUBROUTINE********************************************************
subroutine smenuu2_handler(selection)
integer(kind=glint), value :: selection
real*8::mine,maxe
integer::chogen

sphereselection=selection


call get_color_max_min(colorselection,mine,maxe,chogen)
print*,"minvalue",mine,"maxvalue",maxe

end subroutine


!********************************************************************************************************

subroutine menud_handler(selection)
integer(kind=glint), value :: selection
integer::i

if (flag(selection)==1) then ; flag(selection)=0 ; else ; flag(selection)=1 ; end if

if (selection==5)then
  if(flag(selection)==1 .and. flag(6)==1) flag(6)=0
elseif(selection==6)then
  if(flag(selection)==1 .and. flag(5)==1) flag(5)=0
end if

return
if(selection==22)then
  if(flag(22)==1)then
    flag(8)=1
  else
    flag(8)=0
  end if
end if

if(selection==8)then
  if(flag(8)==1)then
    flag(22)=1
  else
    flag(22)=0
  end if
end if

end subroutine

!*******************************************SUBROUTINE********************************************************
!subroutine menut_handler(selection)
!integer(kind=glint), value :: selection
!integer i,j,k
!
!	select case(selection)
!	case(1)
!		print*,"cell division: which cell?"
!		read(*,*) i
!		call division(i)
!	case(2)
!		if(ffu(3)==0)then
!			print*,"turning growth ON"
!			ffu(3)=1
!		else
!			print*,"turning growth OFF"
!			ffu(3)=0
!		end if
!	case(3)
!		do while (i/=-1)
!			print*,"invaginate: which cell? (-1 to stop)"
!			read(*,*) i
!			if(i>ncels)then
!				print*,"there's no such cell. actual number of cells:",ncels
!			elseif(i>0)then
!!				call invagination(i)
!			end if
!		end do
!	case(4)
!!		i=int(ran2(idum)*nd)+1
!		print*,"duplicate a node, in which cell?"
!		read(*,*) i
!	case(5)
!!		i=int(ran2(idum)*nd)+1
!		print*,"add an ECM node, from which cell? apical (1) or basal (2)?"
!		read(*,*) i,j
!	end select
!
!
!
!end subroutine

!*******************************************SUBROUTINE********************************************************
subroutine menuq_handler(selection)
integer(kind=glint), value :: selection
integer i

if (selection==-1) then
  print *,"how many?"
  read(*,*) i
else
  i=selection
end if
call go_iteration_back(i)
end subroutine

!*******************************************SUBROUTINE******************************************************** !>>>>>>Tommi, editor main menu handler 11.9.2013

subroutine menueditor_handler(selection)
integer(kind=glint), value :: selection
!call view_from_front
integer qgen,altro

select case(selection)
case(1)
         choice=1
         flag(40)=1
         call addnode(cursx,cursy,cursz)
case(2)	
	choice=2
        flag(40)=1
        left_button_func = CURPAN
case(3)
     	mouseswitch = .false.
        switch=.false.
        flag(40)=0
case(4)
case(5)
case(6)
      	choice=6
        flag(40)=1
        target=0
        call selectnode(cursx,cursy,cursz)
case(7)
        if(nki>0 .and. flag(3)==1)then
          choice=7
          flag(40)=1
          call pastenode(cursx,cursy,cursz)
        else
          print*,"*************************************************"
          print*,"no node selected, please select one or more nodes"
          print*,"*************************************************"
        end if
case(8)
	choice=8
        flag(40)=1
        call addcell(cursx,cursy,cursz)
case(9)
        colorselection=24 !Tommi 10.1.2014
	choice=9
        flag(40)=1
        call selectcell(cursx,cursy,cursz)
case(10)
        if(nki>0.and.flag(3)==1)then
          altro=node(nodeindex)%altre
          call deletenode
          flag(40)=1
          !call pastecell(cursx,cursy,cursz)

        !print*,"oopp",oopp(1:nki)
        !print*,"nodeindex",nodeindex,"altre",altro

          oopp(nki)=0

          nki=nki-1
          if(nki>0)then 
            do i=1,nki
              if(node(nodeindex)%tipus>=3)then
                if(oopp(i)>nodeindex) oopp(i)=oopp(i)-1
              else
                if(oopp(i)>nodeindex .and. oopp(i)>altro)then
                  oopp(i)=oopp(i)-2 !; print*,"minus 2"
                elseif(oopp(i)>nodeindex .and. oopp(i)<altro)then
                  oopp(i)=oopp(i)-1 !; print*,"minus 1.0"
                elseif(oopp(i)<nodeindex .and. oopp(i)>altro)then
                  oopp(i)=oopp(i)-1 !; print*,"minus 1.1"
                end if
              end if
            end do
            nodeindex=oopp(nki)
          else
            flag(3)=0 ; flag(7)=0
            nodeindex=0
          end if
        else
          print*,"*************************************************"
          print*,"no node selected, please select one or more nodes"
          print*,"*************************************************"
        end if


case(11)
        if(nki>0.and.flag(4)==1)then
          call deletecell(cellid)
          flag(40)=1

          oopp(nki)=0
          nki=nki-1
          if(nki>0)then
            do i=1,nki
              if(oopp(i)>cellid) oopp(i)=oopp(i)-1
            end do
            cellid=oopp(nki)
          else
            flag(4)=0 ; flag(7)=0
            cellid=0
          end if
        else
          print*,"*************************************************"
          print*,"no cell selected, please select one or more cells"
          print*,"*************************************************"
        end if


case(12)
         choice=12
         flag(40)=1
         call addmescell(cursx,cursy,cursz)
case(13)
         choice=13
         flag(40)=1
         call addecm(cursx,cursy,cursz)
case(14)
      	choice=14
        target=1
        call selectnode(cursx,cursy,cursz)
        flag(40)=1
case(15)
        left_button_func = CURPANZ
case(16)
        if(nki>0 .and. flag(4)==1)then
          flag(40)=1
          call pastecell(cursx,cursy,cursz)
        else
          print*,"*************************************************"
          print*,"no node selected, please select one or more nodes"
          print*,"*************************************************"
        end if

case(17)
        flag(40)=1
        print *,"which gene"
        read(*,*) qgen
        print *,"the old value was",gex(nodeindex,qgen),"give new value"
        read(*,*) a
        gex(nodeindex,qgen)=a


end select

end subroutine !>>>>>Tommi

!*******************************************SUBROUTINE********************************************************

subroutine sel_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection
!call view_from_front
!integer qgen
real*8 :: guardaene

  select case(selection)
  case(1)
    !if (flag(7)==1) then 
    !  flag(7)=0
    !else 
      nki=0
      if (allocated(oopp)) then ; deallocate(oopp) ; endif
      allocate(oopp(nd))
56	  print*,"which node? (-1 to stop)"
	  read(*,*) ki
      if (ki>0) then
        guardaene=node(ki)%e
        call energia(ki) 
        print*,"position",node(ki)%x,node(ki)%y,node(ki)%z
        print*,"down node",node(ki)%altre
	    print*,"energy",node(ki)%e
        print*,"altre",node(ki)%altre
        print*,"tipus",node(ki)%tipus
        print*,"cell",node(ki)%icel
        print *,node(ki)%altre,"altre"
        node(ki)%e=guardaene
        print*,"gene expression"
        print*,gex(ki,1:ng)
        !print *,"node properties"
        !print *,node(ki)
      end if
      nki=nki+1
      oopp(nki)=ki
      nodeindex=ki
      if (ki/=-1) goto 56
      nki=nki-1
      flag(7)=1 ; flag(3)=1 ; flag(4)=0
    !end if
	call draw_func
  case(2)

    print*, "     " !!>> HC 16-11-2020
    print*,"CURSOR MODE ON: DRAG WITH LEFT MOUSE-BUTTON TO MOVE CURSOR ON X-Y PLANE"
    print*,"                CLICK MIDDLE MOUSE TO SELECT NODE"
    print*,"                PRESS DOWN ARROW KEY TO SWITCH BETWEEN MOVING THE CURSOR"
    print*,"                ALONG X-Y PLANE AND Z AXIS"
    print*, "     " !!>> HC 16-11-2020

    if(nki==0)then
      if (allocated(oopp)) then ; deallocate(oopp) ; endif
      allocate(oopp(nd))
    end if
    flag(7)=1 ; flag(3)=1 ;flag(4)=0

    flag(40)=1
    left_button_func=CURPAN
    middle_press=SELNO
	call draw_func


  case (3)
!    if (flag(7)==1) then 
!      flag(7)=0
!    else 
      if(nki==0)then
        if (allocated(oopp)) then ; deallocate(oopp) ; endif
        allocate(oopp(nd))
      end if
566	  print*,"which cell? (-1 to stop)"
	  read(*,*) k
      if (k>0) then
        !do kk=1,cels(k)%nunodes
          nki=nki+1
          !ki=cels(k)%node(kk)
          !print *,kk,"th node in cell",k,"is",ki,"its tipus is",node(ki)%tipus
          oopp(nki)=k
        !end do
        cellid=k
      end if
      if (k/=-1) goto 566
      flag(7)=1 ; flag(4)=1 ; flag(3)=0
!    end if
    call draw_func
  case (4)

    print*, "     " !!>> HC 16-11-2020
    print*,"CURSOR MODE ON: DRAG WITH LEFT MOUSE-BUTTON TO MOVE CURSOR ON X-Y PLANE"
    print*,"                CLICK MIDDLE MOUSE TO SELECT NODE"
    print*,"                PRESS DOWN ARROW KEY TO SWITCH BETWEEN MOVING THE CURSOR"
    print*,"                ALONG X-Y PLANE AND Z AXIS"
    print*, "     " !!>> HC 16-11-2020

    if(nki==0)then
      if (allocated(oopp)) then ; deallocate(oopp) ; endif
      allocate(oopp(nd))
    end if
    flag(7)=1; flag(4)=1 ; flag(3)=0

    flag(40)=1
    left_button_func=CURPAN
    middle_press=SELCEL
	call draw_func

  case (5)
      flag(7)=0 ; flag(3)=0 ; flag(4)=0
      nodeindex=0

      nki=0
      if (allocated(oopp)) then ; deallocate(oopp) ; endif

      flag(40)=0
      if(left_button_func==CURPAN.or.left_button_func==CURPANZ) &
                                           & left_button_func=ROTATE
      if(middle_button_func==SELNO.or.middle_button_func==SELCEL) &
                                           & left_button_func=ZOOM
      middle_press=0

	call draw_func

  end select

end subroutine !>>>>>Miquel

subroutine cursor_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection
!call view_from_front
!integer qgen
real*8 :: guardaene

  select case(selection)
  case(1)
    if(flag(40)==0)then
      flag(40)=1
      left_button_func=CURPAN
      middle_button_func=CURPANZ
    else
      flag(40)=0
      left_button_func=ROTATE
      middle_button_func=ZOOM
    end if
	call draw_func
  case(2)
    flag(40)=1
    left_button_func=CURPAN
	call draw_func
  case (3)
    flag(40)=1
    left_button_func=CURPANZ
    call draw_func
  case (4)
    flag(40)=1
    middle_button_func=CURPANZ
    call draw_func
  end select

end subroutine !>>>>>Miquel

subroutine color_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection
!call view_from_front
!integer qgen
!real*8 :: guardaene

  select case(selection)
  !case(1)
  !  print*,"changing the color palette used to paint nodes"
  !  print*,"actual color palette is",conf_rainbow
  !  print*,"choose new color palette"
  !  print*,"1 = multicolor palette"
  !  print*,"2 = yellow-red palette"
  !  print*,"3 = blue palette"
  !  print*,"4 = green palette"
  !  read(*,*) conf_rainbow
  !	call draw_func
  case(2)
    print*,"changing the node color code by setting other maximum and minimum values"
    print*,"actual colorselection is",colorselection
    print*,"enter new minimum and maximum values:"
    read(*,*) custom_minval,custom_maxval
    custom_colorselection=colorselection
  	call draw_func
  case (3)
    print*,"changing the node color code by a new minimum value"
    print*,"actual colorselection is",colorselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=COLORMIN

    middle_press=0

    if(custom_colorselection/=colorselection)then
      custom_colorselection=colorselection
      call get_color_max_min(custom_colorselection,custom_minval,custom_maxval,chogen)
    end if

    call draw_func
  case (4)
    print*,"changing the node color code by a new maximum value"
    print*,"actual colorselection is",colorselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=COLORMAX
    middle_press=0

    if(custom_colorselection/=colorselection)then
      custom_colorselection=colorselection
      call get_color_max_min(custom_colorselection,custom_minval,custom_maxval,chogen)
    end if

    call draw_func
  case (5)
    print*,"changing the node color code by a new minimum value"
    print*,"actual colorselection is",colorselection
    print*,"drag the mouse with middle button pressed vertically to change the value"

    middle_button_func=COLORMIN
    middle_press=0

    if(custom_colorselection/=colorselection)then
      custom_colorselection=colorselection
      call get_color_max_min(custom_colorselection,custom_minval,custom_maxval,chogen)
    end if

    call draw_func
  case (6)
    print*,"changing the node color code by a new maximum value"
    print*,"actual colorselection is",colorselection
    print*,"drag the mouse with middle button pressed vertically to change the value"

    middle_button_func=COLORMAX
    middle_press=0

    if(custom_colorselection/=colorselection)then
      custom_colorselection=colorselection
      call get_color_max_min(custom_colorselection,custom_minval,custom_maxval,chogen)
    end if

    call draw_func
  end select

end subroutine !>>>>>Miquel
!*************************************************

subroutine palette_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection

 conf_rainbow=selection
 call draw_func


end subroutine !>>>>>Miquel


subroutine arrow_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection
!call view_from_front
!integer qgen
!real*8 :: guardaene

  select case(selection)
  case(1)
    print*,"drag the mouse with left button pressed vertically to change the scale"
    left_button_func=ARROWSCALE
      middle_press=0
  case(2)
    print*,"changing the node arrow length by setting other maximum and minimum values"
    print*,"actual arrowselection is",arrowselection
    print*,"enter new minimum and maximum values:"
    read(*,*) custom_aminval,custom_amaxval
    custom_arrowselection=arrowselection
  	call draw_func
  case (3)
    print*,"changing the node arrow length by a new minimum value"
    print*,"actual arrowselection is",arrowselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=ARROWMIN
    middle_press=0
    if(custom_arrowselection/=arrowselection)then
      custom_arrowselection=arrowselection
      call get_color_max_min(custom_arrowselection,custom_aminval,custom_amaxval,achogen)
    end if

    call draw_func
  case (4)
    print*,"changing the node arrow length by a new maximum value"
    print*,"actual arrowselection is",arrowselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=ARROWMAX
    middle_press=0
    if(custom_arrowselection/=arrowselection)then
      custom_arrowselection=arrowselection
      call get_color_max_min(custom_arrowselection,custom_aminval,custom_amaxval,achogen)
    end if

    call draw_func
  case (5)
    print*,"changing the node arrow length by a new minimum value"
    print*,"actual arrowselection is",arrowselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    middle_button_func=ARROWMIN
    middle_press=0
    if(custom_arrowselection/=arrowselection)then
      custom_arrowselection=arrowselection
      call get_color_max_min(custom_arrowselection,custom_aminval,custom_amaxval,achogen)
    end if

    call draw_func
  case (6)
    print*,"changing the node arrow length by a new maximum value"
    print*,"actual arrowselection is",arrowselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    middle_button_func=ARROWMAX
    middle_press=0
    if(custom_arrowselection/=arrowselection)then
      custom_arrowselection=arrowselection
      call get_color_max_min(custom_arrowselection,custom_aminval,custom_amaxval,achogen)
    end if

    call draw_func
  case (7)
    arrowselection=0
    custom_arrowselection=0
    custom_amaxval=0 ;custom_aminval=0
    left_button_func=ROTATE
    middle_button_func=ZOOM
    middle_press=0
    call draw_func
  end select

end subroutine !>>>>>Miquel

subroutine sphere_menu_handler(selection)  !>>Miquel24-10-14
integer(kind=glint), value :: selection
!call view_from_front
!integer qgen
!real*8 :: guardaene

  select case(selection)
  case(1)
    print*,"drag the mouse with left button pressed vertically to change the scale"
    left_button_func=SPHERESCALE
      middle_press=0
  case(2)
    print*,"changing the node sphere radius by setting other maximum and minimum values"
    print*,"actual sphereselection is",sphereselection
    print*,"enter new minimum and maximum values:"
    read(*,*) custom_sminval,custom_smaxval
    custom_sphereselection=sphereselection
  	call draw_func
  case (3)
    print*,"changing the node sphere radius by a new minimum value"
    print*,"actual sphereselection is",sphereselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=SPHEREMIN
    middle_press=0
    if(custom_sphereselection/=sphereselection)then
      custom_sphereselection=sphereselection
      call get_color_max_min(custom_sphereselection,custom_sminval,custom_smaxval,schogen)
    end if

    call draw_func
  case (4)
    print*,"changing the node sphere radius by a new maximum value"
    print*,"actual sphereselection is",sphereselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    left_button_func=SPHEREMAX
    middle_press=0
    if(custom_sphereselection/=sphereselection)then
      custom_sphereselection=sphereselection
      call get_color_max_min(custom_sphereselection,custom_sminval,custom_smaxval,schogen)
    end if

    call draw_func
  case (5)
    print*,"changing the node sphere radius by a new minimum value"
    print*,"actual sphereselection is",sphereselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    middle_button_func=SPHEREMIN
    middle_press=0
    if(custom_sphereselection/=sphereselection)then
      custom_sphereselection=sphereselection
      call get_color_max_min(custom_sphereselection,custom_sminval,custom_smaxval,schogen)
    end if

    call draw_func
  case (6)
    print*,"changing the node sphere radius by a new maximum value"
    print*,"actual sphereselection is",sphereselection
    print*,"drag the mouse with left button pressed vertically to change the value"

    middle_button_func=SPHEREMAX
    middle_press=0
    if(custom_sphereselection/=sphereselection)then
      custom_sphereselection=sphereselection
      call get_color_max_min(custom_sphereselection,custom_sminval,custom_smaxval,schogen)
    end if

    call draw_func
  case (7)
    sphereselection=0
    custom_sphereselection=0
    custom_smaxval=0 ;custom_sminval=0
    left_button_func=ROTATE
    middle_button_func=ZOOM
    call draw_func
  end select

end subroutine !>>>>>Miquel


!************************************************* !>>>Miguel 8-10-14 made subroutine

subroutine get_color_max_min(selection,mine,maxe,chogen)
 integer, intent(in):: selection
 real*8, intent(out)::mine,maxe
 integer,intent(out)::chogen
real*8         :: disx(nd),disy(nd),disz(nd),distot(nd), discentch(nd)

  do i=1,nd                                                                      !!>> HC 12-2-2024
     discentch(i)=(sqrt( (node(i)%x)**2 + (node(i)%y)**2 + (node(i)%z)**2))**2   !!>> HC 12-2-2024
  enddo                                                                          !!>> HC 12-2-2024
  
  !BALL COLORS MAX AND MINS

    select case(selection)	
    case(1) ; mine=minval(node(:nd)%x) ; maxe=maxval(node(:nd)%x)
    case(2) ; mine=minval(node(:nd)%y) ; maxe=maxval(node(:nd)%y)
    case(3) ; mine=minval(node(:nd)%z) ; maxe=maxval(node(:nd)%z)
    case(4) ; mine=minval(node(:nd)%e) ; maxe=maxval(node(:nd)%e)
    case(5) ; mine=minval(node(:nd)%eqd) ; maxe=maxval(node(:nd)%eqd) 
    case(6) ; mine=minval(node(:nd)%add)  ; maxe=maxval(node(:nd)%add)
    case(7) ; mine=minval(node(:nd)%you) ; maxe=maxval(node(:nd)%you)
    case(8) ; mine=minval(node(:nd)%adh) ; maxe=maxval(node(:nd)%adh)
    case(9) ; mine=minval(node(:nd)%rep) ; maxe=maxval(node(:nd)%rep)	
    case(10) ; mine=minval(node(:nd)%rec) ; maxe=maxval(node(:nd)%rec)
    case(11) ; mine=minval(node(:nd)%erp) ; maxe=maxval(node(:nd)%erp)
    case(12) ; mine=minval(node(:nd)%est) ; maxe=maxval(node(:nd)%est)	
    case(13) ; mine=minval(node(:nd)%eqs) ; maxe=maxval(node(:nd)%eqs)
    case(14) ; mine=minval(node(:nd)%hoo) ; maxe=maxval(node(:nd)%hoo) 
    case(15) ; mine=minval(node(:nd)%mov) ; maxe=maxval(node(:nd)%mov)
    case(16) ; mine=minval(node(:nd)%dmo) ; maxe=maxval(node(:nd)%dmo)
    case(17) ; mine=minval(node(:nd)%orix) ; maxe=maxval(node(:nd)%orix) 
    case(18) ; mine=minval(node(:nd)%oriy) ; maxe=maxval(node(:nd)%oriy)
    case(19) ; mine=minval(node(:nd)%oriz) ; maxe=maxval(node(:nd)%oriz)
    case(20) ; mine=minval(node(:nd)%ecm) ; maxe=maxval(node(:nd)%ecm)
    case(21) ; mine=minval(node(:nd)%cod) ; maxe=maxval(node(:nd)%cod)
    case(22) ; mine=minval(node(:nd)%grd) ; maxe=maxval(node(:nd)%grd)
    case(23) ; mine=minval(node(:nd)%pld) ; maxe=maxval(node(:nd)%pld)
    case(24) ; mine=minval(node(:nd)%vod) ; maxe=maxval(node(:nd)%vod)
    case(25) ; mine=minval(node(:nd)%dif) ; maxe=maxval(node(:nd)%dif)
    case(26) ; mine=minval(node(:nd)%kfi) ; maxe=maxval(node(:nd)%kfi)
    case(27) ; mine=minval(node(:nd)%pla) ; maxe=maxval(node(:nd)%pla)
    case(28) ; mine=minval(node(:nd)%kvol) ; maxe=maxval(node(:nd)%kvol)
    case(29) ; mine=1 ;maxe=7 !mine=minval(node(:nd)%tipus) ; maxe=maxval(node(:nd)%tipus) !this is tipus
    case(30) ; mine=minval(node(:nd)%icel) ; maxe=maxval(node(:nd)%icel)
    case(31) ; mine=minval(node(:nd)%altre) ; maxe=maxval(node(:nd)%altre)
    case(32) ; mine=minval(node(:nd)%marge) ; maxe=maxval(node(:nd)%marge) 
    case(33) ; mine=minval(node(:nd)%talone) ; maxe=maxval(node(:nd)%talone) 
    case(34) ; mine=0;maxe=1!minval(node(:nd)%fix) ; maxe=maxval(node(:nd)%fix) 
    case(35) ; mine=minval(node(:nd)%ndiv) ; maxe=maxval(node(:nd)%ndiv)     !! HC 16-3-2022
    case(36) ; mine=minval(discentch(:nd)) ; maxe=maxval(discentch(:nd))     !! HC 14-2-2024

    !case(36) ; if (chogen==0) then
    !           chogen=1
    !           mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !           print*,"maxim",maxe,"minim",mine; end if
    !           mine=minval(gex(:nd,1)) ; maxe=maxval(gex(:nd,1))
    !case(37) ; if (chogen==0) then
    !             if (ng>1) then
    !               chogen=2
    !                mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !               print*,"maxim",maxe,"minim",mine; 
    !             else
    !               print *,"WARNING, not such a gene"
    !             end if
    !           else
    !             mine=minval(gex(:nd,2)) ; maxe=maxval(gex(:nd,2))
    !           end if
    !case(38) ; if (chogen==0) then
    !           print *,"which gene" ; read(*,*) chogen ; print *,chogen,"chosen regulatory molecule"
    !           if (chogen>ng) then 
    !             print *,"WARNING, not such a gene"
    !           else
    !             mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !             print*,"maxim",maxe,"minim",mine; 
    !           end if
    !             mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
    !           end if
    case(37) ; mine=0d0 ; maxe=1d0
    case(38) ; do i=1,nd ; disx(i)=node(i)%x-node(i)%orix ; end do 
               mine=minval(disx) ; maxe=maxval(disx)
    case(39) ; do i=1,nd ; disy(i)=node(i)%y-node(i)%oriy ; end do 
               mine=minval(disy) ; maxe=maxval(disy)
    case(40) ; do i=1,nd ; disz(i)=node(i)%z-node(i)%oriz ; end do 
               mine=minval(disz) ; maxe=maxval(disz)
    case(41) ; do i=1,nd 
               distot(i)=sqrt((node(i)%x-node(i)%orix)**2+(node(i)%y-node(i)%oriy)**2+(node(i)%z-node(i)%oriz)**2) 
               end do 
               mine=minval(distot) ; maxe=maxval(distot)
    case(42) ; mine=1 ; maxe=nd
    case(43) ; mine=1 ; maxe=nd
    case(44) ; mine=1; maxe=maxval(nodeo(:nd)%icel)  !!>> HC 11-1-2021
    end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !>>>Miguel30-10-14 
    if(selection.gt.44)then                          !!>> HC 11-1-2021
      chogen=selection-44                            !!>> HC 11-1-2021
      mine=minval(gex(:nd,chogen)) ; maxe=maxval(gex(:nd,chogen))
      !print*,"chosen gene",chogen,"maxim",maxe,"minim",mine
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine move_menu_handler(selection)  !>>Miquel31-10-14
integer(kind=glint), value :: selection
integer:: ki
real*8::x,y,z,dx,dy,dz
!call view_from_front
!integer qgen
!real*8 :: guardaene

  select case(selection)
  case(1)
    if(flag(3)==1)then
      print*,"node already selected,",nodeindex
      ki=nodeindex
    else
      print*,"no node selected, enter node index:"
      read(*,*)ki
    end if
    print *,"to which x (from x=",node(ki)%x,") -1 keep the same value"
    read(*,*) a
    if (a/=-1) node(ki)%x=a
    print *,"to which y (from y=",node(ki)%y,") -1  keep the same value"
    read(*,*) a
    if (a/=-1) node(ki)%y=a
    print *,"to which z (from z=",node(ki)%z,") -1 keep the same value"
    read(*,*) a
    if (a/=-1) node(ki)%z=a
    call iniboxes
	call draw_func
  case(2)
    if(flag(3)==1)then
      print*,"node already selected,",nodeindex
      ki=nodeindex
    else
      print*,"no node selected, enter node index:"
      read(*,*)ki
      nodeindex=ki
      flag(7)=1 ; flag(3)=1
      if(nki==0)then
        nki=1
        if(allocated(oopp)) deallocate (oopp)
        allocate(oopp(nki))
        oopp(1)=nodeindex
      end if
    end if
    left_button_func=NOPAN
    middle_press=0
    call draw_func
  case(3)
    if(flag(4)==1)then
      print*,"cell already selected,",cellid
      ki=cellid
    else
      print*,"no cell selected, enter cell index:"
      read(*,*)ki
    end if
    print *,"Move the cell by setting the new coordinates of its centroid"
    print *,"Current centroid position is:",cels(ki)%cex,",",cels(ki)%cey,",",cels(ki)%cez
    print *,"to which x (from x=",node(ki)%x,") -1 keep the same value"
    read(*,*) a
    if (a==-1) x=cels(ki)%cex
    print *,"to which y (from y=",node(ki)%y,") -1  keep the same value"
    read(*,*) a
    if (a==-1) y=cels(ki)%cey
    print *,"to which z (from z=",node(ki)%z,") -1 keep the same value"
    read(*,*) a
    if (a==-1) z=cels(ki)%cez

    dx=x-cels(ki)%cex
    dy=y-cels(ki)%cey
    dz=z-cels(ki)%cez

    do i=1,cels(ki)%nunodes
      j=cels(ki)%node(i)
      node(j)%x=node(j)%x+dx
      node(j)%y=node(j)%y+dy
      node(j)%z=node(j)%z+dz
    end do
    cels(ki)%cex=x
    cels(ki)%cey=y
    cels(ki)%cez=z

    call iniboxes
	call draw_func
  case(4)
    if(flag(4)==1)then
      print*,"cell already selected,",cellid
      ki=cellid
    else
      print*,"no cell selected, enter cell index:"
      read(*,*)ki
      cellid=ki
      flag(7)=1 ; flag(4)=1
      if(nki==0)then
        nki=cels(ki)%nunodes
        if(allocated(oopp)) deallocate (oopp)
        allocate(oopp(nki))
        do i=1,nki
          j=cels(ki)%node(i)
          oopp(i)=j
        end do
      end if
    end if
    left_button_func=CEPAN
    middle_press=0
    call draw_func
  case(5)
    if(left_button_func==NOPAN.or.left_button_func==CEPAN) &
                                           & left_button_func=ROTATE
    if(left_button_func==NOPANZ.or.left_button_func==CEPANZ) &
                                           & left_button_func=ROTATE
    middle_press=0
  end select

end subroutine !>>>>>Miquel
!*************************************************


subroutine pintagnufor(nodoss,geness)
 integer ::ki,kii,el,ell,i,j,ii,jj,im,imm
 character(len=40)  ::                  arxiv  
 integer :: nodoss(10),geness(10)
 integer :: ret

ret=SYSTEM('rm _N_______*')
 
el=getott!int(real(getot)/real(fprint))+1
imm=0
do i=1,10
  ki=nodoss(i)   
  if(ki.ne.0)then
    do j=1,10
      kii=geness(j)
      if(kii.ne.0)then   
        imm=imm+1
        write(arxiv,*)"N",ki,"G",kii
        do im=1,40
          if (arxiv(im:im)==" ")then; arxiv(im:im)="_";end if    
          arxiv(im:im)=arxiv(im:im)
        end do
        open(11111,file=arxiv,status='unknown',action='write')   
        do ii=1,el
          write(11111,*)ejex1(ii), gext1(ki,kii,ii)
          !write(11111,*)ii*real(fprint)-fprint, gext1(ki,kii,ii) 
        end do
        close(11111)
      end if
    end do
  end if
end do
if(imm.gt.0)then 
   call f2gp(nodoss,geness,imm)
end if
  
end subroutine pintagnufor


!*******************************************SUBROUTINE********************************************************

subroutine planemenu_handler(selection) !<<<<<Tommi, planemenu handler 11.9.2013
integer(kind=glint), value :: selection
select case(selection)
case(1)
         plane=1
         call view_from_front
case(2)
         plane=2
         call view_from_front
case(3)
         plane=3
	 call view_from_front
case(4)
         plane=4
	 call view_from_front
case(5)
         plane=5
	 call view_from_front
case(6)
         plane=6
	 call view_from_front
end select
end subroutine!>>>>>Tommi

!*******************************************SUBROUTINE********************************************************
subroutine nodecopymenu_handler(selection) !<<<<<<Tommi 23.9.2013 Node copying menu
integer(kind=glint), value :: selection

if (nodeindex<=0) then
print*,"No node selected! Please select a node"
else
select case(selection)
case(1)
	 print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%x
	 tempnod%x=node(nodeindex)%x
         print*,"Now, go and chose Paste"
print *,nodeindex,tempnod%x,"uuuuuuuuuu"
         
case(2)
	 print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%y
	 tempnod%y=node(nodeindex)%y
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=2
	 call menueditor_handler(6)
      
case(3)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%z
	 tempnod%z=node(nodeindex)%z
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=3
         call menueditor_handler(6)
case(4)
	 print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%e
	 tempnod%e=node(nodeindex)%e
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=4
         call menueditor_handler(6)
case(5)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%eqd
	 tempnod%eqd=node(nodeindex)%eqd
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=5
	 call menueditor_handler(6)
!case(6)
!         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%codel
!	 tempnod%codel=node(nodeindex)%codel
!         print*,"Now, please select a new node where you want to copy this property and select paste property"
!         property=6
!	 call menueditor_handler(6)
case(6)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%add
	  tempnod%add=node(nodeindex)%add
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=6
	 call menueditor_handler(6)
case(7)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%you
	  tempnod%you=node(nodeindex)%you
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=7
	 call menueditor_handler(6)
case(8)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%adh
	  tempnod%adh=node(nodeindex)%adh
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=8
	 call menueditor_handler(6)
case(9)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%rep
	 tempnod%rep=node(nodeindex)%rep
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=9
	 call menueditor_handler(6)
case(10)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%rec
	  tempnod%rec=node(nodeindex)%rec
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=10
	 call menueditor_handler(6)
case(11)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%erp
	  tempnod%erp=node(nodeindex)%erp
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=11
	 call menueditor_handler(6)
!case(12)
!         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%est
!	  tempnod%est=node(nodeindex)%est
!         print*,"Now, please select a new node where you want to copy this property and select paste property"
!         property=12
!	 call menueditor_handler(6)
case(14)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%eqs
	  tempnod%eqs=node(nodeindex)%eqs
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=12
	 call menueditor_handler(6)
case(15)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%hoo
	 tempnod%hoo=node(nodeindex)%hoo
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=13
	 call menueditor_handler(6)
case(16)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%mov
	  tempnod%mov=node(nodeindex)%mov
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=14
	 call menueditor_handler(6)
case(17)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%dmo
	  tempnod%dmo=node(nodeindex)%dmo
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=15
	 call menueditor_handler(6)
case(18) 
          print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%orix
	  tempnod%orix=node(nodeindex)%orix
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=16
	 call menueditor_handler(6)
case(19) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%oriy
	  tempnod%oriy=node(nodeindex)%oriy
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=17
	 call menueditor_handler(6)
case(20) 
          print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%oriz
	  tempnod%oriz=node(nodeindex)%oriz
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=18
	 call menueditor_handler(6)
case(21) 
          print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%ecm
	  tempnod%ecm=node(nodeindex)%ecm
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=19
	 call menueditor_handler(6)
case(22)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%tipus
	 tempnod%tipus=node(nodeindex)%tipus
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=20
	 call menueditor_handler(6)
case(23)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%icel
	tempnod%icel=node(nodeindex)%icel
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=21
	 call menueditor_handler(6)
case(24)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%altre
	tempnod%altre=node(nodeindex)%altre
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=22
	 call menueditor_handler(6)
case(25)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%marge
	tempnod%marge=node(nodeindex)%marge
         print*,"Now, please select a new node where you want to copy this property and select paste property"
         property=23
	 call menueditor_handler(6)
case(26)
         if(property<=0) then
            print*,"No property of a node selected! Please select a property"
         else
	 	 select case (property)
		 case(1)
	         colorselection=Z
		 node(nodeindex)%x=tempnod%x
		 print*,"Property changed"
		 case(2)
                 colorselection=2
		 node(nodeindex)%y=tempnod%y
		 print*,"Property changed"
		 case(3) 
                 colorselection=3
		 node(nodeindex)%z=tempnod%z
		 print*,"Property changed"
	 	 case(4)
	         colorselection=4
     		 node(nodeindex)%e=tempnod%e
		 print*,"Property changed"
		 case(5)
                 colorselection=5
		 node(nodeindex)%eqd=tempnod%eqd
		 print*,"Property changed"
!		 case(6) 
!                 colorselection=6
!		 node(nodeindex)%codel=tempnod%codel
!		 print*,"Property changed"
                 case(6)
	         colorselection=6
		 node(nodeindex)%add=tempnod%add
		 print*,"Property changed"
		 case(7)
                 colorselection=7
		 node(nodeindex)%you=tempnod%you
		 print*,"Property changed"
		 case(8) 
                 colorselection=8
		 node(nodeindex)%adh=tempnod%adh
		 print*,"Property changed"
	 	 case(9)
	         colorselection=9
		 node(nodeindex)%rep=tempnod%rep
		 print*,"Property changed"
		 case(10)
                 colorselection=10
		 node(nodeindex)%rec=tempnod%rec
		 print*,"Property changed"
		 case(11) 
                 colorselection=11
		 node(nodeindex)%erp=tempnod%erp
		 print*,"Property changed"
		 !case(12)
	         !colorselection=12
		 !node(nodeindex)%est=tempnod%est
		 !print*,"Property changed"
		 case(12)
                 colorselection=12
		 node(nodeindex)%eqs=tempnod%eqs
		 print*,"Property changed"
		 case(13) 
                 colorselection=13
		 node(nodeindex)%hoo=tempnod%hoo
		 print*,"Property changed"
	 	 case(14)
	         colorselection=14
     		 node(nodeindex)%mov=tempnod%mov
		 print*,"Property changed"
		 case(15)
                 colorselection=15
		 node(nodeindex)%dmo=tempnod%dmo
		 print*,"Property changed"
		 case(16) 
                 colorselection=16
		 node(nodeindex)%orix=tempnod%orix
		 print*,"Property changed"
                 case(17)
	         colorselection=17
		 node(nodeindex)%oriy=tempnod%oriy
		 print*,"Property changed"
		 case(18)
                 colorselection=18
		 node(nodeindex)%oriz=tempnod%oriz
		 print*,"Property changed"
		 case(19) 
                 colorselection=19
		 node(nodeindex)%ecm=tempnod%ecm
		 print*,"Property changed"
	 	 case(20)
	         colorselection=20
		 node(nodeindex)%tipus=tempnod%tipus
		 print*,"Property changed"
		 case(21)
                 colorselection=21
		 node(nodeindex)%icel=tempnod%icel
		 print*,"Property changed"
		 case(22) 
                 colorselection=22
		 node(nodeindex)%altre=tempnod%altre
		 case(23) 
                 colorselection=23
		 node(nodeindex)%marge=tempnod%marge
		 print*,"Property changed"
		 end select
         endif
end select

endif

end subroutine

!*******************************************SUBROUTINE********************************************************
!*******************************************SUBROUTINE********************************************************
subroutine nodepastemenu_handler(selection)
integer(kind=glint), value :: selection

if (nodeindex<=0) then
print*,"No node selected! Please select a node"
else
print*,"We copy ",nodeparams(selection),"from",nodeindex,"to",nodetarget
select case(selection)
case(1) ; node(nodetarget)%x=tempnod%x
case(2) ; node(nodetarget)%y=tempnod%y      
case(3) ; node(nodetarget)%z=tempnod%z
case(4) ; node(nodetarget)%e=tempnod%e 
case(5) ; node(nodetarget)%eqd=tempnod%eqd
case(6) ; node(nodetarget)%add=tempnod%add
case(7) ; node(nodetarget)%you=tempnod%you
case(8) ; node(nodetarget)%adh=tempnod%adh
case(9) ; node(nodetarget)%rep=tempnod%rep
case(10); node(nodetarget)%rec=tempnod%rec
case(11); node(nodetarget)%erp=tempnod%erp
!case(12); node(nodetarget)%est=tempnod%est
case(12); node(nodetarget)%eqs=tempnod%eqs
case(13); node(nodetarget)%hoo=tempnod%hoo
case(14); node(nodetarget)%mov=tempnod%mov
case(15); node(nodetarget)%dmo=tempnod%dmo
case(16); node(nodetarget)%orix=tempnod%orix
case(17); node(nodetarget)%oriy=tempnod%oriy
case(18); node(nodetarget)%oriz=tempnod%oriz
case(19); node(nodetarget)%ecm=tempnod%ecm
case(20); node(nodetarget)%cod=tempnod%cod
case(21); node(nodetarget)%grd=tempnod%grd
case(22); node(nodetarget)%pld=tempnod%pld
case(23); node(nodetarget)%vod=tempnod%vod
case(24); node(nodetarget)%dif=tempnod%dif
case(25); node(nodetarget)%kfi=tempnod%kfi
case(26); node(nodetarget)%pla=tempnod%pla
case(27); node(nodetarget)%kvol=tempnod%kvol
!case(29); node(nodetarget)%temt=tempnod%temt
case(28); node(nodetarget)%tipus=tempnod%tipus
case(29); node(nodetarget)%icel=tempnod%icel
case(30); node(nodetarget)%altre=tempnod%altre
case(31); node(nodetarget)%marge=tempnod%marge
case(32); node(nodetarget)%talone=tempnod%talone
case(33); node(nodetarget)%fix=tempnod%fix
case(35); node(nodetarget)%ndiv=tempnod%ndiv  !!>> HC 16-3-2022
end select

endif

end subroutine


!*******************************************SUBROUTINE********************************************************
subroutine cellcopymenu_handler(selection) !<<<<<<Tommi 23.9.2013 Cell copying menu
integer(kind=glint), value :: selection

if (cellid<=0) then
print*,"No cell selected! Please select a cell"
else
select case(selection)
case(1)
         print*,"The number of nodes in the chosen cell is:", cels(cellid)%nunodes
         tempcel%nunodes=cels(cellid)%nunodes
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         print*," Warning! Changing number of nodes in a a cell may cause problems."
         property=1 
         call menueditor_handler(9)
	 
case(2)
	 print*,"The fase of the chosen cell is", cels(cellid)%fase
         tempcel%fase=cels(cellid)%fase
         print*,"Now, please select a new cell where you want to copy this property and select paste property"
         property=2 
         call menueditor_handler(9)
case(3)
	 print*,"The number of node matrix in the chosen cell is:", cels(cellid)%nodela
         tempcel%nodela=cels(cellid)%nodela
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         print*," Warning! Changing number of nodes in a a cell may cause problems."
         property=3 
         call menueditor_handler(9)
case(4)
	 print*,"The cex of the chosen cell is", cels(cellid)%cex
         tempcel%cex=cels(cellid)%cex
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=4 
         call menueditor_handler(9)
case(5)
	 print*,"The cey of the chosen cell is", cels(cellid)%cey
          tempcel%cey=cels(cellid)%cey
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=5 
         call menueditor_handler(9)
case(6)
	 print*,"The cez of the chosen cell is", cels(cellid)%cez
         tempcel%cez=cels(cellid)%cez
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=6 
         call menueditor_handler(9)
case(7) 
	 print*,"The ctipus of the chosen cell is", cels(cellid)%ctipus
         tempcel%ctipus=cels(cellid)%ctipus
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=7 
         call menueditor_handler(9)
case(8)
	 print*,"The polx of the chosen cell is", cels(cellid)%polx
         tempcel%polx=cels(cellid)%polx
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=8
         call menueditor_handler(9)
case(9)
	 print*,"The poly of the chosen cell is", cels(cellid)%poly
  	 tempcel%poly=cels(cellid)%poly
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=9 
         call menueditor_handler(9)
case(10)
	 print*,"The polz of the chosen cell is", cels(cellid)%polz
         tempcel%polz=cels(cellid)%polz
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=10 
         call menueditor_handler(9)
case(11)
case(12)
	 print*,"The minsize_for_div of the chosen cell is", cels(cellid)%minsize_for_div
         tempcel%minsize_for_div=cels(cellid)%minsize_for_div
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         property=12 
         call menueditor_handler(9)
case(13)
         k=cels(cellid)%nunodes
         print*,"The cell has ",k," nodes"
         print*,"which are:"
         do i = 1,k
           print*,"Node ",i, "is", cels(cellid)%node(i)
         end do
         do i = 1,k
           tempcel%node(i)=cels(cellid)%node(i)
         end do
         print*,"Now, please select a new cell where you want to copy this property and select paste property."
         print*," Warning! Changing the nodes in a a cell may cause problems."
         property=13 
         call menueditor_handler(9)
case(14)
         if(property<=0) then
            print*,"No property of a cell selected! Please select a property"
         else
	 	 select case (property)
		 case(1)
	         colorselection=23
		 cels(cellid)%nunodes=tempcel%nunodes
		 print*,"Property changed"
		 case(2)
                 colorselection=23
		 cels(cellid)%fase=tempcel%fase
		 print*,"Property changed"
		 case(3) 
                 colorselection=23
		 cels(cellid)%nodela=tempcel%nodela
		 print*,"Property changed"
	 	 case(4)
	         colorselection=23
		 cels(cellid)%cex=tempcel%cex
		 print*,"Property changed"
		 case(5)
                 colorselection=23
		 cels(cellid)%cey=tempcel%cey
		 print*,"Property changed"
		 case(6) 
                 colorselection=23
		 cels(cellid)%cez=tempcel%cez
		 print*,"Property changed"
                 case(7)
	         colorselection=23
		 cels(cellid)%ctipus=tempcel%ctipus
		 print*,"Property changed"
		 case(8)
                 colorselection=23
		 cels(cellid)%polx=tempcel%polx
		 print*,"Property changed"
		 case(9) 
                 colorselection=23
		 cels(cellid)%poly=tempcel%poly
		 print*,"Property changed"
	 	 case(10)
	         colorselection=23
		 cels(cellid)%polz=tempcel%polz
		 print*,"Property changed"
		 case(11)
                 colorselection=23
	         print*,"Property changed"
	         case(12)
                 colorselection=23
		 cels(cellid)%minsize_for_div=tempcel%minsize_for_div
		 print*,"Property changed"
	         case(13)
                 colorselection=23
		 k=tempcel%nunodes
		 do i = 1,k
           		cels(cellid)%node(i)=tempcel%node(i)
        	 end do
		 print*,"Property changed"
		 end select
        endif
    end select	
endif

end subroutine

!*******************************************SUBROUTINE******************************************************** !!!!>>>>>>Tommi 23.9.2013


subroutine nodemenu_handler(selection) !<<<<<Tommi, menu for the node properties 11.9.2013
integer(kind=glint), value :: selection
real*8 :: temp

if (nodeindex<=0) then
print*,"No node selected! Please select a node"
else

select case(selection)
case(1)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%x
	 print*,"Please give new x for the selected node"
         read(*,*) temp
         node(nodeindex)%x=temp
case(2)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%y
	 print*,"Please give new y for the selected node"
         read(*,*) temp
         node(nodeindex)%y=temp
case(3)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%z
	 print*,"Please give new z for the selected node"
         read(*,*) temp
         node(nodeindex)%z=temp
case(5)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%eqd
	 print*,"Please give new req for the selected node"
         read(*,*) temp
         node(nodeindex)%eqd=temp
         nodeo(nodeindex)%eqd=temp
!case(6)
!         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%codel
!	 print*,"Please give new reqcel for the selected node"
!         read(*,*) temp
!         node(nodeindex)%codel=temp
case(6)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%add
	 print*,"Please give new da for the selected node"
         read(*,*) temp
         node(nodeindex)%add=temp
         nodeo(nodeindex)%add=temp
case(7)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%you
	 print*,"Please give new you for the selected node"
         read(*,*) temp
         node(nodeindex)%you=temp
         nodeo(nodeindex)%you=temp
case(8)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%adh
	 print*,"Please give new adh for the selected node"
         read(*,*) temp
         node(nodeindex)%adh=temp
         nodeo(nodeindex)%adh=temp
case(9)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%rep
	 print*,"Please give new rep for the selected node"
         read(*,*) temp
         node(nodeindex)%rep=temp
         nodeo(nodeindex)%rep=temp
case(10)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%rec
	 print*,"Please give new repcel for the selected node"
         read(*,*) temp
         node(nodeindex)%rec=temp
         nodeo(nodeindex)%rec=temp
case(11)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%erp
	 print*,"Please give new tor for the selected node"
         read(*,*) temp
         node(nodeindex)%erp=temp
         nodeo(nodeindex)%erp=temp
case(12)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%est
	 print*,"Please give new stor for the selected node"
         read(*,*) temp
         node(nodeindex)%est=temp
         nodeo(nodeindex)%est=temp
case(13)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%eqs
	 print*,"Please give new reqs for the selected node"
         read(*,*) temp
         node(nodeindex)%eqs=temp
         nodeo(nodeindex)%eqs=temp
case(14)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%hoo
	 print*,"Please give new ke for the selected node"
         read(*,*) temp
         node(nodeindex)%hoo=temp
         nodeo(nodeindex)%hoo=temp
case(15)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%mov
	 print*,"Please give new mo for the selected node"
         read(*,*) temp
         node(nodeindex)%mov=temp
         nodeo(nodeindex)%mov=temp
case(16)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%dmo
	 print*,"Please give new dmo for the selected node"
         read(*,*) temp
         node(nodeindex)%dmo=temp
         nodeo(nodeindex)%dmo=temp
case(17) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%orix
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%orix=temp
case(18) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%oriy
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%oriy=temp
case(19) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%oriz
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%oriz=temp
case(20) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%ecm
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%ecm=temp
case(21) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%cod
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%cod=temp
         nodeo(nodeindex)%cod=temp
case(22) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%grd
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%grd=temp
         nodeo(nodeindex)%grd=temp
case(23) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%pld
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%pld=temp
         nodeo(nodeindex)%pld=temp
case(24) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%vod
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%vod=temp
         nodeo(nodeindex)%vod=temp
case(25) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%dif
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%dif=temp
         nodeo(nodeindex)%dif=temp
case(26) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%kfi
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%kfi=temp
         nodeo(nodeindex)%kfi=temp
case(27) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%pla
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%pla=temp
         nodeo(nodeindex)%pla=temp
case(28) 
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%kvol
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%kvol=temp
         nodeo(nodeindex)%kvol=temp
case(29)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%tipus
	 print*,"Please give new tipus for the selected node"
         read(*,*) temp
         node(nodeindex)%tipus=temp
case(30)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%icel
	 print*,"Please give new icel for the selected node"
         read(*,*) temp
         node(nodeindex)%icel=temp
case(31)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%altre
	 print*,"Please give new altre for the selected node"
         read(*,*) temp
         node(nodeindex)%altre=temp
case(32)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%marge
	 print*,"Please give new marge for the selected node"
         read(*,*) temp
         node(nodeindex)%marge=temp
case(33)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%talone
	 print*,"Please give new marge for the selected node"
         read(*,*) temp
         node(nodeindex)%talone=temp
case(34)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%fix
	 print*,"Please give new marge for the selected node"
         read(*,*) temp
         node(nodeindex)%fix=temp
case(35)
         print*,"The ",nodeparams(selection)," of the node is", node(nodeindex)%ndiv !!>> HC 16-3-2022
	 print*,"Please give new marge for the selected node"
         read(*,*) temp
         node(nodeindex)%ndiv=temp   !!>> HC 16-3-2022
end select

print*,"Property of the node changed"
end if
end subroutine !>>>>>Tommi
!*******************************************SUBROUTINE********************************************************
subroutine cellmenu_handler(selection) !<<<<<Tommi, cell properties menu handler 11.9.2013
integer(kind=glint), value :: selection
real*8 :: temp
integer :: temp2, temp3
if (cellid<=0) then
print*,"No cell selected! Please select a cell"
else
select case(selection)
case(1)
         print*,"The number of nodes in the chosen cell is", cels(cellid)%nunodes
         print*,"Please give new number of nodes (nunodes) for the selected cell. Warning! If you change this the program may crash"
         read(*,*) temp
         cels(cellid)%nunodes=temp
         print*,"Property of the cell changed"
case(2)
	 print*,"The fase of the chosen cell is", cels(cellid)%fase
         print*,"Please give new fase for the selected cell."
         read(*,*) temp
         cels(cellid)%fase=temp
         print*,"Property of the cell changed"
case(3)
	 print*,"The size of the node matrix (nodela) in the chosen cell is", cels(cellid)%nodela
  	 print*,"Please give new number of total nodes (nodela) for the selected cell. Warning! If you change this the program may crash"
         read(*,*) temp
         cels(cellid)%nodela=temp
         print*,"Property of the cell changed"
case(4)
	 print*,"The cex of the chosen cell is", cels(cellid)%cex
         print*,"Please give new cex for the selected cell." 
         read(*,*) temp
         cels(cellid)%cex=temp
         print*,"Property of the cell changed"
case(5)
	 print*,"The cey of the chosen cell is", cels(cellid)%cey
         print*,"Please give new cey for the selected cell. "
         read(*,*) temp
         cels(cellid)%cey=temp
         print*,"Property of the cell changed"
case(6)
	 print*,"The cez of the chosen cell is", cels(cellid)%cez
         print*,"Please give new cez for the selected cell. "
         read(*,*) temp
         cels(cellid)%cez=temp
         print*,"Property of the cell changed"
case(7) 
	 print*,"The ctipus of the chosen cell is", cels(cellid)%ctipus
         print*,"Please give new ctipus for the selected cell." 
         read(*,*) temp
         cels(cellid)%ctipus=temp
         print*,"Property of the cell changed"
case(8)
	 print*,"The polx of the chosen cell is", cels(cellid)%polx
         print*,"Please give new polx for the selected cell."
         read(*,*) temp
         cels(cellid)%polx=temp
         print*,"Property of the cell changed"
case(9)
	 print*,"The poly of the chosen cell is", cels(cellid)%poly
  	 print*,"Please give new poly for the selected cell."
         read(*,*) temp
         cels(cellid)%cex=temp
         print*,"Property of the cell changed"
case(10)
	 print*,"The polz of the chosen cell is", cels(cellid)%polz
         print*,"Please give new polz for the selected cell."
         read (*,*) temp
         cels(cellid)%polz=temp
         print*,"Property of the cell changed"
case(11)
case(12)
	 print*,"The minsize_for_div of the chosen cell is", cels(cellid)%minsize_for_div
         print*,"Please give new minsize_for_div for the selected cell."
         read(*,*) temp
         cels(cellid)%minsize_for_div=temp
         print*,"Property of the cell changed"
case(13)
         k=cels(cellid)%nunodes
         print*,"The cell has ",k," nodes"
         print*,"which are:"
         do i = 1,k
           print*,"Node ",i, "is", cels(cellid)%node(i)
         end do
         print*,"Change which node? Give the number of node to be changed from 1 to", cels(cellid)%nunodes
         read(*,*) temp2
         print*, "Please give in the new index number of the node (in the nodes matrix) in that position."
         read(*,*) temp3
         cels(cellid)%node(temp2)=temp3
         print*,"Property of the cell changed"
end select
endif
end subroutine

!*******************************************SUBROUTINE******************************************************** !>>>>>Tommi

subroutine menus_handler(selection)
integer(kind=glint), value :: selection
character*300 nofi
character*300 crida

select case(selection)
case(1) ! save present time
  call writesnap
case(2) ! make snaps
  if (fsnap==1) then ; fsnap=0 ; else ; fsnap=1 ; end if
case(3) ! change snapshot frequency
  print *,"the snapshots are every",freqsnap,"iterations"
  print *,"every how many iterations do you want it now?"
  read (*,*) freqsnap
  fsnap=1
case(4) ! read
  print *,'name of the file, IT SHOULD BE WRITTEN BETWEEN ""s '
  read (*,*) nofi
  call readsnap(nofi)
case(5) !add label to output file
  print *,"print the text that will be added in the front of each output file"
  read(*,*) label
  flabel=1
case(6)
  if (ffinal==1) then ; ffinal=0 ; else ; ffinal=1 ; endif
!case(7)
!  if (fappend==1) then ; fappend=0 ; else ; fappend=1 ; endif
!case(8)
!  if (frappend==1) then ; frappend=0 ; else ; frappend=1 ; endif
!  print *,"name of the file"
!  read (*,*) nofi
!  call readsnap(nofi)
!case(9)
!  call readsnap(nofi)
case(10) !>>> 24-10-14
  if (fmovie==1) then ; fmovie=0 ; else ; fmovie=1 ; endif
  passed=1
  print *,""
  print *,"WARNING the window needs to be visible" 
  print *,""
  call system("mkdir images/"//trim(caa))
  !Is >>> 24-10-14
!case(11)
!  call system("mkdir images/")
!  call system("mkdir images/"//caa)
!  call savegif
case(12)
  if (fmovie==1) then ; fmovie=0 ; else ; fmovie=1 ; endif
  passed=1
!>>> Is 24-10-14
  call system("mkdir images/"//trim(caa))
  print *,"how often to take images (in iterations)?"
  read(*,*) freqsnap
  print *,"until which iteration?"
  read(*,*) fiite
  print *,""
  print *,"WARNING THE MAIN WINDOW NEEDS TO BE VISIBLE" 
  print *,""
!>>> Is 24-10-14
end select
end subroutine

subroutine menuse_handler(selection)
integer(kind=glint), value :: selection

select case(selection)
case(1)
  if (flag(14)==0) then ; mixi=mix ; call twodplot_tallx(mix) ; flag(14)=1 ;else ; flag(14)=0 ; endif
case(2)
  if (flag(14)==0) then ; mixi=mx ; call twodplot_tallx(mx) ; flag(14)=1 ;else ; flag(14)=0 ; endif
case(3)
  if (flag(15)==0) then ; miyi=miy ; call twodplot_tally(miy) ; flag(15)=1 ;else ; flag(15)=0 ; endif
case(4)
  if (flag(15)==0) then ; miyi=my ; call twodplot_tally(my) ; flag(15)=1 ;else ; flag(15)=0 ; endif
case(5)
  if (flag(16)==0) then ; mizi=miz ; call twodplot_tallz(miz) ; flag(16)=1 ;else ; flag(16)=0 ; endif
case(6)
  if (flag(16)==0) then ; mizi=mz ; call twodplot_tallz(mz) ; flag(16)=1 ;else ; flag(16)=0 ; endif
case(7)
  print *,"decoy is presently",decoy,"this is the node from which we calculate the field energy"
  read(*,*) decoy
case(8)
  print *,"the number of points in each direction per field is:",nq,"enter the new one"
  read (*,*) nq
case(9) 
  if (flag(17)==0) then ; flag(17)=1 ;else ; flag(17)=0 ; endif
case default
  if (flag(selection)==0) then ; flag(19:23)=0; flag(selection)=1 ; else ; flag(selection)=0 ; end if
end select
end subroutine


!******************************************************!>>>Miguel29-10-14 made subroutine
subroutine menuplot_handler(selection)                  
integer(kind=glint), value :: selection
integer:: i,nose,gese
integer,dimension(10)::single_gene
select case(selection)
case(1)
        lock=abs(lock-1)
        selection=0     
case(2)                      ! >>>Miguel8-10-14  made "case" 
        if(lock.eq.1)then
          print*,"Sorry, plotting option is disabled" 
          return
	endif  !>>>Miguel29-10-14        
        nodoss=0 ; kii=0 ; ki=0
     if((nki.ge.1).and.(nki.le.10))then ; nodoss(1:nki)=oopp(1:nki);nose=nki ; goto 1771 ; endif !>>>Miguel30-10-14
        if(nki.gt.10)then ; nodoss(1:10)=oopp(1:10) ; print*, 'too many nodes chosen(max=10)' ; goto 1771 ; endif !>>>Miguel30-10-14
        do while(ki.ne.-1)                                    ! delete the old one !!
          kii=kii+1         
          print*,"which node(s)? (-1 to skip select node)"        ! >>>Miguel8-10-14  
          read(*,*)ki
          if(kii.gt.10)then;write(*,*)'too many nodes chosen(max=10)';exit;endif  
          if(ki.gt.nd)then;write(*,*)'wrong input:there are only',nd,'nodes';kii=kii-1;cycle;endif 
          if((ki.eq.0).or.(ki.lt.-1))then;write(*,*)'wrong input: please choose a positive value';kii=kii-1;cycle;endif	                   
          nodoss(kii)=ki          
        end do
        nose=kii
        if(ki.eq.-1)then;nodoss(kii)=0;endif
1771    geness=0 ;kii=0 ; ki=0
        do while(ki.ne.-1)
          kii=kii+1
          print*,"which genes(s)? (-1 to skip select gene)"        ! >>>Miguel8-10-14  
          read(*,*)ki
          if(kii.gt.10)then;write(*,*)'too many genes chosen(max=10)';exit;endif  
          if(ki.gt.ng)then;write(*,*)'wrong input: there are only',ng,'genes';kii=kii-1;cycle;endif 	                   
          if((ki.eq.0).or.(ki.lt.-1))then;write(*,*)'wrong input: please choose a positive value';kii=kii-1;cycle;endif 
          geness(kii)=ki          
        end do
        if(ki.eq.-1)then;geness(kii)=0;endif
        gese=kii-1

        if(gese==1)then !one node selected, we plot all the genes in the same window !>>Miquel6-11-14
          call pintagnufor(nodoss,geness)    ! >>>Miguel8-10-14
        elseif(nose==1)then !one node selected, we plot all the genes in the same window !>>Miquel6-11-14
          call pintagnufor(nodoss,geness)    ! >>>Miguel8-10-14
        else
          do i=1,kii
            single_gene=0
            single_gene(1)=geness(i)
            call pintagnufor(nodoss,single_gene)
          end do
        end if


end select
end subroutine



!*******************************************SUBROUTINE********************************************************
subroutine contour_color_handler(selection)
integer(kind=glint), value :: selection

contour_color = selection
call draw_func

end subroutine contour_color_handler

!*************************************************************************************************************
subroutine twodplot_tallx(slice)  !makes a raster plot of energy in the slice 
  real*8 slice
  integer i,j,k
  real*8 a,b,dy,dz
  real*8 dmy,dmz
  real*8 oldx,oldy,oldz

  if (allocated(theslice)) deallocate(theslice)
  allocate(theslice(0:nq,0:nq))
  dmy=my-miy ; dmz=mz-miz

  oldx=node(decoy)%x ; oldy=node(decoy)%y ; oldz=node(decoy)%z

  !here we calculate the slice
  node(decoy)%x=slice      
  do i=0,nq
    do j=0,nq
      dy=dmy/nq
      dz=dmz/nq
      a=miy+i*dy ; b=miz+j*dz
      node(decoy)%y=a ; node(decoy)%z=b
      call energia(decoy)      
      if (flag(19)==1) then ;  theslice(i,j)=node(decoy)%e     ;cycle ;end if
      if (flag(20)==1) then ;  theslice(i,j)=eadh(decoy)  ;cycle ;end if
      if (flag(21)==1) then ;  theslice(i,j)=eyou(decoy)  ;cycle ;end if
      if (flag(22)==1) then ;  theslice(i,j)=erep(decoy)  ;cycle ;end if
      if (flag(23)==1) then ;  theslice(i,j)=erepcel(decoy)  ;cycle ;end if
    end do
  end do  

  node(decoy)%x=oldx ; node(decoy)%y=oldy ; node(decoy)%z=oldz

end subroutine
!*************************************************************************************************************
subroutine twodplot_tally(slice)  !makes a raster plot of energy in the slice 
  real*8 slice
  integer i,j,k
  real*8 a,b,dx,dz
  real*8 dmx,dmz
  real*8 oldx,oldy,oldz

  if (allocated(theslice)) deallocate(theslice)
  allocate(theslice(0:nq,0:nq))
  dmx=mx-mix ; dmz=mz-miz

  oldx=node(decoy)%x ; oldy=node(decoy)%y ; oldz=node(decoy)%z

  !here we calculate the slice
  node(decoy)%y=slice      
  do i=0,nq
    do j=0,nq
      dx=dmx/nq
      dz=dmz/nq
      a=mix+i*dx ; b=miz+j*dz
      node(decoy)%x=a ; node(decoy)%z=b
      call energia(decoy)      
      if (flag(19)==1) then ;  theslice(i,j)=node(decoy)%e     ;cycle ;end if
      if (flag(20)==1) then ;  theslice(i,j)=eadh(decoy)  ;cycle ;end if
      if (flag(21)==1) then ;  theslice(i,j)=eyou(decoy)  ;cycle ;end if
      if (flag(22)==1) then ;  theslice(i,j)=erep(decoy)  ;cycle ;end if
      if (flag(23)==1) then ;  theslice(i,j)=erepcel(decoy)  ;cycle ;end if
    end do
  end do  

  node(decoy)%x=oldx ; node(decoy)%y=oldy ; node(decoy)%z=oldz

end subroutine

!*************************************************************************************************************
subroutine twodplot_tallz(slice)  !makes a raster plot of energy in the slice 
  real*8 slice
  integer i,j,k
  real*8 a,b,dx,dy
  real*8 dmx,dmy
  real*8 oldx,oldy,oldz

  if (allocated(theslice)) deallocate(theslice)
  allocate(theslice(0:nq,0:nq))
  dmx=mx-mix ; dmy=my-miy

  oldx=node(decoy)%x ; oldy=node(decoy)%y ; oldz=node(decoy)%z

  !here we calculate the slice
  node(decoy)%z=slice      
  do i=0,nq
    do j=0,nq
      dx=dmx/nq
      dy=dmy/nq
      a=mix+i*dx ; b=miy+j*dy
      node(decoy)%x=a ; node(decoy)%y=b
      call energia(decoy)      
      if (flag(19)==1) then ;  theslice(i,j)=node(decoy)%e     ;cycle ;end if
      if (flag(20)==1) then ;  theslice(i,j)=eadh(decoy)  ;cycle ;end if
      if (flag(21)==1) then ;  theslice(i,j)=eyou(decoy)  ;cycle ;end if
      if (flag(22)==1) then ;  theslice(i,j)=erep(decoy)  ;cycle ;end if
      if (flag(23)==1) then ;  theslice(i,j)=erepcel(decoy)  ;cycle ;end if
    end do
  end do  

  node(decoy)%x=oldx ; node(decoy)%y=oldy ; node(decoy)%z=oldz

end subroutine

!*******************************************SUBROUTINE********************************************************

subroutine surface_color_handler(selection)
integer(kind=glint), value :: selection

surface_color = selection
call draw_func

end subroutine surface_color_handler

!***************************************************************************************

subroutine make_menu(submenuid)
integer, intent(in) :: submenuid
integer submenuidd,menuu,menuu2,menud,menut,menuq,menus,menuse, menueditor, planemenu,nodemenu,cellmenu,nodecopymenu,cellcopymenu, &
          & nodepastemenu,pastemenu,sel_menu,cursor_menu,color_menu,menuplot,& !Tommi 27.9.2013 !>>Miquel29-10-14
          cmenuu,cmenuu2,amenuu,amenuu2,smenuu,smenuu2,palette_menu,arrow_menu,sphere_menu,move_menu
integer :: menuid, param_id, contour_color_menu, surface_color_menu
integer :: nu !Tommi 1.8.2013
 character(70) :: genetitle !>>>Miguel30-10-14

!contour_color_menu = glutCreateMenu(contour_color_handler)
!call glutAddMenuEntry("black"//char(0),black_contour)
!call glutAddMenuEntry("contour value"//char(0),rainbow_contour)

!surface_color_menu = glutCreateMenu(surface_color_handler)
!call glutAddMenuEntry("red"//char(0),red_surface)
!call glutAddMenuEntry("white"//char(0),white_surface)
!call glutAddMenuEntry("surface height"//char(0),rainbow_surface)

!param_id = glutCreateMenu(param_handler)
!call glutAddMenuEntry("reset to initial parameters"//char(0),reset_params)
!call glutAddMenuEntry("point size"//char(0),p_size)

nu=1 !Tommi 1.8.2013

!!!! COLOR MENUS  !>>Miquel30-10-14
  cmenuu = glutCreateMenu(cmenuu_handler)
  do nu=1,27 !nparam_per_node !>>>>>Tommi 1.8.2013
    if (nu==4) cycle
    if (nu>16.and.nu<19) cycle
    call glutAddMenuEntry(nodeparams(nu)//char(0),nu)
  end do
  
  cmenuu2 = glutCreateMenu(cmenuu2_handler)
  call glutAddMenuEntry("node type"//char(0),29)
  call glutAddMenuEntry("cells"//char(0),30)
  call glutAddMenuEntry("cell nucleus as blue"//char(0),32)
  call glutAddMenuEntry("nodes fixed as yellow"//char(0),34)
  call glutAddMenuEntry("Number of divisions"//char(0),35)
  call glutAddMenuEntry("sq distance center"//char(0),36)                              !!>> HC 15-12-2021

  nu=nparam_per_node+1
  call glutAddMenuEntry("as cell cycle"//char(0),nu+1)
  call glutAddMenuEntry("as dx"//char(0),nu+2)
  call glutAddMenuEntry("as dy"//char(0),nu+3)
  call glutAddMenuEntry("as dz"//char(0),nu+4)
  call glutAddMenuEntry("as dtotal"//char(0),nu+5)
  call glutAddMenuEntry("as boxes"//char(0),nu+6)
  call glutAddMenuEntry("as node index"//char(0),nu+7)!>>>>>Tommi 5.8.2013
  call glutAddMenuEntry("cell lineage"//char(0),nu+8)                              !!>> HC 15-12-2021
  do i=1,ng						!>>>Miguel30-10-14
    write(genetitle,*)'Amount of regulatory molecule',i  !>>>Miguel30-10-14
    genetitle=adjustl(genetitle)				!>>>Miguel30-10-14
    call glutAddMenuEntry(genetitle(1:44)//char(0),nu+8+i) !>>>Miguel30-10-14      !!>> HC 15-12-2021
  end do                                                 !>>>Miguel30-10-14 
!!!!!!!!!!!

!!!! ARROW MENUS  !>>Miquel30-10-14
  amenuu = glutCreateMenu(amenuu_handler)
  do nu=1,27 !nparam_per_node !>>>>>Tommi 1.8.2013
    if (nu==4) cycle
    if (nu>16.and.nu<19) cycle
    call glutAddMenuEntry(nodeparams(nu)//char(0),nu)
  end do

  amenuu2 = glutCreateMenu(amenuu2_handler)
  call glutAddMenuEntry("node type"//char(0),29)
  call glutAddMenuEntry("cells"//char(0),30)
  call glutAddMenuEntry("cell nucleus as blue"//char(0),32)
  call glutAddMenuEntry("nodes fixed as yellow"//char(0),34)

  nu=nparam_per_node
  call glutAddMenuEntry("as cell cycle"//char(0),nu+1)
  call glutAddMenuEntry("as dx"//char(0),nu+2)
  call glutAddMenuEntry("as dy"//char(0),nu+3)
  call glutAddMenuEntry("as dz"//char(0),nu+4)
  call glutAddMenuEntry("as dtotal"//char(0),nu+5)
  call glutAddMenuEntry("as boxes"//char(0),nu+6)
  call glutAddMenuEntry("as node index"//char(0),nu+7)!>>>>>Tommi 5.8.2013
  do i=1,ng						!>>>Miguel30-10-14
    write(genetitle,*)'Amount of regulatory molecule',i  !>>>Miguel30-10-14
    genetitle=adjustl(genetitle)				!>>>Miguel30-10-14
    call glutAddMenuEntry(genetitle(1:43)//char(0),nu+7+i) !>>>Miguel30-10-14
  end do                                                 !>>>Miguel30-10-14 
!!!!!!!!!!!!!!

!!!! TRANSPARENT SPHERE MENUS  !>>Miquel30-10-14
  smenuu = glutCreateMenu(smenuu_handler)
  do nu=1,27 !nparam_per_node !>>>>>Tommi 1.8.2013
    if (nu==4) cycle
    if (nu>16.and.nu<19) cycle
    call glutAddMenuEntry(nodeparams(nu)//char(0),nu)
  end do

  smenuu2 = glutCreateMenu(smenuu2_handler)
  call glutAddMenuEntry("node type"//char(0),29)
  call glutAddMenuEntry("cell"//char(0),30)
  call glutAddMenuEntry("cell nucleus as blue"//char(0),32)
  call glutAddMenuEntry("nodes fixed as yellow"//char(0),34)

  nu=nparam_per_node
  call glutAddMenuEntry("as cell cycle"//char(0),nu+1)
  call glutAddMenuEntry("as dx"//char(0),nu+2)
  call glutAddMenuEntry("as dy"//char(0),nu+3)
  call glutAddMenuEntry("as dz"//char(0),nu+4)
  call glutAddMenuEntry("as dtotal"//char(0),nu+5)
  call glutAddMenuEntry("as boxes"//char(0),nu+6)
  call glutAddMenuEntry("as node index"//char(0),nu+7)!>>>>>Tommi 5.8.2013
  do i=1,ng						!>>>Miguel30-10-14
    write(genetitle,*)'Amount of regulatory molecule',i  !>>>Miguel30-10-14
    genetitle=adjustl(genetitle)				!>>>Miguel30-10-14
    call glutAddMenuEntry(genetitle(1:43)//char(0),nu+7+i) !>>>Miguel30-10-14
  end do                                                 !>>>Miguel30-10-14 
!!!!!!!!!!!!!!!


!PALETTE MENU
 palette_menu=glutCreateMenu(palette_menu_handler)
call glutAddMenuEntry("multicolor palette"//char(0),1) !>>Miquel29-10-14
call glutAddMenuEntry("yellow-red palette"//char(0),2)!>>Miquel29-10-14
call glutAddMenuEntry("blue palette"//char(0),3)!>>Miquel29-10-14
call glutAddMenuEntry("green palette"//char(0),4)!>>Miquel29-10-14

!COLOR OPTIONS MENU
 color_menu = glutCreateMenu(color_menu_handler)    !>>Miquel29-10-14
call glutAddSubMenu("change color palette"//char(0),palette_menu) !>>Miquel29-10-14
call glutAddMenuEntry("select color max. and min. from terminal"//char(0),2) !>>Miquel29-10-14
call glutAddMenuEntry("select color min. with left mouse button"//char(0),3)!>>Miquel29-10-14
call glutAddMenuEntry("select color max. with left mouse button"//char(0),4)!>>Miquel29-10-14
call glutAddMenuEntry("select color min. with middle mouse button"//char(0),5)!>>Miquel29-10-14
call glutAddMenuEntry("select color max. with middle mouse button"//char(0),6)!>>Miquel29-10-14
!!!!!!!!!!!!!!!

!ARROW OPTIONS MENU
 arrow_menu = glutCreateMenu(arrow_menu_handler)    !>>Miquel29-10-14
call glutAddMenuEntry("change arrow scale with left mouse button"//char(0),1)
call glutAddMenuEntry("select arrow max. and min. from terminal"//char(0),2) !>>Miquel29-10-14
call glutAddMenuEntry("select arrow min. with left mouse button"//char(0),3)!>>Miquel29-10-14
call glutAddMenuEntry("select arrow max. with left mouse button"//char(0),4)!>>Miquel29-10-14
call glutAddMenuEntry("select arrow min. with middle mouse button"//char(0),5)!>>Miquel29-10-14
call glutAddMenuEntry("select arrow max. with middle mouse button"//char(0),6)!>>Miquel29-10-14
call glutAddMenuEntry("disable arrows"//char(0),7)
!!!!!!!!!!!!!!!


!SPHERE OPTIONS MENU
 sphere_menu = glutCreateMenu(sphere_menu_handler)    !>>Miquel29-10-14
call glutAddMenuEntry("change sphere scale with left mouse button"//char(0),1)
call glutAddMenuEntry("select sphere max. and min. from terminal"//char(0),2) !>>Miquel29-10-14
call glutAddMenuEntry("select sphere min. with left mouse button"//char(0),3)!>>Miquel29-10-14
call glutAddMenuEntry("select sphere max. with left mouse button"//char(0),4)!>>Miquel29-10-14
call glutAddMenuEntry("select sphere min. with middle mouse button"//char(0),5)!>>Miquel29-10-14
call glutAddMenuEntry("select sphere max. with middle mouse button"//char(0),6)!>>Miquel29-10-14
call glutAddMenuEntry("disable spheres"//char(0),7)
!!!!!!!!!!!!!!!


menuu=glutCreateMenu(menuu_handler)
call glutAddSubMenu("with colors"//char(0),cmenuu)
call glutAddSubMenu("with arrows (only epithelia)"//char(0),amenuu)
call glutAddSubMenu("with semitransparent spheres"//char(0),smenuu)
call glutAddSubMenu("color options"//char(0),color_menu)
call glutAddSubMenu("arrow options"//char(0),arrow_menu)
call glutAddSubMenu("sphere options"//char(0),sphere_menu)

!call glutAddMenuEntry("change arrow scale"//char(0),3)
!call glutAddMenuEntry("change sphere scale"//char(0),4)
!call glutAddMenuEntry("dissable arrows"//char(0),1)
!call glutAddMenuEntry("dissable transparent spheres"//char(0),2)



menuu2=glutCreateMenu(menuu2_handler)
call glutAddSubMenu("with colors"//char(0),cmenuu2)
call glutAddSubMenu("with arrows (only epithelia)"//char(0),amenuu2)
call glutAddSubMenu("with semitransparent spheres"//char(0),smenuu2)
call glutAddSubMenu("color options"//char(0),color_menu)
call glutAddSubMenu("arrow options"//char(0),arrow_menu)
call glutAddSubMenu("sphere options"//char(0),sphere_menu)

!call glutAddMenuEntry("change arrow scale"//char(0),3)
!call glutAddMenuEntry("change sphere scale"//char(0),4)
!call glutAddMenuEntry("dissable arrows"//char(0),1)
!call glutAddMenuEntry("dissable transparent spheres"//char(0),2)



menud = glutCreateMenu(menud_handler)
call glutAddMenuEntry("epithelial springs"//char(0),1)
call glutAddMenuEntry("box grid"//char(0),2)
!call glutAddMenuEntry("polygons"//char(0),3)
!call glutAddMenuEntry("full polygons"//char(0),33)
!call glutAddMenuEntry("polygonal sides"//char(0),4)
!call glutAddMenuEntry("polyhedral cells mesenchyme"//char(0),34) ! miguel 25-9-13
!call glutAddMenuEntry("eggshell"//char(0),35) ! miguel4-11-13
call glutAddMenuEntry("connexions between cells"//char(0),12)
call glutAddMenuEntry("connexions between nodes"//char(0),13)
call glutAddMenuEntry("draw cell contour"//char(0),37)
call glutAddMenuEntry("draw intercellular contour"//char(0),38)
call glutAddMenuEntry("balls as radius=EQD"//char(0),5)
call glutAddMenuEntry("balls as radius=ADD"//char(0),6)
call glutAddMenuEntry("small balls"//char(0),11)
call glutAddMenuEntry("no balls"//char(0),8)
call glutAddMenuEntry("cylinders"//char(0),22)
call glutAddMenuEntry("only upper balls"//char(0),9)
call glutAddMenuEntry("only lower balls"//char(0),10)
call glutAddMenuEntry("display box boundaries"//char(0),18)
!call glutAddMenuEntry("dynamic/fixed display box"//char(0),41)
call glutAddMenuEntry("polarization vectors"//char(0),19)
!call glutAddMenuEntry("Node property as transparent balls"//char(0),28)
!call glutAddMenuEntry("Node property as arrows"//char(0),29)
!call glutAddMenuEntry("Change scale of node property"//char(0),33)
call glutAddMenuEntry("centroids"//char(0),21) ! miguel 30-5-13
call glutAddMenuEntry("fixed nodes"//char(0),36) !>> Miquel 15-1-14
call glutAddMenuEntry("show displacement of nodes from origin"//char(0),39) !>> Miquel 15-1-14
call glutAddMenuEntry("movement unitary vectors"//char(0),23)     !>>>>Miquel 30-8-13
call glutAddMenuEntry("force component: repulson-adhesion"//char(0),24)  !
call glutAddMenuEntry("force component: epi. surface tension lateral"//char(0),25)        !
call glutAddMenuEntry("force component: epi. surface tension apical/basal"//char(0),26)       !
call glutAddMenuEntry("movement module vectors"//char(0),27)      !
call glutAddMenuEntry("do not show/show epithelium"//char(0),30)      !
call glutAddMenuEntry("do not show/show mesenchyme"//char(0),31)      !
call glutAddMenuEntry("do not show/show extracellular matrix"//char(0),32)      !  !
call glutAddMenuEntry("Ellongated cells as ellipses"//char(0),33)               !!>> HC 8-4-2022

!menut = glutCreateMenu(menut_handler)
!call glutAddMenuEntry("cell division"//char(0),1)
!call glutAddMenuEntry("cell growth"//char(0),2)
!call glutAddMenuEntry("invaginate cell"//char(0),3)
!call glutAddMenuEntry("add node"//char(0),4)
!call glutAddMenuEntry("add ECM"//char(0),5)

menuq = glutCreateMenu(menuq_handler)
call glutAddMenuEntry("go 1 iterations back"//char(0),1)
call glutAddMenuEntry("go 10 iterations back"//char(0),10)
call glutAddMenuEntry("go x iterations back"//char(0),-1)

menus = glutCreateMenu(menus_handler)
call glutAddMenuEntry("save present time"//char(0),1)
!call glutAddMenuEntry("save only the last time file"//char(0),6)
!call glutAddMenuEntry("save: all times in the same file: append"//char(0),7)
call glutAddMenuEntry("save snaps periodically"//char(0),2)
call glutAddMenuEntry("change frequency of snapshots"//char(0),3)
call glutAddMenuEntry("read from file"//char(0),4)
!call glutAddMenuEntry("read from file successively: append"//char(0),8)
!call glutAddMenuEntry("read again: append"//char(0),9)
call glutAddMenuEntry("add label to the name of the output file"//char(0),5)
!call glutAddMenuEntry("save into gif or alike (it can be changed)"//char(0),11)
call glutAddMenuEntry("save images automatically (the window must be open and uncovered)"//char(0),10) !>>> is 24-10-14
call glutAddMenuEntry("save movie, check terminal (the window must be open and uncovered)"//char(0),12)

!menuse = glutCreateMenu(menuse_handler)
!call glutAddMenuEntry("make 2D plot in the minimal x plane"//char(0),1)
!call glutAddMenuEntry("make 2D plot in the maximal x plane"//char(0),2)
!call glutAddMenuEntry("make 2D plot in the minimal y plane"//char(0),3)
!call glutAddMenuEntry("make 2D plot in the maximal y plane"//char(0),4)
!call glutAddMenuEntry("make 2D plot in the minimal z plane"//char(0),5)
!call glutAddMenuEntry("make 2D plot in the maximal z plane"//char(0),6)
!call glutAddMenuEntry("chose the decoy node: default is node 1"//char(0),7)
!call glutAddMenuEntry("sampling points per interval"//char(0),8)
!call glutAddMenuEntry("make plot instead of raster"//char(0),9)
!call glutAddMenuEntry("show energy"//char(0),19)
!call glutAddMenuEntry("show adhesion energy"//char(0),20)
!call glutAddMenuEntry("show you energy"//char(0),21)
!call glutAddMenuEntry("show rep energy"//char(0),22)
!call glutAddMenuEntry("show repcel energy"//char(0),23)

!planemenu = glutCreateMenu(planemenu_handler) !>>>>Tommi 11.9.2013
!call glutAddMenuEntry("Maximal x plane"//char(0),2)!The menu to select the plane of adding/selecting
!call glutAddMenuEntry("Minimal x plane"//char(0),3)
!call glutAddMenuEntry("Maximal y plane"//char(0),4)
!call glutAddMenuEntry("Minimal y plane"//char(0),1)
!call glutAddMenuEntry("Maximal z plane"//char(0),5)
!call glutAddMenuEntry("Minimal z plane"//char(0),6)!<<<<<Tommi

cellcopymenu = glutCreateMenu(cellcopymenu_handler)  !>>>>>Tommi 23.9.2013
call glutAddMenuEntry("nunodes"//char(0),1)
call glutAddMenuEntry("fase"//char(0),2)!The menu to change the properties of a selected cell
call glutAddMenuEntry("nodela"//char(0),3)
call glutAddMenuEntry("cex"//char(0),4)
call glutAddMenuEntry("cey"//char(0),5)
call glutAddMenuEntry("cez"//char(0),6)
call glutAddMenuEntry("ctipus"//char(0),7)
call glutAddMenuEntry("polx"//char(0),8)
call glutAddMenuEntry("poly"//char(0),9)
call glutAddMenuEntry("polz"//char(0),10) 
call glutAddMenuEntry("minsize_for_div"//char(0),12) 
call glutAddMenuEntry("Copy nodes in the selected cell"//char(0),13) 
call glutAddMenuEntry("Paste the property to the selected cell"//char(0),14)

nodecopymenu= glutCreateMenu(nodecopymenu_handler)
do i=1,nparam_per_node
  call glutAddMenuEntry(nodeparams(i)(1:10)//char(0),i)  
end do

nodepastemenu= glutCreateMenu(nodepastemenu_handler)
do i=1,nparam_per_node
  call glutAddMenuEntry(nodeparams(i)(1:10)//char(0),i)  
end do


nodemenu = glutCreateMenu(nodemenu_handler) !>>>>Tommi 11.9.2013
do i=1,nparam_per_node
  call glutAddMenuEntry(nodeparams(i)(1:10)//char(0),i)  !Is 11-9-13 
end do

cellmenu = glutCreateMenu(cellmenu_handler)  !>>>>>Tommi 11.9.2013
call glutAddMenuEntry("nunodes"//char(0),1)
call glutAddMenuEntry("fase"//char(0),2)!The menu to change the properties of a selected cell
call glutAddMenuEntry("nodela"//char(0),3)
call glutAddMenuEntry("cex"//char(0),4)
call glutAddMenuEntry("cey"//char(0),5)
call glutAddMenuEntry("cez"//char(0),6)
call glutAddMenuEntry("ctipus"//char(0),7)
call glutAddMenuEntry("polx"//char(0),8)
call glutAddMenuEntry("poly"//char(0),9)
call glutAddMenuEntry("polz"//char(0),10) 
call glutAddMenuEntry("minsize_for_div"//char(0),12) 
call glutAddMenuEntry("Change nodes in the cell"//char(0),13) !<<<<<<Tommi

menueditor = glutCreateMenu(menueditor_handler) !>>>>>Tommi !The editor main menu 27.9.2013
!call glutAddMenuEntry("Move the cursor in the x-y plane"//char(0),2) 
!call glutAddMenuEntry("Move the cursor in z"//char(0),15) 
!call glutAddMenuEntry("Select a node"//char(0),6)
!call glutAddMenuEntry("Select a cell"//char(0),9)
call glutAddMenuEntry("Add basic epithelial cell"//char(0),8)
call glutAddMenuEntry("Add basic mesenchymal cell"//char(0),12)
call glutAddMenuEntry("Add node to a selected cell"//char(0),1)
call glutAddMenuEntry("Add extracellular node"//char(0),13)
call glutAddMenuEntry("Paste a node from selection"//char(0),7)
call glutAddMenuEntry("Paste a cell from selection"//char(0),16)
call glutAddMenuEntry("Delete the selected node"//char(0),10)
call glutAddMenuEntry("Delete the selected cell"//char(0),11)
call glutAddSubMenu("Change properties of the selected node"//char(0),nodemenu)
call glutAddSubMenu("Change properties of the selected cell"//char(0),cellmenu)
call glutAddMenuEntry("Change gene expression in the selected node"//char(0),17)
call glutAddMenuEntry("Choose a node in which to Paste properties of the selected node"//char(0),14)
call glutAddSubMenu("Paste a property into it"//char(0),nodepastemenu)
call glutAddMenuEntry("STOP the editor and go back to running mode"//char(0),3)

sel_menu = glutCreateMenu(sel_menu_handler)    !>>Miquel24-10-14
call glutAddMenuEntry("select node by index"//char(0),1)!>>Miquel24-10-14
call glutAddMenuEntry("select node with cursor"//char(0),2)!>>Miquel24-10-14
call glutAddMenuEntry("select cell by index"//char(0),3)!>>Miquel24-10-14
call glutAddMenuEntry("select cell with cursor"//char(0),4)!>>Miquel24-10-14
call glutAddMenuEntry("undo selections"//char(0),5)!>>Miquel24-10-14


 cursor_menu = glutCreateMenu(cursor_menu_handler)    !>>Miquel24-10-14
call glutAddMenuEntry("cursor ON/OFF"//char(0),1)!>>Miquel24-10-14
call glutAddMenuEntry("set left mouse-button for x-y move"//char(0),2)!>>Miquel24-10-14
call glutAddMenuEntry("set left mouse-button for z move"//char(0),3)!>>Miquel24-10-14
call glutAddMenuEntry("set middle mouse-button for z move"//char(0),4)!>>Miquel24-10-14



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
menuplot = glutCreateMenu(menuplot_handler)                        !>>>Miguel29-10-14
  call glutAddMenuEntry("plot [gen] vs time in a node"//char(0),2) !>>>Miguel8-10-14
  call glutAddMenuEntry("Temporal plots: ON/OFF"//char(0),1)        !>>>Miguel8-10-14
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 move_menu = glutCreateMenu(move_menu_handler)    !>>Miquel31-10-14
call glutAddMenuEntry("move node from terminal"//char(0),1)!>>Miquel31-10-14
call glutAddMenuEntry("move node with left mouse button"//char(0),2)!>>Miquel31-10-14
call glutAddMenuEntry("move cell from terminal"//char(0),3)!>>Miquel31-10-14
call glutAddMenuEntry("move cell with left mouse button"//char(0),4)!>>Miquel31-10-14
call glutAddMenuEntry("Stop movement (for cursor movement)"//char(0),5)!>>Miquel31-10-14


menuid = glutCreateMenu(menu_handler)
call glutAddSubMenu("Basic View Modifier"//char(0),submenuid)
call glutAddSubMenu("show node mechanical properties"//char(0),menuu)
call glutAddSubMenu("show other node properties"//char(0),menuu2)
call glutAddSubMenu("what to draw"//char(0),menud)
call glutAddSubMenu("sections"//char(0),menut)
call glutAddSubMenu("input output"//char(0),menus)
call glutAddSubMenu("2D plots"//char(0),menuse)
call glutAddSubMenu("Editor"//char(0),menueditor) !Tommi 5.8.2013
call glutAddSubMenu("selection menu"//char(0),sel_menu)  !>>Miquel24-10-14
call glutAddSubMenu("cursor menu"//char(0),cursor_menu)  !>>Miquel24-10-14
call glutAddSubMenu("gene plotting options"//char(0),menuplot) !>>>Miguel29-10-14
call glutAddSubMenu("move menu"//char(0),move_menu) !>>>Miguel29-10-14

!call glutAddMenuEntry("select node"//char(0),1)      !X>>Miquel24-10-14
!call glutAddMenuEntry("select cell"//char(0),13)     !X>>Miquel24-10-14
!call glutAddMenuEntry("select color max. and min."//char(0),15) !>>Miquel 7-7-14
!call glutAddMenuEntry("change color palette"//char(0),16) !>>Miquel 7-7-14
!call glutAddMenuEntry("move node"//char(0),2)
call glutAddMenuEntry("1 iteration"//char(0),3)
call glutAddMenuEntry("100 iterations"//char(0),4)
call glutAddMenuEntry("1000 iterations"//char(0),5)
call glutAddMenuEntry("50000 iterations"//char(0),6)
call glutAddMenuEntry("x iterations"//char(0),7)
call glutAddMenuEntry("fitness"//char(0),15)
call glutAddMenuEntry("quit"//char(0),12)
!call glutAddMenuEntry("plot [gen] vs time in a node"//char(0),18) !>>>Miguel8-10-14

call glutAttachMenu(GLUT_RIGHT_BUTTON)

end subroutine make_menu

end module function_plotter
