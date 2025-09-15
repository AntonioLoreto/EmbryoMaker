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




module automaticon
use general
use model
use io

contains

!**************************************************

subroutine auto
integer irf
integer iit_step,inu_it
character*10 cf
character*400 cx
integer*4 pid

iit_step=10
inu_it=100
call get_its(iit_step,inu_it)

do irf=1,inu_it
  call iteracio(iit_step)
  call writesnap
end do

call exit(231)

end subroutine

end module
