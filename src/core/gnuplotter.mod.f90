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

module gnuplotter
implicit none
contains
subroutine  f2gp(nodoss,geness,imm)

integer :: n1,n2,j,ii,jj,im,ki,kii,rat                                     
integer :: nodoss(10),geness(10)
integer :: plot_type                                 ! 1 for linear plot, 2 for log plot, 3 for log-log plot
 character(len=20) :: xlabel,ylabel                  ! plot axis labels and title
 character(len=40)  ::title1,title2,arxiv  
!---------------

integer :: i,imm                                     ! number of lines i the plot points
integer :: ret

!---------------




xlabel='TIME'
ylabel='GENE CONCENTRATION'
plot_type=1                                         ! 1 for linear plot, 2 for log plot, 3 for log-log plot

!!!!!!!! it comes back if gnuplot is not installed
rat=0
call checkgnuplot(rat)
if(rat.eq.1)then
 write(*,*)'sorry, gene concentration can not be plotted because'
 write(*,*)'gnuplot is not installed'
 return
end if


! create gnuplot command file

OPEN(10,ACCESS='SEQUENTIAL',FILE='gp.txt')
write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
write(10,*) 'set autoscale'
if (plot_type==2) write(10,*) 'set log y'
if (plot_type==3) then
   write(10,*) 'set log x' ;  write(10,*) 'set log y'
endif
write(10,*) 'set style line 1  lt 2 lw 1 pt 4 ps 2'

jj=0 
do i=1,10
  ki=nodoss(i)   
  if(ki.ne.0)then
    do j=1,10
      kii=geness(j)
      if(kii.ne.0)then        
        jj=jj+1
        write(arxiv,*)"N",ki,"G",kii
        do im=1,40
          if (arxiv(im:im)==" ")then; arxiv(im:im)="_";end if    
          arxiv(im:im)=arxiv(im:im)
        end do
        open(11111+jj,file=arxiv,status='unknown',action='read') 
        write(title1,*)'node',ki,'gen',kii             
        
        if(jj.lt.imm)then  
          if(jj.eq.1)then
            write(10,*) 'plot "'//arxiv//'" using 1:2 with lines title "'//TRIM(title1)//'",\'
          else
            write(10,*) '"'//arxiv//'" using 1:2 with lines title "'//TRIM(title1)//'",\'
          end if
        else
          if(imm.gt.1)then
            write(10,*) '"'//arxiv//'" using 1:2 with lines title "'//TRIM(title1)//'"'   
          else  ! only one gene in a node has been selected
            write(10,*) 'plot "'//arxiv//'" using 1:2 with lines title "'//TRIM(title1)//'"'
          end if
        end if     
        close(11111+jj)
      end if
    end do
  end if
end do

 CLOSE(10,STATUS='KEEP')

! plot curve with gnuplot and cleanup files
ret=SYSTEM('gnuplot -persist gp.txt')
ret=SYSTEM('rm gp.txt')


end subroutine

!************************************************************************************************

subroutine checkgnuplot(rat)
 integer :: rat,ret,i
 character(len=10) :: output

 open(9999,file='check',status='unknown')
 ret=SYSTEM('dpkg -l | grep gnuplot > check')
 open(9959,file='check1',status='unknown')
 ret=SYSTEM('dpkg -l | du -bsh check > check1')
 read(9959,'(A)')output
 if(output(1:1).eq.'0')then
    rat=1
 end if
 ret=SYSTEM('rm check*')
! ret=SYSTEM('rm check1*')
end subroutine checkgnuplot

end module

