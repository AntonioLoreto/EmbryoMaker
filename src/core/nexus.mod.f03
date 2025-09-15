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




module nexus

  use general
  use genetic
  use growth  !>>>> Miquel 17-6-13
  use death    !>>>> Miquel 18-6-13
  use ecm !>>>> Is 31-8-13 
  use single_node
  use mitosis
  use neighboring !>>>> TT 22-05-2020


contains

!******************************************************************

! Declare the interface for POSIX fsync function

subroutine nexe

!interface
!  function fsync (fd) bind(c,name="fsync")
!  use iso_c_binding, only: c_int
!     integer(c_int), value :: fd
!     integer(c_int) :: fsync
!   end function fsync
!end interface


  implicit none                                                                !!>> HC 18-9-2020
  integer i,j,k,ii,jj,kk,ik,ikk,ick,ret
  real*8 a, fitelli                                                            !!>> HC 2-7-2020
  real*8 differ
  real*8 kplast,kvol,dvol
  character*140 nofifit, tarfit                                                !!>> HC 2-7-2020
  integer,save::truenodeos=0
  real*8 :: icr,icrr=0.005 !RZ 30-05-17

! GENETIC REGULATION OF CELL BEHAVIOURS

  ! extracellular matrix secretion

  if (npag(nparam_per_node+4)>0) then ; call should_I_secrete   ; end if

  if(ffu(1)==1)then

    ! cell polarization  THIS ONE SHOULD BE THE FIRST IN HERE
    if (npag(nparam_per_node+8)>0) then ; call polarization ; end if
    !nparam_per_node+9 is to tell that growth is polarized 
    ! cell growth
    if (npag(nparam_per_node+1)>0) then ; call should_I_grow        ; end if  !It considers also polar growth in it, with nparam_per_node+9

    !cell division
    !if (npag(nparam_per_node+2)>0) then ; 
    call should_I_divide ; 


    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if

    ! change the size of the cell required for division
    if (npag(nparam_per_node+10)>0) then; call change_minsize_for_div ; end if

    !nparam_per_node+11 is to orient division according to the chemical polarization and not according to the physical (hertwig) one

    !nparam_per_node+12 the larger the most asymetric (in mass) is the plane of division

    ! change the maximal number of nodes per cell before the cell divides
    if (npag(nparam_per_node+15)>0) then; call change_maxsize_for_div ; end if  !>>> Is 23-3-14

  else
    if (ffu(18)==0) then ; if (npag(nparam_per_node+1)>0) then ; call should_I_grow ; end if  ; end if  !>>> Is 18-4-15 !It considers also polar growth in it, with nparam_per_node+9 ACHTUNG
    if (npag(nparam_per_node+8)>0) then ; call polarization_single_new ; end if !!>> HC 4-2-2021 The old 2016 version of polarization included some small bugs solved by pfh in this new version
    if (npag(nparam_per_node+2)>0) then ; call should_I_divide_single ; end if
    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if  !>>Miquel20-3-14
    if (npag(nparam_per_node+13)>0) then; call emt ; end if   ! >>> Is 4-10-14
   !call emt_single
  end if

  ! epithelial-mesenchymal transition
  if (npag(nparam_per_node+13)>0) then; call emt ; end if


  if (ffu(4)==1.and.nd>2) then   ! IS 23-4-13 this eliminates the nodes that get alone for too long
    ik=1
    do while(ik<=nd) !;print*,"node(",ik,")%talone=",node(ik)%talone
      if (node(ik)%talone>ttalone) then !;print*,ik,"entra mort",node(ik)%tipus
        if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
          ikk=node(ik)%icel !;print*,"entra mesenq"
          call apoptosis(ik)
          if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
            if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
          end if
        else
          if (node(node(ik)%altre)%talone>ttalone) then
            ikk=node(ik)%icel
            call apoptosis(ik)
            if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
              if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
            end if
          end if
        end if
        ik=ik+1 ! Is it? >>> Is 16-1-14
      else
        ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
      end if
    end do
!    call cellbreak ! to see if the cell is split in two 
  end if

! GENETIC REGULATION OF NODE PROPERTIES
  do i=1,nd         ! we update that parameter in each cell that expresses the gene
                    ! WE ONLY UPDATES DE NODES IN WHICH THE GENE IS EXPRESSED, OTHERWISE WE LEAVE IT IS AS IT WAS 
    if(node(i)%fix==1) cycle  !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
    differ=1-node(i)%dif !;print*,"differ",differ

    ! DIFFERENTIATION
    if (npag(25)>0) then
      ii=25;a=0.0; 
      do k=1,npag(ii) 
        kk=whonpag(ii,k) ; 
        if (gex(i,kk)>0.0d0) then ; 
          a=a+gex(i,kk)*gen(kk)%e(ii)  
        endif
      end do
      node(i)%dif=node(i)%dif+a*delta
      if (node(i)%dif>1.0) node(i)%dif=1.0 
      if (node(i)%dif<0.0d0) node(i)%dif=0.0

      nodeo(i)%dif=node(i)%dif !>>Miquel17-9-14
    end if

    if(node(i)%tipus<2)then  !this for epithelial nodes  !>>Miquel12-5-14
      j=node(i)%altre
      if(ffu(6)==1) then        !plastic deformation   !>>Miquel5-2-14
        !kplast=node(i)%pla
        !if(fmeanl(i)<epsilod)then;ki=1/node(i)%rep;else;ki=1/node(i)%you;endif
        !if(fmeanl(j)<epsilod)then;kj=1/node(j)%rep;else;kj=1/node(j)%you;endif
        node(i)%pld=node(i)%pld+(node(i)%pla*fmeanl(i))*delta
        node(j)%pld=node(j)%pld+(node(j)%pla*fmeanl(j))*delta !; print*,"fmeanl",fmeanl(i),fmeanl(j)

        nodeo(i)%pld=node(i)%pld  !>>Miquel17-9-14
        nodeo(j)%pld=node(j)%pld  !>>Miquel17-9-14

      else
        node(i)%pld=0
        node(j)%pld=0
      end if
    
      if (npag(21)>0) then  ! Contraction by genes 
        ii=21 ; a=0 ; b=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) !;print*,"ij",i,j,"gex",gex(i,kk),gex(j,kk),"kk",kk
          if (gex(i,kk)>0.0d0.or.gex(j,kk)>0.0d0) then ; 
            a=a+gex(i,kk)*gen(kk)%e(ii)
            b=b+gex(j,kk)*gen(kk)%e(ii)
          endif 
        enddo
        node(i)%cod=nodeo(i)%cod+a*differ  !;print*,"a",a,"b",b
        node(j)%cod=nodeo(j)%cod+b*differ
        !;print*,"reqc",node(i)%cod,node(j)%cod
      else
        node(i)%cod=0 ;node(j)%cod=0
      end if

      if(ffu(10)==0) then        !volume conservation   !>>Miquel6-5-14                                  ! >>> Is 17-2-21
        !!kvol=node(i)%kvol                                                                              ! >>> Is 17-2-21
        !dvol=0.5*(node(i)%grd+node(j)%grd-(node(i)%eqd+node(j)%eqd))                                    ! >>> Is 17-2-21
        !node(i)%vod=node(i)%vod+node(i)%kvol*dvol*delta                                                 ! >>> Is 17-2-21
        !node(j)%vod=node(j)%vod+node(j)%kvol*dvol*delta                                                 ! >>> Is 17-2-21
        !nodeo(i)%vod=node(i)%vod  !>>Miquel17-9-14                                                      ! >>> Is 17-2-21
        !nodeo(j)%vod=node(j)%vod  !>>Miquel17-9-14                                                      ! >>> Is 17-2-21

        k=node(i)%altre                                                                                  ! >>> Is 17-2-21
        a=node(i)%cod+node(k)%cod  ! if that is not zero if means that there  is no volume conservation  ! >>> Is 17-2-21
        if (abs(a)>0.1d-4) then  !totally arbitrary                                                      ! >>> Is 17-2-21
          ! we check which side is bigger in absolute value                                              ! >>> Is 17-2-21
          if (abs(node(i)%cod)>abs(node(k)%cod)) then                                                    ! >>> Is 17-2-21
            node(k)%cod=-node(i)%cod                                                                     ! >>> Is 17-2-21
          else                                                                                           ! >>> Is 17-2-21
            node(i)%cod=-node(k)%cod                                                                     ! >>> Is 17-2-21
          end if                                                                                         ! >>> Is 17-2-21
        end if                                                                                           ! >>> Is 17-2-21
      else                                                                                               ! >>> Is 17-2-21
        node(i)%vod=0 ; node(j)%vod=0                                                                    ! >>> Is 17-2-21
      end if                                                                                             ! >>> Is 17-2-21

      if (ffu(11)==1) call diffusion_of_reqcr       !diffusion of reqcr ! >>> Is 25-5-14

      a=node(i)%add-node(i)%eqd
      node(i)%eqd=node(i)%grd+node(i)%cod+node(i)%pld+node(i)%vod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(i)%eqd>df_reqmax) node(i)%eqd=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%eqd<reqmin) node(i)%eqd=reqmin !the req can be deformed
      node(i)%add=node(i)%eqd+a

      b=node(j)%add-node(j)%eqd
      node(j)%eqd=node(j)%grd+node(j)%cod+node(j)%pld+node(j)%vod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(j)%eqd>df_reqmax) node(j)%eqd=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(j)%eqd<reqmin) node(j)%eqd=reqmin !the req can be deformed
      node(j)%add=node(j)%eqd+b
      if (npag(6)>0) then;ii=6;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%add=nodeo(i)%add+a*differ;node(j)%add=nodeo(j)%add+b*differ;
      if (node(i)%add<0) node(i)%add=0.0;if (node(j)%add<0) node(j)%add=0.0;end if
      if (npag(7)>0) then;ii=7;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%you=nodeo(i)%you+a*differ;node(j)%you=nodeo(j)%you+b*differ;
      if (node(i)%you<0) node(i)%you=0.0;if (node(j)%you<0) node(j)%you=0.0;end if
      if (npag(8)>0) then;ii=8;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%adh=nodeo(i)%adh+a*differ;node(j)%adh=nodeo(j)%adh+b*differ;
      if (node(i)%adh<0) node(i)%adh=0.0;if (node(j)%adh<0) node(j)%adh=0.0;end if
      if (npag(9)>0) then;ii=9;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%rep=nodeo(i)%rep+a*differ;node(j)%rep=nodeo(j)%rep+b*differ;
      if (node(i)%rep<0) node(i)%rep=0.0;if (node(j)%rep<0) node(j)%rep=0.0;end if
      if (npag(10)>0) then;ii=10;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%rec=nodeo(i)%rec+a*differ;node(j)%rec=nodeo(j)%rec+b*differ;
      if (node(i)%rec<0) node(i)%rec=0.0;if (node(j)%rec<0) node(j)%rec=0.0;end if
      if (npag(11)>0) then;ii=11;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%erp=nodeo(i)%erp+a*differ;node(j)%erp=nodeo(j)%erp+b*differ;
      if (node(i)%erp<0) node(i)%erp=0.0;if (node(j)%erp<0) node(j)%erp=0.0;end if
      if (npag(12)>0) then;ii=12;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%est=nodeo(i)%est+a*differ;node(j)%est=nodeo(j)%est+b*differ;
      if (node(i)%est<0) node(i)%est=0.0;if (node(j)%est<0) node(j)%est=0.0;end if
      if (npag(13)>0) then;ii=13;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%eqs=nodeo(i)%eqs+a*differ;node(j)%eqs=nodeo(j)%eqs+b*differ;
      if (node(i)%eqs<0) node(i)%eqs=0.0;if (node(j)%eqs<0) node(j)%eqs=0.0;end if
      if (npag(14)>0) then;ii=14;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%hoo=nodeo(i)%hoo+a*differ;node(j)%hoo=nodeo(j)%hoo+b*differ;
      if (node(i)%hoo<0) node(i)%hoo=0.0;if (node(j)%hoo<0) node(j)%hoo=0.0;end if
      !if (npag(15)>0) then;ii=15;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units      !!>> HC 28-1-2021
      !if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                     !!>> HC 28-1-2021 This has been commented out because
      !if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                               !!>> HC 28-1-2021 there is no way that a gene can affect
      !node(i)%mov=nodeo(i)%mov+a*differ;node(j)%mov=nodeo(j)%mov+b*differ;                                !!>> HC 28-1-2021 epitelial cell movement due to noise
      !if (node(i)%mov<0) node(i)%mov=0.0;if (node(j)%mov<0) node(j)%mov=0.0;end if                        !!>> HC 28-1-2021
      !if (npag(16)>0) then;ii=16;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  !!>> HC 28-1-2021 It is posible in mesenchyme because of pseudopodia.
      !if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                     !!>> HC 28-1-2021
      !if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                               !!>> HC 28-1-2021
      !node(i)%dmo=nodeo(i)%dmo+a*differ;node(j)%dmo=nodeo(j)%dmo+b*differ;                                !!>> HC 28-1-2021
      !if (node(i)%dmo<0) node(i)%dmo=0.0;if (node(j)%dmo<0) node(j)%dmo=0.0;end if                        !!>> HC 28-1-2021
      if (npag(27)>0) then; ii=27;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%pla=nodeo(i)%pla+a*differ;node(j)%pla=nodeo(j)%pla+b*differ;
      if (node(i)%pla<0) node(i)%pla=0.0;if (node(j)%pla<0) node(j)%pla=0.0;end if
      if (npag(28)>0) then;ii=28;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;
      node(i)%kvol=nodeo(i)%kvol+a*differ;node(j)%kvol=nodeo(j)%kvol+b*differ;
      if (node(i)%kvol<0) node(i)%kvol=0.0;if (node(j)%kvol<0) node(j)%kvol=0.0;end if


    else if (node(i)%tipus>2)then !this for mesenchyme and ECM

      if (npag(21)>0) then 
        ii=21 ; a=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) 
          if (gex(i,kk)>0.0d0) then ; 
            a=a+gex(i,kk)*gen(kk)%e(ii)
          endif 
        enddo
        node(i)%cod=nodeo(i)%cod+a*differ*delta !wa in req-space units
      else
        node(i)%cod=0
      end if

      a=node(i)%add-node(i)%eqd
      node(i)%eqd=node(i)%grd+node(i)%cod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(i)%eqd>df_reqmax) node(i)%eqd=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%eqd<reqmin) node(i)%eqd=reqmin !the req can be deformed
      node(i)%add=node(i)%eqd+a

      if (npag(6)>0) then;ii=6;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%add=nodeo(i)%add+a*differ;if (node(i)%add<0) node(i)%add=0.0;end if
      if (npag(7)>0) then;ii=7;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%you=nodeo(i)%you+a*differ;if (node(i)%you<0) node(i)%you=0.0;end if
      if (npag(8)>0) then;ii=8;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%adh=nodeo(i)%adh+a*differ;if (node(i)%adh<0) node(i)%adh=0.0;end if
      if (npag(9)>0) then;ii=9;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%rep=nodeo(i)%rep+a*differ;if (node(i)%rep<0) node(i)%rep=0.0;end if
      if (npag(10)>0) then;ii=10;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%rec=nodeo(i)%rec+a;if (node(i)%rec<0) node(i)%rec=0.0;end if
      if (npag(11)>0) then;ii=11;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%erp=nodeo(i)%erp+a*differ;if (node(i)%erp<0) node(i)%erp=0.0;end if
      if (npag(12)>0) then;ii=12;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%est=nodeo(i)%est+a*differ;if (node(i)%est<0) node(i)%est=0.0;end if
      if (npag(13)>0) then;ii=13;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%eqs=nodeo(i)%eqs+a*differ;if (node(i)%eqs<0) node(i)%eqs=0.0;end if
      if (npag(14)>0) then;ii=14;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%hoo=nodeo(i)%hoo+a*differ;if (node(i)%hoo<0) node(i)%hoo=0.0;end if
      if (npag(15)>0) then;ii=15;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo;node(i)%mov=nodeo(i)%mov+a*differ;if (node(i)%mov<0) node(i)%mov=0.0;end if
      if (npag(16)>0) then;ii=16;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%e(ii);endif;enddo;node(i)%dmo=nodeo(i)%dmo+a*differ;if (node(i)%dmo<0) node(i)%dmo=0.0;end if


    end if
  end do

  if (ffu(21)==1)then       !!>> HC 17-11-2020  THIS FFU CONTROLS WHETHER WE ARE APPLYING FILTERS OR NOT   
     call filter            !!>> HC 12-2-2021                                                                
  endif                     !!>> HC 12-2-2021

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   NEXE GRADUAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nexe_gradual !!>> HC 12-2-2021 This subrutine includes the code written by Rz and pfh 
                        !!>> HC 12-2-2021 to apply genetic changes in node properties gradually (5% of the previous value)
                        !!>> HC 12-2-2021 it also limits the proportion ADD/EQD
                        !!>> HC 12-2-2021 it was formerly included in nexe2
  implicit none                            !!>> HC 12-2-2021           
  integer i,j,k,ii,jj,kk,ik,ikk,ick        !!>> HC 12-2-2021
  real*8 a                                 !!>> HC 12-2-2021
  real*8 differ                            !!>> HC 12-2-2021
  real*8 dvol                              !!>> HC 12-2-2021
  character*300 cx                         !!>> HC 12-2-2021
  real*8 mz3,mz4                           !!>> HC 12-2-2021

  real*8 :: bbb                            !!>> HC 12-2-2021
  real*8 :: scui,maxa,maxs(3),mins(3)      !!>> HC 12-2-2021
  integer :: nd_a,nd_t,ndm,first=0         !!>> HC 12-2-2021

  integer,save::truenodeos=0               !!>> HC 12-2-2021
  real*8 :: icr,icrr=0.005 !RZ 30-05-17    !!>> HC 12-2-2021
 

        
  ! GENETIC REGULATION OF CELL BEHAVIOURS

  if (ndo==0)then                    !!>> HC 12-2-2021
     ndoo=nd                         !!>> HC 12-2-2021
     ndo=1                           !!>> HC 12-2-2021
!     print*,"ndo",ndo,"ndoo",nd      !!>> HC 12-2-2021
  end if                             !!>> HC 12-2-2021

  !extracellular matrix secretion
  if (npag(nparam_per_node+4)>0) then ; call should_I_secrete   ; end if                !!>> HC 12-2-2021

  if (ffu(1)==1)then                                                                    !!>> HC 12-2-2021
     ! cell polarization  THIS ONE SHOULD BE THE FIRST IN HERE                          !!>> HC 12-2-2021
     if (npag(nparam_per_node+8)>0) then ; call polarization ; end if                   !!>> HC 12-2-2021

     !nparam_per_node+9 is to tell that growth is polarized                             !!>> HC 12-2-2021
     ! cell growth                                                                      !!>> HC 12-2-2021
     if (npag(nparam_per_node+1)>0) then ; call should_I_grow ; end if                  !!>> HC 12-2-2021 !It considers also polar growth in it, with nparam_per_node+9

     !cell division                                                                     !!>> HC 12-2-2021
     call should_I_divide ;                                                             !!>> HC 12-2-2021

     !cell apoptosis                                                                                     !!>> HC 12-2-2021
     if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if     !!>> HC 12-2-2021

     !change the size of the cell required for division                                                  !!>> HC 12-2-2021
     if (npag(nparam_per_node+10)>0) then; call change_minsize_for_div ; end if                          !!>> HC 12-2-2021

     !nparam_per_node+11 is to orient division according to the chemical polarization and not according to the physical (hertwig) one  !!>> HC 12-2-2021
     !nparam_per_node+12 the larger the most asymetric (in mass) is the plane of division                                              !!>> HC 12-2-2021
     ! change the maximal number of nodes per cell before the cell divides                                                             !!>> HC 12-2-2021
     if (npag(nparam_per_node+15)>0) then; call change_maxsize_for_div ; end if  !>>> Is 23-3-14                                       !!>> HC 12-2-2021

  else

     if (ffu(18)==0) then ; if (npag(nparam_per_node+1)>0) then ; call should_I_grow ; end if  ; end if  !>>> Is 18-4-15 !!>> HC 12-2-2021 !It considers also polar growth in it, with nparam_per_node+9 ACHTUNG
     if (npag(nparam_per_node+8)>0) then ; call polarization_single ; end if                                                  !!>> HC 12-2-2021
     if (npag(nparam_per_node+2)>0) then ; call should_I_divide_single ; end if                                               !!>> HC 12-2-2021
     ! cell apoptosis                                                                                                         !!>> HC 12-2-2021
     if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if  !>>Miquel20-3-14        !!>> HC 12-2-2021
     if (npag(nparam_per_node+13)>0) then; call emt ; end if   ! >>> Is 4-10-14                                               !!>> HC 12-2-2021
  endif                                                                                                                       !!>> HC 12-2-2021

  !epithelial-mesenchymal transition                                                                                          !!>> HC 12-2-2021
  if (npag(nparam_per_node+13)>0) then; call emt ; end if                                                                     !!>> HC 12-2-2021

  if (ffu(4)==1.and.nd>2) then   ! IS 23-4-13 this eliminates the nodes that get alone for too long                           !!>> HC 12-2-2021
      ik=1                                                                                                                    !!>> HC 12-2-2021
      do while(ik<=nd) !;print*,"node(",ik,")%talone=",node(ik)%talone                                                        !!>> HC 12-2-2021
         if (node(ik)%talone>ttalone) then !;print*,ik,"entra mort",node(ik)%tipus                                            !!>> HC 12-2-2021
            if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely                        !!>> HC 12-2-2021
               ikk=node(ik)%icel !;print*,"entra mesenq"                                                                      !!>> HC 12-2-2021
               call apoptosis(ik)                                                                                             !!>> HC 12-2-2021
               if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node                          !!>> HC 12-2-2021
                  if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left           !!>> HC 12-2-2021
               end if                                                                                                         !!>> HC 12-2-2021
            else                                                                                                              !!>> HC 12-2-2021
               if (node(node(ik)%altre)%talone>ttalone) then                                                                  !!>> HC 12-2-2021
                  ikk=node(ik)%icel                                                                                           !!>> HC 12-2-2021
                  call apoptosis(ik)                                                                                          !!>> HC 12-2-2021
                  if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node                       !!>> HC 12-2-2021
                     if (cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left       !!>> HC 12-2-2021
                  end if                                                                                                      !!>> HC 12-2-2021
               end if                                                                                                         !!>> HC 12-2-2021
            end if                                                                                                            !!>> HC 12-2-2021
            ik=ik+1 ! Is it? >>> Is 16-1-14                                                                                   !!>> HC 12-2-2021
         else                                                                                                                 !!>> HC 12-2-2021
            ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis !!>> HC 12-2-2021
         end if                                                                                                               !!>> HC 12-2-2021
      end do                                                                                                                  !!>> HC 12-2-2021
      !    call cellbreak ! to see if the cell is split in two                                                                !!>> HC 12-2-2021
  end if                                                                                                                      !!>> HC 12-2-2021

  

  ! GENETIC REGULATION OF NODE PROPERTIES                                                                                  !!>> HC 12-2-2021
  do i=1,nd         ! we update that parameter in each cell that expresses the gene                                        !!>> HC 12-2-2021
                    ! WE ONLY UPDATES DE NODES IN WHICH THE GENE IS EXPRESSED, OTHERWISE WE LEAVE IT IS AS IT WAS          !!>> HC 12-2-2021

     if(node(i)%fix==1) cycle  !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14  !!>> HC 12-2-2021
     if(node(i)%fix==2) cycle                                                                                              !!>> HC 12-2-2021
     differ=1-node(i)%dif !;print*,"differ",differ                                                                         !!>> HC 12-2-2021

     ! DIFFERENTIATION                                                                                                     !!>> HC 12-2-2021
     if (npag(25)>0) then                                                                                                  !!>> HC 12-2-2021
        ii=25;a=0.0;                                                                                                       !!>> HC 12-2-2021
        do k=1,npag(ii)                                                                                                    !!>> HC 12-2-2021
           kk=whonpag(ii,k) ;                                                                                              !!>> HC 12-2-2021
           if (gex(i,kk)>0.0d0) then ;                                                                                     !!>> HC 12-2-2021
              a=a+gex(i,kk)*gen(kk)%e(ii)                                                                                  !!>> HC 12-2-2021
           endif                                                                                                           !!>> HC 12-2-2021
        end do                                                                                                             !!>> HC 12-2-2021
        ! NEW INDEPENDENT DIFFE !RZ 22-3-16                                                                                !!>> HC 12-2-2021
        if(ng.ge.11)  a=gen(11)%e(25)*0.1 ! This makes that differ=0 in 10.000-100.000 iterations                          !!>> HC 12-2-2021
        node(i)%dif=node(i)%dif+a*delta                                                                                    !!>> HC 12-2-2021
        if (node(i)%dif>1.0) node(i)%dif=1.0                                                                               !!>> HC 12-2-2021
        if (node(i)%dif<0.0d0) node(i)%dif=0.0                                                                             !!>> HC 12-2-2021
        nodeo(i)%dif=node(i)%dif !>>Miquel17-9-14                                                                          !!>> HC 12-2-2021
     end if                                                                                                                !!>> HC 12-2-2021

     if (node(i)%tipus<2)then  !this for epithelial nodes  !>>Miquel12-5-14                                                !!>> HC 12-2-2021   
        j=node(i)%altre                                                                                                    !!>> HC 12-2-2021
        if (ffu(6)==1)then        !plastic deformation   !>>Miquel5-2-14                                                  !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually?????                                                     !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                                           !!>> HC 12-2-2021
           aa=node(i)%pld; bb=node(j)%pld  ! previous values RZ 30-5-17                                                    !!>> HC 12-2-2021
                                              
           node(i)%pld=node(i)%pld+(node(i)%pla*fmeanl(i))*delta                                                           !!>> HC 12-2-2021
           node(j)%pld=node(j)%pld+(node(j)%pla*fmeanl(j))*delta !; print*,"fmeanl",fmeanl(i),fmeanl(j)                    !!>> HC 12-2-2021
           ! RZ gradual 30-5-17                                                                                            !!>> HC 12-2-2021
           !if (node(i)%pld>aa*(1+icr).and.aa.ne.0) node(i)%pld=aa*(1+icr)                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (node(j)%pld>bb*(1+icr).and.bb.ne.0) node(j)%pld=bb*(1+icr)                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (node(i)%pld<aa*(1-icr).and.aa.ne.0) node(i)%pld=aa*(1-icr)                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (node(j)%pld<bb*(1-icr).and.bb.ne.0) node(j)%pld=bb*(1-icr)                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (aa==0.and.node(i)%pld>icr) node(i)%pld=icr                                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (bb==0.and.node(j)%pld>icr) node(j)%pld=icr                                                                  !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (aa==0.and.node(i)%pld<icr*(-1)) node(i)%pld=icr*(-1)                                                        !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (bb==0.and.node(j)%pld<icr*(-1)) node(j)%pld=icr*(-1)                                                        !!>> HC 12-2-2021 !!>> TT 2-9-2021
           ! RZ end changes                                                                                                !!>> HC 12-2-2021
           nodeo(i)%pld=node(i)%pld  !>>Miquel17-9-14                                                                      !!>> HC 12-2-2021
           nodeo(j)%pld=node(j)%pld  !>>Miquel17-9-14                                                                      !!>> HC 12-2-2021
        else                                                                                                               !!>> HC 12-2-2021
           node(i)%pld=0                                                                                                   !!>> HC 12-2-2021
           node(j)%pld=0                                                                                                   !!>> HC 12-2-2021
        end if                                                                                                             !!>> HC 12-2-2021
    
        if (npag(21)>0) then  ! Contraction by genes                                                                       !!>> HC 12-2-2021
           ii=21 ; a=0 ; b=0                                                                                               !!>> HC 12-2-2021
           do k=1,npag(ii) ;                                                                                               !!>> HC 12-2-2021
              kk=whonpag(ii,k) !;print*,"ij",i,j,"gex",gex(i,kk),gex(j,kk),"kk",kk                                         !!>> HC 12-2-2021
              if (gex(i,kk)>0.0d0.or.gex(j,kk)>0.0d0) then                                                                 !!>> HC 12-2-2021
                 a=a+gex(i,kk)*gen(kk)%e(ii)                                                                               !!>> HC 12-2-2021
                 b=b+gex(j,kk)*gen(kk)%e(ii)                                                                               !!>> HC 12-2-2021
              endif                                                                                                        !!>> HC 12-2-2021
           enddo                                                                                                           !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually?????                                                     !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                                           !!>> HC 12-2-2021
           aaa=node(i)%cod; bbb=node(j)%cod  ! previous values RZ 30-5-17                                                  !!>> HC 12-2-2021
           aa=0.250d0; bb=0.250d0!nodeoini(5,1); bb=nodeoini(5,1)  ! previous values RZ 30-5-17                            !!>> HC 2-7-2021 COD is not a property, so its nodeo is =0
           node(i)%cod=nodeo(i)%cod+a!*differ  !pfh no differ 26-11-15                                                     !!>> HC 12-2-2021 we have to set up an arbitrary nodeo value
           node(j)%cod=nodeo(j)%cod+b!*differ  !pfh no differ 26-11-15                                                     !!>> HC 12-2-2021 to make changes gradual. This value was in the
           ! RZ gradual 30-5-17                                                                                            !!>> HC 12-2-2021 original code by RZ and pfh and it works, feel free
           !if (node(i)%cod>aaa+aa*(icr).and.aa.ne.0) node(i)%cod=aaa+aa*(icr)                                              !!>> HC 12-2-2021 to change it  !!>> TT 2-9-2021
           !if (node(j)%cod>bbb+bb*(icr).and.bb.ne.0) node(j)%cod=bbb+bb*(icr)                                              !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (node(i)%cod<aaa-aa*(icr).and.aa.ne.0) node(i)%cod=aaa-aa*(icr)                                              !!>> HC 12-2-2021 !!>> TT 2-9-2021
           !if (node(j)%cod<bbb-bb*(icr).and.bb.ne.0) node(j)%cod=bbb-bb*(icr)                                              !!>> HC 12-2-2021 !!>> TT 2-9-2021
        else                                                                                                               !!>> HC 12-2-2021
           node(i)%cod=0 ;node(j)%cod=0                                                                                    !!>> HC 12-2-2021
        end if                                                                                                             !!>> HC 12-2-2021
        
        
        if (npag(22)>0) then  ! growth by genes                                                                            !!>> HC 7-10-2021
           ii=22 ; a=0.0d0 ; b=0.0d0                                                                                       !!>> HC 7-10-2021
           do k=1,npag(ii) ;                                                                                               !!>> HC 7-10-2021
              kk=whonpag(ii,k) !;print*,"ij",i,j,"gex",gex(i,kk),gex(j,kk),"kk",kk                                         !!>> HC 7-10-2021
              if (gex(i,kk)>0.0d0.or.gex(j,kk)>0.0d0) then                                                                 !!>> HC 7-10-2021
                 a=a+gex(i,kk)*gen(kk)%e(ii)                                                                               !!>> HC 7-10-2021
                 b=b+gex(j,kk)*gen(kk)%e(ii)                                                                               !!>> HC 7-10-2021
              endif                                                                                                        !!>> HC 7-10-2021
           enddo                                                                                                           !!>> HC 7-10-2021
           ! Rz 18-5-17 makes node properties to change gradually?????                                                     !!>> HC 7-10-2021
           icr=icrr ! RZ 30-5-17                                                                                           !!>> HC 7-10-2021 
           node(i)%grd=nodeo(i)%grd+a!*differ  !pfh no differ 26-11-15                                                     !!>> HC 7-10-2021 
           node(j)%grd=nodeo(j)%grd+b!*differ  !pfh no differ 26-11-15                                                     !!>> HC 7-10-2021 
        end if                                                                                                             !!>> HC 7-10-2021


      if(ffu(10)==0) then        !volume conservation   !>>Miquel6-5-14                                  ! >>> Is 17-2-21
        !!kvol=node(i)%kvol                                                                              ! >>> Is 17-2-21
        !dvol=0.5*(node(i)%grd+node(j)%grd-(node(i)%eqd+node(j)%eqd))                                    ! >>> Is 17-2-21
        !node(i)%vod=node(i)%vod+node(i)%kvol*dvol*delta                                                 ! >>> Is 17-2-21
        !node(j)%vod=node(j)%vod+node(j)%kvol*dvol*delta                                                 ! >>> Is 17-2-21
        !nodeo(i)%vod=node(i)%vod  !>>Miquel17-9-14                                                      ! >>> Is 17-2-21
        !nodeo(j)%vod=node(j)%vod  !>>Miquel17-9-14                                                      ! >>> Is 17-2-21

        k=node(i)%altre                                                                                  ! >>> Is 17-2-21
        a=node(i)%cod+node(k)%cod  ! if that is not zero if means that there  is no volume conservation  ! >>> Is 17-2-21
        if (abs(a)>0.1d-4) then  !totally arbitrary                                                      ! >>> Is 17-2-21
          ! we check which side is bigger in absolute value                                              ! >>> Is 17-2-21
          if (abs(node(i)%cod)>abs(node(k)%cod)) then                                                    ! >>> Is 17-2-21
            node(k)%vod=-node(i)%cod*node(i)%kvol                                                        ! >>> Is 17-2-21 !!>> HC 6-7-2021 kvol is the percentage of volume conserved
          else                                                                                           ! >>> Is 17-2-21
            node(i)%vod=-node(k)%cod*node(i)%kvol                                                        ! >>> Is 17-2-21 !!>> HC 6-7-2021 kvol is the percentage of volume conserved
          end if                                                                                         ! >>> Is 17-2-21
        end if                                                                                           ! >>> Is 17-2-21
      else                                                                                               ! >>> Is 17-2-21
        node(i)%vod=0 ; node(j)%vod=0                                                                    ! >>> Is 17-2-21
      end if                                                                                             ! >>> Is 17-2-21
            
!        if (ffu(10)==0) then        !volume conservation   !>>Miquel6-5-14 ! >>> Is 24-5-14                                !!>> HC 12-2-2021
!           dvol=0.5*(node(i)%grd+node(j)%grd-(node(i)%eqd+node(j)%eqd))                                                    !!>> HC 12-2-2021
!           icr=icrr ! RZ 30-5-17                                                                                           !!>> HC 12-2-2021
!           aa=node(i)%vod; bb=node(j)%vod  ! previous values RZ 30-5-17                                                    !!>> HC 12-2-2021
!           node(i)%vod=node(i)%vod+node(i)%kvol*dvol*delta                                                                 !!>> HC 12-2-2021
!           node(j)%vod=node(j)%vod+node(j)%kvol*dvol*delta                                                                 !!>> HC 12-2-2021
!           ! RZ gradual 30-5-17                                                                                            !!>> HC 12-2-2021
!           if (node(i)%vod>aa*(1+icr).and.aa.ne.0) node(i)%vod=aa*(1+icr)                                                  !!>> HC 12-2-2021
!           if (node(j)%vod>bb*(1+icr).and.bb.ne.0) node(j)%vod=bb*(1+icr)                                                  !!>> HC 12-2-2021
!           if (node(i)%vod<aa*(1-icr).and.aa.ne.0) node(i)%vod=aa*(1-icr)                                                  !!>> HC 12-2-2021
!           if (node(j)%vod<bb*(1-icr).and.bb.ne.0) node(j)%vod=bb*(1-icr)                                                  !!>> HC 12-2-2021
!           if (aa==0.and.node(i)%vod>icr) node(i)%vod=icr                                                                  !!>> HC 12-2-2021
!           if (bb==0.and.node(j)%vod>icr) node(j)%vod=icr                                                                  !!>> HC 12-2-2021
!           if (aa==0.and.node(i)%vod<icr*(-1)) node(i)%vod=icr*(-1)                                                        !!>> HC 12-2-2021
!           if (bb==0.and.node(j)%vod<icr*(-1)) node(j)%vod=icr*(-1)                                                        !!>> HC 12-2-2021
!           ! RZ end changes                                                                                                !!>> HC 12-2-2021
!           nodeo(i)%vod=node(i)%vod  !>>Miquel17-9-14                                                                      !!>> HC 12-2-2021
!           nodeo(j)%vod=node(j)%vod  !>>Miquel17-9-14                                                                      !!>> HC 12-2-2021
!        else                                                                                                               !!>> HC 12-2-2021
!           node(i)%vod=0 ; node(j)%vod=0                                                                                   !!>> HC 12-2-2021
!        end if                                                                                                             !!>> HC 12-2-2021

        if (ffu(11)==1) call diffusion_of_reqcr       !diffusion of reqcr ! >>> Is 25-5-14                                 !!>> HC 12-2-2021
        
        a=node(i)%add-node(i)%eqd !>>>Miquel9-9-15                                                                         !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually                                                             !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                                              !!>> HC 12-2-2021
        aa=node(i)%eqd; aaa=node(i)%eqd; bb=node(j)%eqd; bbb=node(j)%eqd  ! previous values RZ 30-5-17                     !!>> HC 12-2-2021 !!>> TT 2-9-2021
        
      
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        !!! HERE WE ACTUALLY CHANGE THE EQD AND ADD VALUES
        !!>> TT 2-9-2021 No gradual change of COD, GRD and PLD, but a gradual change of EQD
        node(i)%eqd=node(i)%grd+node(i)%cod+node(i)%pld+node(i)%vod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation !!>> HC 12-2-2021
        if (node(i)%eqd>aaa+aa*icr) node(i)%eqd=aaa+aa*icr !!>> TT 2-9-2021
        if (node(i)%eqd<aaa-aa*icr) node(i)%eqd=aaa-aa*icr !!>> TT 2-9-2021
        if(node(i)%eqd>df_reqmax) node(i)%eqd=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14      !!>> HC 12-2-2021
        if(node(i)%eqd<reqmin) node(i)%eqd=reqmin !the req can be deformed                                                 !!>> HC 12-2-2021
        node(i)%add=node(i)%eqd+a !>>>Miquel9-9-15                                                                         !!>> HC 15-2-2021 it is node not nodeo  
        if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15                                        !!>> HC 15-2-2021 it is node not nodeo          
        b=node(j)%add-node(j)%eqd !>>>Miquel9-9-15                                                                         !!>> HC 12-2-2021
        node(j)%eqd=node(j)%grd+node(j)%cod+node(j)%pld+node(j)%vod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation !!>> HC 12-2-2021
        if (node(j)%eqd>bbb+bb*icr) node(j)%eqd=bbb+bb*icr !!>> TT 2-9-2021
        if (node(j)%eqd<bbb-bb*icr) node(j)%eqd=bbb-bb*icr !!>> TT 2-9-2021
        if(node(j)%eqd>df_reqmax) node(j)%eqd=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14                                        !!>> HC 12-2-2021
        if(node(j)%eqd<reqmin) node(j)%eqd=reqmin !the req can be deformed                                                                                   !!>> HC 12-2-2021
        node(j)%add=node(j)%eqd+b !>>>Miquel9-9-15                                                                                   !!>> HC 15-2-2021 it is node not nodeo                          
!        if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15 1.5                                              !!>> HC 12-2-2021
!        if(node(i)%add.ge.node(i)%eqd*1.2) node(i)%add=node(i)%eqd*1.2  !pfh 9-6-15 2.5                                              !!>> HC 12-2-2021
!        if(node(j)%add.le.node(j)%eqd*1.0) node(j)%add=node(j)%eqd*1.0  !pfh 9-6-15                                                  !!>> HC 12-2-2021
!        if(node(j)%add.ge.node(j)%eqd*1.2) node(j)%add=node(j)%eqd*1.2  !pfh 9-6-15                                                  !!>> HC 12-2-2021
        !pfh no differ 26-11-15                                                                                                      !!>> HC 12-2-2021  
  


     
        if (npag(6)>0) then                                                                                             
           ii=6;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units                !!>> HC 12-2-2021
           if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                             !!>> HC 12-2-2021
           if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                       !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually?????                                 !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                       !!>> HC 12-2-2021
           aaa=node(i)%add; bbb=node(j)%add  ! previous values RZ 30-5-17                              !!>> HC 12-2-2021
           aa=nodeo(i)%add; bb=nodeo(j)%add  ! previous values RZ 30-5-17                              !!>>HC 15-2-2021
           node(i)%add=nodeo(i)%add+a;node(j)%add=nodeo(j)%add+b                                       !!>>HC 15-2-2021
           ! RZ gradual 30-5-17                                                                        !!>>HC 15-2-2021
           if (node(i)%add>aaa+aa*(icr).and.aa.ne.0) node(i)%add=aaa+aa*(icr)                          !!>>HC 15-2-2021
           if (node(j)%add>bbb+bb*(icr).and.bb.ne.0) node(j)%add=bbb+bb*(icr)                          !!>>HC 15-2-2021
           if (node(i)%add<aaa-aa*(icr).and.aa.ne.0) node(i)%add=aaa-aa*(icr)                          !!>>HC 15-2-2021
           if (node(j)%add<bbb-bb*(icr).and.bb.ne.0) node(j)%add=bbb-bb*(icr)                          !!>>HC 15-2-2021
           if (node(i)%add<0) node(i)%add=0.0;if (node(j)%add<0) node(j)%add=0.0                       !!>>HC 15-2-2021

!           if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15                 !!>>HC 15-2-2021
!           if(node(i)%add.ge.node(i)%eqd*1.2) node(i)%add=node(i)%eqd*1.2  !pfh 9-6-15                 !!>>HC 15-2-2021
!           if(node(j)%add.le.node(j)%eqd*1.0) node(j)%add=node(j)%eqd*1.0  !pfh 9-6-15                 !!>>HC 15-2-2021
!           if(node(j)%add.ge.node(j)%eqd*1.2) node(j)%add=node(j)%eqd*1.2  !pfh 9-6-15                 !!>>HC 15-2-2021
           
!        else                                                                                           !!>>HC 15-2-2021
           
!           icr=icrr ! RZ 30-5-17                                                                       !!>>HC 15-2-2021
!           aaa=node(i)%add; bbb=node(j)%add  ! previous values RZ 30-5-17                              !!>>HC 15-2-2021
!           aa=nodeo(i)%add; bb=nodeo(j)%add  ! previous values RZ 30-5-17                              !!>>HC 15-2-2021
!           node(i)%add=nodeo(i)%add ; node(j)%add=nodeo(j)%add !pfh 9-6-15                             !!>>HC 15-2-2021
          ! RZ gradual 30-5-17                                                                         !!>>HC 15-2-2021
!           if (node(i)%add>aaa+aa*(icr).and.aa.ne.0) node(i)%add=aaa+aa*(icr)                          !!>>HC 15-2-2021
!           if (node(j)%add>bbb+bb*(icr).and.bb.ne.0) node(j)%add=bbb+bb*(icr)                          !!>>HC 15-2-2021
!           if (node(i)%add<aaa-aa*(icr).and.aa.ne.0) node(i)%add=aaa-aa*(icr)                          !!>>HC 15-2-2021
!           if (node(j)%add<bbb-bb*(icr).and.bb.ne.0) node(j)%add=bbb-bb*(icr)                          !!>>HC 15-2-2021
                                     
        end if

  !      if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15 1.5                !!>>HC 15-2-2021
  !      if(node(i)%add.ge.node(i)%eqd*1.2) node(i)%add=node(i)%eqd*1.2  !pfh 9-6-15 2.5                !!>>HC 15-2-2021
  !      if(node(j)%add.le.node(j)%eqd*1.0) node(j)%add=node(j)%eqd*1.0  !pfh 9-6-15                    !!>>HC 15-2-2021
  !      if(node(j)%add.ge.node(j)%eqd*1.2) node(j)%add=node(j)%eqd*1.2  !pfh 9-6-15                    !!>>HC 15-2-2021
        if (npag(7)>0) then;ii=7;a=0 ; b=0 ;                                                           !!>>HC 15-2-2021
        do k=1,npag(ii)                                                                                !!>>HC 15-2-2021
           kk=whonpag(ii,k) ; !wa in req-space units                                                   !!>>HC 15-2-2021
           if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                             !!>>HC 15-2-2021
           if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif                              !!>>HC 15-2-2021
        enddo                                                                                          !!>>HC 15-2-2021
        
        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021 
        aaa=node(i)%you; bbb=node(j)%you  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%you; bb=nodeo(j)%you  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%you=nodeo(i)%you+a;node(j)%you=nodeo(j)%you+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%you>aaa+aa*(icr).and.aa.ne.0) node(i)%you=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%you>bbb+bb*(icr).and.bb.ne.0) node(j)%you=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%you<aaa-aa*(icr).and.aa.ne.0) node(i)%you=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%you<bbb-bb*(icr).and.bb.ne.0) node(j)%you=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%you<0) node(i)%you=0.0;if (node(j)%you<0) node(j)%you=0.0;end if                   !!>> HC 12-2-2021
        if (npag(8)>0) then;ii=8;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units   !!>> HC 12-2-2021
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021

        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%adh; bbb=node(j)%adh  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%adh; bb=nodeo(j)%adh  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        
        node(i)%adh=nodeo(i)%adh+a;node(j)%adh=nodeo(j)%adh+b                                          !!>> HC 12-2-2021
         
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%adh>aaa+aa*(icr).and.aa.ne.0) node(i)%adh=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%adh>bbb+bb*(icr).and.bb.ne.0) node(j)%adh=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%adh<aaa-aa*(icr).and.aa.ne.0) node(i)%adh=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%adh<bbb-bb*(icr).and.bb.ne.0) node(j)%adh=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%adh<0) node(i)%adh=0.0;if (node(j)%adh<0) node(j)%adh=0.0;end if                   !!>> HC 12-2-2021
        if (npag(9)>0) then;ii=9;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units   !!>> HC 12-2-2021
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%rep; bbb=node(j)%rep  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%rep; bb=nodeo(j)%rep  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%rep=nodeo(i)%rep+a;node(j)%rep=nodeo(j)%rep+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%rep>aaa+aa*(icr).and.aa.ne.0) node(i)%rep=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%rep>bbb+bb*(icr).and.bb.ne.0) node(j)%rep=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%rep<aaa-aa*(icr).and.aa.ne.0) node(i)%rep=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%rep<bbb-bb*(icr).and.bb.ne.0) node(j)%rep=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%rep<0) node(i)%rep=0.0;if (node(j)%rep<0) node(j)%rep=0.0;end if                   !!>> HC 12-2-2021
        if (npag(10)>0) then;ii=10;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units !!>> HC 12-2-2021
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%rec; bbb=node(j)%rec  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%rec; bb=nodeo(j)%rec  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%rec=nodeo(i)%rec+a;node(j)%rec=nodeo(j)%rec+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%rec>aaa+aa*(icr).and.aa.ne.0) node(i)%rec=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%rec>bbb+bb*(icr).and.bb.ne.0) node(j)%rec=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%rec<aaa-aa*(icr).and.aa.ne.0) node(i)%rec=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%rec<bbb-bb*(icr).and.bb.ne.0) node(j)%rec=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%rec<0) node(i)%rec=0.0;if (node(j)%rec<0) node(j)%rec=0.0;end if                   !!>> HC 12-2-2021
        if (npag(11)>0) then;ii=11;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units !!>> HC 12-2-2021
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%erp; bbb=node(j)%erp  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%erp; bb=nodeo(j)%erp  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%erp=nodeo(i)%erp+a;node(j)%erp=nodeo(j)%erp+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%erp>aaa+aa*(icr).and.aa.ne.0) node(i)%erp=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%erp>bbb+bb*(icr).and.bb.ne.0) node(j)%erp=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%erp<aaa-aa*(icr).and.aa.ne.0) node(i)%erp=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%erp<bbb-bb*(icr).and.bb.ne.0) node(j)%erp=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%erp<0) node(i)%erp=0.0;if (node(j)%erp<0) node(j)%erp=0.0;end if                   !!>> HC 12-2-2021
        if (npag(12)>0) then;ii=12;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units !!>> HC 12-2-2021
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%est; bbb=node(j)%est  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%est; bb=nodeo(j)%est  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%est=nodeo(i)%est+a;node(j)%est=nodeo(j)%est+b                                          !!>> HC 12-2-2021
       ! RZ gradual 30-5-17                                                                            !!>> HC 12-2-2021
        if (node(i)%est>aaa+aa*(icr).and.aa.ne.0) node(i)%est=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%est>bbb+bb*(icr).and.bb.ne.0) node(j)%est=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%est<aaa-aa*(icr).and.aa.ne.0) node(i)%est=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%est<bbb-bb*(icr).and.bb.ne.0) node(j)%est=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%est<0) node(i)%est=0.0;if (node(j)%est<0) node(j)%est=0.0;end if                   !!>> HC 12-2-2021
        if (npag(13)>0) then;ii=13;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units !!>> HC 12-2-2021 
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually?????                                    !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%eqs; bbb=node(j)%eqs  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%eqs; bb=nodeo(j)%eqs  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%eqs=nodeo(i)%eqs+a;node(j)%eqs=nodeo(j)%eqs+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%eqs>aaa+aa*(icr).and.aa.ne.0) node(i)%eqs=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%eqs>bbb+bb*(icr).and.bb.ne.0) node(j)%eqs=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%eqs<aaa-aa*(icr).and.aa.ne.0) node(i)%eqs=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%eqs<bbb-bb*(icr).and.bb.ne.0) node(j)%eqs=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%eqs<0) node(i)%eqs=0.0;if (node(j)%eqs<0) node(j)%eqs=0.0;end if                   !!>> HC 12-2-2021
        if (npag(14)>0) then;ii=14;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units !!>> HC 12-2-2021 
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                                !!>> HC 12-2-2021
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                          !!>> HC 12-2-2021
        ! Rz 18-5-17 makes node properties to change gradually?????                                    !!>> HC 12-2-2021
        icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
        aaa=node(i)%hoo; bbb=node(j)%hoo  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        aa=nodeo(i)%hoo; bb=nodeo(j)%hoo  ! previous values RZ 30-5-17                                 !!>> HC 12-2-2021
        node(i)%hoo=nodeo(i)%hoo+a;node(j)%hoo=nodeo(j)%hoo+b                                          !!>> HC 12-2-2021
        ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
        if (node(i)%hoo>aaa+aa*(icr).and.aa.ne.0) node(i)%hoo=aaa+aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%hoo>bbb+bb*(icr).and.bb.ne.0) node(j)%hoo=bbb+bb*(icr)                             !!>> HC 12-2-2021
        if (node(i)%hoo<aaa-aa*(icr).and.aa.ne.0) node(i)%hoo=aaa-aa*(icr)                             !!>> HC 12-2-2021
        if (node(j)%hoo<bbb-bb*(icr).and.bb.ne.0) node(j)%hoo=bbb-bb*(icr)                             !!>> HC 12-2-2021
        ! RZ end changes                                                                               !!>> HC 12-2-2021
        if (node(i)%hoo<0) node(i)%hoo=0.0;if (node(j)%hoo<0) node(j)%hoo=0.0;end if                   !!>> HC 12-2-2021
        
        !if (npag(15)>0) then;ii=15;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units!!>> HC 15-2-2021 Epithelial cells have no pseudopodia
        !  if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                             !!>> HC 15-2-2021 so noise related params cannot change
        !  if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                       !!>> HC 15-2-2021 due to genes in these cells
        !  ! Rz 18-5-17 makes node properties to change gradually?????                                 !!>> HC 15-2-2021
        !  icr=icrr ! RZ 30-5-17                                                                       !!>> HC 15-2-2021
        !  aaa=node(i)%mov; bbb=node(j)%mov  ! previous values RZ 30-5-17                              !!>> HC 15-2-2021
        !  aa=nodeo(i)%mov; bb=nodeo(j)%mov  ! previous values RZ 30-5-17                              !!>> HC 15-2-2021
        !  node(i)%mov=nodeo(i)%mov+a;node(j)%mov=nodeo(j)%mov+b                                       !!>> HC 15-2-2021
        !  ! RZ gradual 30-5-17                                                                        !!>> HC 15-2-2021
        !  if (node(i)%mov>aaa+aa*(icr).and.aa.ne.0) node(i)%mov=aaa+aa*(icr)                          !!>> HC 15-2-2021
        !  if (node(j)%mov>bbb+bb*(icr).and.bb.ne.0) node(j)%mov=bbb+bb*(icr)                          !!>> HC 15-2-2021
        !  if (node(i)%mov<aaa-aa*(icr).and.aa.ne.0) node(i)%mov=aaa-aa*(icr)                          !!>> HC 15-2-2021
        !  if (node(j)%mov<bbb-bb*(icr).and.bb.ne.0) node(j)%mov=bbb-bb*(icr)                          !!>> HC 15-2-2021
        !  if (node(i)%mov<0) node(i)%mov=0.0 ; if (node(j)%mov<0) node(j)%mov=0.0                     !!>> HC 15-2-2021
        !end if                                                                                        !!>> HC 15-2-2021

        !if (npag(16)>0) then;ii=16;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units !!>> HC 15-2-2021 
        ! if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                              !!>> HC 15-2-2021
        ! if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif                               !!>> HC 15-2-2021
        ! enddo                                                                                        !!>> HC 15-2-2021
        ! ! Rz 18-5-17 makes node properties to change gradually                                       !!>> HC 15-2-2021
        ! icr=icrr ! RZ 30-5-17                                                                        !!>> HC 15-2-2021
        ! aaa=node(i)%dmo; bbb=node(j)%dmo  ! previous values RZ 30-5-17                               !!>> HC 15-2-2021
        ! aa=nodeo(i)%dmo; bb=nodeo(j)%dmo  ! previous values RZ 30-5-17                               !!>> HC 15-2-2021
        ! node(i)%dmo=nodeo(i)%dmo+a;node(j)%dmo=nodeo(j)%dmo+b                                        !!>> HC 15-2-2021
        ! ! RZ gradual 30-5-17                                                                         !!>> HC 15-2-2021
        ! if (node(i)%dmo>aaa+aa*(icr).and.aa.ne.0) node(i)%dmo=aaa+aa*(icr)                           !!>> HC 15-2-2021
        ! if (node(j)%dmo>bbb+bb*(icr).and.bb.ne.0) node(j)%dmo=bbb+bb*(icr)                           !!>> HC 15-2-2021
        ! if (node(i)%dmo<aaa-aa*(icr).and.aa.ne.0) node(i)%dmo=aaa-aa*(icr)                           !!>> HC 15-2-2021
        ! if (node(j)%dmo<bbb-bb*(icr).and.bb.ne.0) node(j)%dmo=bbb-bb*(icr)                           !!>> HC 15-2-2021
        ! if (node(i)%dmo<0) node(i)%dmo=0.0;if (node(j)%dmo<0) node(j)%dmo=0.0                        !!>> HC 15-2-2021
        !end if                                                                                        !!>> HC 15-2-2021

        if (npag(27)>0) then; ii=27;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units      !!>> HC 12-2-2021
          if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                              !!>> HC 12-2-2021
          if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif;enddo;                        !!>> HC 12-2-2021
          ! Rz 18-5-17 makes node properties to change gradually?????                                  !!>> HC 12-2-2021
          icr=icrr ! RZ 30-5-17                                                                        !!>> HC 12-2-2021
          aaa=node(i)%pla; bbb=node(j)%pla  ! previous values RZ 30-5-17                               !!>> HC 12-2-2021
          aa=nodeo(i)%pla; bb=nodeo(j)%pla  ! previous values RZ 30-5-17                               !!>> HC 12-2-2021
          node(i)%pla=nodeo(i)%pla+a;node(j)%pla=nodeo(j)%pla+b                                        !!>> HC 12-2-2021
          ! RZ gradual 30-5-17                                                                         !!>> HC 12-2-2021
          if (node(i)%pla>aaa+aa*(icr).and.aa.ne.0) node(i)%pla=aaa+aa*(icr)                           !!>> HC 12-2-2021
          if (node(j)%pla>bbb+bb*(icr).and.bb.ne.0) node(j)%pla=bbb+bb*(icr)                           !!>> HC 12-2-2021
          if (node(i)%pla<aaa-aa*(icr).and.aa.ne.0) node(i)%pla=aaa-aa*(icr)                           !!>> HC 12-2-2021
          if (node(j)%pla<bbb-bb*(icr).and.bb.ne.0) node(j)%pla=bbb-bb*(icr)                           !!>> HC 12-2-2021
          if (node(i)%pla<0) node(i)%pla=0.0;if (node(j)%pla<0) node(j)%pla=0.0                        !!>> HC 12-2-2021
        end if                                                                                         !!>> HC 12-2-2021
          
        if (npag(28)>0) then;ii=28;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  !!>> HC 12-2-2021
          if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;                              !!>> HC 12-2-2021
          if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%e(ii) ;endif                               !!>> HC 12-2-2021
          enddo                                                                                        !!>> HC 12-2-2021
          ! Rz 18-5-17 makes node properties to change gradually                                       !!>> HC 12-2-2021
          icr=icrr ! RZ 30-5-17                                                                        !!>> HC 12-2-2021
          aaa=node(i)%kvol; bbb=node(j)%kvol  ! previous values RZ 30-5-17                             !!>> HC 12-2-2021
          aa=nodeo(i)%kvol; bb=nodeo(j)%kvol  ! previous values RZ 30-5-17                             !!>> HC 12-2-2021
          node(i)%kvol=nodeo(i)%kvol+a;node(j)%kvol=nodeo(j)%kvol+b                                    !!>> HC 12-2-2021
          ! RZ gradual 30-5-17                                                                         !!>> HC 12-2-2021
          if (node(i)%kvol>aaa+aa*(icr).and.aa.ne.0) node(i)%kvol=aaa+aa*(icr)                         !!>> HC 12-2-2021
          if (node(j)%kvol>bbb+bb*(icr).and.bb.ne.0) node(j)%kvol=bbb+bb*(icr)                         !!>> HC 12-2-2021
          if (node(i)%kvol<aaa-aa*(icr).and.aa.ne.0) node(i)%kvol=aaa-aa*(icr)                         !!>> HC 12-2-2021
          if (node(j)%kvol<bbb-bb*(icr).and.bb.ne.0) node(j)%kvol=bbb-bb*(icr)                         !!>> HC 12-2-2021
          if (node(i)%kvol<0) node(i)%kvol=0.0;if (node(j)%kvol<0) node(j)%kvol=0.0                    !!>> HC 12-2-2021
        end if                                                                                         !!>> HC 12-2-2021
          
        if(node(i)%erp>10*node(i)%est) node(i)%erp=10*node(i)%est                                      !!>> HC 12-2-2021
        if(node(i)%est>500) node(i)%est=500                                                            !!>> HC 12-2-2021
        if(node(i)%erp>500) node(i)%erp=500                                                            !!>> HC 12-2-2021
        if(node(i)%est<1) node(i)%est=1                                                                !!>> HC 12-2-2021
        if(node(i)%erp<1) node(i)%erp=1                                                                !!>> HC 12-2-2021
        if(node(j)%erp>10*node(j)%est) node(j)%erp=10*node(j)%est                                      !!>> HC 12-2-2021
        if(node(j)%est>500) node(j)%est=500                                                            !!>> HC 12-2-2021
        if(node(j)%erp>500) node(j)%erp=500                                                            !!>> HC 12-2-2021
        if(node(j)%est<1) node(j)%est=1                                                                !!>> HC 12-2-2021
        if(node(j)%erp<1) node(j)%erp=1                                                                !!>> HC 12-2-2021
          
     else if (node(i)%tipus>2)then !this for mesenchyme and ECM ! NEW NEW NEW                          !!>> HC 12-2-2021

        if (npag(21)>0) then                                                                           !!>> HC 12-2-2021
           ii=21 ; a=0                                                                                 !!>> HC 12-2-2021
           do k=1,npag(ii) ;                                                                           !!>> HC 12-2-2021
              kk=whonpag(ii,k)                                                                         !!>> HC 12-2-2021
              if (gex(i,kk)>0.0d0) then ;                                                              !!>> HC 12-2-2021
                 a=a+gex(i,kk)*gen(kk)%e(ii)                                                           !!>> HC 12-2-2021
              endif                                                                                    !!>> HC 12-2-2021
           enddo                                                                                       !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                      !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                       !!>> HC 12-2-2021
           aaa=node(i)%cod                                                                             !!>> HC 12-2-2021
           aa=nodeo(i)%cod  ! previous values RZ 30-5-17                                               !!>> HC 12-2-2021
           node(i)%cod=nodeo(i)%cod+a !>>>Miquel 9-9-15 !pfh no differ 26-11-15                        !!>> HC 12-2-2021
           ! RZ gradual 30-5-17                                                                        !!>> HC 12-2-2021
           if (node(i)%cod>aaa+aa*(icr).and.aa.ne.0) node(i)%cod=aaa+aa*(icr)                          !!>> HC 12-2-2021
           if (node(i)%cod<aaa-aa*(icr).and.aa.ne.0) node(i)%cod=aaa-aa*(icr)                          !!>> HC 12-2-2021
        else                                                                                           !!>> HC 12-2-2021
           node(i)%cod=0                                                                               !!>> HC 12-2-2021
        end if                                                                                         !!>> HC 12-2-2021

        a=nodeo(i)%add-node(i)%eqd !>>>Miquel9-9-15                                                   !!>> HC 12-2-2021
        node(i)%eqd=node(i)%grd+node(i)%cod  !now req is the sum of the req components: growth/apoptosis and contraction/deformation  !!>> HC 12-2-2021
        if(node(i)%eqd>df_reqmax) node(i)%eqd=df_reqmax !put an upper an lower boundary on how much    !>>Miquel28-7-14                 !!>> HC 12-2-2021
        if(node(i)%eqd<reqmin) node(i)%eqd=reqmin !the req can be deformed !NEW *5 deleted                !!>> HC 12-2-2021
        node(i)%add=node(i)%eqd+a !>>>Miquel9-9-15                                                        !!>> HC 15-2-2021 it is node not nodeo
!        if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15                       !!>> HC 15-2-2021
!        if(node(i)%add.ge.node(i)%eqd*1.2) node(i)%add=node(i)%eqd*1.2  !pfh 9-6-15                       !!>> HC 15-2-2021
        !pfh no differ 26-11-15                                                                           !!>> HC 15-2-2021
        if (npag(6)>0) then;ii=6;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ;   !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%add                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%add  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%add=nodeo(i)%add+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%add>aaa+aa*(icr).and.aa.ne.0) node(i)%add=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%add<aaa-aa*(icr).and.aa.ne.0) node(i)%add=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%add<0) node(i)%add=0.0                                                             !!>> HC 15-2-2021
!           if(node(i)%add.le.node(i)%eqd*1.0) node(i)%add=node(i)%eqd*1.0  !pfh 9-6-15                    !!>> HC 15-2-2021
!           if(node(i)%add.ge.node(i)%eqd*1.2) node(i)%add=node(i)%eqd*1.2  !pfh 9-6-15                    !!>> HC 15-2-2021
        else                                                                                              !!>> HC 15-2-2021
           node(i)%add=nodeo(i)%add !pfh 9-6-15                                                           !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(7)>0) then;ii=7;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ;   !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually?????                                    !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%you                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%you  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%you=nodeo(i)%you+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%you>aaa+aa*(icr).and.aa.ne.0) node(i)%you=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%you<aaa-aa*(icr).and.aa.ne.0) node(i)%you=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%you<0) node(i)%you=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(8)>0) then;ii=8;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ;   !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%adh                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%adh  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%adh=nodeo(i)%adh+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%adh>aaa+aa*(icr).and.aa.ne.0) node(i)%adh=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%adh<aaa-aa*(icr).and.aa.ne.0) node(i)%adh=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%adh<0) node(i)%adh=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(9)>0) then;ii=9;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ;   !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%rep                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%rep  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%rep=nodeo(i)%rep+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%rep>aaa+aa*(icr).and.aa.ne.0) node(i)%rep=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%rep<aaa-aa*(icr).and.aa.ne.0) node(i)%rep=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%rep<0) node(i)%rep=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
        
        if (npag(10)>0) then;ii=10;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then;  !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%rec                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%rec  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%rec=nodeo(i)%rec+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%rec>aaa+aa*(icr).and.aa.ne.0) node(i)%rec=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%rec<aaa-aa*(icr).and.aa.ne.0) node(i)%rec=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%rec<0) node(i)%rec=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(11)>0) then;ii=11;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%erp                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%erp  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%erp=nodeo(i)%erp+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%erp>aaa+aa*(icr).and.aa.ne.0) node(i)%erp=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%erp<aaa-aa*(icr).and.aa.ne.0) node(i)%erp=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%erp<0) node(i)%erp=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
   
        if (npag(12)>0) then;ii=12;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%est                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%est  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%est=nodeo(i)%est+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%est>aaa+aa*(icr).and.aa.ne.0) node(i)%est=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%est<aaa-aa*(icr).and.aa.ne.0) node(i)%est=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%est<0) node(i)%est=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(13)>0) then;ii=13;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%eqs                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%eqs  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%eqs=nodeo(i)%eqs+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%eqs>aaa+aa*(delta).and.aa.ne.0) node(i)%eqs=aaa+aa*(delta)                         !!>> HC 15-2-2021
           if (node(i)%eqs<aaa-aa*(delta).and.aa.ne.0) node(i)%eqs=aaa-aa*(delta)                         !!>> HC 15-2-2021
           if (node(i)%eqs<0) node(i)%eqs=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
      
        if (npag(14)>0) then;ii=14;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 15-2-2021 
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 15-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 15-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 15-2-2021
           aaa=node(i)%hoo                                                                                !!>> HC 15-2-2021
           aa=nodeo(i)%hoo  ! previous values RZ 30-5-17                                                  !!>> HC 15-2-2021
           node(i)%hoo=nodeo(i)%hoo+a                                                                     !!>> HC 15-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 15-2-2021
           if (node(i)%hoo>aaa+aa*(icr).and.aa.ne.0) node(i)%hoo=aaa+aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%hoo<aaa-aa*(icr).and.aa.ne.0) node(i)%hoo=aaa-aa*(icr)                             !!>> HC 15-2-2021
           if (node(i)%hoo<0) node(i)%hoo=0.0                                                             !!>> HC 15-2-2021
        end if                                                                                            !!>> HC 15-2-2021
  
        if (npag(15)>0) then;ii=15;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 12-2-2021
           a=a+gex(i,kk)*gen(kk)%e(ii) ;endif;enddo                                                       !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
           aaa=node(i)%mov                                                                                !!>> HC 12-2-2021
           aa=nodeo(i)%mov  ! previous values RZ 30-5-17                                                  !!>> HC 12-2-2021
           node(i)%mov=nodeo(i)%mov+a                                                                     !!>> HC 12-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
           if (node(i)%mov>aaa+aa*(icr).and.aa.ne.0) node(i)%mov=aaa+aa*(icr)                             !!>> HC 12-2-2021
           if (node(i)%mov<aaa-aa*(icr).and.aa.ne.0) node(i)%mov=aaa-aa*(icr)                             !!>> HC 12-2-2021
           if (node(i)%mov<0) node(i)%mov=0.0                                                             !!>> HC 12-2-2021
           else; if(node(i)%mov==10) node(i)%mov=1000 !!!!!!                                              !!>> HC 12-2-2021
        end if                                                                                            !!>> HC 12-2-2021
      
        if (npag(16)>0) then;ii=16;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; !!>> HC 12-2-2021
           a=a+gex(i,kk)*gen(kk)%e(ii);endif;enddo                                                        !!>> HC 12-2-2021
           ! Rz 18-5-17 makes node properties to change gradually                                         !!>> HC 12-2-2021
           icr=icrr ! RZ 30-5-17                                                                          !!>> HC 12-2-2021
           aaa=node(i)%dmo                                                                                !!>> HC 12-2-2021
           aa=nodeo(i)%dmo  ! previous values RZ 30-5-17                                                  !!>> HC 12-2-2021
           node(i)%dmo=nodeo(i)%dmo+a                                                                     !!>> HC 12-2-2021
           ! RZ gradual 30-5-17                                                                           !!>> HC 12-2-2021
           if (node(i)%dmo>aaa+aa*(icr).and.aa.ne.0) node(i)%dmo=aaa+aa*(icr)                             !!>> HC 12-2-2021
           if (node(i)%dmo<aaa-aa*(icr).and.aa.ne.0) node(i)%dmo=aaa-aa*(icr)                             !!>> HC 12-2-2021
           if (node(i)%dmo<0) node(i)%dmo=0.0                                                             !!>> HC 12-2-2021
        end if                                                                                            !!>> HC 12-2-2021
     end if                                                                                               !!>> HC 12-2-2021
  end do                                                                                                  !!>> HC 12-2-2021

  if (ffu(21)==1)then       !!>> HC 17-11-2020  THIS FFU CONTROLS WHETHER WE ARE APPLYING FILTERS OR NOT   
     call filter            !!>> HC 12-2-2021                                                                
  endif                     !!>> HC 12-2-2021

end subroutine nexe_gradual




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FILTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine filter  !!>> HC 12-2-2021 This subroutine applies filters if the corresponding ffu(21) is ON
  implicit none    !!>> HC 12-2-2021
  real*8 mx1,my1,my2,mz1,mz3,mx2,my,mz2,mz4                                    !!>> HC 17-9-2020
  real*8 :: vvv, uuu, all_dmov, all_dmov_max, alld, cputime, difk              !!>> HC 18-9-2020
  integer ::  lch, iech, occupiedch                                            !!>> HC 28-1-2021 
  real*8 :: nodes_per_box                                                      !!>> HC 28-1-2021 
  real*8 :: scui,maxa,maxs(3),mins(3)                                          !!>> HC 18-9-2020
  real*8 :: addoch, minaddch, maxaddch                                         !!>> HC 14-1-2021
  integer :: ich, l, jch, nd_a, nd_t, ndm, first=0                             !!>> HC 26-2-2021
  character*400 cx,rangfile,onetwin                                            !!>> HC 12-2-2021
  integer :: ntotbroken, nepi, ic, ica                                         !!>> TT 5-3-2021
  real*8 :: uvx, uvy, uvz, mvecx, mvecy, mvecz, sumdotp, di, threshdotp, dc    !!>> TT 5-3-2021
  real*8, allocatable :: sumdotps(:)                                           !!>> HC 11-11-2021
  integer :: limit                                                             !!>> HC 11-11-2021
  real*8, dimension(1:nga) :: max_elim, min_elim                               !!>> HC 11-11-2020 These vectors store the limits of e matrix 
  integer, dimension(1:nga) :: rembeh                                          !!>> HC 11-11-2021 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  real*8 :: wrongval, maxl, minl,ratio                                         !!>> HC 11-11-2021
  real*8, dimension (1:5) :: min_glim, max_glim                                !!>> HC 11-11-2021
  integer :: sumneigh, many, oripich                                           !!>> HC 14-12-2021
  real*8 :: percent                                                            !!>> HC 14-12-2021

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Calculating displacements!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020

     uuu=0d0 ; vvv=1d0                                                                                                   !!>> HC 17-9-2020
     if(tellme2==0)then ! store at the first iteration. This could also be done in a different way. !tellme2 en general  !!>> HC 17-9-2020
       all_dmov=1d0; all_dmov_max=0d0                                                                                    !!>> HC 17-9-2020
       tellme2=getot                                                                                                     !!>> HC 17-9-2020
       if(allocated(nodeu)) deallocate(nodeu)                                                                            !!>> HC 17-9-2020
       allocate(nodeu(nd,3))                                                                                             !!>> HC 17-9-2020
       do i=1,nd                                                                                                         !!>> HC 17-9-2020
         nodeu(i,1)=nodeo(i)%x;nodeu(i,2)=nodeo(i)%y;nodeu(i,3)=nodeo(i)%z                                               !!>> HC 17-9-2020
       end do                                                                                                            !!>> HC 17-9-2020
       ndu=nd                                                                                                            !!>> HC 17-9-2020
     end if                                                                                                              !!>> HC 17-9-2020

     if((mod(int(rtime/deltamax),10)==0).and.(getot>tellme2))then  ! check every 10 iterations ARBITRARY                 !!>> HC 17-9-2020
       all_dmov=0d0 ; difk=0 ; alld=0d0                                                                                  !!>> HC 17-9-2020
       do l=1,ndu                                                                                                        !!>> HC 17-9-2020
  
         difk=difk+node(l)%dif                                                                                           !!>> HC 17-9-2020
         if(node(l)%dif<vvv) vvv=node(l)%dif                                                                             !!>> HC 17-9-2020
         if(node(l)%dif>uuu) uuu=node(l)%dif                                                                             !!>> HC 17-9-2020
         scui=sqrt((node(l)%x-nodeo(l)%x)**2+(node(l)%y-nodeo(l)%y)**2+(node(l)%z-nodeo(l)%z)**2)+scui ! total displacement !!>> HC 17-9-2020
         if(node(l)%x<mins(1)) mins(1)=node(l)%x  ! store extreme points                                                 !!>> HC 17-9-2020
         if(node(l)%x>maxs(1)) maxs(1)=node(l)%x                                                                         !!>> HC 17-9-2020
         if(node(l)%y<mins(2)) mins(2)=node(l)%y                                                                         !!>> HC 17-9-2020
         if(node(l)%y>maxs(2)) maxs(2)=node(l)%y                                                                         !!>> HC 17-9-2020
         if(node(l)%z<mins(3)) mins(3)=node(l)%z                                                                         !!>> HC 17-9-2020
         if(node(l)%z>maxs(3)) maxs(3)=node(l)%z                                                                         !!>> HC 17-9-2020

         if(getot>10)then                                                                                                !!>> HC 17-9-2020
           all_dmov=all_dmov+sqrt((node(l)%x-nodeu(l,1))**2+(node(l)%y-nodeu(l,2))**2+(node(l)%z-nodeu(l,3))**2)         !!>> HC 17-9-2020
           if(all_dmov>all_dmov_max) all_dmov_max=all_dmov                                                               !!>> HC 17-9-2020
           if(alld<sqrt((node(l)%x-nodeu(l,1))**2+(node(l)%y-nodeu(l,2))**2+(node(l)%z-nodeu(l,3))**2))then              !!>> HC 17-9-2020
             alld=sqrt((node(l)%x-nodeu(l,1))**2+(node(l)%y-nodeu(l,2))**2+(node(l)%z-nodeu(l,3))**2)                    !!>> HC 17-9-2020
           endif                                                                                                         !!>> HC 17-9-2020
         else                                                                                                            !!>> HC 17-9-2020
           all_dmov=0.0001 ! we never end up there                                                                       !!>> HC 17-9-2020
         endif                                                                                                           !!>> HC 17-9-2020
       end do

       if(allocated(nodeu)) deallocate(nodeu)                                                                            !!>> HC 17-9-2020
       allocate(nodeu(nd,3))                                                                                             !!>> HC 17-9-2020
       do i=1,nd                                                                                                         !!>> HC 17-9-2020
         nodeu(i,1)=node(i)%x;nodeu(i,2)=node(i)%y;nodeu(i,3)=node(i)%z  ! update positions saved                        !!>> HC 17-9-2020
       end do                                                                                                            !!>> HC 17-9-2020
       ndu=nd                                                                                                            !!>> HC 17-9-2020

       maxa=maxs(1)-mins(1)                                                                                              !!>> HC 17-9-2020
       if(maxs(2)-mins(2)>maxa) maxa=maxs(2)-mins(2)                                                                     !!>> HC 17-9-2020
       if(maxs(3)-mins(3)>maxa) maxa=maxs(3)-mins(3)                                                                     !!>> HC 17-9-2020
       scui=scui/(maxa*nd)  ! total displacement, not in use                                                             !!>> HC 17-9-2020

     end if                                                                                                              !!>> HC 17-9-2020
 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FILTERING WRONG MORPHOLOGIES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     
     ! WHICHEND 1: Check cpu_time                    !!>> HC 2-7-2020  THIS COMES FROM RENSKE
     if (ffufi(1,1)==1.and.ffufi(1,3)>0)then         !!>> HC 20-2-2021
        if (mod(getot,ffufi(1,3)).eq.0)then          !!>> HC 20-2-2021
           call cpu_time(cputime)                    !!>> HC 2-7-2020
           if ((cputime-start_time)>14400)then       !!>> HC 8-1-2021
              whichend=1                             !!>> HC 2-7-2020
              goto 171                               !!>> HC 14-8-2020
           end if                                    !!>> HC 2-7-2020
        endif                                        !!>> HC 20-2-2021
     endif                                           !!>> HC 20-2-2021
     
     ! WHICHEND 2: TOO MANY NEIGHBORS IT IS IN BIOMECHANIC.MOD.F90

     ! WHICHEND 3: Check total expansion of morphology: a diameter of 50LE is arbitrary, but filters out explosions         !!>> HC 17-9-2020
     if (ffufi(3,1)==1.and.ffufi(3,3)>0)then                                                                                !!>> HC 20-2-2021
        if (mod(getot,ffufi(3,3)).eq.0)then                                                                                 !!>> HC 20-2-2021
           mx1=1000d0;my1=1000d0;mz1=1000d0;mz3=1000d0;mx2=-1000d0;my2=-1000d0;mz2=-1000d0;mz4=-1000d0;ndm=0;k=1;l=1        !!>> HC 17-9-2020
           do i=1,nd                                                                                                        !!>> HC 17-9-2020
              if (node(i)%tipus == 2 .or. node(i)%tipus == 3 .or. node(i)%tipus == 4) then  !!>>> HC 2-6-2020               !!>> HC 17-9-2020
                 if((node(i)%z.lt.mz3).and.(node(i)%tipus<3))then; mz3=node(i)%z; k=i; endif                                !!>> HC 17-9-2020
                 if((node(i)%z.gt.mz4).and.(node(i)%tipus<3))then; mz4=node(i)%z; l=i; endif                                !!>> HC 17-9-2020
                 if(node(i)%x.lt.mx1) mx1=node(i)%x                                                                         !!>> HC 17-9-2020
                 if(node(i)%y.lt.my1) my1=node(i)%y                                                                         !!>> HC 17-9-2020
                 if(node(i)%z.lt.mz1) mz1=node(i)%z                                                                         !!>> HC 17-9-2020
                 if(node(i)%x.gt.mx2) mx2=node(i)%x                                                                         !!>> HC 17-9-2020
                 if(node(i)%y.gt.my2) my2=node(i)%y                                                                         !!>> HC 17-9-2020
                 if(node(i)%z.gt.mz2) mz2=node(i)%z                                                                         !!>> HC 17-9-2020
                 if(node(i)%tipus==3) ndm=ndm+1                                                                             !!>> HC 17-9-2020
              endif !!>>> HC 2-6-2020                                                                                       !!>> HC 17-9-2020
           end do                                                                                                           !!>> HC 17-9-2020
           if ((abs(mx1-mx2).gt.50d0).or.(abs(my1-my2).gt.50d0).or.(abs(mz1-mz2).gt.50d0))then                              !!>> HC 17-9-2020
              whichend=3                                                                                                    !!>> HC 17-9-2020
              goto 171 ! OUT: too large                                                                                     !!>> HC 17-9-2020
           end if                                                                                                           !!>> HC 17-9-2020
        endif                                                                                                               !!>> HC 20-2-2021
     endif                                                                                                                  !!>> HC 20-2-2021


     ! WHICHEND 4 and 5: TOO MANY OR NO NODES IN MODEL.MOD.F90

     ! WHICHEND 6: NO MOVEMENT
     ! Check movement at a late time                                                                                        !!>> HC 17-9-2020
     if (ffufi(6,1)==1.and.ffufi(6,3)>0)then                                                                                !!>> HC 20-2-2021
        if (mod(getot,ffufi(6,3)).eq.0)then                                                                                 !!>> HC 20-2-2021
           if (rtime>100)then ! New: check if movement has gone down                                                        !!>> HC 17-9-2020
              if (all_dmov<0.2*all_dmov_max)then                                                                            !!>> HC 17-9-2020
                 whichend=6                                                                                                 !!>> HC 17-9-2020
                 goto 171  ! THIS ONE IS FUNKY, BUT RARELY APPLIES                                                          !!>> HC 17-9-2020
              endif                                                                                                         !!>> HC 17-9-2020
           end if                                                                                                           !!>> HC 17-9-2020
        endif                                                                                                               !!>> HC 20-2-2021
     endif                                                                                                                  !!>> HC 20-2-2021
     
     
     
     ! WHICHEND 7: TIME IS OVER (GETOT)
     if (ffufi(7,1)==1.and.ffufi(7,3)>0)then                                                                                !!>> HC 20-2-2021
        if (mod(getot,ffufi(7,3)).eq.0)then                                                                                 !!>> HC 20-2-2021     
           if (getot>5000000)then ! TIME OVER!!! Make sure automaticon.mod doesnt stop it before that                       !!>> HC 17-9-2020
              whichend=7                                                                                                    !!>> HC 17-9-2020
              goto 171  ! WE NEVER GET THERE ACTUALLY                                                                       !!>> HC 17-9-2020
           end if                                                                                                           !!>> HC 17-9-2020
        endif                                                                                                               !!>> HC 20-2-2021
     endif                                                                                                                  !!>> HC 20-2-2021

     ! WHICHEND 8: TIME IS OVER (realtime)     
     if (ffufi(8,1)==1.and.ffufi(8,3)>0)then                                                                                !!>> HC 20-2-2021
        if (mod(getot,ffufi(8,3)).eq.0)then                                                                                 !!>> HC 20-2-2021         
           if (rtime>5000)then ! TIME OVER!!! Make sure automaticon.mod doesnt stop it before that                          !!>> HC 17-9-2020
              whichend=8                                                                                                    !!>> HC 17-9-2020
              goto 171  ! OUT:SAVE NEVER FINISHED, should be regarded failure: ca.2%                                        !!>> HC 17-9-2020
           end if                                                                                                           !!>> HC 17-9-2020
        endif                                                                                                               !!>> HC 20-2-2021
     endif                                                                                                                  !!>> HC 20-2-2021


     ! WHICHEND 9: 
     ! if a node gets 3 times the original size of node 1     !!>>HC 14-1-2021
     ! and it is also 3 times larger than the smallest node   !!>>HC 14-1-2021
     ! we fucking kill the program                            !!>>HC 14-1-2021
     if (ffufi(9,1)==1.and.ffufi(9,3)>0)then                  !!>> HC 20-2-2021
        if (mod(getot,ffufi(9,3)).eq.0)then                   !!>> HC 20-2-2021
           addoch=nodeo(1)%add                                !!>>HC 14-1-2021 original value (we assume all nodes had more of less the same eqd)
           maxaddch=maxval(node(:nd)%add)                     !!>>HC 14-1-2021 maximum current value
           minaddch=minval(node(:nd)%add)                     !!>>HC 14-1-2021 minimum current value
           if (3*maxaddch>addoch)then                         !!>>HC 14-1-2021 if the node is more than 3 times the original value
              if (3*minaddch<maxaddch)then                    !!>>HC 14-1-2021 if the node is more than 3 times larger than the smallest node
                 whichend=9                                   !!>>HC 14-1-2021 we kill the program
                 goto 171                                     !!>>HC 14-1-2021
              endif                                           !!>>HC 14-1-2021
           endif                                              !!>>HC 14-1-2021
        endif                                                 !!>> HC 20-2-2021
     endif                                                    !!>> HC 20-2-2021

     ! WHICHEND 10: AVERAGE NUMBER OF NODES PER BOX
     if (ffufi(10,1)==1.and.ffufi(10,3)>0)then                !!>> HC 20-2-2021
        if (mod(getot,ffufi(10,3)).eq.0)then                  !!>> HC 20-2-2021 This filter calculates the average number of nodes 
           nodes_per_box=0.0d0; occupiedch=0                  !!>> HC 28-1-2021 per occupied box
           do ich=-nboxes, nboxes                             !!>> HC 28-1-2021
              do jch=-nboxes, nboxes                          !!>> HC 28-1-2021
                 do lch=-nboxes, nboxes                       !!>> HC 28-1-2021
                    iech=boxes(ich,jch,lch)                   !!>> HC 28-1-2021                     
                    if (iech>0) occupiedch=occupiedch+1       !!>> HC 28-1-2021 This box is occupied
                 enddo                                        !!>> HC 28-1-2021
              enddo                                           !!>> HC 28-1-2021
           enddo                                              !!>> HC 28-1-2021
           nodes_per_box=nd/real(occupiedch)                  !!>> HC 28-1-2021 This is the average number of nodes per box
           if (nodes_per_box.gt.maxbox)then                   !!>> HC 28-1-2021 If it is larger then maxbox, we kill the program
              whichend=10                                     !!>> HC 28-1-2021
              if (ffu(22)==0) print*, "TOO CROWDED: THE AVERAGE NUMBER OF NODES PER BOX IS TOO BIG"  !!>> HC 28-1-2021
              goto 171                                        !!>> HC 28-1-2021
           endif                                              !!>> HC 28-1-2021
        endif                                                 !!>> HC 28-1-2021
     endif                                                    !!>> HC 20-2-2021

     ! WHICHEND 11: NO MOVEMENT
     if (ffufi(11,1)==1.and.ffufi(11,3)>0)then                                                                                     !!>> HC 20-2-2021
        if (mod(getot,ffufi(11,3)).eq.0)then                                                                                       !!>> HC 20-2-2021   
           if ((rtime>5).and.(rtime<6).and.(all_dmov<nd*delta))then ! ca.500 iterations before discarding ! OUT: NO MORE MOVEMENTS !!>> HC 17-9-2020
              whichend=11                                                                                                          !!>> HC 17-9-2020
              print*, "NO MORE MOVEMENT", all_dmov; goto 171 ! most morphologies stop here                                         !!>> HC 17-9-2020
           endif
        endif
     endif       



     !!!!!!!!!!!!!!!!!!!!!!!THIS COMES FROM ENSEMBLE MODE (HAGOLANI 2018)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>> HC 17-9-2020
     ! WHICHEND 12: BROKEN MORPHOLOGIES   !!!!!!!!!WE HAVE TO CHECK THESE!!!!!!!!!!!!!!!!!!
     if (ffufi(12,1)==1.and.ffufi(12,3)>0)then                        !!>> HC 14-12-2021
        if (mod(getot,ffufi(12,3)).eq.0)then                          !!>> HC 14-12-2021
           oripich=0; percent=0.0d0; many=0                           !!>> HC 24-9-2021 Here we calculate the original number of epithelial nodes
           do ich=1,nd                                                !!>> HC 24-9-2021
              if (nodeo(ich)%tipus.ge.3)cycle                         !!>> HC 24-9-2021 WARNING 1 we are assuming that the epitelial cells are listed BEFORE the mesenchyme
              if (nodeo(ich)%icel>oripich) oripich=nodeo(ich)%icel    !!>> HC 24-9-2021 the original number of epithelial nodes is then the maximum number of epitelial cells
           enddo                                                      !!>> HC 24-9-2021 in the nodeo list (WARNING 2 assuming single_node cells)
           oripich=oripich*2                                          !!>> HC 14-12-2021 original number of epitelial nodes        
           do ich=1,nd                                                !!>> HC 14-12-2021 We llok all the epithelial nodes
              if(node(ich)%tipus.ge.3)cycle                           !!>> HC 14-12-2021
              sumneigh=0                                              !!>> HC 14-12-2021
              do jch=1,nneigh(ich)                                    !!>> HC 14-12-2021 Count the neighbors that are epithelial too
                 if(node(neigh(ich,jch))%tipus.ge.3)cycle             !!>> HC 14-12-2021
                 sumneigh=sumneigh+1                                  !!>> HC 14-12-2021
              enddo                                                   !!>> HC 14-12-2021
              if(sumneigh<3) many=many+1                              !!>> HC 14-12-2021 We consider the nodes with less than 3 neighbors broken
           enddo                                                      !!>> HC 14-12-2021 (empirically calculated)
           percent=real(many)/real(oripich)                           !!>> HC 14-12-2021
           if (percent>0.20d0)then                                    !!>> HC 14-12-2021 This is an empirically calculated limit
              whichend=12                                             !!>> HC 14-12-2021 Limit surpased, we kill this simulation
              goto 171                                                !!>> HC 14-12-2021
           endif                                                      !!>> HC 14-12-2021
        endif                                                         !!>> HC 14-12-2021       
     endif                                                            !!>> HC 14-12-2021
!     if (ffufi(12,1)==1.and.ffufi(12,3)>0)then                                                                  !!>> HC 20-2-2021
!        if (mod(getot,ffufi(12,3)).eq.0)then                                                                    !!>> HC 12-7-2021
!           threshdotp=4                                                                                         !!>> TT 5-3-2021
!           ntotbroken=0; nepi=0                                                                                 !!>> TT 5-3-2021
!           do i=1, nd                                                                                           !!>> TT 5-3-2021
!               if(node(i)%tipus .gt. 2)cycle                                                                    !!>> TT 5-3-2021
!               nepi=nepi+1                                                                                      !!>> TT 5-3-2021
!               sumdotp=0                                                                                        !!>> TT 5-3-2021
!               uvx=node(node(i)%altre)%x-node(i)%x                                                              !!>> TT 5-3-2021
!               uvy=node(node(i)%altre)%y-node(i)%y                                                              !!>> TT 5-3-2021
!               uvz=node(node(i)%altre)%z-node(i)%z                                                              !!>> TT 5-3-2021
!               di=sqrt(uvx**2+uvy**2+uvz**2)                                                                    !!>> TT 5-3-2021
!               do j=1, nneigh(i)                                                                                !!>> HC 12-7-2021
!                   ic=neigh(i,j)                                                                                !!>> HC 5-3-2021
!                   if(node(ic)%tipus .ne. node(i)%tipus)cycle                                                   !!>> TT 5-3-2021
!                   mvecx=node(ic)%x-node(i)%x                                                                   !!>> TT 5-3-2021
!                   mvecy=node(ic)%y-node(i)%y                                                                   !!>> TT 5-3-2021
!                   mvecz=node(ic)%z-node(i)%z                                                                   !!>> TT 5-3-2021
!                   sumdotp=sumdotp+(mvecx*uvx+mvecy*uvy+mvecz*uvz)/di                                           !!>> TT 5-3-2021
!                   ica=node(ic)%altre                                                                           !!>> TT 5-3-2021
!                   mvecx=node(ica)%x-node(i)%x                                                                  !!>> TT 5-3-2021
!                   mvecy=node(ica)%y-node(i)%y                                                                  !!>> TT 5-3-2021
!                   mvecz=node(ica)%z-node(i)%z                                                                  !!>> TT 5-3-2021
!                   sumdotp=sumdotp+(mvecx*uvx+mvecy*uvy+mvecz*uvz)/di                                           !!>> TT 5-3-2021
!               end do                                                                                           !!>> TT 5-3-2021
!               if(abs(sumdotp).gt.threshdotp)ntotbroken=ntotbroken+1                                            !!>> TT 5-3-2021
!           end do                                                                                               !!>> TT 5-3-2021
!           if (ntotbroken > 0.05*nepi)then                                                                      !!>> TT 5-3-2021
!              whichend=12                                                                                            !!>> HC 17-9-2020
!              goto 171                                                                                               !!>> HC 17-9-2020
!           endif                                                                                                    !!>> HC 17-9-2020
!        endif
!     endif                                                                                                      !!>> HC 17-9-2020
  
     !!WHICHEND 13: BLACK HOLES
     if (ffufi(13,1)==1.and.ffufi(13,3)>0)then                                                                  !!>> HC 20-2-2021
        if (mod(getot,ffufi(13,3)).eq.0)then                                                                    !!>> HC 20-2-2021 
           c=0; cc=0                                                                                                !!>> HC 17-9-2020
           do i=1,nd                                                                                                !!>> HC 17-9-2020
              if(node(i)%tipus.ne.1) cycle                                                                           !!>> HC 17-9-2020
              b=0                                                                                                    !!>> HC 17-9-2020
              l=0; aa=1000000d0                                                                                      !!>> HC 17-9-2020
              do j=1,nneigh(i)                                                                                       !!>> HC 17-9-2020
                 k=neigh(i,j)                                                                                         !!>> HC 17-9-2020
                if(k==i) cycle                                                                                       !!>> HC 17-9-2020
                a=sqrt((node(i)%x-node(k)%x)**2+(node(i)%y-node(k)%y)**2+(node(i)%z-node(k)%z)**2)                   !!>> HC 17-9-2020
                if(a<0.05*node(i)%eqd) b=b+1                                                                         !!>> HC 17-9-2020
                if((a<aa).and.(node(i)%altre.ne.k).and.(node(k)%tipus<3))then; aa=a; l=k; endif                      !!>> HC 17-9-2020
              end do                                                                                                 !!>> HC 17-9-2020
              if(l==0) cycle                                                                                         !!>> HC 17-9-2020
              if(node(l)%tipus.ne.node(i)%tipus) cc=1+cc ! broken epithelium                                         !!>> HC 17-9-2020
              if(b>1) c=c+1 ! black hole                                                                             !!>> HC 17-9-2020
           end do                                                                                                   !!>> HC 17-9-2020
           !if(c>0.01*ndoo)then                                                                                     !!>> HC 17-9-2020
           if (cc>0.8*ndoo)then !>>> TT 05-06-2020 Originally 0.1*ndoo, made a crash in initial condition 5            !!>> HC 17-9-2020
              whichend=13                                                                                            !!>> HC 17-9-2020
              goto 171                                                                                               !!>> HC 17-9-2020
           endif                                                                                                    !!>> HC 17-9-2020
        endif
     endif                                                                                                      !!>> HC 17-9-2020
     
     !!WHICHEND 14: BLACK HOLES
     if (ffufi(14,1)==1.and.ffufi(14,3)>0)then                                                                  !!>> HC 20-2-2021
        if (mod(getot,ffufi(14,3)).eq.0)then                                                                    !!>> HC 20-2-2021 
           rangfile="ranges.dat"
           call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh) !!>> HC 6-10-2021 Read ranges
             
           limit=0                                               !!>> HC 27-11-2020 if limit stays in 0 no limits have been surpased if =0, we will the embryo
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
           do ich=1,nd                                                                                !!>> HC 28-9-2021 Here we check whether any node has surpased the limits
              !! 6 ADD                                                                                !!>> HC 28-9-2021 of node properties that correspond to the maximum/minimum activation
              ratio=node(ich)%add/node(ich)%eqd                                                       !!>> HC 28-9-2021 in a gene with concentration = 1
              if (ratio>2.30d0)then; limit=1; lch=6; wrongval=ratio;  exit; endif                        !!>> HC 28-9-2021
              if (ratio<1.10d0)then; limit=2; lch=6; wrongval=ratio;  exit; endif                        !!>> HC 28-9-2021 In ADD however, we check the ADD/EHD ration
                                                                                             !!>> HC 28-9-2021
              !! 8 ADH                                                                                !!>> HC 28-9-2021
              maxl=max_elim(8)+nodeo(ich)%adh                                                         !!>> HC 28-9-2021 This is the maximum realistic value of adh
              minl=min_elim(8)+nodeo(ich)%adh                                                         !!>> HC 28-9-2021 This is the minimum realistic value of adh
              if (node(ich)%adh>maxl)then; limit=1; lch=8; wrongval=node(ich)%adh;  exit; endif       !!>> HC 28-9-2021
              if (node(ich)%adh<minl)then; limit=2; lch=8; wrongval=node(ich)%adh;  exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 10 REC                                                                               !!>> HC 28-9-2021 And the same for all the properties
              maxl=max_elim(10)+nodeo(ich)%rec                                                        !!>> HC 28-9-2021
              minl=min_elim(10)+nodeo(ich)%rec                                                        !!>> HC 28-9-2021
              if (node(ich)%rec>maxl)then; limit=1; lch=10; wrongval=node(ich)%rec; exit; endif       !!>> HC 28-9-2021
              if (node(ich)%rec<minl)then; limit=2; lch=10; wrongval=node(ich)%rec; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 11 ERP                                                                               !!>> HC 28-9-2021
              maxl=max_elim(11)+nodeo(ich)%erp                                                        !!>> HC 28-9-2021
              minl=min_elim(11)+nodeo(ich)%erp                                                        !!>> HC 28-9-2021
              if (node(ich)%erp>maxl)then; limit=1; lch=11; wrongval=node(ich)%erp; exit; endif       !!>> HC 28-9-2021
              if (node(ich)%erp<minl)then; limit=2; lch=11; wrongval=node(ich)%erp; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 12 EST                                                                               !!>> HC 28-9-2021
              maxl=max_elim(12)+nodeo(ich)%est                                                        !!>> HC 28-9-2021
              minl=min_elim(12)+nodeo(ich)%est                                                        !!>> HC 28-9-2021
              if (node(ich)%est>maxl)then; limit=1; lch=12; wrongval=node(ich)%est; exit; endif       !!>> HC 28-9-2021
              if (node(ich)%est<minl)then; limit=2; lch=12; wrongval=node(ich)%est; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 13 EQS                                                                               !!>> HC 28-9-2021
              maxl=max_elim(13)+nodeo(ich)%eqs                                                        !!>> HC 28-9-2021
              minl=min_elim(13)+nodeo(ich)%eqs                                                        !!>> HC 28-9-2021
              if (node(ich)%eqs>maxl)then; limit=1; lch=13; wrongval=node(ich)%eqs; exit; endif       !!>> HC 28-9-2021
              if (node(ich)%eqs<minl)then; limit=2; lch=13; wrongval=node(ich)%eqs; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 14 HOO                                                                               !!>> HC 28-9-2021
              maxl=max_elim(14)+nodeo(ich)%hoo                                                        !!>> HC 28-9-2021
              minl=min_elim(14)+nodeo(ich)%hoo                                                        !!>> HC 28-9-2021
              if (node(ich)%hoo>maxl)then; limit=1; lch=14; wrongval=node(ich)%hoo; exit; endif       !!>> HC 28-9-2021
              if (node(ich)%hoo<minl)then; limit=2; lch=14; wrongval=node(ich)%hoo; exit; endif       !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 21 COD                                                                               !!>> HC 28-9-2021
              maxl=max_elim(21)+nodeo(ich)%cod                                                        !!>> HC 28-9-2021
              minl=min_elim(21)+nodeo(ich)%cod                                                        !!>> HC 28-9-2021
              if (node(ich)%cod>maxl)then; limit=1; lch=21; wrongval=node(ich)%cod; exit;  endif      !!>> HC 28-9-2021
              if (node(ich)%cod<minl)then; limit=2; lch=21; wrongval=node(ich)%cod; exit;  endif      !!>> HC 28-9-2021
                                                                                             !!>> HC 28-9-2021
              !! 22 GRD                                                                               !!>> HC 28-9-2021
              maxl=max_elim(22)+nodeo(ich)%grd                                                        !!>> HC 28-9-2021
              minl=min_elim(22)+nodeo(ich)%grd                                                        !!>> HC 28-9-2021
              if (node(ich)%grd>maxl)then; limit=1; lch=22; wrongval=node(ich)%grd; exit;  endif      !!>> HC 28-9-2021
              if (node(ich)%grd<minl)then; limit=2; lch=22; wrongval=node(ich)%grd; exit;  endif      !!>> HC 28-9-2021
           enddo                                                                                      !!>> HC 28-9-2021
  
  
           if (limit>0)then                            !!>> HC 28-9-2021 If limits have been surpased, we kill this individual
              whichend=14                              !!>> HC 11-11-2021
              goto 171                                 !!>> HC 11-11-2021
           endif                                       !!>> HC 28-9-2021
        endif
     endif                                                                                                      !!>> HC 17-9-2020


                                                                                                              !!>> HC 17-9-2020

return !X! THIS IS THE END         !!>> HC 2-7-2020



! REASONS TO STOP A SIMULATION: WHICHEND VALUES          !!>> HC 18-9-2020
!  1 = OVER  4 HOURS RUNNING                             !!>> HC  8-1-2021
!  2 = TOO MANY NEIGHBORS                                !!>> HC  7-1-2021
!  3 = TOO LARGE                                         !!>> HC 18-9-2020
!  4 = TOO MANY NODES                                    !!>> HC  7-1-2021
!  5 = NOT ENOUGH NODES                                  !!>> HC  7-1-2021
!  6 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
!  7 = TIME HAS ENDED (GETOT)                            !!>> HC 18-9-2020
!  8 = TIME HAS ENDED (REALTIME)                         !!>> HC 18-9-2020
!  9 = SOME NODES ARE DISPROPORTIONALY LARGE             !!>> HC 14-1-2020
! 10 = THERE ARE TOO MANY NODES IN AVERAGE PER BOX       !!>> HC 28-1-2021
! 11 = NO MORE MOVEMENT                                  !!>> HC 18-9-2020
! 12 = BROKEN EPITHELIUM                                 !!>> HC 18-9-2020
! 13 = BLACK HOLES                                       !!>> HC 18-9-2020
! 14 = CHECK RANGES                                      !!>> HC 11-11-2021


  
  
171  print *,""                                                                                                 !!>> HC 18-9-2020
     if (ffu(22)==0)then                                                                                           !!>> HC 2-12-2020  Silent mode does not print this
         print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)          !!>> HC 17-11-2020
         print *,""                                                                                                !!>> HC 17-11-2020
         print*,"reason to end:",whichend                                                                          !!>> HC 17-11-2020
     endif                                                                                                         !!>> HC 2-12-2020
     cx=trim(carg)//".loga"                                                                                        !!>> HC 19-2-2021
     open(616, file=trim(cx))                                                                                      !!>> HC 19-2-2021 We want to know
         if (whichend==14)then                                                                                     !!>> HC 11-11-2021
            write(616,*) "reason to end:",whichend, lch, wrongval                                                  !!>> HC 11-11-2021 why we killed this indv
         else                                                                                                      !!>> HC 11-11-2021
            write(616,*) "reason to end:",whichend                                                                 !!>> HC 11-11-2021
         endif                                                                                                     !!>> HC 11-11-2021
     close(616)                                                                                                    !!>> HC 15-10-2020
     cx="pwd >> "//trim(carg)//".logb"                                                                             !!>> HC 29-9-2021 Save the name of the individual (evolution)
     call system(cx)                                                                                               !!>> HC 29-9-2021
     cx="paste "//trim(carg)//".logb "//trim(carg)//".loga > "//trim(carg)//".log"                                 !!>> HC 29-9-2021 paste the name and the reason to end
     call system(cx)                                                                                               !!>> HC 29-9-2021 Remove the temporal files
     cx="rm "//trim(carg)//".logb"                                                                                 !!>> HC 29-9-2021
     call system(cx)                                                                                               !!>> HC 29-9-2021
     cx="rm "//trim(carg)//".loga"                                                                                 !!>> HC 29-9-2021
     call system(cx)                                                                                               !!>> HC 29-9-2021
     if (ffufi(whichend,2)==1)then                                                                                 !!>> HC 29-8-2021  If the end is lethal 
        cx='echo "0.00" > '//trim(carg)//'fitness'                                                                 !!>> HC 17-11-2020 We kill the embryo 
        call system(cx)                                                                                            !!>> HC 17-11-2020 and it has no fitness 
     endif                                                                                                         !!>> HC 17-11-2020
                                                             

     call writesnap                                                                                                !!>> HC 18-9-2020
  
     stop                                                                                                          !!>> HC 18-9-2020
45 print *,"end of file error"                                                                                     !!>> HC 18-9-2020
  stop                                                                                                             !!>> HC 18-9-2020
46 print *,"error in writing",trim(carg)//"t"                                                                      !!>> HC 18-9-2020
  stop                                                                                                             !!>> HC 18-9-2020
48 print *,"other end"                                                                                             !!>> HC 18-9-2020
  stop                                                                                                             !!>> HC 18-9-2020


end subroutine filter 




!*******************************************************************

subroutine cellbreak   !>>> 17-1-14
  real*8  ax,ay,az,ida
  integer doapo(nd)
  integer nocon(nd,nd)
  integer i,j,k,ii,jj,kk,iii,jjj,kkk,ik

  ! now we check that the cell is not split in two or more parts
  doapo=0
  nocon=0

  do i=1,ncels

    kkk=0
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%marge==0) then ; kkk=ii ; ida=node(kkk)%add ; exit ; end if
    end do

    if (kkk==0) then !it means that the cell has no nucleus and then IT MUST DIE!!!!   !>>> Is 11-6-14
      do j=1,cels(i)%nunodes  !>>> Is 11-6-14
        ii=cels(i)%node(j)    !>>> Is 29-6-14
!print *,ii,i,j,cels(i)%node(:cels(i)%nunodes),ii,"ii",nd
        doapo(ii)=1           !>>> Is 11-6-14
      end do                  !>>> Is 11-6-14
      kkk=1  ! >>> Is 11-6-14
    else
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
        if (sqrt((ax-node(kkk)%x)**2+(ay-node(kkk)%y)**2+(az-node(ii)%z)**2)<node(ii)%add+ida) then
          doapo(ii)=1
        end if
      end do
    end if   ! >>> Is 11-6-14

    doapo(kkk)=1
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
      ida=node(ii)%add
      do jj=1,cels(i)%nunodes
        if (j==jj) cycle
        iii=cels(i)%node(jj)
        if (sqrt((ax-node(iii)%x)**2+(ay-node(iii)%y)**2+(az-node(iii)%z)**2)<node(iii)%add+ida) then
          nocon(ii,iii)=1
          nocon(iii,ii)=1
          if (doapo(ii)==1) then
            doapo(iii)=1
          else
            if (doapo(iii)==1) doapo(ii)=1
          end if
        end if
      end do
    end do

    do k=1,cels(i)%nunodes/2
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        do jj=1,cels(i)%nunodes
          if (j==jj) cycle
          iii=cels(i)%node(jj)
          if (nocon(ii,iii)==1) then
            if (doapo(ii)==1) then
              doapo(iii)=1
            else
              if (doapo(iii)==1) doapo(ii)=1
            end if
          end if
        end do
      end do
    end do
  end do

  ik=1
  do while(ik<=nd)
    if (doapo(ik)==0) then
      if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
        call apoptosis(ik)
      else
        if (doapo(node(ik)%altre)==0) then
          call apoptosis(ik)
        end if
      end if
      ik=ik+1
    else
      ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
    end if
  end do


end subroutine

!********************************************************************

subroutine polarization
integer:: celd,nnod,tipi,ggr,ccen
real*8::a,b,c,d,e,ax,ay,az,bx,by,bz,cx,cy,cz,ix,iy,iz,alfa,s
      do celd=1,ncels
        tipi=cels(celd)%ctipus
        nnod=cels(celd)%nunodes        
        if (nnod==0) cycle      ! >>> Is 10-5-14
        iy=1d10 ; cx=0d0 ; cy=0d0 ; cz=0d0
	a=cels(celd)%cex ; b=cels(celd)%cey ; c=cels(celd)%cez   

        do i=1,nnod                                                     ! (gen) in the centroid (in the closest node)
          j=cels(celd)%node(i)
          if(node(j)%tipus==1.or.node(j)%tipus==3)then !in epithelial cells, polarity is planar, so we only take one layer of nodes >>>Miquel22-10-13
            d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
          end if
          if(d.le.iy)then;iy=d;ccen=j;endif             
        end do   

        alfa=0.0d0                                                 ! concentration in the central node
        do k=1,npag(nparam_per_node+8)    
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ccen,kk)>0.0d0) then
            alfa=alfa+gex(ccen,kk)*gen(kk)%e(nparam_per_node+8)   ! wa in units of probability such that it makes things to go from 0 to 1
          end if
        end do  

        ix=0d0 ; iy=0d0 ; iz=0d0                                        ! vector of the gradient within a cell
        do i=1,nnod                                                     
            j=cels(celd)%node(i)
            if(node(j)%tipus==1.or.node(j)%tipus==3)then          
              d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
              if (d<epsilod) cycle
              d=1d0/d                                                   ! module of radial vectors to get unitary vectors     
              s=0.0d0
              do k=1,npag(nparam_per_node+8)
                kk=whonpag(nparam_per_node+8,k)
                if (gex(j,kk)>0.0d0) then
                  s=s+gex(j,kk)*gen(kk)%e(nparam_per_node+8)
                end if
              end do
              ix=ix+((node(j)%x-a)*d)*(s-alfa)                   ! and ignore shape/size effects
              iy=iy+((node(j)%y-b)*d)*(s-alfa)
              iz=iz+((node(j)%z-c)*d)*(s-alfa)
            end if
        end do

        if((ix.eq.0).and.(iy.eq.0).and.(iz.eq.0))then            ! if the gene has uniform expresion, the vector is random ! >>>Miguel1-7-14
          call random_number(a)                                  ! >>>Miguel1-7-14
          k=int(a*nvaloq)+1                                      ! >>>Miguel1-7-14
          cels(celd)%polx=particions_esfera(k,1)                 ! >>>Miguel1-7-14
          cels(celd)%poly=particions_esfera(k,2)                 ! >>>Miguel1-7-14
          cels(celd)%polz=particions_esfera(k,3)                 ! >>>Miguel1-7-14
        else                                                     ! >>>Miguel1-7-14
          a=ix**2+iy**2+iz**2 
          if(a==0)then
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0	! unitary resultant vector (gradient polarization)
          else
            d=1d0/sqrt(a)
            cels(celd)%polx=ix*d ; cels(celd)%poly=iy*d ; cels(celd)%polz=iz*d	! unitary resultant vector (gradient polarization)
          end if
          if((ix.eq.0d0).and.(iy.eq.0d0).and.(iz.eq.0d0))then                     ! miguel27-11-13
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0
          endif   ! miguel27-11-13
        endif                                                    ! >>>Miguel1-7-14
      end do
end subroutine

!*******************************************************************

subroutine polarizationisaac
integer i,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d,aa,bb,cc 

do i=1,ncels
  a=cels(i)%cex ; b=cels(i)%cey ; c=cels(i)%cez
  sx=0.0d0      ; sy=0.0d0      ; sz=0.0d0
  if (node(cels(i)%node(1))%tipus<3) then
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%tipus==1) then
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%e(nparam_per_node+8)*delta
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          d=s/sqrt(aa**2+bb**2+cc**2)
          sx=sx+d*aa ; sy=sy+d*bb ; sz=sz+d*cc
        end if
      end if
    end do
    aa=sqrt(sx**2+sy**2+sz**2)
    if (aa>0.0d0) then
      d=1d0/aa
      cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
    else
      cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
    end if
  else
    if (node(cels(i)%node(1))%tipus==3) then
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%e(nparam_per_node+8)
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          sx=sx+s*aa ; sy=sy+s*bb ; sz=sz+s*cc
        end if
      end do
      aa=sqrt(sx**2+sy**2+sz**2)
      if (aa>0.0d0) then  !epsilod) then
        d=1d0/aa
        cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
      else
        cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
      end if
    end if
  end if
end do
end subroutine

!*************************************************************************************

subroutine change_minsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+10)
      kk=whonpag(nparam_per_node+10,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%e(nparam_per_node+10)  !wa in units of number of nodes but it can be roughly understood as space-req units
      end if                                           ! THIS WA CAN BE NEGATIVE
    end do
  end do
  s=s*delta
!print *,delta,"ACHTUNG"
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%minsize_for_div=cels(ick)%minsize_for_div+s  !this can make a SUDDEN CHANGE
    if (cels(ick)%minsize_for_div<1) cels(ick)%minsize_for_div=1
!  end if
end do

end subroutine

!*************************************************************************************

!>>> Is 5-2-14
subroutine change_maxsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+15)
      kk=whonpag(nparam_per_node+15,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%e(nparam_per_node+15)  ! THIS WA MAY BE NEGATIVE
      end if
    end do
  end do
  s=s*delta
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%maxsize_for_div=cels(ick)%maxsize_for_div+s  !wa in units of space-req roughly this can make a SUDDEN CHANGE
    if (cels(ick)%maxsize_for_div<1) cels(ick)%maxsize_for_div=1
!  end if
end do

!>>> Is 5-2-14

end subroutine


!**************************************************************************************

subroutine emt !epithelial-mensenchymal transitions

               ! this simply makes the tipus equal to 3 but it needs to approach both sides of the originally epithelial cell so that 
               ! the two parts do not drift away since they are more far away than da (they are at reqs normally)
               ! this transition has to be sudden otherwise other existing forces may impede the nodes from getting close enough to each other
               ! that wouldnt be realistic since it would be as if the cell would explode.

integer ick,j,k,ii,kk,iii,jjj,kkk,iv
real*8 a,b,c,d,aa,bb,cc,dd,s

e: do ick=1,ncels
  if (cels(ick)%ctipus==3) cycle
  s=0.0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+13)
      kk=whonpag(nparam_per_node+13,k)
      a=gex(ii,kk) !don't change this, weird stuff will happen if you do !>>Miquel18-12-14
      s=s+a*gen(kk)%e(nparam_per_node+13)
   !if(gex(ii,kk)>epsilod) print*,"emt gene",kk,"node",ii,"cell",ick,"conc",gex(ii,kk)
    end do
  end do
  cels(ick)%temt=cels(ick)%temt+s*delta/real(cels(ick)%nunodes)
  
  if (cels(ick)%temt>=1.0d0) then       ! wa in units of probability: this is an arbitrary value but it is fine since the rest can be re-scaled 
    !so this cell goes to emt
 !print*,"EMT HAPPENS CELL",ick
    cels(ick)%ctipus=3  ! >>> Is 4-1-14
    do jjj=1,cels(ick)%nunodes
      iii=cels(ick)%node(jjj)          
      if (node(iii)%tipus==1) then
        iv=node(iii)%altre
        a=node(iv)%x-node(iii)%x ; b=node(iv)%y-node(iii)%y ; c=node(iv)%z-node(iii)%z 
        d=1d0/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
        aa=0.5d0*(node(iv)%x+node(iii)%x) ; bb=0.5d0*(node(iv)%y+node(iii)%y) ; cc=0.5d0*(node(iv)%z+node(iii)%z) 
        dd=(node(iii)%add+node(iv)%add)*0.1d0
        node(iii)%x=aa-a*dd 
        node(iii)%y=bb-b*dd 
        node(iii)%z=cc-c*dd
        node(iv)%x=aa+a*dd 
        node(iv)%y=bb+b*dd 
        node(iv)%z=cc+c*dd
        node(iii)%altre=0
        node(iv)%altre=0
        node(iii)%tipus=3
        node(iv)%tipus=3
        !if only one node per cell both nodes of the epithelial cell get a nucleus
        if (ffu(1)==0) then ; node(iii)%marge=0 ; node(iv)%marge=0 ; end if ! >>> Is 10-10-14
      end if
    end do
  end if
end do e

end subroutine

!*******************************************************************************************************

subroutine diffusion_of_reqcr
  real*8 hreqcr(nd),hreqp(nd),hreqc(nd)
  integer i,j,k,ii,jj,kk
  real*8 a,b,c,d

  do i=1,nd
    a=0.0d0 ; b=0.0 ; c=0.0d0
    do ii=1,nneigh(i)
      k=neigh(i,ii)
      if (node(i)%icel/=node(k)%icel) cycle    ! only within the same cell
      if (node(i)%tipus/=node(k)%tipus) cycle  ! only within the same side of the cell
      if (node(i)%tipus>2) cycle               ! only for epithelial cells
      d=dneigh(i,ii)
      a=a+(node(k)%grd-node(i)%grd) !/(d+1d0)
      !b=b+(node(k)%cod-node(i)%cod) !/(d+1d0)
      !c=c+(node(k)%pld-node(i)%pld) !/(d+1d0)
    end do
    hreqcr(i)=a*dif_req ! >>> 11-6-14
    !hreqc(i)=b*dif_req  ! >>> 11-6-14
    !hreqp(i)=c*dif_req  ! >>> 11-6-14
  end do  
  do i=1,nd
    node(i)%grd=node(i)%grd+delta*hreqcr(i)
    !node(i)%cod=node(i)%cod+delta*hreqc(i)
    !node(i)%pld=node(i)%pld+delta*hreqp(i)
  end do
end subroutine

end module nexus
