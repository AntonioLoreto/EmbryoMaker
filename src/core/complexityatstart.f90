program complexity

    use io
    use general
    use genetic
    use neighboring

    integer                 :: ihc,jhc,khc,surface,volume,outside,total,differ,shared,ord1,changes,nhc,ijk,puntos_fill,coco,iter,cb
    integer                 :: neich,neichi,newdots,cuantos,w_sect,compa_tot,secciones,phc,comparaciones,size_box,numrowspersect,ca
    integer                 :: guardar
    real*8                  :: ahc,bhc,chc,ahc2,bhc2,chc2,ahc3,bhc3,chc3,maxx,minx,maxy,miny,maxz,maxadd,totmax,crompio,conta_n_nodos
    real*8                  :: u1,u2,pershared,pershared_z,pershared_x,pershared_y,per,alfa,fitelli,voll,cociente,crompiopokito
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

    call getarg(1,input)
    call iniread
    call readsnap(input)
    call iniboxes
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
       volumenes(1) = 1
       volumenes(2) = 1
       volumenes(3) = 1
    end if
 
    open (123, file="fitnessatstart.dat")
          write(123,*) fitelli
    close(123)
    
end program complexity