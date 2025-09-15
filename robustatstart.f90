program assymetryatstart

    use general
    use genetic
    use neighboring
    use io
    use conservative_R ! module created by Renske, contains simple EMD and PF H measures of similarity to ancestor
  
  
    implicit none
    real*8 :: fithc, ifithc, mefithc, emd                                                       !!>> HC 3-7-2023
    real*8 :: osfithc, otfithc, offithc, stfithc, sffithc, tffithc, symch                       !!>> HC 26-11-2020
    character*140 :: input                                         !!>> HC 26-11-2020    
    integer ::  ihc, jhc, khc, nnodes, naltech, surface, volume, outside, total, inside, differ, centrow
    integer ::  iihc,jjhc,kkhc, ord1, ord2, ord3, lhc, prev, changes, mhc, nhc, neich, neichi, newdots, vhc, phc
    real*8, dimension(1:3) ::  u
    real*8 :: upv, sumd,  modu, ahc, bhc, chc, ahc2, bhc2, chc2, ahc3, bhc3, chc3, u1, u2, pershared
    real*8 :: maxx, minx, maxy, miny, maxz, minz, maxadd, minadd, anchx, anchy, anchz, anchadd, totmax, totmin, totanch
    real*8, allocatable, dimension(:,:) :: ncoords
    integer, allocatable, dimension(:,:,:) :: bfillz
    fithc=0.0d0; ifithc=0.0d0 ; mefithc=0.0d0                                                   !!>> HC 26-11-2020
    osfithc=0.0d0; otfithc=0.0d0; offithc=0.0d0; stfithc=0.0d0; sffithc=0.0d0; tffithc=0.0d0    !!>> HC 26-11-2020
    distfitscale=2.250d0; distfitmag=1.0d0

    
    
    call getarg(1,input)
    call iniread
    call readsnap(input)
    call iniboxes
    call neighbor_build
    
  maxx=0.0d0; minx=10000; maxy=0.0d0; miny=10000; maxz=0.0d0; minz=10000                    !!>> HC 4-3-2024
  anchx=0.0d0; anchy=0.0d0; anchz=0.0d0; totmax=0.0d0; totmin=0.0d0; totanch=0.0d0          !!>> HC 4-3-2024
  maxadd=0.0d0                                                                              !!>> HC 4-3-2024
  
  
  newdots=0                                                                                 !!>> HC 4-3-2024 NUMBER OF NEW POINTS
  do ihc=1,nd                                                                               !!>> HC 4-3-2024 placed randomly to fill in gaps in the epithelium
     if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
     newdots=newdots+1                                                                      !!>> HC 4-3-2024
     do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
        neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
        if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
        do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
           neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
           if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
           if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
              newdots=newdots+10                                                            !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  if(allocated(ncoords))deallocate(ncoords)                                                 !!>> HC 4-3-2024 This stores the coordinates of the morphology
  allocate(ncoords(1:newdots,1:3))                                                          !!>> HC 4-3-2024  and the newly added random points
  ncoords=0.0d0                                                                             !!>> HC 4-3-2024
  
  ord1=0                                                                                    !!>> HC 4-3-2024 STORE THE COORDINATES OF NODES
  do ihc=1,nd                                                                               !!>> HC 4-3-2024
     if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
     ord1=ord1+1                                                                            !!>> HC 4-3-2024
     ncoords(ord1,1)=node(ihc)%x                                                            !!>> HC 4-3-2024
     ncoords(ord1,2)=node(ihc)%y                                                            !!>> HC 4-3-2024
     ncoords(ord1,3)=node(ihc)%z                                                            !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  do ihc=1,nd                                                                               !!>> HC 4-3-2024 ADD THE RANDOM POINTS
     if(node(ihc)%tipus.ne.1)cycle                                                          !!>> HC 4-3-2024
     do jhc=1,nneigh(ihc)                                                                   !!>> HC 4-3-2024
        neich=neigh(ihc,jhc)                                                                !!>> HC 4-3-2024
        if(node(neich)%tipus.ne.1)cycle                                                     !!>> HC 4-3-2024
        do khc=1,nneigh(neich)                                                              !!>> HC 4-3-2024
           neichi=neigh(neich,khc)                                                          !!>> HC 4-3-2024
           if(node(neichi)%tipus.ne.1)cycle                                                 !!>> HC 4-3-2024
           if(neichi==ihc)cycle                                                             !!>> HC 4-3-2024
              do phc=1,10                                                                  !!>> HC 4-3-2024
                  ahc = node(neich)%x - node(ihc)%x                                         !!>> HC 4-3-2024
                  bhc = node(neich)%y - node(ihc)%y                                         !!>> HC 4-3-2024
                  chc = node(neich)%z - node(ihc)%z                                         !!>> HC 4-3-2024
                  
                  ahc2 = node(neichi)%x - node(ihc)%x                                       !!>> HC 4-3-2024
                  bhc2 = node(neichi)%y - node(ihc)%y                                       !!>> HC 4-3-2024
                  chc2 = node(neichi)%z - node(ihc)%z                                       !!>> HC 4-3-2024
                  
                  call random_number(u1); call random_number(u2);                           !!>> HC 4-3-2024
                  if ( (u1+u2) > 1.0d0)then                                                 !!>> HC 4-3-2024
                     u1=1-u1                                                                !!>> HC 4-3-2024
                     u2=1-u2                                                                !!>> HC 4-3-2024
                  endif                                                                     !!>> HC 4-3-2024
                  
                  ahc3 = u1*ahc + u2*ahc2                                                   !!>> HC 4-3-2024
                  bhc3 = u1*bhc + u2*bhc2                                                   !!>> HC 4-3-2024
                  chc3 = u1*chc + u2*chc2                                                   !!>> HC 4-3-2024
                  ord1=ord1+1                                                               !!>> HC 4-3-2024
                  ncoords(ord1,1)= node(ihc)%x + ahc3                                       !!>> HC 4-3-2024
                  ncoords(ord1,2)= node(ihc)%y + bhc3                                       !!>> HC 4-3-2024
                  ncoords(ord1,3)= node(ihc)%z + chc3                                       !!>> HC 4-3-2024
              enddo                                                                         !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  do ihc=1,nd                                                                               !!>> HC 4-3-2024 CALCULATE THE BIGGEST COORDINATE IN THE MORPH
     if(abs(node(ihc)%x)>maxx) maxx=abs(node(ihc)%x)                                        !!>> HC 4-3-2024
     if(abs(node(ihc)%y)>maxy) maxy=abs(node(ihc)%y)                                        !!>> HC 4-3-2024
     if(abs(node(ihc)%z)>maxz) maxz=abs(node(ihc)%z)                                        !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  totmax=maxx                                                                               !!>> HC 4-3-2024
  if (maxy>totmax) totmax=maxy                                                              !!>> HC 4-3-2024
  if (maxz>totmax) totmax=maxz                                                              !!>> HC 4-3-2024 The biggest coordinate
  
  maxadd=nodeo(1)%add*2                                                                     !!>> HC 4-3-2024
  urv=1/maxadd                                                                              !!>> HC 4-3-2024
  nboxes=(urv*totmax) +2                                                                    !!>> HC 4-3-2024  THE NUMBER OF BOXES
  
  if(allocated(bfillz))deallocate(bfillz)                                                   !!>> HC 4-3-2024 THIS ARRAY STORES THE BOXES
  allocate(bfillz(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))                            !!>> HC 4-3-2024
  bfillz=0                                                                                  !!>> HC 4-3-2024
  
  do ihc=1,ord1                                                                             !!>> HC 4-3-2024 INTRODUCE THE COORDINATES IN THEIR CORRESPONDING BOXES
     jhc=nint(ncoords(ihc,1)*urv)                                                           !!>> HC 4-3-2024 
     khc=nint(ncoords(ihc,2)*urv)                                                           !!>> HC 4-3-2024
     nhc=nint(ncoords(ihc,3)*urv)                                                           !!>> HC 4-3-2024
     bfillz(jhc,khc,nhc)=1                                                                  !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  bfillz(nboxes,nboxes,nboxes)=2                                                            !!>> HC 4-3-2024 THE FILLING ALGORITHM
  changes=11; ord1=0                                                                        !!>> HC 4-3-2024 Fills with either 2 or 3 the boxes OUTSIDE the morphology
  do while(changes>0)                                                                       !!>> HC 4-3-2024 2-3= OUT
     changes=0                                                                              !!>> HC 4-3-2024 0= IN
     do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024 1= SURFACE
        do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
           do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
              if( bfillz(ihc,jhc,khc).ne.2)cycle                                            !!>> HC 4-3-2024
              do nhc=-1,1                                                                   !!>> HC 4-3-2024
                 ord1=ihc+nhc                                                               !!>> HC 4-3-2024
                 if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                 if(bfillz(ord1,jhc,khc).ne.0)cycle                                         !!>> HC 4-3-2024
                 bfillz(ord1,jhc,khc)=2                                                     !!>> HC 4-3-2024
                 changes=changes+1                                                          !!>> HC 4-3-2024
              enddo                                                                         !!>> HC 4-3-2024
              do nhc=-1,1                                                                   !!>> HC 4-3-2024
                 ord1=jhc+nhc                                                               !!>> HC 4-3-2024
                 if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                 if(bfillz(ihc,ord1,khc).ne.0)cycle                                         !!>> HC 4-3-2024
                 bfillz(ihc,ord1,khc)=2                                                     !!>> HC 4-3-2024
                 changes=changes+1                                                          !!>> HC 4-3-2024
              enddo                                                                         !!>> HC 4-3-2024
              do nhc=-1,1                                                                   !!>> HC 4-3-2024
                 ord1=khc+nhc                                                               !!>> HC 4-3-2024
                 if(ord1>nboxes.or.ord1<(-nboxes))cycle                                     !!>> HC 4-3-2024
                 if(bfillz(ihc,jhc,ord1).ne.0)cycle                                         !!>> HC 4-3-2024
                 bfillz(ihc,jhc,ord1)=2                                                     !!>> HC 4-3-2024
                 changes=changes+1                                                          !!>> HC 4-3-2024
              enddo                                                                         !!>> HC 4-3-2024
              bfillz(ihc,jhc,khc)=3                                                         !!>> HC 4-3-2024
           enddo                                                                            !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  
  surface=0; volume=0; outside=0; total=0                                                   !!>> HC 4-3-2024  Counting outside and inside boxes
  do ihc=-nboxes,nboxes                                                                     !!>> HC 4-3-2024
     do jhc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
        do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
           if (bfillz(ihc,jhc,khc)==1)then                                                  !!>> HC 4-3-2024
              surface=surface+1                                                             !!>> HC 4-3-2024
           elseif(bfillz(ihc,jhc,khc)==0)then                                               !!>> HC 4-3-2024
              volume=volume+1                                                               !!>> HC 4-3-2024
           else                                                                             !!>> HC 4-3-2024
              outside=outside+1                                                             !!>> HC 4-3-2024
           endif                                                                            !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  enddo                                                                                     !!>> HC 4-3-2024
  
  if (volume>0)then                                                                         !!>> HC 4-3-2024 If the morphology is not broken
     total=surface+volume                                                                   !!>> HC 4-3-2024 CALCULATE THE ASYMMETRY BY SHARED VOLUME
     do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
        do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
           do khc=-nboxes,nboxes                                                            !!>> HC 4-3-2024 Transform all the boxes into 1=IN 0=out
              if ( bfillz(ihc,jhc,khc).le.1)then                                            !!>> HC 4-3-2024
                 bfillz(ihc,jhc,khc)=1                                                      !!>> HC 4-3-2024
              else                                                                          !!>> HC 4-3-2024
                 bfillz(ihc,jhc,khc)=0                                                      !!>> HC 4-3-2024
              endif                                                                         !!>> HC 4-3-2024
           enddo                                                                            !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  
  
      centrow=0                                                                             !!>> HC 4-3-2024 NON-SHARED VOLUME
     do ihc=-nboxes,nboxes                                                                  !!>> HC 4-3-2024
        do jhc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
        if ( bfillz(ihc,0,jhc)==0)cycle                                                     !!>> HC 4-3-2024
        centrow=centrow+1                                                                   !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
  
     differ=0                                                                               !!>> HC 4-3-2024
     do ihc=1,nboxes                                                                        !!>> HC 4-3-2024
        jhc=-ihc                                                                            !!>> HC 4-3-2024
        do khc=-nboxes,nboxes                                                               !!>> HC 4-3-2024
           do nhc=-nboxes,nboxes                                                            !!>> HC 4-3-2024
              if( bfillz(ihc,khc,nhc) == bfillz(jhc,khc,nhc))cycle                          !!>> HC 4-3-2024
              differ=differ+1                                                               !!>> HC 4-3-2024 !!»»AL 18-10-24: each half is with respect to x axis so righ side is +x and left -x 
           enddo                                                                            !!>> HC 4-3-2024
        enddo                                                                               !!>> HC 4-3-2024
     enddo                                                                                  !!>> HC 4-3-2024
     symch=real(differ)/real(total-centrow)                                                 !!>> HC 4-3-2024
  else                                                                                      !!>> HC 4-3-2024 if the morphology is broken we will discard it
     symch=666.0d0                                                                          !!>> HC 4-3-2024
  endif                                                                                     !!>> HC 4-3-2024
    
  
    open(777,file="assymetryatstart.dat")
    write(777,*) symch                     !!>> HC 26-11-12020
    close(777)
  
  end program assymetryatstart
  