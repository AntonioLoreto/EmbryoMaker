#!/bin/bash


cd src/core


echo 'descending'
pwd

aleas=aleas.mod.f90

#making object files
echo 'making object files'

gfortran -S -w -fallow-invalid-boz  OpenGL_gl.f90 OpenGL_glu.f90 OpenGL_glut.f90
gfortran -w -c $aleas
gfortran -w -c general.mod.f90
gfortran -w -c geompack3.f90 gnuplotter.mod.f90 
gfortran -w -c shell.mod.f90 neighboring.mod.f90 
gfortran -w -c genetic.mod.f90 
gfortran -w -c energy.mod.f90 io.mod.f90  
gfortran -w -c polarization.mod.f90
gfortran -w -c biomechanic_pola.mod.f90 
gfortran -w -c biomechanic.mod.f90 death.mod.f90 pola.mod.f90 ecm.mod.f90 growth.mod.f90 
gfortran -w -c ic.mod.f90 mitosis.mod.f90 
gfortran -w -c single_node.mod.f90 
gfortran -w -c conservative_R.mod.f90
gfortran -w -c nexus.mod.f03 
gfortran -w -c inicial.mod.f90 
gfortran -w -c editor.mod.f90 model.mod.f90 
gfortran -w -c drawer.mod.f90 automaticon.mod.f90 
gfortran -w -c elli.f90 


#linking
echo 'linking'

gfortran -w -O2 -fbounds-check -fexceptions -fno-underscoring -fcheck=all gnuplotter.mod.o polarization.mod.o biomechanic_pola.mod.o aleas.mod.o general.mod.o neighboring.mod.f90 genetic.mod.o energy.mod.o shell.mod.o io.mod.o pola.mod.o mitosis.mod.o growth.mod.o death.mod.o single_node.mod.o ic.mod.o ecm.mod.o nexus.mod.o biomechanic.mod.o conservative_R.mod.o model.mod.o inicial.mod.o editor.mod.o drawer.mod.f90 automaticon.mod.f90 geompack3.f90 elli.o -o EMaker -lGL -lGLU -lglut

#cleaning
echo 'cleaning'
rm *.s *.o *.mod

cd ../..
mv src/core/EMaker bin

echo 'executables installed in bin/'


