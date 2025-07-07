!fortran general.mod.f90 neighboring.mod.f90 conservative_R.mod.f90  geompack3.f90 io.mod.f90 genetic.mod.f90 prueba.f90 -o p.e
program EMD_comp

    use general
    use genetic
    use neighboring
    use io
    use conservative_R ! module created by Renske, contains simple EMD and PF H measures of similarity to ancestor

implicit none

    character*140 :: m1, m2
    integer        :: SV
    call getarg(1,m1)
    call getarg(2,m2)

    call shared_volume_two_morphs(SV,m1,m2)

end program EMD_comp
