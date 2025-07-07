!gfortran -O2 -Wall -fbounds-check  general.mod.f90  io.mod.f90 genetic.mod.f90 recombination.mod.f90 reco-test.f90 -o rec.e
program recombination2

use io
use general
use genetic
use recombination

implicit none

character*700 :: parent,string2
call getarg(1,parent)

string2 = trim(parent)// CHAR(0)

call do_recombination_4_evopt2(string2) !the last is treated as "donor"

end program recombination2