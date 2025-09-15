!gfortran -O2 -Wall -fbounds-check  general.mod.f90  io.mod.f90 genetic.mod.f90 recombination.mod.f90 reco-test.f90 -o rec.e
program recombination1

use io
use general
use genetic
use recombination

implicit none

character*700 :: rechosen,string1
call getarg(1,rechosen)

string1 = trim(rechosen)// CHAR(0)

call do_recombination_4_evopt1(string1) !the last is treated as "donor"

end program recombination1