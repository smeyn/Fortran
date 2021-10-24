REM  compile all files
REM gfortran -c fpback.f fpbspl.f fpdisc.f fpgivs.f fporde.f fprank.f fprati.f fprota.f fpsurf.f surfit.f
cd fitpack\Surfit\test\fortran
dir *.f
gfortran -g -c ..\..\fortran\*.f
gfortran -g -c  ..\..\..\BISPEV\*.f
gfortran -g mnsurf.f *.o -o mnsurf.exe -fdump-tree-original