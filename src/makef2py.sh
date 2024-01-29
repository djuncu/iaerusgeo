# this export statement is the reason for the separate shell script....
export LDFLAGS="-Wl,-rpath,src/"

F2PYFLAGS="-L. -lAerusL3"
S="C_F90"

f2py -c -m AerusL3 $S/py_ifc.f90 $S/start_mod.f90 $S/start_exe.f90 $S/smoothing.f90 --f90flags="-fallow-invalid-boz -fallow-argument-mismatch -ffree-line-length-none -m$1 -fcheck=all -fmax-array-constructor=137000" $F2PYFLAGS
#f2py -c -m AerusL3 $S/py_ifc.f90 $S/start_mod.f90 $S/start_exe.f90 $S/smoothing.f90 --f90flags="-fallow-invalid-boz -fallow-argument-mismatch -ffree-line-length-none -m$1 -fcheck=all -fmax-array-constructor=137000" $F2PYFLAGS
