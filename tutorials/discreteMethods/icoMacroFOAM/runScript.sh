NEVRY=$(awk '/^noOfMicroSteps/' system/immDict | awk '{print substr($2, 1, length($2)-1)}')
NITER=$(awk '/^nIter/' system/immDict | awk '{print substr($2, 1, length($2)-1)}')
NRUN=$(awk -v nevery=$NEVRY -v nIter=$NITER 'BEGIN{printf("%d",nevery*nIter)}')
echo $NEVRY $NITER $NRUN
STRAVG=100000
VAR_LIST=$(echo "-var NEVRYAVG $NEVRY -var NRUN $NRUN -var STRAVG $STRAVG")
echo $VAR_LIST

#edit this line to include the correct path to MUI-enabled LAMMPS solver
LMPS_PATH=~/lammps/src/lmp_openmpi_c++11

mpirun -np 1 icoMacroFOAM : -np 4 ${LMPS_PATH} ${VAR_LIST} -in lammps/micro_1/lammps.in : -np 4 ${LMPS_PATH} ${VAR_LIST} -in lammps/micro_2/lammps.in : -np 4 ${LMPS_PATH} ${VAR_LIST} -in lammps/micro_3/lammps.in : -np 4 ${LMPS_PATH} ${VAR_LIST} -in lammps/micro_4/lammps.in : -np 4 ${LMPS_PATH} ${VAR_LIST} -in lammps/micro_5/lammps.in
