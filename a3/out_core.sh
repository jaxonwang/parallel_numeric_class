#------pjsub option --------#
#PJM -L rscgrp=lecture
#PJM -L node=2
#PJM --mpi proc=112
#PJM -L elapse=00:15:00
#PJM -g gt54
#PJM -j
#-------Program execution -------#
echo "mpiexec.hydra -n ${PJM_MPI_PROC} ./ldlt_rightlooking 24288 32"
mpiexec.hydra -n ${PJM_MPI_PROC} ./ldlt_rightlooking 24288 32
