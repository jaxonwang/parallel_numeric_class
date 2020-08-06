#------pjsub option --------#
#PJM -L rscgrp=lecture
#PJM -L node=1
#PJM --mpi proc=56
#PJM -L elapse=00:15:00
#PJM -g gt54
#PJM -j
#-------Program execution -------#
#mpiexec.hydra -n ${PJM_MPI_PROC} ./ldlt_rightlooking 12800 32

#echo 1
#mpiexec.hydra -n 1 ./ldlt_rightlooking 16000 32
#echo 2
#mpiexec.hydra -n 2 ./ldlt_rightlooking 16000 32
END=56
for ((i=52;i<=END;i+=4)); do
    echo $i
    mpiexec.hydra -n ${i} ./ldlt_rightlooking 16000 32
done
