#------pjsub option --------#
#PJM -L rscgrp=lecture
#PJM -L node=1
#PJM --mpi proc=56
#PJM -L elapse=00:15:00
#PJM -g gt54
#PJM -j
#-------Program execution -------#
#mpiexec.hydra -n ${PJM_MPI_PROC} ./ldlt_rightlooking 12800 32


args=(
8000	2
10080	4
11552	6
12704	8
13664	10
14528	12
16000	16
17248	20
18304	24
19296	28
20160	32
20960	36
21728	40
22432	44
23072	48
23712	52
24288	56
)

np=()
size=()
for i in ${!args[@]}; do
if [ `echo "$i % 2"| bc` -eq 0 ];then
	size+=(${args[$i]})
else
	np+=(${args[$i]})
fi
done

for i in ${!np[@]}; do
tmp_np=${np[$i]}
tmp_size=${size[$i]}
echo "mpiexec.hydra -n ${tmp_np} ./ldlt_leftlooking $tmp_size 32"
mpiexec.hydra -n ${tmp_np} ./ldlt_leftlooking $tmp_size 32
done
