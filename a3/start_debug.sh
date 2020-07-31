#!/bin/bash

pid=$(pgrep ldlt_mpi|head -n 1);gdb -q \
-ex "attach ${pid}" \
-ex "set variable ifl=1" \
-ex "finish"
