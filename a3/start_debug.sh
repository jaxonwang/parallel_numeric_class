#!/bin/bash

pid=$(pgrep $1|head -n 1);gdb -q \
-ex "attach ${pid}" \
-ex "set variable ifl=1" \
-ex "finish"
