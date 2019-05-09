#!/bin/bash
# creates C library

libname=bjutils

# move .c into source/
# move .h into include/

cd source/
LB=" -lgsl -lgslcblas"

# create static library
gcc -c *.c
ar rs lib${libname}.a *.o

# create shared library
gcc -c -fpic *.c
gcc -shared ${LB} -o lib${libname}.so *.o

mkdir ../lib/
mv lib${libname}.a ../lib/
mv lib${libname}.so ../lib/

