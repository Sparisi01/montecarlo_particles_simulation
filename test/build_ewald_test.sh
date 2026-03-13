#!/bin/bash

mkdir -p build
mkdir -p output

clang -Wall -lm -O2 -std=gnu17 ewald_test.c -o ./build/ewald_test

if [ $? -eq 0 ]; then
    echo "Compilation done succesfully"
    echo "Execution"
    ./build/ewald_test
else
    echo "Compilation Error"
fi
