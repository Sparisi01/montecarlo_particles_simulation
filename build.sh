#!/bin/bash

mkdir -p build
mkdir -p output

clang -Wall -lm -O2 -std=gnu17 main.c -o ./build/main 

if [ $? -eq 0 ]; then
    echo "Compilation done succesfully"
    echo "Execution"
    ./build/main
else
    echo "Compilation Error"
fi
