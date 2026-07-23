#!/bin/bash

mkdir -p build
mkdir -p output

clang -Wall -lm -O2 -std=gnu17 main.c -o ./build/main 
