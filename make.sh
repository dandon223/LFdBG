#!/bin/bash
clear
rm -rf build > /dev/null
mkdir -p build > /dev/null
cd build

cmake ../
make

