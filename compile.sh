#!/bin/sh
# set to executable using: chmod 755 example.sh

# compile source
gcc -g -Wall -c src/*.c
# copy headers
cp src/*.h include/
# move binaries
mv *.o bin/