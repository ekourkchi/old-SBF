#!/bin/bash

export pwd=$PWD

cd sbfsrc

make all

mv ./likenew6.so ../lib/.

make clean

cd $pwd


