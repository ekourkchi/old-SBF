#!/bin/bash

export command=$1
export pwd=$PWD


if [ !  -d "./bin" ] 
then
    mkdir $pwd"/bin"
fi

if [ !  -d "./lib" ] 
then
    mkdir $pwd"/lib"
fi



function setup_xpa() {
    cd xpa
    bash configure
    make 
    make all
    cp xpaget ../bin/.
    cp xpaset ../bin/.
    cp xpainfo ../bin/.
    cp xpaaccess ../bin/.
    cp xpans ../bin/.
    cp xpamb ../bin/.
    make clean
    cd $pwd
}



function setup_likenew() {
    cd sbfsrc
    make all
    mv ./likenew6.so ../lib/.
    make clean
    cd $pwd
}


function setup_likenew() {
    cd sbfsrc
    make all
    mv ./likenew6.so ../lib/.
    make clean
    cd $pwd
}


function setup_dophot() {
    cd dophot_jag
    make clean 
    make all
    cp dophot ../bin/.
    cd $pwd
}


case $command in 
 xpa)
    setup_xpa
;;
 likenew) 
    setup_likenew
;;
dophot)
    setup_dophot
;;
*)
    setup_likenew
    setup_dophot
    setup_xpa
;;
esac




