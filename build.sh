#! /bin/bash

mkdir -p ./build
cd ./build
cmakefile_timestamp=`stat -c %Y ../CMakeLists.txt`
if [ ! -f "./timestamp" ]; then
    touch ./timestamp
    # echo "touch and cmake"
    cmake ..
else
    cmakefile_lasttimestamp=`stat -c %Y ./timestamp`
    if [ $cmakefile_timestamp -gt $cmakefile_lasttimestamp ]; then
        touch ./timestamp
        # echo "cmake"
        cmake ..
    fi
fi

make -j10