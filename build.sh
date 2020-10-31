cd build
export CXX=/act/gcc-4.7.2/bin/g++
export CC=/act/gcc-4.7.2/bin/gcc
cmake ..
make
cd ..
cd script
sbatch submit.sh
