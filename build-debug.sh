rm -rf build-debug
mkdir build-debug
cd build-debug
cmake ..  -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON -DCMAKE_TOOLCHAIN_FILE=/home/orebas/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake  --build .
