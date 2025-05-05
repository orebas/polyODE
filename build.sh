rm -rf build
mkdir build
cd build
cmake ..  -DCMAKE_BUILD_TYPE=Release -DENABLE_SANITIZERS=OFF -DCMAKE_TOOLCHAIN_FILE=/home/orebas/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake  --build .
