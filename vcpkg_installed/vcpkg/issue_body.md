Package: msolve:x64-linux@0.7.5

**Host Environment**

- Host: x64-linux
- Compiler: GNU 14.2.0
-    vcpkg-tool version: 2024-10-18-e392d7347fe72dff56e7857f7571c22301237ae6
    vcpkg-scripts version: d221c5d2c 2024-10-31 (7 months ago)

**To Reproduce**

`vcpkg install `

**Failure logs**

```
CMake Error at scripts/cmake/vcpkg_download_distfile.cmake:152 (message):
  Invalid SHA512:
  1fc8c1cc61d8d67ce0e0b881422b6c2a62b4a9b309105d3e4eaedb9974e7a66699eacc4baf1cdbcc5262a60728b0d42868ef4d30ae3b6ec6a91.


      If you do not know the file's SHA512, set this to "0".
Call Stack (most recent call first):
  scripts/cmake/vcpkg_from_github.cmake:106 (vcpkg_download_distfile)
  /home/orebas/cpp/polyODE/ports/msolve/portfile.cmake:1 (vcpkg_from_github)
  scripts/ports.cmake:192 (include)



```

**Additional context**

<details><summary>vcpkg.json</summary>

```
{
  "name": "polyode",
  "version-string": "0.1.0",
  "dependencies": [
    "boost-odeint",
    "ceres",
    "gtest",
    "msolve"
  ]
}

```
</details>
