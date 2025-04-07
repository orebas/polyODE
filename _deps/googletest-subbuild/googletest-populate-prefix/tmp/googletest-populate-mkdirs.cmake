# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/orebas/cpp/polyODE/_deps/googletest-src"
  "/home/orebas/cpp/polyODE/_deps/googletest-build"
  "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix"
  "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/orebas/cpp/polyODE/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
