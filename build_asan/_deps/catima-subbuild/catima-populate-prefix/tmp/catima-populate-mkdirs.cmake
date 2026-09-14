# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-src"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-build"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/tmp"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/src/catima-populate-stamp"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/src"
  "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/src/catima-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/src/catima-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/ribll2026/ribll2026_www/github_code/brill2/build_asan/_deps/catima-subbuild/catima-populate-prefix/src/catima-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
