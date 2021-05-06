cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
SET(CMAKE_CXX_FLAGS   "-w -g -O3")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_STANDARD 17)

#Elasticity library
project(Elas)
file(                  GLOB_RECURSE source_files sources/*)
file(                  GLOB_RECURSE header_files sources/*.h)
ADD_LIBRARY(           Elas SHARED ${source_files})
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")
INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
LINK_DIRECTORIES(     "$ENV{HOME}/lib")
target_link_libraries( Elas ${Commons})
INSTALL(
  TARGETS Elas
  LIBRARY  DESTINATION "$ENV{HOME}/lib"
  RUNTIME  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")

#Executable
project(main)
add_executable(        elastic sources/elastic.c)
target_link_libraries( elastic ${Commons} Elas)
INSTALL(               TARGETS elastic RUNTIME DESTINATION "$ENV{HOME}/bin")


 



