cmake_minimum_required(VERSION 2.8)

project(Elas)

set(CMAKE_BUILD_TYPE "Release")


#Options C++ (-Wno-deprecated: remove glut deprecated)

#set (CMAKE_C_FLAGS "-std=-03")

#Include files
#include_directories("include/")
include_directories("sources/" "common/")

#Compilation
#Recuperation de la liste des fichiers .c et .h
file(GLOB_RECURSE source_files sources/* common/*)
add_library(Elas SHARED ${source_files})

project(main)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/bin")
add_executable(sulastic elastic.c)
target_link_libraries(sulastic Elas)
