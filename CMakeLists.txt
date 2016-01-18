cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE "Release")

#Variables
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "$ENV{HOME}/lib")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/bin")
SET(ELASTIC_SOLVER_LIB             Elas)
SET(COMMONS_LIB                    Commons)
LINK_DIRECTORIES(                  "$ENV{HOME}/lib")

#Suffixes
IF(APPLE)
set_target_properties(${ELASTIC_SOLVER_LIB}  ${COMMONS_LIB} PROPERTIES SUFFIX .dylib)
ELSEIF("${CMAKE_SYSTEM}" MATCHES "Linux")
set_target_properties(${ELASTIC_SOLVER_LIB}  ${COMMONS_LIB} PROPERTIES SUFFIX .so)
ENDIF()

#Commons library
include(ExternalProject)
ExternalProject_Add(
  commons
  PREFIX "./commons/"
  GIT_REPOSITORY https://github.com/ICStoolbox/Commons.git
  BUILD_COMMAND "make"
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  LOG_INSTALL ON
)

#Elasticity library
project(${ELASTIC_SOLVER_LIB})
file(GLOB_RECURSE source_files sources/*)
file(GLOB_RECURSE header_files sources/*.h)
ADD_LIBRARY(          ${ELASTIC_SOLVER_LIB} SHARED ${source_files})
ADD_DEPENDENCIES(     ${ELASTIC_SOLVER_LIB} ${COMMONS_LIB})
INCLUDE_DIRECTORIES(  build/commons/src/commons/sources)

#Executable
project(main)
add_executable(elastic sources/elastic.c)
IF(APPLE)
target_link_libraries(elastic "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libCommons.dylib" ${ELASTIC_SOLVER_LIB})
ELSEIF("${CMAKE_SYSTEM}" MATCHES "Linux")
target_link_libraries(elastic "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libCommons.so" ${ELASTIC_SOLVER_LIB})
ENDIF()


 



