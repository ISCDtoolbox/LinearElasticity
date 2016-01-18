cmake_minimum_required(VERSION 2.8)
project(Elas)
set(CMAKE_BUILD_TYPE "Release")

SET( SOURCES_DIR sources )

#Compilation
file(GLOB_RECURSE source_files ${SOURCES_DIR}/*)
file(GLOB_RECURSE header_files ${SOURCES_DIR}/*.h)

#Adds commons
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

#Find math and pthread libraries
IF (NOT WIN32)
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${M_LIB})
ENDIF()
SET( LIBRARIES "-lpthread" ${LIBRARIES})


#Build Elasticity library
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "$ENV{HOME}/lib")
add_library(${PROJECT_NAME} SHARED ${source_files})

add_dependencies(${PROJECT_NAME} commons)

#target_link_libraries(elastic )

include_directories(build/commons/src/commons/sources)

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES})

#Build main executable
project(main)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/bin")
add_executable(elastic sources/elastic.c)
target_link_libraries(elastic Elas ${LIBRARIES} "$ENV{HOME}/lib/libCommons.so")
