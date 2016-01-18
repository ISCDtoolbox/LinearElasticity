cmake_minimum_required(VERSION 2.8)
project(Elas)
set(CMAKE_BUILD_TYPE "Release")

SET ( COMMON_DIR common )
SET ( SOURCES_DIR sources )

include_directories(${COMMON_DIR} ${SOURCES_DIR})

#Compilation
file(GLOB_RECURSE source_files ${SOURCES_DIR}/* ${COMMON_DIR}/*)
file(GLOB_RECURSE header_files ${SOURCES_DIR}/*.h ${COMMON_DIR}/*.h)

#Find math and pthread libraries
IF (NOT WIN32)
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${M_LIB})
ENDIF()
SET( LIBRARIES "-lpthread" ${LIBRARIES})

#Build Elasticity library
add_library(${PROJECT_NAME} SHARED ${source_files})
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION 5.0.a SOVERSION 5)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES})

#Build main executable
project(main)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/bin")
add_executable(elastic elastic.c)
target_link_libraries(elastic Elas ${LIBRARIES})
