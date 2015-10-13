cmake_minimum_required(VERSION 2.8)

project(Elas)

set(CMAKE_BUILD_TYPE "Release")


#Options C++ (-Wno-deprecated: remove glut deprecated)

#set (CMAKE_C_FLAGS "-std=-03")

#Include files
#include_directories("include/")
SET ( COMMON_DIR common )
SET ( SOURCES_DIR sources )
include_directories(${COMMON_DIR} ${SOURCES_DIR})

#Compilation
#Recuperation de la liste des fichiers .c et .h
file(GLOB_RECURSE source_files ${SOURCES_DIR}/* ${COMMON_DIR}/*)
file(GLOB_RECURSE header_files ${SOURCES_DIR}/*.h ${COMMON_DIR}/*.h)


#Build SUscElas library
add_library(Elas SHARED ${source_files})
SET(LIBRARY_OUTPUT_PATH $ENV{HOME}/lib)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION 5.0.a SOVERSION 5)

# Install library in /usr/local or equivalent
INSTALL(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION LIBRARY
  DESTINATION lib)
INSTALL(FILES ${header_files} DESTINATION include/SUscElas)
# Install library in project directory
INSTALL(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION LIBRARY
  DESTINATION ${CMAKE_SOURCE_DIR}/lib)
INSTALL(FILES ${header_files} DESTINATION ${CMAKE_SOURCE_DIR}/include/SUscElas)

#Build main executable
project(main)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/bin")
add_executable(sulastic elastic.c)
target_link_libraries(sulastic Elas)
