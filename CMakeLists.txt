# pytrt (mclib, imclib, etc.) CMake file
# Seth R. Johnson
# 

cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

project(BFTC)

#------------------------------ VERSION INFO ------------------------------#
set( ${PROJECT_NAME}_MAJOR_VERSION 0 ) 
set( ${PROJECT_NAME}_MINOR_VERSION 0 ) 
set( ${PROJECT_NAME}_PATCH_LEVEL   1 ) 

find_package(Subversion)
if(Subversion_FOUND)
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} MCGSVN)
else(Subversion_FOUND)
  set( MCGSVN_WC_REVISION "---" )
endif(Subversion_FOUND)

# these variables are used to set the correct numbering of shared libraries
set( ${PROJECT_NAME}_SOVERSION
  "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
set( ${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL} (r${MCGSVN_WC_REVISION})" )

#--------------------------- INCLUDE AND PATHS --------------------------#
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools")

#==================== 
# modify the include paths for my own purposes (after checking that leaving
# these out will allow manual specification)
if(APPLE)
  list(APPEND CMAKE_MODULE_PATH "/Users/Shared/_local/share/cmake")
  list(APPEND CMAKE_INCLUDE_PATH "/Users/Shared/_local/include")
  list(APPEND CMAKE_LIBRARY_PATH "/Users/Shared/_local/lib")
else(APPLE)
  list(APPEND CMAKE_MODULE_PATH "$ENV{HOME}/local/share/cmake")
  list(APPEND CMAKE_INCLUDE_PATH "$ENV{HOME}/local/include")
  list(APPEND CMAKE_LIBRARY_PATH "$ENV{HOME}/local/lib")
endif(APPLE)
#==================== 
##########
# Set this directory as an include directory so we correctly get
# mclib/blah.hpp and imclib/blah.hpp
include_directories(${PROJECT_SOURCE_DIR})

# Set default output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#------------------------------ OPTIONS ------------------------------#

# Set the default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

#------------------------------ TARGETS ------------------------------#

#------------------------------  RECURSE ------------------------------#
#recurse into subdirectories
add_subdirectory(src)

#------------------------------ OUTPUT ------------------------------#
message(STATUS "${PROJECT_NAME} version  ${${PROJECT_NAME}_VERSION}")
message(STATUS "System:             ${CMAKE_SYSTEM_NAME} "
                                    "${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Build type:         ${CMAKE_BUILD_TYPE}")
