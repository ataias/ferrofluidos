# Author: Ataias Pereira Reis
# The files in this project require boost and eigen libraries
# The compilation is made with CMake that creates a makefile
# to guide the compilation

cmake_minimum_required(VERSION 2.6)
project(ff)
#The version number
set (EDP_VERSION_MAJOR 1)
set (EDP_VERSION_MINOR 0)

#CONDITIONALS
#Check of uncheck each as you wish
SET(COMPILE_UNIT_TEST_MATRIX_FILES_IO 1)
SET(COMPILE_UNIT_TEST_POISSON 1)
SET(COMPILE_OTHERS 1)
SET(COMPILE_MAIN_FILE 0) #not finished yet
SET(USING_PYTHON 1)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

#Searches for boost
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    FIND_PACKAGE(Boost COMPONENTS filesystem system date_time python)
    message("Include dirs of boost: " ${Boost_INCLUDE_DIRS} )
    message("Libs of boost: " ${Boost_LIBRARIES} )
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ELSE(Boost_FOUND)
MESSAGE(STATUS "Could not find boost")
ENDIF()

#Searches for Eigen
find_package (Eigen3 3.1.0)
IF (EIGEN3_FOUND) 
    INCLUDE_DIRECTORIES(${Eigen3_INCLUDE_DIRS})
    message(STATUS "Eigen ${EIGEN3_VERSION} found" )
ELSE(EIGEN3_FOUND) 
    message(STATUS "Eigen3 not found" )
ENDIF()

#Python and wrapping
IF(USING_PYTHON)
find_package(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
message("Include dirs of Python: " ${PYTHON_INCLUDE_DIRS} )
message("Libs of Python: " ${PYTHON_LIBRARIES} )
ENDIF()

add_subdirectory (src)
