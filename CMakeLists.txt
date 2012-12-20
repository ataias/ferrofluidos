cmake_minimum_required(VERSION 2.6)
project(ff)
#The version number
set (EDP_VERSION_MAJOR 1)
set (EDP_VERSION_MINOR 0)

add_subdirectory (src)
#find_package (eigen)
#if(EIGEN_FOUND) 
#    include_directories(${EIGEN_INCLUDE_DIRS})
#    target_link_libraries(edp eigen)
#else(EIGEN_FOUND) 
#    message(STATUS "Eigen not found" )
#endif(EIGEN_FOUND)
