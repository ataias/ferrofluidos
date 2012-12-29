######################################################################################
# Packages needed for PIBIC files
######################################################################################

set(CMAKE_C_FLAGS "-Wall -fopenmp")
INCLUDE_DIRECTORIES(/home/ataias/workspace/ff/include/Eigen)
INCLUDE_DIRECTORIES(/home/ataias/workspace/ff/include/OwnMath)

add_executable(edp_01 edp_01.cpp)
add_executable(edp_02 edp_02.cpp)
add_executable(edp_03 edp_03.cpp)
add_executable(edp_04 edp_04.cpp)
add_executable(edp_05 edp_05.cpp)
add_executable(edp_06 edp_06.cpp)
install (TARGETS edp_01 edp_02 edp_03 edp_04 edp_05 edp_06
  DESTINATION bin)
