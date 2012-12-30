######################################################################################
# Packages needed for PIBIC files
######################################################################################

set(CMAKE_C_FLAGS "-Wall -g -fopenmp ${CMAKE_CPP_FLAGS}")
INCLUDE_DIRECTORIES(/home/ataias/workspace/ff/include)

add_executable(edp_06 edp_06.cpp)
target_link_libraries(edp_06 ${Boost_PROGRAM_OPTIONS_LIBRARY} )
install (TARGETS edp_06
  DESTINATION bin)
