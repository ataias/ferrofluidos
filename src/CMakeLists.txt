######################################################################################
# Packages needed for PIBIC files
######################################################################################

set(CMAKE_C_FLAGS "-fopenmp -Wcpp ${CMAKE_CPP_FLAGS}")
INCLUDE_DIRECTORIES(/home/ataias/workspace/ff/include)

IF(USING_PYTHON)
SET(BUILD_SHARED_LIBS 1)

add_library(navierstokes SHARED NavierStokes.cpp)
target_link_libraries(navierstokes ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
install (TARGETS navierstokes DESTINATION bin)

add_library(poisson SHARED Poisson.cpp)
target_link_libraries(poisson ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
install (TARGETS poisson DESTINATION bin)
ENDIF(USING_PYTHON)
