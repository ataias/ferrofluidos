######################################################################################
# Packages needed for PIBIC files
######################################################################################

set(CMAKE_C_FLAGS "-Wall -g -fopenmp ${CMAKE_CPP_FLAGS}")
INCLUDE_DIRECTORIES(/home/ataias/workspace/ff/include)

IF(COMPILE_UNIT_TEST_MATRIX_FILES_IO)
add_executable(MatrixFilesIO_test MatrixFilesIO.cpp MatrixFilesIO_test.cpp)
target_link_libraries(MatrixFilesIO_test ${Boost_PROGRAM_OPTIONS_LIBRARY} )
install (TARGETS MatrixFilesIO_test DESTINATION bin)
ENDIF(COMPILE_UNIT_TEST_MATRIX_FILES_IO)

IF(COMPILE_UNIT_TEST_POISSON)
add_executable(Poisson_test MatrixFilesIO.cpp Poisson.cpp Poisson_test.cpp)
target_link_libraries(Poisson_test ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
install (TARGETS Poisson_test DESTINATION bin)
ENDIF(COMPILE_UNIT_TEST_POISSON)

IF(COMPILE_OTHERS)
add_executable(SimpleWriteAndRead SimpleTestReadAndWrite.cpp)
target_link_libraries(SimpleWriteAndRead ${Boost_PROGRAM_OPTIONS_LIBRARY} )
install (TARGETS SimpleWriteAndRead DESTINATION bin)
ENDIF(COMPILE_OTHERS)

IF(USING_PYTHON)
SET(BUILD_SHARED_LIBS 1)
add_library(fooclass SHARED poisson2py.cpp)
target_link_libraries(fooclass ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
install (TARGETS fooclass DESTINATION bin)

add_library(poisson SHARED Poisson.cpp)
target_link_libraries(poisson ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
install (TARGETS poisson DESTINATION bin)
ENDIF(USING_PYTHON)