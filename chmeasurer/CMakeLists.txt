# USING cmake:
#   mkdir build
#   cd build
#   cmake ..
#   make -j

# RUN unit tests:
#   cd build
#   ctest



# OBJECT library requires 2.8.8
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(ch_measurer CXX)
#project(run_tests CXX)

#GLPK
FIND_LIBRARY( GLPK_LIBRARIES
 glpk
 PATHS /usr/lib $ENV{GLPK_DIR}/lib
 DOC "Directory where GLPK lib file is stored"
)
FIND_PATH(GLPK_INCLUDE_DIR
  glpk.h
  PATHS /usr/include/ /usr/include/glpk $ENV{GLPK_DIR}/include
  DOC "Directory where GLPK header files are stored" )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK "Could not find GLPK " GLPK_INCLUDE_DIR GLPK_LIBRARIES)
include_directories(${GLPK_INCLUDE_DIR})
link_directories(${GLPK_LIBRARIES})

#GMP
FIND_LIBRARY( GMP_LIBRARIES
 gmp
 PATHS /usr/lib $ENV{GMP_DIR}/lib
 DOC "Directory where GMP lib file is stored"
)
FIND_PATH(GMP_INCLUDE_DIR
  gmp.h
  PATHS /usr/include/ /usr/include/GMP $ENV{GLPK_DIR}/include
  DOC "Directory where GMP header files are stored" )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP "Could not find GMP " GMP_INCLUDE_DIR GMP_LIBRARIES)
include_directories(${GMP_INCLUDE_DIR})
link_directories(${GMP_LIBRARIES})


find_package(Boost REQUIRED thread)

find_package(CGAL REQUIRED)

find_package(PkgConfig REQUIRED)



include(CMakeDetermineCXXCompiler)

# default to RelWithDebInfo
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(EXTRA_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -pedantic" CACHE STRING "Extra flags used by the compiler during all build types.")
	set(EXTRA_EXE_LINKER_FLAGS "-Wl,--as-needed" CACHE STRING "Extra flags used by the linker.")
	set(EXTRA_EXE_LINKER_FLAGS_RELEASE "-flto" CACHE STRING "Extra flags used by the linker for the Release build type.")
	set(EXTRA_EXE_LINKER_FLAGS_RELWITHDEBINFO "-flto" CACHE STRING "Extra flags used by the linker for the RelWithDebInfo build type.")
else()
	set(EXTRA_CXX_FLAGS "" CACHE STRING "Extra flags used by the compiler during all build types.")
	set(EXTRA_EXE_LINKER_FLAGS "" CACHE STRING "Extra flags used by the linker.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS} -std=c++11 -fPIE -frounding-math -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EXTRA_EXE_LINKER_FLAGS} -pie -frounding-math -fopenmp")   
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${EXTRA_EXE_LINKER_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} ${EXTRA_EXE_LINKER_FLAGS_RELWITHDEBINFO}")

option(VERBOSE "Verbose logging" ON)

if(NOT VERBOSE)
	add_definitions(-DNVERBOSE)
endif()

# compile shared sources only once, and reuse object files in both,
# as they are compiled with the same options anyway

add_library(common OBJECT	
	src/nodes_and_edges.h
)

add_executable(ch_measurer
	src/ch_measurer.cpp	
	$<TARGET_OBJECTS:common>
)
target_link_libraries(ch_measurer ${GLPK_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${Boost_LIBRARIES} )

set(TESTING ON)

if(TESTING)
	enable_testing()

	add_executable(run_tests
		src/run_tests.cpp
		src/unit_tests.cpp
		$<TARGET_OBJECTS:common>
	)
	target_link_libraries(run_tests ${GLPK_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${Boost_LIBRARIES})

	add_test(NAME unit-test
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/src"
		COMMAND $<TARGET_FILE:run_tests>
	)
endif()

