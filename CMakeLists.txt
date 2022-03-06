cmake_minimum_required(VERSION 3.16)
project(RGMesh)

set(CMAKE_CXX_STANDARD 17)

#add_subdirectory(extern/tetgen)


#find_package(PCL REQUIRED)
find_package(OpenMP)



#include_directories( ${PCL_INCLUDE_DIRS} )
## [C/C++]>[Preprocessor]>[Preprocessor Definitions]
#add_definitions( ${PCL_DEFINITIONS} )
## For Use Not PreCompiled Features
#add_definitions( -DPCL_NO_PRECOMPILE )
## [Linker]>[General]>[Additional Library Directories]
#link_directories( ${PCL_LIBRARY_DIRS} )


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/phtree-cpp)

add_executable(RGMesh RiemanGraph/RiemannGraph.cpp apps/Smoothing.cpp)
target_link_libraries(RGMesh OpenMP::OpenMP_CXX phtree)

#add_executable(MegTriRecons apps/reconstrcut.cpp)
#target_link_libraries(MegTriRecons OpenMP::OpenMP_CXX )
#
#add_executable(MegTriResampling apps/resampling.cpp)
#target_link_libraries(MegTriRecons OpenMP::OpenMP_CXX )

