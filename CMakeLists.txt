cmake_minimum_required(VERSION 3.12.0 FATAL_ERROR)
project(A_Finite_Element_project LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_STANDARD_REQUIRED "ON")
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")



# Set default build type to Release.
if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions(-DBUILD_TYPE_DEBUG)
endif()


# Locate vtk
find_package(VTK REQUIRED HINTS $ENV{mkVtkPrefix})

# Locate eigen
find_package (Eigen3 3.3 REQUIRED HINTS $ENV{mkEigenPrefix})

# Locate boost
find_package(Boost REQUIRED HINTS $ENV{mkBoostPrefix})
find_package(Boost REQUIRED COMPONENTS program_options )


# Add useful compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wnon-virtual-dtor -fopenmp -fconcepts-diagnostics-depth=3")


set(TARGET src/main.cpp  src/quadrature.cpp src/solver.cpp src/psolver.cpp src/psolver2.cpp src/psolverColoring.cpp ) 
add_executable(${CMAKE_PROJECT_NAME} ${TARGET})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
target_link_libraries (${CMAKE_PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
vtk_module_autoinit(
  TARGETS ${CMAKE_PROJECT_NAME}
  MODULES ${VTK_LIBRARIES}
)