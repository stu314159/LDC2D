cmake_minimum_required(VERSION 3.5)

enable_language(CUDA)

project(Lid-Driven-Cavity-2D VERSION 0.1 LANGUAGES CUDA CXX)
find_package(CUDA)

add_library(vtk_lib STATIC src/vtk_lib.cxx include/vtk_lib.h)
add_library(lbm_lib STATIC src/lbm_lib.cu include/lbm_lib.h)

add_executable(ldc2D src/ldc2D.cxx)

set_property(TARGET ldc2D PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(ldc2D PUBLIC include PRIVATE ${CUDA_INCLUDE_DIRS})
target_include_directories(vtk_lib PUBLIC include)
target_include_directories(lbm_lib PUBLIC include)

target_link_libraries(ldc2D ${CUDA_LIBRARIES} vtk_lib lbm_lib)
