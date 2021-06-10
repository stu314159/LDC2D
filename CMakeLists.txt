cmake_minimum_required(VERSION 3.5)

enable_language(CUDA)

project(Lid-Driven-Cavity-2D VERSION 0.1 LANGUAGES CUDA CXX)
find_package(CUDA)

add_library(vtk_lib STATIC vtk_lib.cxx vtk_lib.h)
add_library(lbm_lib STATIC lbm_lib.cu)


include_directories("${CUDA_INCLUDE_DIRS}")

add_executable(ldc2D ldc2D.cxx)

set_property(TARGET ldc2D PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(ldc2D ${CUDA_LIBRARIES} vtk_lib lbm_lib)
