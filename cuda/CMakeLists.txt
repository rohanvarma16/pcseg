cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(pcl_seg)

find_package(CUDA REQUIRED)
set(CUDA_SEPARABLE_COMPILATION ON)
IF(CUDA_FOUND)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11" )
INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
SET(LIBRARIES ${CUDA_LIBRARIES})
cuda_compile(cuda_methods src/cuda_methods.cu)
cuda_add_library(cuda_methods_lib src/cuda_methods.cu)
ENDIF()

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



set(CMAKE_PREFIX_PATH "~/install_dir/pcl/")
find_package(PCL 1.3 REQUIRED COMPONENTS common io kdtree)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")
add_executable(pcd_cuda src/preprocess.cpp src/segmentation.cpp src/helper.hpp src/pcd_cuda.cpp)
target_link_libraries(pcd_cuda cuda_methods_lib ${CUDA_LIBRARIES} ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${PCL_LIBRARIES})

