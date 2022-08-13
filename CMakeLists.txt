cmake_minimum_required(VERSION 2.8)

project(TIM)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# Compiler options
add_compile_options(-std=c++11 -pthread -g)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# Armadillo
find_package (Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# MPI
find_package (MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# base library
file(GLOB base_sources ${PROJECT_SOURCE_DIR}/base/src/*.cpp)
add_library(base STATIC ${base_sources})
target_include_directories(base PUBLIC ${PROJECT_SOURCE_DIR}/base/include)
# target_link_libraries(base PRIVATE ${MPI_CXX_LIBRARIES})

# models library
file(GLOB models_sources ${PROJECT_SOURCE_DIR}/models/src/*.cpp)
add_library(models STATIC ${models_sources})
target_include_directories(models PUBLIC ${PROJECT_SOURCE_DIR}/models/include)
target_link_libraries(models PRIVATE base)

# create special needed directories
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)
#file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/aux)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/plots)

# compile main functions
file(GLOB app_sources ${PROJECT_SOURCE_DIR}/app/*.cpp)
foreach(app ${app_sources})
  get_filename_component(exe ${app} NAME_WE)
  add_executable(${exe} ${app})
  target_link_libraries(${exe} PRIVATE ${ARMADILLO_LIBRARIES} PRIVATE base PRIVATE models PRIVATE ${MPI_CXX_LIBRARIES})
endforeach(app ${app_sources})
