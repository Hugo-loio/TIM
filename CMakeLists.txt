cmake_minimum_required(VERSION 3.0)

project(TIM)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

find_package (Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/plots)

include_directories(SYSTEM ${matplotlib_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)

file(GLOB app_sources ${PROJECT_SOURCE_DIR}/app/*.cpp)

foreach(app ${app_sources})
  get_filename_component(exe ${app} NAME_WE)

  add_executable(${exe} ${app} ${sources})

  target_link_libraries(${exe} PRIVATE ${ARMADILLO_LIBRARIES})

endforeach(app ${app_sources})
