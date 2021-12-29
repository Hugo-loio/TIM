cmake_minimum_required(VERSION 3.0)

project(TIM)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

find_package (Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

file(GLOB base_sources ${PROJECT_SOURCE_DIR}/base/src/*.cpp)
file(GLOB base_headers ${PROJECT_SOURCE_DIR}/base/include/*.h)

file(GLOB models_sources ${PROJECT_SOURCE_DIR}/models/src/*.cpp)
file(GLOB models_headers ${PROJECT_SOURCE_DIR}/models/include/*.h)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/plots)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/base/include ${PROJECT_SOURCE_DIR}/models/include)

file(GLOB app_sources ${PROJECT_SOURCE_DIR}/app/*.cpp)

foreach(app ${app_sources})
  get_filename_component(exe ${app} NAME_WE)

  add_executable(${exe} ${app} ${base_sources} ${models_sources})

  target_link_libraries(${exe} PRIVATE ${ARMADILLO_LIBRARIES})

endforeach(app ${app_sources})
