cmake_minimum_required(VERSION 3.0)

project(TIM)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

find_package (Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# find python libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
  matplotlib
  GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
  GIT_TAG        master
  )
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
  FetchContent_Populate(matplotlib)
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

include_directories(SYSTEM ${matplotlib_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)

file(GLOB app_sources ${PROJECT_SOURCE_DIR}/app/*.cpp)

foreach(app ${app_sources})
  get_filename_component(exe ${app} NAME_WE)

  add_executable(${exe} ${app} ${sources})

  target_link_libraries(${exe} PRIVATE ${ARMADILLO_LIBRARIES} ${PYTHON_LIBRARIES} Python3::NumPy)

endforeach(app ${app_sources})
