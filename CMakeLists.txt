cmake_minimum_required(VERSION 3.16)
project(blackhole CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Apple clang has no bundled OpenMP runtime; point CMake at Homebrew's libomp.
if(APPLE AND EXISTS /opt/homebrew/opt/libomp)
  set(OpenMP_ROOT /opt/homebrew/opt/libomp)
endif()
find_package(OpenMP REQUIRED)

add_executable(blackhole src/main.cpp)
target_include_directories(blackhole PRIVATE src third_party)
target_link_libraries(blackhole PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(blackhole PRIVATE -O3 -Wall -Wextra)
