cmake_minimum_required(VERSION 3.15)
project(ECMC VERSION 1.0.0 LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler options
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Optimization options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        add_compile_options(/O2)
    else()
        add_compile_options(-O3 -march=native)
    endif()
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
set(ECMC_SOURCES
    src/system.cpp
    src/cell.cpp
    src/geometry.cpp
    src/collision.cpp
    src/ecmc.cpp
    src/para.cpp
)

# Create library
add_library(ecmc_lib STATIC ${ECMC_SOURCES})
target_include_directories(ecmc_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Main program
add_executable(ecmc_main main.cpp)
target_link_libraries(ecmc_main ecmc_lib)

# Installation
install(TARGETS ecmc_main DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# Optional: Build tests
option(BUILD_TESTS "Build test programs" ON)
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

# Print configuration information
message(STATUS "ECMC Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
