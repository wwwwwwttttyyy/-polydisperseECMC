cmake_minimum_required(VERSION 3.10)
project(ECMC_Tests)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../include)

# Source files
set(ECMC_SOURCES
    ../src/system.cpp
    ../src/cell.cpp
    ../src/geometry.cpp
    ../src/collision.cpp
    ../src/ecmc.cpp
    ../src/para.cpp
)

# Test executable
add_executable(test_ecmc
    test_ecmc.cpp
    ${ECMC_SOURCES}
)

# IO test executable
add_executable(test_io
    test_io.cpp
    ${ECMC_SOURCES}
)

# Link math library on Linux
if(UNIX AND NOT APPLE)
    target_link_libraries(test_ecmc m)
    target_link_libraries(test_io m)
endif()
