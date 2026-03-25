cmake_minimum_required(VERSION 3.15)

if(NOT WIN32)
    find_program(INTEL_CXX_COMPILER "icpc")
    
    if(INTEL_CXX_COMPILER)
        message(STATUS "Found icpc at ${INTEL_CXX_COMPILER}, forcing it as the compiler.")
        set(CMAKE_CXX_COMPILER "${INTEL_CXX_COMPILER}")
        find_program(INTEL_C_COMPILER "icc")
        if(INTEL_C_COMPILER)
            set(CMAKE_C_COMPILER "${INTEL_C_COMPILER}")
        endif()
    else()
        message(STATUS "icpc not found, falling back to default compiler (g++).")
    endif()
endif()


project(ECMC VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

set(ECMC_SOURCES
    src/system.cpp
    src/cell.cpp
    src/collision.cpp
    src/ecmc.cpp
    src/para.cpp
)

add_library(ecmc_lib STATIC ${ECMC_SOURCES})
target_include_directories(ecmc_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(ecmc_main main.cpp)
target_link_libraries(ecmc_main ecmc_lib)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        add_compile_options(/O2)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        message(STATUS "Intel Compiler detected. Using -xHost optimization.")
        add_compile_options(-O3 -xHost)
    else()
        add_compile_options(-O3 -march=native)
    endif()
endif()
