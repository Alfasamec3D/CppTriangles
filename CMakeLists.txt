cmake_minimum_required(VERSION 3.13)

project("Triangles" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost REQUIRED)

add_executable(main main.cpp)

add_library(geom lib/geom.cpp)

target_include_directories(geom PUBLIC "include"  PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(main PRIVATE geom)


option(BUILD_TESTS "Build tests" ON)
if (BUILD_TESTS)
    add_subdirectory(tests)
endif()