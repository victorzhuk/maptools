cmake_minimum_required(VERSION 2.6)
project(maptools C)

include_directories("${PROJECT_SOURCE_DIR}/src")

set(CMAKE_C_FLAGS "-lm")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist/bin)

# lib
add_subdirectory(src)

# tests
add_subdirectory(tests)
