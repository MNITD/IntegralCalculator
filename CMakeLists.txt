cmake_minimum_required(VERSION 3.6)
project(IntegralCalculator)

#SET(CMAKE_CXX_FLAGS  "-pthread")

set(CMAKE_CXX_STANDARD 11)

#set(SOURCE_FILES main.cpp ${GCC_COVERAGE_LINK_FLAGS})
set(SOURCE_FILES main.cpp)
add_executable(IntegralCalculator ${SOURCE_FILES})