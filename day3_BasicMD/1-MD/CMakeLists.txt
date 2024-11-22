cmake_minimum_required(VERSION 3.15)
project(example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11) 

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(md pybindings.cpp md.cpp sample.cpp)
target_include_directories(md PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})