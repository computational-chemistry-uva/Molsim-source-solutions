cmake_minimum_required(VERSION 3.15)
project(molsim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hardDisks
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecularDynamics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monteCarlo
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
)

pybind11_add_module(molsim
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pybindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hardDisks/hardDisks.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecularDynamics/md.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecularDynamics/sample.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecularDynamics/thermostats.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monteCarlo/mc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monteCarlo/moves.cpp
)
target_include_directories(molsim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hardDisks
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecularDynamics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monteCarlo
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
)

install(TARGETS molsim LIBRARY DESTINATION .)

add_subdirectory(docs)