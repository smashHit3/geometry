cmake_minimum_required(VERSION 3.20)

project(geometry
    VERSION 0.1.0
    DESCRIPTION "Geometry library"
    LANGUAGES CXX
)

add_executable(geometry main.cpp)
target_compile_features(geometry PRIVATE cxx_std_17)

include(CTest)

add_subdirectory(common)
add_subdirectory(boolean)
add_subdirectory(skeleton)
add_subdirectory(trapezoid)

if(BUILD_TESTING)
    add_test(NAME geometry-smoke COMMAND geometry)
endif()
