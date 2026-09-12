cmake_minimum_required(VERSION 3.20)

project(geometry
    VERSION 0.1.0
    DESCRIPTION "Geometry library and interactive visualizer"
    LANGUAGES CXX
)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(geometry
    gui/main.cpp
    gui/GeometryCanvas.cpp
    gui/MainWindow.cpp
)
target_compile_features(geometry PRIVATE cxx_std_17)
target_link_libraries(geometry PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)

include(CTest)

if(BUILD_TESTING)
    find_package(GTest QUIET)

    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        FetchContent_MakeAvailable(googletest)
    endif()
endif()

add_subdirectory(common)
add_subdirectory(boolean)
add_subdirectory(skeleton)
add_subdirectory(trapezoid)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/problems/CMakeLists.txt")
    add_subdirectory(problems)
endif()
