cmake_minimum_required(VERSION 3.20)

project(geometry
    VERSION 0.1.0
    DESCRIPTION "Geometry library and interactive visualizer"
    LANGUAGES CXX
)

find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets)
if(QT_FOUND)
    find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(QT_FOUND)
    add_executable(geometry
        gui/main.cpp
        gui/GeometryCanvas.cpp
        gui/MainWindow.cpp
    )
    target_compile_features(geometry PRIVATE cxx_std_17)
    target_link_libraries(geometry PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)
endif()

include(CTest)

if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(common)
add_subdirectory(boolean)
add_subdirectory(skeleton)
add_subdirectory(trapezoid)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/problems/CMakeLists.txt")
    add_subdirectory(problems)
endif()
