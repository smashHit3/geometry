add_library(geometry_trapezoid INTERFACE)
target_include_directories(geometry_trapezoid INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(geometry_trapezoid INTERFACE geometry_common)
target_compile_features(geometry_trapezoid INTERFACE cxx_std_17)

add_subdirectory(detail)
add_subdirectory(unittest)
