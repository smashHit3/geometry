add_library(geometry_skeleton INTERFACE)
target_include_directories(geometry_skeleton INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(geometry_skeleton INTERFACE geometry_common)
target_compile_features(geometry_skeleton INTERFACE cxx_std_17)

add_subdirectory(detail)
add_subdirectory(unittest)
