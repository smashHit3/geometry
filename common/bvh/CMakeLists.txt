add_library(geometry_common_bvh INTERFACE)
target_include_directories(geometry_common_bvh INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../..")
target_link_libraries(geometry_common_bvh INTERFACE geometry_common)
target_compile_features(geometry_common_bvh INTERFACE cxx_std_17)

add_subdirectory(unittest)
