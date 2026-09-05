add_library(geometry_common_rtree INTERFACE)
target_include_directories(geometry_common_rtree INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../..")
target_link_libraries(geometry_common_rtree INTERFACE geometry_common)
target_compile_features(geometry_common_rtree INTERFACE cxx_std_17)

add_subdirectory(unittest)
