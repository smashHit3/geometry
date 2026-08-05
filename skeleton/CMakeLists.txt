add_library(geometry_skeleton detail/skeleton.cpp)
target_include_directories(geometry_skeleton PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(geometry_skeleton PUBLIC geometry_common)
target_compile_features(geometry_skeleton PUBLIC cxx_std_17)

add_subdirectory(detail)
add_subdirectory(unittest)
