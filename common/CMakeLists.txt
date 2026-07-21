add_library(geometry_common INTERFACE)
target_include_directories(geometry_common INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_compile_features(geometry_common INTERFACE cxx_std_17)

add_subdirectory(detail)
add_subdirectory(unittest)
