add_library(geometry_boolean INTERFACE)
target_include_directories(geometry_boolean INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(geometry_boolean INTERFACE geometry_common)
target_compile_features(geometry_boolean INTERFACE cxx_std_17)

add_subdirectory(detail)
add_subdirectory(unittest)
