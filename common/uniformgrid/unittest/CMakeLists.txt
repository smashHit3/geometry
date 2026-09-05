if(BUILD_TESTING)
    add_executable(
        test_geometry_common_uniformgrid
        UniformGridTest.cpp
    )
    target_link_libraries(
        test_geometry_common_uniformgrid
        PRIVATE
            geometry_common_uniformgrid
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_uniformgrid
        DISCOVERY_MODE PRE_TEST
    )
endif()
