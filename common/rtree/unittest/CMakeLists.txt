if(BUILD_TESTING)
    add_executable(
        test_geometry_common_rtree
        RTreeTest.cpp
    )
    target_link_libraries(
        test_geometry_common_rtree
        PRIVATE
            geometry_common_rtree
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_rtree
        DISCOVERY_MODE PRE_TEST
    )
endif()
