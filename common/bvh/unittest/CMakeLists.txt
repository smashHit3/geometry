if(BUILD_TESTING)
    add_executable(
        test_geometry_common_bvh
        BvhTest.cpp
    )
    target_link_libraries(
        test_geometry_common_bvh
        PRIVATE
            geometry_common_bvh
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_bvh
        DISCOVERY_MODE PRE_TEST
    )
endif()
