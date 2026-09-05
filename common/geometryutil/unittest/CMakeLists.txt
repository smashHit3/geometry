if(BUILD_TESTING)
    add_executable(
        test_geometry_common_geometryutil
        BasePointUtilTest.cpp
        CommonUtilTest.cpp
    )
    target_link_libraries(
        test_geometry_common_geometryutil
        PRIVATE
            geometry_common_geometryutil
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_geometryutil
        DISCOVERY_MODE PRE_TEST
    )
endif()
