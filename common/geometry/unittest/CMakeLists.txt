if(BUILD_TESTING)
    add_executable(
        test_geometry_common_geometry
        BasePointTest.cpp
        AabbTest.cpp
    )
    target_link_libraries(
        test_geometry_common_geometry
        PRIVATE
            geometry_common_geometry
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_geometry
        DISCOVERY_MODE PRE_TEST
    )
endif()
