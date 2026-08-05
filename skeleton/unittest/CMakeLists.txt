if(BUILD_TESTING)
    add_executable(
        test_geometry_skeleton
        MedialAxisTest.cpp
    )
    target_link_libraries(
        test_geometry_skeleton
        PRIVATE
            geometry_skeleton
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(test_geometry_skeleton DISCOVERY_MODE PRE_TEST)
endif()
