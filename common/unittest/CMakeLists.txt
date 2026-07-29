if(BUILD_TESTING)
    add_executable(
        test_geometry_common_base_point
        Main.cpp
        BasePointTest.cpp
        BasePointUtilTest.cpp
    )
    target_link_libraries(
        test_geometry_common_base_point
        PRIVATE
            geometry_common
            GTest::gtest
    )

    include(GoogleTest)
    gtest_discover_tests(test_geometry_common_base_point)
endif()
