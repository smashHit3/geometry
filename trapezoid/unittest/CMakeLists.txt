if(BUILD_TESTING)
    add_executable(
        test_geometry_trapezoid
        OctilinearTrapezoidMapTest.cpp
    )
    target_link_libraries(
        test_geometry_trapezoid
        PRIVATE
            geometry_trapezoid
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(test_geometry_trapezoid DISCOVERY_MODE PRE_TEST)
endif()
