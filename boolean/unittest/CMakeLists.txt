if(BUILD_TESTING)
    add_executable(
        test_geometry_boolean_rect_clipper
        RectClipperTest.cpp
    )
    target_link_libraries(
        test_geometry_boolean_rect_clipper
        PRIVATE
            geometry_boolean
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_boolean_rect_clipper
        DISCOVERY_MODE PRE_TEST
    )
endif()
