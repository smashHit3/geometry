if(BUILD_TESTING)
    add_executable(
        test_geometry_common_base_point
        Main.cpp
        BasePointTest.cpp
        BasePointUtilTest.cpp
        CommonUtilTest.cpp
        SpatialIndexTest.cpp
    )
    target_link_libraries(
        test_geometry_common_base_point
        PRIVATE
            geometry_common
            geometry_common_geometry
            geometry_common_geometryutil
            geometry_common_rtree
            geometry_common_bvh
            geometry_common_uniformgrid
            GTest::gtest
    )

    include(GoogleTest)
    gtest_discover_tests(
        test_geometry_common_base_point
        DISCOVERY_MODE PRE_TEST
    )
endif()
