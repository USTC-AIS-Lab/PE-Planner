cmake_minimum_required(VERSION 3.0)

if(APPLE)
    set(MACOS 1)
elseif(UNIX AND NOT APPLE)
    set(LINUX 1)
endif()

project(pe-planner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)

if (LINUX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma -ffast-math -march=native")
    add_definitions(-DLINUX)
    find_package(Python3 COMPONENTS Development NumPy)
    find_package(Eigen3 REQUIRED)
    find_package(catkin REQUIRED COMPONENTS roslib roscpp tf)
    include_directories(./third/ ${EIGEN3_INCLUDE_DIR} /usr/lib/include ${catkin_INCLUDE_DIRS} src/ modules/ /usr/include/pcl-1.10/)
    link_directories(/usr/local/lib)
elseif(MACOS)
    find_package(Python3 COMPONENTS Development NumPy)
    find_package(Eigen3 REQUIRED)
    find_package(nlopt REQUIRED)
    include_directories(./third/ ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIRS} src/ modules/)
    link_directories(${NLOPT_LIBRARY_DIRS})
endif()

if (LINUX)
    add_executable(planner_px4
        src/planner_px4_main.cpp
        modules/mpcc/nominal_mpcc.cpp
        modules/map/map.cpp
        modules/ros_interface/ros_interface.cpp
        modules/kinodynamic_astar/kinodynamic_astar.cpp
        modules/bspline_opt/bspline_optimizer.cpp
        modules/px4_interface/px4_interface.cpp)
    target_link_libraries(planner_px4
        Python3::Python Python3::NumPy nlopt ${catkin_LIBRARIES}
    )

    add_executable(map_viewer
        src/map_viewer.cpp
        modules/map/map.cpp
        modules/ros_interface/ros_interface.cpp)
    target_link_libraries(map_viewer
        ${catkin_LIBRARIES}
    )

    add_executable(map_generator
        src/map_generator.cpp
        modules/map/map.cpp
        modules/ros_interface/ros_interface.cpp)
    target_link_libraries(map_generator
        ${catkin_LIBRARIES}
    )
endif()
