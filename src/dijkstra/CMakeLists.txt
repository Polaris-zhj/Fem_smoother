cmake_minimum_required(VERSION 3.0.2)
project(dijkstra)

SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_definitions(-DPROJECT_PATH="${PROJECT_SOURCE_DIR}")
add_compile_options(-std=c++11)

file(GLOB_RECURSE PYTHON2.7_LIB "/usr/lib/python2.7/config-x86_64-linux-gnu/*.so") 
set(PYTHON2.7_INLCUDE_DIRS "/usr/include/python2.7") #指定头文件路径

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)
find_package(OsqpEigen REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(tinyspline_ros REQUIRED)


catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES reed_sheep
#  CATKIN_DEPENDS roscpp  std_msgs
#  DEPENDS system_lib
)


include_directories(
  ${include}
  ${include}/tinyspline_ros/
  ${catkin_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
  ${PROJECT_SOURCE_DIR}/
)

include_directories(include
        ${PYTHON2.7_INLCUDE_DIRS}       #添加头文件到工程
)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

add_executable(dijkstra src/main.cpp src/dijkstra.cpp src/fem_smoother.cpp)

target_link_libraries(dijkstra ${PYTHON2.7_LIB} tinyspline_ros_cpp)  #添加库文件到工程
target_link_libraries(dijkstra OsqpEigen::OsqpEigen yaml-cpp)

