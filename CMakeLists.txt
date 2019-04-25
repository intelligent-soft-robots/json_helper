cmake_minimum_required(VERSION 2.8.3)

project(json_helper)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wl,--no-as-needed")

find_package(catkin REQUIRED COMPONENTS )

find_package(Armadillo REQUIRED)

catkin_package(
  LIBRARIES json_helper
  INCLUDE_DIRS include
  CATKIN_DEPENDS
)

include_directories(include ${catkin_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})

add_library(json_helper
  src/json_helper.cpp
)

target_link_libraries(json_helper ${ARMADILLO_LIBRARIES})
