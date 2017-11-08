cmake_minimum_required(VERSION 3.5.0)

project(Triangle
   # triangle 1.6 was the last "official" version; i.e. 1.7 should be good to go:
   VERSION 1.7.0
   )

option(TRIANGLE_ENABLE_ACUTE "Enable the aCute extension." ON)

option(BUILD_TESTING "Enable unit tests." ON)
if(BUILD_TESTING)
   enable_testing()
endif()

add_subdirectory(Triangle)
add_subdirectory(examples)
