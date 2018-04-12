cmake_minimum_required(VERSION 3.5.0)

project(Triangle
   # triangle 1.6 was the last "official" version; i.e. 1.7 should be good to go:
   VERSION 1.7.0
   )

option(TRIANGLE_ENABLE_ACUTE "Enable the aCute extension." ON)
option(BUILD_SHARED_LIBS "Build libraries as shared libraries." ON)
option(BUILD_TESTING "Enable unit tests." ON)
if(BUILD_TESTING)
   enable_testing()
endif()
option(BUILD_EXAMPLES "Enable examples." ON)
if(BUILD_EXAMPLES)
   add_subdirectory(examples)
endif()

add_subdirectory(Triangle)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
   "${CMAKE_CURRENT_BINARY_DIR}/Triangle/TriangleConfigVersion.cmake"
   VERSION ${Triangle_VERSION}
   COMPATIBILITY AnyNewerVersion
   )

export(EXPORT TriangleTargets
   FILE "${CMAKE_CURRENT_BINARY_DIR}/Triangle/TriangleTargets.cmake"
   NAMESPACE Triangle::
   )
configure_file(cmake/TriangleConfig.cmake
   "${CMAKE_CURRENT_BINARY_DIR}/Triangle/TriangleConfig.cmake"
   COPYONLY
   )
install(EXPORT TriangleTargets
   FILE TriangleTargets.cmake
   NAMESPACE Triangle::
   DESTINATION lib/cmake/Triangle
   )
install(
   FILES
   cmake/TriangleConfig.cmake
   "${CMAKE_CURRENT_BINARY_DIR}/Triangle/TriangleConfigVersion.cmake"
   DESTINATION lib/cmake/Triangle
   COMPONENT Devel
   )

# uninstall target
if(NOT TARGET uninstall)
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
		"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
		IMMEDIATE @ONLY
		)
	add_custom_target(uninstall
		COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
		)
endif()
