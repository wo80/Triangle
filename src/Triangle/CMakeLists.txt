# the libraries in this directory are not meant for direct consumption
# (e.g. there are no proper export definitions)

add_library(triangle
   acute.c
   eps_writer.c
   predicates.c
   triangle.c
   triangle_io.c
   )
add_library(Triangle::triangle ALIAS triangle)

set_target_properties(triangle
   PROPERTIES
   CXX_VISIBILITY_PRESET hidden
   POSITION_INDEPENDENT_CODE TRUE
   )

target_link_libraries(triangle PRIVATE -lm)
target_include_directories(triangle
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>
   )

if(NOT DEFINED TRIANGLE_ENABLE_ACUTE)
   message(STATUS "Building without aCute...")
   target_compile_options(triangle PRIVATE -DNO_ACUTE)
else()
   message(STATUS "Building with aCute...")
endif()

set_property(TARGET triangle PROPERTY VERSION ${Triangle_VERSION})
set_property(TARGET triangle PROPERTY SOVERSION 1)

configure_file(triangle_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/triangle_version.h @ONLY)

install(
   TARGETS triangle
   EXPORT TriangleTargets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   RUNTIME DESTINATION bin
   INCLUDES DESTINATION include
   )

install(
   FILES
   acute.h
   eps_writer.h
   predicates.h
   triangle_config.h
   triangle.h
   ${CMAKE_CURRENT_BINARY_DIR}/triangle_version.h
   DESTINATION include
   COMPONENT Devel
   )
