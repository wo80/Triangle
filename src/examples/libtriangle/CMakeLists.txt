set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# this is currently needed as a clutch.
# ... a proper solution is to change the directory layout of internal vs. api lib
configure_file(../../Triangle/triangle.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(../../Triangle/triangle_config.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# force shared library for now, because otherwise check_behavior is defined and clashes with another define in triangle-cli
add_library(triangle-api SHARED
   triangle_api.c
   triangle_helper.c
   )
add_library(Triangle::triangle-api ALIAS triangle-api)

target_link_libraries(triangle-api
   PRIVATE
   Triangle::triangle
   )
target_include_directories(triangle-api
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
   $<INSTALL_INTERFACE:include>
   )

#set_property(TARGET triangle-api PROPERTY VERSION ${Triangle_VERSION})
#set_property(TARGET triangle-api PROPERTY SOVERSION 1)
#set_property(TARGET triangle-api PROPERTY INTERFACE_Triangle_MAJOR_VERSION 3)
#set_property(TARGET triangle-api APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Triangle_MAJOR_VERSION)

include(GenerateExportHeader)
generate_export_header(triangle-api BASE_NAME triangle)

install(
   TARGETS triangle-api
   EXPORT TriangleTargets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   RUNTIME DESTINATION bin
   INCLUDES DESTINATION include
   )

install(
   FILES
   triangle_api.h
   triangle_helper.h
   ${CMAKE_CURRENT_BINARY_DIR}/triangle_export.h
   DESTINATION include
   COMPONENT Devel
   )
