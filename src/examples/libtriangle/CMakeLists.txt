set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# this is currently needed as a clutch.
# ... a proper solution is to change the directory layout of internal vs. api lib
configure_file(../../Triangle/triangle.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(../../Triangle/triangle_config.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# force shared library for now, because otherwise check_behavior is defined and clashes with another define in triangle-cli
add_library(libtriangle SHARED
   triangle_api.c
   triangle_helper.c
   )
add_library(triangle::libtriangle ALIAS libtriangle)

target_link_libraries(libtriangle
   PRIVATE
   triangle::triangle
   )
target_include_directories(libtriangle
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
   $<INSTALL_INTERFACE:include>
   )

#set_property(TARGET libtriangle PROPERTY VERSION ${Triangle_VERSION})
#set_property(TARGET libtriangle PROPERTY SOVERSION 1)
#set_property(TARGET libtriangle PROPERTY INTERFACE_Triangle_MAJOR_VERSION 3)
#set_property(TARGET libtriangle APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Triangle_MAJOR_VERSION)

include(GenerateExportHeader)
generate_export_header(libtriangle BASE_NAME triangle)

install(
   FILES
   triangle_api.h
   triangle_helper.h
   ${CMAKE_CURRENT_BINARY_DIR}/triangle_export.h
   ${CMAKE_CURRENT_BINARY_DIR}/triangle.h
   DESTINATION include
   COMPONENT Devel
   )
