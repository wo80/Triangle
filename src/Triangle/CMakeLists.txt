# the libraries in this directory are not meant for direct consumption
# (e.g. there are no proper export definitions)

add_library(triangle_internal STATIC
   acute.c
   eps_writer.c
   predicates.c
   triangle.c
   triangle_io.c
   )
add_library(triangle::triangle_internal ALIAS triangle_internal)

set_target_properties(triangle_internal
   PROPERTIES
   CXX_VISIBILITY_PRESET hidden
   POSITION_INDEPENDENT_CODE TRUE
   )

target_link_libraries(triangle_internal PRIVATE -lm)
target_include_directories(triangle_internal
   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>
   )

if(NOT DEFINED TRIANGLE_ENABLE_ACUTE)
   message(STATUS "Building without aCute...")
   target_compile_options(triangle_internal PRIVATE -DNO_ACUTE)
else()
   message(STATUS "Building with aCute...")
endif()
