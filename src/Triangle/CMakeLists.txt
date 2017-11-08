# the libraries in this directory are not meant for direct consumption
# (e.g. there are no proper export definitions)

add_library(triangle
   acute.c
   eps_writer.c
   predicates.c
   triangle.c
   triangle_io.c
   )
add_library(triangle::triangle ALIAS triangle)

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
