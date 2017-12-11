add_executable(triangle-cli
   main.c
   util.c
   )
set_target_properties(triangle-cli PROPERTIES OUTPUT_NAME triangle)

target_link_libraries(triangle-cli PRIVATE Triangle::triangle-api)

install(
   TARGETS triangle-cli
   EXPORT TriangleTargets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   RUNTIME DESTINATION bin
   INCLUDES DESTINATION include
   COMPONENT bin
   )
