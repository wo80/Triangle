add_executable(triangle-cli
   main.c
   util.c
   )
set_target_properties(triangle-cli PROPERTIES OUTPUT_NAME triangle)

target_link_libraries(triangle-cli PRIVATE triangle::libtriangle)
