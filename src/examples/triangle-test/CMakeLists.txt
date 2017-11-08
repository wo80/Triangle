add_executable(triangle-test
   main.c
   tests.c
   util.c
   )

target_link_libraries(triangle-test PRIVATE triangle::libtriangle)

add_test(
   NAME triangle-test
   COMMAND triangle-test)
