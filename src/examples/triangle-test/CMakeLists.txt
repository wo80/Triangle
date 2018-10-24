add_executable(triangle-test
   main.c
   tests.c
   util.c
   )

target_link_libraries(triangle-test PRIVATE Triangle::triangle-api)

add_test(
   NAME triangle-test
   COMMAND triangle-test)

add_custom_command(TARGET triangle-test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy            
        "$<TARGET_FILE_DIR:triangle-api>/$<TARGET_FILE_NAME:triangle-api>"
        "$<TARGET_FILE_DIR:triangle-test>/$<TARGET_FILE_NAME:triangle-api>"
)