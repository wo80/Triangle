if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_subdirectory(triangle-cli)
add_subdirectory(triangle-test)
