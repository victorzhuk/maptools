enable_language(C)
set(SRCS main.c)
add_executable(tests ${SRCS})
target_link_libraries(tests "maptools")

