enable_language(C)

set(SRCS maptools.h geotasks.c)

add_library(maptools SHARED ${SRCS})
