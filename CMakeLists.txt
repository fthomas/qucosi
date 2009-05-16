cmake_minimum_required(VERSION 2.6)

project(QuCoSi)

#set(CMAKE_VERBOSE_MAKEFILE ON)

#add_definitions(-Wall -ansi -pedantic)
add_definitions(-Wall -ansi)

add_subdirectory(src)
add_subdirectory(tests)

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
endif (DOXYGEN_FOUND)
