include_directories(. .. /usr/include/eigen2)

add_executable(testall testall.cpp)
target_link_libraries(testall cppunit)

add_executable(sandbox sandbox.cpp)

add_custom_target(test ${CMAKE_CURRENT_BINARY_DIR}/testall DEPENDS testall)
