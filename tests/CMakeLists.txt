set(testall_SOURCES
    testall.cpp
)

include_directories(. .. /usr/include/eigen2)
add_executable(testall ${testall_SOURCES})
target_link_libraries(testall cppunit)
