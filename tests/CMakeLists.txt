set(tests_SOURCES
    main.cpp
)

set(tests_HEADERS
    AlgorithmsTest.h
    GateTest.h
    QubitTest.h
    VectorTest.h
)

include_directories(. .. /usr/include/eigen2)
add_executable(tests ${tests_SOURCES})
target_link_libraries(tests cppunit)
