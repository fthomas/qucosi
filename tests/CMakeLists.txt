SET(tests_SOURCES
    main.cpp
)

SET(tests_HEADERS
    AlgorithmsTest.h
    VectorTest.h
)

INCLUDE_DIRECTORIES(. ../src /usr/include/eigen2)
ADD_EXECUTABLE(tests ${tests_SOURCES})
TARGET_LINK_LIBRARIES(tests cppunit)
