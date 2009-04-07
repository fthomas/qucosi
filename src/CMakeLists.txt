SET(QUCOSI_SOURCES
    main.cpp
)

SET(QUCOSI_HEADERS
    Base.h
    MultiQubit.h
    Qubit.h
)

INCLUDE_DIRECTORIES(. /usr/include/eigen2)
ADD_EXECUTABLE(qctest ${QUCOSI_SOURCES})
