SET(QUCOSI_SOURCES
    main.cpp
)

SET(QUCOSI_HEADERS
    Basis.h
    Gate.h
    Qubit.h
    Vector.h
)

INCLUDE_DIRECTORIES(. /usr/include/eigen2)
ADD_EXECUTABLE(qctest ${QUCOSI_SOURCES})
