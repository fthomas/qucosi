set(QUCOSI_SOURCES
    main.cpp
)

set(QUCOSI_HEADERS
    Aux.h
    Basis.h
    Gate.h
    Qubit.h
    Vector.h
)

include_directories(. /usr/include/eigen2)
add_executable(qctest ${QUCOSI_SOURCES})
