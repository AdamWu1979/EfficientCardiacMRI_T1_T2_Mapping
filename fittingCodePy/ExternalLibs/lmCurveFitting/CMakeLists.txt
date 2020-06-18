#-----------------------------------------------------------------------------

add_definitions(-DHAVE_LMCURVEFIT)

set(NAME lmCurveFitting)

add_library(${NAME} _curveFitting.h _curveFitting.cpp)

install(TARGETS ${NAME}
        
RUNTIME DESTINATION bin
        
LIBRARY DESTINATION lib
        
ARCHIVE DESTINATION lib
)

install(FILES _curveFitting.h DESTINATION include COMPONENT Development)


#-----------------------------------------------------------------------------
