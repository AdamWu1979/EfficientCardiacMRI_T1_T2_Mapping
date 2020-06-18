#-----------------------------------------------------------------------------

add_definitions(-DHAVE_PYIMAGE)

set(NAME pyImageDy)

set(NAMEST pyImageSt)

set(BOOST_ROOT "C:/Fani/libraries/boost_1_61_0/")
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:/Fani/libraries/boost_1_61_0/")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Fani/libraries/boost_1_61_0/lib")

find_package (PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories ("${PROJECT_SOURCE_DIR}/ExternalLibs/lmCurveFitting")

FIND_PACKAGE(Boost 1.61.0)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" "C:/Program Files (x86)/Python35-32/include")
  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost COMPONENTS python)

  ADD_LIBRARY(${NAME} SHARED Images.h Images.cxx)
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES} lmCurveFitting)
  INSTALL(TARGETS ${NAME}
        
	RUNTIME DESTINATION bin
        
	LIBRARY DESTINATION lib
        
	ARCHIVE DESTINATION lib
)

	install(FILES Images.h DESTINATION include COMPONENT Development)

  ADD_LIBRARY(${NAMEST} Images.h Images.cxx)
  TARGET_LINK_LIBRARIES(${NAMEST} ${Boost_LIBRARIES} lmCurveFitting)
  INSTALL(TARGETS ${NAMEST}
        
	RUNTIME DESTINATION bin
        
	LIBRARY DESTINATION lib
        
	ARCHIVE DESTINATION lib
)



ELSEIF(NOT BOOST_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()




#-----------------------------------------------------------------------------
