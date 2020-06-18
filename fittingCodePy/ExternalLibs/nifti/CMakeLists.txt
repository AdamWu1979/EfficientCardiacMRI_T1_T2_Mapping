#-----------------------------------------------------------------------------

add_definitions(-DHAVE_ZLIB)


set(NAME extern_nifti)


find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    #target_link_libraries(${NAME} ${ZLIB_LIBRARIES})
endif()


add_library(${NAME} nifti1.h nifti1_io.c nifti1_io.h znzlib.c znzlib.h)

#target_link_libraries(${NAME} "C:\Program Files (x86)\GnuWin32\lib/zlib.lib")
target_link_libraries(${NAME} ${ZLIB_LIBRARIES})


install(TARGETS ${NAME}
        
RUNTIME DESTINATION bin
        
LIBRARY DESTINATION lib
        
ARCHIVE DESTINATION lib
)

install(FILES nifti1_io.h znzlib.h nifti1.h DESTINATION include COMPONENT Development)


#-----------------------------------------------------------------------------
