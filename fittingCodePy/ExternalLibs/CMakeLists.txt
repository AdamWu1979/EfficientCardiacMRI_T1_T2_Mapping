# Build the nifti file format library


subdirs(nifti)

set(LIBRARIES extern_nifti)


# Create the extern_io library

add_library(_extern_ReadWriteImage _extern_ReadWriteImage.h _extern_ReadWriteImage.cpp)

target_link_libraries(_extern_ReadWriteImage ${LIBRARIES})





install(TARGETS _extern_ReadWriteImage        
RUNTIME DESTINATION bin COMPONENT Development
        
LIBRARY DESTINATION lib COMPONENT Development
        
ARCHIVE DESTINATION lib COMPONENT Development
)


install(FILES _extern_ReadWriteImage.h DESTINATION include COMPONENT Development)



