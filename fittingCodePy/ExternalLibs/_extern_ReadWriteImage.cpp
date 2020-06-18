/*
 *  _extern_ReadWriteImage.cpp
 *
 *
 *  Created by Fani Deligianni on 09/04/2016.
 *  Copyright (c) 2016, InfoReach Consultancy. All rights reserved.
 *  InfoReach Consultancy
 *
 */

#ifndef _EXTERN_READWRITEIMAGE_CPP
#define _EXTERN_READWRITEIMAGE_CPP

#include "_extern_ReadWriteImage.h"
//#include "_extern_tools.h"

/* *************************************************************** */
void extern_fix_filename(nifti_image* image, const char *filename)
{
    std::string name(filename);
    name.append("\0");
    // Free the char arrays if already allocated
    if(image->fname) free(image->fname);
    if(image->iname) free(image->iname);
    // Allocate the char arrays
    image->fname = (char *)malloc((name.size()+1)*sizeof(char));
    image->iname = (char *)malloc((name.size()+1)*sizeof(char));
    // Copy the new name in the char arrays
    strcpy(image->fname,name.c_str());
    strcpy(image->iname,name.c_str());
    // Returns at the end of the function
    return;
}
/* *************************************************************** */
int extern_io_checkFileFormat(const char *filename)
{
    // Nifti format is used by default
    // Check the extention of the provided filename
    std::string b(filename);
    if(b.find( ".nii.gz") != std::string::npos)
        return NR_NII_FORMAT;
    else if(b.find( ".nii") != std::string::npos)
        return NR_NII_FORMAT;
    else if(b.find( ".hdr") != std::string::npos)
        return NR_NII_FORMAT;
    else if(b.find( ".img.gz") != std::string::npos)
        return NR_NII_FORMAT;
    else if(b.find( ".img") != std::string::npos)
        return NR_NII_FORMAT;
    else fprintf(stderr, "[Read WARNING]: Filename extension is not compatible - the Nifti library is used by default\n");

    return NR_NII_FORMAT;
}
/* *************************************************************** */
nifti_image *extern_io_ReadImageFile(const char *filename)
{
    // First read the fileformat in order to use the correct library
    int fileFormat=extern_io_checkFileFormat(filename);

    // Create the nifti image pointer
    nifti_image *image=NULL;

    // Read the image and convert it to nifti format if required
    switch(fileFormat){
    case NR_NII_FORMAT:
        image=nifti_image_read(filename,true);
        extern_fix_filename(image,filename);
        break;
    }
    //extern_checkAndCorrectDimension(image);

    // Return the nifti image
    return image;
}
/* *************************************************************** */
nifti_image *extern_io_ReadImageHeader(const char *filename)
{
    // First read the fileformat in order to use the correct library
    int fileFormat=extern_io_checkFileFormat(filename);

    // Create the nifti image pointer
    nifti_image *image=NULL;

    // Read the image and convert it to nifti format if required
    switch(fileFormat){
    case NR_NII_FORMAT:
        image=nifti_image_read(filename,false);
        break;
    }
    //extern_checkAndCorrectDimension(image);

    // Return the nifti image
    return image;
}
/* *************************************************************** */
void extern_io_WriteImageFile(nifti_image *image, const char *filename)
{
    // First read the fileformat in order to use the correct library
    int fileFormat=extern_io_checkFileFormat(filename);

    // Convert the image to the correct format if required, set the filename and save the file
    switch(fileFormat){
    case NR_NII_FORMAT:
        nifti_set_filenames(image,filename,0,0);
        nifti_image_write(image);
        break;
    }

    // Return
    return;
}
/* *************************************************************** */

#endif
