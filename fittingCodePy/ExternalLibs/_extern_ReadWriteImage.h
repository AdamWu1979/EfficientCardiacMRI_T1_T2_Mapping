/**
 * @file _extern_ReadWriteImage.h
 * @author Fani Deligianni
 * @date 09/04/2016
 * @brief IO interface to the nifti images. It uses the nifti library.
 *
 *  Created by Fani Deligianni on 09/04/2016.
 *  Copyright (c) 2016, InfoReach Consultancy. All rights reserved.
 *  See the LICENSE.txt file in the root folder
 *
 */

#ifndef _EXTERN_READWRITEIMAGE_H
#define _EXTERN_READWRITEIMAGE_H

#include "nifti1_io.h"
#include <string>

/** @defgroup NIFTYREG_FILEFORMAT_TYPE
 *  @brief Codes to define the image file format
 *  @{
 */
#define NR_NII_FORMAT 0
/* @} */

/* *************************************************************** */
/** The function checks the file format using the provided filename
  * Nifti is returned by default if no format are specified
  * @param filename Filename of the input images
  * @return Code, NIFTYREG_FILEFORMAT_TYPE,  that encode the file format
  */
int extern_io_checkFileFormat(const char *filename);
/* *************************************************************** */
/** The function expects a filename and returns a nifti_image structure
  * The function will use to correct library and will return a NULL image
  * if the image can not be read
  * @param filename Filename of the input images
  * @return Image as a nifti image
  */
nifti_image *extern_io_ReadImageFile(const char *filename);
/* *************************************************************** */
/** The function expects a filename and returns a nifti_image structure
  * The function will use to correct library and will return a NULL image
  * if the image can not be read
  * Only the header information is read and the actual data is not store
  * @param filename Filename of the input images
  * @return Image as a nifti image
  */
nifti_image *extern_io_ReadImageHeader(const char *filename);
/* *************************************************************** */
/** The function expects a filename and nifti_image structure
  * The image will be converted to the format specified in the
  * filename before being saved
  * @param image Nifti image to be saved
  * @param filename Filename of the output images
  */
void extern_io_WriteImageFile(nifti_image *image, const char *filename);
/* *************************************************************** */
#endif
