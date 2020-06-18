/** \file fitNiftiImages.h
    \A class to fit each voxel from a set of nifti images.
           - Written by Fani Deligianni, InfoReach Consultancy
 */
#ifndef _FITNIFTIIMAGES_HEADER_
#define _FITNIFTIIMAGES_HEADER_

#include "_extern_ReadWriteImage.h"
#include <vector>
using namespace std;

template <class T> class fitNiftiImages{
	
protected:
	std::vector<nifti_image *> niftiImages;
	nifti_image *inputNiftiImage;
	std::vector<T> currentMins;
	std::vector<T> currentMaxs;
    int dim[8];                  /*!< dim[0]=ndim, dim[1]=nx, dim[2]=ny, dim[3]=nz, etc         */
    size_t nvox;                    /*!< number of voxels = nx*ny*nz*...*nw   */
    int nbyper;                  /*!< bytes per voxel, matches datatype    */
    int datatype;                /*!< type of data in voxels: DT_* code    */
	
	vector<nifti_image*> outputNiftis;
	std::vector<double> *inversionTime;
	
	bool verbose;	
public:
	fitNiftiImages();
	~fitNiftiImages();
	void initNiftiImages(nifti_image *image,std::vector<double> *invTime);
	void allocateOutputNifti();
	void clearOutputImages();
	bool checkCompatibility(nifti_image *image);
	void fitImages(char *ouputBaseName);
	bool testingFunction(char *outputFilename, int imgId);
	
};

#include "fitNiftiImages.cxx"
#endif /* _FITNIFTIIMAGES_HEADER_ */