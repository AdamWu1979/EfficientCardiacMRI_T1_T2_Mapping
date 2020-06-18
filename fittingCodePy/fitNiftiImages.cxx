#include "fitNiftiImages.h"

#ifndef _FITNIFTIIMAGES_CXX_
#define _FITNIFTIIMAGES_CXX_

#include <iostream>
#include <ctime>
#include <cmath>
#include "_curveFitting.h"

template <class T> 
fitNiftiImages<T>::fitNiftiImages(){
	nvox = 0;                    /*!< number of voxels = nx*ny*nz*...*nw   */
	nbyper = 0;                  /*!< bytes per voxel, matches datatype    */
	datatype = 0;                /*!< type of data in voxels: DT_* code    */
	
	for(int i=0; i<8; i++)
		dim[i] = 0;
	
	verbose = true;
}

template<class T>
fitNiftiImages<T>::~fitNiftiImages(){
	clearOutputImages();
}

template <class T> 
bool fitNiftiImages<T>::checkCompatibility(nifti_image *image){
	bool flagComp = true;
	for(int i=0; i<8; i++){
		if(dim[i]!=image->dim[i] && dim[i]!=0){
			flagComp = false;
			break;
		}
	}
	
	if(nvox!=0 && nvox!=image->nvox)
		flagComp = false;

	if(nbyper!=0 && nbyper!=image->nbyper)
		flagComp = false;

	if(datatype!=0 && datatype!=image->datatype)
		flagComp = false;
	
	if(flagComp && nvox==0){
		for(int i=0; i<8; i++)
			dim[i] = image->dim[i];
		nvox = image->nvox;
		nbyper = image->nbyper;
		datatype = image->datatype;
	}
	
	return flagComp;
}

template <class T> 
void fitNiftiImages<T>::initNiftiImages(nifti_image *image, std::vector<double> *invTime){
	if(!checkCompatibility(image)){
        fprintf(stderr,"[fitNiftiImages ERROR] fitNiftiImages::initNiftiImages()\n");
        fprintf(stderr,"[fitNiftiImages ERROR] Compatibility check failed\n");
		exit(1);		
	}
	
	if( invTime->size() != image->dim[3]){
        fprintf(stderr,"[fitNiftiImages ERROR] fitNiftiImages::initNiftiImages()\n");
        fprintf(stderr,"[fitNiftiImages ERROR] InversionTime dimensions do not match image z dimension\n");
		exit(1);		
	}
	
	inversionTime = invTime;
		
    T *imagePtr = static_cast<T *>(image->data);
    unsigned int voxelNumber = image->nx*image->ny*image->nz;
	
    T currentMin=0;
    T currentMax=0;
    switch(image->datatype){
    case NIFTI_TYPE_UINT8:
        currentMin=(T)std::numeric_limits<unsigned char>::max();
        currentMax=0;
        break;
    case NIFTI_TYPE_INT8:
        currentMin=(T)std::numeric_limits<char>::max();
        currentMax=(T)-std::numeric_limits<char>::max();
        break;
    case NIFTI_TYPE_UINT16:
        currentMin=(T)std::numeric_limits<unsigned short>::max();
        currentMax=0;
        break;
    case NIFTI_TYPE_INT16:
        currentMin=(T)std::numeric_limits<char>::max();
        currentMax=-(T)std::numeric_limits<char>::max();
        break;
    case NIFTI_TYPE_UINT32:
        currentMin=(T)std::numeric_limits<unsigned int>::max();
        currentMax=0;
        break;
    case NIFTI_TYPE_INT32:
        currentMin=(T)std::numeric_limits<int>::max();
        currentMax=-(T)std::numeric_limits<int>::max();
        break;
    case NIFTI_TYPE_FLOAT32:
        currentMin=(T)std::numeric_limits<float>::max();
        currentMax=-(T)std::numeric_limits<float>::max();
        break;
    case NIFTI_TYPE_FLOAT64:
        currentMin=(T)std::numeric_limits<double>::max();
        currentMax=-(T)std::numeric_limits<double>::max();
        break;
    }

    for(int t=0;t<image->nt;t++){
        T *volumePtr = &imagePtr[t*voxelNumber];
        // Extract the minimal and maximal values from the current volume
        if(image->scl_slope==0) image->scl_slope=1.0f;
        for(unsigned int index=0; index<voxelNumber; index++){
            T value = (T)(*volumePtr++ * image->scl_slope + image->scl_inter);
			currentMin=(currentMin<value)?currentMin:value;
            currentMax=(currentMax>value)?currentMax:value;
        }
	}
		currentMins.push_back(currentMin);
		currentMaxs.push_back(currentMax);
		niftiImages.push_back(image);
	
	
}

template<class T>
void fitNiftiImages<T>::allocateOutputNifti(){
    if(niftiImages.size()<1){
        fprintf(stderr,"[fitNiftiImages ERROR] fitNiftiImages::allocateNiftiOutputImage()\n");
        fprintf(stderr,"[fitNiftiImages ERROR] Reference and floating images are not defined. Exit.\n");
        exit(1);
    }
	nifti_image *outputImage = NULL;
    outputImage = nifti_copy_nim_info(niftiImages[0]);
    outputImage->dim[0]=outputImage->ndim=2;
    outputImage->dim[4]=outputImage->nt=niftiImages[0]->nt;
    outputImage->pixdim[4]=outputImage->dt=1.0;
    outputImage->nvox = outputImage->nx * outputImage->ny;
    outputImage->datatype = NIFTI_TYPE_FLOAT32;
    outputImage->nbyper = sizeof(float);
    outputImage->data = (void *)calloc(outputImage->nvox, outputImage->nbyper);
	
	outputNiftis.push_back(outputImage);
}

template<class T>
void fitNiftiImages<T>::clearOutputImages(){
	for(int i=0; i<outputNiftis.size(); i++){
		nifti_image_free(outputNiftis[i]);
	}
	outputNiftis.clear();
	
}

template<class T> 
void fitNiftiImages<T>::fitImages(char *outputBaseName){
	T *dataPtr=static_cast<T *>(niftiImages[0]->data);
	double *pntsT = new double[dim[3]];
	int count = 0;
	for(int i=0; i<dim[3]; i++){
		pntsT[i] = inversionTime[0][i];
		count++;
	}
	
	//sort out time points
	int *indexT = new int[dim[3]];
	for(int i=0; i<dim[3]; i++)
		indexT[i] = i;
	
	for(int i=0; i<dim[3]; i++){
		double small = pntsT[i];
		int index0 = i;
		for(int j=i+1;j<dim[3];j++){
			if(pntsT[j]<small){
				small = pntsT[j];
				index0 = j;
			}
		}
		pntsT[index0] = pntsT[i];
		pntsT[i] = small;
		int tmpInd = indexT[i];
		indexT[i] = indexT[index0];
		indexT[index0] = tmpInd;
	}
	
	allocateOutputNifti();
	float *dataT1star = static_cast<float *>(outputNiftis[0]->data);
	allocateOutputNifti();
	float *dataXsquare = static_cast<float *>(outputNiftis[1]->data);

	allocateOutputNifti();
	float *A = static_cast<float *>(outputNiftis[2]->data);
	allocateOutputNifti();
	float *B = static_cast<float *>(outputNiftis[3]->data);
	allocateOutputNifti();
	float *TI = static_cast<float *>(outputNiftis[4]->data);

	allocateOutputNifti();
	float *statusImg = static_cast<float *>(outputNiftis[5]->data);

	allocateOutputNifti();
	float *numIterImg = static_cast<float *>(outputNiftis[6]->data);

	allocateOutputNifti();
	float *minIndexImg = static_cast<float *>(outputNiftis[7]->data);

	
	double xsquare;
	int minindex;

	lm_opts opts;
	opts.prnt = 0;
	opts.maxIter = 150;
	opts.epsilon_1 = 1e-5;
	opts.epsilon_2 = 1e-5;
	opts.epsilon_3 = 1e-5;
	opts.epsilon_4 = 1e-1;
	opts.lambda_0 = 1e-2;
	opts.lambda_UP_fac = 11;
	opts.lambda_DN_fac = 9;
	opts.Update_Type = 1;

	
    std::time_t result = std::time(nullptr);
    std::cout << std::asctime(std::localtime(&result))
              << result << " seconds since the Epoch\n";

	double defaultVal = 1000.0;
	double pars[3] = { defaultVal, defaultVal, defaultVal };
	void(*func)(VectorXd&, Map<VectorXd>&, VectorXd&) = curveFitting::fexpcurve;
	int numpar = 3;
	int status;
	int verbose0 = 0;
	double *pnts = new double[dim[3]];
	for(int ix=0; ix<dim[1]; ix++){
		for(int iy=0; iy<dim[2]; iy++){
			int countNonZeros = 0;
			for(int iz=0; iz<dim[3]; iz++){
				//int coord = ix+iy*dim[1]+dim[1]*dim[2]*iz;
				//pnts[iz] = (double)(dataPtr[coord] * niftiImages[0]->scl_slope + niftiImages[0]->scl_inter);
				pnts[iz] = (double)(dataPtr[ix + iy*dim[1] + dim[1] * dim[2] * iz]);

				if(abs(pnts[iz]) > 0.5)
					countNonZeros++;
			}
		
			if(countNonZeros>numpar){
				double *pntsS = new double[dim[3]];
				for(int ii=0; ii<dim[3]; ii++){
					pntsS[ii] = pnts[indexT[ii]];
				}
				double *tmp = pnts;
				pnts = pntsS;
				delete [] tmp;
				
				//double pars[3]={6.0978, 6.0978, 6.0978};

				//pnts[0]=203.0;pnts[1]=149.0;pnts[2]=557.0;pnts[3]=541.0;pnts[4]=525.0;pnts[5]=494.0;pnts[6]=504.0;pnts[7]=482.0;
				//pnts[0]=-203.0;pnts[1]=-149.0;pnts[2]=557.0;pnts[3]=541.0;pnts[4]=525.0;pnts[5]=494.0;pnts[6]=504.0;pnts[7]=482.0;
				//pntsT[0]=129.0;pntsT[1]=209.0;pntsT[2]=857.0;pntsT[3]=1464.0;pntsT[4]=1504.0;pntsT[5]=2127.0;pntsT[6]=2792.0;pntsT[7]=3449.0;
				//extern_expCurveFitting(pnts, pntsT, dim[3], pars, &xsquare, verbose, &minindex, &status);
				std::vector<double> pnts0(pnts, pnts + dim[3]);

				std::vector<std::vector<double> > parAll;
				std::vector<double> Xsquares;
				int status0[5];
				
				double pars0[3] = { defaultVal, defaultVal, defaultVal };
/*				if (ix > 0 && ix < (dim[1] - 1) && iy>0 && iy < (dim[2] - 1)) {
					//if fitting of previous pnts was succesful
					if (statusImg[ix-1 + iy*dim[1]] && statusImg[ix + (iy-1)*dim[1]]) { 
						//if parameters are close then average
						if ((abs(A[ix - 1 + iy*dim[1]] - A[ix + (iy - 1)*dim[1]]) < 0.1* abs(A[ix - 1 + iy*dim[1]] + A[ix + (iy - 1)*dim[1]])/2.0)
							&& (abs(B[ix - 1 + iy*dim[1]] - B[ix + (iy - 1)*dim[1]]) < 0.1* abs(B[ix - 1 + iy*dim[1]] + B[ix + (iy - 1)*dim[1]])/2.0)
							&& (abs(TI[ix - 1 + iy*dim[1]] - TI[ix + (iy - 1)*dim[1]]) < 0.1* abs(TI[ix - 1 + iy*dim[1]] + TI[ix + (iy - 1)*dim[1]])/2.0)) {

							pars0[0] = (A[ix - 1 + iy*dim[1]] + A[ix + (iy - 1)*dim[1]]) / 2.0;
							pars0[1] = (B[ix - 1 + iy*dim[1]] + B[ix + (iy - 1)*dim[1]]) / 2.0;
							pars0[2] = (TI[ix - 1 + iy*dim[1]] + TI[ix + (iy - 1)*dim[1]]) / 2.0;
						}
					}
				}
*/
				int numIter[5];
				for (int j = 0; j<5; j++) {
					for (int i = 0; i<numpar; i++)
						pars[i] = pars0[i];

					for (int pp = 0; pp<dim[3]; pp++) {
						if (j>pp)
							pnts[pp] = -pnts0[pp];
					}

					double chiSquare;
					curveFitting myfitting(func, numpar, pars, dim[3], pntsT, pnts, verbose0,&opts);
					myfitting.lm_fitting(pars, &chiSquare, &(status0[j]) );
					myfitting.lm_info(numIter[j]);
					//std::cout << pars[0] << " " << pars[1] << " " << pars[2] << std::endl;

					std::vector<double> parTMP(pars,pars+numpar);
					parAll.push_back(parTMP);

					Xsquares.push_back(chiSquare);
				}

				std::vector<double>::iterator minXsquare = std::min_element(Xsquares.begin(), Xsquares.end());
				int minIndex = std::distance(Xsquares.begin(), minXsquare);
				double T1star = parAll[minIndex][2]*(parAll[minIndex][1]/ parAll[minIndex][0]-1.0);

				/*if (status0[minIndex] == 0) {
					for (int i = 0; i<numpar; i++)
						pars[i] = defaultVal;
				}
				else {
					for (int i = 0; i<numpar; i++)
						pars[i] = parAll[minIndex][i];
				}
				*/
					

				dataT1star[ix+iy*dim[1]] = (float)T1star;
				dataXsquare[ix+iy*dim[1]] = (float)(Xsquares[minIndex]);
				A[ix+iy*dim[1]] = (float)parAll[minIndex][0];;
				B[ix+iy*dim[1]] = (float)parAll[minIndex][1];;
				TI[ix+iy*dim[1]] = (float)parAll[minIndex][2];;
				minIndexImg[ix+iy*dim[1]] = (float)(minIndex);
				statusImg[ix + iy*dim[1]] = (float)(status0[minIndex]);
				numIterImg[ix + iy*dim[1]] = (float)(numIter[minIndex]);
				//fprintf(stdout, "fitting:%f %f %f %f %f %d %d\n", pars[0],pars[1],pars[2],T1star,xsquare, minIndex, status0[minIndex]);
				//std::cout << Xsquares[0]<< " "<< Xsquares[1] <<" "<< Xsquares[2] <<" " << Xsquares[3] << std::endl;
				}
				
		}
		
	}
	delete [] pntsT;
	delete [] indexT;
	delete [] pnts;

	result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch\n";
	
	int newSize = strlen(outputBaseName)  + strlen("_T1star.nii.gz") + 1; 
	char * outputFileName = new char[newSize];
    strcpy(outputFileName,outputBaseName);
    strcat(outputFileName,"_T1star.nii.gz"); 
	extern_io_WriteImageFile(outputNiftis[0],outputFileName);
	delete [] outputFileName;
	
	newSize = strlen(outputBaseName)  + strlen("_Xsquare.nii.gz") + 1; 
	outputFileName = new char[newSize];
    strcpy(outputFileName,outputBaseName);
    strcat(outputFileName,"_Xsquare.nii.gz"); 
	extern_io_WriteImageFile(outputNiftis[1],outputFileName);
	delete [] outputFileName;

	newSize = strlen(outputBaseName)  + strlen("_A.nii.gz") + 1; 
	outputFileName = new char[newSize];
    strcpy(outputFileName,outputBaseName);
    strcat(outputFileName,"_A.nii.gz"); 
	extern_io_WriteImageFile(outputNiftis[2],outputFileName);
	delete [] outputFileName;
	
	newSize = strlen(outputBaseName)  + strlen("_B.nii.gz") + 1; 
	outputFileName = new char[newSize];
    strcpy(outputFileName,outputBaseName);
    strcat(outputFileName,"_B.nii.gz"); 
	extern_io_WriteImageFile(outputNiftis[3],outputFileName);
	delete [] outputFileName;
	
	newSize = strlen(outputBaseName)  + strlen("_TI.nii.gz") + 1; 
	outputFileName = new char[newSize];
    strcpy(outputFileName,outputBaseName);
    strcat(outputFileName,"_TI.nii.gz"); 
	extern_io_WriteImageFile(outputNiftis[4],outputFileName);
	delete [] outputFileName;

	newSize = strlen(outputBaseName) + strlen("_status.nii.gz") + 1;
	outputFileName = new char[newSize];
	strcpy(outputFileName, outputBaseName);
	strcat(outputFileName, "_status.nii.gz");
	extern_io_WriteImageFile(outputNiftis[5], outputFileName);
	delete[] outputFileName;

	newSize = strlen(outputBaseName) + strlen("_numIter.nii.gz") + 1;
	outputFileName = new char[newSize];
	strcpy(outputFileName, outputBaseName);
	strcat(outputFileName, "_numIter.nii.gz");
	extern_io_WriteImageFile(outputNiftis[6], outputFileName);
	delete[] outputFileName;

	newSize = strlen(outputBaseName)  + strlen("_minIndex.nii.gz") + 1;
	outputFileName = new char[newSize];
	strcpy(outputFileName,outputBaseName);
	strcat(outputFileName,"_minIndex.nii.gz");
	extern_io_WriteImageFile(outputNiftis[7],outputFileName);
	delete [] outputFileName;
	



}


template <class T> 
bool fitNiftiImages<T>::testingFunction(char *outputFilename, int imgId){
	if(niftiImages.size()<=imgId)
		return false;

	
	int width = dim[1];
	int height = dim[2];
	T *dataPtr=static_cast<T *>(niftiImages[imgId]->data);
	
	int cMin = (int)currentMins[imgId];
	int cMax = (int)currentMaxs[imgId];

	int maxImgInt = 255;
	FILE	*pfile = fopen( outputFilename, "wb" );
	if( pfile == NULL )	//display error - cannot open file for reading
		return false;
	
	//write header
	fprintf( pfile, "P6\n" );
	//write width - height - maximum color value
	fprintf( pfile, "%d %d %d\n", width, height, maxImgInt );
	fprintf(stdout, "slope,inter:%f %f\n", niftiImages[imgId]->scl_slope,niftiImages[imgId]->scl_inter);
	for( int i=0; i<width*height; i++ )
	{
		
		T value = (T)(dataPtr[i] * niftiImages[imgId]->scl_slope + niftiImages[imgId]->scl_inter);
		float valueNew = round((float)maxImgInt*((float)value-(float)cMin)/((float)cMax-(float)cMin));
		fprintf(stdout,"%f,%d,%d,%d\n",valueNew, (short)value,cMin,cMax);

		unsigned char val = (unsigned char)valueNew;
		fwrite( &val, sizeof(unsigned char), 1, pfile );
		fwrite( &val, sizeof(unsigned char), 1, pfile );
		fwrite( &val, sizeof(unsigned char), 1, pfile );
	}

	fclose( pfile );
		
	return true;	
}

#endif /* _FITNIFTIIMAGES_CXX_ */