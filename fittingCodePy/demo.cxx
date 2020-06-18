// A simple program that computes the square root of a number
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "_extern_ReadWriteImage.h"
#include "fitNiftiImages.h"
#include "demoConfig.h"
using namespace std;

void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("This is a demo function.\n");
    printf("Use --h for help and -v for version information.\n");
	
}

 
int main (int argc, char *argv[])
{
    std::vector<char *> referenceImagesName;
	std::vector<nifti_image *> niftiImages;
	std::vector<double> inverseTime;
	double inputValue;
	bool inputValueFlag = false;
	char *ouputPPMFile = NULL;
	bool ouputPPMFileFlag = false;
	char *ouputNiftiImage = NULL;
	
    /* read the input parameter */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }
		else if(strcmp(argv[i], "-version")==0 ||
	            strcmp(argv[i], "-Version")==0 ||
	            strcmp(argv[i], "-V")==0 ||
	            strcmp(argv[i], "-v")==0 ||
	            strcmp(argv[i], "--v")==0 ||
	            strcmp(argv[i], "--version")==0){
				    fprintf(stdout,"%s Version %d.%d\n",
				            argv[0],
				            demo_VERSION_MAJOR,
				            demo_VERSION_MINOR);
				    fprintf(stdout,"Usage: %s number\n",argv[0]);
		            return 0;
			}
			else if(strcmp(argv[i], "-input")==0){
				inputValue = atof(argv[++i]);
				inputValueFlag = true;
			}
			else if(strcmp(argv[i], "-writePPM")==0){
				ouputPPMFile = argv[++i];
				ouputPPMFileFlag = true;
			}
			else if(strcmp(argv[i], "-ouputNameImage")==0){
				ouputNiftiImage = argv[++i];
			}					
			else if(strcmp(argv[i], "-inputImage")==0){
				referenceImagesName.push_back(argv[++i]);
			}
			else{
				inverseTime.push_back(atof(argv[i]));
			}
	} //end of for
	
    /* Read the images and check their dimensions */	
	for(int i=0; i<referenceImagesName.size(); i++){
	    nifti_image *referenceHeader = extern_io_ReadImageFile(referenceImagesName[i]);
	    if(referenceHeader == NULL){
	        fprintf(stderr,"* ERROR Error when reading the reference  image: %s\n",referenceImagesName[i]);
	        return 1;
	    }else{
			niftiImages.push_back(referenceHeader);
	    	fprintf(stdout,"ndim:%d nx:%d ny:%d nz:%d nt:%d\n",niftiImages[i]->ndim,referenceHeader->nx,referenceHeader->ny,referenceHeader->nz,referenceHeader->nt);
			fprintf(stdout,"SUCCESS\n");
	    }
	}
	
	
	if(niftiImages.size()>0){
		fprintf(stdout,"SUCCESS2\n");
		fprintf(stdout,"%d\n",niftiImages[0]->datatype);
		
		if(ouputNiftiImage==NULL){
	        fprintf(stderr,"MAIN ERROR: Error output nifti image name is not specified\n");
	        return 1;
		}
		
        switch(niftiImages[0]->datatype){
          case NIFTI_TYPE_UINT8:{
  			fitNiftiImages<unsigned char> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);					
  			}
			
              break;        	
          }
          case NIFTI_TYPE_INT8:{
  			fitNiftiImages<char> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
          }
          case NIFTI_TYPE_UINT16:{
  			fitNiftiImages<unsigned short> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;        	
          }
          case NIFTI_TYPE_INT16:{
  			fitNiftiImages<short> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
          }
          case NIFTI_TYPE_UINT32:{
  			fitNiftiImages<unsigned int> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
          }
          case NIFTI_TYPE_INT32:{
  			fitNiftiImages<int> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
          }
          case NIFTI_TYPE_FLOAT32:{
  			fitNiftiImages<float> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
          }
          case NIFTI_TYPE_FLOAT64:{
  			fitNiftiImages<double> cardiacMRI;
  			for(int i=0; i<niftiImages.size(); i++){
  				cardiacMRI.initNiftiImages(niftiImages[i], &inverseTime);
  				if(ouputPPMFileFlag)
  					cardiacMRI.testingFunction(ouputPPMFile, i);
  				cardiacMRI.fitImages(ouputNiftiImage);
  			}
              break;
        	
          }
        }
		
  	}
	
    if(inputValueFlag){
	    double outputValue = sqrt(inputValue);
	    fprintf(stdout,"The square root of %g is %g\n",
	            inputValue, outputValue);    	
    }
    return 0;
}

