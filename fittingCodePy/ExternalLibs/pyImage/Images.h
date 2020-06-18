/** \file Images.h
    \A class that encodes images.
           - Written by Fani Deligianni
 */
#ifndef _IMAGES_HEADER_
#define _IMAGES_HEADER_

#define BOOST_PYTHON_STATIC_LIB

#include <boost/python.hpp>
#include <vector>

using namespace boost::python;
using namespace std;

template <class T> class Images{
	
protected:
	std::vector<vector<T>> inputImages;
	std::vector<T> currentMins;
	std::vector<T> currentMaxs;
    std::vector<int> imgDims; 

	std::vector<vector<float>> outputImages;
	std::vector<double> inversionTime;

	bool verbose;	
	bool optim;
public:
	Images();
	~Images();

	void fitImages();
	void fitImagesT2();
	void fitImagesT2_3();

	void initImagePy(boost::python::list& img, boost::python::list& dims, float invTime)
	{
		for (int i = 0; i < len(dims); ++i)
		{
			if(inputImages.size()<1)
				imgDims.push_back(boost::python::extract<int>(dims[i]));
			else if (boost::python::extract<int>(dims[i]) != imgDims[i]) {
				std::cout << "Image dimension should match original image." << std::endl;
				return;
			}

		}

		std::vector<T> tmpimg;
		for (int i = 0; i < len(img); ++i)
		{
			tmpimg.push_back(boost::python::extract<T>(img[i]));
		}
		inputImages.push_back(tmpimg);


		inversionTime.push_back(invTime);
	}

	boost::python::list getImagePy(int index)
	{
		boost::python::list pyList;

		if (outputImages.size() < 1) {
			pyList.append(0);
			return pyList;
		}

		for (unsigned int i = 0; i < outputImages[index].size(); i++)
			pyList.append( (outputImages[index])[i] );

		return pyList;
	}


};


#include "Images.cxx"
#endif /* _IMAGES_HEADER_ */