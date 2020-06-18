#include "Images.h"

#ifndef _IMAGES_CXX_
#define _IMAGES_CXX_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include "_curveFitting.h"

template <class T> 
Images<T>::Images(){
	
}

template<class T>
Images<T>::~Images(){

}

template<class T>
void Images<T>::fitImages() {
	
	optim = true;

	double *pntsT = new double[inputImages.size()];
	for (int i = 0; i < inputImages.size(); i++) {
		pntsT[i] = inversionTime[i];
	}

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

	vector<float> A((imgDims[1] * imgDims[2]),0.0);
	//A.reserve(imgDims[1]* imgDims[2]);
	vector<float> B((imgDims[1] * imgDims[2]), 0.0);
	vector<float> TI((imgDims[1] * imgDims[2]), 0.0);
	vector<float> dataT1star((imgDims[1] * imgDims[2]), 0.0);
	vector<float> dataXsquare((imgDims[1] * imgDims[2]), 0.0);
	vector<float> minIndexImg((imgDims[1] * imgDims[2]), 0.0);
	vector<float> statusImg((imgDims[1] * imgDims[2]), 0.0);
	vector<float> numIterImg((imgDims[1] * imgDims[2]), 0.0);


	std::ofstream out("logfile.txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf());				 //redirect std::cout to out.txt!

	std::time_t result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch\n" <<std::endl;

	double defaultVal = 1000.0;
	double pars[3] = { defaultVal, defaultVal, defaultVal };
	void(*func)(VectorXd&, Map<VectorXd>&, VectorXd&) = curveFitting::fexpcurve;
	int numpar = 3;
	int verbose0 = 0;
	long count = 0;
	double *pnts = new double[inputImages.size()];
	for (int ix = 0; ix<imgDims[1]; ix++) {
		for (int iy = 0; iy<imgDims[2]; iy++) {
			count = (long)(ix + iy*imgDims[1]);
			int countNonZeros = 0;
			for (int iz = 0; iz < inputImages.size(); iz++) {
				//int coord = ix+iy*dim[1]+dim[1]*dim[2]*iz;
				//pnts[iz] = (double)(dataPtr[coord] * niftiImages[0]->scl_slope + niftiImages[0]->scl_inter);
				pnts[iz] = (double) (inputImages[iz][ix + iy*imgDims[1]] );

				if (abs(pnts[iz]) > 0.5)
					countNonZeros++;
			}

			if (countNonZeros>numpar) {

				std::vector<double> pnts0(pnts, pnts + inputImages.size());

				std::vector<std::vector<double> > parAll;
				std::vector<double> Xsquares;
				int status0[5];

				double pars0[3] = { defaultVal, defaultVal, defaultVal };

				if (optim) {
					if (ix > 0 && ix < (imgDims[1] - 1) && iy>0 && iy < (imgDims[2] - 1)) {
				
						//if fitting of previous pnts was succesful
						int nx = ix - 1 + iy*imgDims[1];
						int ny = ix + (iy - 1)*imgDims[1];
						if (statusImg[nx] && statusImg[ny]) {
							//if parameters are close then average
							if ((abs(A[nx] - A[ny]) < 0.1* abs(A[nx] + A[ny])/2.0)
							&& (abs(B[nx] - B[ny]) < 0.1* abs(B[nx] + B[ny])/2.0)
							&& (abs(TI[nx] - TI[ny]) < 0.1* abs(TI[nx] + TI[ny])/2.0)) {

								pars0[0] = (A[nx] + A[ny]) / 2.0;
								pars0[1] = (B[nx] + B[ny]) / 2.0;
								pars0[2] = (TI[nx] + TI[ny]) / 2.0;
							}
						}
					}
					
				} //end of if(optim)

				int numIter[5];
				for (int j = 0; j<5; j++) {
					for (int i = 0; i<numpar; i++)
						pars[i] = pars0[i];

					for (int pp = 0; pp<inputImages.size(); pp++) {
						if (j>pp)
							pnts[pp] = -pnts0[pp];
					}

					double chiSquare;
					curveFitting myfitting(func, numpar, pars, inputImages.size(), pntsT, pnts, verbose0, &opts);
					myfitting.lm_fitting(pars, &chiSquare, &(status0[j]));
					myfitting.lm_info(numIter[j]);
					//std::cout << pars[0] << " " << pars[1] << " " << pars[2] << std::endl;

					std::vector<double> parTMP(pars, pars + numpar);
					parAll.push_back(parTMP);

					Xsquares.push_back(chiSquare);
				}

				std::vector<double>::iterator minXsquare = std::min_element(Xsquares.begin(), Xsquares.end());
				int minIndex = std::distance(Xsquares.begin(), minXsquare);
				double T1star = parAll[minIndex][2] * (parAll[minIndex][1] / parAll[minIndex][0] - 1.0);

				dataT1star[count] = (float)T1star;
				dataXsquare[count] = (float)(Xsquares[minIndex]);
				A[count] = (float)parAll[minIndex][0];;
				B[count] = (float)parAll[minIndex][1];;
				TI[count] = (float)parAll[minIndex][2];;
				minIndexImg[count] = (float)(minIndex);
				statusImg[count] = (float)(status0[minIndex]);
				numIterImg[count] = (float)(numIter[minIndex]);

				//fprintf(stdout, "fitting:%f %f %f %f %f %d %d\n", pars[0],pars[1],pars[2],T1star,xsquare, minIndex, status0[minIndex]);
				//std::cout << Xsquares[0]<< " "<< Xsquares[1] <<" "<< Xsquares[2] <<" " << Xsquares[3] << std::endl;

				//if (count % 1000 == 0)
				//	std::cout << "iter:" << ix << " " << iy << std::endl;

			}

		}

	}
	delete[] pntsT;
	delete[] pnts;

	result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch" <<std::endl;

	std::cout.rdbuf(coutbuf); //reset to standard output again


	outputImages.push_back(A);
	outputImages.push_back(B);
	outputImages.push_back(TI);
	outputImages.push_back(dataT1star);
	outputImages.push_back(dataXsquare);
	outputImages.push_back(minIndexImg);
	outputImages.push_back(statusImg);
	outputImages.push_back(numIterImg);

}

template<class T>
void Images<T>::fitImagesT2() {

	optim = false;

	double *pntsT = new double[inputImages.size()];
	for (int i = 0; i < inputImages.size(); i++) {
		pntsT[i] = inversionTime[i];
	}

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

	vector<float> A((imgDims[1] * imgDims[2]), 0.0);
	vector<float> T2((imgDims[1] * imgDims[2]), 0.0);
	vector<float> dataXsquare((imgDims[1] * imgDims[2]), 0.0);
	vector<float> statusImg((imgDims[1] * imgDims[2]), 0.0);
	vector<float> numIterImg((imgDims[1] * imgDims[2]), 0.0);

	std::time_t result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch\n" << std::endl;

	double defaultVal[2] = { 500.0, 40.0 };
	double pars[2] = { defaultVal[0], defaultVal[1] };
	void(*func)(VectorXd&, Map<VectorXd>&, VectorXd&) = curveFitting::fexpcurveT2;
	int numpar = 2;
	int verbose0 = 0;
	long count = 0;
	double *pnts = new double[inputImages.size()];
	for (int ix = 0; ix<imgDims[1]; ix++) {
		for (int iy = 0; iy<imgDims[2]; iy++) {
			count = (long)(ix + iy*imgDims[1]);
			int countNonZeros = 0;
			for (int iz = 0; iz < inputImages.size(); iz++) {
				//int coord = ix+iy*dim[1]+dim[1]*dim[2]*iz;
				//pnts[iz] = (double)(dataPtr[coord] * niftiImages[0]->scl_slope + niftiImages[0]->scl_inter);
				pnts[iz] = (double)(inputImages[iz][ix + iy*imgDims[1]]);

				if (abs(pnts[iz]) > 0.5)
					countNonZeros++;
			}

			if (countNonZeros>=numpar) {

				std::vector<double> pnts0(pnts, pnts + inputImages.size());

				std::vector<std::vector<double> > parAll;
				std::vector<double> Xsquares;
				int status0;
				int numIter;

				//double pars0[2] = { defaultVal[0], defaultVal[1]};
				for (int i = 0; i < numpar; i++)
					pars[i] = defaultVal[i];

				if (optim) {
					if (ix > 0 && ix < (imgDims[1] - 1) && iy>0 && iy < (imgDims[2] - 1)) {

						//if fitting of previous pnts was succesful
						int nx = ix - 1 + iy*imgDims[1];
						int ny = ix + (iy - 1)*imgDims[1];
						if (statusImg[nx] && statusImg[ny]) {
							//if parameters are close then average
							if ((abs(A[nx] - A[ny]) < 0.1* abs(A[nx] + A[ny]) / 2.0)
								&& (abs(T2[nx] - T2[ny]) < 0.1* abs(T2[nx] + T2[ny]) / 2.0)) {

								pars[0] = (A[nx] + A[ny]) / 2.0;
								pars[1] = (T2[nx] + T2[ny]) / 2.0;
							}
						}
					}

				} //end of if(optim)




				double chiSquare;
				curveFitting myfitting(func, numpar, pars, inputImages.size(), pntsT, pnts, verbose0, &opts);
				myfitting.lm_fitting(pars, &chiSquare, &status0);
				myfitting.lm_info(numIter);


				dataXsquare[count] = (float)(chiSquare);
				A[count] = (float)pars[0];
				T2[count] = (float)pars[1];
				statusImg[count] = (float)(status0);
				numIterImg[count] = (float)(numIter);

				//fprintf(stdout, "fitting:%f %f %f %f %f %d %d\n", pars[0],pars[1],pars[2],T1star,xsquare, minIndex, status0[minIndex]);
				//std::cout << Xsquares[0]<< " "<< Xsquares[1] <<" "<< Xsquares[2] <<" " << Xsquares[3] << std::endl;
				if (count % 1000 == 0)
					std::cout << "iter:" << ix << " " << iy << std::endl;

			}

		}

	}
	delete[] pntsT;
	delete[] pnts;

	result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch" << std::endl;

	outputImages.push_back(A);
	outputImages.push_back(T2);
	outputImages.push_back(dataXsquare);
	outputImages.push_back(statusImg);
	outputImages.push_back(numIterImg);

}

template<class T>
void Images<T>::fitImagesT2_3() {

	optim = false;

	double *pntsT = new double[inputImages.size()];
	for (int i = 0; i < inputImages.size(); i++) {
		pntsT[i] = inversionTime[i];
	}

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

	vector<float> A((imgDims[1] * imgDims[2]), 0.0);
	vector<float> B((imgDims[1] * imgDims[2]), 0.0);
	vector<float> T2((imgDims[1] * imgDims[2]), 0.0);
	vector<float> T2s((imgDims[1] * imgDims[2]), 0.0);
	vector<float> dataXsquare((imgDims[1] * imgDims[2]), 0.0);
	vector<float> statusImg((imgDims[1] * imgDims[2]), 0.0);
	vector<float> numIterImg((imgDims[1] * imgDims[2]), 0.0);

	std::time_t result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch\n" << std::endl;

	double defaultVal[3] = { 200.0, 400.0, 15.0 };
	double pars[3] = { defaultVal[0], defaultVal[1], defaultVal[2] };
	void(*func)(VectorXd&, Map<VectorXd>&, VectorXd&) = curveFitting::fexpcurveT2_3;
	int numpar = 3;
	int verbose0 = 0;
	long count = 0;
	double *pnts = new double[inputImages.size()];
	for (int ix = 0; ix<imgDims[1]; ix++) {
		for (int iy = 0; iy<imgDims[2]; iy++) {
			count = (long)(ix + iy*imgDims[1]);
			int countNonZeros = 0;
			for (int iz = 0; iz < inputImages.size(); iz++) {
				//int coord = ix+iy*dim[1]+dim[1]*dim[2]*iz;
				//pnts[iz] = (double)(dataPtr[coord] * niftiImages[0]->scl_slope + niftiImages[0]->scl_inter);
				pnts[iz] = (double)(inputImages[iz][ix + iy*imgDims[1]]);

				if (abs(pnts[iz]) > 0.5)
					countNonZeros++;
			}

			if (countNonZeros>=numpar) {

				std::vector<double> pnts0(pnts, pnts + inputImages.size());

				std::vector<std::vector<double> > parAll;
				std::vector<double> Xsquares;
				int status0;
				int numIter;

				//double pars0[3] = { defaultVal[0], defaultVal[1], defaultVal[2] };
				for (int i = 0; i < numpar; i++)
					pars[i] = defaultVal[i];


				if (optim) {
					if (ix > 0 && ix < (imgDims[1] - 1) && iy>0 && iy < (imgDims[2] - 1)) {

						//if fitting of previous pnts was succesful
						int nx = ix - 1 + iy*imgDims[1];
						int ny = ix + (iy - 1)*imgDims[1];
						if (statusImg[nx] && statusImg[ny]) {
							//if parameters are close then average
							if ((abs(A[nx] - A[ny]) < 0.1* abs(A[nx] + A[ny]) / 2.0)
								&& (abs(B[nx] - B[ny]) < 0.1* abs(B[nx] + B[ny]) / 2.0)
								&& (abs(T2[nx] - T2[ny]) < 0.1* abs(T2[nx] + T2[ny]) / 2.0)) {

								pars[0] = (A[nx] + A[ny]) / 2.0;
								pars[1] = (B[nx] + B[ny]) / 2.0;
								pars[2] = (T2[nx] + T2[ny]) / 2.0;
							}
						}
					}

				} //end of if(optim)

				double chiSquare;
				curveFitting myfitting(func, numpar, pars, inputImages.size(), pntsT, pnts, verbose0, &opts);
				myfitting.lm_fitting(pars, &chiSquare, &status0);
				myfitting.lm_info(numIter);


				dataXsquare[count] = (float)(chiSquare);
				A[count] = (float)pars[1];
				B[count] = (float)pars[0];
				T2[count] = (float)pars[2];
				T2s[count] = (float)(pars[2]*(pars[0]/pars[1]-1.0));
				statusImg[count] = (float)(status0);
				numIterImg[count] = (float)(numIter);

				//fprintf(stdout, "fitting:%f %f %f %f %f %d %d\n", pars[0],pars[1],pars[2],T1star,xsquare, minIndex, status0[minIndex]);
				//std::cout << Xsquares[0]<< " "<< Xsquares[1] <<" "<< Xsquares[2] <<" " << Xsquares[3] << std::endl;
				if (count % 1000 == 0)
					std::cout << "iter:" << ix << " " << iy << std::endl;

			}

		}

	}
	delete[] pntsT;
	delete[] pnts;

	result = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&result))
		<< result << " seconds since the Epoch" << std::endl;

	outputImages.push_back(A);
	outputImages.push_back(B);
	outputImages.push_back(T2);
	outputImages.push_back(T2s);
	outputImages.push_back(dataXsquare);
	outputImages.push_back(statusImg);
	outputImages.push_back(numIterImg);

}



BOOST_PYTHON_MODULE(pyImageDy)
{

	class_<Images<float>>("Images")
		.def("initImagePy", &Images<float>::initImagePy)
		.def("getImagePy", &Images<float>::getImagePy)
		.def("fitImagesPy", &Images<float>::fitImages)
		.def("fitImagesT2Py", &Images<float>::fitImagesT2)
		.def("fitImagesT2_3Py", &Images<float>::fitImagesT2_3)
		;
}



#endif /* _IMAGES_CXX_ */