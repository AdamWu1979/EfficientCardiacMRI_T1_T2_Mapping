/**
 * @file _extern_curveFitting.h
 * @author Fani Deligianni
 * @date 03/05/2016
 * @Curve fitting library - .
 *
 *  Created by Fani Deligianni on 03/05/2016.
 *  Copyright (c) 2016, InfoReach Consultancy. All rights reserved.
 *
 * *************************************************************** 
 ** This class performs curve fitting based on the Levenberg-Marquardt algorithm
 */

#ifndef _CURVEFITTING_H
#define _CURVEFITTING_H

#include <Eigen/Dense>

using namespace Eigen;

typedef struct {
    int prnt;          		// >1 intermediate results; >2 plots
    int maxIter;            // maximum number of iterations
    double epsilon_1;     	// convergence tolerance for gradient
    double epsilon_2;     	// convergence tolerance for parameters
    double epsilon_3;     	// convergence tolerance for Chi-square
    double epsilon_4;     	// determines acceptance of a L-M step
    double lambda_0;      	// initial value of damping paramter, lambda
    double lambda_UP_fac; 	// factor for increasing lambda
    double lambda_DN_fac; 	// factor for decreasing lambda
    int Update_Type;   		// 1: Levenberg-Marquardt lambda update
} lm_opts;


class curveFitting{
protected:
	bool init_flag;
	int verbose;
	int iteration;
		
	int npar;
	int npnts;
	lm_opts *opts;
	
	VectorXd par;

	Map<VectorXd> t;
	Map<VectorXd> y;

	MatrixXd J;
	MatrixXd JtWJ; //(npar,npar);
	VectorXd JtWdy; //(npar);
	
	VectorXd weight_sq;
	double lambda;
	VectorXd dp; 
		
	void (*func)(VectorXd&, Map<VectorXd>&, VectorXd&);
	
	void initParameters();
	void lm_matx(VectorXd& p_old, VectorXd& y_old, double dx2, double& x2, VectorXd& y_hat);
	void lm_FD_J(VectorXd& p, VectorXd& y_hat); 
	void lm_Broyden_J(VectorXd& p_old, VectorXd& y_old, VectorXd& p, VectorXd& y_hat);

public:
	curveFitting(void(*pfn)(VectorXd& , Map<VectorXd>&,VectorXd&), int numpar, double *p, int numpnts, double *pntsT, double *pnts, int verbose0, lm_opts *opts0);
	~curveFitting();

	int lm_fitting(double *pr, double *chiSq, int *status0);
	void lm_info(int &numIter);
	
	static void fexpcurve(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat);
	static void fexpcurveT2(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat);
	static void fexpcurveT2_3(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat);
	
};


#endif
