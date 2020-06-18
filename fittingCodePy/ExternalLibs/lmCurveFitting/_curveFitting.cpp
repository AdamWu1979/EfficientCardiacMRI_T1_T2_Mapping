/*
 *  _curveFitting.cpp
 *
 *
 *  Created by Fani Deligianni on 15/04/2016.
 *  Copyright (c) 2016, InfoReach Consultancy. All rights reserved.
 *  InfoReach Consultancy
 *
 */

#ifndef _CURVEFITTING_CPP
#define _CURVEFITTING_CPP

#include <iostream>
#include <limits>
#include "_curveFitting.h"

using namespace Eigen;
using namespace std;


curveFitting::curveFitting(void(*pfn)(VectorXd&, Map<VectorXd>&,VectorXd&), int numpar, double *p, int numpnts, double *pntsT, double *pnts, int verbose0, lm_opts *opts0):
	t(pntsT,numpnts),
	y(pnts,numpnts)
{
	try{

		init_flag = false;

		npar = numpar;
		npnts = numpnts;
		
		verbose = verbose0;
		
		J = MatrixXd::Zero(npnts,npar);
		JtWJ = MatrixXd::Zero(npar,npar);
		JtWdy = VectorXd::Zero(npar);
		weight_sq = VectorXd::Zero(npnts);
		func 		= pfn;
		
		par = Map<VectorXd>(p,npar);
		
		double tmp = npnts/ std::sqrt(y.transpose()*(y));
		weight_sq = (tmp*VectorXd::Ones(npnts)).array().square();
		
		dp = -0.01*(VectorXd::Ones(npar));
		
		opts = opts0;
				
		if(verbose){
			cout<< "t:" << t.transpose() <<endl;
			cout<< "y:" << y.transpose() <<endl;
			cout<< "p:" << par.transpose() <<endl;
			cout << "curveFitting::init SUCCESS"<<endl;
		}
		
	}catch(const char* msg){
		std::cout << msg << std::endl;	
		std::cout << "curveFitting::init FAILED" << std::endl;	
	}
	
}

curveFitting::~curveFitting(){
	try{
		iteration = 0;

		if (init_flag)
			delete opts;
		
		if(verbose)
			cout << "curveFitting::~curveFitting -success"<<endl;
		
	}catch(const char* msg){
		std::cout << msg << std::endl;		
	}
	
}

void curveFitting::initParameters(){
	try{
		opts = new lm_opts;
		opts->prnt = 0;
		opts->maxIter = 150;
		opts->epsilon_1 = 1e-5;
		opts->epsilon_2 = 1e-5;
		opts->epsilon_3 = 1e-5;
		opts->epsilon_4 = 1e-1;
		opts->lambda_0 =  1e-2;
		opts->lambda_UP_fac = 11;
		opts->lambda_DN_fac = 9;
		opts->Update_Type = 1;
		init_flag = true;

		if(verbose)
			cout << "curveFitting::initParameters -success"<<endl;
		
	}catch(const char* msg){
		std::cout << msg << std::endl;		
	}
 
}

void curveFitting::lm_info(int &numIter) {
	numIter = iteration;
}


void curveFitting::fexpcurve(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat){
	y_hat = p(0) - p(1)*(-t_hat/(p(2))).array().exp(); 
}

void curveFitting::fexpcurveT2(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat) {
	y_hat = p(0)*(-t_hat / (p(1))).array().exp();
}

void curveFitting::fexpcurveT2_3(VectorXd& p, Map<VectorXd>& t_hat, VectorXd& y_hat) {
	y_hat = p(0) + p(1)*(-t_hat / (p(2))).array().exp();
}


/*
* lm_FD_J(VectorXd& p0, VectorXd& y_hat)
*
* partial derivatives (Jacobian) dy/dp for use with lm.m
* computed via Finite Differences
*/
void curveFitting::lm_FD_J(VectorXd& p0, VectorXd& y_hat){

	try{
		
		VectorXd del = VectorXd::Zero(npar);
		VectorXd p = p0;
		VectorXd ps = p0;
		
		for(int j=0; j<npar; j++){
			del(j) = dp(j)*(1.0+abs( p(j)) );
			p(j) = ps(j) + del(j); 
		
			if( std::abs(del(j)) > 0.0000001 ){
				VectorXd y1(npnts); 
				func(p, t, y1);

				if(dp(j)<0){  //backwards difference
					J.col(j) = (y1-y_hat) / del(j);			
				}else{
					p(j) = ps(j) - del(j);
					VectorXd y2(npnts); 
					func(p, t, y2);

					J.col(j) = (y1-y2) / (2.0*del(j));										
				}
			}
			
			p(j) = ps(j);
		
		} //end of for
		
	}catch(const char* msg){
		std::cout << msg << std::endl;		
	}
	
} 

/*
* lm_Broyden_J(VectorXd& p_old, VectorXd& y_old, VectorXd& p, VectorXd& y_hat)
*
* carry out a rank-1 update to the Jacobian matrix using Broyden's equation
*/
void curveFitting::lm_Broyden_J(VectorXd& p_old, VectorXd& y_old, VectorXd& p, VectorXd& y_hat){
	VectorXd h = p-p_old;
	J = J + (y_hat-y_old - J*h) * h.transpose() / (h.transpose()*h);
}

/*
* lm_matxlm_matx(VectorXd& p_old, VectorXd& y_old, double dx2, double& chi_sq, VectorXd& y_hat)
*
* Evaluate the linearized fitting matrix, JtWJ, and vector JtWdy, 
* and calculate the Chi-squared error function, Chi_sq 
* Used by Levenberg-Marquard algorithm  
*/
void curveFitting::lm_matx(VectorXd& p_old, VectorXd& y_old, double dx2, double& chi_sq, VectorXd& y_hat){
	try{
		func(par, t, y_hat);

		VectorXd p = par;
		if ( (iteration%(2*npar))==0 || dx2 > 0 ){
		    lm_FD_J(p,y_hat);		// finite difference
	    } else{
	    	lm_Broyden_J(p_old,y_old,p,y_hat);
	    } 
				
		VectorXd delta_y(npnts);
	    delta_y = y - y_hat;	// residual error between model and data

	    chi_sq = (delta_y.transpose())*( delta_y.array() * weight_sq.array() ).matrix(); 	// Chi-squared error criteria
		
		VectorXd onesV = VectorXd::Ones(npar);
		JtWJ = J.transpose() * ( J.array() * ( weight_sq * onesV.transpose() ).array() ).matrix();  		
	    JtWdy = J.transpose() * ( weight_sq.array() * delta_y.array() ).matrix();
	}catch(const char* msg){
		std::cout << msg << std::endl;		
	}
}



int curveFitting::lm_fitting(double *pr, double *chiSq, int *status0){
	try{
		//initialise parameters
		int stop = 0;
		int converged = 0;
		iteration = 0;

		VectorXd p_old = par;
		VectorXd y_old = y;

		double x2 = 1e-3/std::numeric_limits<double>::epsilon();
		double x2_old = 1e-3/std::numeric_limits<double>::epsilon();
		double dx2 = 1.0;
		
		//initialize Jacobian with finite difference calculation
		VectorXd y_hat = VectorXd::Zero(npnts);
		lm_matx(p_old,y_old,dx2, x2,y_hat);
		
		double lambda = opts->lambda_0;
		x2_old = x2;
				
	    if ( (JtWdy.array().abs()).maxCoeff() < opts->epsilon_1 ){
			if (verbose) {
				std::cout << " *** Your Initial Guess is Extremely Close to Optimal ***" << std::endl;
				std::cout << " *** epsilon_1 = " << opts->epsilon_1 << std::endl;
			}
		   	stop = 1;
			*status0 = 1;
			*chiSq = x2;
	    }
		
	    if(verbose){
			cout << "----------------------------------" << endl;
	    	cout << "par=" << par.transpose() << endl;
			cout << "J=" << J << endl;
			cout << "JtWdy=" << JtWdy.transpose() <<endl;
			cout << "JtWJ=" << JtWJ <<endl;
	    }
		//main loop
		while(!stop && iteration<=opts->maxIter){
			iteration += 1;

			
			VectorXd h(npar);

			//inverse based on the LU decomposition
			//h = ( JtWJ + lambda* MatrixXd(JtWJ.diagonal().asDiagonal()) ).inverse() * JtWdy;
			//cout << "LU=" << (JtWJ + lambda*MatrixXd(JtWJ.diagonal().asDiagonal())).inverse() * JtWdy << endl;
			
			//inverse based on cholesky decomposition LDL
			h = (JtWJ + lambda* MatrixXd(JtWJ.diagonal().asDiagonal())).ldlt().solve(JtWdy);
			//cout << "LDL=" << (JtWJ + lambda*MatrixXd(JtWJ.diagonal().asDiagonal())).ldlt().solve(JtWdy) << endl;

			VectorXd p_try = par + h;
			VectorXd y_try(npnts);
			func( p_try, t, y_try);

			VectorXd delta_y(npnts);
		    delta_y = y - y_try;			// residual error between model and data
			
			double x2_try;					//Chi-squared error criteria
		    x2_try = delta_y.transpose() * ( delta_y.array() * weight_sq.array() ).matrix();
						
            double rho = (x2 - x2_try) / ( h.transpose() * (lambda*MatrixXd(JtWJ.diagonal().asDiagonal())*h + JtWdy) );
			
			if(rho > opts->epsilon_4){ 		//it is significant better
				double dX2 = x2-x2_old;
				x2_old = x2;
				p_old = par;
				y_old = y_hat;
				par = p_try; 				//accept p_try
				
		        lm_matx(p_old,y_old,dX2,x2,y_hat);
				
				//decrease lambda ==> Gauss-Newton method
				lambda = std::max(lambda/opts->lambda_DN_fac,1.0e-7);
				
			} else{  						//it is not better -- do not accept p_try
				x2 = x2_old;
				
				if ( (iteration%(2*npar))==0 ){ //rank-1 update of Jacobian
					double tmp = x2;
					lm_matx(p_old,y_old,-1.0,tmp,y_hat);
				}
				
				//increase lambda  ==> gradient descent method
				lambda = std::min(lambda*opts->lambda_UP_fac,1.0e7);
			}
			
			//check convergence
			if( (JtWdy.array().abs()).maxCoeff() < opts->epsilon_1 &&  iteration > 2 ){
				if(verbose){
					std::cout <<  "**** Convergence in r.h.s. ('JtWdy')  ****"<< std::endl;
					std::cout << "**** epsilon_1 = "<< opts->epsilon_1 << std::endl;
				}
				stop = 1;
				converged += 2;
			}
			if( (abs(h.array()/par.array()).abs()).maxCoeff() < opts->epsilon_2 &&  iteration > 2 ){
				if(verbose){
					std::cout <<  "**** Convergence in Parameters  ****"<< std::endl;
					std::cout << "**** epsilon_2 = "<< opts->epsilon_2 << std::endl;
				}
				converged += 3;
				stop = 1;
			}
			if( x2/(npnts-npar+1.0) < opts->epsilon_3 &&  iteration > 2 ){
				if(verbose){
					std::cout <<  "**** Convergence in Chi-square  ****" << std::endl; 
					std::cout << "**** epsilon_3 = "<< opts->epsilon_3 << std::endl;
				}
				converged += 4;
				stop = 1;
			}
			if(iteration == opts->maxIter && converged==0 ){
				if(verbose){
					std::cout << " !! Maximum Number of Iterations Reached Without Convergence !!" << std::endl;
				}
				stop = 1;
			}

		    if(verbose>1){
				std::cout << "it=" <<iteration<<" chiSq="<<x2<<"  lambda:"<<lambda<<endl;
		    	std::cout << "par=" << par.transpose() << endl;
				std::cout << "h/p=" << h.array()/par.array() << endl;
				std::cout << "J=" << J << endl;
				std::cout << "JtWdy=" << JtWdy.transpose() <<endl;
				std::cout << "JtWJ=" << JtWJ <<endl;
		    }
		
		}   //End of main loop
		*chiSq = x2;
		*status0 = converged;
		Map<VectorXd>(pr, npar) = par;
		return converged;
		
	}catch(const char* msg){
		std::cout << msg << std::endl;	
		return 0;	
	}
}



#endif