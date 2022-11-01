#pragma once

#include <iostream>

#define NUMSTAGES 5
#define FINITEDIFFRES 200

#pragma warning(disable:4996)  //needed in VS2022 to stop security check error

class LorenzDiffusion
{
public:
	LorenzDiffusion(double rho, double sigma, double beta, double diffusionCoeff, double timeStep, int steps);
	~LorenzDiffusion();
	void EvaluateTimestep();
	void EvaluateTimeDerivs();
	double OutputNusselt();

	//these were originally in diffusion solver
	//this version diffuses three time series independently, the third one is the temperature
	void EvaluateTimestepDiffusion(double LHS_val0, double LHS_val1, double LHS_val2, double timeStepScaleFac);
	void OutputDataDiffusion();

private:
	//model parameters
	//all the following model parameters are >0
	//see https://www2.physics.ox.ac.uk/sites/default/files/profiles/read/lect6-43147.pdf
	double mRho;    //Rayleigh number
	double mSigma;  //Prandtl number
	double mBeta;   //coupling strength

	//initial data
	double mX0;
	double mY0;
	double mZ0;

	//computational parameters
	//time step
	double mDt;

	//model state storage appropriate to LSRK
	double mX;
	double mY;
	double mZ;
	double mXDot;
	double mYDot;
	double mZDot;
	double mXS;
	double mYS;
	double mZS;

	//storage for time series
	int mSteps;
	int mCurrentStep;
	double* mpTimeSeries;

	//model parameters for diffusion
	double mD;  //diffusion coefficient

	//computational parameters for diffusion
	int mResolution;
	double mDx;

	//dimensionless group appearing in matrix used for time-evolution
	double mAlpha;

	double mTime;

	//coefficients used to solve tridiagonal matrix
	double* mpCdash;
	double* mpDdash0;
	double* mpDdash1;
	double* mpDdash2;

	//storage for model state
	double* mpData0;
	double* mpData1;
	double* mpData2;

	//I/O
	FILE * mpOutFile;

};