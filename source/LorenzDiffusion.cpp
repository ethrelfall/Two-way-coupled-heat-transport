#include "LorenzDiffusion.h"
#include <algorithm>

LorenzDiffusion::LorenzDiffusion(double rho, double sigma, double beta, double diffusionCoeff, double timeStep, int steps) : mRho(rho), mSigma(sigma), mBeta(beta), mD(diffusionCoeff), mDt(timeStep), mSteps(steps)
{
    mX0 = 0.0;
    mY0 = 1.0;
    mZ0 = 0.0;

    mX = mX0;
    mY = mY0;
    mZ = mZ0;

    mXS = mXDot = 0.0;
    mYS = mYDot = 0.0;
    mZS = mZDot = 0.0;

    mpTimeSeries = new double[mSteps];
    mCurrentStep = 0;

    mResolution = FINITEDIFFRES;
    mDx = 1.0 / mResolution;

    mAlpha = 0.5 * mD * mDt / (mDx * mDx);  //Crank-Nicolson method gives the factor of 0.5

    mTime = 0.0;

    mpCdash = new double[mResolution];
    mpDdash0 = new double[mResolution];
    mpDdash1 = new double[mResolution];
    mpDdash2 = new double[mResolution];
    mpData0 = new double[mResolution];
    mpData1 = new double[mResolution];
    mpData2 = new double[mResolution];

    for (int i = 0; i < mResolution; i++)
    {
        mpCdash[i] = 0.0;
        mpDdash0[i] = 0.0;
        mpDdash1[i] = 0.0;
        mpDdash2[i] = 0.0;
        mpData0[i] = 0.0;
        mpData1[i] = 0.0;
        mpData2[i] = 0.0;
    }

    mpOutFile = fopen("diffusion_heat_flux_vs_t.txt", "w");
    fprintf(mpOutFile, "t,x,y,z\n");
}

LorenzDiffusion::~LorenzDiffusion()
{
    delete[] mpCdash;
    delete[] mpDdash0;
    delete[] mpDdash1;
    delete[] mpDdash2;
    delete[] mpData0;
    delete[] mpData1;
    delete[] mpData2;

    fclose(mpOutFile);

    delete[] mpTimeSeries;
}

void LorenzDiffusion::EvaluateTimestep()
{
    double A[NUMSTAGES];
    double B[NUMSTAGES];
    double C[NUMSTAGES];

    //scheme 2 Carpenter - Kennedy (NUMSTAGES must be 5)
    //see Fourth-Order 2N-Storage Runge-Kutta Schemes, NASA Technical Memorandum 109112
    A[0] = 0.0;
    A[1] = -0.4801594388478;
    A[2] = -1.4042471952;
    A[3] = -2.016477077503;
    A[4] = -1.056444269767;

    B[0] = 0.1028639988105;
    B[1] = 0.7408540575767;
    B[2] = 0.7426530946684;
    B[3] = 0.4694937902358;
    B[4] = 0.1881733382888;

    C[0] = 0.0;
    C[1] = 0.1028639988105;
    C[2] = 0.487989987833;
    C[3] = 0.6885177231562;
    C[4] = 0.9023816453077;

    for (int iStage = 0; iStage < NUMSTAGES; iStage++)
    {
        EvaluateTimestepDiffusion(mX, mY, -OutputNusselt() / mD, C[iStage]);

        //now update the Rayleigh number based on the current temperature
        //this part of the model is somewhat experimental
        mRho = 1.0+27.0 * exp(-mpData2[0]/1.); 

        EvaluateTimeDerivs();

        mXS = A[iStage] * mXS + mDt * mXDot;
        mYS = A[iStage] * mYS + mDt * mYDot;
        mZS = A[iStage] * mZS + mDt * mZDot;
        mX += B[iStage] * mXS;
        mY += B[iStage] * mYS;
        mZ += B[iStage] * mZS;
    }
    mCurrentStep++;
}

void LorenzDiffusion::EvaluateTimeDerivs()
{
    mXDot = mSigma * (mY - mX);
    mYDot = mX * (mRho - mZ) - mY;
    mZDot = mX * mY - mBeta * mZ;
}

double LorenzDiffusion::OutputNusselt()
{
    mpTimeSeries[mCurrentStep] = 1.0 + 2.0 * mZ / mRho;  //record time series

    return 1.0 + 2.0 * mZ / mRho;  //actual Nusselt number
}

void LorenzDiffusion::EvaluateTimestepDiffusion(double LHS_val0, double LHS_val1, double LHS_val2, double timeStepScaleFac)
{
    //Gauss elimination to solve tridiagonal system
    //see e.g. The Nature of Mathematical Modelling, N. Gershenfeld, p.83
    //uses Crank-Nicolson to improve time accuracy to second-order

    //timestep scaling
    double alpha = timeStepScaleFac * mAlpha;

    //mpData0[0] = LHS_val;  //for Dirichlet BC LHS; not needed if Neumann LHS 
    mpData0[mResolution - 1] = 0.0;              //Dirichlet BC RHS (zero-temperature)
    mpData1[mResolution - 1] = 0.0;              //Dirichlet BC RHS (zero-temperature)
    mpData2[mResolution - 1] = 0.0;              //Dirichlet BC RHS (zero-temperature)

    //For Crank-Nicolson data needs to be transformed into data + alpha/2 * matrix * data
    double* pDataCopy0 = new double[mResolution];
    double* pDataCopy1 = new double[mResolution];
    double* pDataCopy2 = new double[mResolution];
#if 1
    for (int i = 0; i < mResolution; i++)
    {
        pDataCopy0[i] = mpData0[i];
        pDataCopy1[i] = mpData1[i];
        pDataCopy2[i] = mpData2[i];
    }

    //Neumann BC LHS.  Was derived using central difference at point zero
    //and T_{1}-T_{-1} = dx * gradient to derive value for putative T_{-1}
    mpData0[0] = pDataCopy0[0] + 2.0 * alpha * (-pDataCopy0[0] + pDataCopy0[1] - 2.0 * mDx * LHS_val0);
    mpData1[0] = pDataCopy1[0] + 2.0 * alpha * (-pDataCopy1[0] + pDataCopy1[1] - 2.0 * mDx * LHS_val1);
    mpData2[0] = pDataCopy2[0] + 2.0 * alpha * (-pDataCopy2[0] + pDataCopy2[1] - 2.0 * mDx * LHS_val2);

    for (int i = 1; i < mResolution - 1; i++)
    {
        mpData0[i] = mpData0[i] + alpha * (pDataCopy0[i + 1] - 2.0 * pDataCopy0[i] + pDataCopy0[i - 1]);
        mpData1[i] = mpData1[i] + alpha * (pDataCopy1[i + 1] - 2.0 * pDataCopy1[i] + pDataCopy1[i - 1]);
        mpData2[i] = mpData2[i] + alpha * (pDataCopy2[i + 1] - 2.0 * pDataCopy2[i] + pDataCopy2[i - 1]);
    }
#endif

    //for Neumann BC LHS
    mpCdash[0] = -2.0 * alpha / (1.0 + 2.0 * alpha);
    mpDdash0[0] = mpData0[0] / (1.0 + 2.0 * alpha);
    mpDdash1[0] = mpData1[0] / (1.0 + 2.0 * alpha);
    mpDdash2[0] = mpData2[0] / (1.0 + 2.0 * alpha);

    //for Dirichlet BC LHS
    //mpCdash[0] = 0.0;
    //mpDdash0[0] = mpData0[0];

    for (int i = 0; i < mResolution - 2; i++)
    {
        double recipDenom = 1.0 / (1.0 + 2.0 * alpha + alpha * mpCdash[i]);
        mpCdash[i + 1] = -alpha * recipDenom;
        mpDdash0[i + 1] = (mpData0[i + 1] + alpha * mpDdash0[i]) * recipDenom;
        mpDdash1[i + 1] = (mpData1[i + 1] + alpha * mpDdash1[i]) * recipDenom;
        mpDdash2[i + 1] = (mpData2[i + 1] + alpha * mpDdash2[i]) * recipDenom;
    }

    mpCdash[mResolution - 1] = 0;
    mpDdash0[mResolution - 1] = mpData0[mResolution - 1];
    mpDdash1[mResolution - 1] = mpData1[mResolution - 1];
    mpDdash2[mResolution - 1] = mpData2[mResolution - 1];

    for (int i = 0; i < mResolution - 1; i++)
    {
        int idx = mResolution - i - 2;
        mpData0[idx] = mpDdash0[idx] - mpCdash[idx] * mpData0[idx + 1];
        mpData1[idx] = mpDdash1[idx] - mpCdash[idx] * mpData1[idx + 1];
        mpData2[idx] = mpDdash2[idx] - mpCdash[idx] * mpData2[idx + 1];
    }

    mTime += mDt*timeStepScaleFac;

    delete[] pDataCopy0;
    delete[] pDataCopy1;
    delete[] pDataCopy2;
}

void LorenzDiffusion::OutputDataDiffusion()
{
    //output gradients at LHS and RHS
    //fprintf(mpOutFile, "%.6e, %.6e, %.6e\n", mTime, mZ); 
    //fprintf(mpOutFile, "%.6e, %.6e, %.6e\n", mTime, OutputNusselt());
     
    //temperatures at endpoints
    //fprintf(mpOutFile, "%.6e, %.6e, %.6e\n", mTime, mpData2[0], mpData2[mResolution - 1]);

    //output heat fluxes
    fprintf(mpOutFile, "%.6e, %.6e, %.6e\n", mTime, mD*(mpData2[0] - mpData2[1]) / mDx, mD*(mpData2[mResolution - 2] - mpData2[mResolution - 1]) / mDx);
}