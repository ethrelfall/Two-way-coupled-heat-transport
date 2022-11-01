// TightCoupledTurbulenceDiffusion1D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "LorenzDiffusion.h"

int main()
{
    std::cout << "Hello World!\n";

    //read parameters
    double rho, sigma, beta;  //Lorenz model parameters
    double diffusionCoeff;

    //if want to read in diffusivity from file
    //FILE* paramFile = fopen("inputs.txt", "r");
    //fscanf(paramFile, "%lf", &diffusionCoeff);
    //fclose(paramFile);

    rho = 28.0;
    sigma = 10.0;
    beta = 8.0 / 3.0;
    diffusionCoeff = 100.0;

    double timeStep = 0.001;
    int numSteps = 50000;
    int outputFrequency = 1;

    LorenzDiffusion lorenzDiffusion = LorenzDiffusion(rho, sigma, beta, diffusionCoeff, timeStep, numSteps);

    for (int i = 0; i < numSteps; i++)
    {
        if (!(i % outputFrequency))
        {
            lorenzDiffusion.OutputDataDiffusion();
        }
        lorenzDiffusion.EvaluateTimestep();
    }
}

