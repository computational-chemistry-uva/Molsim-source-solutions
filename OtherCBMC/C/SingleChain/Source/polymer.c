//Samples configurations of a single polymer chain; the polymer chain is generated either using a "normal" MC  scheme or CBMC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

int main(void)
{
  FILE *FilePtr;

  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  //  Read input file (see file ../Source/run)
  FilePtr=fopen("input","r"); 
  fscanf(FilePtr,"%d %d %d %d\n",&NumberOfSteps,&NumberOfInitializationSteps,
         &ChainLength,&NumberOfTrialPositions);
  fscanf(FilePtr,"%d\n",&Lcbmc);
  fscanf(FilePtr,"%lf %lf %lf %lf %lf\n",&Temperature,&Rcut,&A,&Kt,&Theta0);
  Beta=1.0/Temperature;

  printf("Cbmc program of a single chain\n\n");

  //Output input parameters  
  if(Lcbmc)
    printf("Cbmc is used\n");
  else
    printf("Random chains are grown\n");

  printf("\n");
  printf("Number of cycles  : %d\n",NumberOfSteps);
  printf("Number of init.   : %d\n",NumberOfInitializationSteps);
  printf("\n");
  printf("Temperature       : %f\n",Temperature);
  printf("Beta              : %f\n",Beta);
  printf("Chain length      : %d\n",ChainLength);
  printf("Number of trials  : %d\n",NumberOfTrialPositions);
  printf("\n");
  printf("Rcut              : %f\n",Rcut);
  printf("A                 : %f\n",A);
  printf("Kt                : %f\n",Kt);
  printf("Theta0 [radians]   : %f\n",Theta0);
  printf("\n");

  //If certain input parameters are no good, warn the user and kill the program
  if(ChainLength>MAX_CHAIN_LENGTH)
  {
    printf("Chain is too long !!!\n");
    exit(0);
  }

  if(NumberOfTrialPositions>MAX_TRIALS)
  {
    printf("Too many trial positions !!!\n"); 
    exit(0);
  }
  //Let the action begin! 
  Mcloop();

  return 0;
} 
