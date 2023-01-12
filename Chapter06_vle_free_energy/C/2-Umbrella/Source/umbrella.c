#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

#define SQR(x) ((x)*(x))

#define Max 100000
#define Width 0.001
#define CycleMultiplication 1000

int main(void)
{
  int i,j,k;
  int NumberOfCycles,NumberOfAcceptedMoves,NumberOfAttempts;
  double LeftBoundary,RightBoundary,Xold,Xnew,Uold,Unew;
  double MaximumDisplacement,PositionCounter[2*Max+1];
  FILE *FilePtr;

  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  printf("How many cycles (x %d) ? ",CycleMultiplication);
  fscanf(stdin,"%d",&NumberOfCycles);

  printf("Maximum Displacement  ? ");
  fscanf(stdin,"%lf",&MaximumDisplacement);

  printf("Minimum Value Of X    ?  ");
  fscanf(stdin,"%lf",&LeftBoundary);

  printf("Maximum Value Of X    ?  ");
  fscanf(stdin,"%lf",&RightBoundary);
 
  if(LeftBoundary>=RightBoundary)
  {
    printf("LeftBoundary should be smaller than RightBoundary\n");
    exit(-1);
  }

  // generate initial position and energy
  Xold=0.5*(RightBoundary+LeftBoundary);
  Uold=SQR(Xold);
  Xnew=0.0;
  Unew=0.0;
  NumberOfAcceptedMoves=0;
  NumberOfAttempts=0;
      
  for(i=-Max;i<Max;i++)
    PositionCounter[i+Max]=0.0;
         
  printf("\n");
  printf("Calculating......\n");
  printf("\n");

  // loop over all cycles
  for(i=0;i<NumberOfCycles;i++)
  {
    for(j=0;j<CycleMultiplication;j++)
    {
      Xnew=Xold+(RandomNumber()-0.5)*MaximumDisplacement;
      Unew=SQR(Xnew);
      NumberOfAttempts++;
      k=(int)(Xold/Width);

      if(k<-Max||k>Max)
      {
        printf("Out of range : j !!!\n");
        exit(-1);
      }
      else
      {
        if(k==0) // Bin 0 is double the size of the other bins.
          PositionCounter[k+Max]+=0.5;
        else
          PositionCounter[k+Max]+=1.0;
      }

      // reject when outside slice
      if(Xnew>LeftBoundary&&Xnew<RightBoundary)
      {
        if(RandomNumber()<exp(Uold-Unew))
        {
           Xold=Xnew;
           Uold=Unew;
           NumberOfAcceptedMoves++;
        }
      }
    }
  }

  // write results
  FilePtr=fopen("Umbrella.dat","w");
  for(i=-Max;i<Max;i++)
     if(PositionCounter[i+Max]>0.5)
       fprintf(FilePtr,"%f %f\n",i*Width,PositionCounter[i+Max]/NumberOfAttempts);
  fclose(FilePtr);

  printf("\n");
  printf("Fraction of accepted displacements : %f\n",(double)NumberOfAcceptedMoves/NumberOfAttempts);
  printf("\n");
  printf("Histogram Written To Umbrella.dat\n");

  return 0;
}
