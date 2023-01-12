#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

#define SQR(x) ((x)*(x))

#define Maxtrial 100
#define CycleMultiplication 1000

int main(void)
{
  int Nchoi,i,j,k,Nstep,Ichoi;
      
  double Uold,CumW,Ws,Sumnew,Sumold,CumEnergySq;
  double X1,X2,X3,X1n,X2n,X3n;
  double X1t[Maxtrial],X2t[Maxtrial],X3t[Maxtrial];
  double Ut[Maxtrial],Deltax,CumEnergy,CumCycle,CumAccepted;
 
  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  // read the input parameters
  printf("How Many Cycles (x %d)    ? ",CycleMultiplication);
  fscanf(stdin,"%d",&Nstep);

  printf("How Many Trial Directions ? ");
  fscanf(stdin,"%d",&Nchoi);

  if(Nchoi>Maxtrial)
  {
    printf("Error Nchoi>Maxtrial\n");
    exit(-1);
  }

  printf("Maximum Displacement      ? ");
  fscanf(stdin,"%lf",&Deltax);

  X1=0.0;
  X2=0.0;
  X3=0.0;
  CumEnergy=0.0;
  CumCycle=0.0;
  CumEnergySq=0.0;
  CumAccepted=0.0;

  // start the simulation
  // calculate initial energy

  Uold=SQR(X1)+SQR(X2)+SQR(X3);

  for(i=1;i<=Nstep;i++)
  {
    for(j=1;j<=CycleMultiplication;j++)
    {
      // Cbmc; generate Nchoi Displacements
      // Selectone According To Its Boltzmann Factor...
      // Old Coordinates : X1, X2, X3
      // Trial Coordinates : X1t, Y1t, Z1t
      // Nchoi           : Number Of Trials
      // Sumnew          : Weight New Config.

      Sumnew=0.0;
      for(Ichoi=1;Ichoi<=Nchoi;Ichoi++)
      {
        X1t[Ichoi]=X1+(RandomNumber()-0.5)*Deltax;
        X2t[Ichoi]=X2+(RandomNumber()-0.5)*Deltax;
        X3t[Ichoi]=X3+(RandomNumber()-0.5)*Deltax;
 
        Ut[Ichoi]=exp(-(SQR(X1t[Ichoi])+SQR(X2t[Ichoi])+SQR(X3t[Ichoi])));
        Sumnew+=Ut[Ichoi];
      }

      // select one of the trial directions...
      Ws=RandomNumber()*Sumnew;
      CumW=0;
      Ichoi=0;

      while(CumW<Ws)
      {
        Ichoi++;
        CumW+=Ut[Ichoi];
      }

      // old configuration 
      // generate Nchoi-1 positions around 
      // the new config; the first one is the 
      // old configuration... 
      // sumold = weight of old configuration 

      Sumold=0.0;

      for(k=1;k<=Nchoi;k++)
      {
        if(k==1)
        {
          X1n=X1;
          X2n=X2;
          X3n=X3;
        }
        else
      // start modification
        {
          X1n=X1t[Ichoi]+(RandomNumber()-0.5)*Deltax;
          X2n=X2t[Ichoi]+(RandomNumber()-0.5)*Deltax;
          X3n=X3t[Ichoi]+(RandomNumber()-0.5)*Deltax;
        }
        Sumold+=exp(-(SQR(X1n)+SQR(X2n)+SQR(X3n)));
      // end modification
      }

      // accept or reject this configuration

      if(RandomNumber()<(Sumnew/Sumold))
      {
        X1=X1t[Ichoi];
        X2=X2t[Ichoi];
        X3=X3t[Ichoi];
               
        Uold=SQR(X1)+SQR(X2)+SQR(X3);
        CumAccepted=CumAccepted+1.0;
      }

      CumEnergy+=Uold;
      CumEnergySq+=SQR(Uold);
      CumCycle+=1.0;
    }
  }

  // write results
  printf("Average Energy          : %f\n",CumEnergy/CumCycle);
  printf("Sigma <E>               : %f\n",sqrt((CumEnergySq/CumCycle)-SQR(CumEnergy/CumCycle)));
  printf("Fraction Accepted Moves : %f\n",CumAccepted/CumCycle);

  return 0;
}
