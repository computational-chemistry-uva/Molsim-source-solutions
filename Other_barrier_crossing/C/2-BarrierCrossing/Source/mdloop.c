#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

#define CycleMultiplication 1000

void MdLoop(void)
{
  int i,j;
  double EnergyDrift,CumEnergyDrift,CumCycles;

  Sample(1,0);

  CumEnergyDrift=0.0;
  CumCycles=0.0;

  for(i=0;i<NumberOfCycles;i++)
    for(j=0;j<CycleMultiplication;j++)
    {
      // generate initial coordinates
      //
      // perform 1000*Ncycle MD simulations with 
      // different initial conditions
      //
      // Xpos   = starting position
      // Vpos   = starting velocity
      // Theta  = storage for initial velocity
      // Qstar  = place of the dividing surface
      //
      // To Program:
      //
      // -Generate Initial Position/Velocity
      // -Integrate The Equations Of Motion By A 
      //  Function Call To Subroutine Integrate 
      // -Beware That In Integrate The Subroutine
      //  Sample Is Called !!
      //
      // The Averaged Energy Drift (Over All
      // Simulations) Equals CumEnergyDrift/CumCycles

      // start modification
        Xpos=Qstar;
        Vpos=RandomVelocity(Temperature);
        Theta=Vpos;

        Integrate(&EnergyDrift);

        CumEnergyDrift+=EnergyDrift;
        CumCycles+=1.0;
      // end modification

    }

  Sample(3,0);
  printf("Av. Energy Drift      : %f\n",CumEnergyDrift/CumCycles);
}
