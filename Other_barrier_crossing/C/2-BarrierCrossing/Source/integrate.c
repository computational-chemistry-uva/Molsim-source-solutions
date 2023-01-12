#include <stdlib.h> 
#include <stdio.h> 
#include <math.h> 
#include "system.h"
 
// integrate the equations of motion for an nve system
// velocity verlet integrator
// DeltaUtotal is the energy drift
void Integrate(double *DeltaUdrift)
{
  int i;
  double U,F,Consv,Utotal0,DeltaUtotal,DeltaUcounter;
 
  Force(Xpos,&U,&F);

  DeltaUtotal=0.0;
  DeltaUcounter=0.0;
  Utotal0=0.0;

  // integrate the E.O.M. for Nstep steps
  for(i=0;i<NumberOfSteps;i++)
  {
    Xpos+=Vpos*Tstep+0.5*F*Tstep*Tstep;
    Vpos+=0.5*Tstep*F;

    Force(Xpos,&U,&F);

    Vpos+=0.5*Tstep*F;
    Consv=0.5*Vpos*Vpos+U;

    // used for calculating the energy drift

    if(i==0)
      Utotal0=Consv;
    else
    {
      if(!(i>5&&((Xpos<-1.0&&Vpos<0.0)||(Xpos>1.0&&Vpos>0.0))))
      {
        DeltaUtotal+=fabs((Consv-Utotal0)/Utotal0);
        DeltaUcounter+=1.0;
      }
    }

    Sample(2,i);
  }

  if(DeltaUcounter>0.5)
    *DeltaUdrift=DeltaUtotal/DeltaUcounter;
  else
    *DeltaUdrift=0.0;
}
