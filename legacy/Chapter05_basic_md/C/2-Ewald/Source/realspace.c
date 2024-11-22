#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

// real Part; also direct calculation 
// loop over all particle pairs
void Realspace(double *Ureal)
{
  int i,j;
  VECTOR dr;
  double U,r,r2;

  U=0.0;

  // start modification
  // 1. For all particle pairs calculate the distance in x, y, and z.
  // 2. Apply periodic boundary conditions where necessary.
  // 3. Calculate the real-space contribution to the energy. 
  // (Use the "ErrorFunctionComplement(x)" function.) 

    for(i=0;i<NumberOfParticles-1;i++)
    {
      for(j=i+1;j<NumberOfParticles;j++)
      {
        dr.x=Positions[i].x-Positions[j].x;
        dr.y=Positions[i].y-Positions[j].y;
        dr.z=Positions[i].z-Positions[j].z;

        // apply boundary conditions
        dr.x-=Box*rint(dr.x/Box);
        dr.y-=Box*rint(dr.y/Box);
        dr.z-=Box*rint(dr.z/Box);

        r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(r2);
        if(r<0.5*Box)
        {
          U+=Charges[i]*Charges[j]*ErrorFunctionComplement(Alpha*r)/r;
        }
      }
    }

  // end  modification
  *Ureal=U;
}
