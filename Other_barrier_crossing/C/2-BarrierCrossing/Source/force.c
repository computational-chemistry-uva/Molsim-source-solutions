#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// calculate the forces and potential energy
void Force(double x, double *U, double *F)
{
  *U=0.0;
  *F=0.0;

  if(x>=-1.0&&x<=1.0)
  {
    *U=0.5*(1.0+cos(M_PI*x));
    *F=0.5*M_PI*sin(M_PI*x);
  }
} 
