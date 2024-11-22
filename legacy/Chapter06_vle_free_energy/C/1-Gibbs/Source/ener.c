#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "potential.h"

int Ener(double *En, double *Vir, double R2, int BoxID)
{
  // ---Calculates Energy (En) And Virial (Vir) For Given
  //    Distance Squared Between (R2) Two Particles
  
  double R2i, R6i;
  
  if(R2<=Rc2[BoxID])
    {
      R2i = Sig2/R2;
      R6i = R2i*R2i*R2i;
      *En  = Eps4*(R6i*R6i-R6i);
      *Vir = Eps48*(R6i*R6i-0.5*R6i);
    }
  else
    {
      *En  = 0.0;
      *Vir = 0.0;
    }
  return 0;
}      
