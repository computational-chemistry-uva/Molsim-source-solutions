#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"

void Toterg(double *Ener,double *Vir,int BoxID)
{
  //     ---Calculates Total Energy
  
  double Xi, Yi, Zi, Eni, Viri;
  int I, Jb;
  
  *Ener = 0.0;
  *Vir = 0.0;
  for(I = 0;I<Npart-1;I++)
    {
      if(Id[I]==BoxID) 
	{
	  Xi = X[I];
	  Yi = Y[I];
	  Zi = Z[I];
	  Jb = I + 1;
	  Eneri(Xi, Yi, Zi, I, Jb, &Eni,&Viri, BoxID);
	  *Ener = *Ener + Eni;
	  *Vir = *Vir + Viri;
	}
    }
  return;
}
