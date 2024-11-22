#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"   

int Eneri(double Xi, double Yi, double Zi, int I, int Jb, double *En, double *Vir, int BoxID)
{
  double Dx, Dy, Dz, R2, Virij, Enij;
  int  J;
  
  *En = 0.0;
  *Vir = 0.0;
  for(J=Jb;J<Npart;J++)
    {
      if(Id[J]==BoxID)
	{
	  if(J!=I)
	    {
	      Dx = Xi - X[J];
	      Dy = Yi - Y[J];
	      Dz = Zi - Z[J]; 
	      if(Dx>Hbox[BoxID])
		{
		  Dx = Dx - Box[BoxID];
		}
	      else if(Dx<-Hbox[BoxID])
		{
		  Dx = Dx + Box[BoxID];
		}
	      if(Dy>Hbox[BoxID])
		{
		  Dy = Dy - Box[BoxID];
		}
	      else if(Dy<-Hbox[BoxID])
		{
		  Dy = Dy + Box[BoxID];
		}
	      if(Dz>Hbox[BoxID])
		{
		  Dz = Dz - Box[BoxID];
		}
	      else if(Dz<-Hbox[BoxID]) 
		{
		  Dz = Dz + Box[BoxID];
		}
	      R2 = Dx*Dx + Dy*Dy + Dz*Dz;
	      Ener(&Enij, &Virij, R2, BoxID);
	      *En+= Enij;
	      *Vir+= Virij;
	    }
	}
    }
  
  return 0;
}
