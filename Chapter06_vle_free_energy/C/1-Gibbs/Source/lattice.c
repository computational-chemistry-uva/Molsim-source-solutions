#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"   

int Lattice(void)
{
  // ---Place `Npart' Particles On A Lattice With Density 'Rho'
  // --Half The Number In Box 1 And The Other Half In Box 2
  
  int I, J, K, Itel, N, BoxID;
  double Del;
  
  Del      = pow((pow(Box[0],3)),(1.0/3.0));
  Npbox[0] = Npart/2;
  Npbox[1] = Npbox[0];
  
  if(Npbox[0]+Npbox[1]!=Npart)
    {
      printf("Error: Npart\n");
      exit(0);
    }
  
  printf("Generate Simple Cubic Lattice\n");
  
  N = (int)(pow(Npart,(1.0/3.0)) + 1);
  if(N==0) N = 1;
  Del = Del/(double)N;
  Itel = 0;
  for(I=0;I<N;I++)
    {
      for(J=0;J<N;J++)
	{
	  for(K=0;K<N;K++)
	    {
	      for(BoxID=0;BoxID<2;BoxID++)
		{
		  if(Itel<Npart)
		    {
		      X[Itel] = (double)K*Del;
		      Y[Itel] = (double)J*Del;
		      Z[Itel] = (double)I*Del;
		      Id[Itel] = BoxID;
		      Itel++;
		    }
		} 
	    }
	}
    }
  
  printf("Initialisation on Lattice: \n"); 
  printf("%d Particles Placed on a Lattice \n",Itel); 
  
  return(0);
}
