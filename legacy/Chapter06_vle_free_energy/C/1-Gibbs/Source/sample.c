#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "conf.h"
#include "system.h"

void Sample(int I,double *En,double *Vir)
{
  //C     Writes Quantities To File
  int BoxID;
  double Enp[2], Press[2], Vol, Rho[2];
  FILE* fileptr;
  
  for(BoxID = 0; BoxID<2;BoxID++)
    {
      Vol = pow(Box[BoxID],3);
      Rho[BoxID] = (double)(Npbox[BoxID])/Vol;
      Press[BoxID] = Rho[BoxID]/Beta + Vir[BoxID]/(3.0*Vol);
         
      if (Npbox[BoxID]!=0)
	{
	  Enp[BoxID] = En[BoxID]/(double)(Npbox[BoxID]);
	}
      else
	{
	  Enp[BoxID] = 0.;
	}
    }
  
  fileptr=fopen("fort.66","a");
  fprintf(fileptr,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",I, Enp[0],Enp[1],Press[0],Press[1],Rho[0],Rho[1]);
  fclose(fileptr);
  fileptr=fopen("fort.44","a");
  fprintf(fileptr,"%d\t%lf\t%d\t%lf\n",Npbox[0],pow(Box[0],3),Npbox[1],pow(Box[1],3));
  fclose(fileptr);
  fileptr=fopen("fort.45","a");
  fprintf(fileptr,"%d\t%lf\t%lf\n",I,Npbox[0]/pow(Box[0],3),Npbox[1]/pow(Box[1],3));
  fclose(fileptr);
  return;
}
