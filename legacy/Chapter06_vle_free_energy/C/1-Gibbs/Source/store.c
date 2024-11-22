#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "parameter.h"
#include "conf.h"
#include "system.h"

void Store(FILE* Iout,double Dr,double Vmax)
{
  //C     Writes Configuration To Disk
  int I;
 
  fprintf(Iout,"%lf\t%lf\t%lf\t%lf\n",Box[0],Hbox[0],Box[1],Hbox[1]);
  fprintf(Iout,"%d\t%d\t%d\n",Npart,Npbox[0],Npbox[1]);
  fprintf(Iout,"%lf\t%lf\n",Dr,Vmax);
  for(I=0;I<Npart;I++)
    {
      fprintf(Iout,"%lf\t%lf\t%lf\t%d\n", X[I], Y[I], Z[I], Id[I]);
    }
  //   Rewind (Iout)
  return;
}
