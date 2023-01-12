#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..
void WritePdb(FILE *FilePtr)
{
  int i;
  static int Countmodel=0;

  Countmodel++;
  fprintf(FilePtr,"%s %9d\n","MODEL",Countmodel);

  fprintf(FilePtr,"CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1         1\n");
  for(i=0;i<NumberOfParticles;i++)
  {
    fprintf(FilePtr,"%s%7d%s%12d    %8.3lf%8.3lf%8.3lf\n",
      "ATOM",i,"  H",i,Positions[i].x*2.0,Positions[i].y*2.0,Positions[i].z*2.0);
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}
