#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..
void WritePdb(FILE *FilePtr)
{
  int i;
  static int CountModel=0;

  CountModel++;
  fprintf(FilePtr,"%s %9d\n","MODEL",CountModel);

  fprintf(FilePtr,"CRYST1   %6.3f   %6.3f   %6.3f  90.00  90.00  90.00 P 1         1\n", 2.0*BOXSIZE, 2.0*BOXSIZE, 2.0*BOXSIZE);
  for(i=0;i<NumberOfParticles;i++)
  {
    fprintf(FilePtr,"%s%7d%s%12d    %8.3lf%8.3lf%8.3lf                      %2s\n",
      "ATOM",i,"  H",i,Positions[i].x*2.0,Positions[i].y*2.0,0.0," H");
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}
