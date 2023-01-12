#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..
void WritePdb(FILE *FilePtr)
{
  int i;
  static int CountModel=0,CountAtom=0;

  CountModel++;
  fprintf(FilePtr,"%s %9d\n","MODEL",CountModel);

  for(i=0;i<NumberOfParticles;i++)
  {
    CountAtom++;
    fprintf(FilePtr,"%s%7d%s%12d%8.3lf%8.3lf%8.3lf\n",
      "ATOM",CountAtom,"  H",CountAtom,Positions[i].x*2.0,Positions[i].y*2.0,0.0);
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}
