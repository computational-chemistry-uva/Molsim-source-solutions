#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
 
// read in system information
void ReadData(void)
{      
  int i;
  double U,F,Suma;
  FILE *FilePtr;

  FilePtr=fopen("input","r");
  fscanf(FilePtr,"%d %d %lf %lf %lf\n",&NumberOfCycles,&NumberOfSteps,&Tstep,&Temperature,&Qstar);

  printf("Barrier Crossing Of A Particle Over A Energy Barrier\n");
  printf("\n");
  printf("Number Of Cycles      : %d\n",NumberOfCycles);
  printf("Number Of Md Steps    : %d\n",NumberOfSteps);
  printf("Temperature           : %f\n",Temperature);
  printf("Position Of Q(Star)   : %f\n",Qstar);
  printf("Timestep              : %f\n",Tstep);
  printf("\n");

  if(Qstar<-1.7||Qstar>1.7)
  {
    printf("Wrong Value Qstar\n");
    printf("-0.7<Qstar<0.7\n");
    exit(-1);
  }

  // Suma = Integral Over A
  Suma=0.0;
  for(i=0;i<=100000;i++)
  {
    Xpos=(i-1)*0.0001-4.0;
    Force(Xpos,&U,&F);
    if(Xpos<=Qstar) 
      Suma+=exp(-U/Temperature)*0.0001;
  }
  Force(Qstar,&U,&F);
  printf("P(Qstar)/P(A)         : %18.10f\n",exp(-U/Temperature)/Suma);
}
