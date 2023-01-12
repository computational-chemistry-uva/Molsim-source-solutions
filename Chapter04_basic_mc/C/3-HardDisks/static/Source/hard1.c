#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

#define CycleMultiplication 1000

// hard disks on a square; random placement

int main(void)
{
  int NumberOfCycles,i,j,k,l,overlap;
  int NumberOfAttemptedMoves,NumberOfAcceptedMoves;
  double r2;
  VECTOR dr;
  FILE *FilePtr;

  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  FilePtr=fopen("movie.pdb","w");

  // read the input parameters
  printf("How many cycles (x %d)? ",CycleMultiplication);
  fscanf(stdin,"%d",&NumberOfCycles);

  printf("How many particles    ?        (always 2 < i < 80 ) ");
  fscanf(stdin,"%d",&NumberOfParticles);

  if((NumberOfParticles>80)||(NumberOfParticles<2))
  {
    printf("input error\n");
    exit(0);
  }

  Sample(INITIALIZE);

  NumberOfAttemptedMoves=0;
  NumberOfAcceptedMoves=0;


  for(l=0;l<NumberOfCycles;l++)
  {
    for(k=0;k<CycleMultiplication;k++)
    {
      NumberOfAttemptedMoves++;

      overlap=FALSE;
      for(i=0;(i<NumberOfParticles)&&(overlap==FALSE);i++)
      {

        // new positions
        Positions[i].x=RandomNumber()*BOXSIZE;
        Positions[i].y=RandomNumber()*BOXSIZE;

        // check distances with previous placed particles
        for(j=0;j<i&&overlap==FALSE;j++)
        {
          dr.x=Positions[i].x-Positions[j].x;
          dr.y=Positions[i].y-Positions[j].y;

          // apply boundary conditions
          if(PBC)
          {
            dr.x-=BOXSIZE*rint(dr.x/BOXSIZE);
            dr.y-=BOXSIZE*rint(dr.y/BOXSIZE);
          }

          r2=SQR(dr.x)+SQR(dr.y);
          if(r2<1.0) overlap=TRUE;
        }
      }

      if(overlap==FALSE)
      {
        NumberOfAcceptedMoves++;
        Sample(SAMPLE);
      }
    }
  }
  if(NumberOfAcceptedMoves>0) Sample(WRITE_RESULTS);
  printf("Fraction succes : %f\n",(double)NumberOfAcceptedMoves/NumberOfAttemptedMoves);

  return 0;
}
