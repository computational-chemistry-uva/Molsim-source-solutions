#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

#define CycleMultiplication 1000

int PBC;
// hard disks on a square; random placement

int main(void)
{
  int NumberOfCycles,NumberOfInitializationCycles;
  int i,j,k,SelectedParticle,NumberOfSamples,overlap;
  int NumberOfAttemptedMoves,NumberOfAcceptedMoves;
  double r2,MaximumDisplacement;
  VECTOR PosNew,dr;
  FILE *FilePtr;
      
  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  FilePtr=fopen("movie.pdb","w");

  // read the input parameters
  printf("How many cycles (x %d)      ? ",CycleMultiplication);
  fscanf(stdin,"%d",&NumberOfCycles);

  printf("How many init. cycles ?        (example: 100      ) ");
  fscanf(stdin,"%d",&NumberOfInitializationCycles);


  printf("How many particles    ?        (always 2 < i < 80 ) ");
  fscanf(stdin,"%d",&NumberOfParticles);

  printf("Maximum displacement  ?        (disp > 0          ) ");
  fscanf(stdin,"%lf",&MaximumDisplacement);

  printf("Period boundary conditions ?   (0 or 1) ");
  fscanf(stdin,"%d",&PBC);

  if((NumberOfParticles>MAX_PARTICLES)||(NumberOfParticles<2))
  {
    printf("input error\n");
    exit(0);
  }

  Sample(INITIALIZE);

  NumberOfAttemptedMoves=0;
  NumberOfAcceptedMoves=0.0;
  NumberOfSamples=0;

  // put particles on a lattice 

  k=0;
  for(i=0;i<(int)(BOXSIZE-1.0);i++)
  {
    for(j=0;j<(int)(BOXSIZE-1.0);j++)
    {
      Positions[k].x=i*BOXSIZE/(BOXSIZE-1.0);
      Positions[k].y=j*BOXSIZE/(BOXSIZE-1.0);
      k++;
    }
  }

  for(i=0;i<NumberOfCycles;i++)
  {
    for(j=0;j<CycleMultiplication;j++)
    {
      NumberOfAttemptedMoves++;

      // select particle at random 
      SelectedParticle=RandomNumber()*NumberOfParticles;

      PosNew.x=Positions[SelectedParticle].x+(RandomNumber()-0.5)*MaximumDisplacement;
      PosNew.y=Positions[SelectedParticle].y+(RandomNumber()-0.5)*MaximumDisplacement;

      // always reject when particle is out of the box
      //if(PBC||(MIN(PosNew.x,PosNew.y)>=0.0&&MAX(PosNew.x,PosNew.y)<BOXSIZE))
      if(PBC||((!PBC)&&(MIN(PosNew.x,PosNew.y)>=0.0&&MAX(PosNew.x,PosNew.y)<BOXSIZE)))
      {
        // see if there is an overlap with any other particle
        overlap=FALSE;
        for(k=0;((k<NumberOfParticles));k++)
        {
          if(k!=SelectedParticle)
          {
            dr.x=PosNew.x-Positions[k].x;
            dr.y=PosNew.y-Positions[k].y;

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
          Positions[SelectedParticle].x=PosNew.x;
          Positions[SelectedParticle].y=PosNew.y;
        }
      }
      if(((j%5)==0)&&(i>NumberOfInitializationCycles))
      {
        Sample(SAMPLE);
        NumberOfSamples++;
      }
    }
    if((i%5)==0) WritePdb(FilePtr);
  }

  if(NumberOfSamples>0) Sample(WRITE_RESULTS);
  printf("Fraction Succes : %f\n",(double)NumberOfAcceptedMoves/NumberOfAttemptedMoves);

  fclose(FilePtr);
  return 0;
}
