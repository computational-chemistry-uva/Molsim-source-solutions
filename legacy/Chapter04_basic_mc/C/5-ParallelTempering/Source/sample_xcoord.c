#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

#define NUMBER_OF_BINS 1000

// sample the distribution of the x-coordinate 
// of the first particle
void SampleXcoord(int Switch)
{
  int i,j;
  static double Gg1[NUMBER_OF_BINS][MAXNUMBEROFSYSTEMS],Gg2;
  FILE *FilePtr;

  switch(Switch)
  {
    case INITIALIZE:
      for(i=0;i<NUMBER_OF_BINS;i++)
        for(j=0;j<MAXNUMBEROFSYSTEMS;j++)
           Gg1[i][j]=0.0;
      Gg2=0.0;
      break;
    case SAMPLE:
      Gg2+=1.0;
      for(j=0;j<NumberOfSystems;j++)
      {
        i=(int)(NUMBER_OF_BINS*Positions[0][j].x/Box);
        if(i<NUMBER_OF_BINS) Gg1[i][j]+=1.0;
      }
      break;
    case WRITE_RESULTS:
      FilePtr=fopen("position_distribution.dat","w");
      for(j=0;j<NumberOfSystems;j++)
      {
        for(i=0;i<NUMBER_OF_BINS;i++)
          if(Gg1[i][j]>0.5) 
            fprintf(FilePtr,"%f %f\n",(i+1.0)/NUMBER_OF_BINS,Gg1[i][j]/Gg2);
        fprintf(FilePtr,"\n");
      }
      fclose(FilePtr);
      break;
  }
}



/* We can't believe that you're reading this!! You should have seen this code 
before we went through it! Unbelievable!! Are you going to change it again? 
Good luck! And let us know please -- send us an e-mail!
Patrick and Edith:
patrick.charbonneau@post.harvard.edu 
beerdsen@science.uva.nl.
*/
