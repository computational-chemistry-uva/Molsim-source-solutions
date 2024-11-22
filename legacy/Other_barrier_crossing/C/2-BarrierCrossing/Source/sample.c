#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

#define Maxx 50000

void Sample(int Switch, int Ttt)
{
  int i;
  static double Kt[Maxx],Ks[Maxx];
  FILE *FilePtr;

  switch(Switch)
  {
    case 1:
      // initialize everything
      for(i=0;i<Maxx;i++)
      {
        Kt[i]=0.0;
        Ks[i]=0.0;
      }
      break;
    case 2:
      // Compute the Transmission Coeffient
      // see also Case Study 23 on page 440/441 (equation c)
      // the maximum time equals Maxx

      // begin modification
        if(Ttt<Maxx)
        {
          if(Xpos>=Qstar)
             Kt[Ttt]+=Theta;
          Ks[Ttt]+=fabs(Theta);
        }
      // end modification
      break;
    case 3:
      // write results
      FilePtr=fopen("results.dat","w");
      for(i=0;i<Maxx;i++)
        if(i<NumberOfSteps&&Ks[i]>0.5)
          fprintf(FilePtr,"%lf %lf\n",i*Tstep,2.0*Kt[i]/Ks[i]);
      fclose(FilePtr);
      break;
  }
}
