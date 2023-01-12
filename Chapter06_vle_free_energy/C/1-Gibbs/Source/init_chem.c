#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "chem.h"

int Init_Chem(int Switch)
{
  //  ---Initialize And Calculate Chemical Potentials
  int BoxID;
  
  if(Switch==0)
    {
      //    ---Initialize
      for(BoxID=0;BoxID<2;BoxID++)
	{   
	  Chp[BoxID] = 0.0;
	  Ichp[BoxID] = 0;
	}
    }
  else if(Switch==2)
    {
      //    ---Print Final Results
      for(BoxID=0;BoxID<2;BoxID++)
	{
	  if(Ichp[BoxID]!=0)
	    {
	      Chp[BoxID] = -log(Chp[BoxID]/(double)(Ichp[BoxID]))/Beta;
	    }
	}
      printf("------------------- \n");
      printf("Chemical Potentials \n");
      printf("Number of Samples: %d\n", ((Ichp[0]+Ichp[1])/2));
      printf("Box 1: %g\n",Chp[0]);
      printf("Box 2: %g\n",Chp[1]);
    }
   else
     {
       printf("Error: Init_Chem\n");
       exit(0);
     } 
  return 0;
}
