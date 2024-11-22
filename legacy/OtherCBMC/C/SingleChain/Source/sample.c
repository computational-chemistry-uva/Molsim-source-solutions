#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

#define MAX_NUMBER_OF_BINS 501
//Samples the end-to-end distance of the chain
void Sample(int Choice,double r,double W)
{
  int i;
  static double Normalize;
  static double histogram[MAX_NUMBER_OF_BINS];
  static double dR;
  static double averageR2;
  FILE *FilePtr;
//    Choice: parameter that specifies what the subroutine will do:  initiallize,    sample or  output final data 
//         R:      end-to-end distance of the current chain configuration
//         W:      weight of the current configuration.  The sampling scheme used in subroutine mcloop is called a dynamic scheme; in dynamics schems W is always 1.  The subroutine Sample is written so it can also be used in static schemes where W is not always 1.

//       normalize:   the variable name says it all...
//       histogram:   the variable name says it all...
//       dR:          inverse bin width in the histogram (converted to bin width when writing the output)
//       averageR2:   use it to calculate <R^2>

  switch(Choice)
  {
    case INITIALIZE:      // initialize variables for this subroutine
      dR=(MAX_NUMBER_OF_BINS-1.0)/MAX_CHAIN_LENGTH;
//        If the number of bins in the histogram is MAX_NUMBER_OF_BINS, why is dR not
//        MAX_NUMBER_OF_BINS/MAX_CHAIN_LENGTH ?
//        Answer: As we use function int to populate the histogram, 
//        if dR=MAX_NUMBER_OF_BINS/MAX_CHAIN_LENGTH then a fully 
//        stretched chain of length  MAX_CHAIN_LENGTH  will fall outside the 
//        histogram and the program will crash.
      Normalize=0.0;
      averageR2=0.0;
      for(i=0;i<MAX_NUMBER_OF_BINS;i++)
        histogram[i]=0.0;
      break;
    case SAMPLE:          // acumulate statistics on the chain end-to-end distances
      i=(int)(r*dR);
      if(i<MAX_NUMBER_OF_BINS)
        histogram[i]+=W;
      Normalize+=W;
//    Start Modification to calculate R2: sum over all  R*R (with the correct weight). 
//      Use variable averageR2 to store the sum
         averageR2+=r*r*W;
//    End Modification to calculate R2      
      break;
    case WRITE_RESULTS:   // output the final probability density of the chain end-to-end distances
      FilePtr=fopen("results.dat","w");
      fprintf(FilePtr, "%s", "# R (sigma)              P(R)\n");
      dR=1.0/dR;
      Normalize=1.0/Normalize;
      for(i=0;i<MAX_NUMBER_OF_BINS;i++)
        if(histogram[i]>0.5)
          fprintf(FilePtr,"%14f %14f\n",i*dR+dR/2.,histogram[i]*Normalize/dR);
      fclose(FilePtr);
//    Start Modification to calculate R2
//      Normalize and print R2
      printf ("averageR2                   : %4.2f \n", averageR2*Normalize);
//    End Modification to calculate R2
      break;
  }
}

