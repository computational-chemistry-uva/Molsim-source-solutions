#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

// compute the Ewald sum of a cubic lattice 
// (for example Nacl)

int main(void)
{
  double UReal,UFourier,USelf;
  FILE *FilePtr;

  InitializeRandomNumberGenerator(time(0l));

  FilePtr=fopen("input","r");
  fscanf(FilePtr,"%d %lf %d\n",&NumberOfCells,&Alpha,&NumberOfWavevectors);
  fclose(FilePtr);

  if(NumberOfWavevectors>MAX_WAVEVECTOR)
    {
      printf("%d is larger than the allowed number of wavevectors %d\n",NumberOfWavevectors,MAX_WAVEVECTOR);
      return -1;
    }

  Lattice();
  Realspace(&UReal);
  Fourierspace(&UFourier,&USelf);

  UReal/=NumberOfParticles;
  UFourier/=NumberOfParticles;
  USelf/=NumberOfParticles;

  printf("Real space energy          : %18.10f\n",UReal);
  printf("Fourier space energy       : %18.10f\n",UFourier);
  printf("Self interaction energy    : %18.10f\n",USelf);
  printf("Total Coulomb energy       : %18.10f\n",UReal+UFourier+USelf);
  printf("\n");
// The energy is calculated as a sum over all particle pairs. Since the
// energy needs to be added to the energy of each particle in the pair, 
// we need a factor two: 
  printf("Madelung constant          : %18.10f\n",-2.0*(UReal+UFourier+USelf));

  return 0;
}
