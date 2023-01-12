// Generate multiple configurations of a polymer chain and calculate its averaged properties (energy, end-to-end distance)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Mcloop(void)
{
  int i,j,k;
  int NumberAccepted,NumberTrials;
  double RosenbluthFactor,RosenbluthFactorOld,r;
  double BondedEnergyOld,BondedEnergyNew,NonBondedEnergyNew,NonBondedEnergyOld;
  double BondedEnergySum,NonBondedEnergySum,EnergyNormalize,BondedEnergy,NonBondedEnergy;
  FILE *FilePtr;

  FilePtr=fopen("movie.pdb","w");
  InitializePdb(FilePtr);

//   NumberTrials:                               number of attempts to create a chain
//   NumberAccepted:                             number of accepted chains 
//   BondedEnergySum:                            average bonded energy
//   NonBondedEnergySum:                         average non bonded energy
//   EnergyNormalize:                            number of configurations used to compute averages 
//   BondedEnergyOld, BondedEnergyNew:           bonded energy of the old and trial chains
//   NonBondedEnergyOld, NonBondedEnergyNew:    non bonded energy of the old and trial chains
//   RosenbluthFactorOld, RosenbluthFactor:      statistical weight of the old and trial chains.  If random chains are generated, both values are 1.
//   r:                                          end-to-end distance
//   BondedEnergy,NonBondedEnergy:                Bonded and non bonded energy of the  chain 

  // initiallize variables
  RosenbluthFactor=1.0;
  BondedEnergyNew=0.0;
  NonBondedEnergyNew=0.0;
  RosenbluthFactorOld=1.0;
  BondedEnergyOld=0.0;
  NonBondedEnergyOld=0.0;
  NumberTrials=0.0;
  NumberAccepted=0.0;
  BondedEnergySum=0.0;
  NonBondedEnergySum=0.0;
  EnergyNormalize=0.0;
      
  // generate a first chain
  Grow(FALSE,&RosenbluthFactorOld,&BondedEnergyOld,&NonBondedEnergyOld);


  for(j=0;j<ChainLength;j++)
  {
    Positions[j].x=TrialPositions[j].x;
    Positions[j].y=TrialPositions[j].y;
    Positions[j].z=TrialPositions[j].z;
  }

  BondedEnergy=BondedEnergyOld;
  NonBondedEnergy=NonBondedEnergyOld;

  // Initialize the subroutine "Sample" (the subroutine that will sample the average end-to-end distance)
  Sample(INITIALIZE,1.0,1.0);

  // Generate the Markov chain of states
  for(i=0;i<NumberOfSteps;i++)
  {
    if(fmod(i,20) == 0) {
      printf("cycle: %d \n", i);
    }
    for(k=0;k<1000;k++)
    {
      NumberTrials++;
      // retrace old chain
      Grow(TRUE,&RosenbluthFactorOld,&BondedEnergyOld,&NonBondedEnergyOld);
      // create a new chain
      Grow(FALSE,&RosenbluthFactor,&BondedEnergyNew,&NonBondedEnergyNew);

      // If the new chain is accepted, then copy coordinates
      if(RandomNumber()<(RosenbluthFactor/RosenbluthFactorOld))
      {
        NumberAccepted++;
        NonBondedEnergy=NonBondedEnergyNew;
        BondedEnergy=BondedEnergyNew;

        for(j=0;j<ChainLength;j++)
        {
          Positions[j].x=TrialPositions[j].x;
          Positions[j].y=TrialPositions[j].y;
          Positions[j].z=TrialPositions[j].z;
        }
      }
      RosenbluthFactor=1.0;

      // If the simulation is beyond the set number of equilibration steps,
      // sample End-To-End Distance And Average Energies
      if(i>=NumberOfInitializationSteps) 
      {
        r=sqrt(SQR(Positions[ChainLength-1].x-Positions[0].x)+
              SQR(Positions[ChainLength-1].y-Positions[0].y)+
              SQR(Positions[ChainLength-1].z-Positions[0].z));
        Sample(SAMPLE,r,RosenbluthFactor);
        BondedEnergySum+=BondedEnergy;
        NonBondedEnergySum+=NonBondedEnergy;

        // update number of configurations sampled
        EnergyNormalize+=1.0;

        // Periodically write a pdb file containing the particle coordinates
        if(k==0&&(i%5)==0) 
          WritePdb(FilePtr);
      }
    }
  }

  // Call subroutine "Sample" to output the histogram of the end-to-end distance
  Sample(WRITE_RESULTS,1.0,1.0);

  // Output average energies and fraction of accepted chains
  printf("Fraction Accepted           : %f\n", (double)NumberAccepted/(double)NumberTrials ); 
  printf("Average bonded energy       : %f\n",BondedEnergySum/EnergyNormalize);
  printf("Average nonbonded energy    : %f\n",NonBondedEnergySum/EnergyNormalize);

  fclose(FilePtr);
}
