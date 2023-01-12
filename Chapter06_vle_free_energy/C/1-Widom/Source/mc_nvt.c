#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"


// equation of state of the Lennard-Jones fluid
 
int main(void)
{
  int i,j,k;
  int NumberOfCycles;
  double PressureSum,PressureCount,Pressure;
  double ChemicalPotentialSum,ChemicalPotentialCount,Dummy1,Dummy2;
  double EnergySquaredSum,EnergySum,EnergyCount;
  VECTOR pos;
  FILE *FilePtrMovie,*FilePtr;

  printf("**************** Mc_Nvt ***************\n");

  // initialize system
  ReadInputData();

  FilePtr=fopen("results.dat","w");
  FilePtrMovie=fopen("movie.pdb","w");

  NumberOfCycles=0;
  PressureSum=0.0;
  PressureCount=0.0;
  ChemicalPotentialSum=0.0;
  ChemicalPotentialCount=0.0;
  EnergySquaredSum=0.0;
  EnergySum=0.0;
  EnergyCount=0;

  // total energy of the system
  EnergySystem();
  printf(" Total Energy Initial Configuration: %lf\n",TotalEnergy);
  printf(" Total Virial Initial Configuration: %lf\n",TotalVirial);

  RunningEnergy=TotalEnergy;
  RunningVirial=TotalVirial;

  // start MC-cycles

  for(i=1;i<=2;i++)
  {
    // i=1 equilibration
    // i=2 production

    if (i==EQUILIBRATION)
    {
      NumberOfCycles=NumberOfEquilibrationCycles;
      if (NumberOfCycles!=0) printf(" Start Equilibration\n");
    }
    else
    {
      if(NumberOfCycles!=0) printf(" Start Production\n");
      NumberOfCycles=NumberOfProductionCycles;
    }
        
    NumberOfAttempts=0;
    NumberOfAcceptedMoves=0;

    // intialize the subroutine that adjust the maximum displacement
    Adjust();

    for(j=0;j<NumberOfCycles;j++)
    {
      // attempt to displace a particle
      for(k=0;k<NumberOfDisplacementsPerCycle;k++)
        Mcmove();

      if(i==PRODUCTION)
      {
        // sample averages
        if((j%SamplingFrequency)==0)
        {
          EnergySquaredSum+=SQR(RunningEnergy);
          EnergySum+=RunningEnergy;
          EnergyCount+=1.0;

          Sample(j,RunningEnergy,RunningVirial,&Pressure,FilePtr);
          PressureSum+=Pressure;
          PressureCount+=1.0;

          // calculate the chemical potential:
          // do 10 trials
          // calculate the average of [exp(-Beta*Energy)]
          // you can use the subroutine EnergyParticle for this

          for(k=0;k<10;k++)
          {
            // start modification

            pos.x = RandomNumber()*Box;
            pos.y = RandomNumber()*Box;
            pos.z = RandomNumber()*Box;
            EnergyParticle(pos,NumberOfParticles,0,&Dummy1,&Dummy2);
            ChemicalPotentialSum += exp(-Beta*Dummy1);
            ChemicalPotentialCount += 1.0;

            // end modification
          }
        }
      }

      if((j%20)==0) WritePdb(FilePtrMovie);

      if((j%(NumberOfCycles/5))==0)
      {
        printf("======>> Done %d out of %d\n",j,NumberOfCycles);

        // write intermediate configuration to file
        Store();

        // adjust maximum displacements
        Adjust();
      }
    }
    if(NumberOfCycles!=0)
    {
      if (NumberOfAttempts!=0) 
      {
        printf("Number Of Att. To Displ. A Part.  : %d\n",NumberOfAttempts);
        printf("Success: %d (%lf)\n",NumberOfAcceptedMoves,100.0*NumberOfAcceptedMoves/NumberOfAttempts);
      }

      // test total energy
      EnergySystem();

      if(fabs(TotalEnergy-RunningEnergy)>1.0e-6)
        printf("######### Problems Energy ################\n");
      if(fabs(TotalVirial-RunningVirial)>1.0e-6)
        printf("######### Problems Virial ################\n");

      printf(" Total Energy End Of Simulation    : %lf\n",TotalEnergy);
      printf("       Running Energy              : %lf\n",RunningEnergy);
      printf("       Difference                  : %lf\n",TotalEnergy-RunningEnergy);
      printf(" Total Virial End Of Simulation    : %lf\n",TotalVirial);
      printf("       Running Virial              : %lf\n",RunningVirial);
      printf("       Difference                  : %lf\n",TotalVirial-RunningVirial);

      // print chemical potential and pressure
      if(i==PRODUCTION)
      {
        printf("Average Pressure                  : %lf\n",
          PressureSum/PressureCount);
        printf("Chemical Potential                : %lf\n",
          -log((ChemicalPotentialSum/ChemicalPotentialCount)*(CUBE(Box)/NumberOfParticles))/Beta);
        printf("Heat capacity                     : %lf\n",
          (EnergySquaredSum/EnergyCount)-SQR(EnergySum/EnergyCount));
      }
    }
  }
  Store();

  fclose(FilePtrMovie);
  fclose(FilePtr);

  return 0;
}
