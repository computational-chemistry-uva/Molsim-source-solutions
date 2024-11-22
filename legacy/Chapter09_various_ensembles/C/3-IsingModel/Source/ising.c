#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

#define SQR(x) ((x)*(x))
#define CycleMultiplication 100

// 2d Ising Model
// No Periodic boundary conditions
// 4 Neighbors 

#define MAX_LATTICE_SIZE 32
#define Maxmag 1024
#define MAX_MAGNETIZATION 1024

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

int main(void)
{
  int i,j,k,l;
  int NumberOfCycles,NumberOfInitializationSteps;
  int Magnetization,LatticeSize;
  int Inew,Jnew,Lnew,Lold,Diff,Ilat,Jlat;
  int Mnew,Mold;
  double Energy,Beta,Av1,Dist2,Weight,Move1,Move2,Dist3,Av2;
  double W[MAX_MAGNETIZATION+1],Dist1[2*MAX_MAGNETIZATION+1];
  int Lattice[MAX_LATTICE_SIZE+2][MAX_LATTICE_SIZE+2];
  int Iup[4]={1,-1,0,0};
  int Jup[4]={0,0,1,-1};
  FILE *FilePtr;

  InitializeRandomNumberGenerator(time(0l)) ;
 
  // read info 
  printf("2d Ising Model, 4 Neightbors, No Pbc\n\n");

  printf("Lattice Size                     ? ");
  fscanf(stdin,"%d",&LatticeSize);

  if(LatticeSize%2!=0||MAX_MAGNETIZATION%2!=0||MAX_LATTICE_SIZE%2!=0)
  {
    printf("error in input parameters\n");
    exit(1);
  }
      
  printf("Beta                             ? ");
  fscanf(stdin,"%lf",&Beta);

  printf("Number Of Cycles                 ? ");
  fscanf(stdin,"%d",&NumberOfCycles);

  if(LatticeSize<3||LatticeSize>MAX_LATTICE_SIZE||NumberOfCycles<10)
  {
    printf("error in input parameters\n");
    exit(1);
  }

  for(i=0;i<MAX_MAGNETIZATION;i++)
    W[i]=1.0;
      
  // multicanonical algorithm 
  FilePtr=fopen("w.dat","r");
  for(i=0;i<MAX_MAGNETIZATION;i+=2)
  {
    fscanf(FilePtr,"%d %lf\n",&j,&Av1);
    if(abs(j)>MAX_MAGNETIZATION||j%2!=0)
    {
      printf("j out of range\n");
      exit(1);
    }
    W[j]=Av1;
  }
  fclose(FilePtr);
      
  NumberOfInitializationSteps=NumberOfCycles/3;
  Av1=0.0;
  Av2=0.0;
  Move1=0.0;
  Move2=0.0;
      
  if(NumberOfInitializationSteps>100) NumberOfInitializationSteps=100;

  for(i=0;i<2*MAX_MAGNETIZATION+1;i++)
    Dist1[i]=0.0;

  // Initialize Lattice
  // Random Lattice
  // 
  // add one boundary layer with spin=0
  // this is a way to avoid if-statements
  // in the computation of the energy

  for(i=0;i<=MAX_LATTICE_SIZE+1;i++)
    for(j=0;j<=MAX_LATTICE_SIZE+1;j++)
      Lattice[i][j]=0;

  for(i=1;i<=LatticeSize;i++)
    for(j=1;j<=LatticeSize;j++)
    {
      if(RandomNumber()<0.5)
        Lattice[i][j]=1;
      else
        Lattice[i][j]=-1;
    }

  // calculate initial energy
  Energy=0.0;
  Magnetization=0.0;

  for(i=1;i<=LatticeSize;i++)
    for(j=1;j<=LatticeSize;j++)
    {
      Magnetization+=Lattice[i][j];

      for(k=0;k<4;k++)
      {
        Inew=i+Iup[k];
        Jnew=j+Jup[k];
        Energy-=0.5*Lattice[i][j]*Lattice[Inew][Jnew];
      }
    }

    printf("\n");
    printf("Lattice Size               : %d\n",LatticeSize);
    printf("Beta                       : %f\n",Beta);
    printf("Number Of Cycles           : %d\n",NumberOfCycles);
    printf("Number Of Init             : %d\n",NumberOfInitializationSteps);
    printf("\n");
    printf("\n");
    printf("Initial Energy             : %f\n",Energy);
    printf("Initial Magnetization      : %d\n",Magnetization);
    printf("\n");
      
    if(fabs(Magnetization)>MAX_MAGNETIZATION)
    {
      printf("|Magnetization|>MAX_MAGNETIZATION\n");
      exit(0);
    }

  // loop over all cycles
  for(l=0;l<NumberOfCycles;l++)
  {
    for(k=0;k<CycleMultiplication*SQR(LatticeSize);k++)
    {
      // flip a single spin
      // Metropolis algorithm
      Ilat=1+(int)(RandomNumber()*LatticeSize);
      Jlat=1+(int)(RandomNumber()*LatticeSize);
            
      Mold=Magnetization;
      Lold=Lattice[Ilat][Jlat];
      Lnew=-Lold;
      Mnew=Mold+Lnew-Lold;
      Diff=0;
      Move2=Move2+1.0;
            
      // calculate the energy difference
      // between new and old

      // start modification
        for(i=0;i<4;i++)
        {
          Inew=Ilat+Iup[i];
          Jnew=Jlat+Jup[i];
          Diff=Diff+Lattice[Inew][Jnew]*(Lold-Lnew);
        }
      // end modification

      // acceptance/rejection rule
      // use weight function W 
      if(RandomNumber()<exp(-Beta*Diff+W[abs(Mnew)]-W[abs(Mold)]))
      {
        // update the lattice/energy/magnetisation
        Move1=Move1+1.0;

        // start modification
          Lattice[Ilat][Jlat]=Lnew;
          Energy=Energy+Diff;
          Magnetization=Mnew;
        // end modification
               
        if(fabs(Magnetization)>MAX_MAGNETIZATION)
        {
          printf("Increase the value of maxmag\n");
          exit(0);
        }
      }
                  
      if(l>NumberOfInitializationSteps)
      {
        Weight=exp(-W[abs(Magnetization)]);
        Av1=Av1+Weight*Energy;
        Av2=Av2+Weight;
        Dist1[Magnetization+MAX_MAGNETIZATION]+=1.0;
      }
    }
  }

  printf("\n");
  printf("Average Energy             : %f\n",Av1/Av2);
  printf("Fraction Accepted Swaps    : %f\n",Move1/Move2);
  printf("\n");
  printf("Final Energy        (Simu) : %f\n",Energy);
  printf("Final Magnetization (Simu) : %d\n",Magnetization);
  printf("\n");

  // calculate distributions magnetization
  Dist2=0.0;
  Dist3=0.0;
            
  for(i=-MAX_MAGNETIZATION;i<=MAX_MAGNETIZATION;i+=2)
  {
    Weight=exp(-W[abs(i)]);
    Dist2+=Dist1[i+MAX_MAGNETIZATION];
    Dist3+=Dist1[i+MAX_MAGNETIZATION]*Weight;
  }

  FilePtr=fopen("magnetic.dat","w");
  for(i=-MAX_MAGNETIZATION;i<=MAX_MAGNETIZATION;i+=2)
  {
    if(Dist1[i+MAX_MAGNETIZATION]>0.5)
    {
      Weight=exp(-W[abs(i)]);
      fprintf(FilePtr,"%d %f %f\n",i,Dist1[i+MAX_MAGNETIZATION]/Dist2,Dist1[i+MAX_MAGNETIZATION]*Weight/Dist3);
    }
  }
  fclose(FilePtr);

  // weight distribution
  Dist2=0.0;
  for(i=0;i<=MAX_MAGNETIZATION;i+=2)
  {
    Dist1[i+MAX_MAGNETIZATION]=-log(MAX(1.0,Dist1[i+MAX_MAGNETIZATION]+Dist1[-i+MAX_MAGNETIZATION])*exp(-W[abs(i)]));
    Dist2=MIN(Dist2,Dist1[i+MAX_MAGNETIZATION]);
  }

  for(i=0;i<=MAX_MAGNETIZATION;i+=2)
  {
    Dist1[i+MAX_MAGNETIZATION]=Dist1[i+MAX_MAGNETIZATION]-Dist2;
    Dist1[i+MAX_MAGNETIZATION]=MIN(10.0,1.05*Dist1[i+MAX_MAGNETIZATION]);
  }

  FilePtr=fopen("w.dat","w");
  for(i=0;i<=MAX_MAGNETIZATION;i+=2)
    fprintf(FilePtr,"%d %f\n",i,Dist1[i+MAX_MAGNETIZATION]);
  fclose(FilePtr);
      
  // calculate final energy
  // used to check the code
  Energy=0;
  Magnetization=0;
  for(i=1;i<=LatticeSize;i++)
    for(j=1;j<=LatticeSize;j++)
    {
      Magnetization+=Lattice[i][j];
      for(k=0;k<4;k++)
      {
        Inew=i+Iup[k];
        Jnew=j+Jup[k];
        Energy-=0.5*Lattice[i][j]*Lattice[Inew][Jnew];
      }
    }

  printf("Final Energy        (Calc) : %f\n",Energy);
  printf("Final Magnetization (Calc) : %d\n",Magnetization);

  return 0;
}
