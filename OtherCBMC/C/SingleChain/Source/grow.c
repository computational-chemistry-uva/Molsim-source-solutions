#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

// Calculate angular potential energy  associated with 3 consecutive beads:  I-2 .....  I-1 ..... I
double Angle(VECTOR dr1,VECTOR dr2)
{
  double a,U;
  // dr1:          vector from I-1 to I
  // dr2:          vector from I-1 to I-2 
  // U:            angular potential energy
  // a:            cosine of the angle formed by the three beads

  //      Calculate the cosine a of the angle formed by the  vectors dr1 and dr2 by calculating their dot product: 
  //         dr1.dr2 = |dr1| |dr2| cos(dr1^dr2) <=> 
  //     <=> cos(dr1^dr2) = dr1.dr2/(|dr1|*|dr2|)

  a=(dr1.x*dr2.x+dr1.y*dr2.y+dr1.z*dr2.z)/
      sqrt((SQR(dr1.x)+SQR(dr1.y)+SQR(dr1.z))*
           (SQR(dr2.x)+SQR(dr2.y)+SQR(dr2.z)));

  // Calculate the potential energy associated with that angle
  //    Note: the function acos(a) returns the angle between 0 and Pi associated
  //    with a cosine value of R
  U=0.5*Kt*SQR(acos(a)-Theta0);
  return U;
}

// Calculate repulsive potential energy between two beads
double Repulsion(VECTOR dr)
{
//       Dr:       vector connecting the two beads
//       U:        repulsive energy between the two beads
//       R:        distance between two beads
  double R,U;

  R=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  U=0.0;
  // If the distance between the two beads is less than cutoff, calculate the repulsive energy.  
  //   Otherwise, the repulsive energy is 0.

  if(R<Rcut) U=A*SQR(R-Rcut)/SQR(Rcut);
  return U;
}

// grow a chain using cbmc or random insertion
void Grow(int LOldChain,double *Weight,double *Ubonded,double *Unonb)
{
  int LValidTrial;
  int i,j,k,Ntrial,Iwalk;
  VECTOR dr1,dr2,vec,Trial[MAX_TRIALS];
  double BendEnergy[MAX_TRIALS],ExternalEnergy[MAX_TRIALS],BoltzmannFactor[MAX_TRIALS];
  double Ws,Cumw,Cmmm;

  *Ubonded=0.0;
  *Unonb=0.0;
  *Weight=1.0;
  Ws=0.0;

  if(Lcbmc)
    Ntrial=NumberOfTrialPositions;
  else
    Ntrial=1;

  // the first bead of the trial chain is always at (0,0,0)
  TrialPositions[0].x=0.0;
  TrialPositions[0].y=0.0;
  TrialPositions[0].z=0.0;

  // Grow The Chain; Loop Over All Segments Exept The First One 
  // 
  // Positions            : Position Of The Old Chain 
  // TrialPositions       : Position Of The Trial Chain
  // Trial                : Position Of A Trial Segment
  // Lcbmc                : Do We Use Cbmc (.True. Or .False.)
  // LOldChain            : Do We Have The Old Configuration
  //                        (Only For Dynamic Schemes)

  // Loop over all the particles in the polymer chain
  for(i=1;i<ChainLength;i++)
  {
    // loop over all trial positions
    for(j=0;j<Ntrial;j++)
    { 
      do
      {
        BendEnergy[j]=0.0;
        ExternalEnergy[j]=0.0;
        LValidTrial=FALSE;
            
        // generate trial position
        // j=0 And LOldChain : old position
        // Otherwise    : new position
        // old coordinates are stored in Positions
        // the new chain is stored in TrialPositions
        // trial positions are stored in Trial
        // function RanSphere(&vec): creates a random unit
        // vector on a sphere centered at (0,0,0)
         
          // if an old chain exists, we must retrace it to obtain its
          // statistical weight; in that case the position of bead I in the old chain is also used as the first trial position of bead I

          // Note: retracing the old chain is not indispensable in our example
          // because we only have one chain; instead we could have just saved its statistical
          // weight from the previous step.    
          // However, in a real CBMC implementation with more than one chain, the
          // environment of the chain may have changed since the previous time its
          // statistical weight was calculated; in that case, retracing is necessary. To
          // keep the algorithm general, CBMC is always implemented with retracing of
          // the old chain.
 
          if(j==0&&LOldChain)
          {
            Trial[j].x=Positions[i].x;
            Trial[j].y=Positions[i].y;
            Trial[j].z=Positions[i].z;
          }
          // Otherwise, generate a new trial position for bead I
          else
          {
        // start modification
            RanSphere(&vec);
            Trial[j].x=TrialPositions[i-1].x+vec.x;
            Trial[j].y=TrialPositions[i-1].y+vec.y;
            Trial[j].z=TrialPositions[i-1].z+vec.z;
        // end modification
          }

        // calculate bond-bending potential
        if(i>1)
        {
          // Vector from bead I-1 to trial position J
          dr1.x=Trial[j].x-TrialPositions[i-1].x;
          dr1.y=Trial[j].y-TrialPositions[i-1].y;
          dr1.z=Trial[j].z-TrialPositions[i-1].z;
          // Vector from bead I-1 to bead I-2
          dr2.x=TrialPositions[i-2].x-TrialPositions[i-1].x;
          dr2.y=TrialPositions[i-2].y-TrialPositions[i-1].y;
          dr2.z=TrialPositions[i-2].z-TrialPositions[i-1].z;
          // Calculate angular potential energy  associated with 3 consecutive beads:  I-2 .....  I-1 ..... J
          BendEnergy[j]=Angle(dr1,dr2);
        }

        // accept or reject the chosen configuration
        // this is for cmbc only, because only in cbmc
        // a trial position according to a bond-bending
        // potential has to be generated... 
        // when j=0&&LOldChain : always accepted !! (why?)

        if(Lcbmc)
        {
            if(j==0&&LOldChain) // old chain is always accepted
              LValidTrial=TRUE;
            else
            {
          // start modification
          // Accept or reject the other trial positions according to the Metropolis algorithm
              if(RandomNumber()<exp(-Beta*BendEnergy[j]))
                LValidTrial=TRUE;
          // end modification
            }
        }
        else  // If CMBC is not used, always accept  the single trial position that was randomly generated.
          LValidTrial=TRUE;
      } while(!LValidTrial);

      // calculate the external energy 
      // only when there are at least 4 beads...
      if(i>2)
      {
        for(k=0;k<(i-2);k++)
        {
          dr1.x=Trial[j].x-TrialPositions[k].x;
          dr1.y=Trial[j].y-TrialPositions[k].y;
          dr1.z=Trial[j].z-TrialPositions[k].z;
          ExternalEnergy[j]+=Repulsion(dr1);
        }
      }
    }

    // end loop over all trial positions

    // choose a trial position for cbmc 
    // for cbmc and an old config.: 1st config.
    Ws=0.0;
    Iwalk=0;
    if(Lcbmc)
    {
      for(k=0;k<Ntrial;k++)
      {
        BoltzmannFactor[k]=exp(-Beta*ExternalEnergy[k]);
        Ws+=BoltzmannFactor[k];
      }
      Cumw=BoltzmannFactor[0];
      Cmmm=Ws*RandomNumber();
      if(!LOldChain)
      {
        while(Cumw<Cmmm)
        {
          Iwalk++;
          Cumw+=BoltzmannFactor[Iwalk];
        }
      }
    }
    else
      Ws=exp(-Beta*(ExternalEnergy[0]+BendEnergy[0]));

    // store the chosen configuration
    // update energies/Rosenbluth weight
    TrialPositions[i].x=Trial[Iwalk].x;
    TrialPositions[i].y=Trial[Iwalk].y;
    TrialPositions[i].z=Trial[Iwalk].z;

    (*Ubonded)+=BendEnergy[Iwalk];
    (*Unonb)+=ExternalEnergy[Iwalk];
    (*Weight)*=(Ws/Ntrial);
  } // END loop over all particles in chain
}
