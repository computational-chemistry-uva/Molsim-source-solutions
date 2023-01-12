#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

// attempts to displace a randomly selected particle
void Mcmove(void)
{
  double EnergyNew,VirialNew,EnergyOld,VirialOld;
  VECTOR NewPosition;
  int i;
 
  NumberOfAttempts++;

  // choose a random particle
  i=NumberOfParticles*RandomNumber();

  // calculate old energy
  EnergyParticle(Positions[i],i,0,&EnergyOld,&VirialOld);

  // give a random displacement
  NewPosition.x=Positions[i].x+(RandomNumber()-0.5)*MaximumDisplacement;
  NewPosition.y=Positions[i].y+(RandomNumber()-0.5)*MaximumDisplacement;
  NewPosition.z=Positions[i].z+(RandomNumber()-0.5)*MaximumDisplacement;

  // calculate new energy
  EnergyParticle(NewPosition,i,0,&EnergyNew,&VirialNew);
         
  if(RandomNumber()<exp(-Beta*(EnergyNew-EnergyOld)))
  {
    // accept
    NumberOfAcceptedMoves++;
    RunningEnergy+=(EnergyNew-EnergyOld);
    RunningVirial+=(VirialNew-VirialOld);

    // put particle in simulation box
    if(Positions[i].x<0.0)
      Positions[i].x+=Box;
    else if(Positions[i].x>=Box)
      Positions[i].x-=Box;
            
    if(Positions[i].y<0.0)
      Positions[i].y+=Box;
    else if(Positions[i].y>=Box)
      Positions[i].y-=Box;
          
    if(Positions[i].z<0.0)
      Positions[i].z+=Box;
    else if(Positions[i].z>=Box)
      Positions[i].z-=Box;

   // update new position
   Positions[i].x=NewPosition.x;
   Positions[i].y=NewPosition.y;
   Positions[i].z=NewPosition.z;
  }
}
