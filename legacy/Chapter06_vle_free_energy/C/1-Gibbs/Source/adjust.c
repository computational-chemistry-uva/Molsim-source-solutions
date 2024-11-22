#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

static int Naccp,Attempp,AttemptVolumep,AcceptVolumep;

int Adjust(int Attemp, int Nacc, double Dr, int AttemptVolume, int AcceptVolume, double Vmax, double Succ)
{
  // Sets Maximum Displacement And Maximum Volume Change
  // Such That 50% Of The Move Will Be Accepted
  
  double Dro, Frac, Vmaxo;
  
  // ---Displacement:
  if(Attemp==0||Attempp>=Attemp)
    { 
      Naccp = Nacc;
      Attempp = Attemp;
    }
  else
    {
      Frac = (double)(Nacc-Naccp)/(double)(Attemp-Attempp);
      Dro = Dr;
      Dr = Dr*fabs(Frac/(Succ/100.0));
      
      //    ---Limit The Change:
      
      if(Dr/Dro>1.5) Dr = Dro*1.5;
      if(Dr/Dro<0.5) Dr = Dro*0.5;
      if(Dr>Hbox[1]/2.0) Dr = Hbox[1]/2.0;
      
      printf("Max. Displ. set to: %g (Old: %g)\n", Dr, Dro);
      printf("Frac. Acc.: %g Attempts: %d Success: %d\n",
	     Frac, (Attemp - Attempp), (Nacc - Naccp));
      
      //    ---Store Nacc And Attemp For Next Use
      
      Naccp = Nacc;
      Attempp = Attemp;
    }
  
  // ---Volume:
  if(AttemptVolume==0||AttemptVolumep>=AttemptVolume)
    {
      AcceptVolumep = AcceptVolume;
      AttemptVolumep = AttemptVolume;
    }
  else
    {
      Frac = (double)(AcceptVolume-AcceptVolumep)/
	(double)(AttemptVolume-AttemptVolumep);
      Vmaxo = Vmax;
      Vmax = Vmax*fabs(Frac/(Succ/100.0));
      
      //    ---Limit The Change:
      
      if(Vmax/Vmaxo>1.5) Vmax = Vmaxo*1.5;
      if(Vmax/Vmaxo<0.5) Vmax = Vmaxo*0.5;
      printf("Max. Vol. Change set to: %g (Old: %g)\n", Vmax, Vmaxo);
      printf("Frac. Acc.: %g Attempts: %d Success: %d\n",
	     Frac, (AttemptVolume - AttemptVolumep), 
	     (AcceptVolume - AcceptVolumep));
      
      //    ---Store Nacc And Attemp For Next Use

      AcceptVolumep = AcceptVolume;
      AttemptVolumep = AttemptVolume;
    }
  return 0;
}
