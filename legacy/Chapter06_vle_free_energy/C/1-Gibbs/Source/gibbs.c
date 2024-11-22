#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system_2.h"
#include "parameter.h"
#include "ran_uniform.h"

double X[Npmax],Y[Npmax],Z[Npmax];
int Id[Npmax],Npart,Npbox[2];
double Box[2],Hbox[2],Temp,Beta;
double  Eps4,Eps48,Sig2,Mass,Rc[2],Rc2[2];
double Chp[2];
int Ichp[2]; 

int main()
{
  //    ---Gibbs-Ensemble Simulation Of The Lennard-Joned Fluid
  int  Equil, Prod, Nsamp, I, Icycl, Ndispl, Attempt; 
  int  Nacc, Ncycl, Nmoves, Imove, Nvol, AcceptVolume; 
  int  AttemptVolume, BoxID, Nswap, AcceptSwap, AttemptSwap;
  double En[2], Ent[2], Vir[2], Virt[2], Vmax, Dr; 
  double Ran, Succ;
  FILE* Fileptr;
  
  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));
  
  printf("**************** GIBBS ***************\n");
  
  //    ---Initialize System
  Readdat(&Equil, &Prod, &Nsamp, &Ndispl, &Dr, &Nvol, &Vmax, &Nswap, &Succ);
  Nmoves = Ndispl + Nvol + Nswap;
  
  //    ---Total Energy Of The System
  for(BoxID=0;BoxID<2;BoxID++)
    {
      Toterg(&(En[BoxID]),&(Vir[BoxID]),BoxID);
      printf("Box: %d\n",BoxID);
      printf("Total Energy Initial Configuration: %lf\n",En[BoxID]);
      printf("Total Virial Initial Configuration: %lf\n",Vir[BoxID]);
    }

  //    ---Start MC Cycle
  
  for(I=1;I<3;I++)
    {
      //       --- I=1 Equilibration
      //       --- I=2 Production
      
      if(I==1)
	{
	  Ncycl = Equil;
	  if(Ncycl!=0) printf("Start Equilibration\n"); 
	}
      else
	{
	  Ncycl = Prod;
	  if(Ncycl!=0) printf("Start Production\n"); 
	}
      Attempt = 0;
      Nacc = 0;
      AttemptVolume = 0;
      AcceptVolume = 0;
      AttemptSwap = 0;
      AcceptSwap = 0;
      
      //       ---Initialize Calculation Chemical Potential
      
      Init_Chem(0);
      
      //       ---Intialize The Subroutine That Adjust The Maximum Displacement

      Adjust(Attempt, Nacc, Dr, AttemptVolume, AcceptVolume, Vmax, Succ);
      
      for(Icycl=0;Icycl<Ncycl;Icycl++)
	{
	  for(Imove=0;Imove<Nmoves;Imove++)
            {
	      Ran = RandomNumber()*(Ndispl+Nvol+Nswap);
	      if(Ran<(double)(Ndispl)) 
		{ 
		  //                ---Attempt To Displace A Particle
                  Mcmove(En, Vir, &Attempt, &Nacc, Dr);
		}
	      else if (Ran<(double)(Ndispl+Nvol))
		{
		  //                ---Attempt To Change The Volume
                  Mcvol(En, Vir, &AttemptVolume,&AcceptVolume, Vmax);
		}
	      else
		{
		  //                ---Attemp To Exchange Particles
		  Mcswap(En, Vir, &AttemptSwap, &AcceptSwap);
		}
            }
	  if(I==2)
            {  
	      //             ---Sample Averages
	      if((Icycl%Nsamp)==0) Sample(Icycl, En, Vir);
            }
	  if((Icycl%(Ncycl/10))==0)
            {
	      printf("======>> Done %d out of %d %d %d\n", Icycl, Ncycl,
		     Npbox[0], Npbox[1]);
	      
	      //             ---Write Intermediate Configuration To File
	      Fileptr=fopen("fort.8","w");
	      Store(Fileptr, Dr, Vmax);
	      fclose(Fileptr);
	      for(BoxID=0;BoxID<2;BoxID++)
		{
		  printf("Density: %g\n",(double)Npbox[BoxID]/pow(Box[BoxID],3.));
		}
	      //             ---Adjust Maximum Displacements
	      
	      Adjust(Attempt, Nacc, Dr, AttemptVolume, AcceptVolume, Vmax, Succ);
            }
	}
      if(Ncycl!=0)
	{
	  if(Attempt!=0) 
	    {
	      printf("Number of Att. to Displ. a Part.: %d\n", Attempt);
	      printf("Success: %d (= %lf %%)\n",Nacc,
		     100.*(double)(Nacc)/(double)(Attempt));
            }
	  if(AttemptVolume!=0) 
            {
	      printf("Number of Att. to Change Volume: %d\n", AttemptVolume);
	      printf("Success: %d (= %lf %%)\n",AcceptVolume,
		     100.*(double)(AcceptVolume)/(double)(AttemptVolume));
            }
	  if(AttemptSwap!=0) 
            {
	      printf("Number of Att. to Exchange Part.: %d\n", AttemptSwap);
	      printf("Success: %d (= %lf %%)\n",AcceptSwap,
		     100.*(double)(AcceptSwap)/(double)(AttemptSwap));
            }
	  for(BoxID=0;BoxID<2;BoxID++)
            { 
	      //             ---Test Total Energy
	      Toterg(&(Ent[BoxID]), &(Virt[BoxID]), BoxID);
	      if(fabs(Ent[BoxID]-En[BoxID])>1e-6) 
		{
                  printf(" ######### Problems Energy ################ \n");
		}
	      if(fabs(Virt[BoxID]-Vir[BoxID])>1e-6)
		{
                  printf(" ######### Problems Virial ################ \n");
		}
	      printf("Box: %d\n",BoxID);
	      printf("Total Energy End of Simulation: %g\n",Ent[BoxID]);
	      printf("Running Energy: %g\n",En[BoxID]);
	      printf("Difference: %g\n",(Ent[BoxID] - En[BoxID]));
	      printf("Total Virial End of Simulation: %g\n",Virt[BoxID]);
	      printf("Running Virial: %g\n",Vir[BoxID]);
	      printf("Difference: %g\n",(Virt[BoxID] - Vir[BoxID]));
	      printf("Number of particles: %d\n",Npbox[BoxID]);
	      printf("Volume of box: %g\n",pow(Box[BoxID],3.));
	      printf("Density: %g\n",(double)Npbox[BoxID]/pow(Box[BoxID],3.));
            }
	  //          ---Calculation Chemical Potential
	  Init_Chem(2);
	}
    }
  Fileptr=fopen("fort.21","w");
  Store(Fileptr, Dr, Vmax);
  fclose(Fileptr);
  return(0);
}
