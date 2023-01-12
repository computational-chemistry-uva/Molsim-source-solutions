#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

int main(void)
{
  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  // read data from disk
  ReadData();
 
  // Finally, Perform An Md Simulation        C
  MdLoop();

  return 0;
} 
