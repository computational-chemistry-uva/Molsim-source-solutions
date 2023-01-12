#include "system.h"

int NumberOfSteps;
int NumberOfInitializationSteps;
int NumberOfTrialPositions;
int ChainLength;
int Lstatic;
int Lcbmc;

VECTOR Positions[MAX_CHAIN_LENGTH];
VECTOR TrialPositions[MAX_CHAIN_LENGTH];

double Rcut;
double Kt;
double Theta0;
double A;
double Temperature;
double Beta;
