#include <stdio.h>

#define SQR(x) ((x)*(x))

#define FALSE 0
#define TRUE 1

enum{INITIALIZE,SAMPLE,WRITE_RESULTS};

extern int NumberOfSteps;
extern int NumberOfInitializationSteps;
extern int NumberOfTrialPositions;
extern int ChainLength;
extern int Lstatic;
extern int Lcbmc;

#define MAX_CHAIN_LENGTH 100
#define MAX_TRIALS 20

typedef struct
{
  double x;
  double y;
  double z;
} VECTOR;

extern VECTOR Positions[MAX_CHAIN_LENGTH];
extern VECTOR TrialPositions[MAX_CHAIN_LENGTH];
extern double Rcut;
extern double Kt;
extern double Theta0;
extern double A;
extern double Temperature;
extern double Beta;

void Mcloop(void);
void WritePdb(FILE *FilePtr);
void InitializePdb(FILE *FilePtr);
void Sample(int Choise,double R2,double W);
void Grow(int Lold,double *Weight,double *Ubonded,double *Unonb);
