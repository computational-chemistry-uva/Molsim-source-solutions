extern int NumberOfCycles;
extern int NumberOfSteps;
extern double Tstep;
extern double Temperature;
extern double Xpos;
extern double Vpos;
extern double Theta;
extern double Qstar;

void Force(double X, double *U, double *F);
void Sample(int Switch, int Ttt);
void ReadData(void);
void MdLoop(void);
void Integrate(double *Cnn);
