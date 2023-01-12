#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..

int nr_frame=1;
char RecordName[7]="ATOM  ";   // Atom
int SerialNumber;              // Atom serial number
char AtomName[4]="  C";        // Atom Name
char AltLoc[2]=" ";            // Alternate location indicator
char ResIdueName[4]="   ";     // ResIdue Name
char ChainId[2]=" ";           // Chain Identifier
int ResSeq=0;                  // ResIdue sequence number
char iCode[2]=" ";             // code for insertion of resIdues
double Occupancy=0.0;          // Occupancy
double Temp=533.0;             // Temperature factor
char SegID[4]="   ";           // Segment Identifier, left-justified
char Element[3]="  ";          // Element symbol, right-justified
char Charge[3]="  ";           // Charge on the Atom


void InitializePdb(FILE *FilePtr)
{
  int j;

  SerialNumber=0;
  fprintf(FilePtr,"MODEL%9d\n",nr_frame++);
  ChainId[0]='A';
  for(j=0;j<ChainLength;j++)
  {
    SerialNumber++;
    fprintf(FilePtr,"%s %4d %s %s %s %s %3d %s %8.3f%8.3f%8.3f%6.2f%6.2f %s %s %s\n",
        RecordName,SerialNumber,AtomName,AltLoc,ResIdueName,
        ChainId,ResSeq,iCode,1.5*j,0.0,0.0,Occupancy,Temp,SegID,Element,Charge);
  }
  fprintf(FilePtr,"ENDMDL\n");
}


void WritePdb(FILE *FilePtr)
{
  int j;

  SerialNumber=0;
  fprintf(FilePtr,"MODEL%9d\n",nr_frame++);
  ChainId[0]='A';
  
  for(j=0;j<ChainLength;j++)
  {
    SerialNumber++;
    fprintf(FilePtr,"%s %4d %s %s %s %s %3d %s %8.3f%8.3f%8.3f%6.2f%6.2f %s %s %s\n",
        RecordName,SerialNumber,AtomName,AltLoc,ResIdueName,
        ChainId,ResSeq,iCode,Positions[j].x,Positions[j].y,Positions[j].z,Occupancy,Temp,SegID,Element,Charge);
  }
  fprintf(FilePtr,"ENDMDL\n");
}
