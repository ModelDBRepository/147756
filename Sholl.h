#ifndef __SHOLL__H_
#define __SHOLL__H_

#include"SWC.h"

#define STLENGTH 50.0  // Sholl TRACK width For hippocampal CA3 neurons

#define MAXTRACK 20	// Maximum possible Sholl tracks

//#define RESLN 0.04237288135593220338 
#define RESLN 1
	// How many um does each point in morphology correspond to

#define MAXPTS 4500
	// How many points can a typical swc segment have with the above resolution

 
 
struct LineStruct{
	int num;			// Number of points in the line
	double * x;			// X coordinates of the points
	double * y;			// Y coordinates of the points
	double * z;			// Z coordinates of the points
};

struct ShollData {
    int nb;   // number of basal tracks
    int na;   // number of apical tracks
    double * basalcount ;    // Count in basal tracks
    double * apicalcount ;    // Count in apical tracks
};


class Sholl : public SWC
{
protected:
	double somax;	// CoG of Soma (x,y,z)
	double somay;
	double somaz;
	double totlen; // Total dendritic length
	ShollData sdata;
	LineStruct * ls;

public:
	Sholl();
	ShollData & ShollAnalysis ();
	void WriteSholl(char * basefilename=NULL);     

protected:
	void createLine (SWCData, SWCData);
	void isDone ();			// To see if Sholl's analysis has been done.
	void ComputeTotalLength();
};	

#endif 
