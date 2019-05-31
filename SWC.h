
#ifndef __SWC__H_
#define __SWC__H_

#include<math.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<fstream.h>
#include<iomanip.h>

struct SWCData{
	double pi;
    int tc;
    double x;
    double y;
    double z;
    double r;
    int ppi;
};

ostream & operator << (ostream & tmp, SWCData swcd);
istream & operator >> (istream & tmp, SWCData & swcd);

class SWC
{
protected:
	int ncomp;				// number of SWC compartments in swc file
	SWCData * Comptmt ; 	// The actual SWC compartments

public:

	SWC();
	void ReadSWC(char * filename=NULL);
	void WriteSWC(char * filename=NULL);
	void SetComptmt(SWCData * Cpts, int npts);

protected:

	void isRead();
};	
#endif 
