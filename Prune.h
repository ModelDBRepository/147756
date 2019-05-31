
#ifndef __PRUNE__H_
#define __PRUNE__H_

#include"Sholl.h"
#include<sys/time.h>


#define BASAL 0
#define APICAL 1

struct Subset{
	int npts;		// Number of points in the subset
	int * Pts;		// The indices to the points
};

struct twodouble{
	double x1;
	double x2;
};

class Prune: public Sholl
{
protected:

	Subset TP;		// Terminal Points;
	Subset BP;		// Branch Points;
	Subset Stems;	// Stems;
	Subset Soma;	// Soma;

	Subset * Children;

	ShollData BPStats;	// BP Stats in sholl fashion
	
	int npc;              // number of SWC compartments in pruned neuron
	SWCData * PruCom ;     // The Pruned Compartments

	int antra;		// number of sholl tracks on apical side
	int bntra;		// number of sholl tracks on basal side

	double * abpzeta;	// Actual BP pruning values along apical side
	double * bbpzeta;	// Actual BP pruning values along basal side;

	double * adlzeta;	// Actual DL pruning values along apical side
	double * bdlzeta;	// Actual DL pruning values along basal side;
	

	double * abpP;			// Probability of apical BP Pruning
	double * adlP;			// Probability of apical DL Pruning

	double * bbpP;			// Probability of basal BP Pruning
	double * bdlP;			// Probability of basal DL Pruning

	double * bpPruneP;

	double maxdlP;

	double * abpprune;	// Completed BP pruning along apical side
	double * bbpprune;	// Completed BP pruning along basal side;

	double * adlprune;	// Completed DL pruning along apical side
	double * bdlprune;	// Completed DL pruning along basal side;

	char * outfilename;
	char * basefilename;

	int * TPBP;
	double * TPBPdist;

	int cmr;
	double dlr;
	double dlsum;
	double dlred;	// Specification on the amount of reduction


public:

	Prune();
	void LoadParams(char * = NULL);
	void PruneTree();
	void PrintFinal();
	void WriteBPStats(char*);

protected:

	void printSubsets();
	void printStatus();
	void generateSpecs();
	void saveSWC(char * = NULL);
	void checkTPs();

// ******************* Base Functions *****************************//

	void getBP();		// Get Branch Points
	void getTP();		// Get Terminal Points
	void getStems();	// Get Stems
	void getSoma();		// Get Soma

	int isParentSoma(int);	// To find if a pt's parent is soma
	int isParentBP(int);	// To find if a pt's parent is BP

	void findTPBP();		// Find the BP's corres. to TP's
	void findChildren();	// Update/form Children list
	void getBPHist();		// Per segment BP count
	void getStemHist();		// Per segment Stem count

	double disttoParent(int); 	// Distance of a cmpt to its parent     

// ******************* Specification Functions **********************//

	twodouble normal(double, double, double, double); // Normal distribution
	double pnorm(double);	// Gaussian CDF
	double dnorm(double);	// Gaussian PDF

	double ratio(double,double,double);	// PDF of a+x/b+y
	double maxratio (double,double,double, double);
	double maxxratio (double,double,double, double);
	double finratio (double,double,double, double);

	void setZeta(char *);
	void setMaxZeta(char *);
	void setApicalMaxZeta(char *);
	void setBasalMaxZeta(char *);
	double generateZeta(double, double, double, double, double);
	void setDLProbs();
	void setBPProbs();
	void setProbabilities();	// Set the initial probabilities

// ******************* TPPrune Functions *****************************//

	void pruneTP(int, int);
	void tpPrune(int);
	void rcParentBP(int, int);
	void rcParentSoma(int);
	void rcParentDendrite(int);
	void rmcmpt(int);
	void removeCompartment(int, int);

// ******************* BPPrune Functions *****************************//

	void pruneBP(int);
	void bpPrune(int);
	int removeCmprtmnt(int);

// ******************* Prune Functions *******************************//

	void initialize();
	void prune();	// Prune
	double convmeasure();
};	
#endif 
