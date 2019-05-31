
#include"Prune.h"


/*****************************************************************/

Prune :: Prune () :Sholl()
{
	TP.Pts=NULL;
	BP.Pts=NULL;
	Soma.Pts=NULL;
	Stems.Pts=NULL;
	Soma.npts=0;
	TP.npts=0;
	BP.npts=0;
	Stems.npts=0;
}

/*****************************************************************/

void Prune :: LoadParams(char * filename)
{
    char * flname=new char [35];
    if(!filename){
        cout << "\nGive input PRN filename : ";
        cin >> flname ;
    }
    else strcpy(flname,filename);
 
    ifstream prnfile (flname);
	basefilename=new char[35];
	outfilename=new char[35];

    if(!prnfile){
        cerr << "\nFile " << flname << " does not exist.\n";
        exit(1);
    }

	prnfile >> flname;	// SWC filename.
	ReadSWC(flname);
	prnfile >> basefilename; 	
	prnfile >> antra ;
	prnfile >> bntra ;
	prnfile >> dlred ;
	prnfile >> outfilename; // File to output SWC.
	prnfile.close();
	//delete(flname);
}
/*****************************************************************/
// Dendritic points which are not parents to 
// anything else are terminal points
/*****************************************************************/

void Prune :: getTP ()
{
	isRead();

	int i;
	int * temp=new int[ncomp];

	for(i=0; i<ncomp; i++)
		temp[i]=0;

	for(i=0; i<ncomp; i++)		// Finds all the points which are parents
		temp[PruCom[i].ppi-1]++;

	TP.npts=0;
	for(i=0; i<ncomp; i++)		// If it was not a parent and is dendrite
		if((!temp[i]) && ((PruCom[i].tc==3)||(PruCom[i].tc==4)))
				TP.npts ++;

	if(TP.Pts) delete TP.Pts;
	TP.Pts=new int[TP.npts];

	TP.npts=0;
	for(i=0; i<ncomp; i++)
		if((!temp[i]) && ((PruCom[i].tc==3)||(PruCom[i].tc==4)))
				TP.Pts[TP.npts++]=i;
	
	//cerr << "No of TP's: " << TP.npts << endl ;
	delete(temp);
}	

/*****************************************************************/
// If a segment has the soma as the parent and is a dendrite,
//  it is representative of a stem
/*****************************************************************/

void Prune :: getStems ()
{
	isRead();

	int i;

	Stems.npts=0;
	for(i=0; i<ncomp; i++)
		if((isParentSoma(i)) && (PruCom[i].tc==3 || PruCom[i].tc==4))
			Stems.npts++;

	if(Stems.Pts) delete Stems.Pts;
	Stems.Pts=new int[Stems.npts];

	Stems.npts=0;
	for(i=0; i<ncomp; i++)
		if((isParentSoma(i)) && (PruCom[i].tc==3 || PruCom[i].tc==4))
			Stems.Pts[Stems.npts++]=i;

	cerr << "No of Stems: " << Stems.npts << endl ;;
	getStemHist();
}	

/*****************************************************************/

void Prune :: getSoma()
{
	isRead();

	int i;
	Soma.npts=0;
	for(i=0; i<ncomp; i++)
		if(PruCom[i].tc==1) Soma.npts++;

	if(Soma.Pts) delete Soma.Pts;
	Soma.Pts=new int[Soma.npts];

	Soma.npts=0;
	for(i=0; i<ncomp; i++)
		if(PruCom[i].tc==1)
			Soma.Pts[Soma.npts++]=i;

	cerr << "No of Soma pts: " << Soma.npts << endl ;;
}			

/*****************************************************************/
// Dendritic points which have more than one children are branch 
// points. We are not bothered about soma.
/*****************************************************************/

void Prune :: getBP ()
{
	isRead();
	int i;
	int * temp=new int[ncomp];

	for(i=0; i<ncomp; i++)
		temp[i]=0;

	for(i=0; i<ncomp; i++)		// Finds all the points which are parents
		temp[PruCom[i].ppi-1]++;

	BP.npts=0;

// If a point is parent to more than one child and is a dendrite

	for(i=0; i<ncomp; i++)		
		if(temp[i]>1 && (PruCom[i].tc==3 || PruCom[i].tc==4)) 
			BP.npts ++;

	if(BP.Pts) delete BP.Pts;
	BP.Pts=new int[BP.npts];

	BP.npts=0;
	for(i=0; i<ncomp; i++){
		if(temp[i]>1 && (PruCom[i].tc==3 || PruCom[i].tc==4)){
			BP.Pts[BP.npts++]=i;
			if(temp[i]==3) 
				cerr << "Warning: Trifurcation at " << PruCom[i] ;
		}		
	}
	
	getBPHist();

	delete(temp);
}	

/*****************************************************************/
// To find whether the parent of a given pt is the soma
// If the parent is any of the segments part of the soma, then 
// return is 1
/*****************************************************************/

int Prune :: isParentSoma (int in)
{
	int ret=0;	
	int i;
	
	if(!Soma.npts) getSoma();	//Need to know which pts are soma

	for(i=0; i<Soma.npts; i++){
		if (PruCom[in].ppi==(Soma.Pts[i]+1)){
			ret=1; 
			break;
		}
	}	
	return ret;
}	

/*****************************************************************/

int Prune :: isParentBP (int in)
{
	int ret=0;	
	int i;

	if(!BP.npts) getBP();	//Need to know which pts are soma

	for(i=0; i<BP.npts; i++)
		if (PruCom[in].ppi==(BP.Pts[i]+1)){
			ret=1; 
			break;
		}
	return ret;
}	

/*****************************************************************/

void Prune :: getBPHist()
{
	int i;
	int cn,sn;
	double dist;

	BPStats.na=sdata.na;
	BPStats.nb=sdata.nb;

	BPStats.apicalcount=new double [BPStats.na];
	BPStats.basalcount=new double [BPStats.nb];

	for(i=0;i <BPStats.na; i++)
		BPStats.apicalcount[i]=0.0;

	for(i=0;i <BPStats.nb; i++)
		BPStats.basalcount[i]=0.0;

	for(i=0; i<BP.npts; i++){
		cn    = BP.Pts[i];
		dist  = 0.0;
		dist += (somax-PruCom[cn].x)*(somax-PruCom[cn].x);
		dist += (somay-PruCom[cn].y)*(somay-PruCom[cn].y);
		dist += (somaz-PruCom[cn].z)*(somaz-PruCom[cn].z);
		dist  = sqrt(dist);
		sn    = (int)(dist/(double)STLENGTH); // Segment where BP[i] belongs
		if (PruCom[cn].tc==4)	// Apical
			BPStats.apicalcount[sn]++;
		if (PruCom[cn].tc==3)	// Basal
			BPStats.basalcount[sn]++;
	}	
}

/*****************************************************************/
// Find distance between a compartment and its parent.

double Prune :: disttoParent(int in)
{
	double dist;
	int cn=PruCom[in].ppi-1;

	dist  = 0.0;
	dist += (PruCom[in].x-PruCom[cn].x)*(PruCom[in].x-PruCom[cn].x);
	dist += (PruCom[in].y-PruCom[cn].y)*(PruCom[in].y-PruCom[cn].y);
	dist += (PruCom[in].z-PruCom[cn].z)*(PruCom[in].z-PruCom[cn].z);
	dist  = sqrt(dist);
	
	return dist;
}

/*****************************************************************/

void Prune :: findChildren ()
{
	int tn;
	int i;

	for(i=0; i<npc; i++)	
		Children[i].npts=0;

	for(i=1; i<npc; i++){	// 1st point doesnot have a valid parent
		tn=PruCom[i].ppi-1;
		Children[tn].Pts[Children[tn].npts++]=i;	
	}	
}			

/*****************************************************************/

void Prune :: WriteBPStats(char * flname)
{
	char * filename=new char [35];
	sprintf(filename,"%s.asl",flname);
	ofstream aoutfile (filename,ios::app);
	int i;
	npc=ncomp;
	PruCom=new SWCData[npc];
	for(i=0; i<npc; i++)
		PruCom[i]=Comptmt[i];
	getBP();
	
	aoutfile << endl << endl ;
	for(i=0; i<BPStats.na; i++)
		aoutfile << BPStats.apicalcount[i] << endl ;

	sprintf(filename,"%s.bsl",flname);
	ofstream boutfile (filename,ios::app);

	boutfile << endl << endl ;
	for(i=0; i<BPStats.nb; i++)
		boutfile << BPStats.basalcount[i] << endl ;

	cerr << "Total number of BP is: " << BP.npts << endl ;

}	
/*****************************************************************/

void Prune :: getStemHist()
{
	int i;
	int cn,sn;
	double dist;

	ShollData StemsStats;	// Stems Stats in sholl fashion
	StemsStats.na=sdata.na;
	StemsStats.nb=sdata.nb;

	StemsStats.apicalcount=new double [StemsStats.na];
	StemsStats.basalcount=new double [StemsStats.nb];

	for(i=0;i <StemsStats.na; i++)
		StemsStats.apicalcount[i]=0.0;

	for(i=0;i <StemsStats.nb; i++)
		StemsStats.basalcount[i]=0.0;

	for(i=0; i<Stems.npts; i++){
		cn    = Stems.Pts[i];
		dist  = 0.0;
		dist += (somax-PruCom[cn].x)*(somax-PruCom[cn].x);
		dist += (somay-PruCom[cn].y)*(somay-PruCom[cn].y);
		dist += (somaz-PruCom[cn].z)*(somaz-PruCom[cn].z);
		dist  = sqrt(dist);
		sn    = (int)(dist/(double)STLENGTH); // Segment where Stems[i] belongs
		if (PruCom[cn].tc==4)	// Apical
			StemsStats.apicalcount[sn]++;
		if (PruCom[cn].tc==3)	// Basal
			StemsStats.basalcount[sn]++;
	}	

/*
	cerr << "\n\nSTEM STATISTICS\n\n";

	cerr << "APICAL\n\n";
	for(i=0; i<StemsStats.na; i++)
		cerr  << StemsStats.apicalcount[i] << endl ;

	cerr << "\n\nBASAL\n\n";
	for(i=0; i<StemsStats.nb; i++)
		cerr  << StemsStats.basalcount[i] << endl ;
*/		
}

/*****************************************************************/
