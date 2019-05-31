
#include"Prune.h"

/*****************************************************************/

void Prune :: initialize()
{
	int i;
	isRead();
	ShollAnalysis();

	npc=ncomp;
	PruCom=new SWCData[npc];
	for(i=0; i<npc; i++)
		PruCom[i]=Comptmt[i];


	delete(Comptmt);

	getSoma();
	getStems();
	getBP();
	getTP();



	Children=new Subset[npc] ;
	for(i=0; i<npc; i++){
		Children[i].npts=0;
		Children[i].Pts=new int [3];	// Maximum is trifurcation
	}	

	findChildren();

 	//generateSpecs(); exit(1);

	abpprune=new double[sdata.na]; bbpprune=new double[sdata.nb];
	adlprune=new double[sdata.na]; bdlprune=new double[sdata.nb];
	abpP    =new double[sdata.na]; bbpP    =new double[sdata.nb];
	adlP    =new double[sdata.na]; bdlP    =new double[sdata.nb];

	for(i=0; i<sdata.na; i++)		// Amount of initial pruning is zero
		adlprune[i]=abpprune[i]=0;
	for(i=0; i<sdata.nb; i++)
		bdlprune[i]=bbpprune[i]=0;

	cmr=0;	// Reduction in number of compartments
	dlr=0;	// Reduction in dl.

	setZeta(basefilename);
	//setBasalMaxZeta(basefilename);
	setProbabilities();

/* findTPBP requires probabilities to be initialized and hence is called
 * after setProbabilities() */

	TPBP=new int[TP.npts];
	bpPruneP=new double[TP.npts];
	findTPBP();
}

/*****************************************************************/

void Prune :: PruneTree()
{
	int i;
	initialize();
	prune();
//	cerr << "\nNo. of Soma Pts: " << Soma.npts ;
//	cerr << "\nNo. of Stems: " << Stems.npts ;
//	cerr << "\nNo. of BP's: " << BP.npts ;
//	cerr << "\nNo. of TP's: " << TP.npts << endl ;
	char  * filename=new char[35];
	sprintf(filename,"%s_final.swc",outfilename);
	cerr << "\nFinal save in " << filename << " after pruning " 
		<< dlr << " micron of length" << endl ;
	saveSWC(filename);
	
	//printStatus();
}

/*****************************************************************/

void Prune :: prune()
{
	int i;
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);
	int SEED=tv.tv_usec;
	srand48(SEED);

	int count =0;
	int scount = 1000; // Save for every 1000 um reduction in DL.
	char  * filename=new char[35];

	//while(convmeasure()>0.05){
	while(dlr<dlred){
	//while(dlr<dlsum){
	//while (count < 5000000){
		if(TP.npts==0){
			cerr << "\nNo more to Prune...\n";
			break;
		}	
		for(i=0; i<TP.npts; i++){
			bpPrune(i);
			tpPrune(i);
		}	
		//cerr << count << "\t" << cmr << "\t" << dlr << "\t" 
		//	 << convmeasure() << endl ;
		count ++;
		if(dlr>scount){
			sprintf(filename,"%s_%d.swc",outfilename,scount);
			cerr << "\nSaving file " << filename << " after pruning " 
				<< dlr << " micron of length" << endl ;
			saveSWC(filename);
			scount += 1000;
		}
	}	

	//cerr << "\nReduction in DL: " << dlr ;
	//cerr << "\nReduction in Cmptmt: " << cmr << endl ;
}

/*****************************************************************/
// Prunes TP[i]

void Prune :: tpPrune(int i)
{
	double dist;
	int l;
	int sn;	// Segment number
	int cn;		// Compartment number
	int PRUNE;

	/* BP has been pruned and the number of TP's has been reduced.
	 * So, the value of i will be greater than what is allowed. We
	 * just return back without bothering much as it does not mean
	 * anything */ 

	if(i>(TP.npts-1)){
		return;
	}	

	cn    = TP.Pts[i];
	dist  = 0.0;
	dist += (somax-PruCom[cn].x)*(somax-PruCom[cn].x);
	dist += (somay-PruCom[cn].y)*(somay-PruCom[cn].y);
	dist += (somaz-PruCom[cn].z)*(somaz-PruCom[cn].z);
	dist  = sqrt(dist);
	sn   = (int)(dist/(double)STLENGTH); // Segment where TP[i] belongs
	PRUNE=0;
	if (PruCom[cn].tc==4){	// Apical
		if(drand48()*maxdlP<adlP[sn]) PRUNE=1;
	}
	else if (PruCom[cn].tc==3){	// Basal
		if(drand48()*maxdlP<bdlP[sn]) PRUNE=1;
	}
	else{
		cerr << "\ntpP :: Terminal point is not dendrite!!!\n\n";
		checkTPs();
		cerr <<  i << endl << PruCom[cn] ;
		cerr <<  Children[cn].npts << endl  ;
		for(l=0; l<Children[cn].npts; l++){
			cerr << PruCom[Children[cn].Pts[l]];
		}
		exit(1);
	}	

	if(PRUNE){
		pruneTP(cn,sn);
		//setProbabilities();		// Update probability
		findChildren();			// Update children list
		findTPBP();				// Update TPBP list
	}	
}	

/*****************************************************************/
// Prunes TPBP[i]

void Prune :: bpPrune(int i)
{
	double dist;
	int sn;	// Segment number
	int cn;		// Compartment number
	int PRUNE;
	cn    = TPBP[i];
	
	if(cn==-2) return; // Case where the TP does not end on a BP.
	
	dist  = 0.0;
	dist += (somax-PruCom[cn].x)*(somax-PruCom[cn].x);
	dist += (somay-PruCom[cn].y)*(somay-PruCom[cn].y);
	dist += (somaz-PruCom[cn].z)*(somaz-PruCom[cn].z);
	dist  = sqrt(dist);
	sn   = (int)(dist/(double)STLENGTH); // Segment where TP[i] belongs
	PRUNE=0;
	if (PruCom[cn].tc==4){	// Apical
		//if(drand48()*maxdlP<(bpPruneP[i]*abpP[sn])) PRUNE=1;
		if(drand48()*maxdlP<bpPruneP[i]) PRUNE=1;
	}
	else if (PruCom[cn].tc==3){	// Basal
		//if(drand48()*maxdlP<(bpPruneP[i]*bbpP[sn])) PRUNE=1;
		if(drand48()*maxdlP<bpPruneP[i]) PRUNE=1;
	}
	else{
		cerr << "\nBranching point is not dendrite!!!\n\n";
		//cerr << j << " " << i << " "<< PruCom[cn] ;
		exit(1);
	}	
	if(PRUNE){
		pruneBP(i);
		//setProbabilities();		// Update probability
		findChildren();
		findTPBP();
	}	
}

/*****************************************************************/
// Convergence measure.

double Prune :: convmeasure()
{
	double cm=0.0;
	int i;
	
/* If the percentage reduction is negative it is not considered. Can be
 * reviewed at a later stage. TEST */

	for(i=0; i<sdata.na; i++){
		if((adlzeta[i]-adlprune[i])/adlzeta[i]>cm)
			cm=fabs((adlzeta[i]-adlprune[i])/adlzeta[i]);
		if((abpzeta[i]-abpprune[i])/abpzeta[i]>cm)
			cm=fabs((abpzeta[i]-abpprune[i])/abpzeta[i]);
	}		

	for(i=0; i<sdata.nb; i++){
		if((bdlzeta[i]-bdlprune[i])/bdlzeta[i]>cm)
			cm=fabs((bdlzeta[i]-bdlprune[i])/bdlzeta[i]);
		if((bbpzeta[i]-bbpprune[i])/bbpzeta[i]>cm)
			cm=fabs((bbpzeta[i]-bbpprune[i])/bbpzeta[i]);
	}		

	return cm;
}

/*****************************************************************/

void Prune :: saveSWC(char * filename)
{
	char fname[35];
	if(!filename) strcpy (fname,outfilename);
	else strcpy(fname,filename);

   	ofstream outfile (fname);
   	for(int i=0; i<npc; i++)
       	outfile << PruCom[i] ;
   	outfile.close();
}

/*****************************************************************/

void Prune :: checkTPs()
{
	int cn;
	for (int i=0; i<TP.npts; i++){
		cn=TP.Pts[i];
		if((PruCom[cn].tc != 3) && (PruCom[cn].tc != 4)){
    		cerr << "\ncheckTPs :: TP not dendrite\n";
			exit(1);
		}
	}
}

/*****************************************************************/

