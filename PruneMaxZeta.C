#include"Prune.h"

/*****************************************************************/

void Prune :: setMaxZeta(char * basefilename)
{
	int i,k;
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);
	int SEED=tv.tv_usec;
	srand48(SEED);

	char * filename=new char[35]; 

	cerr << basefilename << endl ;

	sprintf(filename,"%s.cb",basefilename);
    ifstream cbfile (filename);		// Control, BP
	sprintf(filename,"%s.cd",basefilename);
    ifstream cdfile (filename);		// Control, DL
	sprintf(filename,"%s.sb",basefilename);
    ifstream sbfile (filename);		// Stress, BP
	sprintf(filename,"%s.sd",basefilename);
    ifstream sdfile (filename);		// Stress, DL

    if(!cbfile || !cdfile || !sbfile || !sdfile){
        cerr << "\nOne or More of the Stats files not found.\n";
        exit(1);
	}	

// ************* CONVENTION *************/
// xyzp; x \in {c,s} - control and stress
//       y \in {b,d} - B.P. or D.L.
// 		 z \in {a,b} - apical or basal
// 		 p \in {m,v} - mean or variance
// *************************************/

	double * cdam; double * cdav; double * cdbm; double * cdbv;	
	double * cbam; double * cbav; double * cbbm; double * cbbv;

	double * sdam; double * sdav; double * sdbm; double * sdbv;
	double * sbam; double * sbav; double * sbbm; double * sbbv;

	cdam=new double[antra]; cdav=new double[antra];
	cbam=new double[antra]; cbav=new double[antra];
	cdbm=new double[bntra]; cdbv=new double[bntra];
	cbbm=new double[bntra]; cbbv=new double[bntra];
	sdam=new double[antra]; sdav=new double[antra];
	sbam=new double[antra]; sbav=new double[antra];
	sdbm=new double[bntra]; sdbv=new double[bntra];
	sbbm=new double[bntra]; sbbv=new double[bntra];

	for (i=0; i<antra; i++){	// Apical means
		cdfile >> cdam[i]; cbfile >> cbam[i];
		sdfile >> sdam[i]; sbfile >> sbam[i];
	}	
	for (i=0; i<antra; i++){	// Apical variances
		cdfile >> cdav[i]; cbfile >> cbav[i];
		sdfile >> sdav[i]; sbfile >> sbav[i];
	}	
	for (i=0; i<bntra; i++){	// Basal means
		cdfile >> cdbm[i]; cbfile >> cbbm[i];
		sdfile >> sdbm[i]; sbfile >> sbbm[i];
	}	
	for (i=0; i<bntra; i++){	// Basal variances
		cdfile >> cdbv[i]; cbfile >> cbbv[i];
		sdfile >> sdbv[i]; sbfile >> sbbv[i];
	}	
	cdfile.close(); cbfile.close();
	sdfile.close(); sbfile.close();

	double a,b,max;

	abpzeta=new double[sdata.na];
	bbpzeta=new double[sdata.nb];
	adlzeta=new double[sdata.na];
	bdlzeta=new double[sdata.nb];

	double dlsum, bpsum;
	double dlcnt, bpcnt;
	dlsum=bpsum=0.0;
	dlcnt=bpcnt=0.0;
	
	double temp;

//  Generate Mean Zeta on Apical Side for DL and BP //

	for(i=0; i<sdata.na; i++){
		if(i>(antra-1)) k=antra-1;
		else k=i;
		temp=finratio(sdam[k],sdav[k],cdam[k],cdav[k]);
		if(temp==-1) {
			if(adlzeta[i-1])
				adlzeta[i]=adlzeta[i-1]/sdata.apicalcount[i-1];
			else adlzeta[i]=drand48();
		}	
		else
			adlzeta[i]=temp;
		
		cerr << adlzeta[i] << "\t\t\t" ;

		temp=finratio(sbam[k],sbav[k],cbam[k],cbav[k]);
		if(temp==-1) {
			if(abpzeta[i-1])
				abpzeta[i]=abpzeta[i-1]/BPStats.apicalcount[i-1];
			else abpzeta[i]=drand48();
		}	
		else
			abpzeta[i]=temp;
		
		cerr << abpzeta[i] << endl ;

		adlzeta[i] *= sdata.apicalcount[i];
		abpzeta[i] *= BPStats.apicalcount[i];
		
		dlsum += adlzeta[i] ;
		bpsum += abpzeta[i];
		dlcnt += sdata.apicalcount[i];
		bpcnt += BPStats.apicalcount[i];

		//cerr << adlzeta[i] << "  " << sdata.apicalcount[i] 
			 //<< "\t\t\t" << abpzeta[i] << "  " 
			 //<< BPStats.apicalcount[i] << endl ;
	}

	cerr << endl << endl ;

//  Generate Mean Zeta on Basal Side for DL and BP //

	for(i=0; i<sdata.nb; i++){	
		if(i>(bntra-1)) k=bntra-1;
		else k=i;
		temp=finratio(sdbm[k],sdbv[k],cdbm[k],cdbv[k]);
		if(temp==-1) {
			if(bdlzeta[i-1])
				bdlzeta[i]=bdlzeta[i-1]/sdata.basalcount[i-1];
			else bdlzeta[i]=drand48();	
		}	
		else
			bdlzeta[i]=temp;
		
		cerr << bdlzeta[i] << "\t\t\t" ;

		temp=finratio(sbbm[k],sbbv[k],cbbm[k],cbbv[k]);
		if(temp==-1) {
			if(bbpzeta[i-1])
				bbpzeta[i]=bbpzeta[i-1]/BPStats.basalcount[i-1];
			else bbpzeta[i]=drand48();	
		}	
		else
			bbpzeta[i]=temp;
		
		cerr << bbpzeta[i] << endl ;

		bdlzeta[i] *= sdata.basalcount[i];
		bbpzeta[i] *= BPStats.basalcount[i];

		dlsum += bdlzeta[i] ;
		bpsum += bbpzeta[i];
		dlcnt += sdata.basalcount[i];
		bpcnt += BPStats.basalcount[i];

		//cerr << bdlzeta[i] << "  " << sdata.basalcount[i] 
			 //<< "\t\t\t" << bbpzeta[i] << "  " 
			 //<< BPStats.basalcount[i] << endl ;
	}

	cerr << "\nReduction in DL: " << dlsum ;
	cerr << "\nTotal DL:        " << dlcnt ;
	cerr << "\nReduction in BP: " << bpsum ;
	cerr << "\nTotal BP:        " << bpcnt ;
	cerr << endl ;
	
}	
/*****************************************************************/


void Prune :: setApicalMaxZeta(char * basefilename)
{
	int i;
	setMaxZeta(basefilename);

//  Set Zeta on Basal Side for DL and BP to zero//

	for(i=0; i<sdata.nb; i++){	
		bdlzeta[i]=0;
		bbpzeta[i]=0;
	}

	cerr << "\nSetting basal specifications to zero ....";

}	
/*****************************************************************/

void Prune :: setBasalMaxZeta(char * basefilename)
{
	int i;
	setMaxZeta(basefilename);

//  Set Zeta on Apical Side for DL and BP to zero //

	for(i=0; i<sdata.na; i++){
			adlzeta[i]=0;
			abpzeta[i]=0;
	}
	cerr << "\nSetting apical specifications to zero ....";
}	

/*****************************************************************/

double Prune :: finratio (double sm, double sv, double cm, double cv)
{
    double max=maxxratio(sm,sv,cm,cv);
    if(!max) return -1;

    if (!cv && !sv) return 1-(sm/cm);
    else if(!sv) return 1-max/cv;
    else if(!cv) return 1-max*sv;
    else return 1-(max*sv)/cv;
}

/*****************************************************************/

double Prune :: maxxratio (double sm, double sv, double cm, double cv)
{
	double a,b;
	if (sv) a=sm/sv; 			// Avoid division by zero
	else if (!sm) return 0;		// Variance is zero and mean is zero
	else a=sm;					// Variance is zero and mean is nonzero

	if (cv) b=cm/cv;
	else if (!cm) return 0; 
	else b=cm;

	double max=0.0;
	double t;
	double maxt;
	for(t=-3; t<5; t+=0.001){
		if(max<ratio(a,b,t)) {
			max=ratio(a,b,t);
			maxt=t;
		} 
	}	
	return maxt;		
}	

/*****************************************************************/

