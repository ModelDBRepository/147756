#include"Prune.h"

/*****************************************************************/
// Polar Method (Box, Muller and Marsaglia) for two independent
// normal distributions
/*****************************************************************/

twodouble Prune :: normal (double mu1, double var1, double mu2,
										double var2)
{
	double v1, v2, s;
	do {
		v1 = 2 * drand48() - 1;
		v2 = 2 * drand48() - 1;
		s = v1*v1 + v2*v2;
	}
	while (s >= 1 || s < 1e-30);
	s = sqrt((-2*log(s))/s);
	v1 *= s;
	v2 *= s;
	twodouble x;
	x.x1=v1*sqrt(var1)+mu1;
	x.x2=v2*sqrt(var2)+mu2;
	return (x);
}


/*****************************************************************/
// Taylor Series expansion to find area under Normal distribution
// P(z\le x)=0.5 + {{1}\over{\sqrt{2\pi}}} \sum_{k=0}^\infty 
// {{(-1)^k x^{2*k+1}} \over {(2k+1) 2^k k!}}
// z is standard normal.
/*****************************************************************/

double Prune :: pnorm(double x) // Gaussian CDF
{
	if(x<-6.5) return 0;
	if(x>6.5) return 1;
	double factK=1;
	double sum=0.0;
	double term=1.0;
	int k=0;
	while(fabs(term)>exp(-23)){
		term=(1.0/sqrt(2.0*M_PI))*pow(-1,k)*pow(x,2*k+1)/
						((2*k+1)*pow(2,k)*factK);
		sum+=term;
		k++;
		factK*=k;
	}
	sum+=0.5;
	if(sum<1e-10) sum=0;
	return sum;
}

/*****************************************************************/
// Returns pdf of (a+x)/(b+y), where x and y are independent 
// standard normal distributed.
/*****************************************************************/

double Prune :: ratio(double a, double b, double t)
{
	double q,tf;
	q=(b+a*t)/sqrt(1+t*t);
	tf=(exp(-0.5*(a*a+b*b))/(M_PI*(1+t*t)))*(1+(q/dnorm(q))*
					(pnorm(q)-pnorm(0)));
	return tf;				
}

/*****************************************************************/

double Prune :: dnorm(double x) 	// Gaussian PDF
{
	return ((1.0/sqrt(2.0*M_PI))*exp(-x*x/2));
}

/*****************************************************************/
// The actual values of pruning are generated here
/*****************************************************************/


void Prune :: setZeta(char * basefilename)
{
	int i,k;
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);
	int SEED=tv.tv_usec;
	srand48(SEED);
	
	char * filename=new char[35]; 

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

	double bpsum;
	double dlcnt, bpcnt;
	dlsum=bpsum=0.0;
	dlcnt=bpcnt=0.0;

	abpzeta=new double[sdata.na];
	bbpzeta=new double[sdata.nb];
	adlzeta=new double[sdata.na];
	bdlzeta=new double[sdata.nb];

//  Find maximum ratio on Apical side for D.L. and B.P. //
//  Generate Zeta on Apical Side for DL and BP //

	for(i=0; i<sdata.na; i++){
		if(i>(antra-1)) k=antra-1;
		else k=i;
		max=maxratio(sdam[k], sdav[k], cdam[k], cdav[k]);
		adlzeta[i]=generateZeta(sdam[k], sdav[k], cdam[k], cdav[k], max);
		adlzeta[i] *= sdata.apicalcount[i];

		max=maxratio(sbam[k], sbav[k], cbam[k], cbav[k]);
		abpzeta[i]=generateZeta(sbam[k], sbav[k], cbam[k], cbav[k], max);
		abpzeta[i] *= BPStats.apicalcount[i];

		dlsum += adlzeta[i] ;
		bpsum += abpzeta[i];
		dlcnt += sdata.apicalcount[i];
		bpcnt += BPStats.apicalcount[i];

	}	
	

//  Find maximum ratio on Basal side for DL and BP  //
//  Generate Zeta on Basal Side for DL and BP //

	for(i=0; i<sdata.nb; i++){
		if(i>(bntra-1)) k=bntra-1;
		else k=i;
		max=maxratio(sdbm[k], sdbv[k], cdbm[k], cdbv[k]);
		bdlzeta[i]=generateZeta(sdbm[k], sdbv[k], cdbm[k], cdbv[k], max);
		bdlzeta[i] *= sdata.basalcount[i];

		max=maxratio(sbbm[k], sbbv[k], cbbm[k], cbbv[k]);
		bbpzeta[i]=generateZeta(sbbm[k], sbbv[k], cbbm[k], cbbv[k], max);
		bbpzeta[i] *= BPStats.basalcount[i];

		dlsum += bdlzeta[i] ;
		bpsum += bbpzeta[i];
		dlcnt += sdata.basalcount[i];
		bpcnt += BPStats.basalcount[i];

	}	

/*
	cerr << "\nReduction in DL: " << dlsum ;
	cerr << "\nTotal DL:        " << dlcnt ;
	cerr << "\nReduction in BP: " << bpsum ;
	cerr << "\nTotal BP:        " << bpcnt ;
	cerr << endl ;
*/

	delete(cdam); delete(cdav); 
	delete(cbam); delete(cbav); 
	delete(cdbm); delete(cdbv); 
	delete(cbbm); delete(cbbv); 
	delete(sdam); delete(sdav); 
	delete(sbam); delete(sbav); 
	delete(sdbm); delete(sdbv); 
	delete(sbbm); delete(sbbv); 
}

/*****************************************************************/

double Prune :: generateZeta(double sm, double sv, double cm, 
							double cv, double maxratio)
{							 
	if(!maxratio) return drand48();

	double vara=1,varb=1;
	double a,b;
	if (sv) a=sm/sv; 			// Avoid division by zero
	else if (!sm) return 0;		// Variance is zero and mean is zero
	else {a=sm; vara=0;}		// Variance is zero and mean is nonzero
		
	if (cv) b=cm/cv;
	else if (!cm) return 0; 
	else {b=cm; varb=0;}
		
	twodouble zeta;
	double t;
	double dzeta;
	double min;
	if(sm>cm) min=1.0;
	else min=sm/cm;
	do{
		if(!vara && !varb) {
			return 1-sm/cm;
		}
		else if(!vara){
			zeta=normal(sm,0,b,1);
			t=zeta.x1/zeta.x2;
			dzeta=t/cv;  
		}	
		else if(!varb){
			zeta=normal(a,1,cm,0);
			t=zeta.x1/zeta.x2;
			dzeta=t*sv;
		}	
		else {
			zeta=normal(a,1,b,1);	// This is ORIG.
			//zeta=normal(sm,sv,cm,cv);			// TEST
			t=zeta.x1/zeta.x2;
			dzeta=(t*sv)/(cv);
		}	
	}
	while(ratio(a,b,t)<(0.1*maxratio) || dzeta > 1.0 || dzeta < 0.0);
	//while(dzeta > 1.0 || dzeta < (min*0.88));
	//while(ratio(a,b,t)<(0.65*maxratio) || dzeta > 1.0 || dzeta < (min*0.85));
	return (1-dzeta);
}	

/*****************************************************************/

double Prune :: maxratio (double sm, double sv, double cm, double cv)
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
	for(t=-3; t<5; t+=0.01)
		if(max<ratio(a,b,t)) 
			max=ratio(a,b,t);
	return max;		
}	

/*****************************************************************/

void Prune :: setProbabilities()
{
	setDLProbs();
	setBPProbs();
}

/*****************************************************************/

void Prune :: setBPProbs()
{
	int i;
	double P;

	// Apical

	for(i=0; i<sdata.na; i++){
		if(abpzeta[i]==0 || adlzeta[i]==0)
			abpP[i]=0; 
		else{
			if(abpprune[i]==0) abpprune[i]=1;
			if(adlprune[i]==0) adlprune[i]=RESLN;
			P=(adlzeta[i]*abpprune[i])/(abpzeta[i]*adlprune[i]);
			abpP[i]=pow(adlP[i],P);
		}
	}	

	// Basal

	for(i=0; i<sdata.nb; i++){
		if(bbpzeta[i]==0 || bdlzeta[i]==0)
			bbpP[i]=0; 
		else{
			if(bbpprune[i]==0) bbpprune[i]=1;
			if(bdlprune[i]==0) bdlprune[i]=RESLN;
			P=(bdlzeta[i]*bbpprune[i])/(bbpzeta[i]*bdlprune[i]);
			bbpP[i]=pow(bdlP[i],P);
		}
	}	
}

/*****************************************************************/

void Prune :: setDLProbs()
{
/* Set Probabilities for DL Pruning by just setting the sum of  *
 * probabilities to 1 by normalizing (zeta-prune) values.		*/

	double dlSum=0.0;

	maxdlP=0;

	int i;
	for(i=0; i<sdata.na; i++){
		adlP[i]=adlzeta[i]-adlprune[i];
		dlSum+=adlP[i];
	}	
	
	for(i=0; i<sdata.nb; i++){
		bdlP[i]=bdlzeta[i]-bdlprune[i];
		dlSum+=bdlP[i];
	}	

	for(i=0; i<sdata.na; i++){
		adlP[i]/=dlSum;
		if(maxdlP<adlP[i]) maxdlP=adlP[i];
		//cerr << "ADLP: " << i << "   "  << adlP[i] << endl;
	}	

	for(i=0; i<sdata.nb; i++){
		bdlP[i]/=dlSum;
		if(maxdlP<bdlP[i]) maxdlP=bdlP[i];
		//cerr << "BDLP: " << i << "   "  << bdlP[i] << endl;
	}	

	maxdlP += 0.1*maxdlP;
}
/*****************************************************************/
