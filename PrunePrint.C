
#include"Prune.h"

/*****************************************************************/

void Prune :: generateSpecs()
{
	ofstream adlfile("CA3b4p.adl");
	ofstream bdlfile("CA3b4p.bdl");
	ofstream abpfile("CA3b4p.abp");
	ofstream bbpfile("CA3b4p.bbp");
	ofstream tdlfile("CA3b4p.tdl");
	ofstream tbpfile("CA3b4p.tbp");
	int i,j;

	double adlsum=0.0;
	double bdlsum=0.0;
	double abpsum=0.0;
	double bbpsum=0.0;

	double adlorg=0.0;
	double bdlorg=0.0;
	double abporg=0.0;
	double bbporg=0.0;

	for(i=0; i<sdata.na; i++){
		adlorg+=sdata.apicalcount[i];
		abporg+=BPStats.apicalcount[i];
	}	

	for(i=0; i<sdata.nb; i++){
		bdlorg+=sdata.basalcount[i];
		bbporg+=BPStats.basalcount[i];
	}	

	for(j=0; j<10000; j++){
		setZeta(basefilename);
		adlsum=abpsum=bdlsum=bbpsum=0.0;

		for(i=0; i<sdata.na; i++){
			//adlfile << sdata.apicalcount[i]-adlzeta[i] << " ";
			adlfile << adlzeta[i]/sdata.apicalcount[i] << " ";
			adlsum += sdata.apicalcount[i]-adlzeta[i];
		}	
		adlfile << endl ;
	
		for(i=0; i<sdata.nb; i++){
			//bdlfile << sdata.basalcount[i]-bdlzeta[i] << " ";
			bdlfile << bdlzeta[i]/sdata.basalcount[i] << " ";
			bdlsum += sdata.basalcount[i]-bdlzeta[i];
		}	
		bdlfile << endl ;
			
		for(i=0; i<sdata.na; i++){	
			//abpfile << BPStats.apicalcount[i]-abpzeta[i] << " " ;
			abpfile << abpzeta[i]/BPStats.apicalcount[i] << " " ;
			abpsum += BPStats.apicalcount[i]-abpzeta[i];
		}	
		abpfile << endl ;
			
		for(i=0; i<sdata.nb; i++){
			//bbpfile << BPStats.basalcount[i]-bbpzeta[i] << " ";
			bbpfile << bbpzeta[i]/BPStats.basalcount[i] << " ";
			bbpsum += BPStats.basalcount[i]-bbpzeta[i];
		}
		bbpfile << endl ;
		
		tdlfile << adlorg << " " << adlsum << " " 
				<< ((adlorg-adlsum)*100)/adlorg  << " "
				<< bdlorg << " " << bdlsum << " "
				<< ((bdlorg-bdlsum)*100)/bdlorg  << " "
				<< ((adlorg+bdlorg-adlsum-bdlsum)*100)/(adlorg+bdlorg) 
				<< endl ;
		tbpfile	<< abporg << " " << abpsum << " "
				<< ((abporg-abpsum)*100)/abporg  << " "
				<< bbporg << " " << bbpsum << " "
				<< ((bbporg-bbpsum)*100)/bbporg  << " " 
				<< ((abporg+bbporg-abpsum-bbpsum)*100)/(abporg+bbporg) 
				<< endl ;

		cerr << j << " .....\n";
	}
}

/*****************************************************************/

void Prune:: printSubsets()
{
	int i;
	
	cerr << "\nThe set of TP's are:  ";
	for(i=0; i<TP.npts; i++)
		cerr << TP.Pts[i] << "   ";

	cerr << "\nThe set of BP's are:  ";
	for(i=0; i<BP.npts; i++)
		cerr << BP.Pts[i] << "   ";


	cerr << "\nThe set of Stems are:  ";
	for(i=0; i<Stems.npts; i++)
		cerr << Stems.Pts[i] << "   ";

	cerr << "\nThe set of Soma are:  ";
	for(i=0; i<Soma.npts; i++)
		cerr << Soma.Pts[i] << "   ";

	cerr << endl ;

}
/*****************************************************************/

void Prune :: printStatus()
{
	int i;
	Sholl test;
	test.SetComptmt(PruCom,npc);
	ShollData stest=test.ShollAnalysis();
	cerr << "\nADLPrune:";
	for (i=0; i<sdata.na; i++)
		cerr << (int)adlprune[i] << " ";
	cerr << "\nADLZeta: ";
	for (i=0; i<sdata.na; i++)
		cerr << (int)adlzeta[i] << " ";
	cerr << "\nADLOrig:";
	for (i=0; i<sdata.na; i++)
		cerr << (int)sdata.apicalcount[i] << " ";
	cerr << "\nADLDiff: ";
	for (i=0; i<sdata.na; i++)
		cerr << (int)(sdata.apicalcount[i] - adlzeta[i]) << " ";
	cerr << "\nADLFinal:";
	for(i=0; i<stest.na; i++)
		cerr << (int) stest.apicalcount[i] << " " ;
	cerr << endl ;	


	cerr << endl << endl << "ABPPrune: ";
	for (i=0; i<sdata.na; i++)
		cerr << (int)abpprune[i] << " ";
	cerr << endl << "ABPZeta: ";
	for (i=0; i<sdata.na; i++)
		cerr << (int)abpzeta[i] << " ";
	cerr << endl << "ABPOrig: ";
	for (i=0; i<sdata.na; i++)
		cerr << (int)BPStats.apicalcount[i] << " ";

	cerr << endl << endl << "BDLPrune: ";
	for (i=0; i<sdata.nb; i++)
		cerr << (int)bdlprune[i] << " ";
	cerr << endl << "BDLZeta: ";
	for (i=0; i<sdata.nb; i++)
		cerr << (int)bdlzeta[i] << " ";
	cerr << endl << "BDLOrig: ";
	for (i=0; i<sdata.nb; i++)
		cerr << (int)sdata.basalcount[i] << " ";

	cerr << endl << endl << "BBPPrune: " ;
	for (i=0; i<sdata.nb; i++)
		cerr << (int)bbpprune[i] << " ";
	cerr << endl << "BBPZeta: " ;
	for (i=0; i<sdata.nb; i++)
		cerr << (int)bbpzeta[i] << " ";
	cerr << endl << "BBPOrig: " ;
	for (i=0; i<sdata.nb; i++)
		cerr << (int)BPStats.basalcount[i] << " ";
	cerr << endl ;
}	

/*****************************************************************/

void  Prune :: PrintFinal()
{
	int i;

	ofstream totreddl;

	ofstream specabp;
	ofstream specbbp;
	ofstream specadl;
	ofstream specbdl;

	ofstream shollabp;
	ofstream shollbbp;
	ofstream sholladl;
	ofstream shollbdl;

	totreddl.open("4_totred.dl",ios::app);

	totreddl << dlsum << endl ;


	specabp.open("4_spe.abp",ios::app);
	for (i=0; i<sdata.na; i++)
		specabp << BPStats.apicalcount[i]-abpzeta[i] << " ";
	specabp << endl ;	

	specbbp.open("4_spe.bbp",ios::app);
	for (i=0; i<sdata.nb; i++)
		specbbp << BPStats.basalcount[i]-bbpzeta[i] << " ";
	specbbp << endl ;	


	specadl.open("4_spe.adl",ios::app);
	for (i=0; i<sdata.na; i++)
		specadl << sdata.apicalcount[i]-adlzeta[i] << " ";
	specadl << endl ;	

	specbdl.open("4_spe.bdl",ios::app);
	for (i=0; i<sdata.nb; i++)
		specbdl << sdata.basalcount[i]-bdlzeta[i] << " ";
	specbdl << endl ;	

	Sholl test;
	test.SetComptmt(PruCom,npc);
	ShollData stest=test.ShollAnalysis();
	cerr << stest.na << "\t" << stest.nb << "\n";
	getBPHist();
	

	shollabp.open("4_sho.abp",ios::app);
	for (i=0; i<stest.na; i++)
		shollabp << (int)BPStats.apicalcount[i] << " ";
	shollabp << endl ;	

	shollbbp.open("4_sho.bbp",ios::app);
	for (i=0; i<stest.nb; i++)
		shollbbp << (int)BPStats.basalcount[i] << " ";
	shollbbp << endl ;	

	sholladl.open("4_sho.adl",ios::app);
	for(i=0; i<stest.na; i++)
		sholladl << (int) stest.apicalcount[i] << " " ;
	sholladl << endl ;	

	shollbdl.open("4_sho.bdl",ios::app);
	for(i=0; i<stest.nb; i++)
		shollbdl << (int) stest.basalcount[i] << " " ;
	shollbdl << endl ;	
}
/*****************************************************************/


