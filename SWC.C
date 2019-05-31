#include"SWC.h"


/*****************************************************************/
// Overloading SWCData for I/0
/*****************************************************************/

ostream & operator << (ostream & tmp, SWCData swcd)
{
	tmp << swcd.pi << " " ;
	tmp << swcd.tc << " " ;
	tmp << swcd.x << " " ;
	tmp << swcd.y << " " ;
	tmp << swcd.z << " " ;
	tmp << swcd.r << " " ;
	tmp << swcd.ppi << endl ;
	return tmp;
}	

/*****************************************************************/

istream & operator >> (istream & tmp, SWCData & swcd)
{
	tmp >> swcd.pi ;
	tmp >> swcd.tc ;
	tmp >> swcd.x  ;
	tmp >> swcd.y  ;
	tmp >> swcd.z  ;
	tmp >> swcd.r  ;
	tmp >> swcd.ppi;
	return tmp;
}	


/*****************************************************************/
// SWC Class Files 
/*****************************************************************/

SWC :: SWC ()
{
	ncomp=0;
}

/*****************************************************************/

void SWC :: ReadSWC(char * filename)
{
	char * flname=new char [35];
	if(!filename){
		cout << "\nGive input SWC filename : ";
		cin >> flname ;
	}
	else strcpy(flname,filename);

	ifstream swcfile (flname);
	if(!swcfile){
		cerr << "\nFile " << flname << " does not exist.\n";
		exit(1);
	}	

 	char c;
    char * buffer=new char[80];
	streampos crps;
	int i;


	swcfile.get(c);
	if (c=='#'){

		// Skip through header...

		while (c == '#'){
			swcfile.getline(buffer,80);
			swcfile.get(c);
		}	

		crps=swcfile.tellg();
	}
	else crps=0;

	// Browse through the end of file to get the number of points....

	ncomp=0;
	while(!swcfile.eof()){
		swcfile.getline(buffer,80);
		ncomp ++;
	}

	ncomp --;
	Comptmt=new SWCData[ncomp];
	
	// Read all points...

	swcfile.close();

	ifstream swc1file;
	swc1file.open(flname);
	if(!swc1file){
		cerr << "\nFile " << flname << " does not exist.\n";
		exit(1);
	}	
	swc1file.seekg(crps);
	
	for(i=0; i<ncomp; i++){
		swc1file >> Comptmt[i];
		if(Comptmt[i].tc==-1)
			cerr << "\nWarning: Undefined compartment type at " 
				 << Comptmt[i] ;
	}		

	delete(buffer);
	delete(flname);
}
	
/*****************************************************************/

void SWC :: isRead()
{
	if(!ncomp){
		cerr << "\nSWC file has been read.\n";
		exit(1);
	}	
}

/*****************************************************************/

void SWC :: WriteSWC(char * filename)
{
	isRead ();

	char * flname=new char [35];
	if(!filename){
		cout << "\nGive output SWC filename : ";
		cin >> flname ;
	}
	else strcpy(flname,filename);     

	ofstream outfile (flname);
	int i;
	
	for(i=0; i<ncomp; i++)
		outfile << Comptmt[i] ;

	delete(flname);
}	
	
/*****************************************************************/

void SWC :: SetComptmt(SWCData * Cpts, int npts)
{
	int i;
	ncomp=npts;
	Comptmt=new SWCData[ncomp];
	for(i=0; i<ncomp; i++)
		Comptmt[i]=Cpts[i];
}
/*****************************************************************/
