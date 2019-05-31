
#include"Sholl.h"


/*****************************************************************/
// ShollAnalyze Class Files 
/*****************************************************************/

Sholl :: Sholl ():SWC()
{
	ls=new LineStruct;
	ls -> x=new double[MAXPTS];
	ls -> y=new double[MAXPTS];
	ls -> z=new double[MAXPTS];

	sdata.na=0;
	sdata.nb=0;
}

/*****************************************************************/

void Sholl :: createLine (SWCData s1, SWCData s2)
{			
	double x,y,z;
	double x1, y1, z1, x2, y2, z2;
	x1=s1.x; y1=s1.y; z1=s1.z;
	x2=s2.x; y2=s2.y; z2=s2.z;

	double dist=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));

	int number=(int)(ceil(dist/RESLN))+1;


	if(number > MAXPTS){
		cerr << "ERROR: No. of pts in line greater than MAX and is: " 
			 << number << endl;
		exit(1);
	}	


	double b1=(x2-x1)/dist;		// Make it unit length vector
	double b2=(y2-y1)/dist;
	double b3=(z2-z1)/dist;

	x=x1;y=y1;z=z1;

	int i=0;

	for(double t=RESLN; t<=dist; t+=RESLN){
		ls -> x[i]=x1+t*b1;
		ls -> y[i]=y1+t*b2;
		ls -> z[i]=z1+t*b3;
		i++;
	}	

/* We have to add the last point to the Line Structure. However, if the
 * last pt generated is actually the the actual last point, we would not
 * want to do that.
 */

	if (i<=0) {
		ls -> x[0]=x1;	
		ls -> y[0]=y1;
		ls -> z[0]=z1;
	
		ls -> x[1]=x2;
		ls -> y[1]=y2;
		ls -> z[1]=z2;

		ls -> num=2;
		return;
	}	
	double tdist;
	tdist=fabs(x2-ls->x[i-1])+fabs(y2-ls->y[i-1])+fabs(z2-ls->z[i-1]);

	if(tdist!=0){ 	// If the last pt generated and the actual last point differ
		ls -> x[i]=x2;	// The last point is added; but not the first point.
		ls -> y[i]=y2;
		ls -> z[i]=z2;
		ls -> num=i+1;
	}
	else{ 
		ls -> num=i;
	}	
}	


/*****************************************************************/

ShollData & Sholl :: ShollAnalysis ()
{
	somax=0; somay=0; somaz=0;
	int j,i=0;
	int count=1;		// Initial count is 1 for reasons given below.

	isRead();

	for(i=1;i<ncomp;i++){	//0th point does not have valid parent
		if(Comptmt[i].tc==1){	// If the current point is part of soma
			createLine(Comptmt[Comptmt[i].ppi-1],Comptmt[i]);
			for(j=0; j<ls -> num; j++){
				somax+=ls -> x[j];
				somay+=ls -> y[j];
				somaz+=ls -> z[j];
				count ++;
			}
		}	
	}		

	somax += Comptmt[0].x;	// The first point is added;
	somay += Comptmt[0].y;	// Hence initial count above is 1.
	somaz += Comptmt[0].z;

	// Compute CoG of the Soma.

	somax /= count ;
	somay /= count ;
	somaz /= count ;

	cerr << "\nCentroid of Soma is: " << somax << "  " 
		 << somay << "  " << somaz << endl ;

	double dist;
	double tx1, tx2;
	double ty1, ty2;
	double tz1, tz2;
	double bmax=0.0;
	double amax=0.0;

	sdata.basalcount=new double[MAXTRACK];
	sdata.apicalcount=new double[MAXTRACK];
	double ipdist;

	for(i=0; i<MAXTRACK; i++){
		sdata.basalcount[i] =0;
		sdata.apicalcount[i] =0;
	}

	for(i=1;i<ncomp;i++){
		if(Comptmt[i].tc==3  || Comptmt[i].tc==4){ 	//  Sholl's for dendrites.
			createLine(Comptmt[Comptmt[i].ppi-1],Comptmt[i]);
			for(j=0; j<ls -> num; j++){
				dist=0.0;
				dist += (somax-ls->x[j])*(somax-ls->x[j]);
				dist += (somay-ls->y[j])*(somay-ls->y[j]);
				dist += (somaz-ls->z[j])*(somaz-ls->z[j]);
				dist=sqrt(dist);

				if(j == 0){	
					tx1=Comptmt[Comptmt[i].ppi-1].x;
					ty1=Comptmt[Comptmt[i].ppi-1].y;
					tz1=Comptmt[Comptmt[i].ppi-1].z;
				}	
				else{
					tx1=ls->x[j-1];
					ty1=ls->y[j-1];
					tz1=ls->z[j-1];
				}

				tx2=ls -> x[j];
				ty2=ls -> y[j];
				tz2=ls -> z[j];
						
				ipdist=sqrt((tx1-tx2)*(tx1-tx2)+ (ty1-ty2)*(ty1-ty2)+
							(tz1-tz2)*(tz1-tz2));

				if(Comptmt[i].tc==3){
					sdata.basalcount[(int)(dist/(double)STLENGTH)]+=
										ipdist/RESLN;
					if (dist>bmax) bmax=dist;
				}	
				if(Comptmt[i].tc==4){
					sdata.apicalcount[(int)(dist/(double)STLENGTH)]+=
										ipdist/RESLN;
					if (dist>amax) amax=dist;
				}	
			}	
		}	
	}	

	sdata.nb=1+(int)(bmax/(double)STLENGTH);
	sdata.na=1+(int)(amax/(double)STLENGTH);

	ComputeTotalLength();

	return sdata;
}

/*****************************************************************/

void Sholl :: WriteSholl(char * bfilename)
{
	isDone ();
	char * flname=new char [35];
	if(!bfilename){
		cout << "\nGive output SHL basefilename : ";
		cin >> flname ;
	}
	else strcpy(flname,bfilename);     
	
	char * filename=new char [35];
	sprintf(filename,"%s.asl",flname);
	ofstream aoutfile (filename);
	int i;
	
	for(i=0; i<sdata.na; i++){
		aoutfile << sdata.apicalcount[i] << endl ;
	}	

	sprintf(filename,"%s.bsl",flname);
	ofstream boutfile (filename);

	for(i=0; i<sdata.nb; i++){
		boutfile << sdata.basalcount[i] << endl ;
	}	

	//delete(flname);
	//delete(filename);
}	
	
/*****************************************************************/

void Sholl :: ComputeTotalLength()
{
	isDone ();
	totlen=0.0;

	int i;
	
	for(i=0; i<sdata.na; i++){
		totlen += sdata.apicalcount[i];
	}	


	for(i=0; i<sdata.nb; i++){
		totlen += sdata.basalcount[i];
	}	

	cerr << "Total dendritic length is: " << totlen << endl ;

}	
	
/*****************************************************************/

void Sholl :: isDone()
{
	if(!sdata.na){
		cerr << "\nSholl's analysis is not yet done.\n";
		exit(1);
	}	
}

/*****************************************************************/

