
#include"Prune.h"

/*****************************************************************/

void Prune :: pruneTP(int cn, int sn)
{
	double tx1,ty1,tz1;
	double tx2,ty2,tz2;
	double ipdist;
	createLine(PruCom[PruCom[cn].ppi-1],PruCom[cn]);

	if(ls -> num==2){	//	Current pt is last pt in comptmt.
		tx1=PruCom[PruCom[cn].ppi-1].x;
		ty1=PruCom[PruCom[cn].ppi-1].y;
		tz1=PruCom[PruCom[cn].ppi-1].z;
		removeCompartment(cn,sn);
		cmr++;
   	}
	else{
		// Remove Last Pt.
		tx1=ls->x[ls->num-2];
		ty1=ls->y[ls->num-2];
		tz1=ls->z[ls->num-2];
		PruCom[cn].x = ls -> x[ls->num-2];
		PruCom[cn].y = ls -> y[ls->num-2];
		PruCom[cn].z = ls -> z[ls->num-2];
	}	
	tx2=ls -> x[ls->num-1];
	ty2=ls -> y[ls->num-1];
	tz2=ls -> z[ls->num-1];
	ipdist=sqrt((tx1-tx2)*(tx1-tx2)+ (ty1-ty2)*(ty1-ty2)+
				(tz1-tz2)*(tz1-tz2));

	if(PruCom[cn].tc==4)				// Apical
		adlprune[sn] += ipdist/RESLN;
	else if (PruCom[cn].tc==3)			// Basal
		bdlprune[sn] += ipdist/RESLN;
	dlr += ipdist/RESLN;	
}	

/*****************************************************************/
	// If the removed compartment's parent is BP, remove the point 
	// from the TP list, from the BP list (if it is not a 
	// trifurcation), and update the BP prune numbers and
	// probabilities.
/*****************************************************************/
	
void Prune :: rcParentBP(int cn, int sn)
{
	int i,j;

	if (PruCom[cn].tc==3)	// Basal
		bbpprune[sn]++;
	else if (PruCom[cn].tc==4)	// Apical
		abpprune[sn]++;

	// Remove pt from TP list

	int SET=0;

	for(i=0; i<TP.npts; i++)
		if(TP.Pts[i]==cn){
			SET=1;
			break;
		}
	
	if(!SET){
		cerr <<"\nrcParentBP :: Error, TP not found\n\n";
		exit(1);
	}

	for(j=i;j<TP.npts-1;j++)	
		TP.Pts[j]=TP.Pts[j+1];
	TP.npts--;

	// Remove parent from bp list if that point is a bifurcation

	if(Children[PruCom[cn].ppi-1].npts==2){
		for(i=0; i<BP.npts; i++)
			if(BP.Pts[i]==(PruCom[cn].ppi-1)) break;
		for(j=i;j<BP.npts-1;j++)	
			BP.Pts[j]=BP.Pts[j+1];
		BP.npts--;
	}
}

/*****************************************************************/
	// If removed compartment's parent is soma, remove the 
	// compartment from the TP list and from the Stems list
/*****************************************************************/

void Prune :: rcParentSoma(int cn)
{
	int i,j;

	// Remove pt from tp list

	for(i=0; i<TP.npts; i++)
		if(TP.Pts[i]==cn) break;
	for(j=i;j<TP.npts-1;j++)	
		TP.Pts[j]=TP.Pts[j+1];
	TP.npts--;

	// Remove pt from stems list

	for(i=0; i<Stems.npts; i++)
		if(Stems.Pts[i]==cn) break;
	for(j=i;j<Stems.npts-1;j++)	
		Stems.Pts[j]=Stems.Pts[j+1];
	Stems.npts--;
}

/*****************************************************************/
// The case where the parent is neither soma nor BP
// Implies parent is also dendrite. So, update the TP to the current
// compartment's parent.
/*****************************************************************/

void Prune :: rcParentDendrite(int cn)
{
	int i,j;

	for(i=0; i<TP.npts; i++)
		if(TP.Pts[i]==cn) break;

	TP.Pts[i]=PruCom[cn].ppi-1;	

// If parent (which is the current TP) has an undefined TC (-1),
// remove that from TP list.

	if(PruCom[TP.Pts[i]].tc!=3 && PruCom[TP.Pts[i]].tc!=4){	
		for(j=i;j<TP.npts-1;j++)	
			TP.Pts[j]=TP.Pts[j+1];
		TP.npts--;
		cerr << "\nWarning: Undefined compartment removed from TP list.\n\n";
	}	
}

/*****************************************************************/
	// Remove comptmt cn from PruCom.
	// Update the indices of the SWC file to accomodate the 
	// removal of the compartment number cn. Also updates the parents'
	// indices, the TP, BP, Soma list.
/*****************************************************************/

void Prune :: rmcmpt (int cn)
{
	int i;

	if((PruCom[cn].tc != 3) && (PruCom[cn].tc != 4)){
		cerr << "\nrmcmpt :: GRAVE ERROR\n";
		exit(1);
	}	


	for(i=cn; i<npc-1; i++){
		PruCom[i]=PruCom[i+1];
		PruCom[i].pi --;		// One comptmt is removed.
	}	
	npc--;

// Update Parent Indices to accomodate above update.

	for(i=0; i<npc; i++){
		if((PruCom[i].ppi-1)>cn)
			PruCom[i].ppi--;	// One Comptmt is removed.
	}	

// Update lists of TP, BP, Stems, Soma to accomodate above.

	for(i=0; i<TP.npts; i++)
		if(TP.Pts[i]>cn) TP.Pts[i]--;

	for(i=0; i<BP.npts; i++)
		if(BP.Pts[i]>cn) BP.Pts[i]--;

	for(i=0; i<Stems.npts; i++)
		if(Stems.Pts[i]>cn) Stems.Pts[i]--;

	for(i=0; i<Soma.npts; i++)
		if(Soma.Pts[i]>cn) Soma.Pts[i]--;

}


/*****************************************************************/
	// Remove the current compartment totally.
	// If the parent to this compartment were a BP, then, 
	// remove that BP from the list of BP's and update BP 
	// probability; sn stands for the sholl's segment number
/*****************************************************************/

void Prune :: removeCompartment(int cn, int sn)
{
	if(isParentBP(cn))			//	If parent is BP
		rcParentBP(cn,sn);
	else if(isParentSoma(cn))	//	If parent is Soma
		rcParentSoma(cn);
	else 						// If parent is dendrite
		rcParentDendrite(cn);

	rmcmpt(cn);
}


/*****************************************************************/
//****************  BPPrune starts here **************************/
/*****************************************************************/

void Prune :: findTPBP()
{
	int i;
	double dist;
	int sn;

	int trgt;
	for(i=0; i<TP.npts; i++){
		trgt=TP.Pts[i];
		bpPruneP[i]=1.0;
		while((isParentBP(trgt)==0) && (isParentSoma(trgt)==0)){
			dist  = 0.0;
			dist += (somax-PruCom[trgt].x)*(somax-PruCom[trgt].x);
			dist += (somay-PruCom[trgt].y)*(somay-PruCom[trgt].y);
			dist += (somaz-PruCom[trgt].z)*(somaz-PruCom[trgt].z);
			dist  = sqrt(dist);
			sn   = (int)(dist/(double)STLENGTH); 

			if(PruCom[trgt].tc==3)			// Basal
				bpPruneP[i]*=pow(bdlP[sn], disttoParent(trgt));
			else if(PruCom[trgt].tc==4)		// Apical
				bpPruneP[i]*=pow(adlP[sn], disttoParent(trgt));
			trgt=PruCom[trgt].ppi-1;
		}
		if(isParentBP(trgt))
			TPBP[i]=PruCom[trgt].ppi-1;
		else if (isParentSoma(trgt))
			TPBP[i]=-2;		
		else{
			cerr << "\nTPBP :: There exists some error...\n";
			cerr << PruCom[trgt] ;
			cerr << PruCom[PruCom[trgt].ppi-1] ;
			exit(1);
		}	
	}	
}

/*****************************************************************/
// Remove BP which has its TP as ntp'th TP. This is accomplished by
// continually pruning TP[ntp] until a BP is encountered.
/*****************************************************************/

void Prune :: pruneBP(int ntp)
{
	int trgt;

	// Actually remove the part connecting the TP to the BP by
	// recursively pruning TP[ntp] until a branching point is
	// encountered.

	int ret=0;
	int cmcnt=0;

	while(!ret){
		trgt=TP.Pts[ntp];
		ret=removeCmprtmnt(trgt);
		cmcnt++;
	}
}

/*****************************************************************/

int Prune :: removeCmprtmnt(int cn)
{
	double dist  = 0.0;
    dist += (somax-PruCom[cn].x)*(somax-PruCom[cn].x);
    dist += (somay-PruCom[cn].y)*(somay-PruCom[cn].y);
    dist += (somaz-PruCom[cn].z)*(somaz-PruCom[cn].z);
    dist = sqrt(dist);

	// Segment where current compartment belongs to

    int sn = (int)(dist/(double)STLENGTH); 

	if (PruCom[cn].tc==3) 
		bdlprune[sn]+=disttoParent(cn);
	else if (PruCom[cn].tc==4) 
		adlprune[sn]+=disttoParent(cn);
	else{
		cerr << "\nGravest error\n";
		exit(1);
	}	
	int ret=0;

	if(isParentBP(cn)){
		ret=1;
		rcParentBP(cn, sn);
	}	
	else if (isParentSoma(cn)){
		cerr << "\nrC :: Grave Error!\n";
	}	
	else
		rcParentDendrite(cn);

	rmcmpt(cn);  	// Removing comptmt cn from PruCom.

	return ret;
}

/*****************************************************************/

