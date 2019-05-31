
#include "Prune.h"
#include <iomanip.h>

main(int argc, char **argv)
{
	if(argc !=3){
		cerr << "\nProgram to perform analysis on a SWC file\n\n";
		cerr << "Usage: " << argv[0] 
			 << " <Input SWC File> <Output Filename> \n\n";
		exit(1);
	}	
	Prune sholl;
	sholl.ReadSWC(argv[1]);
	sholl.ShollAnalysis();
	sholl.WriteSholl(argv[2]);
	sholl.WriteBPStats(argv[2]);
	return 0;
}
