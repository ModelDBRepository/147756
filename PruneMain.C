
#include"Prune.h"

main(int, char ** argv)
{
	Prune prune;
	prune.LoadParams(argv[1]);
	prune.PruneTree();
	//prune.PrintFinal();
	return 0;
}
