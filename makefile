OBJS=SWC.o Sholl.o PruneBase.o PruneZeta.o PruneTPBP.o Prune.o\
	PruneMaxZeta.o PruneMain.o PrunePrint.o
SBJS=SWC.o Sholl.o PruneBase.o ShollMain.o

CC=g++ -g -Wno-deprecated

Prune: $(OBJS)
	$(CC) -o Prune $(OBJS) -lm

Sholl: $(SBJS)
	$(CC) -o Sholl $(SBJS) -lm

$(OBJS): SWC.h Sholl.h Prune.h
$(SBJS): SWC.h Sholl.h 

.C.o:
	$(CC) -c $*.C

clean:
	rm *.o
