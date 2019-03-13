EXEC_PATH=/home/filipo/master_thesis/ngl-beta/build/binsrc/

FILENAME=4SV_B60.cat

# Skeleton parameters
METHOD=BSkeleton
BETA=1.3

FILE_IN=$(FILENAME)
FILE_SKEL=SV_R20_BS60_B_$(BETA).bsk


$(FILE_SKEL): $(FILE_IN)
	echo Input File to NGL: $(FILE_IN)
	$(EXEC_PATH)./getNeighborGraph -i $(FILE_IN) -d 3 -b $(BETA) -m $(METHOD) > $(FILE_SKEL)


clean:
	rm $(FILE_SKEL)
#rm -f *~
