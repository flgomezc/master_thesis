


EXEC_PATH=/home/filipo/master_thesis/ngl-beta/build/binsrc/

# Skeleton parameters
METHOD=BSkeleton
BETA=1.7

FILE_NGL=test_in
FILE_SKEL=test_out


$(FILE_SKEL): $(FILE_NGL)
	$(EXEC_PATH)./getNeighborGraph -i $(FILE_NGL) -d 3 -b $(BETA) -m $(METHOD) > $(FILE_SKEL)


clean:
	rm $(FILE_SKEL)
	rm -f *~
