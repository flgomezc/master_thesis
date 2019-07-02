# Where is NGL
EXEC_PATH=/home/filipo/master_thesis/ngl-beta/build/binsrc/

# Number of Particles
NPARTICLES=5000
# Void Type: SV MSV IV MIV
# (Spherical Void, Multiple Spherical Voids, Irregular Void, Multiple I. V.)
VOID_TYPE=MIV
# Skeleton parameters
METHOD=BSkeleton
BETA=1.0
DIMENSION=3

FILE_IN := $(shell mktemp  XXXXXXXXXX )

FILENAME=FC_N$(NPARTICLES)_$(VOID_TYPE).cat
FILE_SKEL=FC_N$(NPARTICLES)_$(VOID_TYPE)_Beta$(BETA).bsk

$(FILE_SKEL): $(FILE_IN)
	@echo Input File $(FILENAME)
	cp $(FILENAME) $(FILE_IN)
	@echo 
	@echo Calling NGL
	$(EXEC_PATH)./getNeighborGraph -i $(FILE_IN) -d $(DIMENSION) -m $(METHOD) -b $(BETA)c > $(FILE_SKEL)
	@echo Output file: $(FILE_SKEL)
	@echo Deleting auxiliar file
	rm $(FILE_IN)


clean:
	rm $(FILE_SKEL)
