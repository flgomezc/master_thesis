#wget https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/corner_data_cut_0.dat
SHELL := /bin/bash # Use bash syntax

N_OBS = 26


FILENAME = corner_data_cut_
URL_PATH = https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/

OBS_PATH = observed_catalogs/
RND_PATH = random_catalogs/
FLL_PATH = full_catalogs/
MLS_PATH = masterlists/
VEX_PATH = volume_and_excentricity/
FIG_PATH = figures/

all : folders download_cuts random_cats
	@echo Boing


# Create working directories
folders :
	mkdir -p $(OBS_PATH)
	mkdir -p $(RND_PATH)
	mkdir -p $(FLL_PATH)
	mkdir -p $(MLS_PATH)
	mkdir -p $(VEX_PATH)
	mkdir -p $(FIG_PATH)

# Download "SDSS PLANK15 cuts"
download_cuts :
	@number=0 ; while [[ $$number -le $(N_OBS) ]] ; do \
		echo $(URL_PATH)$(FILENAME)$$number.dat ; \
		wget $(URL_PATH)$(FILENAME)$$number.dat ; \
		mv $(FILENAME)$$number.dat $(OBS_PATH) ; \
		(( number = number + 1 )) ; \
	done

# Generate random catalogs
random_cats :
	python 01_calculate_beta_skeleton.py



clean :
	rm -f *~
	rm -f *.LOG
	rm -rf $(OBS_PATH)
	rm -rf $(RND_PATH)
	rm -rf $(FLL_PATH)
	rm -rf $(MLS_PATH)
	rm -rf $(VEX_PATH)
	rm -rf $(FIG_PATH)
