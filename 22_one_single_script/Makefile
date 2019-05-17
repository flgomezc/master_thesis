#wget https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/corner_data_cut_0.dat
SHELL := /bin/bash # Use bash syntax

N_OBS = 26


FILENAME = corner_data_cut_
URL_PATH = https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/

OC_PATH = observed_catalogs/
RC_PATH = random_catalogs/
FC_PATH = full_catalogs/
ML_PATH = masterlists/
VE_PATH = volume_and_excentricity/
FG_PATH = figures/
BS_PATH = xdl_beta_skeleton/


all : folders download_cuts random_cats
	@echo Boing


# Create working directories
folders :
	mkdir -p $(OC_PATH)
	mkdir -p $(RC_PATH)
	mkdir -p $(FC_PATH)
	mkdir -p $(ML_PATH)
	mkdir -p $(VE_PATH)
	mkdir -p $(FG_PATH)
	mkdir -p $(BS_PATH)

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
