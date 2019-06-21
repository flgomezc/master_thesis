#wget https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/corner_data_cut_0.dat
SHELL := /bin/bash # Use bash syntax

N_OBS = 1
FILENUM = 0
FILEPREFIX = corner_data_cut_
FILENAME = $(FILEPREFIX)$(FILENUM).dat

URL_PATH = https://github.com/forero/abacus/raw/master/data/AbacusCosmos_720box_planck_00_0_rockstar_halos/z0.1/cuts/

OC_PATH = ./observed_catalogs/
RC_PATH = ./random_catalogs/
FC_PATH = full_catalogs/
ML_PATH = masterlists/
VE_PATH = volume_and_excentricity/
FG_PATH = figures/
BS_PATH = xdl_beta_skeleton/
AN_PATH = analysis/

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
	mkdir -p $(AN_PATH)

# Download "SDSS PLANK15 cuts"
download_cuts :
	@number=0 ; while [[ $$number -le $(N_OBS) ]] ; do \
		echo $(URL_PATH)$(FILEPREFIX)$$number.dat ; \
		wget $(URL_PATH)$(FILEPREFIX)$$number.dat ; \
		mv $(FILENAME)$$number.dat $(OBS_PATH) ; \
		(( number = number + 1 )) ; \
	done

# Generate random catalogs
main_squence :
	python calculate_beta_skeleton.py $(FILENAME) $(FILENUM)
	python void_identifier-recursion.py $(FILENAME) $(FILENUM)
	python ellipsoid_analysis.py $(FILENAME) $(FILENUM)
	python morphologycal_analysis.py $(FILENAME) $(FILENUM)


clean :
	rm -f *~
	rm -f *.LOG
	rm -rf $(OC_PATH)
	rm -rf $(RC_PATH)
	rm -rf $(FC_PATH)
	rm -rf $(ML_PATH)
	rm -rf $(VE_PATH)
	rm -rf $(FG_PATH)
	rm -rf $(BS_PATH)
	rm -rf $(AN_PATH)
