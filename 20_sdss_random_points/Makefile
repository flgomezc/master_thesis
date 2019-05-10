

data:
	mkdir -p observed_catalogs
	#mkdir -p random_catalogs
	mkdir -p full_catalogs
	#wget https://github.com/forero/abacus/raw/master/data/SDSS/SDSS_data_random_Planck15.txt
	#mv SDSS_data_random_Planck15.txt random_catalogs
	wget https://github.com/forero/abacus/raw/master/data/SDSS/SDSS_data_Planck15.txt
	mv SDSS_data_Planck15.txt observed_catalogs

clean:
	rm *~
	rm -r full_catalogs
	rm -r masterlists
	rm -r volume_and_excentricity
	rm -r xdl_beta_skeleton
