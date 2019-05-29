import argparse

description = "This script generate a Random Catalog, the Full Catalog and calls the Xiao-Dong Li's Fortran Beta-Skeleton Calculator"
epilog = "At the end, the script stores a BetaSkeleton (.bsk) file."

parser = argparse.ArgumentParser(description=description, epilog=epilog)

parser.add_argument('filein', type=str,
                    help='Name of the Observed Catalog file, must be stored in the folder OC_PATH="./observed_catalgos/"')

parser.add_argument('filenumber', type=int,
                     help='The consecutive number of the file, this is necessary to generate the Beta-Skeleton file.')

parser.add_argument('-b', '--beta', type=float,
                    default=1.0,
                    help='Beta Skeleton Value, a float value "b>=1". Default Value = 1.0')

parser.add_argument('-n', '--nrand', type=float,
                    default=1.0,
                    help='The ratio between Number of Random Points and Number of Observational Points (nrand= N_random/N_obs)')
parser.add_argument('-T', '--TEST', 
                    action='store_true',
                    default=False,
                    help='Tests filenames and folders generating empty files, does not runs the hard calculations.')

arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
OC_FILE_IN = arg.filein
FILENUM = arg.filenumber
TEST = arg.TEST

print(arg)

print("output with filename {}, filenumber {}, and optional input parameters beta {} and nrand {}. TEST={}".format(OC_FILE_IN, FILENUM,BETA, nrand, TEST))