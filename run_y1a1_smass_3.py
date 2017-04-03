import smass_only
indir = "/data/des60.b/data/palmese/lambda_star/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="/data/des60.b/data/palmese/stellarmasses_y1a1/splashback/parallel_stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_3.out"
import helperfunctions
inputDataDict=helperfunctions.read_fits_smassonly("/data/des60.b/data/palmese/stellarmasses_y1a1/splashback/stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_3.fits")
smass_only.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
