# BMAStellarMasses

Prepare the universe input seds that are given by FSPS:

import addMinority
import addRest
addMinority.doall(indir)
addRest.doall(indir)

addMinority adds a lower metallicity population to a single metallicity population SSP. addRest adds rest frame colors.
dir is the directory that contains the universe.

Compute stellar masses for single galaxies:

import smass_only
indir = "/data/des60.b/data/palmese/lambda_star/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_1.out"
import helperfunctions
inputDataDict=helperfunctions.read_fits_smassonly("stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_1.fits")
smass_only.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")

