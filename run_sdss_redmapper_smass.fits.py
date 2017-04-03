import smass
indir = "/data/des60.b/data/palmese/lambda_star/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="/data/des60.b/data/palmese/lambda_star/sdss/sdss_redmapper_smass.fits_masses.out"
import helperfunctions
inputDataDict=helperfunctions.read_fits("/data/des60.b/data/palmese/lambda_star/sdss/sdss_redmapper_smass.fits")
smass.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
