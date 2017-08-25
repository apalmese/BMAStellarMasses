import smass_only
import helperfunctions

basedir = "/Users/apalmese/work/scripts/BMAStellarMasses/"
libdir = basedir+"lib/simha_miles_test/"
inputdir = basedir+"in/"
infile = "test.fits"
outfile= basedir+"out/test"
libtype = "test"

inputDataDict=helperfunctions.read_fits_smassonly(inputdir+infile)
smass_only.calc(inputDataDict, outfile=outfile, indir=libdir, lib=libtype)
