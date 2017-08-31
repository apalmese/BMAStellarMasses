import smass_only
import helperfunctions

basedir = "/Users/apalmese/work/scripts/BMAStellarMasses/"
libdir = basedir+"lib/simha_miles_test/"
infile = basedir+"in/redmapper_test.fits"
outfile= basedir+"out/redmapper_test.fits"
outfile_smass= basedir+"out/redmapper_test_smass.fits"
libtype = "test"

inputDataDict=helperfunctions.read_redmapper_smassonly(infile)
smass_only.calc(inputDataDict, outfile=outfile_smass, indir=libdir, lib=libtype)
helperfunctions.redmapper_output(infile,outfile_smass, outfile)