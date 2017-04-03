import smass
indir = "/data/des60.b/data/palmese/VT_clusters/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="partial_redmapper_sva1_00_masses.txt"
import helperfunctions
inputDataDict=helperfunctions.read("/home/s1/palmese/scripts/stellarmass/redmapper_test.in")
smass.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
