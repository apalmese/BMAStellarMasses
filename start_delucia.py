import smass
indir = "/data/des60.b/data/palmese/VT_clusters/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="delucia_masses.txt"
import helperfunctions
inputDataDict=helperfunctions.read_delucia("/data/des60.b/data/palmese/VT_clusters/millennium/DeLucia_smass.in")
smass.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
