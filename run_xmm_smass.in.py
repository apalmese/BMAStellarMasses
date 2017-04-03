import smass
indir = "/data/des60.b/data/palmese/lambda_star/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="/data/des60.b/data/palmese/lambda_star/Xray/xmm_masses.out"
import helperfunctions
inputDataDict=helperfunctions.read_xmm("/data/des60.b/data/palmese/lambda_star/Xray/xmm_smass.in")
smass.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
