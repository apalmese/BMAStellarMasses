# BMAStellarMasses

## Introduction
This code uses Bayesian Model Averaging (BMA) to compute stellar masses and other properties (age, metallicity, k-corrections, SFR...) for galaxies and galaxy clusters. The BMA is run over a set of templates, and here we provide some of those sets computed using Conroy and Gunn 2010 FSPS.
This code can be used to compute galaxy properties of single galaxies, or to produce a mass proxy \mu* for galaxy clusters. This mass proxy is defined and studied in Palmese et al. 2017. Refer to Palmese et al. 2017 for further details about the code. 

This code was developed by Jim Annis and adapted by Antonella Palmese.

## Usage 1: compute stellar mass for single galaxies
Prepare the universe input seds that are given by FSPS:
```
import addMinority
import addRest
addMinority.doall(indir)
addRest.doall(indir)
```
addMinority adds a lower metallicity population to a single metallicity population SSP. addRest adds rest frame colors.
dir is the directory that contains the universe.

Compute stellar masses for single galaxies:
```
import smass_only
indir = "/data/des60.b/data/palmese/lambda_star/fsps_v3.0_modified_Nov16/OUTPUTS/simha_miles_Nov2016/"
outfile="stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_1.out"
import helperfunctions
inputDataDict=helperfunctions.read_fits_smassonly("stellar_masses_mags_y1a1_mof_BPZ_splashback_cuts_1.fits")
smass_only.calc(inputDataDict, outfile=outfile, indir=indir, lib="miles")
```

## Usage 2: compute stellar mass for galaxy clusters and the cluster mass proxy \mu*
To compute clusters stellar masses and \mu*
```
import clusterSMass
clusterSMass.haloStellarMass(filename="/your/path/to/file/from/smass",outfile="/your/path/to/output")
```


## smass.py outputs
* Age [Gyr] is the mass weighted age, currently only available for Simha SFH
* sSFR [(M_{sun}/yr)/M_{sun}]

## clusterSMass.py outputs
