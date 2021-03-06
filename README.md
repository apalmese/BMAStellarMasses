# BMAStellarMasses

## Introduction
This code uses Bayesian Model Averaging (BMA) to compute stellar masses and other properties (age, metallicity, k-corrections, SFR...) for galaxies and galaxy clusters. The BMA is run over a set of templates, and here we provide some of those sets computed using Conroy and Gunn 2010 FSPS.
This code can be used to compute galaxy properties of single galaxies, or to produce a mass proxy \mu* for galaxy clusters. This mass proxy is defined and studied in Palmese et al. 2017. Refer to Palmese et al. 2017 for further details about the code. 

This code was developed by Jim Annis and adapted by Antonella Palmese.

## Usage 1: compute stellar mass for single galaxies

### Step 1

Prepare the universe input seds that are given by FSPS:
```
import addMinority
import addRest
addMinority.doall(indir)
addRest.doall(indir)
```
addMinority adds a lower metallicity population to a single metallicity population SSP. addRest adds rest frame colors.
indir is the directory that contains the universe.

### Step 2

Compute stellar masses for single galaxies by running: ```python start_smass_only.py``` in src/. All that is needed is to modify the paths and filenames in that file:

```
import smass_only
import helperfunctions

basedir = "your-path-to-BMAStellarMasses/"
libdir = basedir+"lib/simha_miles_test/"
inputdir = basedir+"in/"
infile = "test.fits"
outfile= basedir+"out/test"
libtype = "test"

inputDataDict=helperfunctions.read_fits_smassonly(inputdir+infile)
smass_only.calc(inputDataDict, outfile=outfile, indir=libdir, lib=libtype)
```

with libtype="miles" and libdir = basedir+"lib/simha_miles_Nov2016/" for a non-test usage. In that case, also the lib/simha_miles_Nov2016.tar.gz file in the directory needs to be decompressed in order to make use of the full template set available. The test run only takes into account a few templates in order to make the computation quicker.

The galaxies given in input in "infile" expects by default a fits file that includes the following columns:
* 'ID': galaxy ID
* 'MAG_G', 'MAGERR_G', 'MAG_R', 'MAGERR_R, 'MAG_I', 'MAGERR_I', 'MAG_Z', 'MAGERR_Z': DES griz magnitudes and errors
* 'Z': galaxy/cluster redshift

The output is a fits file named as 'outfile' and contains the following columns:
* 'ID': galaxy ID as in input
* 'Z': galaxy/cluster redshift as in input
* 'gr_o', 'gi_o' rest frame g-r and g-i colors, with relative errors
* 'kri', 'kii' k-corrections
* 'distmod' distance modulus
* 'rabs', 'iabs' absolute magnitudes
* 'mass', 'mass_err' stellar mass in Log([M_Sun]) and error
* 'zmet' BMA of the metallicity in Z/Z_Sun
* 'best_zmet' best fitting model metallicity

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
