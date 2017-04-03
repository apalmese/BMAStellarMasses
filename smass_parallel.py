#Antonella Mar 2017
#Modified smass script that parallelizes the jobs


import tempfile
import shutil
import os
import numpy as np
import loadPopColors
import helperfunctions
from joblib import Parallel, delayed, load, dump
from scipy.stats import chi2
import CosmologicalDistance

def gal_smass(id,i,ierr,gr,ri,iz,gre,rie,ize,zed, splines, zmet,galaxy):

    print galaxy,"  at z = ",zed
    cd = CosmologicalDistance.CosmologicalDistance()
    ldistDict = dict()
    rest_gr,rest_gi,weight = [],[],[]
    masslight, sfrs,ages, zmets,kii, kri = [],[],[],[],[],[]
    minChiSq = 999; spIndex = -1
    for sp in range(0,len(splines)) :
        # for speed
        skey = str(sp) + "-" + str(zed)
        sgr = splines[sp][0](zed)
        sri = splines[sp][1](zed)
        siz = splines[sp][2](zed)
        sgrr = splines[sp][4](zed) ;# restframe g-r
        sgir = splines[sp][5](zed) ;# restframe g-i
        skii = splines[sp][6](zed) ;# kcorrection: i_o - i_obs
        skri = splines[sp][7](zed) ;# kcorrection: r_o - i_obs
        sml = splines[sp][8](zed) ;# log(mass/light)  (M_sun/L_sun)
        ssfr = splines[sp][9](zed)
        sage_cosmic = splines[sp][10](zed)
        sage = splines[sp][11](zed)
        szmet = zmet[sp]
                #To be changed if SFH changes

        gr_chisq = pow((gr - sgr)/gre,2)
        ri_chisq = pow((ri - sri)/rie,2)
        iz_chisq = pow((iz - siz)/ize,2)
        rest_gr.append(sgrr)
        rest_gi.append(sgir)
        kii.append(skii)
        kri.append(skri)
        masslight.append(sml)
        sfrs.append(ssfr)
        ages.append(sage)
        zmets.append(szmet)
        chisq = gr_chisq + ri_chisq + iz_chisq
        probability = 1-chi2.cdf(chisq, 3-1) ;# probability of chisq greater than this
        weight.append(probability)
    spIndex = np.argmax(weight)
    rest_gr = np.array(rest_gr)
    rest_gi = np.array(rest_gi)
    kii = np.array(kii)
    kri = np.array(kri)
    masslight = np.array(masslight)
    sfrs = np.array(sfrs)
    ages = np.array(ages)
    weight = np.array(weight)
    gr_weighted = rest_gr * weight
    gi_weighted = rest_gi * weight
    kii_weighted = kii * weight
    kri_weighted = kri * weight
    masslight_weighted = masslight * weight
    sfr_weighted = sfrs * weight
    age_weighted = ages * weight
    zmet_weighted = zmets * weight
    w1 = weight.sum()
    w2 = (weight**2).sum()
    if w1 == 0 : w1 = 1e-10
    if w2 == 0 : w2 = 1e-10
    mean_gr  = gr_weighted.sum()/w1
    mean_gi  = gi_weighted.sum()/w1
    mean_kii = kii_weighted.sum()/w1
    mean_kri = kri_weighted.sum()/w1
    mean_masslight = masslight_weighted.sum()/w1
    mean_sfr = sfr_weighted.sum()/w1
    mean_age = age_weighted.sum()/w1
    mean_zmet = zmet_weighted.sum()/w1
    # unbiased weighted estimator of the sample variance
    w3 = w1**2 - w2
    if w3 == 0 : w3 = 1e-10
    var_gr = ( w1/w3 ) * (weight*(rest_gr - mean_gr)**2).sum()
    var_gi = ( w1/w3 ) * (weight*(rest_gi - mean_gi)**2).sum()
    var_kii = ( w1/w3 ) * (weight*(kii - mean_kii)**2).sum()
    var_kri = ( w1/w3 ) * (weight*(kri - mean_kii)**2).sum()
    var_masslight = ( w1/w3 ) * (weight*(masslight - mean_masslight)**2).sum()
    var_sfr = ( w1/w3 ) * (weight*(sfrs - mean_sfr)**2).sum()
    var_age = ( w1/w3 ) * (weight*(ages - mean_age)**2).sum()
    var_zmet = ( w1/w3 ) * (weight*(zmets - mean_zmet)**2).sum()
    std_gr = var_gr**0.5
    std_gi = var_gi**0.5
    std_kii = var_kii**0.5
    std_kri = var_kri**0.5
    std_masslight = var_masslight**0.5
    std_sfr = var_sfr**0.5
    std_age = var_age**0.5
    std_zmet = var_zmet**0.5
    if std_gr > 99.99 : std_gr = 99.99
    if std_gi > 99.99 : std_gi = 99.99
    if std_kii > 99.99 : std_kii = 99.99
    if std_kri > 99.99 : std_kri = 99.99
    if std_sfr > 99.99 : std_sfr = 99.99
    if std_age > 99.99 : std_age = 99.99
    if std_masslight > 99.99 : std_masslight = 99.99
    if std_zmet > 99.99 : std_zmet = 99.99
    # Comment -distanceModulus out for fsps versions <2.5, as their mags don't include distance modulus
    if zed in ldistDict :
        lumdist = ldistDict[zed]
    else :
        lumdist = cd.luminosity_distance(zed) ;# in Mpc
        ldistDict[zed] = lumdist
    distanceModulus = 5*np.log10(lumdist*1e6/10.)
    iabs = i + mean_kii - distanceModulus
    rabs = i + mean_kri - distanceModulus
    taMass = taylorMass(mean_gi, iabs)
    mcMass = mcintoshMass(mean_gr, rabs)
    fsMass = fspsMass( mean_masslight, iabs )
    # JTA: to make purely distance modulus
    #iabs = i[galaxy] - distanceModulus 
    #fsMass = gstarMass( iabs )

    # saving for output
    # perhaps: out_id[galaxy] = id[galaxy]

    return [id, zed, mean_gr , std_gr, mean_gi, std_gi, mean_kii, std_kii, mean_kri, std_kri , i, distanceModulus , iabs,\
    rabs, mcMass, taMass, fsMass, std_masslight,spIndex,zmets[spIndex],mean_sfr,std_sfr, mean_age, std_age, mean_zmet, std_zmet ]


def calc (inputDataDict, outfile, indir="simha/", lib="miles") :
    from scipy.stats import chi2
    import CosmologicalDistance
    import os

    n_outputs = 26
    Njobs = 3

    cd = CosmologicalDistance.CosmologicalDistance()
    ldistDict = dict()
    splineDict = dict()
    splines, zmet = loadPopColors.doall(indir, lib=lib)
    id  = inputDataDict["id"]
    i  = inputDataDict["i"]
    ierr  = inputDataDict["ierr"]
    gr  = inputDataDict["gr"]
    ri  = inputDataDict["ri"]
    iz  = inputDataDict["iz"]
    grerr  = inputDataDict["grerr"]
    rierr  = inputDataDict["rierr"]
    izerr  = inputDataDict["izerr"]
    zed = inputDataDict["zed"]
    size=i.size
    #size=3

    # protect against too small of errors => values = 0
    ix = np.nonzero(grerr < 0.02)
    grerr[ix] = 0.02
    ix = np.nonzero(rierr < 0.02)
    rierr[ix] = 0.02
    ix = np.nonzero(izerr < 0.02)
    izerr[ix] = 0.02

    size = 10

    #folder = tempfile.mkdtemp()
    #outputs_name = os.path.join(folder, 'out')

    # Pre-allocate a writeable shared memory map as a container for the
    # results of the parallel computation
    #outputs = np.memmap(outputs_name,dtype='float32',shape=((size,26)), mode='w+')

    # put all inputs in one array, it's horrible this way...

    outputs = Parallel(n_jobs=Njobs)(delayed(gal_smass)(id[galaxy],i[galaxy],ierr[galaxy],gr[galaxy],ri[galaxy],iz[galaxy],grerr[galaxy],rierr[galaxy],izerr[galaxy],zed[galaxy], splines,zmet, galaxy) for galaxy in range(0,size))


    #print outputs
    #print outputs.shape
    header = "# id, photoz, gr_o,  std,  gi_o, std, kri, std, kii, std, i,  distmod, "
    header = header + "r_abs, i_abs, mcMass, taMass, mass, std, zed, sfr, sfrstd, age, agestd, bestsp,best_zmet,mean_zmet \n"
    fd = open(outfile,"w")
    fd.write(header)
    fd.close()
    np.savetxt(outfile+".dat", outputs,fmt="%d,%6.3f,%6.4f,%6.3f,%6.4f,%6.3f,%6.4f,%6.3f,%6.4f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%d,%6.3f,%6.4f,%6.3f,%6.3f,%6.3f,%6.3f, %6.3f")
    os.system("cat {} >> {}; rm {}".format(outfile+".dat", outfile, outfile+".dat"))



def taylorMass (gi, iabs) :
    # equation 8, taylor et al 2011 MNRAS, V418, Issue 3, pp. 1587-1620
    mass = -0.68 + 0.70*gi - 0.4*(iabs - 4.58)
    # assumes h=0.7
    out_stdgi = np.array(out_stdgi)
    out_kii = np.array(out_kii)
    out_stdkii = np.array(out_stdkii)
    out_kri = np.array(out_kri)
    out_stdkri  = np.array(out_stdkri )
    out_iobs = np.array(out_iobs)
    out_distmod = np.array(out_distmod)
    out_iabs = np.array(out_iabs)
    out_rabs = np.array(out_rabs)
    out_mass_gr = np.array(out_mass_gr)
    out_mass_gi = np.array(out_mass_gi)
    out_mass = np.array(out_mass)
    out_stdmass = np.array(out_stdmass)
    out_sfr = np.array(out_sfr)
    out_stdsfr = np.array(out_stdsfr)
    out_age = np.array(out_age)
    out_stdage = np.array(out_stdage)
    out_bestsp = np.array(out_bestsp)
    data = np.array([out_id, out_gr, out_stdgr, out_gi, out_stdgi, \
        out_kri, out_stdkri, out_kii, out_stdkii, out_iobs, out_distmod, \
        out_rabs, out_iabs, out_mass_gr, out_mass_gi, out_mass, out_stdmass, \
        out_sfr, out_stdsfr, out_age, out_stdage, out_bestsp])
    header = "# id, gr_o,  std,  gi_o, std, kri, std, kii, std, i,  distmod, "
    header = header + "r_abs, i_abs, mcMass, taMass, mass, std, sfr, sfrstd, age, agestd, bestsp\n"

    data = np.array([out_id, out_gr, out_stdgr, out_gi, out_stdgi, \
        out_kri, out_stdkri, out_kii, out_stdkii, out_iobs, out_distmod, \
        out_rabs, out_iabs, out_mass_gr, out_mass_gi, out_mass, out_stdmass, out_bestsp])
    header = "# id, gr_o,  std,  gi_o, std, kri, std, kii, std, i,  distmod, "
    header = header + "r_abs, i_abs, mcMass, taMass, mass, std, bestsp\n"
    fd = open(outfile,"w")
    fd.write(header)
    fd.close()
    np.savetxt(outfile+".dat", data.T, "%d,%6.3f,%6.4f,%6.3f,%6.4f,%6.3f,%6.4f,%6.3f,%6.4f,\
         %6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.4f,%6.3f,%6.3f,%6.3f,%6.3f,%d")
    os.system("cat {} >> {}; rm {}".format(outfile+".dat", outfile, outfile+".dat"))


def mcintoshMass (gr, rabs, h=0.7) :
    # equation 1, McIntosh et al 2014, arXiv:1308.0054v2
    mass = -0.406 + 1.097*gr - 0.4*(rabs - 5*np.log10(h) - 4.64)
    # log(mass/h^2) 
    return mass

def taylorMass (gi, iabs) :
    # equation 8, taylor et al 2011 MNRAS, V418, Issue 3, pp. 1587-1620
    mass = -0.68 + 0.70*gi - 0.4*(iabs - 4.58)
    # assumes h=0.7
    return mass

def fspsMass( masstolight , iabs) :
    # following the above, which assumes h=0.7
    mass =  masstolight - 0.4*(iabs - 4.58)
    return mass
# how far off is just distance modulus?
def gstarMass(iabs) :
    mass =  - 0.4*(iabs - 4.58)
    return mass

