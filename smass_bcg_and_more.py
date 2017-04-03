import numpy as np
import loadPopColors

#
# read in the targets
#
def readVT (filename="../vtColors_annis.csv") :
    print "\t reading {}".format(filename)
    id,objID,ra,dec,modelMag_i,fracDev_i,devMag_i,expMag_i,\
        extinction_i,petromag_r,model_ug,model_gr,model_ri,\
        model_iz,modelMagErr_i,model_ug_err,model_gr_err,\
        model_ri_err,model_iz_err = \
        np.genfromtxt(filename,unpack=True,delimiter=",",skiprows=1)

    cmodel_i = -2.5*np.log10(fracDev_i *10**(-0.4*devMag_i) + (1-fracDev_i)*10**(-0.4*expMag_i))
    cmodel_i = cmodel_i - extinction_i

    return id,cmodel_i, model_ug, model_gr, model_ri, model_iz, \
        model_ug_err, model_gr_err, model_ri_err, model_iz_err


def readVTmembers( filename="../vtmembers.csv") :
    print "\t reading {}".format(filename)
    ra,dec,z,ze,id,central,hostid = np.genfromtxt(filename,unpack=True,delimiter=",",skiprows=1)
    return ra,dec,z,ze,id,central,hostid

#
# prep the population files
#
def prenike ( dir="simha/") :
    import addMinority
    import addRest
    addMinority.doall(dir)
    addRest.doall(dir)

def nike (output, indir="simha/", vtcFile = "../vtColors_annis.csv") :
    from scipy.stats import chi2
    import CosmologicalDistance
    cd = CosmologicalDistance.CosmologicalDistance()
    ldistDict = dict()
    splineDict = dict()
    splines = loadPopColors.doall(indir)
    memberData = readVTmembers()
    id,i,ug,gr,ri,iz,ugerr,grerr,rierr,izerr = readVT(vtcFile)

    # protect against too small of errors => values = 0
    ix = np.nonzero(ugerr < 0.02)
    ugerr[ix] = 0.02
    ix = np.nonzero(grerr < 0.02)
    grerr[ix] = 0.02
    ix = np.nonzero(rierr < 0.02)
    rierr[ix] = 0.02
    ix = np.nonzero(izerr < 0.02)
    izerr[ix] = 0.02

    gal_id = memberData[4]
    gal_z = memberData[2]

    # prepping for output
    gal_ra = memberData[0]
    gal_dec = memberData[1]
    gal_ze = memberData[3]
    gal_central = memberData[5]
    gal_hostid = memberData[6]

    # prepping for output
    out_ra, out_dec, out_z, out_ze = [],[],[],[]
    out_id, out_central, out_hostid = [],[],[]
    out_gr, out_stdgr, out_gi, out_stdgi = [],[],[],[]
    out_kii, out_stdkii, out_kri, out_stdkri = [],[],[],[]
    out_iobs, out_distmod, out_iabs, out_rabs, out_mass_gr, out_mass_gi= [],[],[],[],[],[]
    out_mass, out_stdmass = [],[]
    size = id.size
    #size = 10
    for galaxy in range(0,size) :
        ix = np.nonzero( gal_id == id[galaxy])
        zed = gal_z[ix][0]
        if zed == 0.0 : continue
        print galaxy,"  of  ",size

        rest_gr,rest_gi,weight = [],[],[]
        masslight, kii, kri = [],[],[]
        for sp in range(0,len(splines)) :
            # for speed
            skey = str(sp) + "-" + str(zed)
            if skey in splineDict :
                sug,sgr,sri,siz,sgrr,sgir,skii,skri,sml = splineDict[skey]
            else :
                sug = splines[sp][0](zed)
                sgr = splines[sp][1](zed)
                sri = splines[sp][2](zed)
                siz = splines[sp][3](zed)
                sgrr = splines[sp][4](zed) ;# restframe g-r
                sgir = splines[sp][5](zed) ;# restframe g-i
                skii = splines[sp][6](zed) ;# kcorrection: i_o - i_obs
                skri = splines[sp][7](zed) ;# kcorrection: r_o - i_obs
                sml = splines[sp][8](zed) ;# log(mass/light)  (M_sun/L_sun)
                splineDict[skey] = sug,sgr,sri,siz,sgrr,sgir,skii,skri,sml 

            uge = ugerr[galaxy]
            gre = grerr[galaxy]
            rie = rierr[galaxy]
            ize = izerr[galaxy]
            ug_chisq = pow((ug[galaxy] - sug)/uge,2)
            gr_chisq = pow((gr[galaxy] - sgr)/gre,2)
            ri_chisq = pow((ri[galaxy] - sri)/rie,2)
            iz_chisq = pow((iz[galaxy] - siz)/ize,2)
            rest_gr.append(sgrr)
            rest_gi.append(sgir)
            kii.append(skii)
            kri.append(skri)
            masslight.append(sml)
            chisq = ug_chisq + gr_chisq + ri_chisq + iz_chisq
            probability = 1-chi2.cdf(chisq, 3) ;# probability of chisq greater than this
            weight.append(probability)
        rest_gr = np.array(rest_gr)
        rest_gi = np.array(rest_gi)
        kii = np.array(kii)
        kri = np.array(kri)
        masslight = np.array(masslight)
        weight = np.array(weight)
        gr_weighted = rest_gr * weight
        gi_weighted = rest_gi * weight
        kii_weighted = kii * weight
        kri_weighted = kri * weight
        masslight_weighted = masslight * weight
        # weighted sum
        w1 = weight.sum()
        w2 = (weight**2).sum()
        if w1 == 0 : w1 = 1e-10
        if w2 == 0 : w2 = 1e-10
        mean_gr  = gr_weighted.sum()/w1
        mean_gi  = gi_weighted.sum()/w1
        mean_kii = kii_weighted.sum()/w1
        mean_kri = kri_weighted.sum()/w1
        mean_masslight = masslight_weighted.sum()/w1
        # unbiased weighted estimator of the sample variance
        w3 = w1**2 - w2
        if w3 == 0 : w3 = 1e-10
        var_gr = ( w1/w3 ) * (weight*(rest_gr - mean_gr)**2).sum()
        var_gi = ( w1/w3 ) * (weight*(rest_gi - mean_gi)**2).sum()
        var_kii = ( w1/w3 ) * (weight*(kii - mean_kii)**2).sum()
        var_kri = ( w1/w3 ) * (weight*(kri - mean_kii)**2).sum()
        var_masslight = ( w1/w3 ) * (weight*(masslight - mean_masslight)**2).sum()
        std_gr = var_gr**0.5
        std_gi = var_gi**0.5
        std_kii = var_kii**0.5
        std_kri = var_kri**0.5
        std_masslight = var_masslight**0.5
        if std_gr > 99.99 : std_gr = 99.99
        if std_gi > 99.99 : std_gi = 99.99
        if std_kii > 99.99 : std_kii = 99.99
        if std_kri > 99.99 : std_kri = 99.99
        if std_masslight > 99.99 : std_masslight = 99.99

        if zed in ldistDict :
            lumdist = ldistDict[zed]
        else :
            lumdist = cd.luminosity_distance(zed) ;# in Mpc
            ldistDict[zed] = lumdist
        distanceModulus = 5*np.log10(lumdist*1e6/10.)
        iabs = i[galaxy] - distanceModulus + mean_kii 
        rabs = i[galaxy] - distanceModulus + mean_kri 
        taMass = taylorMass(mean_gi, iabs) 
        mcMass = mcintoshMass(mean_gr, rabs) 
        fsMass = fspsMass( mean_masslight, iabs )

        # saving for output
        out_ra.append( gal_ra[ix][0] )
        out_dec.append( gal_dec[ix][0] )
        out_z.append( gal_z[ix][0] )
        out_ze.append( gal_ze[ix][0] )
        out_id.append( gal_id[ix][0] )
        out_central.append( gal_central[ix][0] )
        out_hostid.append( gal_hostid[ix][0] )
        out_gr.append( mean_gr )
        out_stdgr.append(std_gr )
        out_gi.append( mean_gi )
        out_stdgi.append( std_gi)
        out_kii.append( mean_kii )
        out_stdkii.append( std_kii )
        out_kri.append( mean_kri )
        out_stdkri .append( std_kri )
        out_iobs.append( i[galaxy] )
        out_distmod.append( distanceModulus )
        out_iabs.append( iabs )
        out_rabs.append( rabs )
        out_mass_gr.append( mcMass )
        out_mass_gi.append( taMass )
        out_mass.append( fsMass )
        out_stdmass.append( std_masslight )

    out_ra = np.array(out_ra)
    out_dec = np.array(out_dec)
    out_z = np.array(out_z)
    out_ze = np.array(out_ze)
    out_id = np.array(out_id).astype(int)
    out_central = np.array(out_central).astype(int)
    out_hostid = np.array(out_hostid).astype(int)
    out_gr = np.array(out_gr)
    out_stdgr = np.array(out_stdgr)
    out_gi = np.array(out_gi)
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
    data = np.array([out_ra, out_dec, out_z, out_ze, out_id, \
        out_central, out_hostid, out_gr, out_stdgr, out_gi, out_stdgi, \
        out_kri, out_stdkri, out_kii, out_stdkii, out_iobs, out_distmod, \
        out_rabs, out_iabs, out_mass_gr, out_mass_gi, out_mass, out_stdmass])
    np.savetxt(output, data.T, "%10.6f %10.5f %4.2f %5.3f %5d %d %10d \
        %6.3f %6.4f %6.3f %6.4f   %6.3f %6.4f %6.3f %6.4f \
         %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.4f")
    fd = open(output + ".hdr", "w")
    fd.write("# ra  dec  z zerr id  central hostid                    g-r_o std g-i_o ")
    fd.write("std   kri std kii std                  i  distmod r_abs i_abs mcMass taMass mass std\n")
    fd.close()


def taylorMass (gi, iabs) :
    # equation 8, taylor et al 2011 MNRAS, V418, Issue 3, pp. 1587-1620
    mass = -0.68 + 0.70*gi - 0.4*(iabs - 4.58)
    # assumes h=0.7
    return mass

def mcintoshMass (gr, rabs, h=0.7) :
    # equation 1, McIntosh et al 2014, arXiv:1308.0054v2
    mass = -0.406 + 1.097*gr - 0.4*(rabs - 5*np.log10(h) - 4.64)
    # log(mass/h^2) 
    return mass

def fspsMass( masstolight , iabs) :
    # following the above, which assumes h=0.7
    mass =  masstolight - 0.4*(iabs - 4.58)
    return mass


def readrmBCG (filename="../rmBCGColors_annis.csv") :
    print "\t reading {}".format(filename)
    id,ra,dec,z,cmodel_i,\
        petromag_r,model_ug,model_gr,model_ri,\
        model_iz,model_ug_err,model_gr_err,\
        model_ri_err,model_iz_err = \
        np.genfromtxt(filename,unpack=True,delimiter=",",skiprows=1)

    return id,ra,dec,z,cmodel_i, model_ug, model_gr, model_ri, model_iz, \
        model_ug_err, model_gr_err, model_ri_err, model_iz_err
def readLrg (filename="../lrgSpecColor_annis.csv") :
    print "\t reading {}".format(filename)
    sid,z,id,objID,ra,dec,cmodel_i,\
        petromag_r,model_ug,model_gr,model_ri,\
        model_iz,model_ug_err,model_gr_err,\
        model_ri_err,model_iz_err = \
        np.genfromtxt(filename,unpack=True,delimiter=",",skiprows=1)

    return id,ra,dec,z,cmodel_i, model_ug, model_gr, model_ri, model_iz, \
        model_ug_err, model_gr_err, model_ri_err, model_iz_err
def readMax (filename="../maxColors.csv") :
    print "\t reading {}".format(filename)
    id,objID,ra,dec,z,objra,objdec,extinction_i,cmodel_i,\
         petromag_r,model_ug,model_gr,model_ri,\
        model_iz,modelMagErr_i,model_ug_err,model_gr_err,\
        model_ri_err,model_iz_err = \
        np.genfromtxt(filename,unpack=True,delimiter=",",skiprows=1)

    return id,ra,dec,z,cmodel_i, model_ug, model_gr, model_ri, model_iz, \
        model_ug_err, model_gr_err, model_ri_err, model_iz_err
#
# a version for maxBCG
#   learning how to generalize the code
#
def nike2 (output, cat="max", indir="simha/") :
    from scipy.stats import chi2
    import CosmologicalDistance
    cd = CosmologicalDistance.CosmologicalDistance()
    ldistDict = dict()
    splineDict = dict()
    splines = loadPopColors.doall(indir)
    if cat == "max" :
        catFile = "../maxColors.csv"
        id,ra,dec,z,i,ug,gr,ri,iz,ugerr,grerr,rierr,izerr = readMax(catFile)
    elif cat == "lrg" :
        catFile = "../lrgSpecColor_annis.csv"
        id,ra,dec,z,i,ug,gr,ri,iz,ugerr,grerr,rierr,izerr = readLrg(catFile)
    elif cat == "rmBCG" :
        catFile = "../rmBCGColors_annis.csv"
        id,ra,dec,z,i,ug,gr,ri,iz,ugerr,grerr,rierr,izerr = readrmBCG(catFile)

    # protect against too small of errors => values = 0
    # 0.03 for maxBcg ?
    minErr = 0.02
    ix = np.nonzero(ugerr < minErr)
    ugerr[ix] = minErr
    ix = np.nonzero(grerr < minErr)
    grerr[ix] = minErr
    ix = np.nonzero(rierr < minErr)
    rierr[ix] = minErr
    ix = np.nonzero(izerr < minErr)
    izerr[ix] = minErr

    # prepping for output
    out_ra, out_dec, out_z = [],[],[]
    out_id = []
    out_gr, out_stdgr, out_gi, out_stdgi = [],[],[],[]
    out_kii, out_stdkii, out_kri, out_stdkri = [],[],[],[]
    out_iobs, out_distmod, out_iabs, out_rabs, out_mass_gr, out_mass_gi= [],[],[],[],[],[]
    out_mass, out_stdmass = [],[]
    size = id.size
    #size = 10
    for galaxy in range(0,size) :
        zed = z[galaxy]
        print galaxy,"  of  ",size

        rest_gr,rest_gi,weight = [],[],[]
        masslight, kii, kri = [],[],[]
        for sp in range(0,len(splines)) :
            # for speed
            skey = str(sp) + "-" + str(zed)
            if skey in splineDict :
                sug,sgr,sri,siz,sgrr,sgir,skii,skri,sml = splineDict[skey]
            else :
                sug = splines[sp][0](zed)
                sgr = splines[sp][1](zed)
                sri = splines[sp][2](zed)
                siz = splines[sp][3](zed)
                sgrr = splines[sp][4](zed) ;# restframe g-r
                sgir = splines[sp][5](zed) ;# restframe g-i
                skii = splines[sp][6](zed) ;# kcorrection: i_o - i_obs
                skri = splines[sp][7](zed) ;# kcorrection: r_o - i_obs
                sml = splines[sp][8](zed) ;# log(mass/light)  (M_sun/L_sun)
                splineDict[skey] = sug,sgr,sri,siz,sgrr,sgir,skii,skri,sml 

            uge = ugerr[galaxy]
            gre = grerr[galaxy]
            rie = rierr[galaxy]
            ize = izerr[galaxy]
            ug_chisq = pow((ug[galaxy] - sug)/uge,2)
            gr_chisq = pow((gr[galaxy] - sgr)/gre,2)
            ri_chisq = pow((ri[galaxy] - sri)/rie,2)
            iz_chisq = pow((iz[galaxy] - siz)/ize,2)
            rest_gr.append(sgrr)
            rest_gi.append(sgir)
            kii.append(skii)
            kri.append(skri)
            masslight.append(sml)
            chisq = ug_chisq + gr_chisq + ri_chisq + iz_chisq
            probability = 1-chi2.cdf(chisq, 3) ;# probability of chisq greater than this
            weight.append(probability)
        rest_gr = np.array(rest_gr)
        rest_gi = np.array(rest_gi)
        kii = np.array(kii)
        kri = np.array(kri)
        masslight = np.array(masslight)
        weight = np.array(weight)
        gr_weighted = rest_gr * weight
        gi_weighted = rest_gi * weight
        kii_weighted = kii * weight
        kri_weighted = kri * weight
        masslight_weighted = masslight * weight
        # weighted sum
        w1 = weight.sum()
        w2 = (weight**2).sum()
        if w1 == 0 : w1 = 1e-10
        if w2 == 0 : w2 = 1e-10
        mean_gr  = gr_weighted.sum()/w1
        mean_gi  = gi_weighted.sum()/w1
        mean_kii = kii_weighted.sum()/w1
        mean_kri = kri_weighted.sum()/w1
        mean_masslight = masslight_weighted.sum()/w1
        # unbiased weighted estimator of the sample variance
        w3 = w1**2 - w2
        if w3 == 0 : w3 = 1e-10
        var_gr = ( w1/w3 ) * (weight*(rest_gr - mean_gr)**2).sum()
        var_gi = ( w1/w3 ) * (weight*(rest_gi - mean_gi)**2).sum()
        var_kii = ( w1/w3 ) * (weight*(kii - mean_kii)**2).sum()
        var_kri = ( w1/w3 ) * (weight*(kri - mean_kii)**2).sum()
        var_masslight = ( w1/w3 ) * (weight*(masslight - mean_masslight)**2).sum()
        std_gr = var_gr**0.5
        std_gi = var_gi**0.5
        std_kii = var_kii**0.5
        std_kri = var_kri**0.5
        std_masslight = var_masslight**0.5
        if std_gr > 99.99 : std_gr = 99.99
        if std_gi > 99.99 : std_gi = 99.99
        if std_kii > 99.99 : std_kii = 99.99
        if std_kri > 99.99 : std_kri = 99.99
        if std_masslight > 99.99 : std_masslight = 99.99

        if zed in ldistDict :
            lumdist = ldistDict[zed]
        else :
            lumdist = cd.luminosity_distance(zed) ;# in Mpc
            ldistDict[zed] = lumdist
        distanceModulus = 5*np.log10(lumdist*1e6/10.)
        iabs = i[galaxy] - distanceModulus + mean_kii 
        rabs = i[galaxy] - distanceModulus + mean_kri 
        taMass = taylorMass(mean_gi, iabs) 
        mcMass = mcintoshMass(mean_gr, rabs) 
        fsMass = fspsMass( mean_masslight, iabs )

        # saving for output
        out_ra.append( ra[galaxy])
        out_dec.append( dec[galaxy])
        out_z.append( z[galaxy])
        out_id.append( id[galaxy])
        out_gr.append( mean_gr )
        out_stdgr.append(std_gr )
        out_gi.append( mean_gi )
        out_stdgi.append( std_gi)
        out_kii.append( mean_kii )
        out_stdkii.append( std_kii )
        out_kri.append( mean_kri )
        out_stdkri .append( std_kri )
        out_iobs.append( i[galaxy] )
        out_distmod.append( distanceModulus )
        out_iabs.append( iabs )
        out_rabs.append( rabs )
        out_mass_gr.append( mcMass )
        out_mass_gi.append( taMass )
        out_mass.append( fsMass )
        out_stdmass.append( std_masslight )

    out_ra = np.array(out_ra)
    out_dec = np.array(out_dec)
    out_z = np.array(out_z)
    out_id = np.array(out_id).astype(int)
    out_gr = np.array(out_gr)
    out_stdgr = np.array(out_stdgr)
    out_gi = np.array(out_gi)
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
    data = np.array([out_ra, out_dec, out_z, out_id, \
        out_gr, out_stdgr, out_gi, out_stdgi, \
        out_kri, out_stdkri, out_kii, out_stdkii, out_iobs, out_distmod, \
        out_rabs, out_iabs, out_mass_gr, out_mass_gi, out_mass, out_stdmass])
    np.savetxt(output, data.T, "%10.6f %10.5f %4.2f %5d \
        %6.3f %6.4f %6.3f %6.4f   %6.3f %6.4f %6.3f %6.4f \
         %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.4f")
    fd = open(output + ".hdr", "w")
    fd.write("# ra  dec  z id                     g-r_o std g-i_o ")
    fd.write("std   kri std kii std                  i  distmod r_abs i_abs mcMass taMass mass std\n")
    fd.close()
