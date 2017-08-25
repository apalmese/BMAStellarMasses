import numpy as np
import loadPopColors
import helperfunctions
#import weightedstats as ws
import pyfits as pf
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#
# inputDataDict is what holds the input photometry
# it expects id, i, ierr, gr,ri,iz, grerr, rierr, izerr
#
#
def calc (inputDataDict, indir="simha/", lib="miles") :
    from scipy.stats import chi2
    import CosmologicalDistance
    import os
    cd = CosmologicalDistance.CosmologicalDistance()
    ldistDict = dict()
    splineDict = dict()
    splines, zmet = loadPopColors.doall(indir, lib=lib)
    id  = inputDataDict["id"]
    haloid  = inputDataDict["haloid"]
    i  = inputDataDict["i"]
    ierr  = inputDataDict["ierr"]
    gr  = inputDataDict["gr"]
    ri  = inputDataDict["ri"]
    iz  = inputDataDict["iz"]
    grerr  = inputDataDict["grerr"]
    rierr  = inputDataDict["rierr"]
    izerr  = inputDataDict["izerr"]
    allzed = inputDataDict["zed"]
    GR_P_COLOR=inputDataDict["GR_P_COLOR"]
    RI_P_COLOR=inputDataDict["RI_P_COLOR"]
    IZ_P_COLOR=inputDataDict["IZ_P_COLOR"]
    GR_P_MEMBER=inputDataDict["GR_P_MEMBER"]
    RI_P_MEMBER=inputDataDict["RI_P_MEMBER"]
    IZ_P_MEMBER=inputDataDict["IZ_P_MEMBER"]
    DIST_TO_CENTER=inputDataDict["DIST_TO_CENTER"]
    GRP_RED=inputDataDict["GRP_RED"]
    GRP_BLUE=inputDataDict["GRP_BLUE"]
    RIP_RED=inputDataDict["RIP_RED"]
    RIP_BLUE=inputDataDict["RIP_BLUE"]
    IZP_RED=inputDataDict["IZP_RED"]
    IZP_BLUE=inputDataDict["IZP_BLUE"]

    # protect against too small of errors => values = 0
    ix = np.nonzero(grerr < 0.02)
    grerr[ix] = 0.02
    ix = np.nonzero(rierr < 0.02)
    rierr[ix] = 0.02
    ix = np.nonzero(izerr < 0.02)
    izerr[ix] = 0.02

    mod_tab = np.genfromtxt('simha_miles_Nov16_models_lookup.tab')

    t_start = mod_tab[:,2]

    # prepping for output
    out_id,out_haloid, out_gr, out_stdgr, out_gi, out_stdgi, \
        out_kri, out_stdkri, out_kii, out_stdkii, out_iobs, out_distmod, \
        out_bestzmet,out_stdzmet, out_rabs, out_iabs, out_mass_gr, out_mass_gi, out_mass, out_stdmass, out_sfr, out_stdsfr, out_age, out_stdage, \
        out_GR_P_COLOR,out_RI_P_COLOR,out_IZ_P_COLOR,out_GR_P_MEMBER,out_RI_P_MEMBER,out_IZ_P_MEMBER,out_DIST_TO_CENTER, \
        out_GRP_RED,out_GRP_BLUE,out_RIP_RED,out_RIP_BLUE,out_IZP_RED,out_IZP_BLUE, out_zmet, out_zed =\
        [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    out_bestsp, out_bestchisq = [], []

    size = id.size
    bins = 10
    colors = iter(cm.rainbow(np.linspace(0, 1, bins)))
    plt.clf()
    #size = 10
    for z in np.linspace(0.1,1.,num=bins) :
        zed = z
        print z #galaxy,"  of  ",size,"    z = ",zed

        rest_gr,rest_gi,gr,ri,iz = [], [],[],[],[]
        masslight, sfrs,ages,ages_cosmic,ages_cut,zmets,kii, kri, models = [],[],[],[],[],[],[],[],[]
        minChiSq = 999; spIndex = -1
        for sp in range(0,len(splines)) :
            # for speed
            skey = str(sp) + "-" + str(zed)
            #if skey in splineDict :
            #    sgr,sri,siz,sgrr,sgir,skii,skri,sml = splineDict[skey]
            #else :
            #    sgr = splines[sp][0](zed)
            #    sri = splines[sp][1](zed)
            #    siz = splines[sp][2](zed)
            #    sgrr = splines[sp][4](zed) ;# restframe g-r
            #    sgir = splines[sp][5](zed) ;# restframe g-i
            #    skii = splines[sp][6](zed) ;# kcorrection: i_o - i_obs
            #    skri = splines[sp][7](zed) ;# kcorrection: r_o - i_obs
            #    sml = splines[sp][8](zed) ;# log(mass/light)  (M_sun/L_sun)
            #    ssfr = splines[sp][9](zed)
            #    sage_cosmic = splines[sp][10](zed)
            #    sage = splines[sp][11](zed)
            #    zmet = splines[12]
            sgr = splines[sp][0](zed)
            sri = splines[sp][1](zed)
            siz = splines[sp][2](zed)
            sgrr = splines[sp][4](zed) ;# restframe g-r
            sgir = splines[sp][5](zed) ;# restframe g-i
            skii = splines[sp][6](zed) ;# kcorrection: i_o - i_obs
            skri = splines[sp][7](zed) ;# kcorrection: r_o - i_obs
            sml = splines[sp][8](zed) ;# log(mass/light)  (M_sun/L_sun)
            ssfr = splines[sp][9](zed)
            if (ssfr<-20.): ssfr=-20.
            sage_cosmic = splines[sp][10](zed)
            sage = splines[sp][11](zed)
            szmet = zmet[sp]
                #To be changed if SFH changes

            #splineDict[skey] = sgr,sri,siz,sgrr,sgir,skii,skri,sml

            #print sage_cosmic
            gr.append(sgr)
            ri.append(sri)
            iz.append(siz)
            models.append(sp)
            rest_gr.append(sgrr)
            rest_gi.append(sgir)
            kii.append(skii)
            kri.append(skri)
            masslight.append(sml)
            sfrs.append(ssfr)
            ages.append(sage)
            ages_cosmic.append(sage_cosmic)
            zmets.append(szmet)
            ages_cut.append(10.**(sage_cosmic-9.)-t_start[sp])
        #plt.scatter(masslight,sfrs,edgecolor='none',alpha=0.7)
        #plt.savefig('ssfr_vs_ml_z'+str(zed)+'.png')
        #plt.clf()
        print ages_cut
        ages_cut = np.array(ages_cut)
        sfrs=np.array(sfrs)
        ages_cosmic=np.array(ages_cosmic)
        idx = (ages_cut > 0.1)
        plt.scatter(sfrs[idx],ages_cosmic[idx],edgecolor='none',alpha=0.7, label = str(zed),color = next(colors))
        #plt.savefig('restgr_vs_ssfr_z'+str(zed)+'.png')
        #plt.clf()
        #plt.scatter(ri,gr,edgecolor='none',alpha=0.7)
        #plt.savefig('gr_vs_ri_z'+str(zed)+'.png')
        #plt.clf()
    plt.legend()
    plt.xlabel('ssfr')
    plt.ylabel('age')
    plt.savefig('age_vs_ssfr_new.png')
