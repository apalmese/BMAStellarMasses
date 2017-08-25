import numpy as np
import cosmolopy.distance as cd
from math import pi

def jackknife(x, func):
    """Jackknife estimate of the estimator func"""
    n = len(x)
    idx = np.arange(n)
    return np.sum(func(x[idx!=i]) for i in range(n))/float(n)

def jackknife_var(x, func):
    """Jackknife estiamte of the variance of the estimator func."""
    n = len(x)
    idx = np.arange(n)
    j_est = jackknife(x, func)
    return (n-1)/(n + 0.0) * np.sum((func(x[idx!=i]) - j_est)**2.0 for i in range(n))

def lambda_star_jk(weightedmass):
    return (weightedmass.sum())/10.**10.


def haloStellarMass ( filename="vt-stellar-masses.txt", outfile="vt-halos-stellar-masses.txt", verbose=1 ) :

    #ra,dec,z,zerr,id,central,hostid, gr,gre, gi,gie, kri,krie, kii,kiie, \
    #    i,distmod,rabs,iabs, mcMass, taMass, maiss, std = np.genfromtxt(filename, unpack=True)
    #id, gr_o,  grostd,  gi_o, giostd, kri, kristd, kii, kiistd, i,  distmod, r_abs, i_abs, mcMass, taMass, mass, std, bestsp,COADD_OBJECTS_ID_1,hostid,ra,dec,MAG_AUTO_G,MAG_AUTO_R,MAG_AUTO_I,MAG_AUTO_Z,P_RADIAL,P_REDSHIFT,GR_P_COLOR,RI_P_COLOR,IZ_P_COLOR,GR_P_MEMBER,RI_P_MEMBER,IZ_P_MEMBER,DIST_TO_CENTER,GRP_ED,GRP_BLUE,GRP_BACKGROUND,RIP_RED,RIP_BLUE,RIP_BACKGROUND,IZP_RED,IZP_BLUE,IZP_BG,MAGERR_AUTO_G,MAGERR_AUTO_R, MAGERR_AUTO_I,MAGERR_AUTO_Z,z,R200,M200,N200,LAMBDA_CHISQ,GR_SEP_FLAG,RI_SEP_FLAG,IZ_SEP_FLAG = np.genfromtxt(filename, unpack=True)

    #For xmm?
    #id,mass,std,COADD_OBJECTS_ID_1,hostid,ZP,ZPE,MAG_AUTO_G,MAG_AUTO_R,MAG_AUTO_I,MAG_AUTO_Z,P_RADIAL,P_REDSHIFT,GR_P_COLOR,RI_P_COLOR,IZ_P_COLOR,GR_P_MEMBER,RI_P_MEMBER, IZ_P_MEMBER, AMAG_R,DIST_TO_CENTER,GRP_RED,GRP_BLUE,GRP_BG,RIP_RED,RIP_BLUE, RIP_BG,IZP_RED,IZP_BLUE,IZP_BG,RESTP_RED,RESTP_BLUE,RESTP_BG,REST_P_COLOR,REST_P_MEMBER,MAGERR_AUTO_G,MAGERR_AUTO_R,MAGERR_AUTO_I,MAGERR_AUTO_Z,z,R200,M200,N200,GR_SEP_FLAG,RI_SEP_FLAG,IZ_SEP_FLAG = np.genfromtxt(filename, unpack=True)

    #For Chandra
    id, hostid, gr_o,  std,  gi_o, std, kri, std, kii, std, i,  distmod, r_abs, i_abs, mcMass, taMass, mass, std, z, sfr, sfrstd, age,agestd,bestsp, best_zmet,mean_zmet, GR_P_COLOR, RI_P_COLOR, IZ_P_COLOR, GR_P_MEMBER, RI_P_MEMBER, IZ_P_MEMBER, DIST_TO_CENTER, GRP_RED, GRP_BLUE, RIP_RED, RIP_BLUE, IZP_RED, IZP_BLUE, R200, M200 = np.genfromtxt(filename, unpack=True, delimiter=",")

# FOR X-ray cat
#    id,gr_o, std,gi_o, std,kri,std,kii,std,i,  distmod,r_abs,  i_abs,  mcMass,taMass,mass,std,sfr, sfrstd,age,agestd,bestsp,COADD_OBJECTS_ID,hostid,RA,  DEC,ZP,ZPE,  DERED_G_1,DERED_R_1,DERED_I_1,DERED_Z_1,P_RADIAL,  P_REDSHIFT, GR_P_COLOR,  RI_P_COLOR,IZ_P_COLOR,  GR_P_MEMBER, RI_P_MEMBER, IZ_P_MEMBER, AMAG_R, DIST_TO_CENTER,GRP_RED,GRP_BLUE, GRP_BG, RIP_RED,RIP_BLUE, RIP_BG, IZP_RED,IZP_BLUE, IZP_BG, RESTP_RED, RESTP_BLUE,RESTP_BG,REST_P_COLOR, REST_P_MEMBER,OBJID,  RA_1,DEC_1, ZPHOT,  ZPHOTERR,DERED_U,  DERED_G_2,DERED_R_2,DERED_I_2,DERED_Z_2,ERR_U,ERR_G, ERR_R, ERR_I, ERR_Z, ID,z,R200,  M200,N200,LAMBDA,GR_SLOPE,GR_INTERCEPT,GRMU_R,GRMU_B, GRSIGMA_R, GRSIGMA_B, GRW_R,GRW_B,  RI_SLOPE, RI_INTERCEPT,RIMU_R, RIMU_B, RISIGMA_R, RISIGMA_B,  RIW_R,  RIW_B,  GRMU_BG,  GRSIGMA_BG,GRW_BG,  RIMU_BG,  RISIGMA_BG,RIW_BG, IZ_SLOPE, IZ_INTERCEPT,IZMU_R, IZMU_B,IZSIGMA_R,  IZSIGMA_B, IZW_R,  IZW_B,  IZMU_BG,IZSIGMA_BG,IZW_BG, GR_SEP_FLAG,RI_SEP_FLAG,IZ_SEP_FLAG,REST_SLOPE,  REST_INTERCEPT,RESTMU_R,  RESTMU_B,  RESTMU_BG,RESTSIGMA_R,RESTSIGMA_B,RESTSIGMA_BG,RESTW_R,  RESTW_B,  RESTW_BG,REST_SEP_FLAG = np.genfromtxt(filename, unpack=True)

    #id,mass,std,COADD_OBJECTS_ID_1,hostid,ra,dec,ZP,ZPE,MAG_AUTO_G,MAG_AUTO_R,MAG_AUTO_I,MAG_AUTO_Z,P_RADIAL,P_REDSHIFT,GR_P_COLOR,RI_P_COLOR,IZ_P_COLOR,GR_P_MEMBER,RI_P_MEMBER,IZ_P_MEMBER,AMAG_R,DIST_TO_CENTER,GRP_RED,GRP_BLUE,GRP_BG,RIP_RED,RIP_BLUE,RIP_BG,IZP_RED,IZP_BLUE,IZP_BG,RESTP_RED,RESTP_BLUE,RESTP_BG,REST_P_COLOR,REST_P_MEMBER,MAGERR_AUTO_G,MAGERR_AUTO_R,MAGERR_AUTO_I,MAGERR_AUTO_Z,z,R200,M200,N200,GR_SEP_FLAG,RI_SEP_FLAG,IZ_SEP_FLAG,RESTW_R,RESTW_B,RESTW_BG,REST_SEP_FLAG = np.genfromtxt(filename, unpack=True)


    cosmo={'omega_M_0' : 0.23, 'omega_lambda_0' : 0.77, 'h' : 1.}
    cosmo=cd.set_omega_k_0(cosmo)

    rmax=2.9     #maximum test radius in mpc. rmax will always be included as a test radius regardless of rmin,step
    rmin=0.3   #minimum test radius in mpc. rmin always included as test radius
    step=0.1  #step size, stepping from rmin to rmax, in mpc

    id = np.array(id).astype(int)
    #central = np.array(central).astype(int)
    hostid = np.array(hostid).astype(int)    

    fd = open(outfile,"w")

    halos = np.unique(hostid)
    halos = np.sort(halos)
    if verbose: print "#halo,rad_cut, ngals, sum_mass, sum_mass_std,lambda_iz,mass_std_iz,lambda_gr_err_jk" 
    #fd.write("# hostid z  median(ra) median(dec) median(z) ngals log(stellar_mass) std lambda_gr std lambda_ri std lambda_iz std M200 izflag lambda_gr_red lambda_ri_red lambda_iz_red std lambda_gr_blue lambda_ri_blue lambda_iz_blue std (h=0.7, Om=0.3, flat)\n")
    fd.write("# halo, ngals, sum_mass, sum_mass_std, lambda_gr,lambda_gr_err_jk, lambda_ri, lambda_ri_err_jk,lambda_iz,lambda_iz_err_jk, lambda_gr_red, lambda_ri_red, lambda_iz_red,lambda_iz_red_err_jk, lambda_gr_blue, lambda_ri_blue, lambda_iz_blue,lambda_iz_blue_err_jk\n")

    

    for halo in halos:
        ix_cl = np.nonzero(hostid == halo)
        z_cl = z[ix_cl][0] 
        r200_cl = R200[ix_cl][0]
        #print DIST_TO_CENTER[ix_cl]
        #ang_diam_dist=cd.angular_diameter_distance(z_cl,z0=0,**cosmo)
        #distmpc = ang_diam_dist*DIST_TO_CENTER[ix_cl]*pi/180.
        distmpc = DIST_TO_CENTER[ix_cl]
    
        mass_cl = mass[ix_cl]
        GR_P_MEMBER_cl = GR_P_MEMBER[ix_cl]
        RI_P_MEMBER_cl = RI_P_MEMBER[ix_cl]
        IZ_P_MEMBER_cl = IZ_P_MEMBER[ix_cl]
        std_cl = std[ix_cl]
        GR_P_COLOR_cl = GR_P_COLOR[ix_cl]
        RI_P_COLOR_cl = RI_P_COLOR[ix_cl]
        IZ_P_COLOR_cl = IZ_P_COLOR[ix_cl]
        GRP_RED_cl = GRP_RED[ix_cl]
        RIP_RED_cl = RIP_RED[ix_cl]
        IZP_RED_cl = IZP_RED[ix_cl]
        GRP_BLUE_cl = GRP_BLUE[ix_cl]
        RIP_BLUE_cl = RIP_BLUE[ix_cl]
        IZP_BLUE_cl = IZP_BLUE[ix_cl]
#        P_RADIAL_cl = P_RADIAL[ix_cl]

        radii=np.r_[rmin:rmax:step,rmax]
        radii = np.append(radii,r200_cl)
        #radii = [0.2,r200_cl]

        for rad_cut in radii:
            ix = np.nonzero(distmpc <=rad_cut)
            
        #zmedian = np.median(z[ix])
        #ramedian = np.median(ra[ix])
        #decmedian = np.median(dec[ix])
            ngals = id[ix].size
            if ngals>0.:
                linear_mass = 10**mass_cl[ix]
                linear_mass_weight_gr = 10**mass_cl[ix]*GR_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                linear_mass_weight_ri = 10**mass_cl[ix]*RI_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                linear_mass_weight_iz = 10**mass_cl[ix]*IZ_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                mass_errors = np.log(10.)*linear_mass*std_cl[ix]
                mass_std = np.sqrt((mass_errors**2).sum())

                mass_err_gr = np.log(10.)*linear_mass*std_cl[ix]*GR_P_MEMBER_cl[ix]
                mass_err_ri = np.log(10.)*linear_mass*std_cl[ix]*RI_P_MEMBER_cl[ix]
                mass_err_iz = np.log(10.)*linear_mass*std_cl[ix]*IZ_P_MEMBER_cl[ix]
                mass_std_gr = 10**(-10)*np.sqrt((mass_err_gr**2).sum())
                mass_std_ri = 10**(-10)*np.sqrt((mass_err_ri**2).sum())
                mass_std_iz = 10**(-10)*np.sqrt((mass_err_iz**2).sum())

                linear_mass_weight_gr_red = 10**mass_cl[ix]*GR_P_MEMBER_cl[ix]/GR_P_COLOR_cl[ix]*GRP_RED_cl[ix]
                linear_mass_weight_ri_red = 10**mass_cl[ix]*RI_P_MEMBER_cl[ix]/RI_P_COLOR_cl[ix]*RIP_RED_cl[ix]
                linear_mass_weight_iz_red = 10**mass_cl[ix]*IZ_P_MEMBER_cl[ix]/IZ_P_COLOR_cl[ix]*IZP_RED_cl[ix]#/P_RADIAL_cl[ix]

                linear_mass_weight_gr_blue = 10**mass_cl[ix]*GR_P_MEMBER_cl[ix]/GR_P_COLOR_cl[ix]*GRP_BLUE_cl[ix]
                linear_mass_weight_ri_blue = 10**mass_cl[ix]*RI_P_MEMBER_cl[ix]/RI_P_COLOR_cl[ix]*RIP_BLUE_cl[ix]
                linear_mass_weight_iz_blue = 10**mass_cl[ix]*IZ_P_MEMBER_cl[ix]/IZ_P_COLOR_cl[ix]*IZP_BLUE_cl[ix]

                #jackknife for errors on lambda_star
                weightmass_gr = 10**mass_cl[ix]*GR_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                lambda_gr_err_jk = (jackknife_var(weightmass_gr, lambda_star_jk))**0.5
                weightmass_ri = 10**mass_cl[ix]*RI_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                lambda_ri_err_jk = (jackknife_var(weightmass_ri, lambda_star_jk))**0.5
                weightmass_iz = 10**mass_cl[ix]*IZ_P_MEMBER_cl[ix] #/P_RADIAL_cl[ix]
                lambda_iz_err_jk = (jackknife_var(weightmass_iz, lambda_star_jk))**0.5
                weightmass_iz_red = 10**mass_cl[ix]*IZ_P_MEMBER_cl[ix]/IZ_P_COLOR_cl[ix]*IZP_RED_cl[ix]
                lambda_iz_red_err_jk = (jackknife_var(weightmass_iz_red, lambda_star_jk))**0.5
                weightmass_iz_blue = 10**mass[ix]*IZ_P_MEMBER_cl[ix]/IZ_P_COLOR_cl[ix]*IZP_BLUE_cl[ix]
                lambda_iz_blue_err_jk = (jackknife_var(weightmass_iz_blue, lambda_star_jk))**0.5

                sum_mass = linear_mass.sum()
                lambda_gr = (linear_mass_weight_gr.sum())/10.**10.
                lambda_ri = (linear_mass_weight_ri.sum())/10.**10.
                lambda_iz = (linear_mass_weight_iz.sum())/10.**10.
                sum_mass_std = mass_std/sum_mass/np.log(10.)
                sum_mass = np.log10(sum_mass)

                lambda_gr_red = (linear_mass_weight_gr_red.sum())/10.**10.
                lambda_ri_red = (linear_mass_weight_ri_red.sum())/10.**10.
                lambda_iz_red = (linear_mass_weight_iz_red.sum())/10.**10.

                lambda_gr_blue = (linear_mass_weight_gr_blue.sum())/10.**10.
                lambda_ri_blue = (linear_mass_weight_ri_blue.sum())/10.**10.
                lambda_iz_blue = (linear_mass_weight_iz_blue.sum())/10.**10.

        #sum_mass_gr = np.log10(sum_mass_gr)
        #sum_mass_ri = np.log10(sum_mass_ri)
        #sum_mass_iz = np.log10(sum_mass_iz)
        #lambda_rm = LAMBDA_CHISQ[ix[0][0]]

        #TO UNCOMMENT!!!!!!!
        #M200_GMM = M200[ix[0][0]]
        #zout = z[ix[0][0]]
        #iz_flag = IZ_SEP_FLAG[ix[0][0]]
                if verbose: 
                    print "{:10d}  {:6.3f} {:4d}      {:6.3f} {:6.4f} {:6.3f} {:6.4f} {:6.4f}".format(halo,rad_cut, ngals, sum_mass, sum_mass_std,lambda_iz,mass_std_iz,lambda_gr_err_jk)

                fd.write("{:10d} {:6.3f} {:4d} {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f}  {:6.3f} {:6.3f} {:6.3f} {:6.3f}  {:6.3f} {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n".format(halo, rad_cut, ngals, sum_mass, sum_mass_std, lambda_gr,lambda_gr_err_jk, lambda_ri, lambda_ri_err_jk,lambda_iz,lambda_iz_err_jk, lambda_gr_red, lambda_ri_red, lambda_iz_red,lambda_iz_red_err_jk, lambda_gr_blue, lambda_ri_blue, lambda_iz_blue,lambda_iz_blue_err_jk))
        #LATEST:
        #fd.write("{:10d} {:6.3f}  {:6.3f}     {:4d}      {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f}  {:6.3f} {:6.3f} {:6.3f} {:6.3f}  {:6.3f} {:6.3f} {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:2.3f}\n".format(halo,zout, zmedian, ngals, sum_mass, sum_mass_std, lambda_gr,lambda_gr_err_jk, lambda_ri, lambda_ri_err_jk,lambda_iz,lambda_iz_err_jk,M200_GMM, lambda_gr_red, lambda_ri_red, lambda_iz_red,lambda_iz_red_err_jk, lambda_gr_blue, lambda_ri_blue, lambda_iz_blue,lambda_iz_blue_err_jk,iz_flag))
        #fd.write("{:10d} {:6.3f}  {:6.3f}     {:4d}      {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}  {:6.3f} {:6.3f} {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:10d}\n".format(halo,zout, zmedian, ngals, sum_mass, sum_mass_std, lambda_gr,mass_err_gr, lambda_ri, mass_err_ri,lambda_iz,mass_err_iz,M200_GMM, lambda_gr_red, lambda_ri_red, lambda_iz_red, lambda_gr_blue, lambda_ri_blue, lambda_iz_blue,iz_flag))
        #fd.write("{:10d} {:6.3f}  {:6.3f}     {:4d}      {:6.3f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:3.0f} {:6.4f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n".format(halo,zout, zmedian, ngals, sum_mass, sum_mass_std, lambda_gr, lambda_ri, lambda_iz,mass_err_iz,M200_GMM,iz_flag, lambda_gr_red, lambda_ri_red, lambda_iz_red, lambda_gr_blue, lambda_ri_blue, lambda_iz_blue))

    fd.close()
# cat vt-halos-stellar-masses.txt | sort -k6 -nr > vtt
