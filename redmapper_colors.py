import numpy as np
import color_tracks
import matplotlib.pyplot as plt; 

# degenerate on g-r: simha-miles
# s-5-2.0-7-1.3-0.175-zm.mags   
# s-4-0.7-7-0.3-0.175-zm.mags playing off sf_start and sf_tau with metals
# on g-r, r-i plot only this one gets turnup right, but isnt good on g-r
# s-4-0.7-7-1.3-0.175-zm.mag 

#file="/Users/annis/Code/jta_outputs_fsps/s-20-0.7-11-0.3-0.785-z.mags"
# file = "/Users/annis/Code/fsps/simha-miles/s-20-0.7-11-0.3-0.785-z.mags"
# file = "/Users/annis/Code/fsps/simha-miles-imf/s-20-0.7-11-0.3-0.785-z.mags"

# rm_ra,rm_dec, rm_zed, rm_zspec, rm_imag,rm_ug,rm_gr,rm_ri,rm_iz,rm_imagerr, rm_ugerr, rm_grerr,rm_rierr,rm_izerr = redmapper_colors.get()
# ix=(rm_gr>0)&(rm_gr<2.0);plt.clf();plt.hexbin(rm_zed[ix], rm_gr[ix])  
# reload(redmapper_colors);redmapper_colors.plot_simha(dir+"/s-3-2.0-7-1.3-0.175-zm.mags"); plt.xlim(0.09,0.46);plt.ylim(0.9,1.95)
#  r-i : plt.xlim(0.09,0.46);plt.ylim(0.38,0.9)
# g-r, r-i: ;plt.xlim(0.9,1.95) ;plt.ylim(0.38,0.9)
def plot_simha(file="/Users/annis/Code/jta_outputs_fsps/s-20-0.7-11-0.3-0.785-z.mags", zlow=0.06, zhi=0.6) :
    zed,age,mass,lbol,sfr,u,g,r,i,z = \
        np.genfromtxt(file, unpack=True,skiprows=9);
    ug=u-g;gr=g-r;ri=r-i;iz=i-z; 
    ix2 = (i < 99)&(zed<5.) & (zed<zhi )&(zed>zlow)

    zlist=np.array([0.1, 0.15,0.2, 0.25, 0.3, 0.35, 0.4, 0.45]); 
    zlist=np.array([0.075, 0.1, 0.15,0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,1.1,1.2]); 
    nzed, nug, ngr, nri, niz = color_tracks.get_sdss_zpoints(
            zlist,zed[ix2], ug[ix2],gr[ix2],ri[ix2],iz[ix2]);
    #return nzed, nug, ngr, nri, niz

    # g-r vs z
    #plt.plot(zed[ix2],gr[ix2],c="w");
    #plt.xlabel("z");plt.ylabel("g-r"); 
    #plt.scatter(nzed,ngr,c=nzed,s=100); 
    #x1 = nzed; x2 = ngr; t1="zed"; t2="gr"
    #plt.xlim(0.09,0.47);plt.ylim(0.9,1.95)
    # r-i vs z
    plt.plot(zed[ix2],ri[ix2],c="w");
    plt.xlabel("z");plt.ylabel("r-i"); 
    plt.scatter(nzed,nri,c=nzed,s=100); 
    x1 = nzed; x2 = nri; t1="zed"; t2="ri"
    plt.xlim(0.09,0.46);plt.ylim(0.38,0.9)
    # g-r vs r-i
    #del_ri = 0.03
    #del_ri=0.0
    #plt.plot(gr[ix2],ri[ix2]+del_ri,c="w");
    #plt.xlabel("g-r");plt.ylabel("r-i"); 
    #plt.scatter(ngr,nri+del_ri,c=nzed,s=100); 
    #x1 = ngr; x2 = nri+del_ri; t1="gr"; t2="ri"
    #plt.xlim(0.9,1.95) ;plt.ylim(0.38,0.9)
    # u-g vs z
    #plt.plot(zed[ix2],ug[ix2],c="w");
    #plt.xlabel("z");plt.ylabel("u-g"); 
    #plt.scatter(nzed,nug,c=nzed,s=100); 
    # u-r vs r-z  (not useful- sn problems?)
    #plt.plot( (ri[ix2]-iz[ix2]), (ug[ix2]-gr[ix2]), c="w")
    #plt.ylabel("u-r");plt.xlabel("r-z"); 
    #plt.scatter( (nri-niz), (nug-ngr) ,c=nzed,s=100); 
    # g-r vs r-z
    #plt.plot( (ri[ix2]+iz[ix2]), gr[ix2], c="w")
    #plt.ylabel("g-r");plt.xlabel("r-z"); 
    #plt.scatter( (nri+niz), ngr ,c=nzed,s=100); 
    #x1 = ngr; x2 = nri+niz; t1="gr"; t2="rz"
    # plt.xlim(0.3,1.7);plt.ylim(0.9,1.9)
    # i vs dperp
    #plt.plot( (i[ix2], ri[ix2]-gr[ix2]/8., c="w")
    #plt.xlabel("i");plt.xlabel("r-i - (g-r)/8."); 
    #plt.scatter( (ni), ngr ,c=nzed,s=100); 
    print "{:s}=".format(t1),  
    for x in x1 :
        print "{:.2f} ".format(x),
    print "\n {:s}=".format(t2),
    for x in x2 :
        print "{:.2f} ".format(x),
    print ""



# ra,dec, zed, zspec, imag,ug,gr,ri,iz,imagerr, ugerr, grerr,rierr,izerr = redmapper_colors.get()
# rm_ra,rm_dec, rm_zed, rm_zspec, rm_imag,rm_ug,rm_gr,rm_ri,rm_iz,rm_imagerr, rm_ugerr, rm_grerr,rm_rierr,rm_izerr = redmapper_colors.get()
def get() :
    import pyfits
    rm_data="/Users/annis/Documents/Annis/Papers/redmapper/"
    rm_data = rm_data + "dr8_run_redmapper_v6.3.1_redmagic_0.5-10_wdr12.fit"
    hdu=pyfits.open(rm_data)
    td=hdu[1].data
    ra = td["ra"]
    dec = td["dec"]
    flag = td["redmagicflag"]
    imag = td["imag"]
    imagerr = td["imag_err"]
    model_mag = td["model_mag"]
    model_magerr = td["model_magerr"]
    u = model_mag[:,0]
    uerr = model_magerr[:,0]
    g = model_mag[:,1]
    gerr = model_magerr[:,1]
    r = model_mag[:,2]
    rerr = model_magerr[:,2]
    i = model_mag[:,3]
    ierr = model_magerr[:,3]
    z = model_mag[:,4]
    zerr = model_magerr[:,4]
    ug = u-g
    gr = g-r
    ri = r-i
    iz = i-z
    ugerr = np.sqrt(uerr**2+gerr**2)
    grerr = np.sqrt(gerr**2+rerr**2)
    rierr = np.sqrt(rerr**2+ierr**2)
    izerr = np.sqrt(ierr**2+zerr**2)
    zed = td["zredmagic"]
    zspec = td["zspec"]
    return ra,dec, zed, zspec, imag,ug,gr,ri,iz,imagerr, ugerr, grerr,rierr,izerr

def get_id (zlow, zhi, filename="") :
    import pyfits
    rm_data="/Users/annis/Documents/Annis/Papers/redmapper/"
    rm_data = rm_data + "dr8_run_redmapper_v6.3.1_redmagic_0.5-10_wdr12.fit"
    hdu=pyfits.open(rm_data)
    td=hdu[1].data
    id = td["id"]
    ra = td["ra"]
    dec = td["dec"]
    flag = td["redmagicflag"]
    photoid = td["photoid"]
    zspec = td["zspec"]
    imag = td["imag"]
    ix = (zspec >= zlow)&(zspec < zhi)&(imag<17.5)
    ra = ra[ix]
    dec = dec[ix]
    id=id[ix]
    flag=flag[ix]
    photoid=photoid[ix]
    zspec=zspec[ix]
    if filename != "" :
        data = ra,dec,id,flag,photoid,zspec
        data = np.array([data]).T
        np.savetxt(filename, data, "%.6f %.5f %d %d %d %.3f")
    return ra,dec, id, flag, photoid, zspec

# so, no reddening correction?
# ra,dec,zed, imag, ug,gr,ri,iz, imagerr, ugerr,grerr,rierr,izerr = redmapper_colors.get()
def get2 () :
    rm_data="/Users/annis/Documents/Annis/Papers/redmapper/maxColors_annis.csv"
    id,objID,ra,dec,z,extinction_i,cmodelMag_i,petromag_r,\
            model_ug,model_gr,model_ri,model_iz,modelMagErr_i, \
            model_ug_err,model_gr_err,model_ri_err,model_iz_err = \
            np.genfromtxt(rm_data,unpack=True,delimiter=",",skiprows=1)
    id,objID,ra,dec,z,extinction_i,cmodelMag_i,petromag_r,\
            model_ug,model_gr,model_ri,model_iz,modelMagErr_i, \
            model_ug_err,model_gr_err,model_ri_err,model_iz_err 
    zed = z
    imag = cmodelMag_i
    ug = model_ug
    gr = model_gr
    ri = model_ri
    iz = model_iz
    imag_err = modelMagErr_i 
    ugerr = model_ug_err
    grerr = model_gr_err
    rierr = model_ri_err
    izerr = model_iz_err
    return  ra,dec,zed, imag, ug,gr,ri,iz, imag_err, ugerr,grerr,rierr,izerr
# rm_ra,rm_dec,rm_zed, rm_imag, rm_ug,rm_gr,rm_ri,rm_iz, rm_imagerr, rm_ugerr,rm_grerr,rm_rierr,rm_izerr = redmapper_colors.get()


