import numpy as np
from astropy.io import fits

def read_fits(filename):
    h=fits.open(filename)
    d=h[1].data
    id = d['COADD_OBJECTS_ID']
    haloid = d['HOST_HALOID']
    GR_P_COLOR=d['GR_P_COLOR']
    RI_P_COLOR=d['RI_P_COLOR']
    IZ_P_COLOR=d['IZ_P_COLOR']
    GR_P_MEMBER=d['GR_P_MEMBER']
    RI_P_MEMBER=d['RI_P_MEMBER']
    IZ_P_MEMBER=d['IZ_P_MEMBER']
    DIST_TO_CENTER=d['DIST_TO_CENTER']
    GRP_RED=d['GRP_RED']
    GRP_BLUE=d['GRP_BLUE']
    RIP_RED=d['RIP_RED']
    RIP_BLUE=d['RIP_BLUE']
    IZP_RED=d['IZP_RED']
    IZP_BLUE=d['IZP_BLUE']
    g = d['MAG_AUTO_G']
    gerr = d['MAGERR_AUTO_G']
    r = d['MAG_AUTO_R']
    rerr = d['MAGERR_AUTO_R']
    i = d['MAG_AUTO_I']
    ierr = d['MAGERR_AUTO_I']
    z = d['MAG_AUTO_Z']
    zerr = d['MAGERR_AUTO_Z']
    zed = d['Z']
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5

    inputDataDict = {'id':id,'haloid':haloid,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed, \
        'GR_P_COLOR':GR_P_COLOR,'RI_P_COLOR':RI_P_COLOR,'IZ_P_COLOR':IZ_P_COLOR,'GR_P_MEMBER':GR_P_MEMBER,'RI_P_MEMBER':RI_P_MEMBER,'IZ_P_MEMBER':IZ_P_MEMBER,\
        'DIST_TO_CENTER':DIST_TO_CENTER,'GRP_RED':GRP_RED,'GRP_BLUE':GRP_BLUE,'RIP_RED':RIP_RED,'RIP_BLUE':RIP_BLUE,'IZP_RED':IZP_RED,'IZP_BLUE':IZP_BLUE}

    return inputDataDict




def read_fits_sdss(filename):
    h=fits.open(filename)
    d=h[1].data
    id = d['COADD_OBJECTS_ID']
    haloid = d['HOST_HALOID']
    GR_P_COLOR=d['GR_P_COLOR']
    RI_P_COLOR=d['RI_P_COLOR']
    IZ_P_COLOR=d['IZ_P_COLOR']
    GR_P_MEMBER=d['GR_P_MEMBER']
    RI_P_MEMBER=d['RI_P_MEMBER']
    IZ_P_MEMBER=d['IZ_P_MEMBER']
    DIST_TO_CENTER=d['DIST_TO_CENTER']
    GRP_RED=d['GRP_RED']
    GRP_BLUE=d['GRP_BLUE']
    RIP_RED=d['RIP_RED']
    RIP_BLUE=d['RIP_BLUE']
    IZP_RED=d['IZP_RED']
    IZP_BLUE=d['IZP_BLUE']
    u = d['DERED_U']
    uerr = d['ERR_U']
    g = d['DERED_G_1']
    gerr = d['ERR_G']
    r = d['DERED_R_1']
    rerr = d['ERR_R']
    i = d['DERED_I_1']
    ierr = d['ERR_I']
    z = d['DERED_Z_1']
    zerr = d['ERR_Z']
    zed = d['Z']
    ugerr = (gerr**2+uerr**2)**0.5
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'haloid':haloid,'ug':u-g,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'ugerr':ugerr,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed, \
        'GR_P_COLOR':GR_P_COLOR,'RI_P_COLOR':RI_P_COLOR,'IZ_P_COLOR':IZ_P_COLOR,'GR_P_MEMBER':GR_P_MEMBER,'RI_P_MEMBER':RI_P_MEMBER,'IZ_P_MEMBER':IZ_P_MEMBER,\
        'DIST_TO_CENTER':DIST_TO_CENTER,'GRP_RED':GRP_RED,'GRP_BLUE':GRP_BLUE,'RIP_RED':RIP_RED,'RIP_BLUE':RIP_BLUE,'IZP_RED':IZP_RED,'IZP_BLUE':IZP_BLUE}
         
    return inputDataDict

def read_fits_smassonly(filename):
    h=fits.open(filename)
    d=h[1].data
    id = d['ID']
    g = d['MAG_G']
    gerr = d['MAGERR_G']
    r = d['MAG_R']
    rerr = d['MAGERR_R']
    i = d['MAG_I']
    ierr = d['MAGERR_I']
    z = d['MAG_Z']
    zerr = d['MAGERR_Z']
    zed = d['Z']
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict


def read_fits_redmagic(filename):
    h=fits.open(filename)
    d=h[1].data
    id = d['COADD_OBJECTS_ID']
    g = d['MODEL_MAG'][:,0]
    gerr =  d['MODEL_MAGERR'][:,0]
    r =  d['MODEL_MAG'][:,1]
    rerr = d['MODEL_MAGERR'][:,1]
    i =  d['MODEL_MAG'][:,2]
    ierr = d['MODEL_MAGERR'][:,2]
    z =  d['MODEL_MAG'][:,3]
    zerr = d['MODEL_MAGERR'][:,3]
    zed = d['ZREDMAGIC']
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict


def read(filename):
    d = np.genfromtxt(filename)
    id = d[:,0]
    g = d[:,1]
    gerr = d[:,2]
    r = d[:,3]
    rerr = d[:,4]
    i = d[:,5]
    ierr = d[:,6]
    z = d[:,7]
    zerr = d[:,8]
    zed = d[:,12]
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict


def read_GMM(filename):
    d = np.genfromtxt(filename)
    id = d[:,0]
    g = d[:,4]
    gerr = d[:,26]
    r = d[:,5]
    rerr = d[:,27]
    i = d[:,6]
    ierr = d[:,28]
    z = d[:,7]
    zerr = d[:,29]
    zed = d[:,30]
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    haloid = d[:,1]
    GR_P_COLOR=d[:,10]
    RI_P_COLOR=d[:,11]
    IZ_P_COLOR=d[:,12]
    GR_P_MEMBER=d[:,13]
    RI_P_MEMBER=d[:,14]
    IZ_P_MEMBER=d[:,15]
    DIST_TO_CENTER=d[:,16]
    GRP_RED=d[:,17]
    GRP_BLUE=d[:,18]
    RIP_RED=d[:,20]
    RIP_BLUE=d[:,21]
    IZP_RED=d[:,23]
    IZP_BLUE=d[:,24]
    inputDataDict = {'id':id,'haloid':haloid,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed, \
        'GR_P_COLOR':GR_P_COLOR,'RI_P_COLOR':RI_P_COLOR,'IZ_P_COLOR':IZ_P_COLOR,'GR_P_MEMBER':GR_P_MEMBER,'RI_P_MEMBER':RI_P_MEMBER,'IZ_P_MEMBER':IZ_P_MEMBER,\
        'DIST_TO_CENTER':DIST_TO_CENTER,'GRP_RED':GRP_RED,'GRP_BLUE':GRP_BLUE,'RIP_RED':RIP_RED,'RIP_BLUE':RIP_BLUE,'IZP_RED':IZP_RED,'IZP_BLUE':IZP_BLUE}

    #inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict

def read_delucia(filename):
    d = np.genfromtxt(filename)
    id = d[:,0]
    g = d[:,8]
    gerr = d[:,12]
    r = d[:,9]
    rerr = d[:,13]
    i = d[:,10]
    ierr = d[:,14]
    z = d[:,11]
    zerr = d[:,15]
    zed = d[:,4]
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict

def read_xmm(filename):
    d = np.genfromtxt(filename)
    id = d[:,0]
    g = d[:,4]
    gerr = d[:,32]
    r = d[:,5]
    rerr = d[:,33]
    i = d[:,6]
    ierr = d[:,34]
    z = d[:,7]
    zerr = d[:,35]
    zed = d[:,36]
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict

	
def read_chandra(filename):
    d = np.genfromtxt(filename)
    id = d[:,0]
    g = d[:,6]
    gerr = d[:,34]
    r = d[:,7]
    rerr = d[:,35]
    i = d[:,8]
    ierr = d[:,36]
    z = d[:,9]
    zerr = d[:,37]
    zed = d[:,38]
    grerr = (gerr**2+rerr**2)**0.5
    rierr = (rerr**2+ierr**2)**0.5
    izerr = (ierr**2+zerr**2)**0.5
    inputDataDict = {'id':id,'i':i,'ierr':ierr,'gr':g-r,'ri':r-i,'iz':i-z,'grerr':grerr,'rierr':rierr,'izerr':izerr,'zed':zed}

    return inputDataDict
