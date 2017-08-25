import CosmologicalDistance
import numpy as np

#  z     log(age) log(mass) Log(lbol) log(SFR)  u      g       r       i       z
# 0.447  9.9500   8.1393   7.8800    -70.0000  -9.667 -11.829 -13.535 -14.207 -14.599
def get_kii(zlimit=2.0) :
    dir = "/Users/annis/Code/fsps/"
    file = "simha-miles-imf/s-5-2.0-7-1.3-0.175-zm.mags"
    file = dir+file
    zed, imag = np.genfromtxt(file,unpack=True,usecols=(0,8))
    #zed, zmag = np.genfromtxt(file,unpack=True,usecols=(0,9))
    #zed, gmag = np.genfromtxt(file,unpack=True,usecols=(0,6))
    #zed, rmag = np.genfromtxt(file,unpack=True,usecols=(0,7))
    #imag = gmag-rmag
    ix = np.argsort(zed)
    zed = zed[ix]
    imag = imag[ix]
    ix = (zed > 0) & (zed <= zlimit)
    zed = zed[ix]
    imag = imag[ix]

    kii = (imag-imag[0])
    return zed,kii



def appMag(zed, iabs,mean_kii) :
    cd = CosmologicalDistance.CosmologicalDistance()
    ap = np.array([])
    for i in range (0,zed.size) :
        lumdist = cd.luminosity_distance(zed[i]) ;# in Mpc
        distanceModulus = 5*np.log10(lumdist*1e6/10.)
        iap = iabs + distanceModulus - mean_kii[i]
        ap = np.append(ap, iap)
    return ap
