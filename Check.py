import numpy as np
def do() :
    z = [1.528, 1.428, 1.331, 1.237, 1.147, 1.059, 0.974, \
        0.892, 0.812, 0.735, 0.660, 0.587, 0.516, 0.447, \
        0.380, 0.314, 0.250, 0.188, 0.127, 0.068, 0.010]
    colors = dict()
    for zed in z:
        ug,gr,ri,iz,index = doOne(zed)
        colors[zed] = ug,gr,ri,iz,index
    return colors
    
def doOne( zed) :
    ugAll, grAll, riAll, izAll = [], [] ,[] , []
    index = []
    counter = 0
    dir = "simha/"
    for metal in [10, 14, 17, 20, 22] :
        for start in [0.7, 1.0, 1.5, 2.0] :
            for trunc in [7, 9, 11, 13] :
                # need to map 0.3 < tau < 2 more carefully
                for tau in [0.3, 0.7, 1.0, 1.3, 2.0, 9.0, 13.0] :
                    for theta in [-0.175, -0.524, -0.785, -1.047, -1.396] :
                        file = "s-" + str(metal) + "-" +str(start) + "-"
                        file = file + str(trunc) + "-" + str(tau) + str(theta)
                        filename = dir + file + ".mags"
                        z, age, mass, lbol, SFR, ug, gr, ri, iz, \
                            gro, gio, kii, kri = np.genfromtxt( \
                            filename, unpack=True, skiprows=10)
                        ix = np.nonzero(zed == z)
                        ugAll.append(ug[ix])
                        grAll.append(gr[ix])
                        riAll.append(ri[ix])
                        izAll.append(iz[ix])
                        #if counter >= 75 and counter <= 100 : print counter, filename
                        index.append(counter)
                        counter += 1

    ug = np.array(ugAll)
    gr = np.array(grAll)
    ri = np.array(riAll)
    iz = np.array(izAll)
    index = np.array(index)
    return ug,gr,ri,iz,index
    
