import numpy as np
def go(mass,z) :
    for zed in range(10,55,5) :
        zed = zed/100.
        zed2 = zed + 0.05
        ix = np.nonzero( (z >= zed) & (z < zed2))
        m = mass[ix].mean()
        ms = mass[ix].std()/np.sqrt(mass[ix].size)
        dm = 11.54 - m
        print "{:0.2f} {:.3f} {:.4f} {:7.3f}".format(zed,m,ms,dm)
