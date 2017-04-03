import numpy as np
import math
import glob, os

# do them all
def doall (dir = "simha-miles-imf-alpha/", \
        alpha=False, no_alpha_dir = "simha-miles-imf/", stellar="miles", inputs="allfiles" ) :
    #Antonella updated new metallicities
    metallicites_miles = [6,12,17,20,21,22]
    metallicites_basel = [10,14,17, 20, 22] # Z=0.002, 0.0049, 0.0096, 0.019,0.03
    metallicites_parsec = [14,15]

    if stellar == "miles" :
        metallicites = metallicites_miles
    elif stellar == "basel" : 
        metallicites = metallicites_basel
    else : raise Exception("only basel and miles implemented")

    if inputs == "allfiles":
        os.chdir(dir)
        for file1 in glob.glob("*-z.mags"):
            majorMinor(dir, file1, dir, stellar)
        for file2 in glob.glob("*-0.mags"):
            majorMinor(dir, file2, dir, stellar)
    else : 
        sf_start = [1.0,3.0,5.0,7.0]
        sf_trunc = [5, 7, 9, 11, 13]
        sf_tau = [0.3, 0.7, 1.0, 1.3, 2.0, 9.0]
        sf_theta = [-1.00,-1.73,-5.66]

        # main work
        for metal in metallicites :
            for start in sf_start :
                for trunc in sf_trunc :
                    for tau in sf_tau :
                        for theta in sf_theta :
                            file = "s-" + str(metal) + "-" +str(start) + "-"
                            file = file + str(trunc) + "-" + str(tau) + str("%4.3f" % theta)
                            file1 = file + "-z.mags"
                            # working the case of alpha, 
                            # where sub-solar metallicty unaffected
                            minority_dir = dir 
                            if alpha :
                                minority_dir = no_alpha_dir
                                if (metallicites == metallicites_miles) and \
                                    (metal < 20) :
                                        dir = no_alpha_dir
                                elif (metallicites == metallicites_basel) and \
                                    (metal < 20) :
                                        dir = no_alpha_dir
                            file2 = file + "-0.mags"
                            majorMinor(dir, file1, minority_dir, stellar)
                            majorMinor(dir, file2, minority_dir, stellar)


#   Log(Z/Zsol):  0.000
#   Fraction of blue HB stars:  0.200; Ratio of BS to HB stars:  0.000
#   Shift to TP-AGB [log(Teff),log(Lbol)]:  0.00  0.00
#   IMF: 1
#   Mag Zero Point: AB (not relevant for spec/indx files)
#   SFH: Tage= 14.96 Gyr, log(tau/Gyr)=  0.954, const=  0.000, fb=  0.000, tb=  11.00 Gyr, sf_start=  1.500, dust=(  0.00,  0.00)
#   SFH 5: tau=  9.000 Gyr, sf_start=  1.500 Gyr, sf_trunc=  9.000 Gyr, sf_theta= -0.785, dust=(  0.00,  0.00)
#
#   z log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)
#20.000  5.5000 -70.0000   0.0000 -70.0000  99.000  99.000  99.000  99.000  99.000
#20.000  5.5250 -70.0000   0.0000 -70.0000  99.000  99.000  99.000  99.000  99.000

def majorMinor ( dir, majorFile, minorityDir, stellarLib="miles" ) :
    scale = .03   # we'll add a three percent by mass minority population
#    scale = .00   # we'll add no minority population
    
    pwords = majorFile.split("/")
    fname = pwords[-1]
    words = fname.split("-")
    #Antonella new low metallicity number
    minorName_miles = "s-6"
    minorName_basel = "s-4"
    if stellarLib == "miles" :
        minorName = minorName_miles
    else :
        minorName = minorName_basel
    for j in range(2,len(words)) :
        minorName = minorName + "-" + words[j]
    minorFile = ""
    for j in range(0,len(pwords)-1) :
        minorFile = minorFile + pwords[j] + "/"
    minorFile = minorFile + minorName
    outputFile = dir + majorFile.replace(".mags","m.mags")
    dummyNum1 = -70.0000
    dummyNum2 = 99.000
    maxNum1 = -60.0
    maxNum2 = 90.00
    
    majorHdr, majorData = read(dir + majorFile)
    minorHdr, minorData = read(minorityDir + minorFile)
    zed, age, mass, lum, sfr, g, r, i, z, Y = lineToNP(majorData) 
    mzed, mage, mmass, mlum, msfr, mg, mr, mi, mz, mY = lineToNP(minorData) 

    ix = np.nonzero((mY < maxNum2) & (Y < maxNum2)) # Add a check for nans & ~np.isnan(mY) & ~np.isnan(Y))
    Y[ix] = 10**(-0.4*Y[ix]) + scale*10**(-0.4*mY[ix])
    Y[ix] = -2.5*np.log10(Y[ix])
    ix = np.nonzero((mg < maxNum2) & (g < maxNum2))
    g[ix] = 10**(-0.4*g[ix]) + scale*10**(-0.4*mg[ix])
    g[ix] = -2.5*np.log10(g[ix])
    ix = np.nonzero((mr < maxNum2) & (r < maxNum2))
    r[ix] = 10**(-0.4*r[ix]) + scale*10**(-0.4*mr[ix])
    r[ix] = -2.5*np.log10(r[ix])
    ix = np.nonzero((mi < maxNum2) & (i < maxNum2))
    i[ix] = 10**(-0.4*i[ix]) + scale*10**(-0.4*mi[ix])
    i[ix] = -2.5*np.log10(i[ix])
    ix = np.nonzero((mz < maxNum2) & (z < maxNum2))
    z[ix] = 10**(-0.4*z[ix]) + scale*10**(-0.4*mz[ix])
    z[ix] = -2.5*np.log10(z[ix])

    ix = np.nonzero((mlum > maxNum1) & (lum > maxNum1))
    lum[ix] = np.log10(10**(lum[ix]) + scale*10**(mlum[ix]))
    ix = np.nonzero((msfr > maxNum1) & (sfr > maxNum1))
    sfr[ix] = np.log10(10**(sfr[ix]) + scale*10**(msfr[ix]))
    ix = np.nonzero((mmass > maxNum1) & (mass > maxNum1))
    mass[ix] = np.log10(10**(mass[ix]) + scale*10**(mmass[ix]))

    fd = open(outputFile,"w")
    fd.write(majorHdr[0])
    fd.write("#   + {:.1f}% by mass Log(Z/Zsol): {:6.3f}\n".format(
        100*scale, float(minorHdr[0][16:23])))
    for j in range(1,len(majorHdr)) :
        fd.write(majorHdr[j])
    for j in range(0,zed.size) :
        fd.write("{:6.3f} ".format( zed[j]))
        fd.write("{:7.4f} ".format(age[j]))
        fd.write("{:8.4f} ".format(mass[j]))
        fd.write("{:8.4f} ".format(lum[j]))
        fd.write("{:8.4f} ".format(sfr[j]))
        fd.write("{:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}\n".format(
            g[j],r[j],i[j],z[j],Y[j]))
    fd.close()

def read (file) :
    dummyNum = -70.0000
    header = []
    data = []
    fd = open(file, "r") 
    for line in fd :
        if line[0] == "#" :
            header.append(line)
        else :
            words = line.split()
            zed = float(words[0])
            age = float(words[1])
            mass = float(words[2])
            lum = float(words[3])
            sfr = float(words[4])
            #u = float(words[5])
            g = float(words[5])
            r = float(words[6])
            i = float(words[7])
            z = float(words[8])
            Y = float(words[9])
            if age != dummyNum :
                newLine = [zed, age, mass, lum, sfr, g, r, i, z, Y]
                data.append(newLine)
    return header, data

# convert to np.arrays
def lineToNP(data) :
    zed, age, mass, lum, sfr = [],[],[],[],[]
    g, r, i, z, Y = [],[],[],[],[]
    for j in range(0,len(data)) :
        zed.append(data[j][0])
        age.append(data[j][1])
        mass.append(data[j][2])
        lum.append(data[j][3])
        sfr.append(data[j][4])
        Y.append(data[j][9])
        #u.append(data[j][5])
        g.append(data[j][5])
        r.append(data[j][6])
        i.append(data[j][7])
        z.append(data[j][8])
    zed = np.array(zed)
    age = np.array(age)
    mass = np.array(mass)
    lum = np.array(lum)
    sfr = np.array(sfr)
    #u = np.array(u)
    g = np.array(g)
    r = np.array(r)
    i = np.array(i)
    z = np.array(z)
    Y = np.array(Y)
    return zed, age, mass, lum, sfr, g, r, i, z, Y
        


# end
