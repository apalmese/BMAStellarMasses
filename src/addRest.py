import numpy as np
import CosmologicalDistance
import glob, os

# do them all
def doall ( dir = "simha-miles-imf-alpha/", stellar="miles", inputs="allfiles", noMinority="False") :
    metallicites_basel = [10,14,17, 20, 22]
    metallicites_miles = [12,17,20,21,22]
    if stellar == "miles" :
        metallicites = metallicites_miles
    elif stellar == "basel" :
        metallicites = metallicites_basel
    else : raise Exception("only basel and miles implemented")

    if inputs == "allfiles":
        os.chdir(dir)
        if noMinority == "True":
            print noMinority
            for file in glob.glob("*-z.mags"):
                file = dir + file
                addRestframeGR (file, noMinority="True" )
        else:
            for file in glob.glob("*-zm.mags"):
                file = dir + file
                addRestframeGR (file )
    else :

        sf_start = [1.0,3.0,5.0, 7.0]
        sf_trunc = [5, 7, 9, 11, 13]
        sf_tau = [0.3, 0.7, 1.0, 1.3, 2.0, 9.0]
        sf_theta = [-1.00,-1.73,-5.66]

        for metal in metallicites :
            for start in sf_start :
                for trunc in sf_trunc :
                    for tau in sf_tau :
                        for theta in sf_theta :
                            file = "s-" + str(metal) + "-" +str(start) + "-"
                            file = file + str(trunc) + "-" + str(tau) + str("%4.3f" % theta)
                            file = dir + file + "-zm.mags"
                            addRestframeGR (file ) 


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

def addRestframeGR ( filename , noMinority="False" ) :
    cd = CosmologicalDistance.CosmologicalDistance()
    # filename = s-4-2.0-7-1.3-0.175-zm.mags
    # rest = s-4-2.0-7-1.3-0.175-0m.mags
    # out = s-4-2.0-7-1.3-0.175.mags
    if noMinority == "True":
        restFilename = filename.replace("-z.mags","-0.mags")
        outputFile =  filename.replace("-z.mags",".mags")
    else:
        restFilename = filename.replace("-zm.mags","-0m.mags")
        outputFile =  filename.replace("-zm.mags",".mags")
    dummyNum1 = -70.0000
    dummyNum2 = 99.000
    
    hdr, data = read(filename)
    restHdr, restData = read(restFilename)
    zed, age, mass, lum, sfr, g, r, i, z, Y = lineToNP(data) 
    rfzed, rfage, rfmass, rflum, rfsfr, rfg, rfr, rfi, rfz, rfY = lineToNP(restData) 

    # mags contain distance modulus and a redshift term in fsps v>2.4 so we subtract them from them

    lumdist = np.array([ cd.luminosity_distance(one_zed) for one_zed in zed])
    redshift_term = np.array([ -2.5*np.log10(1.0+one_zed) for one_zed in zed ])
    
    #lumdist = cd.luminosity_distance(zed) ;# in Mpc
    distanceModulus = 5*np.log10(lumdist*1e6/10.)
    ix = zed > 0
    
    g[ix] = g[ix] - distanceModulus[ix] - redshift_term[ix]
    r[ix] = r[ix] - distanceModulus[ix] - redshift_term[ix]
    i[ix] = i[ix] - distanceModulus[ix] - redshift_term[ix]
    z[ix] = z[ix] - distanceModulus[ix] - redshift_term[ix]
    Y[ix] = Y[ix] - distanceModulus[ix] - redshift_term[ix]

    ix = np.nonzero( (z == dummyNum2) | (Y == dummyNum2) )
    zY = z - Y
    zY[ix] = dummyNum2
    ix = np.nonzero( (g == dummyNum2) | (r == dummyNum2) )
    gr = g - r
    gr[ix] = dummyNum2
    ix = np.nonzero( (r == dummyNum2) | (i == dummyNum2) )
    ri = r - i
    ri[ix] = dummyNum2
    ix = np.nonzero( (i == dummyNum2) | (z == dummyNum2) )
    iz = i - z
    iz[ix] = dummyNum2
    ix = np.nonzero( (rfg == dummyNum2) | (rfr == dummyNum2) )
    rgr = rfg - rfr
    rgr[ix] = dummyNum2
    ix = np.nonzero( (rfg == dummyNum2) | (rfi == dummyNum2) )
    rgi = rfg - rfi
    rgi[ix] = dummyNum2
    ix = np.nonzero( (i == dummyNum2) | (rfi == dummyNum2) )
    rfii = rfi-i
    rfii[ix] = dummyNum2
    ix = np.nonzero( (i == dummyNum2) | (rfr == dummyNum2) )
    rfri = rfr-i
    rfri[ix] = dummyNum2
    
    fd = open(outputFile,"w")
    for j in range(0,len(hdr) - 1) : 
        fd.write(hdr[j])
    fd.write("# z  log(age) log(mass) log(lbol) log(SFR)  g-r 	r-i	i-z	z-Y   (g-r)_o (g-i)_o i_o-i r_o-i\n")
    for j in range(0,zed.size) :
        fd.write("{:6.3f} ".format( zed[j]))
        fd.write("{:7.4f} ".format(age[j]))
        fd.write("{:8.4f} ".format(mass[j]))
        fd.write("{:8.4f} ".format(lum[j]))
        fd.write("{:8.4f} ".format(sfr[j]))
        fd.write("{:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}\n".format(
            gr[j],ri[j],iz[j],zY[j],rgr[j], rgi[j], rfii[j], rfri[j]))
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
        #u.append(data[j][5])
        g.append(data[j][5])
        r.append(data[j][6])
        i.append(data[j][7])
        z.append(data[j][8])
        Y.append(data[j][9])
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
