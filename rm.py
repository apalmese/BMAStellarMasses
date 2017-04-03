# z,id,mass,lam,gro = np.genfromtxt("rmBCG-stellar-masses-lam.txt",unpack=True,usecols=(2,3,18,20,4))
import numpy as np
def do(z,mass,lamda,gro, z1=0.2, z2=0.3) :
    ixz = np.nonzero( (gro >0) & (z >= z1) & (z < z2))
    zed = z[ixz]
    mas = mass[ixz]
    lam = lamda[ixz]

    newmass, newlam, newsig = [],[],[]
    min = int(lam.min() )
    max = int(lam.max())
    max1 = 51
    min2 = max1
    #for i in range(min,max1) :
    #    ix = np.nonzero(lam == i)
    #    if ix[0].size > 2 :
    #        newsig.append(mas[ix].std())
    #        newmass.append( np.median(mas[ix]) )
    #        newlam.append(np.log10(i))
    for i in [1.3,1.4,1.5,1.6,1.7, 1.8, 1.9, 2.0] :
        ix = np.nonzero((np.log10(lam)  >= i)&(np.log10(lam) <=i+0.1))
        if ix[0].size > 2 :
            newsig.append(mas[ix].std()/np.sqrt(mas[ix].size))
            newmass.append( np.median(mas[ix]) )
            newlam.append(i+0.05)
    #for i in [2.0,2.2] :
    #    ix = np.nonzero((np.log10(lam)  >= i)&(np.log10(lam) <=i+0.2))
    #    if ix[0].size > 2 :
    #        newsig.append(mas[ix].std())
    #        newmass.append( np.median(mas[ix]) )
    #        if i != 2.2 :
    #            newlam.append(i+0.1)
    #        else :
    #            newlam.append(i+0.05)
    newmass = np.array(newmass)
    newlam = np.array(newlam)
    newsig = np.array(newsig)
    return newlam, newmass, newsig, ixz

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def func(x,a0,a1): return a0 + a1*x

def plot(z,mass,lamda,newlam,newm,news,ix,z1=0.2,z2=0.3) :

    size = lamda[ix].size 
    deviate = np.random.normal(scale=0.01,size=size)
    loglam = np.log10(lamda[ix])+deviate
    plt.hexbin(loglam,mass[ix])
    #poly,cov=np.polyfit(newlam,newm,1,cov=True);
    poly,cov = curve_fit(func,newlam,newm,sigma=news)
    print poly
    print cov
    print "{:6.3f} {:6.3f} +- {:6.3f}  {:6.3f}".format(
        poly[0],poly[1],np.sqrt(cov[0][0]),np.sqrt(cov[1][1]))
    
    x=np.arange(1.3,2.3,0.01)
    #plt.plot(x,poly[1]+x*poly[0], c="w",lw=2)
    plt.plot(x,poly[0]+x*poly[1], c="w",lw=2)
    plt.errorbar(newlam,newm,yerr=news,fmt="o", ecolor="k",
        mfc='red', mec='k', ms=7, mew=2,lw=2)
    plt.xlim(np.log10(20),np.log10(200-5));plt.ylim(10.9,12.4)
    plt.xlabel("log(Lambda)")
    plt.ylabel("log(stellar mass)")
    plt.title("redmapper {:3.1f} <= z <= {:3.1f}".format(z1,z2))

# we pivot around lambda = 20
def norm(mass, lamda) :
    mass = mass - 0.45*np.log10(lamda/20.)
    return mass
####
def donorm(mass, lamda, z, gro) :
    z1=0.1;z2=0.2; donorma(mass, lamda,z,gro,z1,z2)
    z1=0.2;z2=0.3; donorma(mass, lamda,z,gro,z1,z2)
    z1=0.3;z2=0.4; donorma(mass, lamda,z,gro,z1,z2)
    z1=0.4;z2=0.5; donorma(mass, lamda,z,gro,z1,z2)
    z1=0.5;z2=0.6; donorma(mass, lamda,z,gro,z1,z2)
def donorma(mass, lamda, z, gro,z1,z2) :
    ix = np.nonzero((gro>0)&(z>=z1)&(z<z2))
    m= mass[ix]
    median = np.median(m)
    std =  m.std()
    sdm = std/np.sqrt(m.size)
    n = m.size
    print z1,z2,median, sdm, std, n
def donormb(mass, lamda, z, gro,z1,z2) :
    ix = np.nonzero((gro>0)&(z>=z1)&(z<z2))
    m= norm(mass[ix],lamda[ix])
    median = np.median(m)
    std =  m.std()
    sdm = std/np.sqrt(m.size)
    n = m.size
    print z1,z2,median, sdm, std, n
###

def donorm2(mass, lamda, z, gro) :
    z1=0.1;z2=0.2; donorm2a(mass, lamda,z,gro,z1,z2)
    z1=0.2;z2=0.3; donorm2a(mass, lamda,z,gro,z1,z2)
    z1=0.3;z2=0.4; donorm2a(mass, lamda,z,gro,z1,z2)
    z1=0.4;z2=0.5; donorm2a(mass, lamda,z,gro,z1,z2)
    z1=0.5;z2=0.6; donorm2a(mass, lamda,z,gro,z1,z2)
def donorm2a (mass, lamda, z, gro, z1, z2) :
    ix = np.nonzero((gro>0)&(z>=z1)&(z<z2))
    m= norm2(mass[ix],lamda[ix],z1)
    median = np.median(m)
    std =  m.std()
    sdm = std/np.sqrt(m.size)
    n = m.size
    print z1,z2,median, sdm, std, n
def norm2(mass, lamda, z) :
    #if z == 0.1 : slope = 0.46
    #elif z == 0.2 : slope = 0.45
    #elif z == 0.3 : slope = 0.39
    #elif z == 0.4 : slope = 0.41
    #elif z == 0.5 : slope = 0.43
    if z == 0.1 : slope = 0.39
    elif z == 0.2 : slope = 0.42
    elif z == 0.3 : slope = 0.36
    elif z == 0.4 : slope = 0.34
    elif z == 0.5 : slope = 0.53
    mass = mass - slope*np.log10(lamda/20.)
    return mass
    
