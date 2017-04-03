import numpy as np


in=np.genfromtxt('DeLucia_smass.in')
cin=np.genfromtxt('DeLucia_smass.in')
cin.shape
cout=np.genfromtxt('delucia_masses.txt',delimiter=",")
cout.shape
rabs=cout[:,11]
haloids=cin[:,1]
halos=np.unique(haloids)
halos=np.sort(halos)
outfile='Delucia_halos.txt'

fd = open(outfile,"w")

for halo in halos:
    idx=np.logical_and(haloids==halo,rabs<-19.)
    linearmass=10**cout[idx,15]
    sum_mass = linearmass.sum()
    mass_errors = np.log(10.)*linearmass*cout[idx,16]
    mass_std = np.sqrt((mass_errors**2).sum())
    sum_mass_std = mass_std/sum_mass/np.log(10.)
    sum_mass = np.log10(sum_mass)
    m200=np.log10(10.0**10*cin[idx[0],3])
    linearmass_true=10**10*cin[idx,6]
    sum_mass_true=linearmass_true.sum()
    sum_mass_true = np.log10(sum_mass_true)
    z=cin[idx[0],4]
    fd.write("{:10d} {:6.3f} {:10.6f} {:10.5f}  {:10.6f} {:10.5f}\n".format(halo,z, sum_mass, sum_mass_std, sum_mass_true, m200))

fd.close()


