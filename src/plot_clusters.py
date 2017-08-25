import matplotlib.pyplot as plt
import numpy as np

ind=np.genfromtxt(infile,dtype=str)

grfl=ind[:,51]
rifl=ind[:,52]
izfl=ind[:,53]


rifl_cl=np.full(halos.shape,'true',dtype=str)
i=0
for halo in halos:
    ix = np.nonzero(hostid == halo)
    rifl_cl[i]=rifl[ix[0][0]]
    i=i+1

ri = d[rifl_cl=='t',9]
l = d[rifl_cl=='t',7]

plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.title('r-i, flag_sep_ri==1')
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda_\star$')
plt.scatter(l,ri,s=5,edgecolor='none')
plt.savefig('ri_flag.png')

izfl_cl=np.full(halos.shape,'true',dtype=str)
i=0
for halo in halos:
    ix = np.nonzero(hostid == halo)
    izfl_cl[i]=izfl[ix[0][0]]
    i=i+1

iz = d[izfl_cl=='t',10]
l = d[izfl_cl=='t',7]

plt.clf()
plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.title('i-z, flag_sep_iz==1')
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda_\star$')
plt.scatter(l,iz,s=5,edgecolor='none')
plt.savefig('iz_flag.png')

gr_iz_sel = d[izfl_cl=='t',8]
l = d[izfl_cl=='t',7]

plt.clf()
plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.title('g-r, flag_sep_iz==1')
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda_\star$')
plt.scatter(l,gr_iz_sel,s=5,edgecolor='none')
plt.savefig('gr_iz_flag.png')


in_float=np.genfromtxt(infile)
z=in_float[:,44]
z_cl=np.zeros(halos.shape)

i=0
for halo in halos:
    ix = np.nonzero(hostid == halo)
    z_cl[i]=z[ix[0][0]]
    i=i+1

z_fl= z_cl[izfl_cl=='t']

plt.clf()
plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.title('i-z, flag_sep_iz==1')
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda_\star$')
plt.scatter(l,iz,s=5,edgecolor='none',c=z_fl)
plt.colorbar(label='$z_{ph}$')
plt.savefig('redshift_iz_flag.png')
