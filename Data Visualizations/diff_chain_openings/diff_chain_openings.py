from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

def CustomCmap(from_rgb,to_rgb):

    # from color r,g,b
    r1,g1,b1 = from_rgb

    # to color r,g,b
    r2,g2,b2 = to_rgb

    cdict = {'red': ((0, r1, r1),
                   (1, r2, r2)),
           'green': ((0, g1, g1),
                    (1, g2, g2)),
           'blue': ((0, b1, b1),
                   (1, b2, b2))}

    cmap = LinearSegmentedColormap('custom_cmap', cdict)
    return cmap


mycmap = CustomCmap([1.00, 1.00, 1.00], [0.13333, 0.35294, 0.38824]) # from white to teal

mycol = (0.13333, 0.35294, 0.38824)

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
LARGEST_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGEST_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

######################################################################
#------------ FILES AND SETUP -----------------#

EBDOSYTICKS = [-0.14,0,0.14]
DOSXTICKS = [-0.14, 0, 0.14]

address = 'C:/Users/rigas/Desktop/' # Change this

EBDOS30 = np.loadtxt(address + 'EBDOS30.txt', delimiter=',')
EBDOS31 = np.loadtxt(address + 'EBDOS31.txt', delimiter=',')
EBDOS32 = np.loadtxt(address + 'EBDOS32.txt', delimiter=',')
HOSTDEN = np.loadtxt(address + 'hostdensities.txt', delimiter=',')
IMPDEN30 = np.loadtxt(address + 'impdensities30.txt', delimiter=',')
IMPDEN31 = np.loadtxt(address + 'impdensities31.txt', delimiter=',')
IMPDEN32 = np.loadtxt(address + 'impdensities32.txt', delimiter=',')

BMAGS = EBDOS30[:,0] # FOR EBDOS
ENERGIES = EBDOS30[:,1] # FOR EBDOS
DOSENERGIES = HOSTDEN[:,0] # FOR DOS
HOSTDOS = HOSTDEN[:,1] # FOR DOS - Same for all atoms

DOS30 = EBDOS30[:,2] # FOR EBDOS
DOS31 = EBDOS31[:,2] # FOR EBDOS
DOS32 = EBDOS32[:,2] # FOR EBDOS
IMP30 = IMPDEN30[:,1] # FOR DOS
IMP31 = IMPDEN31[:,1] # FOR DOS
IMP32 = IMPDEN32[:,1] # FOR DOS

#-----------------------------------------------#

rcParams['axes.titlepad'] = 12

fig = plt.figure(figsize=(13, 14))

# 30 ATOMS
ax1 = fig.add_subplot(321)
sc1 = ax1.scatter(BMAGS, ENERGIES, s=0.5, c=DOS30, cmap=mycmap)
ax1.set_title('(a)')
#ax1.set_aspect(1.5)
ax1.set_yticks(EBDOSYTICKS)
ax1.set_xlabel('B', labelpad=14)
ax1.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc1, ax=ax1, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax2 = fig.add_subplot(322)

ax2.plot(DOSENERGIES, IMP30, color=mycol)
ax2.plot(DOSENERGIES, HOSTDOS, color='black', linestyle='dashed')
ax2.set_title('(b)')
#ax2.set_aspect(1.5)
ax2.set_xticks(DOSXTICKS)
ax2.set_xlabel('E', labelpad=12)
ax2.set_ylabel(r'$\mathcal{D}(E)$', labelpad=12)
#ax2.set_ylim([-0.5, 14.0])


# 31 ATOMS
ax3 = fig.add_subplot(323)
sc3 = ax3.scatter(BMAGS, ENERGIES, s=0.5, c=DOS31, cmap=mycmap)
ax3.set_title('(c)')
#ax3.set_aspect(1.5)
ax3.set_yticks(EBDOSYTICKS)
ax3.set_xlabel('B', labelpad=14)
ax3.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc3, ax=ax3, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax4 = fig.add_subplot(324)

ax4.plot(DOSENERGIES, IMP31, color=mycol)
ax4.plot(DOSENERGIES, HOSTDOS, color='black', linestyle='dashed')
ax4.set_title('(d)')
#ax4.set_aspect(1.5)
ax4.set_xticks(DOSXTICKS)
ax4.set_xlabel('E', labelpad=12)
ax4.set_ylabel(r'$\mathcal{D}(E)$', labelpad=12)
#ax4.set_ylim([-0.5, 14.0])


# 32 ATOMS
ax5 = fig.add_subplot(325)
sc5 = ax5.scatter(BMAGS, ENERGIES, s=0.5, c=DOS32, cmap=mycmap)
ax5.set_title('(e)')
#ax5.set_aspect(1.5)
ax5.set_yticks(EBDOSYTICKS)
ax5.set_xlabel('B', labelpad=14)
ax5.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc5, ax=ax5, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax6 = fig.add_subplot(326)

ax6.plot(DOSENERGIES, IMP32, color=mycol)
ax6.plot(DOSENERGIES, HOSTDOS, color='black', linestyle='dashed')
ax6.set_title('(f)')
#ax6.set_aspect(1.5)
ax6.set_xticks(DOSXTICKS)
ax6.set_xlabel('E', labelpad=12)
ax6.set_ylabel(r'$\mathcal{D}(E)$', labelpad=12)
#ax6.set_ylim([-0.5, 14.0])

plt.tight_layout(w_pad=3.5, h_pad=2.0)

#plt.savefig('diffchainopenings.pdf')
plt.savefig('diffchainopenings.png')
#plt.show()