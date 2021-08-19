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

address = 'C:/Users/rigas/Desktop/'

THETA4 = np.loadtxt(address+'THETA4.txt', delimiter=',')
THETA8 = np.loadtxt(address+'THETA8.txt', delimiter=',')
THETA12 = np.loadtxt(address+'THETA12.txt', delimiter=',')
THETA18 = np.loadtxt(address+'THETA18.txt', delimiter=',')

BMAGS = THETA4[:,0] # FOR BMAGS
ENERGIES = THETA4[:,1] # FOR ENERGIES

EBDOS4 = THETA4[:,2]
EBDOS8 = THETA8[:,2]
EBDOS12 = THETA12[:,2]
EBDOS18 = THETA18[:,2]

#-----------------------------------------------#

rcParams['axes.titlepad'] = 12

fig = plt.figure(figsize=(13, 10))

# THETA = PI/4
ax1 = fig.add_subplot(221)
sc1 = ax1.scatter(BMAGS, ENERGIES, s=0.5, c=EBDOS4, cmap=mycmap)
ax1.set_title('(a)')
#ax1.set_aspect(1.5)
ax1.set_yticks(EBDOSYTICKS)
ax1.set_xlabel('B', labelpad=14)
ax1.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc1, ax=ax1, ticks=[5, 10, 15, 20])
cbar.ax.set_title(r'$\mathcal{D}(E)$')


# THETA = PI/8
ax2 = fig.add_subplot(222)
sc2 = ax2.scatter(BMAGS, ENERGIES, s=0.5, c=EBDOS8, cmap=mycmap)
ax2.set_title('(b)')
#ax2.set_aspect(1.5)
ax2.set_yticks(EBDOSYTICKS)
ax2.set_xlabel('B', labelpad=14)
#ax2.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc2, ax=ax2, ticks=[5, 10, 15, 20, 25])
cbar.ax.set_title(r'$\mathcal{D}(E)$')


# THETA = PI/12
ax3 = fig.add_subplot(223)
sc3 = ax3.scatter(BMAGS, ENERGIES, s=0.5, c=EBDOS12, cmap=mycmap)
ax3.set_title('(c)')
#ax3.set_aspect(1.5)
ax3.set_yticks(EBDOSYTICKS)
ax3.set_xlabel('B', labelpad=14)
ax3.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc3, ax=ax3, ticks=[5, 10, 15, 20, 25])
cbar.ax.set_title(r'$\mathcal{D}(E)$')


# THETA = PI/18
ax4 = fig.add_subplot(224)
sc4 = ax4.scatter(BMAGS, ENERGIES, s=0.5, c=EBDOS18, cmap=mycmap)
ax4.set_title('(d)')
#ax4.set_aspect(1.5)
ax4.set_yticks(EBDOSYTICKS)
ax4.set_xlabel('B', labelpad=14)
#ax4.set_ylabel('E', labelpad=12)

cbar = fig.colorbar(sc4, ax=ax4, ticks=[5, 10, 15, 20, 25])
cbar.ax.set_title(r'$\mathcal{D}(E)$')


plt.tight_layout(w_pad=3.5, h_pad=2.0)

#plt.savefig('diffhelixangles-hd.pdf')
plt.savefig('diffhelixangles.png')
#plt.show()