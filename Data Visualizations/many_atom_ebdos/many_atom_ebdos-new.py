from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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

TICKS = [-0.14,0,0.14]

address = 'C:/Users/rigas/Desktop/'

EBDOS31 = np.loadtxt(address+'EBDOS31.txt', delimiter=',')
BMAGS = EBDOS31[:,0]
ENERGIES = EBDOS31[:,1]
DOS31 = EBDOS31[:,2]

EBDOS61 = np.loadtxt(address+'EBDOS61.txt', delimiter=',')
DOS61 = EBDOS61[:,2]

EBDOS91 = np.loadtxt(address+'EBDOS91.txt', delimiter=',')
DOS91 = EBDOS91[:,2]

miniebdos31 = np.loadtxt(address+'miniebdos31.txt', delimiter=',')
minibmags = miniebdos31[:,0]
minienergies = miniebdos31[:,1]
minidos31 = miniebdos31[:,2]

miniebdos61 = np.loadtxt(address+'miniebdos61.txt', delimiter=',')
minidos61 = miniebdos61[:,2]

miniebdos91 = np.loadtxt(address+'miniebdos91.txt', delimiter=',')
minidos91 = miniebdos91[:,2]


rcParams['axes.titlepad'] = 12

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(311)
sc1 = ax1.scatter(BMAGS, ENERGIES, s=0.5, c=DOS31, cmap=mycmap)
ax1.set_title('(a)') #31 impurities
#ax1.set_aspect(1.5)
ax1.set_yticks(TICKS)
ax1.set_xlabel('B', labelpad=14)
ax1.set_ylabel('E', labelpad=12)

axins1 = inset_axes(ax1, width=2.5, height=1.5)
scaxins1 = axins1.scatter(minibmags, minienergies, s=0.5, c=minidos31, cmap=mycmap)
axins1.set_yticks([0])
axins1.set_xticks([1.6,2.1])

cbar = fig.colorbar(sc1, ax=ax1, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax2 = fig.add_subplot(312)
sc2 = ax2.scatter(BMAGS, ENERGIES, s=0.5, c=DOS61, cmap=mycmap)
#ax2.set_aspect(5.5)
ax2.set_title('(b)') #61 impurities
ax2.set_yticks(TICKS)
ax2.set_xlabel('B', labelpad=14)
ax2.set_ylabel('E', labelpad=12)
#ax2.set_ylim([-1.1, 1.1])

axins2 = inset_axes(ax2, width=2.5, height=1.5)
scaxins2 = axins2.scatter(minibmags, minienergies, s=0.5, c=minidos61, cmap=mycmap)
axins2.set_yticks([0])
axins2.set_xticks([1.6,2.1])

cbar = fig.colorbar(sc2, ax=ax2, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax3 = fig.add_subplot(313)
sc3 = ax3.scatter(BMAGS, ENERGIES, s=0.5, c=DOS91, cmap=mycmap)
#ax3.set_aspect(5.5)
ax3.set_title('(c)') #91 impurities
ax3.set_yticks(TICKS)
ax3.set_xlabel('B', labelpad=14)
ax3.set_ylabel('E', labelpad=12)
#ax3.set_ylim([-1.1, 1.1])

axins3 = inset_axes(ax3, width=2.5, height=1.5)
scaxins3 = axins3.scatter(minibmags, minienergies, s=0.5, c=minidos91, cmap=mycmap)
axins3.set_yticks([0])
axins3.set_xticks([1.6,2.1])

cbar = fig.colorbar(sc3, ax=ax3, ticks=[2.5, 7.5, 12.5, 17.5])
cbar.ax.set_title(r'$\mathcal{D}(E)$')

plt.tight_layout(h_pad=3.0)

#plt.savefig('manyatomebdos.pdf')
plt.savefig('manyatomebdos.png')
#plt.show()