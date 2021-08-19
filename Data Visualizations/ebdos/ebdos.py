from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib import rcParams
from matplotlib import colors

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

divnorm = colors.TwoSlopeNorm(vmin=0.0, vcenter=3.0, vmax=55.523)

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
LARGEST_SIZE = 20
YTICK_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

rcParams['axes.titlepad'] = 12

######################################################################

address = 'C:/Users/rigas/Desktop/ebdos_def_data.txt' # Change this
EBDOS = np.loadtxt(address, delimiter=',')

BMAGS = EBDOS[:,0]
ENERGIES = EBDOS[:,1]
DOS = EBDOS[:,2]

YTICKS = [-0.14,0.00,0.14]

fig = plt.figure(figsize=(14, 5))

ax1 = fig.add_subplot(121)

sc1 = ax1.scatter(BMAGS, ENERGIES, s=0.5, c=DOS, cmap=mycmap)
ax1.set_yticks(YTICKS)
ax1.set_xlabel('B', labelpad=14)
ax1.set_ylabel('E', labelpad=12)
ax1.set_title('(a)')

cbar = fig.colorbar(sc1, ax=ax1)
cbar.ax.set_title(r'$\mathcal{D}(E)$')

ax2 = fig.add_subplot(122)

sc2 = ax2.scatter(BMAGS, ENERGIES, s=0.5, c=DOS, norm=divnorm, cmap=mycmap)
ax2.set_yticks(YTICKS)
ax2.set_xlabel('B', labelpad=14)
ax2.set_title('(b)')

cbar = fig.colorbar(sc2, ax=ax2)
cbar.ax.set_title(r'$\mathcal{D}(E)$')

plt.tight_layout(w_pad=3.5, h_pad=2.0)
plt.show()