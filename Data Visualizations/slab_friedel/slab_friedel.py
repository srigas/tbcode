import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import blended_transform_factory

mycol = (0.13333, 0.35294, 0.38824)

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

address = 'C:/Users/rigas/Desktop/'

SCRESULTS = np.loadtxt(address+'SCresults.txt', delimiter=',')
LAYERINDEX = SCRESULTS[:,0]
CHARGESZ = SCRESULTS[:,1]
DELTASZ = SCRESULTS[:,3]
czero = CHARGESZ[51]
deltazero = DELTASZ[51]

CHARGES = (CHARGESZ-czero)/czero
DELTAS = (DELTASZ-deltazero)/deltazero

#XTICKS = []
#DELTATICKS = [0.23, 0.25, 0.27]

fig = plt.figure(figsize=(14, 4))
ax1 = fig.add_subplot(121)
ax1.plot(LAYERINDEX, CHARGES, color=mycol)
ax1.set_title('(a)', pad=12)
#ax1.set_aspect(1.5)
#ax1.set_xticks(XTICKS)
ax1.set_ylabel(r'$\delta n/n$', labelpad=12)
ax1.set_xlabel('Layer index', labelpad=12)


transform = blended_transform_factory(fig.transFigure, ax1.transAxes)
axins1 = inset_axes(ax1, width="50%", height="55%",
                    bbox_to_anchor=(0.03, 0.4, 0.5, 0.5),
                    bbox_transform=transform, loc=8, borderpad=0)

scaxins1 = axins1.plot(LAYERINDEX, CHARGES, color=mycol)
axins1.set_yticks([])
###axins1.set_xticks([])
axins1.set_xlim([3, 98])
axins1.set_ylim([-0.003, 0.003])

ax2 = fig.add_subplot(122)
ax2.plot(LAYERINDEX, DELTAS, color=mycol)
ax2.set_title('(b)', pad=12)
#ax2.set_aspect(1.5)
#ax2.set_yticks(DELTATICKS)
ax2.set_ylabel(r'$\delta \Delta/\Delta$', labelpad=12)
ax2.set_xlabel('Layer index', labelpad=12)

plt.tight_layout(w_pad=2.5)

plt.show()