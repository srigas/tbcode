from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
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

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=BIGGER_SIZE)
plt.rc('axes', labelsize=BIGGER_SIZE)
plt.rc('xtick', labelsize=MEDIUM_SIZE)
plt.rc('ytick', labelsize=MEDIUM_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)

## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

#plt.rcParams['font.sans-serif'] = "CMU Sans Serif"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['mathtext.fontset'] = 'cm'

######################################################################

address = 'C:/Users/rigas/Desktop/2dhist_def_data.txt' # Change this

VALS = np.loadtxt(address, delimiter=',')
XAXIS = VALS[:,0]
YAXIS = VALS[:,1]
WGHTS = VALS[:,2]

# Comment the following line and delete the norm=divnorm from the hist2d to have a
# uniform colormap
divnorm = colors.TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=16)

plt.hist2d(XAXIS, YAXIS, bins=[73,11], weights=WGHTS, norm=divnorm, cmap=mycmap)

plt.xlabel('X')
plt.ylabel('Y')
#plt.title('Change ' r'$\delta\Delta$' ' of the order parameter')
#plt.legend()
clb = plt.colorbar()
clb.ax.set_title(r'$\mathcal{D}$' '(0)')
plt.show()

