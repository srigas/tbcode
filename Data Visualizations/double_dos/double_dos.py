import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

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

HOSTSFULL = np.loadtxt(address+'hostsfull.txt', delimiter=',')
IMPSFULL = np.loadtxt(address+'impsfull.txt', delimiter=',')

FULLXTICKS = [-4.0, 0, 8.0]

ENERGIESFULL = IMPSFULL[:,0]
IMPDENFULL = IMPSFULL[:,1] # Change 1 into any other atom number
HOSTDENFULL = HOSTSFULL[:,1] # Change 1 into any other atom number

######

HOSTSGAP = np.loadtxt(address+'hostsgap.txt', delimiter=',')
IMPSGAP = np.loadtxt(address+'impsgap.txt', delimiter=',')

GAPXTICKS = [-0.14, 0, 0.14]

ENERGIESGAP = IMPSGAP[:,0]
IMPDENGAP = IMPSGAP[:,1] # Change 1 into any other atom number
HOSTDENGAP = HOSTSGAP[:,1] # Change 1 into any other atom number

######

fig = plt.figure(figsize=(10, 4))

# FULL DOS
ax1 = fig.add_subplot(121)

sc1 = ax1.plot(ENERGIESFULL, IMPDENFULL, color=mycol)
sc1 = ax1.plot(ENERGIESFULL, HOSTDENFULL, color='black', linestyle='dotted')
sc1 = ax1.axhline(y=0.0, color='black', linestyle='-')
ax1.set_title('(a)')
#ax1.set_aspect(1.5)
ax1.set_xticks(FULLXTICKS)
ax1.set_xlabel('E', labelpad=12)
ax1.set_ylabel(r'$\mathcal{D}(E)$', labelpad=12)

# GAP DOS
ax2 = fig.add_subplot(122)

sc2 = ax2.plot(ENERGIESGAP, IMPDENGAP, color=mycol)
sc2 = ax2.plot(ENERGIESGAP, HOSTDENGAP, color='black', linestyle='dotted')
sc2 = ax2.axhline(y=0.0, color='black', linestyle='-')
ax2.set_title('(b)')
#ax2.set_aspect(1.5)
ax2.set_xticks(GAPXTICKS)
ax2.set_xlabel('E', labelpad=12)
#ax2.set_ylabel(r'$\mathcal{D}(E)$', labelpad=12)

plt.tight_layout(w_pad=3.5, h_pad=2.0)

plt.show()