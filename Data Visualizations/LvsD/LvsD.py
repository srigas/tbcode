import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

mycol = (0.13333, 0.35294, 0.38824)
#lightcol = (0.55294, 0.67059, 0.67059) # 141 171 171

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
LARGEST_SIZE = 20
YTICK_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGEST_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

######################################################################

rcParams['axes.titlepad'] = 12

#YTICKS = []
address = 'C:/Users/rigas/Desktop/'

SC = np.loadtxt(address+'LvsDsc.txt', delimiter=',')
LAMBDAS = SC[:,0]
DELTASSC = SC[:,1]

FCC = np.loadtxt(address+'LvsDfcc.txt', delimiter=',')
DELTASFCC = FCC[:,1]

BCC = np.loadtxt(address+'LvsDbcc.txt', delimiter=',')
DELTASBCC = BCC[:,1]

fig = plt.figure(figsize=(12, 4))

#scatterSC = ax.scatter(x=LAMBDAS, y=DELTASSC, color=mycol, s=3.0, marker='o')
#scatterFCC = ax.scatter(x=LAMBDAS, y=DELTASFCC, color=mycol, s=3.0, marker='+')
#scatterBCC = ax.scatter(x=LAMBDAS, y=DELTASBCC, color=mycol, s=3.0, marker='s')

ax1 = fig.add_subplot(121)
plotSC = ax1.plot(LAMBDAS, DELTASSC, color=mycol, label='SC')
plotBCC = ax1.plot(LAMBDAS, DELTASBCC, color=mycol, linestyle='dashed', label='BCC')
plotFCC = ax1.plot(LAMBDAS, DELTASFCC, color=mycol, linestyle='dotted', label='FCC')

ax1.set_title('(a)')
#ax.set_aspect(1.5)
#ax.set_yticks(YTICKS)
ax1.set_ylabel(r'$\Delta$', labelpad=12)
ax1.set_xlabel(r'$\Lambda$', labelpad=12)
#ax.set_yticklabels(['120', '90', '72', '60', '45', '40', '36', '24'])

leg = ax1.legend()

ax2 = fig.add_subplot(122)
plotSC = ax2.plot(LAMBDAS, DELTASSC, color=mycol, label='SC')
plotBCC = ax2.plot(LAMBDAS, DELTASBCC, color=mycol, linestyle='dashed', label='BCC')
plotFCC = ax2.plot(LAMBDAS, DELTASFCC, color=mycol, linestyle='dotted', label='FCC')

ax2.set_title('(b)')
#ax.set_aspect(1.5)
#ax.set_yticks(YTICKS)
#ax2.set_ylabel(r'$\Delta$', labelpad=12)
ax2.set_xlabel(r'$\Lambda$', labelpad=12)
ax2.set_yscale('log')
ax2.set_xscale('log')
#ax.set_yticklabels(['120', '90', '72', '60', '45', '40', '36', '24'])

leg = ax2.legend()

plt.tight_layout(w_pad=2.5)
plt.show()