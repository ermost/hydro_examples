import numpy as np
import matplotlib.pyplot as plt
from numpy import log10
import sys

input_string =sys.argv[1]

name_string =input_string
if name_string.endswith('.asc'):
        name_string = name_string[:-4]

x, y, u = np.loadtxt(input_string, usecols=(0,1,2), unpack=True, delimiter='\t')

jet = plt.cm.get_cmap('jet') 

fig, ax1 = plt.subplots(1, 1)

ax1.set_title(name_string)

# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
CS=ax1.tripcolor(x,y,u, vmin=0, vmax=2, cmap=jet,shading='gouraud')

ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_autoscaley_on(True)
ax1.set_autoscalex_on(True)

cb= plt.colorbar(CS, extend='both')
#cb.set_label(label=r'$u$')

xmin = np.amin(x)
ymin = np.amin(y)
xmax = np.amax(x)
ymax = np.amax(y)

ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])

#fig.show()
plt.savefig(name_string + '.png', dpi=150.)


