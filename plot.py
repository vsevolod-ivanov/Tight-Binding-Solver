import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib.patches import Rectangle
import math as math
import numpy as np
import sys

dim1 = int(sys.argv[1])
dim2 = int(sys.argv[2])
filein  = sys.argv[3]
#filein2 = sys.argv[2]
fileout = sys.argv[4]

a = np.loadtxt(filein)
x = np.linspace(0,2*math.pi,dim2)
#z = a[:,3]*1000

#a2=np.loadtxt(filein2)
#x2 = a2[:,0] - 10.2
#y2 = a2[:,1] - 10.2
#z2 = a2[:,2]

plt.xlim([-0,6.28])
plt.ylim([-5,5])

#plt.xlim([-0.2,20.4])
#plt.ylim([-0.2,20.4])
for i in range(0,2*dim1):
    y = a[:,i]
    plt.scatter(x, y, s=1, cmap=cm.Reds, edgecolors='none')
#plt.scatter(x2, y2, c=z2, cmap=cm.Blues, edgecolors='none',alpha=0.4)
#bwr

#plt.colorbar(ticks=[ ])
plt.savefig(fileout, dpi=300)
#plt.show()

