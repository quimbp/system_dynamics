#! /usr/bin/env python3

import sys
#import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'DejaVu Sans',
        'weight'   : 'bold',
        'size'     : 14}

fig = {'figsize': (12,8),
       'dpi'    : 100}

savefig = {'dpi'       : 100,
           'directory' : './plots',
           'format'    : 'pdf'}

matplotlib.rc('font',**font)
matplotlib.rc('figure',**fig)
matplotlib.rc('savefig',**savefig)

print(' ')
print(' =============')
print(' > lorenz.py')

try:
  tfile = sys.argv[1]
except:
  tfile = 'lorenz_simulation.dat'

try:
  title = sys.argv[1]
except:
  title = 'Lorenz model'

print(' > Input file: ',tfile)

F = open(tfile,'r')

time = []
x = []
y = []
z = []
for line in F:
  tmp = line.split()
  time.append(float(tmp[0]))
  x.append(float(tmp[1]))
  y.append(float(tmp[2]))
  z.append(float(tmp[3]))
F.close()

fig = plt.figure('Lorenz')

ax = fig.add_subplot(311)
ax.plot(time,x,'-r',linewidth=2)
ax.set_title(title,fontsize=16)
ax.set_ylabel('x',fontweight='bold',fontsize=14)
ax.grid(True)
ax.get_xaxis().set_ticklabels([])

ay = fig.add_subplot(312)
ay.plot(time,y,'-b',linewidth=2)
ay.set_ylabel('y',fontweight='bold',fontsize=14)
ay.grid(True)
ay.get_xaxis().set_ticklabels([])

az = fig.add_subplot(313)
az.plot(time,z,'-g',linewidth=2)
az.set_xlabel('time',fontweight='bold',fontsize=14)
az.set_ylabel('z',fontweight='bold',fontsize=14)
az.grid(True)

plt.show()

#print(' > Saving file as lorenz.pdf')
#fig.savefig('lorenz.pdf',bbox_inches='tight')

