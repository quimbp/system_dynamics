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
  tfile = 'lorenz.dat'

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

fig = plt.figure()

fig.add_subplot(311)
plt.plot(time,x,'-r',linewidth=2)
plt.ylabel('x',fontweight='bold',fontsize=14)
plt.grid(True)
plt.gca().axes.get_xaxis().set_ticklabels([])

fig.add_subplot(312)
plt.plot(time,y,'-b',linewidth=2)
plt.ylabel('y',fontweight='bold',fontsize=14)
plt.grid(True)
plt.gca().axes.get_xaxis().set_ticklabels([])

fig.add_subplot(313)
plt.plot(time,z,'-g',linewidth=2)
plt.xlabel('time',fontweight='bold',fontsize=14)
plt.ylabel('z',fontweight='bold',fontsize=14)
plt.grid(True)

plt.show()

#print(' > Saving file as lorenz.pdf')
#fig.savefig('lorenz.pdf',bbox_inches='tight')

