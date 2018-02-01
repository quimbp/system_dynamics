import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

try:
  #import tkinter as tk
  #from tkinter import ttk
  from tkinter import messagebox
  from tkinter import filedialog
except:
  #import Tkinter as tk
  #import ttk
  import tkMessageBox as messagebox
  import tkFileDialog as filedialog


# Plot and figure general options
#
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 14}

fig = {'figsize': (12,6),
       'dpi'    : 100}

savefig = {'dpi'       : 100,
           'directory' : './plots',
           'format'    : 'png'}

matplotlib.rc('font',**font)
matplotlib.rc('figure',**fig)
matplotlib.rc('savefig',**savefig)



# ================
def empty(string):
# ================
  '''Logical function that checks if a string is empty:
     empty('     ') = True'''

  __version__ = "1.0"
  __author__  = "Quim Ballabrerera"
  __date__    = "June 2017"

  if bool(string.strip()):
    return False
  else:
    return True

# =================
class DataStruct():
# =================

  def __init__(self):
  # -----------------

    self.time = []
    self.x    = []
    self.y    = []
    self.z    = []

  def Read(self,filename):
  # ----------------------
    print(' > Input file: ',filename)

    self.time = []
    self.x = []
    self.y = []
    self.z = []

    F = open(filename,'r')
    for line in F:
      tmp = line.split()
      self.time.append(float(tmp[0]))
      self.x.append(float(tmp[1]))
      self.y.append(float(tmp[2]))
      self.z.append(float(tmp[3]))
    F.close()

  def Plot(self):
  # ------------------

    self.fig = plt.figure()

    ax = self.fig.add_subplot(311)
    ax.plot(self.time,self.x,'-',color='red',linewidth=1,label='Truth')
    ax.set_ylabel('x',fontweight='bold',fontsize=14)
    ax.set_xlim(self.time[0],self.time[-1])
    #ax.set_ylim([-30,30])
    ax.get_xaxis().set_ticklabels([])
    ax.grid(True)


    ay = self.fig.add_subplot(312)
    ay.plot(self.time,self.y,'-',color='blue',linewidth=1,label='Truth')
    ay.set_ylabel('y',fontweight='bold',fontsize=14)
    ay.set_xlim(self.time[0],self.time[-1])
    #ay.set_ylim([-30,30])
    ay.get_xaxis().set_ticklabels([])
    ay.grid(True)


    az = self.fig.add_subplot(313)
    az.plot(self.time,self.z,'-',color='green',linewidth=1,label='Truth')
    az.set_xlabel('time',fontweight='bold',fontsize=14)
    az.set_ylabel('z',fontweight='bold',fontsize=14)
    az.set_xlim(self.time[0],self.time[-1])
    #az.set_ylim([0,50])
    az.grid(True)

    plt.show()


