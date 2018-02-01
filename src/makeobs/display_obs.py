import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

try:
  import tkinter as tk
  from tkinter import ttk
  from tkinter import messagebox
  from tkinter import filedialog
except:
  import Tkinter as tk
  import ttk
  import tkMessageBox as messagebox
  import tkFileDialog as filedialog


font = {'family' : 'DejaVu Sans',
        'weight'   : 'bold',
        'size'     : 14}

fig = {'figsize': (12,8),
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



# ==============
class WinPlot():
# ==============

  def __init__(self,master,tfile,ofile):

    self.master = master
    self.tfile = tk.StringVar()
    self.ofile = tk.StringVar()

    self.tfile.set(tfile)
    self.ofile.set(ofile)

    F0 = ttk.Frame(master,padding=5)
    ttk.Label(F0,text='Truth trajectory filename').grid(row=0,    \
                                                        column=0, \
                                                        padx=3,   \
                                                        pady=5,   \
                                                        sticky='e')
    ttk.Entry(F0,textvariable=self.tfile,width=40).grid(row=0,        \
                                                        column=1,     \
                                                        columnspan=4, \
                                                        padx=3,       \
                                                        pady=5)
    ttk.Button(F0,text='Select',command=lambda: self.select_in(self.tfile)).   \
                                                   grid(row=0,              \
                                                   column=5,                \
                                                   padx=3,                  \
                                                   pady=5)
  

    ttk.Label(F0,text='Observation trajectory filename').grid(row=1, \
                                                        column=0,    \
                                                        padx=3,      \
                                                        pady=5,      \
                                                        sticky='e')
    ttk.Entry(F0,textvariable=self.ofile,width=40).grid(row=1,        \
                                                        column=1,     \
                                                        columnspan=4, \
                                                        padx=3,       \
                                                        pady=5)
    ttk.Button(F0,text='Select',command=lambda:self.select_in(self.ofile)). \
                                                   grid(row=1,              \
                                                   column=5,                \
                                                   padx=3,                  \
                                                   pady=5)

    ttk.Button(F0,text='Draw',command=self.make_plot).            \
                                                   grid(row=2,    \
                                                   column=4,      \
                                                   padx=3,        \
                                                   pady=5)
    ttk.Button(F0,text='Close',command=quit).           \
                                         grid(row=2,    \
                                         column=5,      \
                                         padx=3,        \
                                         pady=5)
    F0.grid()


  def make_plot(self):
  # ------------------

    tt = DataStruct()
    tt.Read(self.tfile.get())

    oo = DataStruct()
    oo.Read(self.ofile.get())

    self.fig = plt.figure(dpi=100,frameon=True)

    ax = self.fig.add_subplot(311)
    ax.plot(tt.time,tt.x,'--r',linewidth=0.5)
    ax.plot(oo.time,oo.x,'o',color='red')
    ax.set_ylabel('x',fontweight='bold',fontsize=14)
    ax.grid(True)
    ax.axes.get_xaxis().set_ticklabels([])

    ay = self.fig.add_subplot(312)
    ay.plot(tt.time,tt.y,'--b',linewidth=0.5)
    ay.plot(oo.time,oo.y,'o',color='blue')
    ay.set_ylabel('y',fontweight='bold',fontsize=14)
    ay.grid(True)
    ay.axes.get_xaxis().set_ticklabels([])

    az = self.fig.add_subplot(313)
    az.plot(tt.time,tt.z,'--g',linewidth=0.5)
    az.plot(oo.time,oo.z,'o',color='green')
    az.set_xlabel('time',fontweight='bold',fontsize=14)
    az.set_ylabel('z',fontweight='bold',fontsize=14)
    az.grid(True)

    plt.show()


  def select_in(self,FILENAME):
  # ---------------------------
    nn = filedialog.askopenfile()
    if nn is not None:
      FILENAME.set(nn.name)
    else:
      FILENAME.set('')


# =========
def main():
# =========

  print(' ')
  print(' ================')
  print(' > obs_display.py')

  try:
    tfile = sys.argv[1]
    try:
      ofile = sys.argv[2]
    except:
      ofile = ''
      
  except:
    tfile = ''
    ofile = ''

  print(tfile)
  print(ofile)

  root = tk.Tk()
  root.title('OBS_DISPLAY')
  root.grid_rowconfigure(0,weight=0)
  root.resizable(width=False,height=True)
  root.protocol('WM_DELETE_WINDOW',quit)

  WinPlot(root,tfile,ofile)
  root.mainloop()


if __name__ == '__main__':
  main()
