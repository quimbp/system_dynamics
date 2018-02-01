#! /usr/bin/env python3
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



# ==============
class WinPlot():
# ==============

  def __init__(self,master,tfile,experiment):

    methods_list = ['Interpolation', 'Nudging','EnKF Background',
                    'EnKF Analysis', '4DVAR']
    self.mdict = {'Interpolation': 'interp_',   \
                  'Nudging': 'nud_',            \
                  'EnKF Background': 'enkf_b_', \
                  'EnKF Analysis': 'enkf_a_',   \
                  '4DVAR': 'fdvar_' }

    self.master = master

    self.method = tk.StringVar()
    self.method.set('Interpolation')

    self.tfile = tk.StringVar()

    self.exper = tk.StringVar()
    if empty(experiment):
      self.exper.set('noise_02_sampl_010')  
    else:
      self.exper.set(experiment)


    self.ofile = tk.StringVar()
    self.ifile = tk.StringVar()
    self.ofile.set('obs_'+self.exper.get()+'.dat')
    self.ifile.set(self.mdict[self.method.get()]+self.exper.get()+'.dat')

    self.tfile.set(tfile)
    #self.ofile.set(ofile)

    F0 = ttk.Frame(master,padding=5)
    ttk.Label(F0,text='Methodology').grid(row=0,    \
                                          column=0, \
                                          padx=3,   \
                                          pady=3)
    mbox = ttk.Combobox(F0,textvariable=self.method,        \
                    values=methods_list,             \
                    width=20)
    mbox.grid(row=0,        \
              column=1,     \
              columnspan=2, \
              padx=3,       \
              pady=3)
    mbox.bind('<<ComboboxSelected>>',lambda e: self.experiment_select())

    F0.grid()


    F1 = ttk.Frame(master,padding=5)
    ttk.Label(F1,text='Truth trajectory filename').grid(row=0,    \
                                                        column=0, \
                                                        padx=3,   \
                                                        pady=5,   \
                                                        sticky='e')
    ttk.Entry(F1,textvariable=self.tfile,width=40).grid(row=0,        \
                                                        column=1,     \
                                                        columnspan=4, \
                                                        padx=3,       \
                                                        pady=5)
    ttk.Button(F1,text='Select',command=lambda: self.select_in(self.tfile)).   \
                                                   grid(row=0,              \
                                                   column=5,                \
                                                   padx=3,                  \
                                                   pady=5)
  

    ttk.Label(F1,text='Experiment (noise_XX_sampl_XXX)').grid(row=1, \
                                                        column=0,    \
                                                        padx=3,      \
                                                        pady=5,      \
                                                        sticky='e')
    wexp = ttk.Entry(F1,textvariable=self.exper,width=20)
    wexp.grid(row=1,        \
              column=1,     \
              columnspan=4, \
              sticky='w',   \
              padx=3,       \
              pady=5)
    wexp.bind('<Return>',lambda e: self.experiment_select())

    ttk.Label(F1,text='Observation filename').grid(row=2,    \
                                                        column=0, \
                                                        padx=3,   \
                                                        pady=5,   \
                                                        sticky='e')
    ttk.Entry(F1,textvariable=self.ofile,width=40).grid(row=2,        \
                                                        column=1,     \
                                                        columnspan=4, \
                                                        padx=3,       \
                                                        pady=5)
    ttk.Button(F1,text='Select',command=lambda: self.select_in(self.ifile)).   \
                                                   grid(row=2,              \
                                                   column=5,                \
                                                   padx=3,                  \
                                                   pady=5)


    ttk.Label(F1,text='Analysis filename').grid(row=3,    \
                                                        column=0, \
                                                        padx=3,   \
                                                        pady=5,   \
                                                        sticky='e')
    ttk.Entry(F1,textvariable=self.ifile,width=40).grid(row=3,        \
                                                        column=1,     \
                                                        columnspan=4, \
                                                        padx=3,       \
                                                        pady=5)
    ttk.Button(F1,text='Select',command=lambda: self.select_in(self.ifile)).   \
                                                   grid(row=3,              \
                                                   column=5,                \
                                                   padx=3,                  \
                                                   pady=5)



    ttk.Button(F1,text='Draw',command=self.make_plot).            \
                                                   grid(row=4,    \
                                                   column=4,      \
                                                   padx=3,        \
                                                   pady=5)
    ttk.Button(F1,text='Close',command=quit).           \
                                         grid(row=4,    \
                                         column=5,      \
                                         padx=3,        \
                                         pady=5)
    F1.grid()


  def experiment_select(self):
  # --------------------------
    self.ofile.set('obs_'+self.exper.get()+'.dat')
    self.ifile.set(self.mdict[self.method.get()]+self.exper.get()+'.dat')
  

  def make_plot(self):
  # ------------------

    # Read the data:
    tt = DataStruct()
    tt.Read(self.tfile.get())

    oo = DataStruct()
    oo.Read(self.ofile.get())

    ii = DataStruct()
    ii.Read(self.ifile.get())

    # Create three plots: Each variable (x,y,z) has its own plot.

    try:
      a = np.array(tt.x)
      b = np.array(ii.x)
      c = np.corrcoef(a,b)
      r = c[0,1]
      rms = np.sqrt(np.mean(np.square(a-b)))
      s = 'r={: 5.2f}, RMS={: 5.2f}'.format(r,rms)
      print('r   = ', r)
      print('RMS = ', rms)
    except:
      s = None
      print('No correlation calculated')

    tag = self.mdict[self.method.get()]+self.exper.get()+'_x'
    self.figx = plt.figure(tag)
    ax = self.figx.add_subplot(111)
    ax.plot(tt.time,tt.x,'--',color='red',linewidth=1,label='Truth')
    ax.plot(oo.time,oo.x,'o',color='red',ms=6,label='Obs')
    ax.plot(ii.time,ii.x,'-',color='red',linewidth=2,label=self.method.get())
    ax.set_ylabel('x',fontweight='bold',fontsize=14)
    ax.set_xlabel('time',fontweight='bold',fontsize=14)
    ax.set_xlim(tt.time[0],tt.time[-1])
    ax.set_ylim([-30,30])
    ax.legend(ncol=3)
    if s is not None:
      self.figx.text(0.12,0.025,s)
    ax.grid(True)


    try:
      a = np.array(tt.y)
      b = np.array(ii.y)
      c = np.corrcoef(a,b)
      r = c[0,1]
      rms = np.sqrt(np.mean(np.square(a-b)))
      s = 'r={: 5.2f}, RMS={: 5.2f}'.format(r,rms)
      print('r   = ', r)
      print('RMS = ', rms)
    except:
      s = None
      print('No correlation calculated')

    tag = self.mdict[self.method.get()]+self.exper.get()+'_y'
    self.figy = plt.figure(tag)
    ay = self.figy.add_subplot(111)
    ay.plot(tt.time,tt.y,'--',color='blue',linewidth=1,label='Truth')
    ay.plot(oo.time,oo.y,'o',color='blue',ms=6,label='obs')
    ay.plot(ii.time,ii.y,'-',color='blue',linewidth=2,label=self.method.get())
    ay.set_ylabel('y',fontweight='bold',fontsize=14)
    ay.set_xlabel('time',fontweight='bold',fontsize=14)
    ay.set_xlim(tt.time[0],tt.time[-1])
    ay.set_ylim([-30,30])
    ay.legend(ncol=3)
    if s is not None:
      self.figy.text(0.12,0.025,s)
    ay.grid(True)


    try:
      a = np.array(tt.z)
      b = np.array(ii.z)
      c = np.corrcoef(a,b)
      r = c[0,1]
      rms = np.sqrt(np.mean(np.square(a-b)))
      s = 'r={: 5.2f}, RMS={: 5.2f}'.format(r,rms)
      print('r   = ', r)
      print('RMS = ', rms)
    except:
      s = None
      print('No correlation calculated')

    tag = self.mdict[self.method.get()]+self.exper.get()+'_z'
    self.figz = plt.figure(tag)
    az = self.figz.add_subplot(111)
    az.plot(tt.time,tt.z,'--',color='green',linewidth=1,label='Truth')
    az.plot(oo.time,oo.z,'o',color='green',ms=6,label='obs')
    az.plot(ii.time,ii.z,'-',color='green',linewidth=2,label=self.method.get())
    az.set_xlabel('time',fontweight='bold',fontsize=14)
    az.set_ylabel('z',fontweight='bold',fontsize=14)
    az.set_xlim(tt.time[0],tt.time[-1])
    az.set_ylim([0,50])
    az.legend(ncol=3)
    if s is not None:
      self.figz.text(0.12,0.025,s)
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
  print(' ====================')
  print(' > results_display.py')

  try:
    tfile = sys.argv[1]
    try:
      exper = sys.argv[2]
    except:
      exper = ''
      
  except:
    tfile = ''
    exper = ''

  root = tk.Tk()
  root.title('DISPLAY_RESULTS')
  root.grid_rowconfigure(0,weight=0)
  root.resizable(width=False,height=True)
  root.protocol('WM_DELETE_WINDOW',quit)

  WinPlot(root,tfile,exper)
  root.mainloop()


if __name__ == '__main__':
  main()
