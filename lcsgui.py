import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import Tkinter as tk
import h5py as hp
from mpl_toolkits.basemap import Basemap
import scipy.ndimage.filters as filter
timestep = 0
originlat = 36.7964
originlon =-120.822
xres = 2000
yres = 2000
xstep = 0.3
ystep = 0.3
#Button Functions-----------------------------------------------------------
def _nextTime(time):
    global ridges
    for c in ridges.collections:
        c.remove()
    global timestep
    timestep = time + 1
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    quad.set_array(np.ravel(np.transpose(np.squeeze(ftledata[timestep,:,:]))))
    ax.set_title("Time Step {0}, Unfiltered".format(timestep))
    timeentry.delete(0,tk.END)
    timeentry.insert(tk.END,timestep)
    canvas.draw()

def _jumpTime():
    global ridges
    for c in ridges.collections:
        c.remove()
    newtimestep = int(timeentry.get())
    global timestep
    timestep = newtimestep
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    quad.set_array(np.ravel(np.transpose(np.squeeze(ftledata[timestep,:,:]))))
    ax.set_title("Time Step {0}, Unfiltered".format(timestep))
    timeentry.delete(0,tk.END)
    timeentry.insert(tk.END,timestep)
    canvas.draw()

def _runFilter():
    global ridges
    for c in ridges.collections:
        c.remove()
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    filterdata = np.squeeze(ftledata[timestep,:,:])
    maxits = int(iterationentry.get())
    print maxits
    iterations = 0
    stddev = float(deviationentry.get())
    print stddev
    while iterations < maxits:
        print iterations
        filterdata = filter.gaussian_filter(filterdata,sigma=stddev)
        iterations += 1
    quad.set_array(np.ravel(np.transpose(filterdata[:,:])))
    canvas.draw()
    
def _extractLCS():
    global ridges
    for c in ridges.collections:
        c.remove()
    filterdata = np.squeeze(ftledata[timestep,:,:])
    maxits = int(iterationentry.get())
    print maxits
    iterations = 0
    stddev = float(deviationentry.get())
    print stddev
    while iterations < maxits:
        print iterations
        filterdata = filter.gaussian_filter(filterdata,sigma=stddev)
        iterations += 1
    dx, dy = np.gradient(filterdata,xstep,ystep)
    dxdx, dydx = np.gradient(dx,xstep,ystep)
    dxdy, dydy = np.gradient(dy,xstep,ystep) 
    dirdiv = np.empty([dim[1],dim[2]])
    minimumeig = np.empty([dim[1],dim[2]])
    for i in range(dim[1]):
        for j in range(dim[2]):
            eig = np.linalg.eig([[dxdx[i,j],dxdy[i,j]],[dydx[i,j],dydy[i,j]]])
            eigmin =  np.argmin(eig[0])
            dirdiv[i,j] = np.dot(eig[1][:,eigmin],[dx[i,j],dy[i,j]])
            minimumeig[i,j] = eig[0][eigmin]    
    tol = float(thresholdentry.get())
    potridge = np.ma.masked_where(minimumeig>=tol,dirdiv)
    ridges = m.contour(xx, yy, np.transpose(potridge),levels =[0],colors='orange')
    canvas.draw()
    

#Main -----------------------------------------------------------------------
root = tk.Tk()
root.wm_title("Embedding in TK")
#root.filename = tkFileDialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
#print (root.filename)
#Interface Objects-----------------------------------------------------------
timelabel = tk.Label(root, text="Time Step")
timeentry = tk.Entry(root)
timeentry.insert(tk.END,timestep)

deviationlabel = tk.Label(root, text="Filter Standard Deviation")
deviationentry = tk.Entry(root)
deviationentry.insert(tk.END,0)

iterationlabel = tk.Label(root, text="Filter Iteration")
iterationentry = tk.Entry(root)
iterationentry.insert(tk.END,0)

thresholdlabel = tk.Label(root, text="Eigenvalue Threshold")
thresholdentry = tk.Entry(root)
thresholdentry.insert(tk.END,0)

nextbutton = tk.Button(root, text='Next Step', command= lambda: _nextTime(timestep))
jumpbutton = tk.Button(root, text='Jump to Time', command= lambda: _jumpTime())
filterbutton = tk.Button(root, text='Run Filter', command= lambda: _runFilter())
extractbutton = tk.Button(root, text='Extract LCS', command= lambda: _extractLCS())#attdata,aridges))

#Object Placement-----------------------------------------------------------
timelabel.grid(row=4,column=0,sticky="e")
timeentry.grid(row=4,column=1,padx=10)
deviationlabel.grid(row=5,column=0,sticky="e")
deviationentry.grid(row=5,column=1,padx=10)
iterationlabel.grid(row=6,column=0,sticky="e")
iterationentry.grid(row=6,column=1,padx=10)
thresholdlabel.grid(row=7,column=0,sticky="e")
thresholdentry.grid(row=7,column=1,padx=10)
jumpbutton.grid(row=8,column=0,sticky="e")
nextbutton.grid(row=8,column=1,sticky="w",padx=10)
filterbutton.grid(row=9,column=1,sticky="w",padx=10)
extractbutton.grid(row=10,column=1,sticky="w",padx=10)
# Adding graph to Canvas--------------------------------------------
f=hp.File('attFTLEOutput.mat','r')
ftledata = f['F'][:,:,:]
f.close()

dim = ftledata.shape
fig = Figure(figsize=(8, 6), dpi=100)
ax = fig.add_subplot(111)

m = Basemap(width=xres*1000,height=yres*1000,\
    rsphere=(6378137.00,6356752.3142),\
    resolution='i',area_thresh=100.,projection='lcc',\
    lat_1=35.,lat_0=originlat,lon_0=originlon,ax=ax)

m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(20,50,2.5),labels=[True,False,False,False])
m.drawmeridians(np.arange(-135,-105,2.5),labels=[False,False,False,True])
m.drawstates()
x = np.linspace(0, m.urcrnrx, dim[1])
y = np.linspace(0, m.urcrnry, dim[2])
xx, yy = np.meshgrid(x, y)

quad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(ftledata[timestep,:,:])),shading='gouraud')
ax.set_title('Time Step {0}, Unfiltered'.format(timestep))
ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])

canvas = FigureCanvasTkAgg(fig,root)
canvas.show() 
canvas.get_tk_widget().grid(row=1,column=3,rowspan=40)


toolbar_frame = tk.Frame(root)
toolbar_frame.grid(row=0,column=3,sticky="w")
toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame) 
toolbar.update() 

#root.mainloop()
tk.mainloop()