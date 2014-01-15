import numpy
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.mlab as mlab
from matplotlib.ticker import NullFormatter
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.lines as mlines
from pylab import *
from numpy import array
from configobj import ConfigObj
from subprocess import call
import os
from scipy import interpolate
from scipy.interpolate import griddata


# As first thing, I have to manually set the interval for the p, pmax , wich model I am using and how many lines I have in each file
MODELUSED = "Random Splitting"
MODELUSED1= "random_splitting"
PInterval = 1
PMAX= 10.0
NUMLINES=151 # The script will run way faster in this way. An alternative would be to create a matrix by adding line by line using the append method

# Now let's compute how many files I have and create the matrix
dummy = PMAX/PInterval
#NUMFILES = int (math.floor(dummy)+1)
NUMFILES = 29
NUMROWS = NUMFILES*NUMLINES
A = numpy.zeros(shape=(NUMROWS,5)) # This is creating a matrix with NUMROWS rows and 5 columns. The columns will be s,p,<x>, <T> and the checkvalue

#**************** NOW LET'S FILL THE MATRIX ***************************
# Here I will basically just attach the files one after the other adding as second column the value of p!
l=0
for i in range(NUMFILES):
	#p = i*PInterval
	if i==0:
		p=0
	if (i>0 and i<=10):
		p=i*0.01
	if (i>10 and i<19):
		p=(i-9)*0.1
	if i>=19:
		p=i-18
	nameoffile= "p+%s" % p
	b=numpy.loadtxt(nameoffile)
	for j in range(NUMLINES):
		A[l][0] = b[j][0] #That's the s
		A[l][1] = p #That's p
		A[l][2] = b[j][1] # That's <x>
		A[l][3] = b[j][2] # That's <T>
		A[l][4] = b[j][3] # That's the checkvalue
		l += 1
		
	

# Now everything is loaded in the big matrix A shaped in the form s p <x> <T>
# Now I just transpose the matrix
data=A.transpose()

# Now let's try to interpolate the data that I Have to create a better plot.
x = data[0] #That's s
y = data[1] #That's p
z = data[2] #That's <x>

xtemp=numpy.linspace(0,1,1001)
ytemp=numpy.linspace(0,10,1001)
xgrid, ygrid = meshgrid(xtemp, ytemp)

extrapolated = griddata((x,y),z,(xgrid,ygrid),method='cubic')


figure(num=None, dpi=200, facecolor='w', edgecolor='k')
axis([-0.02, 1.02, -0.2, PMAX+0.2])
title("Scatter plot of <x> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
scatter(xgrid,ygrid,s=10,c=extrapolated, marker = 's',  cmap = cm.jet , lw = 0, vmin=0,vmax=1);
plt.colorbar()
#plt.show()
plt.savefig("no_white_phase_diagram"+".png",dpi=200)


# Now let's plot the same but only in the lower left corner

figure(num=None, dpi=200, facecolor='w', edgecolor='k')
axis([-0.02, 0.52, -0.05, 1.05])
title("Scatter plot of <x> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
scatter(xgrid,ygrid,s=10,c=extrapolated, marker = 's',  cmap = cm.jet , lw = 0, vmin=0,vmax=1);
plt.colorbar()
#plt.show()
plt.savefig("no_white_phase_diagram_detailed"+".png",dpi=200)

