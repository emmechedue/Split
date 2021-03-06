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


# As first thing, I have to manually set the interval for the p, pmax , wich model I am using and how many lines I have in each file
MODELUSED = "Random Splitting"
MODELUSED1= "random_splitting"
PInterval = 1
PMAX= 1.0
NUMLINES=121 # The script will run way faster in this way. An alternative would be to create a matrix by adding line by line using the append method

# Now let's compute how many files I have and create the matrix
dummy = PMAX/PInterval
NUMFILES = int (math.floor(dummy)+1)
NUMROWS = NUMFILES*NUMLINES
A = numpy.zeros(shape=(NUMROWS,5)) # This is creating a matrix with NUMROWS rows and 5 columns. The columns will be s,p,<x>, <T> and the checkvalue

#**************** NOW LET'S FILL THE MATRIX ***************************
# Here I will basically just attach the files one after the other adding as second column the value of p!
l=0
for i in range(NUMFILES):
	p = i*PInterval
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

#********************************************************************************
#************************** NOW THE PLOTTING ***************************

# Let's define a new color map
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
		   (0.25,1.0, 1.0),
                   (0.75, 0.5, 0.5),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}
my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)

# First let's plot the plot for <x>
figure(num=None, dpi=260, facecolor='w', edgecolor='k')
title("Scatter plot of <x> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
x = data[0]
y = data[1]
z = data[2]
axis([-0.1, 6.1, -0.5, PMAX+0.5])
scatter(x,y,s=20,c=z, marker = 's',  cmap = cm.jet, lw = 0 );
plt.colorbar()
plt.savefig("x_phase_diagram_"+MODELUSED1+".png",dpi=200)


#Second let's plot the one for <T>
figure(num=None, dpi=260, facecolor='w', edgecolor='k')
title("Scatter plot of <T> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
x = data[0]
y = data[1]
z = data[3]
axis([-0.1, 6.1, -0.5, PMAX+0.5])
scatter(x,y,s=20,c=z, marker = 's',  cmap = cm.jet , lw = 0);
plt.colorbar()
plt.savefig("T_phase_diagram_"+MODELUSED1+".png",dpi=200)

#**********	Now let's create the plot for <x> with the markers to signal that the time was not reached	************
# <x>
figure(num=None, dpi=260, facecolor='w', edgecolor='k')
title("Scatter plot of <x> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
x = data[0]
y = data[1]
z = data[2]
axis([-0.1, 6.1, -0.5, PMAX+0.5])
scatter(x,y,s=20,c=z, marker = 's',  cmap = cm.jet, lw = 0 );
for j in range(NUMROWS):
	if data[4][j]<1:
		plt.plot(data[0][j],data[1][j], 'k_',ms=2.5,mew=0.6)
		plt.plot(data[0][j],data[1][j], 'k|',ms=4,mew=0.6)
	

plt.colorbar()
plt.savefig("x_phase_diagram_with_marks_"+MODELUSED1+".png",dpi=200)

# <T>
figure(num=None, dpi=260, facecolor='w', edgecolor='k')
title("Scatter plot of <T> as a function of s and p \n for the "+MODELUSED+" model",fontsize=14)
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
x = data[0]
y = data[1]
z = data[3]
axis([-0.1, 6.1, -0.5, PMAX+0.5])
scatter(x,y,s=20,c=z, marker = 's',  cmap = cm.jet , lw = 0);
for j in range(NUMROWS):
	if data[4][j]<1:
		plt.plot(data[0][j],data[1][j], 'k_',ms=2.5,mew=0.6)
		plt.plot(data[0][j],data[1][j], 'k|',ms=4,mew=0.6)
	

plt.colorbar()
plt.savefig("T_phase_diagram_with_marks_"+MODELUSED1+".png",dpi=200)
"""cmap = cm.jet"""
"""
#fig = plt.figure(figsize=(6,6))
#ax = fig.add_subplot(111)
#figure(num=None, figsize=(8, 6), dpi=160, facecolor='w', edgecolor='k')
figure(num=None, dpi=160, facecolor='w', edgecolor='k')
title("X vs Y AVG",fontsize=14)
xlabel("XAVG",fontsize=12)
ylabel("YAVG",fontsize=12)
#ax.grid(True,linestyle='-',color='0.75')
x = np.random.random(30)
y = np.random.random(30)
z = np.random.random(30)

# scatter with colormap mapping to z value
scatter(x,y,s=20,c=z, marker = 'o', cmap = cm.jet );
plt.colorbar()
"""
