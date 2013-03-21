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


#****************************************Let's first of all get all the parameters:
config = ConfigObj("./parameters.txt")

dummy = config["N0"]
N0=int(dummy)

dummy = config["x0"]
x0=float(dummy)

dummy = config["T"]
T=float(dummy)

dummy = config["interval"]
interval=float(dummy)

dummy = config["s"]
s=float(dummy)

dummy = config["p"]
p=float(dummy)

dummy = config["K"]
K=float(dummy)

dummy = config["N_max"]
N_max=int(dummy)

dummy = config["M_max"]
M_max=int(dummy)

dummy = config["choice"]
cho=int(dummy)
if cho==1:
	choice="propagule"
	choice2=choice
else:
	choice="random splitting"
	choice2="random_splitting"
	
#***************************************************************************************

#*****************************Loading everything *************************************
# Here I load everything
t, N_mean, x_mean , M= np.loadtxt("./output.txt", usecols=(0,1,2,3), unpack=True)
N_table=np.loadtxt("./ensambleN.txt")
x_table=np.loadtxt("./ensamblex.txt")

TMAX=len(t)

#**************************************************************************************

#***************************Preparing for making the video****************************
#parts of the output filename
name = "scattering_"+choice2+"_"
filetype = '.png'

PicN=20 #How many pictures do I want per second
delta=int(floor(1/(PicN*interval))) #Here I am computing the delta for each step and the number of steps
steps = int(floor((TMAX-1)/delta))

#************************************Now I start creating the pictures:**************

#************************************Now I start creating the pictures:**************

for i in range(steps): #I do all the pictures but the last (the one at T=18)
	nullfmt   = NullFormatter()         # no labels
	# start with a rectangular Figure, i.e. the whole figure
	plt.figure(1, figsize=(10,10))

	subplot(111)
	subplots_adjust(left=0.1, bottom=0.25)

	j=i*delta
	# the scatter plot
	testo1="M= "+str(int(M[j]))
	plt.plot(x_mean[j],N_mean[j], 'k_',ms=16,mew=3)
	plt.plot(x_mean[j],N_mean[j], 'k|',ms=16,mew=3)
	scatter(x_table[j], N_table[j], s = 30,marker='o',  vmin=0, vmax=1)
	plt.xlabel('$x$, Percentage of cooperators', fontsize = 10)
	plt.ylabel('$N$, Number of individuals', fontsize = 10)
	plt.axis([0, 1, 0, N_max*1.1])
	plt.title("Evolution of groups for the "+choice+" model")
	text(0, N_max, testo1, bbox=dict(facecolor='red', alpha=0.8))
	
	
	# create time diagram
	# definitions for the axes, determines position of scatter plot and histograms, general
	left, width = 0.1, 0.65 #left = position left in % of figsize, dito width
	bottom, height = 0.1, 0.05
	rect_time = [left, bottom,width, height]
	timebox =  plt.axes(rect_time)
	timebox.bar(0, 1, t[i], 0, color = '0.7', orientation = 'horizontal')
	timebox.set_xlim( 0, T )
	timebox.yaxis.set_ticks_position("none")
	timebox.yaxis.set_ticklabels("")
	plt.title('Time')
	
	# save plot as png file
	#stringlist.append(str(i))
	if i<10: ending = '000' + str(i)
	if 10 <= i and i<100: ending = '00' + str(i)
	if 100<= i and i<1000: ending = '0' + str(i)
	if 1000<= i and i<10000: ending = str(i)
	filename = name + ending + filetype
	plt.savefig(filename)
	plt.close()
	print i


#*******************************************The last plot***************
i=steps #steps +1 -1
j=TMAX-1 #TMAX 
nullfmt   = NullFormatter()         # no labels
# start with a rectangular Figure, i.e. the whole figure
plt.figure(1, figsize=(10,10))

subplot(111)
subplots_adjust(left=0.1, bottom=0.25)

# the scatter plot
stringa2="M= "+str(int(M[j]))
plt.plot(x_mean[j],N_mean[j], 'k_',ms=16,mew=3)
plt.plot(x_mean[j],N_mean[j], 'k|',ms=16,mew=3)
scatter(x_table[j], N_table[j], s = 30,marker='o',  vmin=0, vmax=1)
plt.xlabel('$x$, Percentage of cooperators', fontsize = 10)
plt.ylabel('$N$, Number of individuals', fontsize = 10)
plt.axis([0, 1, 0, N_max*1.1])
plt.title("Evolution of groups for the "+choice+" model")
text(0, N_max, stringa2, bbox=dict(facecolor='red', alpha=0.8))
	

	
# create time diagram
# definitions for the axes, determines position of scatter plot and histograms, general
left, width = 0.1, 0.65 #left = position left in % of figsize, dito width
bottom, height = 0.1, 0.05
rect_time = [left, bottom,width, height]
timebox =  plt.axes(rect_time)
timebox.bar(0, 1, t[i], 0, color = '0.7', orientation = 'horizontal')
timebox.set_xlim( 0, T )
timebox.yaxis.set_ticks_position("none")
timebox.yaxis.set_ticklabels("")
plt.title('Time')
	
# save plot as png file
#stringlist.append(str(i))
if i<10: ending = '000' + str(i)
if 10 <= i and i<100: ending = '00' + str(i)
if 100<= i and i<1000: ending = '0' + str(i)
if 1000<= i and i<10000: ending = str(i)
filename = name + ending + filetype
plt.savefig(filename)
plt.close()
print i

