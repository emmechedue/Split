import numpy as np
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

# Here I just define the constants!
p=10
s=0.05
K=100
M=300

t, N_mean, x_mean , M= np.loadtxt("./output.txt", usecols=(0,1,2,3), unpack=True)
N_table=np.loadtxt("./ensambleN.txt")
x_table=np.loadtxt("./ensamblex.txt")

#parts of the output filename
name = 'scattering_propagule_'
filetype = '.png'

#Parameters of the simulation:
T = 3
step_size = 0.01
M=4
steps = int(math.floor(T/step_size)) + 1

N_max = np.max(N_table)

#Now I start creating the pictures:

for i in range(steps):
	nullfmt   = NullFormatter()         # no labels
	# start with a rectangular Figure, i.e. the whole figure
	plt.figure(1, figsize=(10,10))

	ax = subplot(111)
	subplots_adjust(left=0.1, bottom=0.25)

	# the scatter plot
	plt.plot(x_mean[i],N_mean[i], 'k_',ms=16,mew=3)
	plt.plot(x_mean[i],N_mean[i], 'k|',ms=16,mew=3)
	scatter(x_table[i], N_table[i], s = 30,marker='o',  vmin=0, vmax=1)
	plt.xlabel('$x$, Percentage of cooperators', fontsize = 16)
	plt.ylabel('$N$, Number of individuals', fontsize = 16)
	plt.axis([0, 1, 0, N_max*1.1])#np.max(t)
	plt.title("Evolution of groups for the propagule model")
	
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



