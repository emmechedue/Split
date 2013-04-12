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
	choice2="propagule"
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
name = name = "./video/scattering_"+choice2+"_"
filetype = '.png'

PicN=20 #How many pictures do I want per second
delta=int(floor(1/(PicN*interval))) #Here I am computing the delta for each step and the number of steps
steps = int(floor((TMAX-1)/delta))

call(["mkdir", "video"]) #I am creating a new directory to not create confusion
#************************************Now I start creating the pictures:**************

for i in range(steps): #I do all the pictures but the last (the one at T=18)
	nullfmt   = NullFormatter()         # no labels
	# start with a rectangular Figure, i.e. the whole figure
	plt.figure(1, figsize=(10,14))

	subplot(211)
	subplots_adjust(left=0.1, bottom=0.25,top=0.93)

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
	
	#The x_avereage plot
	
	stringa = "K= "+str(K)+"\ns= "+str(s)+"\np= "+str(p)+"\nN_max= "+str(N_max)+"\nM_max= "+str(M_max)
	subplot(212)
	subplots_adjust(left=0.1, bottom=0.07)
	a=numpy.loadtxt("./output.txt")
	data=a.transpose()
	plot(data[0],data[2])
	title("<x> vs. t for the "+choice+" model")
	xlabel("t")
	ylabel("<x>")
	text(0, 0.9, stringa, bbox=dict(facecolor='orange', alpha=0.8))
	plt.plot([t[j], t[j]], [-0.2, 1.2], 'k-', lw=2.0) #Adding the line for the time
	
	
	
	
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
plt.figure(1, figsize=(10,16))

subplot(211)
subplots_adjust(left=0.1, bottom=0.25,top=0.93)

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
	
#The x_average plot
	
stringa = "K= "+str(K)+"\ns= "+str(s)+"\np= "+str(p)+"\nN_max= "+str(N_max)+"\nM_max= "+str(M_max)
subplot(212)
subplots_adjust(left=0.1, bottom=0.07)
a=numpy.loadtxt("./output.txt")
data=a.transpose()
plot(data[0],data[2])
title("<x> vs. t for the "+choice+" model")
xlabel("t")
ylabel("<x>")
text(0, 0.9, stringa, bbox=dict(facecolor='orange', alpha=0.8))
plt.plot([t[j], t[j]], [-0.2, 1.2], 'k-', lw=2.0) #Adding the line for the time

	

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

#*********************NOW I ALSO MAKE THE MOVIE*******************************
os.system("mencoder -really-quiet -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \"mf://./video/*.png\" -mf type=png:fps=10 -o movie.avi")


