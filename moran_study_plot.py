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


EMME=1000 # This is M_max
LENPOINT=6 # This is how many points per p I have (i.e. how many n's)


p=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10,0.20,0.50,1.0,1.50]

#****************** Note that to solve the problem with n 1 and 2 that don't have all the points, I just fill the text files with negative bullshit so everything has the same lenght**************

LENGP=len(p) #For how many different p's I want to do that


data=numpy.empty((LENPOINT,LENGP))
points=numpy.arange(1,LENPOINT+1)
moran=numpy.empty((LENPOINT,LENGP)) #That's the matrix with the computable moran values. LENGP is the number of rows and LENPOINT of columns
error=numpy.empty((LENPOINT,LENGP)) #That's the matrix of the error points for the data matrix

#************Filling the matrix*****************
temp=numpy.empty((LENGP,2)) #This is where I charge the n+i.txt file. It has 2 rows cuz x and error!

for j in range(LENPOINT):
	i=j+1
	tempstring='./n+'+str(i)+'.txt'
	temp=numpy.loadtxt(tempstring)
	temp2=numpy.transpose(temp)
	data[j]=temp2[0]
	error[j]=temp2[1]

	
#*******************************************************************************************************************
#fill the moran matrix
for j in range(LENPOINT):
	for i in range(LENGP):
		l=1+j
		r=1+p[i]
		if r<2:
			moran[j][i]=(1-1/(r**l))/(1-1/(r**EMME))
			
		else :
			moran[j][i]=1-1/(r**l)
		

colors = ['b','g','r','c','m','y','k'] #list the colors
#colors = plt.get_cmap('jet')(np.linspace(0, 1.0, LENPOINT))


#****************************HERE I CONSTRUCT THE PLOT THAT HAS P AS VARIABLE AND N AS PARAMETER!! **************************

for j in range(LENPOINT):

	plt.plot(p, moran[j], color=colors[j],marker='^',linestyle='None')
	plt.errorbar(p, data[j], yerr=error[j],color=colors[j],marker='o',linestyle='None')
	#plt.show()



plt.title("Fixation probabilities as a function of r and i")
plt.legend(('Moran process for n=1','Simulation for n=1','Moran process for n=2','Simulation for n=2','Moran process for n=3','Simulation for n=3','Moran process for n=4','Simulation for n=4'))
plt.xlabel('p')
plt.ylabel('fixation probability')
#yscale('log')
#xscale('log')
axis([-0.05, 1.05, -0.05, 1.05])
plt.tight_layout()
plt.show()


#*********************HERE I INSTEAD PLOT A SINGLE N AT THE TIME*********************

for j in range(LENPOINT):
	
	figure(num=None, figsize=(21, 10), dpi=180, facecolor='w', edgecolor='k')
	plt.plot(p, moran[j], color=colors[j],marker='^',linestyle='None')
	plt.errorbar(p, data[j], yerr=error[j],color=colors[j],marker='o',linestyle='None')
	plt.title("Fixation probabilities as a function of r and i",fontsize=20)
	legend1='Moran process for $n_0=$'+str(j+1)
	legend2='Simulation for $n_0=$'+str(j+1)
	legend=plt.legend((legend1,legend2),loc=4)
	legend.get_title().set_fontsize('6') #legend 'Title' fontsize
	plt.setp(gca().get_legend().get_texts(), fontsize='20') #legend 'list' fontsize
	plt.xlabel('p',fontsize=20)
	plt.ylabel('fixation probability',fontsize=20)
	maxi=numpy.amax(data[j])
	mini=numpy.amin(moran[j])
	axis([0.008, 1.6, mini-0.05, maxi + 0.05])
	#plt.ylim([mini-0.05, maxi + 0.05])
	plt.tight_layout()
	#yscale('log')
	xscale('log')
	#plt.show()
	savename='Moran_'+str(j+1)+'.png'
	plt.savefig(savename,dpi=100)


#*****************************************************************************************************************************
#**************************************************************************************************************************************************
#*****************************************************************************************************************************

#*********************************	NOW LET'S TRY TO SEE WHAT HAPPENS IF I PLOT THINGS AS A FUNCTION OF N WITH P AS PARAMETERS	******************
#Let's create all the the matrices

datat=np.transpose(data)
morant=np.transpose(moran)
errort=np.transpose(error)
points=numpy.arange(1,LENPOINT+1)

colors = cm.rainbow(np.linspace(0, 1, len(p)))
for j in range(LENGP):

	figure(num=None, figsize=(21, 10), dpi=180, facecolor='w', edgecolor='k')
	plt.plot(points, morant[j],color=colors[j], marker='^',linestyle='None')
	plt.errorbar(points, datat[j], color=colors[j], yerr=errort[j],marker='o',linestyle='None')
	plt.title("Fixation probabilities as a function of r and i")
	legend1='Moran process for p='+str(p[j])
	legend2='Simulation for p='+str(p[j])
	plt.legend((legend1,legend2),loc=4)
	plt.xlabel('n')
	plt.ylabel('fixation probability')
	maxi=numpy.amax(datat[j])
	mini=numpy.amin(morant[j])
	axis([0.8, LENPOINT+0.2, mini-0.05, maxi + 0.05])
	plt.tight_layout()
	savename='./p_plots/Moran_'+str(p[j])+'.png'
	plt.savefig(savename,dpi=100)
	
	






















