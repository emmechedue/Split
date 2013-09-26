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

#Here I manually have to write the parameter p, the model I am using, the interval for s that I used in the simulations and how many jobs I have to analyze
P=10
MODEL="Random Splitting"
NJOBS=121
intervalofs=0.05

#*******************************************************************************************************************************
# As first thing, I have to create the file where I want to write everything


nameoftheglobalfile="p+%s" % P
filefordiagram=open(nameoftheglobalfile, "w+")
filefordiagram.write("# Data to create the phase diagram for the "+MODEL+" model with p= %s \n# Data are printed in the form s <x> <Tstable> \n" % P)

#***********************************************************************************************************************************
# Now I start the loop where basically I do what I was doing in dotheanalysis but for NJOBS times. Furthermore, inside each loop, I also compute what is the average time to reach a "stable state for each realization in the ensamble.
# The average time is defined in this way: First of all I check if a stable state was reached at all, if so I save this time and then I can give an average at the end of the loop. Note that the code for "No stable state reached" is T=-1


for njobs in range(NJOBS):
#*********************************** I AM INSIDE THE NJOBS LOOP NOW ***********************************************
	# As first thing I compute for which s I am analyzing and I go in the correspondent folder
	esse=njobs*intervalofs
	esse="%.2f" % esse
	cmdstring = "s+%s" % (esse)
	os.chdir(cmdstring)
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

	dummy = config["c"]
	c=float(dummy)

	dummy = config["b"]
	b=float(dummy)

	dummy = config["p"]
	p=float(dummy)

	dummy = config["K"]
	K=float(dummy)

	dummy = config["N_max"]
	N_max=int(dummy)

	dummy = config["M_max"]
	M_max=int(dummy)

	dummy = config["N_loop"]
	N_loop=int(dummy)

	dummy = config["choice"]
	cho=int(dummy)
	if cho==1:
		choice="propagule"
	elif cho==2:
		choice="random splitting"
	else:
		choice="deterministic splitting"
	
	#***************************************************************************************
	
	#*****************************Loading everything and computing all the quantities I need
	# Here I load everything
	a=numpy.loadtxt("./ensambleN.txt")
	ensN=a.transpose()
	a=numpy.loadtxt("./ensamblex.txt")
	ensx=a.transpose()
	time=numpy.loadtxt("./time.txt")

	TMAX=len(time)


	#Here I define all the arrays i will need
	N=range(TMAX)
	x=range(TMAX)
	Nstd=range(TMAX)
	xstd=range(TMAX)
	Nerrorstdplus=range(TMAX)
	xerrorstdplus=range(TMAX)
	Nerrorstdminus=range(TMAX)
	xerrorstdminus=range(TMAX)
	Nerrorabsplus=range(TMAX)
	xerrorabsplus=range(TMAX)
	Nerrorabsminus=range(TMAX)
	xerrorabsminus=range(TMAX)

	# Support arrays
	S=range(N_loop)
	P=range(N_loop)
	Tstablearr=range(N_loop)



	Root=math.sqrt(N_loop)
	A=floor(N_loop*0.05) # A is how many data I have to throw away!
	A=int(A)
	
	# Computing all the arrays
	for i in range(TMAX):
		N[i]=ensN[i].mean()
		x[i]=ensx[i].mean()
		Nstd[i]=ensN[i].std()
		xstd[i]=ensx[i].std()
		Nerrorstdplus[i]=N[i]+1.96*Nstd[i]/Root
		Nerrorstdminus[i]=N[i]-1.96*Nstd[i]/Root
		xerrorstdplus[i]=x[i]+1.96*xstd[i]/Root
		xerrorstdminus[i]=x[i]-1.96*xstd[i]/Root
		for j in range(N_loop):
			S[j]=ensN[i][j]-N[i]
			P[j]=ensx[i][j]-x[i]
		for dummy in range(A):
			l1=0
			Boh1=0
			l2=0
			Boh2=0
			for j in range(N_loop):
				if math.fabs(S[j]) > Boh1:
					Boh1=math.fabs(S[j])
					l1=j
				if math.fabs(P[j]) > Boh2:
					Boh2=math.fabs(S[j])
					l2=j
			S[l1]=0
			P[l2]=0
		S.sort()
		P.sort()
		Nerrorabsplus[i]=S[N_loop-1]+N[i]
		Nerrorabsminus[i]=S[0]+N[i]
		xerrorabsplus[i]=P[N_loop-1]+x[i]
		xerrorabsminus[i]=P[0]+x[i]
	

	# To count how many all cooperators and all defectors I have

	fullc=0
	fulld=0
	others=0

	for i in range(N_loop):
		if ensx[TMAX-1,i] > 0.98: #2% is roughly 6 groups, since I am always using N_loop=300
			fullc=fullc+1
		else:
			if ensx[TMAX-1,i] < 0.02:
				fulld=fulld+1
			else:
				others=others+1
	
	#***************************************************************************
	
	# Now I start to compute the things I need to compute the average time of stabilization, remember that here I have to use a instead of ensx because it's better to have the things organized as x[m][t]
	"""	
	for i in range(N_loop):
		Tstablearr[i]= -1.0 #Here I am setting the value to -1, if the stability is not changed, the value is not going to be changed
		for j in range(TMAX):
			if a[i,j] > 0.98 or a[i,j] < 0.02: #2% is roughly 6 groups, since I am always using N_loop=300
				Tstablearr[i]=j*interval
				break
	
	# Now I have to compute the average, I decided to do that if not all realization reached a stable state then the average will be -1
	Tstableaverage = 0
	for i in range(N_loop):
		if Tstablearr[i] < 0 :
			Tstableaverage = -N_loop #In this way, when I get the average it goes to -1 again!
			break
			
		else : 
			Tstableaverage += Tstablearr[i]
		
	Tstableaverage=Tstableaverage/N_loop
	"""
	checkforreachedstability = 1 #Here I am going to modify the value of this variable in 0 if at least one of the iterations hasn't reached a stable configuration
	for i in range(N_loop):
		Tstablearr[i]= -1.0 #Here I am setting the value to -1, if the stability is not changed, the value is going to be TMAX and checkforreachedstability is going to be set to 1 
		for j in range(TMAX):
			if a[i,j] > 0.98 or a[i,j] < 0.02: #2% is roughly 6 groups, since I am basically always using N_loop=300
				Tstablearr[i]=j*interval
				break
				
		if Tstablearr[i] < 0 : # Here I just check if I have to set Tstablearr[i] to the default value (the maximum) or not
			Tstablearr[i] = TMAX
			checkforreachedstability = 0
			
	
	# Now I have to compute the average, I decided to do that if not all realization reached a stable state then the average will be still computed but I am going to print this information
	Tstableaverage = 0
	for i in range(N_loop):
		Tstableaverage += Tstablearr[i]
		
	Tstableaverage=Tstableaverage/N_loop
		
	#*************************************************************************
	#********************* NOW I PLOT **************************************

	#Preparing the strings for the box
	stringa="K= "+str(K)+"\ns= "+str(s)+"\nc= "+str(c)+"\nb= "+str(b)+"\np= "+str(p)+"\nN_max= "+str(N_max)+"\nM_max= "+str(M_max)+"\nN_loop= "+str(N_loop)
	stringa2="fullc= "+str(fullc)+"\nfulld= "+str(fulld)+"\nothers= "+str(others)


	#x plotting
	figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
	plot(time,x,label="data")
	plot(time,xerrorstdplus,'r--', label="std error +")
	plot(time,xerrorstdminus,'r--', label="std error -")
	plot(time,xerrorabsplus,'k--', label="abs error +")
	plot(time,xerrorabsminus,'k--', label="abs error -")
	axis([0, T, 0, 1])
	title("Ensemble average of <x> vs. t for the "+choice+" model")
	ylabel("<x>")
	xlabel("t")
	legend(("data","std error +","std error -","abs error +","abs error -"))
	text(0, 0.8, stringa, bbox=dict(facecolor='orange', alpha=0.8))
	text(T/9+2, 0.8, stringa2, bbox=dict(facecolor='red', alpha=0.8))
	plt.savefig("x.png",dpi=100)
	plt.close()

	#x plotting with no absolute error
	figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
	plot(time,x,label="data")
	plot(time,xerrorstdplus,'k--', label="std error +")
	plot(time,xerrorstdminus,'k--', label="std error -")
	axis([0, T, 0, 1])
	title("Ensemble average of <x> vs. t for the "+choice+" model")
	ylabel("<x>")
	xlabel("t")
	legend(("data","std error +","std error -"))
	text(T-T/5, x0-0.1, stringa, bbox=dict(facecolor='orange', alpha=0.8))
	text(T-T/20, x0-0.1, stringa2, bbox=dict(facecolor='red', alpha=0.8))
	plt.savefig("x_noabs.png",dpi=100)
	plt.close()


	#N plotting
	figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
	plot(time,N,label="data")
	plot(time,Nerrorstdplus, 'r--',label="std error +")
	plot(time,Nerrorstdminus, 'r--',label="std error -")
	plot(time,Nerrorabsplus,'k--', label="abs error +")
	plot(time,Nerrorabsminus,'k--', label="abs error -")
	title("Ensemble average of <N> vs. t for the "+choice+" model")
	ylabel("<N>")
	xlabel("t")
	legend(("data","std error +","std error -","abs error +","abs error -"))
	text(0, N_max-10, stringa, bbox=dict(facecolor='orange', alpha=0.8))
	plt.savefig("N.png",dpi=100)
	plt.close()


	#N plotting without errors
	figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
	plot(time,N,label="data")
	title("Ensemble average of <N> vs. t for the "+choice+" model")
	ylabel("<N>")
	xlabel("t")
	text(0, N_max-10, stringa, bbox=dict(facecolor='orange', alpha=0.8))
	plt.savefig("N_noerr.png",dpi=100)
	plt.close()
	
	#**********************************************************************************************
	
	#****************** NOW THE LAST DETAILS ********************************
	os.chdir("..") #Let's go back to the general folder
	filefordiagram.write("%.2f" % s)
	filefordiagram.write("%8.3f" % x[TMAX-1])
	filefordiagram.write("%8.3f" % Tstableaverage)
	filefordiagram.write("%8.3f" % checkforreachedstability)
	
	filefordiagram.write("\n")
	
	#********************************* END OF THE SINGLE ITERATION OF THE LOOP

filefordiagram.close()


	
	
	
				
	
	
					

	
	
