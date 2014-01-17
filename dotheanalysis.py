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
	

#To count how many all cooperators and all defectors I have

fullc=0
fulld=0
others=0

for i in range(N_loop):
	if ensx[TMAX-1,i]>=0.99:
		fullc=fullc+1
	else:
		if ensx[TMAX-1,i]<=0.01:
			fulld=fulld+1
		else:
			others=others+1

#***************************************************************************

#********************************************************Now the plotting:

#Preparing the strings for the box
stringa="K= "+str(K)+"\ns= "+str(s)+"\nc= "+str(c)+"\nb= "+str(b)+"\np= "+str(p)+"\nN_max= "+str(N_max)+"\nM_max= "+str(M_max)+"\nN_loop= "+str(N_loop)
stringa2="fullc= "+str(fullc)+"\nfulld= "+str(fulld)+"\nothers= "+str(others)


#x plotting
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
plot(time,x,label="data")
plot(time,xerrorstdplus,'r--', label="confidence interval +")
plot(time,xerrorstdminus,'r--', label="confidence interval -")
plot(time,xerrorabsplus,'k--', label="abs error +")
plot(time,xerrorabsminus,'k--', label="abs error -")
axis([0, T, 0, 1])
title("Ensemble average of <x> vs. t for the "+choice+" model")
ylabel("<x>")
xlabel("t")
legend(("data","confidence interval +","confidence interval -","abs error +","abs error -"))
text(0, 0.8, stringa, bbox=dict(facecolor='orange', alpha=0.8))
text(T/9+2, 0.8, stringa2, bbox=dict(facecolor='red', alpha=0.8))
plt.tight_layout()
plt.savefig("x.png",dpi=100)
plt.close()

#x plotting with no absolute error
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
plot(time,x,label="data")
plot(time,xerrorstdplus,'k--', label="confidence interval +")
plot(time,xerrorstdminus,'k--', label="confidence interval -")
axis([0, T, 0, 1])
title("Ensemble average of <x> vs. t for the "+choice+" model",fontsize=20)
ylabel("<x>",fontsize=20)
xlabel("t",fontsize=20)
legend(("data","confidence interval +","confidence interval -"))
text(0, 0.8, stringa, bbox=dict(facecolor='orange', alpha=0.8))
text(T/9+2, 0.8, stringa2, bbox=dict(facecolor='red', alpha=0.8))
plt.tight_layout()
plt.savefig("x_noabs.png",dpi=100)
plt.close()

#N plotting
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
plot(time,N,label="data")
plot(time,Nerrorstdplus, 'r--',label="confidence interval +")
plot(time,Nerrorstdminus, 'r--',label="confidence interval -")
plot(time,Nerrorabsplus,'k--', label="abs error +")
plot(time,Nerrorabsminus,'k--', label="abs error -")
title("Ensemble average of <N> vs. t for the "+choice+" model")
ylabel("<N>")
xlabel("t")
legend(("data","confidence interval +","confidence interval -","abs error +","abs error -"))
text(0, N_max-10, stringa, bbox=dict(facecolor='orange', alpha=0.8))
plt.tight_layout()
plt.savefig("N.png",dpi=100)
plt.close()


#N plotting without errors
figure(num=None, figsize=(16, 12), dpi=160, facecolor='w', edgecolor='k')
plot(time,N,label="data")
title("Ensemble average of <N> vs. t for the "+choice+" model")
ylabel("<N>")
xlabel("t")
text(0, N_max-10, stringa, bbox=dict(facecolor='orange', alpha=0.8))
plt.tight_layout()
plt.savefig("N_noerr.png",dpi=100)
plt.close()
