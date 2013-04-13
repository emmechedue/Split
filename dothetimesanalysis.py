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
else:
	choice="random splitting"
	
#***************************************************************************************

#*****************************Loading everything *************************************
# Here I load everything
t, N_mean, x_mean , M= np.loadtxt("./output.txt", usecols=(0,1,2,3), unpack=True)
N_table=np.loadtxt("./ensambleN.txt")
x_table=np.loadtxt("./ensamblex.txt")

TMAX=len(t)


#**************************************************************************************

#*****************************Now I will search for all the first times******************
#Searching for Tm
Tm=-1

for i in range(TMAX):
	if M[i]==M_max:
		Tm=t[i]
		break

#Searching for Tg and Tg_one
Tg=-1
Tg_one=-1
boolcheck= False
dummy=0

for i in range(TMAX): #Here I search for Tg_one
	for j in range(M_max):
		if x_table[i][j] in (0,1):
			Tg_one=t[i]
			boolcheck=True
			dummy=i
			break
	if boolcheck==True :
		break



for i in range(dummy,TMAX): #Here I start from dummy to search for Tg
	boolcheck=True
	for j in range(M_max):
		if x_table[i][j] not in (0,1):
			boolcheck=False
			break
	if boolcheck==True:
		Tg=t[i]
		break
	
#Searching for Tp
Tp=-1
check1=False
check0=False

for i in range(TMAX):
	check1=True
	check0=True
	for j in range(M_max):
		if x_table[i][j]!=0 :
			check0=False
		if x_table[i][j]!=1 :
			check1=False
		if check1==False and check0==False :
			break
	if check1==True or check0==True :
		Tp=t[i]
		break


#*************************************************************************************

#***********************Plotting**************************************

stringa = "K= "+str(K)+"\ns= "+str(s)+"\np= "+str(p)+"\nN_max= "+str(N_max)+"\nM_max= "+str(M_max)

#Plot for N
figure(num=None, figsize=(12, 9), dpi=160, facecolor='w', edgecolor='k')
plt.plot(t,N_mean)
plt.title("<N> vs. t for the "+choice+" model")
xlabel("t")
ylabel("<N>")
text(0, N_max-10, stringa, bbox=dict(facecolor='orange', alpha=0.8))
plt.savefig("N.png",dpi=100)
plt.close()

#Plot for x
figure(num=None, figsize=(12, 9), dpi=160, facecolor='w', edgecolor='k')
plot(t,x_mean)
title("<x> vs. t for the "+choice+" model")
xlabel("t")
ylabel("<x>")
text(0, 1.05, stringa, bbox=dict(facecolor='orange', alpha=0.8))
if Tm!=-1:
	plt.plot([Tm, Tm], [-0.2, 1.2], 'r-', lw=1.2) #Adding the line for Tm
	plt.annotate("Tm= "+str(Tm), xy=(Tm, 0), xytext=(Tm-3, -0.05),arrowprops=dict(facecolor="red", shrink=0.05, width=2.5),)
if Tg!=-1:
	plt.plot([Tg, Tg], [-0.2, 1.2], 'g-', lw=1.2) #Adding the line for Tg
	plt.annotate("Tg= "+str(Tg), xy=(Tg, 0), xytext=(Tg-3, -0.15),arrowprops=dict(facecolor="green", shrink=0.05, width=2.5),)
if Tg_one!=-1:
	plt.plot([Tg_one, Tg_one], [-0.2, 1.2], 'y-', lw=1.2) #Adding the line for Tg_one
	plt.annotate("Tg_one= "+str(Tg_one), xy=(Tg_one, 0), xytext=(Tg_one-3, -0.15),arrowprops=dict(facecolor="yellow", shrink=0.05, width=2.5),)
if Tp!=-1:
	plt.plot([Tp, Tp], [-0.2, 1.2], 'm-', lw=1.2) #Adding the line for Tp
	plt.annotate("Tp= "+str(Tp), xy=(Tp, 0), xytext=(Tp+3, -0.1),arrowprops=dict(facecolor="magenta", shrink=0.05, width=2.5),)
plt.savefig("x.png",dpi=100)
plt.close()

#*******************************************************************************

#**********************Making the histogram*****************************
#Counting everything
Number=400 #This is the numbers of BINS in my histogram, unless TMAX is perfectly divisible by Number-1, then i simply shift Number by 1, to get one less bin and it's done
if TMAX%(Number-1)==0:
	Number=Number+1

I=floor(TMAX/(Number-1)) #I'm computing the time steps for my histogram.
I=int(I)
middle=floor(I/2)#Middle is the time point that I take to represnt the interval
middle=int(middle)
Nu=range(Number) #This array will count how many cells get fully cooperative or fully defective in that time interval, the time intervals are 0,I,2I,.....,(Number-1)I,TMAX

for i in range(Number):
	Nu[i]=0

for j in range(M_max): #computing Nu[0] at time t[middle]
	if x_table[middle][j] in (0,1):
		Nu[0]=Nu[0]+1

Cumulative=Nu[0]

for i in range(1,Number-1): #computing from Nu[1] to Nu[Number-2]
	for j in range(M_max):
		if x_table[I*i+middle][j] in (0,1):
			Nu[i]=Nu[i]+1
	Nu[i]=Nu[i]-Cumulative
	if Nu[i]<0:
		Nu[i]=0
	Cumulative=Cumulative+Nu[i]
	

middle2=floor((TMAX-I*Number+I)/2) #TMAX - (Number-1)*I, here I am computing the last middle point (the last interval is shorter)
middle2=int(middle2)
i=(Number-1)*I+middle2
for j in range(M_max): #computing Nu[0] at time t[middle]
	if x_table[i][j] in (0,1):
		Nu[Number-1]=Nu[Number-1]+1

Nu[Number-1]=Nu[Number-1]-Cumulative

#Preparing the array for the histogram
k=0
Tau=range(M_max)
for i in range(M_max):
	Tau[i]=0

for i in range(Number-1): #filling the array for all the middle points except for the last one
	timeindex=I*i+middle
	for j in range(Nu[i]):
		Tau[k]=t[timeindex]
		k=k+1

timeindex=I*(Number-1)+middle2
for i in range(Nu[Number-1]):
	Tau[k]=t[timeindex]
	k=k+1
	
#Plotting
maximus=max(Nu)
figure(num=None, figsize=(12, 9), dpi=160, facecolor='w', edgecolor='k')
stringa="Tm= "+str(Tm)+"\nTg= "+str(Tg)+"\nTg_one= "+str(Tg_one)+"\nTp= "+str(Tp)
plt.hist(Tau, Number, normed=0, facecolor='yellow', alpha=0.75)
text(Tg_one, maximus-maximus/7, stringa, bbox=dict(facecolor='orange', alpha=0.8))
plt.xlabel("t")
plt.title("Histogram of group fixation times")
plt.savefig("T_histogram.png",dpi=100)
plt.close()


#*******************Times legend********************************************
# Tm is when for the first time my colony reaches M_max
# Tg is when all groups reach fixation (i.e. every group has or x=0 or x=1)
# Tg_one is when the first group fixates
# Tp is when the entire population fixates
