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

