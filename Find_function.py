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
NUMLINES=151 # The script will run way faster in this way. An alternative would be to create a matrix by adding line by line using the append method

# Now let's compute how many files I have and create the matrix
#dummy = PMAX/PInterval
#NUMFILES = int (math.floor(dummy)+1)
NUMFILES = 29
NUMROWS = NUMFILES*NUMLINES
A = numpy.zeros(shape=(NUMROWS,5)) # This is creating a matrix with NUMROWS rows and 5 columns. The columns will be s,p,<x>, <T> and the checkvalue
Dx= numpy.zeros(shape=(4,NUMFILES)) #This creates a matrix with p as first row and s as second row, and x as third row,the 4th row is the value of s such that x=0.5 (is the one computed with the linear interpolation). This is the matrix with the first value of s such that x<=0.5
xtoprint=range(NUMFILES-1)
ytoprint=range(NUMFILES-1)

#********************************Let's define the function I am going to use to to check for which value of s I have exactly 0.5
def computeexactx(matrix):
	#print " "
	dum1=numpy.where(matrix[1]<=0.5)
	dum2=numpy.where(matrix[1]>=0.5)
	DUM1=numpy.amin(dum1)
	DUM2=numpy.amax(dum2)
	s1=matrix[0][DUM1]
	s2=matrix[0][DUM2]
	x1=matrix[1][DUM1]
	x2=matrix[1][DUM2]
	mdummy=(s2-s1)/(x2-x1)
	s0dummy=(s1*x2-s2*x1)/(x2-x1)
	s0five=mdummy*0.5+s0dummy
	#print s0five, s1, s2, x1, x2, DUM1, dum1[0][0]
	return s0five


#**************** NOW LET'S FILL THE MATRIX ***************************
# Here I will basically just attach the files one after the other adding as second column the value of p!
l=0
for i in range(NUMFILES):
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
	dummy1=b.transpose() # In this line and in the next for I look for the first value of s such that x is smaller equal than 0.52 (I chose 0.52 because in this way I'm getting values for x that are closer to 0.5)
	dummy2=numpy.where(dummy1[1]<=0.52)
	dummy3=computeexactx(dummy1)
	Dx[0][i]=p
	Dx[1][i]=dummy1[0][dummy2[0][0]]
	Dx[2][i]=dummy1[1][dummy2[0][0]]
	Dx[3][i]=dummy3
	

#Dx[3][0]=0.5 #Here I manually correct for p=0,s=0
#**************************************************************Let's plot it ************************************

scatter(Dx[1],Dx[0],c=Dx[2],cmap=cm.jet)
plt.colorbar()
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_normal"+".png",dpi=200)
exponent1=0.9305
a1=exp(-1.4026)
exponent2=4.3559
a2=exp(4.8443)
x1=arange(0,0.15,0.01)
x2=arange(0.15,0.56,0.01)
y1=a1*(x1**exponent1)
y2=a2*(x2**exponent2)
plot(x1,y1)
plot(x2,y2)
#plt.show()
close()

#Let's make a log plot instead:
yscale('log')
scatter(Dx[1],Dx[0],c=Dx[2],cmap=cm.jet)
plt.colorbar()
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_log"+".png",dpi=200)
#plt.show()
close()


#Let's make a log log plot
yscale('log')
xscale('log')
scatter(Dx[1],Dx[0],c=Dx[2],cmap=cm.jet)
plt.colorbar()
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_log_log"+".png",dpi=200)
#plt.show()
close()

#*******************************************************************************************

#**********************NOW LET'S PLOT THE ONE WITH THE INTERPOLATED DATA************************
scatter(Dx[3],Dx[0])
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_normal_interpolated"+".png",dpi=200)
exponent1=0.9305
a1=exp(-1.4026)
exponent2=4.3559
a2=exp(4.8443)
x1=arange(0,0.15,0.01)
x2=arange(0.15,0.56,0.01)
y1=a1*(x1**exponent1)
y2=a2*(x2**exponent2)
plot(x1,y1)
plot(x2,y2)
#plt.show()
close()

#Let's make a log plot instead:
yscale('log')
scatter(Dx[3],Dx[0])
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_log_interpolated"+".png",dpi=200)
close()

#Let's make a log log plot
yscale('log')
xscale('log')
scatter(Dx[3],Dx[0])
title("Transition line for the random splitting model")
xlabel("s",fontsize=12)
ylabel("p",fontsize=12)
plt.savefig("transition_line_log_log_interpolated"+".png",dpi=200)
exponent1=0.9305
a1=exp(-1.4026)
exponent2=4.3559
a2=exp(4.8443)
x1=arange(0,0.15,0.01)
x2=arange(0.15,0.56,0.01)
y1=a1*(x1**exponent1)
y2=a2*(x2**exponent2)
plot(x1,y1)
plot(x2,y2)
#plt.show()
close()

#***************Now let's save the one for the thesis!*****************************
axes([0,0,1,1])
figure(num=None, figsize=(12, 9), dpi=160, facecolor='w', edgecolor='k')
title("Transition line for the random splitting model",fontsize=20)
scatter(Dx[3],Dx[0])
xlabel("s",fontsize=20)
ylabel("p",fontsize=20)
exponent1=0.9305
a1=exp(-1.4026)
exponent2=4.3559
a2=exp(4.8443)
x1=arange(0,0.15,0.01)
x2=arange(0.15,0.56,0.01)
y1=a1*(x1**exponent1)
y2=a2*(x2**exponent2)
plot(x1,y1)
plot(x2,y2)
axes([0.2,0.43,0.42,0.42])
yscale('log')
xscale('log')
scatter(Dx[3],Dx[0])
plt.tight_layout()
plt.savefig("Double_pic"+".png",dpi=200)
close()

#*************Now let's prepare the data for printing**********************
# I want to print them in log-log scale

for i in range(1,NUMFILES): #I don't print the 0,0 value because in log log scale is going to be -infinity
	xtoprint[i-1]=log(Dx[3][i])
	ytoprint[i-1]=log(Dx[0][i])



print xtoprint
print ytoprint



