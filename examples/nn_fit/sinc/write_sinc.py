#*************************************
# Import Statements
#*************************************
import math
import random as rand

#*************************************
# Global Variables
#*************************************

#*************************************
# Gaussian Function
#*************************************
def sinc(x):
    return x*math.sin(x)

xmin=0.0
xmax=10.0
dx=0.05
nx=int((xmax-xmin)/dx)

# write tranining data
xmin=0.0
xmax=10.0
f=open("sinc_train.dat","w")
f.write("#X Y\n")
for i in range(0,nx):
    x=xmin+i*dx
    y=sinc(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

# write validation data
xmin=0.0
xmax=10.0
f=open("sinc_val.dat","w")
f.write("#X Y\n")
for i in range(0,int(nx/5)):
    x=xmin+(xmax-xmin)*rand.random()
    y=sinc(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

# write test data
xmin=0.0
xmax=10.0
f=open("sinc_test.dat","w")
f.write("#X Y\n")
for i in range(0,nx):
    x=xmin+(xmax-xmin)*rand.random()
    y=sinc(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

