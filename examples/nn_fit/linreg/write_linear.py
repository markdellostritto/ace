#*************************************
# Import Statements
#*************************************
import math
import random as rand

#*************************************
# Global Variables
#*************************************
m=1.0
b=0.0

#*************************************
# Gaussian Function
#*************************************
def linear(x):
    return x*m+b+rand.uniform(-1.0,1.0)

xmin=-3.0
xmax=3.0
dx=0.05
nx=int((xmax-xmin)/dx)

# write tranining data
xmin=-3.0
xmax=3.0
f=open("linear_train.dat","w")
f.write("#X Y\n")
for i in range(0,nx):
    x=xmin+i*dx
    y=linear(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

# write validation data
xmin=-3.0
xmax=3.0
f=open("linear_val.dat","w")
f.write("#X Y\n")
for i in range(0,int(nx/5)):
    x=xmin+(xmax-xmin)*rand.random()
    y=linear(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

# write test data
xmin=-3.0
xmax=3.0
f=open("linear_test.dat","w")
f.write("#X Y\n")
for i in range(0,nx):
    x=xmin+(xmax-xmin)*rand.random()
    y=linear(x)
    f.write(str(x)+" "+str(y)+"\n")
f.close()

