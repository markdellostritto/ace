#******************************************************
# Import Statements
#******************************************************

import sys
import getopt

#******************************************************
# Global Variables
#******************************************************

ifile=''
ofile=''
lat=0

#******************************************************
# Arguments
#******************************************************

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:l:",["ifile=","ofile=","lat="])
except getopt.GetoptError:
    print ('vdos.py -i <ifile> -o <ofile> -t <ts> -w <width>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('test.py -i <ifile> -o <ofile> -t <ts> -w <width>')
        sys.exit()
    elif opt in ("-i", "--ifile"):
        ifile = arg
    elif opt in ("-o", "--ofile"):
        ofile = arg
    elif opt in ("-l","--lat"):
        lat = float(arg)

print("ifile = ",ifile)
print("ofile = ",ofile)
print("lat   = ",lat)

#******************************************************
# Read
#******************************************************

# open xyz file
reader=open(ifile,"r")
# read natoms
natoms=int(reader.readline())
#print('natoms = ',natoms)
# read header
reader.readline()
# read atoms
names=[]
posns=[]
for i in range(natoms):
    sarr=reader.readline().split()
    names.append(sarr[0])
    posns.append([float(sarr[1]),float(sarr[2]),float(sarr[3])])
# set mass
mass=[]
for i in range(natoms):
    if(names[i]=="Cl"): mass.append(35.45)
    if(names[i]=="Na"): mass.append(22.99)
types=[]
for i in range(natoms):
    if(names[i]=="Cl"): types.append(1)
    if(names[i]=="Na"): types.append(2)
chg=[]
for i in range(natoms):
    if(names[i]=="Cl"): chg.append(-1.0)
    if(names[i]=="Na"): chg.append(1.0)
#for i in range(natoms): print(names[i],mass[i],posns[i])
# close reader
reader.close()

#******************************************************
# Write
#******************************************************

# open data file
writer=open(ofile,"w")
# write header
writer.write("#header\n\n")
# write atom data
writer.write(" %i atoms\n"%(natoms))
writer.write(" 2 atom types\n")
writer.write("\n")
# write box
writer.write(" 0.0 %.1f xlo xhi\n"%(lat))
writer.write(" 0.0 %.1f ylo yhi\n"%(lat))
writer.write(" 0.0 %.1f zlo zhi\n"%(lat))
writer.write("\n")
# write masses
writer.write(" Masses\n\n")
writer.write("  1 35.45\n")
writer.write("  2 22.99\n")
writer.write("\n")
# write atoms
writer.write(" Atoms # charge\n\n")
for i in range(natoms):
    writer.write("%i %i %3.1f %19.10f %19.10f %19.10f\n"%(i+1,types[i],chg[i],posns[i][0],posns[i][1],posns[i][2]))
# close writer
writer.close()

