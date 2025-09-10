#!/usr/bin/env python
# coding: utf-8

# In[1]:


#### Thsi is python programme to generate input file for CADSSS ####
####              Framework must be in lammps format              ####
####       Written By Supriyo Naskar, CNRS Montpellier          ####


# In[2]:


from __future__ import division
import math 
import matplotlib.pyplot as plt
import sys
import sympy as sp
import numpy as np
import cmath
import mdtraj as md
import MDAnalysis as mda
import os
import shutil



# In[3]:
os.popen('cp ../final-interface.lmps final.data')


name1    = "Systema"
exten    = ".car"
dirname1 = "/home/adiaz/programs/CADSS/share/"
dirname2 = "/home/adiaz/programs/CADSS/share/structures/car/"
dirname3 = "/home/adiaz/programs/CADSS/share/forcefield/"
dirname4 = dirname3 + name1
dirname5 = dirname3 + name1 + "/"
name2    = dirname2 + name1 + exten
name3    = dirname4 + "/pseudo_atoms.def"
name4    = dirname4 + "/force_field_mixing_rules.def"

#s = "os.popen('cp /home/snaskar/CADSS/V3/share/forcefield/calf/force_field.def 'dirname4' ')"




#Check the directory name exist or not
if os.path.isdir(dirname4) == False:
    #Create the directory
    os.mkdir(dirname4)
    #Print success message
    print("The directory is created.")
else:
    #Print the message if the directory exists
    print("The directory already exists.")


shutil.copy2("/home/adiaz/programs/CADSS/share/forcefield/force_field.def",dirname5)



print("#### Thsi is python programme to generate input file for CADSSS ####")
print("####              Framework must be in lammps format              ####")
print("####       Written By Supriyo Naskar, CNRS Montpellier          ####")



u = mda.Universe("final.data", atom_style="id resid type charge x y z")
pos=u.atoms.positions
charge=u.atoms.charges
resid=u.atoms.resids
mass=u.atoms.masses
no_atoms=len(pos)
no_res=len(resid)
box=u.dimensions   ##xx,yy,zz,alpha,beta,gamma
type_name=u.atoms


# In[4]:


### since lammpsdata file doesnot have type, and covalent radius, you need to put them by hand ###
### Remeber to put the order same as lammps data file

n1=len(set(u.atoms.types))
found_type = False
Mass1 = []


with open('final.data', 'r') as f:
    for line in f:
        if 'Masses' in line:
            found_type = True
            k=0
            continue 
        if found_type:
            if k==0:    
                k= k+1
            elif k==n1+1:
                found_type = False                
            else:    
                t_line = str(line).rstrip('\n')
                t_line=t_line.split()
                t_line=list(np.float_(t_line))
                #a=np.array(t_line)
                if t_line!='':
                    Mass1.append(t_line)
                k=k+1
Mass = np.empty(len(Mass1))
for i in range (len(Mass1)):
        Mass[i]= float(Mass1[i][1])
atom_name = []
atom_type = []
k1=0
k2=0
k3=0
k4=0
k5=0
k6=0
k7=0
k8=0
k9=0
k10=0
k11=0
for i in range (len(Mass)):
    if abs(Mass[i]-12.01) < 0.2:
        k1 += 1
        atom_type.append("C"+str(k1))
        atom_name.append("C")
    if abs(Mass[i]-1.008) < 0.2:          
        k2 = k2 + 1
        atom_type.append("H"+str(k2))
        atom_name.append("H")
    if abs(Mass[i]-16) < 0.2:          
        k3 = k3 + 1
        atom_type.append("O"+str(k3))
        atom_name.append("O")
    if abs(Mass[i]-14.01) < 0.2:          
        k4 = k4 + 1
        atom_type.append("N"+str(k4))
        atom_name.append("N")
    if abs(Mass[i]-28.0855) < 0.2:          
        k5 = k5 + 1
        atom_type.append("Si"+str(k5))
        atom_name.append("Si")
    if abs(Mass[i]-19) < 0.2:          
        k6 = k6 + 1
        atom_type.append("F"+str(k6))
        atom_name.append("F")
    if abs(Mass[i]-65.38) < 0.2:          
        k7 = k7 + 1
        atom_type.append("Zn"+str(k7))
        atom_name.append("Zn")
    if abs(Mass[i]-26.9815) < 0.2:          
        k8 = k8 + 1
        atom_type.append("Al"+str(k8))
        atom_name.append("Al")
    if abs(Mass[i]-58.6934) < 0.2:          
        k9 = k9 + 1
        atom_type.append("Ni"+str(k9))
        atom_name.append("Ni")
    if abs(Mass[i]-91.224) < 0.2:          
        k10 = k10 + 1
        atom_type.append("Zr"+str(k9))
        atom_name.append("Zr")
    if abs(Mass[i]-32.06) < 0.2:          
        k11 = k11 + 1
        atom_type.append("S"+str(k9))
        atom_name.append("S")
                

n2=len (atom_type)
if n1 != n2: 
    print ("*****Warning*******")
    print ("total no of atomtypes not matching with lammps data file... Something is wrong")

print (n1,n2,atom_name,atom_type)
def co_rad (type):
    if type=="Zn":
        out1 =  1.21
    if type=="C":
        out1 =  0.77
    if type=="N":
        out1 =  0.74
    if type=="O":
        out1 =  0.74
    if type=="F":
        out1 =  0.57
    if type=="H":
        out1 =  0.37
    if type=="S":
        out1 =  0.80
    if type=="Si":
        out1 =  1.05
    if type=="Ni":
        out1 =  1.24
    if type=="Al":
        out1 =  1.21
    if type=="Zr":
        out1 =  1.75
    if type=="S":
        out1 = 1.02  
    return out1
    



# In[5]:


#### Read the pair potential /Mass values from lammpsdata as mdanalysis cant read it ####
found_type = False
P_coeff = []
with open('final.data', 'r') as f:
    for line in f:
        if 'Pair Coeffs' in line:
            found_type = True
            k=0
            continue 
        if found_type:
            if k==0:    
                k= k+1
            elif k==n1+1:
                found_type = False                
            else:    
                t_line = str(line).rstrip('\n')
                t_line=t_line.split()
                t_line=list(np.float_(t_line))
                #a=np.array(t_line)
                if t_line!='':
                    P_coeff.append(t_line)
                k=k+1


                
                


# In[6]:
import datetime
from datetime import date
import calendar
now = datetime.datetime.now()
my_date = date.today()
week = calendar.day_name[my_date.weekday()]
time = now.strftime("%H:%M:%S")
Date = now.strftime("%d")
Month= now.strftime("%B")
year = now.strftime("%Y")
time = now.strftime("%H:%M:%S")


#### Now first write the car file
f_car = open(name2,"w+")
f_car.write("!BIOSYM archive 3\nPBC=ON\nS Naskar python code Generated CAR File\n")
f_car.write("!Date %s %s %s %s %s\n"%(week, Month, Date, time, year))
f_car.write("PBC%9.3f%10.3f%11.4f%10.4f%10.4f%10.4f (P1)\n"%(box[0],box[1],box[2],box[3],box[4],box[5]))
for i in range (len(pos)):
    name_ref = atom_type[int(u.atoms.types[i])-1]
    type_ref = atom_name[int(u.atoms.types[i])-1]
    f_car.write("%-5s%15.9f%15.9f%15.9f%5s%2d      %-5s   %-2s%11.6f\n"%(name_ref,pos[i][0],pos[i][1],pos[i][2],'XXXX',1,name_ref,type_ref,charge[i]))
f_car.write("end\nend")
f_car.close()








#### Now Lets write atomic types file 
### select guest molecule type 1 for CO2, 2 for N2, 3 for CH4
guest_molecule_type = 1

f_at = open(name3,"w+")
f_at.write("#Number of pseudo atoms\n")
if guest_molecule_type == 1 :
    print ("Guest molecule is CO2")
    tot_at=len(atom_type)+2
    f_at.write("%d\n"%(tot_at))
if guest_molecule_type == 2 :
    print ("Guest molecule is N2")
    tot_at=len(atom_type)+2
    f_at.write("%d\n"%(tot_at))
if guest_molecule_type == 3 :
    print ("Guest molecule is CH4")
    tot_at=len(atom_type)+1
    f_at.write("%d\n"%(tot_at))
if guest_molecule_type == 4 :
    print ("Guest molecule is all")
    tot_at=len(atom_type)+5
    f_at.write("%d\n"%(tot_at))
f_at.write("#ID TypeName  output   element   type mass         charge(e)  beta covalent_radius(Angstrom)          \n")
for i in range(len(atom_type)):
    rad=co_rad(atom_name[i])
    #print (i+1,atom_type[i],atom_name[i],atom_name[i],Mass[i],rad)
    f_at.write("%-2d   %-5s      yes     %-2s  %-2s%15.6f     0.000   1.0%8.2f  0\n"%(i+1,atom_type[i],atom_name[i],atom_name[i],Mass[i],rad))
if guest_molecule_type == 1 :
    f_at.write("%-2d   CO2_C      yes     C   C       12.011000    0.6512   1.0    0.77  0\n"%(i+2))
    f_at.write("%-2d   CO2_O      yes     O   O       15.999000   -0.3256   1.0    0.74  0\n"%(i+3))    
if guest_molecule_type == 2 :
    f_at.write("%-2d   N2_N       yes     N   N       14.006740   -0.4820   1.0    0.74  0\n"%(i+2))
    f_at.write("%-2d   N2_COM     yes     X   X       00.000000    0.9640   1.0    0.70  0\n"%(i+3))  

if guest_molecule_type == 3 :
    f_at.write("%-2d   CH4_sp3    yes     C   C       16.043000     0.000   1.0    1.00  0\n"%(i+2))
    
if guest_molecule_type == 4 :    
    f_at.write("%-2d   CO2_C      yes     C   C       12.011000    0.6512   1.0    0.77  0\n"%(i+2))
    f_at.write("%-2d   CO2_O      yes     O   O       15.999000   -0.3256   1.0    0.74  0\n"%(i+3)) 
    f_at.write("%-2d   N2_N       yes     N   N       14.006740   -0.4820   1.0    0.74  0\n"%(i+4))
    f_at.write("%-2d   N2_COM     yes     X   X       00.000000    0.9640   1.0    0.70  0\n"%(i+5))  
    f_at.write("%-2d   CH4_sp3    yes     C   C       16.043000     0.000   1.0    1.00  0\n"%(i+6))    
f_at.write("\n\n\n\n\n\n\n#----------------------------------------------------------------------------------\nRadius used for determining the connectivities of all the framework atoms")
f_at.close()


# In[8]:


#### Now Lets write pair-potential file 

f_parm = open(name4,"w+")
f_parm.write("# Calculation method for VDW Interactions (options: shifted or truncated)\ntruncated\n")
f_parm.write("# Tail correction for VDW interactions (options: yes or no; only can be used for truncated method)\nyes\n")
f_parm.write("# Number of atomic types with potential parameters\n%d\n"%(tot_at))
f_parm.write("# Parameter Lists (Atomic_Type, Potential_Type, Epslon_[K], Sigma_[A]); additional Reduced Mass for FH potential (Half molecular mass)  \n")
for i in range(len(atom_type)):
    f_parm.write("%-2d    %-5s   lennard-jones%13.6f%11.3f\n"%(i+1,atom_type[i],P_coeff[i][1]/(0.001987204259),P_coeff[i][2]))
if guest_molecule_type == 1 :
    f_parm.write("%-2d    CO2_C   lennard-jones    28.129         2.757000    \n"%(i+2))
    f_parm.write("%-2d    CO2_O   lennard-jones    80.507         3.033000    \n"%(i+3))
if guest_molecule_type == 2 :                                    
    f_parm.write("%-2d    N2_N    lennard-jones    36.000         3.310000    \n"%(i+2))
    f_parm.write("%-2d    N2_COM  lennard-jones    00.000         0.000000    \n"%(i+3))
if guest_molecule_type == 3 :                                    
    f_parm.write("%-2d    CH4_sp3 lennard-jones   148.000         3.730000    \n"%(i+3))
                                                                 
if guest_molecule_type == 4 :                                    
    f_parm.write("%-2d    CO2_C   lennard-jones    28.129         2.757000    \n"%(i+2))
    f_parm.write("%-2d    CO2_O   lennard-jones    80.507         3.033000    \n"%(i+3))
    f_parm.write("%-2d    N2_N    lennard-jones    36.000         3.310000    \n"%(i+4))
    f_parm.write("%-2d    N2_COM  lennard-jones    00.000         0.000000    \n"%(i+5))
    f_parm.write("%-2d    CH4_sp3 lennard-jones   148.000         3.730000    \n"%(i+6))        

f_parm.write("# Mixing rules for LJ-Type Interactions (Options: Lorentz-Berthelot or Jorgensen)\nLorentz-Berthelot")
f_parm.close()


print ("Hola!!!! Enjoy")




os.remove("final.data")
