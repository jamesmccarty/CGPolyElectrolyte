#!/usr/bin/env python

''' script to generate LAMMPS input for running bead spring polyelectrolyte chains '''

import math
import numpy as np
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.special import erfc
from newsavetxt import savetxthd
from FTSparams import *

M=20
Npoints =2**M
gewald = 0.467289719626

def gauss(r,a,norm):
    ''' Excluded volume soft Gaussian repulsive potential '''
    return norm*np.exp(-r**2/4.0/a**2)

def gaussforce(r,a,norm):
    return (r*norm*np.exp(-r**2/4.0/a**2))/(2*a**2)

def coul(r,lb,a,gewald,norm):
    ''' Smeared Coulomb potential '''
    #print(lb)
    vel=(lb/r)*(erf(r/(2.0*a))-erf(r*gewald))
    excl = gauss(r,a,norm)
    return vel + excl

def eforce(r,lb,a,gewald,norm):
    force=lb*(erf(r/(2.0*a))-erf(r*gewald))/r**2-lb*(np.exp(-r**2/(4.0*a**2))/(a*math.sqrt(math.pi))-2.0*gewald*np.exp(-r**2*gewald**2)/(math.sqrt(math.pi)))/r
    exclf=gaussforce(r,a,norm)
    return force + exclf

def coul12(r,lb,a,gewald,norm):
    ''' Smeared Coulomb potential '''
    #print(lb)
    vel=(lb/r)*(erf(r/(2.0*a))-erf(r*gewald))
    excl = gauss(r,a,norm)
    return -vel + excl

def eforce12(r,lb,a,gewald,norm):
    force=lb*(erf(r/(2.0*a))-erf(r*gewald))/r**2-lb*(np.exp(-r**2/(4.0*a**2))/(a*math.sqrt(math.pi))-2.0*gewald*np.exp(-r**2*gewald**2)/(math.sqrt(math.pi)))/r
    exclf=gaussforce(r,a,norm)
    return -force + exclf

def get_params(vex,lB,atilde,Nref,rho):
    ''' convert between FTS and CG parameters '''
  #  B=0.1
  #  E=10.0
  #  abar=1.0
  #  Nref=1.0

    ### convert to FTS parameters:
    Rg=math.sqrt(Nref/6.0)
    B=vex*Nref*Nref/(Rg*Rg*Rg)   # B=vex*N^2/Rg^3
    abar = atilde/Rg # a/Rg
    E = (4.0*math.pi)*lB*Nref*Nref/Rg # reduced Bjerrum length E=4 pi lB N^2/Rg
    C = rho*Rg*Rg*Rg/Nref

    print('FTS parameters')
    print('excluded volume: B', B)
    print('smearing length: abar', abar)
    print('reduced Bjerrum length: E', E)
    print('chain density: C', C)

    if(scale_Rg):
        norm = B/(8.0*math.pi**(3.0/2.0)*abar**3)
    else:
        norm = vex/(8.0*math.pi**(3.0/2.0)*atilde**3) # prefactor of excluded volume potential
    return (norm)

fout = open('in.soft','w')

norm = get_params(vex,lB,atilde,Nref,rho)

bond_const = 1.5

Rg=math.sqrt(Nref/6.0)
if(scale_Rg):
    print("scaling distances by Rg: ",Rg)
    bond_const = 1.5*Rg*Rg
    print("bond K: ",bond_const)
    atilde = atilde/Rg
    print("smearing length: ", atilde)
    lB = lB/Rg
    print("Bjerrum length: ", lB)
    print("Set reduced density to ",rho*Rg*Rg*Rg/Nref)
else:
    print("scaling distances by b")
    print("bond K: ",bond_const)
    print("smearing length: ", atilde)
    print("Bjerrum length: ",lB)
    print("set reduced density to ",rho)

if not salts:
    NCA=0
    NCC=0


print('reduced Bjerrum length [lB/b]', lB)

dielectric = 1.0/lB   # effective dielectric constant for LAMMPS
print('dielectric', dielectric)
print('Gaussian potential normalization constant',norm)
rfinal = rcut+0.5
xdata1 = np.linspace(0.0,rfinal,Npoints+1)
# drop zero point to avoid LAMMPS error
xdata1 = np.delete(xdata1,0)
xdata = xdata1
index = np.arange(Npoints)
index = index+1
ydata = coul(xdata,lB,atilde,gewald,norm)
fdata = eforce(xdata,lB,atilde,gewald,norm)
dataout = np.column_stack((index,xdata,ydata,fdata))
hdrtxt='COULOMBSMEARAA\n'
hdrtxt += ' N       ' + str(Npoints)+'\n\n'
savetxthd('potentialfileAA.dat',dataout, header=hdrtxt, fmt=('%12.0i   %24.16e %24.16E %24.16E'))
hdrtxt='COULOMBSMEARCC\n'
hdrtxt += ' N       ' + str(Npoints)+'\n\n'
savetxthd('potentialfileCC.dat',dataout, header=hdrtxt, fmt=('%12.0i   %24.16e %24.16E %24.16E'))
hdrtxt='COULOMBSMEARAC\n'
hdrtxt += ' N       ' + str(Npoints)+'\n\n'
ydata = coul12(xdata,lB,atilde,gewald,norm)
fdata = eforce12(xdata,lB,atilde,gewald,norm)
dataout = np.column_stack((index,xdata,ydata,fdata))
savetxthd('potentialfileAC.dat',dataout, header=hdrtxt, fmt=('%12.0i   %24.16e %24.16E %24.16E'))

fout.write('atom_style full\n\n')
fout.write('boundary p p p\n\n')
fout.write('units lj \n\n')
fout.write('dielectric '+str(dielectric)+'\n\n')
fout.write('read_data startconf.dat\n\n')
if(salts):
    fout.write('create_atoms 3 random '+str(NCA)+' 23053 NULL # anion\n')
    fout.write('create_atoms 4 random '+str(NCC)+' 34201 NULL # cation\n\n')
    fout.write('set type 3 charge -1.0\n')
    fout.write('set type 4 charge 1.0\n\n')
fout.write('mass 1 1.0\n')  # mass set equal to 1
fout.write('mass 2 1.0\n\n')
if(salts):
    fout.write('mass 3 1.0\n')
    fout.write('mass 4 1.0\n\n')

fout.write('bond_style harmonic\n')

fout.write('bond_coeff 1 '+str(bond_const)+' 0.0\n\n')     # harmonic bond
fout.write('pair_style table bitmap '+str(M)+' '+kspace+'\n')
#fout.write('pair_style table spline '+str(Npoints-2)+' ewald\n')
fout.write('\n')
txt11 = 'pair_coeff 1 1 potentialfileAA.dat COULOMBSMEARAA '+str(rcut)+'\n'
fout.write(txt11)
txt22 = 'pair_coeff 2 2 potentialfileCC.dat COULOMBSMEARCC '+str(rcut)+'\n'
fout.write(txt22)
txt12 = 'pair_coeff 1 2 potentialfileAC.dat COULOMBSMEARAC '+str(rcut)+'\n'
fout.write(txt12)
if(salts):
    txt11 = 'pair_coeff 3 3 potentialfileAA.dat COULOMBSMEARAA '+str(rcut)+'\n'
    fout.write(txt11)
    txt22 = 'pair_coeff 4 4 potentialfileCC.dat COULOMBSMEARCC '+str(rcut)+'\n'
    fout.write(txt22)
    txt12 = 'pair_coeff 3 4 potentialfileAC.dat COULOMBSMEARAC '+str(rcut)+'\n'
    fout.write(txt12)
    txt11 = 'pair_coeff 1 3 potentialfileAA.dat COULOMBSMEARAA '+str(rcut)+'\n'
    fout.write(txt11)
    txt22 = 'pair_coeff 2 4 potentialfileCC.dat COULOMBSMEARCC '+str(rcut)+'\n'
    fout.write(txt22)
    txt12 = 'pair_coeff 1 4 potentialfileAC.dat COULOMBSMEARAC '+str(rcut)+'\n'
    fout.write(txt12)
    txt12 = 'pair_coeff 2 3 potentialfileAC.dat COULOMBSMEARAC '+str(rcut)+'\n'
    fout.write(txt12)
fout.write('\n')
fout.write('kspace_style '+kspace+' 1e-4\n\n')
fout.write('kspace_modify gewald '+str(gewald)+'\n')
#fout.write('velocity all create 1.0 59683\n\n')
fout.write('pair_write 1 1 1000 r 0.002 '+str(rcut-0.1)+' pair11.txt Vex -1.0 -1.0\n')
fout.write('pair_write 2 2 1000 r 0.002 '+str(rcut-0.1)+' pair22.txt Vex 1.0 1.0\n')
fout.write('pair_write 1 2 1000 r 0.002 '+str(rcut-0.1)+' pair12.txt Vex -1.0 1.0\n')
fout.write('\n')
fout.write('special_bonds lj 1.0 1.0 1.0\n\n')
fout.write('neighbor '+str(neighbin)+' bin\n')
fout.write('neigh_modify every 1 delay 0 check yes\n')
fout.write('neigh_modify one 100000\n')
fout.write('neigh_modify page 1000000\n\n')
fout.write('fix 2 all nve\n')
fout.write('fix 3 all temp/csld 1.0 1.0 0.1 54324\n\n')   # langevin thermostat
fout.write('thermo_style custom step temp ke pe etotal press pxx pyy pzz pxy pxz pyz elong\n\n')
fout.write('thermo_modify lost warn\n\n')
fout.write('thermo_modify flush yes\n')
fout.write('thermo 1000\n\n')
fout.write('restart 1000 restart1.dat restart2.dat\n\n')
fout.write('dump dump1 all custom 5000 coords.lammpstrj id mol type x y z ix iy iz\n')
fout.write('dump_modify dump1 sort id\n\n')
fout.write('timestep '+str(dt)+'\n\n')
# uncomment next line for plumed
#fout.write('fix 22 all plumed plumedfile plumed.dat\n\n')
fout.write('run 10000000\n\n')
fout.write('write_data outputconfig.dat\n')
fout.write('write_restart output.restart\n')
