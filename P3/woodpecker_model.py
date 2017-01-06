# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 14:15:39 2016
@author: Jonathan
"""
import numpy as N

def init_woodpecker():
    q = N.zeros((3,))
    qd = N.zeros((3,))
    qdd = N.zeros((3,))
    lam = N.zeros((2,))
    lamd = N.zeros((2,))
    
#    qdd[0] = -9.81
    qd[1] = -1
    y = N.hstack((q,qd,lam))
    yd = N.hstack((qd,qdd,lamd))
    
    sw0 = [True,False,False]
    
    return y,yd,sw0


def res(t,y,yd,sw):
#    sw = N.array([False,True,False])
#    print(sw)
    # Variable alisases
    q = y[0:3]
    qd = y[3:6]
    lam = y[6:]
    
    v = yd[0:3]
    w = yd[3:6]
    lamd = yd[6:]
    
    # Position coordinates
    z = q[0]
    phiS = q[1]
    phiB = q[2]
    
    # Mechanical and geometric constants of the woodpecker model
    # according to 
    # Ch. Glocker: Dynamik von Starrkörpersystemen mit Reibung und Stössen
    # PhD thesis TU München, July 1995, pp. 162ff
    #
    mS = 3.0e-4 # Mass of sleeve [kg]
    JS = 5.0e-9 # Moment of inertia of the sleeve [kgm]
    mB = 4.5e-3 # Mass of bird [kg]
    masstotal=mS+mB # total mass
    JB = 7.0e-7 # Moment of inertia of bird [kgm]
    r0 = 2.5e-3 # Radius of the bar [m]
    rS = 3.1e-3 # Inner Radius of sleeve [m]
    hS = 5.8e-3 # 1/2 height of sleeve [m]
    lS = 1.0e-2 # verical distance sleeve origin to spring origin [m]
    lG = 1.5e-2 # vertical distance spring origin to bird origin [m]
    hB = 2.0e-2 # y coordinate beak (in bird coordinate system) [m]
    lB = 2.01e-2 # -x coordinate beak (in bird coordinate system) [m]
    cp = 5.6e-3 # rotational spring constant [N/rad]
    g  = 9.81 #  [m/s^2]
    
    # Mass matrix
    M = N.zeros((3,3))
    M[0,0] = mB+mS
    M[0,1] = mB*lS
    M[0,2] = mB*lG
    M[1,0] = M[0,1]
    M[1,1] = JB + mB*(lS**2)
    M[1,2] = mB*lS*lG
    M[2,0] = M[0,2]
    M[2,1] = M[1,2]
    M[2,2] = JB + mB*(lG**2)
    
    # h
    h = N.zeros((3,))
    h[0] = -(mB+mS)*g
    h[1] = -cp*(phiS-phiB) - mB*g*lS
    h[2] = -cp*(phiB-phiS) - mB*g*lG
    
    # Force directions wN and wT
    wN = N.zeros((3,))
    wT = N.zeros((3,))
    
    # In the following it is assumed that
    # len(sw)=3 and only one state is activated at a time
    # (handle_event is responsible for this)
    # TODO if time allows: implement a check for this
    
    # State 1
    if sw[0]:
#        wN[2] = -hB
#        wT[0] = 1
#        wT[1] = lS
#        wT[2] = lG-lB
        wN[1] = -1
        wT[2] = -1
        gc = lam[0]
        f = lam[1]

    # State 2
    if sw[1]:
        wN[1] = -hS
        wT[0] = -1
        wT[1] = -rS
        wT[2] = 0
        gc = (rS-r0) + hS*phiS # Index-3 constraint
        f = v[0] + rS*v[1]

    # State 3
    if sw[2]:
        wN[1] = hS
        wT[0] = -1
        wT[1] = -rS
        wT[2] = 0
        gc = (rS-r0) - hS*phiS # Index-3 constraint
        f = v[0] + rS*v[1]

    # State 4 - never occurs since handle_event immediately sends it back to 3
        
    # Calculate residuals
    res1 = v - qd
    res2 = N.dot(M,w) - h - wN*lam[0] - wT*lam[1]
    res3 = gc
    res4 = f
    
    return N.hstack((res1,res2,res3,res4)) # len(y)