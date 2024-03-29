#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 17:00:54 2016

@author: Jonathan Astermark, Anders Hansson, Tuong Lam and Oskar Smedman

**********************************
    Detects when event occurs
**********************************

"""
import numpy as N

def state_event(t,y,yd,sw):
    '''
    Calculates 4 values from 4 functions. Each function represents an event and changes sign when
    an event occurs.
    '''
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

    phiS=y[1]
    lambda1=y[6]
    phiB=y[2]

    # Crossing function 1, changes sign when transition from State 1 to State 2 occurs
    R1=hS*phiS+(rS-r0)

    # Crossing function 2, changes sign when transition from State 1 to 3 occurs
    R2=hS*phiS-(rS-r0)

    # Crossing function 3, changes sign when transition from State 2 to 1 occurs
    R3=lambda1

    # Crossing function 4, changes sign when transition from State 3 to 4 and back occurs
    R4=hB*phiB-lS-lG+lB+r0

    return N.array([R1,R2,R3,R4])
