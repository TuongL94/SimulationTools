#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 17:00:54 2016

@author: anders
"""
from woodpecker_constants import*


def state_event(t,y,sw):
    '''
    Calculates 5 values from 5 functions. Each function represents an event and changes sign when 
    an event occurs.
    '''
    
    #Function 1, transition from State 1 to State 2
    phiS=y[1]
    R1=hS*phiS+(rS-r0)
    
    #FUnction 2, transition from State 1 to 3
    phiS=y[1]
    R2=hS*phiS-(rS-r0)
    
    #Function 3, transition from State 2 to 1
    lambda1=y[6]
    R3=lambda1
    
    #Function 4, transition from State 3 to 1
    R4=lambda1
    
    #Function 5, transition from State 3 to 4 and back
    phiB=y[2]
    R5=hB*phiB-lS-lG+lB+r0
    
    return N.array([R1,R2,R3,R4,R5])
    