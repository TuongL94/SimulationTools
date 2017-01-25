#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 18:14:23 2016

@author: anders

*****************************

    Handles events

*****************************
"""
import numpy as N


def handle_event(solver,event_info):
    """
    Event handling. This functions is called when Assimulo finds an event as
    specified by the event functions.
    """
    state_info=event_info[0] #We are only interested in state events
#    print(event_info)
    '''
    Event 1 occurs - transition from state 1 to 2
    '''
    if state_info[0] !=0:
        if solver.sw[0] ==True: #Currently in State 1
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]<0: # Bird rotates in the right direction
                change_momentum(solver)
                solver.y[3] = 0
                solver.y[4] = 0
                solver.sw=[False,True,False] #Transistion from state 1 to 2
    '''
    Event 2 occurs - transition from state 1 to 3
    '''
    if state_info[1] !=0:
        if solver.sw[0] ==True: #Currently in state 1
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]>0: #Bird rotates in the right direction
                change_momentum(solver)
                solver.y[3] = 0
                solver.y[4] = 0
                solver.sw=[False,False,True] #Transition from state 1 to 3

    '''
    Event 3 occurs - transition from state 2 to 1
    '''
    if state_info[2] !=0:
        if solver.sw[1] ==True: #Currently in state 2
            #raise StateError('Event occured when it was not allowed')
            solver.sw=[True,False,False] #Transition to state 1

    '''
    Event 4 occurs - transition from state 3 to 1
    '''
    if state_info[3] !=0:
        if solver.sw[2] ==True: #Currently in state 3
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]<0: #Bird rotates in the right direction
                solver.sw=[True,False,False] #Transition to state 1

    '''
    Event 5 occurs - beak hits the bar
    '''
    if state_info[3] !=0:
        if solver.sw[2] ==True: #Currently in state 3
            #raise StateError('Event occured when it was not allowed')
            #if solver.y[5]>0: #Bird rotates in the right direction
                solver.y[5]=-solver.y[5] # Bird hits the bar, changes rotation direction without losses

def change_momentum(solver):
    '''
    Conserves the momentum
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

    I_minus=mB*lG*solver.y[3]+(mB*lG*lS)*solver.y[4]+(JB+mB*lG**2)*solver.y[5]
    solver.y[5]=I_minus/(JB+mB*lG**2)

class StateError(Exception):
    pass
