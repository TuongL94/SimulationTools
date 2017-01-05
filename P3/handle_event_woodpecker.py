#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 18:14:23 2016

@author: anders
"""
import numpy as N


def handle_event(solver,event_info):
    """
    Event handling. This functions is called when Assimulo finds an event as
    specified by the event functions.
    """
    state_info=event_info[0] #We are only interested in state events
    '''
    Event 1 occurs
    '''
    if state_info[0] !=0:
        if solver.sw[0] ==True: #Currently in State 1
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]<0: # Bird rotates in the right direction
                change_momentum()
                solver.sw=N.array([False,True,False]) #Transistion from state 1 to 2
    '''
    Event 2 occurs
    '''
    if state_info[1] !=0:
        if solver.sw[0] ==True: #Currently in state 1
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]>0: #Bird rotates in the right direction
                change_momentum()
                solver.sw=N.array([False,False,True]) #Transition from state 1 to 3
        
    '''
    Event 3 occurs
    '''
    if state_info[2] !=0:
        if solver.sw[1] ==True: #Currently in state 2
            #raise StateError('Event occured when it was not allowed')
            solver.sw=N.array([True,False,False]) #Transition to state 1
        
    '''
    Event 4 occurs
    '''
    if state_info[3] !=0:
        if solver.sw[2] ==True: #Currently in state 3
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]<0: #Bird rotates in the right direction
                solver.sw=N.array([True,False,False]) #Transition to state 1
    
    '''
    Event 5 occurs
    '''
    if state_info[4] !=0:
        if solver.sw[2] ==True: #Currently in state 3
            #raise StateError('Event occured when it was not allowed')
            if solver.y[5]>0: #Bird rotates in the right direction
                solver.y[5]=-solver.y[5] # Bird hits the bar, changes rotation direction without losses
                          
    '''
    Event 6 occurs, i.e. the woodpecker hits the bar. Don't need this?
    '''
    '''
    if state_info[5] !=0: #Woodpecker hits the bar
        if solver.sw[3] !=True: # Must be in state 4 in order for this to happen
            raise StateError('Event occured when it was not allowed')
        solver.yd[2]=-solver.yd[2] #Change sign of bird rotation
        solver.sw=N.array([False,False,True,False]) #Transition to state 3
    '''
        
                          

def change_momentum():
    '''
    Conserves the momentum
    '''
    I_minus=mB*lG*y[3]+(mB*lG*lS)*y[4]+(JB+mB*lG**2)*y[5]
    solver.y[5]=I_minus/(JB+mB*lG**2)
    
    
class StateError(Exception):
    pass