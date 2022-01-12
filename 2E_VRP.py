# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 12:02:34 2021

@author: stijn
"""

import numpy as np
from gurobipy import *
import matplotlib.pyplot as plt
import pandas as pd
import math

#%% ----- Problem -----

model = Model('Two-Echelon Vehicle Routing Problem')

#%% ----- Data -----

DEPOT_ID = [0]
DEPOT_xc = [1]
DEPOT_yc = [1]
DEPOT_S_TIME = [0]
DEPOT_DEMAND = [0]

SATELLITE_ID = [1,2]
SATELLITE_xc = [5,4]
SATELLITE_yc = [2,6]
SATELLITE_S_TIME = [50,50]
SATELLITE_DEMAND = [0,0]

CUSTOMER_ID = [3,4,5,6,7,8]
CUSTOMER_xc = [7,7,9,7,7,8]
CUSTOMER_yc = [9,6,8,8,1,3]
CUSTOMER_S_TIME = [10,10,10,10,10,10]
CUSTOMER_DEMAND = [1,1,1,1,1,1]

LOC_ID = DEPOT_ID + SATELLITE_ID + CUSTOMER_ID
xc = DEPOT_xc + SATELLITE_xc + CUSTOMER_xc
yc = DEPOT_yc + SATELLITE_yc + CUSTOMER_yc
S_TIME = DEPOT_S_TIME + SATELLITE_S_TIME + CUSTOMER_S_TIME
DEMAND = DEPOT_DEMAND + SATELLITE_DEMAND + CUSTOMER_DEMAND

d = len(DEPOT_ID)            #number of depots
h = len(SATELLITE_ID)        #number of satellites
c = len(CUSTOMER_ID)         #number of customers
Q_v = 10                      #capacity vessel
Q_cb = 2                       #capacity cargo bike
MNOVA = 1                    #maximum number of vessels available
MNOCBA = 3                   #maximum number of cargo bikes available
M = 1000000                                    

#%% ----- Sets and indices -----

D = [i for i in range (d)]                                  #set of depots
H = [i for i in range (d, d+h)]                             #set of satellites
C = [i for i in range (d+h, d+h+c)]                         #set of customers
N1 = D + H                                                  #set of nodes for 1st echelon
N2 = H + C                                                  #set of nodes for 2nd echelon
N = D + H + C                                               #set of all nodes
A = [(i, j) for i in N for j in N if i != j]                #set of arcs

K1 = [i for i in range(MNOVA)]                              #set of vessels
K2 = [i for i in range(MNOVA, MNOVA + MNOCBA)]              #set of cargo bikes
K = K1 + K2                                                 #set of vehicles

#%% ----- Parameters -----

c = {(i, j): np.hypot(xc[i]-xc[j],yc[i]-yc[j]) for i, j in A}         #euclidean distance
s = S_TIME                                                            #service time
q = {i: DEMAND[i] for i in N}                                         #demand of the customer

#%% ----- Decision variables -----

# R(i,j) - The route of the truck from customer i to j, 0 = don’t visit, 1 = visit
x = {}
for i in N:
    for j in N:
        for k in K:
            x[i,j,k] = model.addVar (lb = 0, vtype = GRB.BINARY)
        
# T(i) - Time counter of elapsed time at the arrival of customer i	        
T = {}
for i in N:
    T[i] = model.addVar (lb = 0, vtype = GRB.CONTINUOUS)
    
# L(i) - Load counter    
L = {}
for i in N:
    L[i] = model.addVar (lb = 0, vtype = GRB.CONTINUOUS)

#%% ----- Objective function (minimize total distance) -----

Total_distance_travelled = quicksum (c[i,j]*x[i,j,k] for i in N for j in N if i != j for k in K)

model.setObjective (Total_distance_travelled)
model.modelSense = GRB.MINIMIZE
model.update ()

#%% ----- Constraints -----

#first echelon  

#M times number of incoming arcs at hub must be greater or equal to outgoing arcs at hub
con1_1 = {}
for j in H:
    con1_1[i,j,k] = model.addConstr((M * quicksum(x[i,j,k] for i in N1 for k in K1)) >= quicksum(x[j,i,k] for i in N2 for k in K2))

con3_1 = {}
for i in N1:
    for k in K1:
        con3_1[i,i,k] = model.addConstr(x[i,i,k] == 0)        
    
con4_1 = {}
for j in N1:
    for k in K1:
        con4_1[i,j,k] = model.addConstr(quicksum(x[i,j,k] for i in N1) == quicksum(x[j,i,k] for i in N1))
        
con5_1 = {}
for k in K1:
    con5_1[k] = model.addConstr(quicksum(x[h,j,k] for j in H for h in D) <= 1)

    
con6_1 = {}
for i in N1:
    for j in N1:
        if j >= d:
            if i != j:
                for k in K1:
                    con6_1[i,j,k] = model.addConstr(T[i] + c[i,j]*x[i,j,k] + s[i] - M*(1-x[i,j,k]) <= T[j])
                
con7_1 = {}
for i in N1:
    for j in N1:
        if j >= d:
            for k in K1:
                con7_1[i,j,k] = model.addConstr(L[i] - q[j] + M*(1-x[i,j,k]) >= L[j]) 

con8_1 = {}
for i in N1:
    con8_1[i] = model.addConstr(L[i] + q[i] <= Q_v)               

#-------------------------------------------------------

#second echelon
con1_2 = {}
for i in C:
    con1_2[i,j,k] = model.addConstr(quicksum(x[i,j,k] for j in N2 for k in K2) == 1)
    
con2_2 = {}
for j in C:
    con2_2[i,j,k] = model.addConstr(quicksum(x[i,j,k] for i in N2 for k in K2) == 1)
        
con3_2 = {}
for i in N2:
    for k in K2:
        con3_2[i,i,k] = model.addConstr(x[i,i,k] == 0)    

con4_2 = {}
for j in N2:
    for k in K2:
        con4_2[i,j,k] = model.addConstr(quicksum(x[i,j,k] for i in N2) == quicksum(x[j,i,k] for i in N2))

con5_2 = {}
for k in K2:
    con5_2[k] = model.addConstr(quicksum(x[h,j,k] for j in C for h in H) <= 1)    
 
con6_2 = {}
for i in N2:
    for j in N2:
        if j >= (d+h):
            if i != j:
                for k in K2:
                    con6_2[i,j,k] = model.addConstr(T[i] + c[i,j]*x[i,j,k] + s[i] - M*(1-x[i,j,k]) <= T[j])

con7_2 = {}
for i in N2:
    for j in N2:
        if j >= (d+h):
            for k in K2:
                con7_2[i,j,k] = model.addConstr(L[i] - q[j] + M*(1-x[i,j,k]) >= L[j])

con8_2 = {}
for i in N2:
    con8_2[i] = model.addConstr(L[i] + q[i] <= Q_cb)            

#%% ----- Solve -----

model.update ()

model.setParam( 'OutputFlag', True) # silencing gurobi output or not
model.setParam ('MIPGap', 0);       # find the optimal solution (or within a percentage 0,05 = 5%, 0 = optimal)
model.write("output.lp")            # print the model in .lp format file

model.optimize ()

#%% ----- Results -----

print ('\n--------------------------------------------------------------------\n')
if model.status == GRB.Status.OPTIMAL:                          # If optimal solution is found
    print ('Minimal distance : %10.2f ' % model.objVal)
    print('\nFinished\n')
else:
    print ('\nNo feasible solution found\n')

active_arcs = [(i,j,k) for i in N for j in N if i != j for k in K if x[i,j,k].x == 1]
print (('Route : '), sorted(active_arcs))

fig, ax = plt.subplots(figsize=(6,6))
plt.xlim([0, 10])
plt.ylim([0, 10])

for i, j, k in active_arcs:
    plt.plot([xc[i], xc[j]], [yc[i], yc[j]], c='grey', zorder=0)
for d1 in range(d):
    plt.plot(xc[d1], yc[d1], c='black', marker='s', markersize=10, label='depot')
for d2 in range(d,d+h):
    plt.plot(xc[d2], yc[d2], c='grey', marker='s', markersize=10)
plt.plot(xc[d], yc[d], c='grey', marker='s', markersize=10, label='transshipment location')
plt.scatter(xc[d+h:], yc[d+h:], c='purple', label='customer')

plt.legend()

#plot LOC_ID next to points
for i, txt in enumerate(LOC_ID):
    plt.annotate(txt, (xc[i], yc[i]), xytext=(xc[i]+0.12, yc[i]+0.25), bbox=dict(boxstyle="round", alpha=0.1))

#print decision variable T = Time counter (TC)
print ('\nArrival time at location:')
for d1 in range(d):
    print ('%35.0f' % LOC_ID[d1] + '%8.0f' % T[d1].x)
for d2 in range(d,d+h):
    print ('%35.0f' % LOC_ID[d2] + '%8.0f' % T[d2].x)
    
for i in C: 
    AT_C1 = LOC_ID[i]               #arrival time column 1
    AT_C2 = T[i].x                  #arrival time column 2
    AT = '%35.0f' % AT_C1 + '%8.0f' % AT_C2
    print(AT)

#print decision variable L = Load counter (LC)
print ('\nLoad after delivery at location:')
for d1 in range(d):
    print ('%35.0f' % LOC_ID[d1] + '%8.0f' % L[d1].x)
for d2 in range(d,d+h):
    print ('%35.0f' % LOC_ID[d2] + '%8.0f' % L[d2].x)
    
for i in C:
    LC_C1 = LOC_ID[i]               #load counter column 1
    LC_C2 = L[i].x                  #load counter column 2
    LC = '%35.0f' % LC_C1 + '%8.0f' % LC_C2
    print(LC)






