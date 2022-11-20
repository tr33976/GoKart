import math
import numpy as np
import scipy as sp
import json
import ast
import matplotlib.pyplot as plt

# global consts
gravity = 9.81
maxTime = 40
airDens = 1.225  # kg m^3

# read in tyre lookup file
with open('tyres.txt') as f:
    tyres = f.read()
Tyres = ast.literal_eval(tyres)
tyre = 'HiPerf'


# kart consts #
# dims in mm
frontWidth = 1200
wheelHeight = 400
cageHeight = 800
cageWidth = 600
# weights in kg
kartWeight = 150
driverWeight = 88
totalWeight = kartWeight + driverWeight
# coefficients
kartAirDrag = 0.8
wetDry = 'Dry'
wheelFriction = Tyres[tyre]['DryFCoeff' if wetDry == 'Dry' else 'WetFCoeff']
tyrePressure = Tyres[tyre]['Pressure']


# ODE constants #
dT = 0.0001  # constant stepsize
time = 0
nSteps = int(np.ceil(40/dT))

# back front P split for power
mBack = 70
mFront = 100 - mBack
pBack = 50
pFront = 6
pConst = 1
powerKW = pBack + pFront

# calc values
kartFrontArea = (frontWidth * wheelHeight + cageWidth * cageHeight)/1000000  # m^2
maxacc = wheelFriction * gravity * (
        1 - abs(1 - (pBack / (pFront + pBack) / mBack * 100 + pFront / (pFront + pBack) / mFront * 100) / 2)
) * pConst

# data
velocityVals = np.zeros((nSteps), dtype=float)
accelVals = np.zeros((nSteps), dtype=float)
dragVals = np.zeros((nSteps), dtype=float)
resistanceVals = np.zeros((nSteps), dtype=float)


# stepwise calculation functions
def Velocity(preVel, preDrag, preRollRes):
    vcalc = (np.sqrt(2*(powerKW*1000*dT+0.5*totalWeight*(
            preVel**2))/totalWeight) - preDrag/totalWeight*dT - preRollRes/totalWeight*dT-preVel)/dT
    if vcalc > maxacc:
        vel = preVel + maxacc * dT
    else:
        vel = np.sqrt(2*(powerKW*1000*dT+0.5*totalWeight*(
            preVel**2))/totalWeight) - preDrag/totalWeight*dT - preRollRes/totalWeight*dT
    return vel

def Acceleration(vel, prevel):
    acc = (vel-prevel)/dT
    return maxacc if acc > maxacc else acc


# calc t=0 values
velocityVals[0] = 0
accelVals[0] = maxacc
powerForce = accelVals[0] * totalWeight
dragVals[0] = kartAirDrag * airDens * (velocityVals[0]**2 / 2) * kartFrontArea
rollCoeff = 0.005 + (1/tyrePressure) * (0.01+0.0095*(velocityVals[0]/100)**2)
resistanceVals[0] = rollCoeff * totalWeight * gravity

for i in range(1, nSteps):
    velocityVals[i] = Velocity(velocityVals[i-1], dragVals[i-1], resistanceVals[i-1])
    accelVals[i] = Acceleration(velocityVals[i], velocityVals[i-1])
    powerForce = accelVals[i] * totalWeight
    dragVals[i] = kartAirDrag * airDens * (velocityVals[i] ** 2 / 2) * kartFrontArea
    rollCoeff = 0.005 + (1 / tyrePressure) * (0.01 + 0.0095 * (velocityVals[i] / 100) ** 2)
    resistanceVals[i] = rollCoeff * totalWeight * gravity


plt.plot(np.linspace(0, 40, nSteps), velocityVals)
plt.plot(np.linspace(0, 40, nSteps), accelVals)
plt.show()
max(velocityVals)
