import math
import numpy as np
import scipy as sp
import json
import ast


def Calculate(nSteps, dT, powerKW, totalWeight, maxacc, airDens, tyrePressure,
              kartAirDrag, kartFrontArea, gravity):
    # data
    velocityVals = np.zeros((nSteps), dtype=float)
    accelVals = np.zeros((nSteps), dtype=float)
    dragVals = np.zeros((nSteps), dtype=float)
    resistanceVals = np.zeros((nSteps), dtype=float)

    # stepwise calculation functions
    def Velocity(preVel, preDrag, preRollRes):
        vcalc = (np.sqrt(2 * (powerKW * 1000 * dT + 0.5 * totalWeight * (
                preVel ** 2)) / totalWeight) - preDrag / totalWeight * dT - preRollRes / totalWeight * dT - preVel) / dT
        if vcalc > maxacc:
            vel = preVel + maxacc * dT
        else:
            vel = np.sqrt(2 * (powerKW * 1000 * dT + 0.5 * totalWeight * (
                    preVel ** 2)) / totalWeight) - preDrag / totalWeight * dT - preRollRes / totalWeight * dT
        return vel

    def Acceleration(vel, prevel):
        acc = (vel - prevel) / dT
        return maxacc if acc > maxacc else acc

    # calc t=0 values
    velocityVals[0] = 0
    accelVals[0] = maxacc
    powerForce = accelVals[0] * totalWeight
    dragVals[0] = kartAirDrag * airDens * (velocityVals[0] ** 2 / 2) * kartFrontArea
    rollCoeff = 0.005 + (1 / tyrePressure) * (0.01 + 0.0095 * (velocityVals[0] / 100) ** 2)
    resistanceVals[0] = rollCoeff * totalWeight * gravity

    for i in range(1, nSteps):
        velocityVals[i] = Velocity(velocityVals[i - 1], dragVals[i - 1], resistanceVals[i - 1])
        accelVals[i] = Acceleration(velocityVals[i], velocityVals[i - 1])
        powerForce = accelVals[i] * totalWeight
        dragVals[i] = kartAirDrag * airDens * (velocityVals[i] ** 2 / 2) * kartFrontArea
        rollCoeff = 0.005 + (1 / tyrePressure) * (0.01 + 0.0095 * (velocityVals[i] / 100) ** 2)
        resistanceVals[i] = rollCoeff * totalWeight * gravity

    return velocityVals, accelVals, dragVals, resistanceVals
