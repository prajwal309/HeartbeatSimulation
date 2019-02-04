import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

Dist = 5.0
e = 0.5

MeanAngles = np.linspace(0, 2.0*np.pi,200)

def EccentricAngleFunc(E):
    global MeanAngle
    return MeanAngle - E +  e*np.sin(E)

X1 = Dist*np.cos(MeanAngles)
Y1 = Dist*np.sin(MeanAngles)

global MeanAngle
EccentricAngles = []
for MeanAngle in MeanAngles:
    EccentricAngle = fsolve(EccentricAngleFunc,np.pi/4.0)
    EccentricAngles.append(EccentricAngle)

EccentricAngles = np.array(EccentricAngles)
Angle = 2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(EccentricAngles/2.0))

a_R = 1.5*Dist*(1.0-e*e)/(1+e*np.cos(Angle))
X2 = a_R*np.cos(Angle)
Y2 = a_R*np.sin(Angle)


fig,axis = plt.subplots(figsize=(10,10), dpi=200)
axis.plot(X1,Y1, "r.", lw=2.5)
axis.plot(X2+1.25,Y2, "k-", lw=2.5)
axis.plot([0],[0],color='orange', markersize=50, marker="o")
axis.plot([-6.15],[5.65],color='brown', markersize=15, marker="o")
axis.set_xticks([])
axis.set_yticks([])
plt.axis('off')
axis.set_aspect('equal')
plt.savefig("Tidal_Circularization.png")
