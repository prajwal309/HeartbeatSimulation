#Author:Prajwal Niraula
#Institute: EAPS, MIT

import numpy as np
import mayavi.mlab as mlab
import  moviepy.editor as mpy
from scipy.optimize import fsolve
#Create the Surface of the Sphere


#Create a spherical object
RStar = 1.0         #Radius of the star
RPlanet = 2.5       #Radius of the planet
Dist = 7.0           #Semi Scaled angle
FPS = 30            #Frames per Second
MassRatio = 0.3    #Mass ratio of planet and star
TDur = 7.0          #Transit Duration
e = 0.75             #eccentricity
#Generate the light curve

#Phi goes from 0 to pi, Theta goes from 0 to 2 pi
Resolution = 100
Phi, Theta = np.meshgrid(np.linspace(0,np.pi,Resolution),np.linspace(0,2.0*np.pi,Resolution))



RDot = 0.2
XDot = RDot*np.sin(Phi)*np.cos(Theta)
YDot = RDot*np.sin(Phi)*np.sin(Theta)
ZDot = RDot*np.cos(Phi)


#Generate the light curve using batman

def EccentricAngleFunc(E):
    global MeanAngle
    return MeanAngle - E +  e*np.sin(E)


def make_frame(t):

    PhaseShift = np.pi/2.0
    global MeanAngle
    MeanAngle = t/TDur*2.0*np.pi + PhaseShift
    EccentricAngle = fsolve(EccentricAngleFunc,np.pi/4.0)
    Angle = 2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(EccentricAngle/2.0))

    #distance
    a_R = Dist*(1-e*e)/(1+e*np.cos(Angle))

    XLoc = a_R*np.sin(Angle)
    ZLoc = a_R*np.cos(Angle)

    SFact = 15.0*MeanAngle
    XFactor = 8.0/(7.0+np.cos(SFact))
    YFactor = 8.0/(7.0+np.sin(SFact))
    XStar = RStar*np.sin(Phi)*np.cos(Theta)
    YStar = RStar*np.sin(Phi)*np.sin(Theta)
    ZStar = RStar*np.cos(Phi)

    XPlanet = XFactor*RPlanet*np.sin(Phi)*np.cos(Theta)
    YPlanet = YFactor*RPlanet*np.sin(Phi)*np.sin(Theta)
    ZPlanet = RPlanet*np.cos(Phi)

    #Location of the star
    XStar_Loc = MassRatio*a_R*XLoc
    ZStar_Loc = MassRatio*a_R*ZLoc


    #See which angle faces the star vs which is away from it
    PlanetSurface = a_R + np.sqrt((XPlanet-XLoc)**2.0+ YPlanet**2.0+ (ZPlanet-ZLoc)**2.0)

    #Construct Stellar Surface
    StellarSurface = ((PlanetSurface-np.mean(PlanetSurface))**2.0)**0.002

    #Plot the star and the planet
    mlab.clf()
    #Find the corresponding point in the LC
    PhaseCalc = (Angle - PhaseShift)/(2.0*np.pi)
    Shift = -0.050                                               #Shift at the beginning

    mlab.mesh(XStar+XStar_Loc , YStar, ZStar+ZStar_Loc, scalars=StellarSurface, colormap='Reds')

    #showing the mesh of the planet
    mlab.mesh(XPlanet-XLoc, YPlanet, ZPlanet-ZLoc,scalars=StellarSurface,colormap='Blues')
    #Create same length time text
    mlab.view(azimuth=60.0, elevation=15.0, roll=15.0, distance=50, focalpoint=(0,0,0))
    return mlab.screenshot(antialiased=True)


mlab.figure(1, bgcolor=(0, 0, 0), fgcolor=(1, 1, 1), size=(1400, 800))
animation = mpy.VideoClip(make_frame, duration=TDur)

#Generating video and the animation
animation.write_gif("PhaseCurve.gif", fps=FPS)
#animation.write_videofile("PhaseCurve.mp4", fps=FPS)
mlab.close()
