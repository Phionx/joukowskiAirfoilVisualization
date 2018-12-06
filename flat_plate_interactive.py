##############################################################
#                             ASEN 5051
#                           Mini Project 1 - TESTING
#
# Authors: Lucas Calvert, Duncan McGough
##############################################################
# IMPORT PACKAGES

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

##############################################################
# DEFINE COMMON VARIABLES
pi = np.pi
alpha = pi/4
U_inf = 1
b = 1

##############################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Polar coordinates in the Zeta-plane
#
# Note: This avoids plotting inside
#       the cylinder where the
#       Joukowsky map is singular
res = 100
radius = np.linspace(b,4*b,res)
angle = np.linspace(0,2*pi,res)

[RADIUS,ANGLE] = np.meshgrid(radius,angle)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Xi-Eta coordiantes in the Zeta-plane

Xi = RADIUS*np.cos(ANGLE)
Eta = RADIUS*np.sin(ANGLE)

Zeta = Xi + 1j*Eta

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Complex potential in the Zeta-plane

Psi = U_inf*(np.exp(-1j*alpha)*Zeta + np.divide(np.exp(1j*alpha)*b**2,Zeta) + 2*np.sin(alpha)*b*1j*np.log(Zeta/b))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Apply the Joukowsky map
Z = Zeta + np.divide(b**2,Zeta)

X = np.real(Z)
Y = np.imag(Z)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Grab streamfunction

psi = np.imag(Psi)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Plot!
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.3)
p1 = plt.contour(X,Y,psi,50,colors=('k',),linewidths=(0.3,))
contour_axis = plt.gca()
plt.contourf(X,Y,psi,50)  #plotting contour
nLine = 1000
plt.plot(np.linspace(-2*b,2*b,nLine),np.zeros(nLine),linewidth=3,color='k')
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Flat Plate Streamlines',fontsize=16)

axcolor = [200.0/255.0,255.0/255.0,255.0/255.0] #'lightgoldenrodyellow'
axlev = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
axaoa = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

slev = Slider(axlev, 'N Contour', 1, 100, valinit=50)
saoa = Slider(axaoa, 'AoA', 0.1, pi/2, valinit=pi/4)


def update(val):
    contour_axis.clear()
    lev = slev.val
    aoa = saoa.val

    [RADIUS,ANGLE] = np.meshgrid(radius,angle)
    Xi = RADIUS*np.cos(ANGLE)
    Eta = RADIUS*np.sin(ANGLE)
    Zeta = Xi + 1j*Eta
    Psi = U_inf*(np.exp(-1j*aoa)*Zeta + np.divide(np.exp(1j*aoa)*b**2,Zeta) + 2*np.sin(aoa)*b*1j*np.log(Zeta/b))
    Z = Zeta + np.divide(b**2,Zeta)
    X = np.real(Z)
    Y = np.imag(Z)
    psi = np.imag(Psi)

    contour_axis.contour(X,Y,psi,int(lev),colors=('k',),linewidths=(0.3,))
    contour_axis.contourf(X,Y,psi,int(lev))  #plotting contour
    contour_axis.plot(np.linspace(-2*b,2*b,nLine),np.zeros(nLine),linewidth=3,color='k')
    ax.set_xlabel('x',fontsize=16)
    ax.set_ylabel('y',fontsize=16)
    ax.set_title('Flat Plate Streamlines',fontsize=16)
    fig.canvas.draw_idle()
slev.on_changed(update)
saoa.on_changed(update)

plt.show()


##############################################################
