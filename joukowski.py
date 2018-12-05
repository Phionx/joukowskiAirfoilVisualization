#####################################################################
#               JOUKOWSKI PLOTTING SCRIPT
#
# AUTHORS: Lucas Calvert and Duncan McGough
#
#####################################################################
# IMPORT PACKAGES
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons

##########################################
# DEFINE THE DRAGGABLE OBJECT class
##########################################
class DraggablePoints(object):

    ###########################################
    #  INITIALIZATION OF OBJECTS
    def __init__(self, artists, tolerance=1):
        for artist in artists:
            artist.set_picker(tolerance)
        self.artists = artists

        self.currently_dragging = False     #Var to keep track of dragging
        self.current_artist = None       #Var to keep track of which object is being interacted with
        self.offset = (0, 0)

        ###### This sets connections between actions '__' and functions.
        for canvas in set(artist.figure.canvas for artist in self.artists):
            canvas.mpl_connect('button_press_event', self.on_press)
            canvas.mpl_connect('button_release_event', self.on_release)
            canvas.mpl_connect('pick_event', self.on_pick)
            canvas.mpl_connect('motion_notify_event', self.on_motion)
    ######################################################################


    #If clicked, make an object "currently dragging"
    def on_press(self, event):
        self.currently_dragging = True

    #When released, forget the astist and set it to not dragging. also print the value of the new center
    def on_release(self, event):
        self.currently_dragging = False
        self.current_artist = None


    # Object movement
    def on_pick(self, event):
        if self.current_artist is None:
            self.current_artist = event.artist
            x0, y0 = event.artist.center
            x1, y1 = event.mouseevent.xdata, event.mouseevent.ydata
            self.offset = (x0 - x1), (y0 - y1)

    def on_motion(self, event):
        if not self.currently_dragging:
            return
        if self.current_artist is None:
            return
        dx, dy = self.offset
        self.current_artist.center = event.xdata + dx, event.ydata + dy
        #self.current_artist.radius = self.current_artist.radius + np.sqrt(dx**2+dy**2)/2
        global centerX
        global centerY
        centerX,centerY = self.current_artist.center
        #reDrawAirfoil(centerX,centerY)
        #updateCenter(centerXnew,centerYnew)
        reDrawAirfoil()
        self.current_artist.figure.canvas.draw()

#################################################################################################



#########################################
# FUNCTION TO REDRAW THE AIRFOIL
#########################################
def reDrawAirfoil():

    ##### DEFINE NEW VALUES
    global centerX
    global centerY
    #centerX, centerY = dr.current_artist.center
    alpha = s_aoa.val
    res = s_res.val
    nContour = int(s_numCont.val)
    b = s_b.val
    U_inf = s_uinf.val
    P_inf = s_pinf.val
    rho_inf = s_rinf.val
    global plotOption
    #print(plotOption)
    ##########################################


    R = np.sqrt((centerX-b)**2 + (centerY)**2)
    angle1 = np.arccos((-centerX+b)/R)
    angle2 = (alpha+angle1)

    radius = np.linspace(R,4*b,res)
    angle = np.linspace(0,2*pi,res)
    [RADIUS,ANGLE] = np.meshgrid(radius,angle)

    zetaPrime = centerX + centerY*1j

    Xi = RADIUS*np.cos(ANGLE) + centerX
    Eta = RADIUS*np.sin(ANGLE) + centerY

    zeta = (Xi + 1j*Eta);

    Gamma = -2*pi*R*2*U_inf*np.sin(angle2)
    Psi = U_inf*(np.exp(-1j*alpha)*(zeta-zetaPrime) + np.divide(np.exp(1j*alpha)*R**2,(zeta-zetaPrime))) + np.divide(Gamma,(2*pi*1j))*np.log((zeta-zetaPrime)/R)

    Z = (zeta) + np.divide((b)**2,(zeta))

    X = np.real(Z);
    Y = np.imag(Z);

    psi = np.imag(Psi);


    ################################
    # CALCULATE PRESSURE AND VELOCTITIES  UNCOMMENT THIS!
    ################################
    dPsidJ = np.divide(((-U_inf*np.exp(-1j*alpha)) + (np.divide(Gamma*1j, 2*pi*(zeta-zetaPrime))) + (np.divide(R**2*U_inf*np.exp(1j*alpha),(zeta-zetaPrime)**2)) ) , (np.divide(b**2,zeta**2)-1))

    u = np.real(dPsidJ)
    v = -np.imag(dPsidJ)

    magVel = np.sqrt(u**2 + v**2)

    const = P_inf + 0.5*rho_inf*U_inf**2

    P = const - 0.5*rho_inf*magVel**2




    ################################
    # CHANGE PLOT DEPENDING ON RADIO BUTTON
    ################################
    # clear axis
    ax2.clear()

    # set title
    ax1.set_title('Cylinder in Complex Plane',fontsize=16)
    ax2.set_title('Airfoil in Z-Plane',fontsize=16)

    if plotOption == 0: #streamline case
        ax2.contour(X,Y,psi,nContour,colors=('k',),linewidths=(0.3,))
        ax2.contourf(X,Y,psi,nContour)
        ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")
        ax2.axis('equal')
    elif plotOption == 1: #Pressure case
        ax2.contour(X,Y,P,nContour,colors=('k',),linewidths=(0.3,))
        ax2.contourf(X,Y,P,nContour)
        ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")
        ax2.axis('equal')
    elif plotOption == 2: #u-velocity case
        ax2.contour(X,Y,u,nContour,colors=('k',),linewidths=(0.3,))
        ax2.contourf(X,Y,u,nContour)
        ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")
        ax2.axis('equal')
    elif plotOption == 3: #v-velocity case
        ax2.contour(X,Y,v,nContour,colors=('k',),linewidths=(0.3,))
        ax2.contourf(X,Y,v,nContour)
        ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")
        ax2.axis('equal')
    elif plotOption == 4: #velMag case
        ax2.contour(X,Y,magVel,nContour,colors=('k',),linewidths=(0.3,))
        ax2.contourf(X,Y,magVel,nContour)
        ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")
        ax2.axis('equal')


    # show the plot
    plt.show()

    # print pressure and velocity magnitudes
    print('\n Velocity Magnitude Max:\n')
    print(np.max(magVel))
    print('\n Velocity Magnitude Min:\n')
    print(np.min(magVel))
    print('\n Pressure Magnitude Max:\n')
    print(np.max(P))
    print('\n Pressure Magnitude Min:\n')
    print(np.min(P))


#########################################
# FUNCTION TO INITIALLY DRAW THE AIRFOIL
def drawAirfoil():
    ax2.clear()
    R = np.sqrt((centerX-b)**2 + (centerY)**2)
    angle1 = np.arccos((-centerX+b)/R)
    angle2 = (alpha+angle1)

    radius = np.linspace(R,4*b,res)
    angle = np.linspace(0,2*pi,res)
    [RADIUS,ANGLE] = np.meshgrid(radius,angle)

    zetaPrime = centerX + centerY*1j

    Xi = RADIUS*np.cos(ANGLE) + centerX
    Eta = RADIUS*np.sin(ANGLE) + centerY

    zeta = (Xi + 1j*Eta);

    Gamma = -2*pi*R*2*U_inf*np.sin(angle2)
    Psi = U_inf*(np.exp(-1j*alpha)*(zeta-zetaPrime) + np.divide(np.exp(1j*alpha)*R**2,(zeta-zetaPrime))) + np.divide(Gamma,(2*pi*1j))*np.log((zeta-zetaPrime)/R)

    Z = (zeta) + np.divide((b)**2,(zeta))

    X = np.real(Z);
    Y = np.imag(Z);

    psi = np.imag(Psi);

    ax2.contour(X,Y,psi,nContour,colors=('k',),linewidths=(0.3,))
    contour_axis = plt.gca()
    ax2.contourf(X,Y,psi,nContour)
    nLine = 1000
    ax2.plot(np.linspace(-2*b,2*b,nLine),np.zeros(nLine),linewidth=3,color='k')
    ax2.set_title('Airfoil in Z-Plane',fontsize=16)
    ax1.plot(np.linspace(np.min(X),np.max(X),nLine),np.zeros(nLine),linewidth=1,color='k')
    ax1.plot(np.zeros(nLine),np.linspace(np.min(Y),np.max(Y),nLine),linewidth=1,color='k')
    ax1.set_title('Cylinder in Complex Plane',fontsize=16)
    ax2.axis('equal')
    ax1.plot(centerX,centerY,marker='o', markersize=3, color="red")

    plt.show()
#########################################




##########################################
# MAIN
    ##########################################
#if __name__ == '__main__':

################################################
# Define global variable and common constants
pi = np.pi
global alpha
alpha = pi/16
U_inf = 100
P_inf = 1000000
rho_inf = 1
b = 1
res = 100
nContour = 50
centerX = 0
centerY = 0

#Define the subplots and the size of the window
fig, (ax1,ax2) = plt.subplots(1,2,False,False,figsize=(15,9))

#Adjust to "push" subplots up and to the right
fig.subplots_adjust(left=0.2, bottom=0.35)

#Set axis and make then equal, etc.
ax1.axis('equal')
ax2.axis('equal')
ax1.axis([-1,1,-1,1])

################################################
# Define the sliders
################################################

aoa = plt.axes([0.25,0.03,0.65,0.03])
bax = plt.axes([0.25, 0.03+1*(0.04),0.65,0.03])
uinfax = plt.axes([0.25, 0.03+2*(0.04),0.65,0.03])
pinfax = plt.axes([0.25, 0.03+3*(0.04),0.65,0.03])
rinfax = plt.axes([0.25, 0.03+4*(0.04),0.65,0.03])
resolution = plt.axes([0.25,0.03+5*(0.04),0.65,0.03])
numCont = plt.axes([0.25,0.03+6*(0.04),0.65,0.03])

s_aoa = Slider(aoa, 'Angle of Attack',-pi/16,pi/2, valinit=alpha)
s_b = Slider(bax,'b',0,5,valinit=1)
s_uinf = Slider(uinfax,'U_free',0,500,valinit=100)
s_pinf = Slider(pinfax,'P_free',0,5e6,valinit=1e6)
s_rinf = Slider(rinfax,'rho_free',0,5,valinit=1)
s_res = Slider(resolution, 'Resolution',1,500, valinit=res)
s_numCont = Slider(numCont, 'Number of Contour Lines',1,150, valinit=nContour)

################################################
#Define the circle object, add it to subplot, and make it draggable
################################################
circle = [patches.Circle((0, 0), 0.3, fc='r', alpha=0.5)]
for circ in circle:
    ax1.add_patch(circ)
    dr = DraggablePoints(circle)

################################################
#Update after slider movement
###############################################
def update(val):
    reDrawAirfoil()
s_aoa.on_changed(update)
s_b.on_changed(update)
s_uinf.on_changed(update)
s_pinf.on_changed(update)
s_rinf.on_changed(update)
s_res.on_changed(update)
s_numCont.on_changed(update)



#################################################
## Create Radio Button Set For Selecting Plot Options
#################################################
pOp = plt.axes([0.025, 0.5, 0.1, 0.15])
pOpRadio = RadioButtons(pOp, ('Streamlines', 'Pressure', 'X-Velociy','Y-Velocity','|U|'), active=0)

global plotOption
plotOption = 0

def setPlot(label):

    global plotOption
    global plotOption
    if label == 'Streamlines':
        plotOption = 0
    elif label == 'Pressure':
        plotOption = 1
    elif label == 'X-Velociy':
        plotOption = 2
    elif label == 'Y-Velocity':
        plotOption = 3
    elif label == '|U|':
        plotOption = 4

    reDrawAirfoil()
    #print(plotOption)
pOpRadio.on_clicked(setPlot)



#################################################
## Create Reset button
#################################################
#resetax = plt.axes([0.01, 0.025, 0.1, 0.04])
#resetButton = Button(resetax, 'Reset', hovercolor='0.975')
#
#def reset(event):
#    ax1.clear()
#    ax2.clear()
#    s_aoa.reset()
#    s_res.reset()
#    s_numCont.reset()
#    global centerX
#    global centerY
#    centerX = 0
#    centery = 0
#    circle = [patches.Circle((0, 0), 0.3, fc='r', alpha=0.5)]
#    for circ in circle:
#        ax1.add_patch(circ)
#        dr = DraggablePoints(circle)
#    ax1.axis([-1,1,-1,1])
#    drawAirfoil()
#resetButton.on_clicked(reset)


drawAirfoil()
