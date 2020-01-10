from __future__ import division #this allows use of integers in fractions

#James McFarlin, T-handle intermediate axis theorem, 10/11/2015

from visual import * #visual library
from visual.graph import * #graphing library
import math #math library


"""scene specs"""

scene.width = 1200
scene.height = 450
scene.range = 14
scene.title = 'Rotating T Handle in Zero G'
scene.center = (5,0,0)
scene.background = color.white


"""creating objects for t-handle, "filler" is simply for cleanliness"""

tbar = frame(omega=vector(.1,.1,4)) #"frame" makes a compound object of the objects below
cross = cylinder(frame=tbar, pos=(0,3,0), axis=(0,-6,0), radius=1, mass=4) #this is the crossbar of the handle
perp = cylinder(frame=tbar, pos=(-1,0,0), axis=(-3,0,0), radius=.75, mass=1) #perpindicular portion
filler = cylinder(frame=tbar, pos=(-1,0,0), axis=(.5,0,0), radius=.75, mass=0) #fills space btwn two cylinders (is not involved in calculations)

#let's make it all yellow
for obj in tbar.objects:
    obj.color = color.yellow


"""parameters and constants"""

q = -(perp.length/6 + cross.radius/3) #distance between center of mass of cross and center of mass of total object

#moments of intertia, Izz lies along center of mass, Ixx & Iyy need to use parallel axis theorem to find I about that axis
Ixx = (1/2)*(cross.mass)*(cross.radius**2) + (cross.mass)*q**2 + (perp.mass/4)*(perp.radius**2) + (perp.mass/12)*(perp.length**2) + perp.mass*(perp.length/2 + cross.radius - q)
Iyy = (1/4)*(cross.mass)*(cross.radius**2) + (1/12)*(cross.mass)*(cross.length**2) + cross.mass*(q**2) + (1/4)*perp.mass*(perp.radius**2) + (1/12)*(perp.mass)*(perp.length**2) + perp.mass*(perp.length/2 + cross.radius - q)
Izz = (1/2)*(perp.mass)*(perp.radius**2) + (1/4)*(cross.mass)*(cross.radius**2) + (1/12)*(cross.mass)*(cross.length**2)

#form a basis
basis_vectors = []
basis_vectors.append(vector(1.,0,0))
basis_vectors.append(vector(0,1.,0))
basis_vectors.append(vector(0,0,1.))


"""Create axis objects"""

#these objects indicate the initial orientation of the object. All share the same origin, the length 4 is for visibility
xaxis = arrow(pos=(0,0,0), axis=(0,4,0), shaftwidth=.1, color=color.red, opacity=.5)
yaxis = arrow(pos=(0,0,0), axis=(0,0,4), shaftwidth=.1, color=color.green, opacity=.5)
zaxis = arrow(pos=(0,0,0), axis=(4,0,0), shaftwidth=.1, color=color.blue, opacity=.5)

axis_List = []    #putting the axes in a list for simplicity
axis_List.append(xaxis)
axis_List.append(yaxis)
axis_List.append(zaxis)


"""Create vector objects"""

#these vectors will indicate the magnitudes of angular momenta along the basis vectors, as well as a total w.
xvec = arrow(pos=(10,0,0), length=(Ixx*tbar.omega.x)/5, axis=basis_vectors[0], shaftwidth=.1, color=color.red, opacity=.7)  
yvec = arrow(pos=(10+xvec.length,0,0), length=(2*Ixx*tbar.omega.y)/5, axis=basis_vectors[1], shaftwidth=.1, color=color.green, opacity=.7)
zvec = arrow(pos=(10+xvec.length,yvec.length,0), length=(Ixx*tbar.omega.z)/5, axis=basis_vectors[2], shaftwidth=.1, color=color.blue, opacity=.7)

totvec = arrow(pos=(10,0,0), length=sqrt(xvec.length*xvec.length+yvec.length*yvec.length+zvec.length*zvec.length), axis=norm((xvec.length,yvec.length,zvec.length)), shaftwidth=.2, color=color.yellow)

vector_List = [] #putting vectors in list for simplicity
vector_List.append(xvec)
vector_List.append(yvec)
vector_List.append(zvec)
vector_List.append(totvec)


"""Create plot"""

#this line sets a title and makes the graph not lay up over the top of the animation window.
graph1 = gdisplay(title = 'Angular Momenta and Energy vs. Time',x=200,y=450,height = 325)

#angular momenta in x,y,z axes, total momentum (conserved) and total energy (conserved)
p_x = gcurve(gdisplay = graph1,color=color.red)
p_y = gcurve(gdisplay = graph1,color=color.green)
p_z = gcurve(gdisplay = graph1,color=color.blue)
p_tot = gcurve(gdisplay = graph1,color=color.yellow)
E_tot = gcurve(gdisplay = graph1,color=color.white)

plot_List = [] #putting plot curves in list for simplicity
plot_List.append(p_x)
plot_List.append(p_y)
plot_List.append(p_z)
plot_List.append(p_tot)
plot_List.append(E_tot)


"""Define functions"""

def find_omega_dot(w):  #function for finding omegadot using Euler's equations
    wx_dot = (Iyy-Izz)*w.y*w.z/Ixx
    wy_dot = (Izz-Ixx)*w.z*w.x/Iyy
    wz_dot = (Ixx-Iyy)*w.x*w.y/Izz
    w_dot = vector(wx_dot, wy_dot, wz_dot)
    return(w_dot)

def update_omega(w,wdot): #updating angular momentum
    w.x = w.x + dt*wdot.x
    w.y = w.y + dt*wdot.y
    w.z = w.z + dt*wdot.z
    
def rotate_tbar(obj,axes):   #rotating the tbar in the graphics
    for i in range(len(axes)):
        obj.rotate(angle=obj.omega[i]*dt, axis=norm(axes[i].axis), origin=obj.pos)

def update_axis(List,w): #rotating the vectors and adjusting length to be proportional to angular momentum
    old = List
    for i in range(len(List)):
        k = 0
        while k<3:
            List[i].rotate(angle=w[k]*dt, axis=norm(old[k].axis), origin=List[i].pos)
            k += 1

def update_vectors(List,w): #updates the position of the vectors to show change in ang momenta
    k = vector(0,0,0)
    for i in range(len(List)-1):
        k[i] = ((Ixx*w[i])/5)/abs((Ixx*w[i])/5)
        List[i].axis = k[i]*basis_vectors[i]
        List[i].length = abs((Ixx*w[i])/5)
        if i>0:
            List[i].pos.x = 10+(List[0].length)*k[0]
            if i>1:
                List[i].pos.y = (List[1].length)*k[1]
    List[3].axis = norm((k[0]*List[0].length, k[1]*List[1].length, k[2]*List[2].length))
    List[3].length = sqrt(List[0].length*List[0].length+List[1].length*List[1].length+List[2].length*List[2].length)

def update_plot(pList, w, time): #sets equations for plots above
    pList[0].plot(pos=(time, w.x*Ixx))
    pList[1].plot(pos=(time, w.y*Iyy))
    pList[2].plot(pos=(time, w.z*Izz))
    pList[3].plot(pos=(time, sqrt(w.x*w.x*Ixx*Ixx+w.y*w.y*Iyy*Iyy+w.z*w.z*Izz*Izz)))
    pList[4].plot(pos=(time, .5*(w.x*w.x*Ixx+w.y*w.y*Iyy+w.z*w.z*Izz)))

       
t = 0.0000
dt = .001

while t<60:
    rate(1000)
    omega_dot = find_omega_dot(tbar.omega)
    update_omega(tbar.omega, omega_dot)
    rotate_tbar(tbar,axis_List)
    update_axis(axis_List, tbar.omega)
    update_vectors(vector_List, tbar.omega)
    update_plot(plot_List, tbar.omega, t)
    t += dt
