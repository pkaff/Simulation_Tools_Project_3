from __future__ import division
from assimulo.solvers.sundials import IDA
from assimulo.problem import Implicit_Problem
import numpy as N
import pylab as P
from scipy import *
 
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
     
def woodpeckertoy(t,y, yp, sw):
     if sw[0]:
		res_1 = masstotal*g + masstotal*yp[3] + mB*lS*yp[4] + mB*lG*yp[5]
		res_2 = -cp*(y[2]-y[1])+mB*lS*g + y[6] +mB*lS*yp[3] + (JS + mB*lS**2) * yp[4] + mB*lS*lG*yp[5]
		res_3 = -cp*(y[1]-y[2])+mB*lG*g + y[7] +mB*lG*yp[3] + mB*lS*lG*yp[4] + (JB+mB*lG**2)*yp[5]
		res_4 = -y[6]
		res_5 = -y[7]
		res_6 = y[3:6] - yp[0:3]
		
     if sw[1]:
		res_1 = masstotal*g + masstotal*yp[3] + mB*lS*yp[4] + mB*lG*yp[5] + y[7]
		res_2 = -cp*(y[2]-y[1])+mB*lS*g + y[6]*hS + rS*y[7] +mB*lS*yp[3] + (JS + mB*lS**2) * yp[4] + mB*lS*lG*yp[5]
		res_3 = -cp*(y[1]-y[2])+mB*lG*g + mB*lG*yp[3] + mB*lS*lG*yp[4] + (JB+mB*lG**2)*yp[5]
		res_4 = -hS*yp[1]
		res_5 = -yp[0] - rS*yp[1]
		res_6 = y[3:6] - yp[0:3]
		
     if sw[2]:
		res_1 = masstotal*g + masstotal*yp[3] + mB*lS*yp[4] + mB*lG*yp[5] + y[7]
		res_2 = -cp*(y[2]-y[1])+mB*lS*g + y[6]*hS + rS*y[7] +mB*lS*yp[3] + (JS + mB*lS**2) * yp[4] + mB*lS*lG*yp[5]
		res_3 = -cp*(y[1]-y[2])+mB*lG*g + mB*lG*yp[3] + mB*lS*lG*yp[4] + (JB+mB*lG**2)*yp[5]
		res_4 = hS*yp[1]
		res_5 = -yp[0] - rS*yp[1]
		res_6 = y[3:6] - yp[0:3]
		
		
     return hstack([res_1, res_2, res_3, res_4, res_5, res_6])
 
def state_events(t, y, yp, sw):

	e = ones(4,)

	if sw[0]:
		e[0] = -(rS - r0) - y[1]*hS
		e[1] = -(rS - r0) + y[1]*hS
	if sw[1] or sw[2]:
		 e[2] = y[6]
	if sw[2]:
		e[3] = lS + lG - lB - r0 - hB*y[2]
	return e
 
def handle_event(solver, event_info):

	yp = solver.yd
	y = solver.y
	
	state_info = event_info[0]
	
	if state_info[0] != 0:
		if solver.sw[0] == 1 and yp[1] < 0:
			solver.sw[0] = False
			solver.sw[1] = True
			print("I to II")
			y[5]=yp[2] = (mB*lG*yp[0] + (mB*lS*lG)*yp[1] + (JB + mB*lG**2)*yp[2] )/(JB + mB*lG**2) # !!!!!!!!!!
			y[3]=yp[0] = 0 
			y[4]=yp[1] = 0
			
	if state_info[1] != 0:
		if solver.sw[0] == 1 and yp[1] > 0:
			solver.sw[0] = False
			solver.sw[2] = True
			print("I to III")
			y[5]=yp[2] = (mB*lG*yp[0] + (mB*lS*lG)*yp[1] + (JB + mB*lG**2)*yp[2] )/(JB + mB*lG**2) # !!!!!!!!!!
			y[3]=yp[0] = 0 
			y[4]=yp[1] = 0 
	
	elif state_info[2] != 0:
		if solver.sw[2] == 1 and yp[2] < 0:
			solver.sw[2] = False
			solver.sw[0] = True
			print("III to I")
		
		elif solver.sw[1] == 1:
			solver.sw[1] = False
			solver.sw[0] = True
			print("II to I")
		
	elif state_info[3] != 0:
		if solver.sw[2] == 1 and yp[2] > 0:
			yp[2] = -yp[2]
			print("DING")
	solver.yd = yp
	solver.y = y	
 
t0 = 0
tfinal = 1
ncp = 500
y0 = [0., -0.10344, -0.65, 0. ,0. , 0., -0.6911, -0.14161]
yp0 = [0., 0., 0., 0., 0., 1.40059e2, 0., 0.]
switches0 = [False, True, False]
mod = Implicit_Problem(woodpeckertoy, y0, yp0, t0, sw0 = switches0)
mod.state_events = state_events
mod.handle_event = handle_event

sim = IDA(mod)
sim.suppress_alg = True
sim.rtol=1.e-6
sim.atol[3:8] = 1e6
sim.algvar[3:8] = 0.
t, y, yp = sim.simulate(tfinal, ncp)

P.plot(t,y[:,0:3])
P.legend(["z","phiS", "phiB", "z vel", "phiS vel", "phiB vel"])
P.ylabel('y')
P.xlabel('x')
fontlabel_size = 20
tick_size = 14
params = {'lines.markersize' : 0, 'axes.labelsize': fontlabel_size, 'text.fontsize': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size}
P.rcParams.update(params)
P.show()