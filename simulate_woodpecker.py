import numpy as np
from assimulo.problem import Implicit_Problem #Imports the problem formulation from Assimulo
from assimulo.solvers import IDA              #Imports the solver IDA from Assimulo
from woodpecker import *
import matplotlib.pyplot as P
#sleeve constants
m_s = 3.0e-4 #mass sleeve
r_s = 3.1e-3 #inner radius of the sleeve
h_s = 2.0e-2 #half height of the sleeve
J_s = 5.0e-9 #moment of inertia for sleeve
l_s = 1.0e-2 #vertical distance between sleeve and rotation axis

#bird constants
l_b = 2.01e-2 #beak x-coords in the birds coordinate system
h_b = 2.1e-2 # beak y-coords in the birds coordinate system
J_b = 7.0e-7 #moment of inertia for bird
m_b = 4.5e-3 #mass bird

#other constants
l_g = 1.5e-2 #vertical distance between bird's origin and rotation axis of the spring
c_p = 5.6e-3 #spring constant
g = 9.81 #graviy constant
r0 = 2.5e-3 #radius of the bar
def run_example():
	def residual(t, y, yd, sw):
		return pecker(t, y, yd, sw)

	#initial values
	t0 = 0
	y0 = np.array([10, -0.1, -0.01, 0, 0, 0, 0, 0])
	yd0 = np.array([0, 0, 0, 0, 0, 0, 0, 0])
	sw = np.array([True, False, False])
	
	#problem
	model = Implicit_Problem(residual, y0, yd0, t0, sw0=sw)
	model.state_events = state_events #from woodpecker.py
	model.handle_event = handle_event #from woodpecker.py
	model.name = 'Woodpeckermodel'
	sim = IDA(model) #create IDA solver
	tfinal = 2.0 #final time
	ncp = 500 #number control points	
	t, y, yd = sim.simulate(tfinal, ncp) #simulate
	
	#plot
	fig, ax = P.subplots()
	ax.plot(t, y[:, 0], label='z')
	legend = ax.legend(loc='upper center', shadow=True)
	P.grid()
	
	P.figure(1)
	fig2, ax2 = P.subplots()
	ax2.plot(t, y[:, 1], label='phi_s')
	ax2.plot(t, y[:, 2], label='phi_b')
	legend = ax2.legend(loc='upper center', shadow=True)
	P.grid()
	
	#event data
	sim.print_event_data()
	
	#show plots
	P.show()
	print("...")
	P.show()
	

if __name__=='__main__':
	run_example()
