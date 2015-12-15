import numpy as np
from assimulo.problem import Implicit_Problem #Imports the problem formulation from Assimulo
from assimulo.solvers import IDA              #Imports the solver IDA from Assimulo
from woodpecker import *
import matplotlib.pyplot as P

def run_example():

	#initial values
	t0 = 0
	y0 = np.array([0., -0.10344, -0.65, 0., 0. , 0., -0.628993, 0.047088]) #|phi_s| <= 0.1034 rad, |phi_b| <= 0.12 rad
	yd0 = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
	sw = [False, True, False]
	
	#problem
	model = Implicit_Problem(pecker, y0, yd0, t0, sw0=sw)
	model.state_events = state_events #from woodpecker.py
	model.handle_event = handle_event #from woodpecker.py
	model.name = 'Woodpeckermodel'
	sim = IDA(model) #create IDA solver
	tfinal = 2.0 #final time
	ncp = 500 #number control points	
	sim.suppress_alg = True
	sim.rtol=1.e-6
	sim.atol[3:8] = 1e6
	sim.algvar[3:8] = 0
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
	
	P.figure(2)
	fig3, ax3 = P.subplots()
	ax3.plot(t, yd[:, 4], label='phi_sp')
	ax3.plot(t, yd[:, 5], label='phi_bp')
	legend = ax3.legend(loc='upper center', shadow=True)
	P.grid()
	
	#event data
	sim.print_event_data()
	
	#show plots
	P.show()
	print("...")
	P.show()
	

if __name__=='__main__':
	run_example()
