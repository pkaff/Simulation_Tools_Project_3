import numpy as np
from assimulo.problem import Implicit_Problem #Imports the problem formulation from Assimulo
from assimulo.solvers import IDA              #Imports the solver IDA from Assimulo
from woodpecker import *
import matplotlib.pyplot as P

def run_example():
	def residual(t, y, yd, sw):
		return pecker(t, y, yd, sw)

	#initial values
	t0 = 0
	y0 = np.array([1, 0, 0.02, 0, 0, 0, 0, 0])
	yd0 = np.array([0, 0, 0, 0, 0, 0, 0, 0])
	sw0 = np.array([True, False, False, False])
	
	#problem
	model = Implicit_Problem(residual, y0, yd0, t0, sw0)
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
