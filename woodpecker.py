from scipy import *
from scipy.linalg import *
import numpy as np

def pecker(t, y, yd, sw):
	#y = [z, phi_s, phi_b, z', phi_s', phi_b', lambda_1, lambda_2]
	#yd = [z', phi_s', phi_b', z'', phi_s'', phi_b'', lambda_1', lambda_2']
	
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
	
	#initial computations and assignments
	lamb = y[6:8]
	z = y[0]
	phi_s = y[1]
	phi_b = y[2]
	zp = y[3]
	phi_sp = y[4]
	phi_bp = yd[5]
	
	
	#Mass matrix
	m = np.zeros((3, 3))
	m[0, 0] = m_s + m_b
	m[1, 0] = m_b * l_s
	m[2, 0] = m_b * l_g
	m[0, 1] = m_b * l_s
	m[1, 1] = J_s + m_b * l_s**2
	m[2, 1] = m_b * l_s * l_g
	m[0, 2] = m_b * l_g
	m[1, 2] = m_b * l_s * l_g
	m[2, 2] = J_b + m_b * l_g**2

			
	#Applied forces (f matrix)
	ff = np.array([-g * (m_s + m_b), 
		c_p * (phi_b - phi_s) - m_b * l_s * g, 
		c_p * (phi_s - phi_b) - m_b * l_g * g])
	
	#Constraint matrix G
	gp = np.zeros((2, 3))

	#index 1 constraints
	gyy = np.zeros((2,))
	
	if sw[0]: #state 1
		gp[0, 0] = 0
		gp[1, 0] = 0
		gp[0, 1] = 1
		gp[1, 1] = 0
		gp[0, 2] = 0
		gp[1, 2] = 1
		
		gyy[0] = 0
		gyy[1] = 0
			
	elif sw[1]: #state 2
		gp[0, 0] = 0
		gp[1, 0] = 1
		gp[0, 1] = h_s
		gp[1, 1] = r_s
		gp[0, 2] = 0
		gp[1, 2] = 0
	
		gyy[0] = yd[4]
		gyy[1] = yd[3] + r_s * yd[4]
			
	else: #state 3
		gp[0, 0] = 0
		gp[1, 0] = 1
		gp[0, 1] = -h_s
		gp[1, 1] = r_s
		gp[0, 2] = 0
		gp[1, 2] = 0
	
		gyy[0] = yd[4]
		gyy[1] = yd[3] + r_s * yd[4]
		
	res_1 = yd[0:3] - y[3:6]
	res_2 = dot(m, yd[3:6]) - ff[0:3] + dot(gp.T, lamb)
	res_3 = lamb - gyy
	
	return hstack((res_1, res_2, res_3))
	
def state_events(t, y, sw):
	'''
	This is the function that keeps track of events. When the sign of any of the functions
	changed, we have an event.
	'''
	if sw[0]:
		#transition 1: State 1 and phi_b' < 0 switch to state 2 when h_s*phi_s = -(r_s - r0)
	else if sw[1]:
		#transition 2: State 1 and phi_b' > 0 switch to state 3 when h_s*phi_s = (r_s - r0)
	else if sw[2]:
		#transition 3: State 2 switch to state 1 if lambda_1 changes sign
	else if sw[3]:
		#transition 4: State 3 and phi_b' < 0 switch to state 1 if lambda_1 changes sign
	else if sw[4]:
		#transition 5: State 3 and phi_b' < 0 switch to state 4 (beak hit, switch to state 3 and change sign of phi_s') if h_b * phi_b = l_s + l_g - l_b - r0
		
def handle_event(solver, event_info):
	'''
	Event handling. This functions is called when Assimulo finds an event as
	specified by the event functions
	'''
	state_info = event_info[0]
	
	if state_info[0] != 0: #Check if the first event function has been triggered
	
		if solver.sw[0]: #transition 1
			