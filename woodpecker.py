from scipy import *
from scipy.linalg import *
import numpy as np

def pecker(t, y, yd, sw): #?
	#y = [z, phi_s, phi_b]
	
	#?
	yn = y
	ydn = yd
	
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

	#Constraint matrix G
	gp = np.zeros((2, 3))

	#index 1 constraints
	gyy = np.zeros((2,))
	
	if sw[0]: #state 1
		gp[0, 0] = 0
		gp[1, 0] = 0
		gp[0, 1] = 1
		gp[1, 1] = 0
		gp[2, 0] = 0
		gp[2, 1] = 1
		
		gyy[0] = 0
		gyy[1] = 0
		
		#Applied forces (f matrix)
		ff = np.array([-g * (m_s + m_b), 
			c_p * (y[2] - y[1]) - m_b * l_s * g, 
			c_p * (y[1] - y[2]) - m_b * l_g * g])
			
	elif sw[1]: #state 2
		gp[0, 0] = 0
		gp[1, 0] = 1
		gp[0, 1] = h_s
		gp[1, 1] = r_s
		gp[2, 0] = 0
		gp[2, 1] = 0
	
		gyy[0] = yd[1]
		gyy[1] = yd[0] + r_s * yd[1]
		
		#Applied forces (f matrix)
		ff = np.array([-g * (m_s + m_b), 
			c_p * (y[2] - y[1]) - m_b * l_s * g, 
			c_p * (y[1] - y[2]) - m_b * l_g * g])
			
		#?
		yn[1] = (r0 - r_s)/h_s
		ydn[1] = 0
		
	else: #state 3
		gp[0, 0] = 0
		gp[1, 0] = 1
		gp[0, 1] = -h_s
		gp[1, 1] = r_s
		gp[2, 0] = 0
		gp[2, 1] = 0
	
		gyy[0] = yd[1]
		gyy[1] = yd[0] + r_s * yd[1]
		
		#Applied forces (f matrix)
		ff = np.array([-g * (m_s + m_b), 
			c_p * (y[2] - y[1]) - m_b * l_s * g, 
			c_p * (y[1] - y[2]) - m_b * l_g * g])
		
		#?
		yn[1] = (r_s - r0)/h_s
		ydn[1] = 0
			
	#? compute residual
	A = np.bmat([m, ydd], [-1*ff, gp.T])
	b = zeros(A)
	res_1 = yn
	res_2 = ydn
	res_3 = solve(A,b)
	
	#compute residual
	A = np.bmat([[m, gp.T], [gp, zeros([2,2])])
	b = np.bmat([[ff[]])