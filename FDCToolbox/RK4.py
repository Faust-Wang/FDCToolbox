from numpy import array

# Ref [1] page 179-180
def RK4(f, time, dt, xx, u, xcg):
	"""
	fourth-order Runge-Kutta Algorithm
	"""
	# k1
	xd = f(xx, u, Xcg=xcg)
	xa = xd*dt
	# k2
	x = xx + 0.5*xa
	t = time + 0.5*dt
	xd = f(x, u, Xcg=xcg)
	q = xd*dt
	# k3
	x = xx + 0.5*q
	xa = xa + 2.0*q
	xd = f(x, u, Xcg=xcg)
	q = xd*dt
	# k4
	x = xx + q
	xa = xa + 2.0*q
	time = time + dt
	xd = f(x, u, Xcg=xcg)
	
	# x_new	
	xnew = xx + (xa + xd*dt)/6.0
	return xnew

