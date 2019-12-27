'''
This program numerically calculates the equilibrium state and control vectors of an F-16 model given
certain parameters.

states:                                              controls:
	x1 = Vt		x4 = phi	x7 = p	  x10 = pn			u1 = throttle
	x2 = alpha	x5 = theta	x8 = q	  x11 = pe			u2 = elevator
	x3 = beta	x6 = psi    x9 = r	  x12 = alt		    u3 = aileron
                                      x13 = pow         u4 = rudder
'''

import numpy as np
from numpy import sin, cos, arcsin, sqrt, tan, arctan

from scipy.optimize import fmin

from atmospheric_model import adc
from F16_Engine_data import tgear
from F16_Model import f16_model

# from https://github.com/AeroPython/PyFME/blob/master/src/pyfme/utils/trimmer.py
def turn_coord_cons(turn_rate, alpha, beta, TAS, gamma=0):
    """Calculates phi for coordinated turn.
    """

    g0 = 32.17
    G = turn_rate * TAS / g0

    if abs(gamma) < 1e-8:
        phi = G * cos(beta) / (cos(alpha) - G * sin(alpha) * sin(beta))
        phi = arctan(phi)
    else:
        a = 1 - G * tan(alpha) * sin(beta)
        b = sin(gamma) / cos(beta)
        c = 1 + G ** 2 * cos(beta) ** 2

        sq = sqrt(c * (1 - b ** 2) + G ** 2 * sin(beta) ** 2)

        num = (a - b ** 2) + b * tan(alpha) * sq
        den = a ** 2 - b ** 2 * (1 + c * tan(alpha) ** 2)

        phi = arctan(G * cos(beta) / cos(alpha) * num / den)
    return phi

# from https://github.com/AeroPython/PyFME/blob/master/src/pyfme/utils/trimmer.py
## ROC constraints
def rate_of_climb_cons(gamma, alpha, beta, phi):
    """Calculates theta for the given ROC, wind angles, and roll angle.
    """
    a = cos(alpha) * cos(beta)
    b = sin(phi) * sin(beta) + cos(phi) * sin(alpha) * cos(beta)
    sq = sqrt(a ** 2 - sin(gamma) ** 2 + b ** 2)
    theta = (a * b + sin(gamma) * sq) / (a ** 2 - sin(gamma) ** 2)
    theta = arctan(theta)
    return theta


# Ref [1] page 645
# Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/conf16.py
## constraints
def conf16(x, u, const):
    """
    Apply constraints to x variable
    used when finding trim conditions
    """

    radgam, singam, rr, pr, tr, phi, cphi, sphi, thetadot, coord, stab, orient = const
    gamm = arcsin(singam)

    #
    # Steady Level Flight
    #
    if orient == 1:
        x[3] = phi          # Phi
        x[4] = x[1]         # Theta
        x[6] = rr           # Roll Rate
        x[7] = pr           # Pitch Rate
        x[8] = 0.0          # Yaw Rate

    #
    # Steady Climb
    #
    if orient == 2:
        x[3] = phi          # Phi
        x[4] = x[1] + radgam  # Theta
        x[6] = rr           # Roll Rate
        x[7] = pr           # Pitch Rate
        x[8] = 0.0          # Yaw Rate

    #
    # orient=3 implies coordinated turn
    #
    if orient == 3: # tr = turn rate (page 342 book)
        x[3] = phi # turn_coord_cons(tr, x[1], x[2], x[0], gamma=gamm)
        x[4] = rate_of_climb_cons(gamm, x[1], x[2], phi) # theta
        x[6] = -tr * sin(x[4])            # Roll Rate
        x[7] = tr * cos(x[4]) * sin(x[3])  # Pitch Rate
        x[8] = tr * cos(x[4]) * cos(x[3])  # Yaw Rate

    #
    # Pitch Pull Up
    #
    if orient == 4:
        x[4] = x[1]         # Theta = alpha
        x[3] = phi          # Phi
        x[6] = rr           # Roll Rate
        x[7] = thetadot     # Pitch Rate
        x[8] = 0.0          # Yaw Rate

    x[12] = tgear(u[0])

    return x, u


# Ref [1] page 645
# Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/clf16.py
def clf16(s, x, u, xcg, const):
    '''
    This is the objective function for finding the trim condition of the initial states
    objective function of the optimization to find the trim conditions

    x and u get modified in-place
    returns the cost
    '''

    _, singam, _, _, tr, _, _, _, thetadot, _, _, orient = const
    gamm = arcsin(singam)

    if len(s) == 3:
        u[0] = s[0]
        u[1] = s[1]
        x[1] = s[2]
    else:
        x[1] = s[0]
        u[0] = s[1]
        u[1] = s[2]
        u[2] = s[3]
        u[3] = s[4]

    #
    # Get the current power and constraints
    #
    x[12] = tgear(u[0])
    x, u = conf16(x, u, const)

    # we just want the derivative
    subf16 = lambda x, u: f16_model(x, u, Xcg=xcg)
    
    xd = subf16(x, u)

    #
    # Steady Level flight
    #
    if orient == 1:
        r = 100.0*(xd[0]**2 + xd[1]**2 + xd[2]**2 + xd[6]**2 + xd[7]**2 + xd[8]**2)

    #
    # Steady Climb
    #
    if orient == 2:
        r = 500.0*(xd[11]-x[0]*sin(gamm))**2 + xd[0]**2 + 100.0*(xd[1]**2 + xd[2]**2) + \
            10.0*(xd[6]**2 + xd[7]**2 + xd[8]**2)



    #
    # Coord Turn
    #
    if orient == 3:
        r = xd[0]*xd[0] + 100.0 * (xd[1] * xd[1] + xd[2]*xd[2] + xd[11]*xd[11]) + 10.0*(xd[6]*xd[6] + \
            xd[7]*xd[7]+xd[8]*xd[8]) + 500.0*(xd[5] - tr)**2

    #
    # Pitch Pull Up
    #


    if orient == 4:
        r = 500.0*(xd[4]-thetadot)**2 + xd[0]**2 + 100.0*(xd[1]**2 + xd[2]**2) + 10.0*(xd[6]**2 + xd[7]**2 + xd[8]**2)

    #
    # Scale r if it is less than 1
    #
    if r < 1.0:
        r = r**0.5

    return r


# Ref [1] page 644
# Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/trim.py
def trim(orient, inputs, Xcg=0.35, printOn=False):
    """
    calculate equilibrium state'
    """

    Xguess = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    Uguess = np.array([0.0, 0.0, 0.0, 0.0])

    x = Xguess.copy()
    u = Uguess.copy()

    xcg = Xcg

    if printOn:
        print('------------------------------------------------------------')
        print('Running trimmerFun.py')

    # gamma singam rr  pr   tr  phi cphi sphi thetadot coord stab  orient
    const = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1]
    rtod = 57.29577951

    # orient: 'Wings Level (gamma = 0)','Wings Level (gamma <> 0)','Steady Turn','Steady Pull Up'
    const[11] = orient

    # inputs: [Vt, h, gamm, psidot, thetadot]
    x[0] = inputs[0]
    x[11] = inputs[1]

    if orient == 2:
        gamm = inputs[2]
        const[0] = gamm/rtod
        const[1] = sin(const[0])

    elif orient == 3:
        psidot = inputs[3]
        const[4] = psidot/rtod  # tr = turn 
        tr = const[4]

        gamm = inputs[2]
        const[0] = gamm

        phi = turn_coord_cons(tr, x[1], x[2], x[0], gamma=gamm)
        const[5] = phi
        const[6] = sin(phi)
        const[7] = cos(phi)
        x[4] = rate_of_climb_cons(gamm, x[1], x[2], phi)

    elif orient == 4:
        thetadot = inputs[4]
        const[8] = thetadot/rtod

    if orient == 3:
        s = np.zeros(shape=(5,))
        s[0] = x[1]
        s[1] = u[0]
        s[2] = u[1]
        s[3] = u[2]
        s[4] = u[3]
    else:               # for orient 1, 2, 4
        s = np.zeros(shape=(3,))
        s[0] = u[0]
        s[1] = u[1]
        s[2] = x[1]

    if printOn:
        print(f"initial cost = {clf16(s, x, u, xcg, const)}")


    ##=== FMIN Algorithm ================================================================
    s = fmin(clf16, s, args=(x, u, xcg, const), xtol=1e-12, maxiter=2000, disp=printOn)

    J = clf16(s, x, u, xcg, const)
    if printOn:
        print(f"cost = {J}")

    if orient != 3:
        x[1] = s[2]
        u[0] = s[0]
        u[1] = s[1]
    else:
        x[1] = s[0]
        u[0] = s[1]
        u[1] = s[2]
        u[2] = s[3]
        u[3] = s[4]
    
    ##===================================================================================

    if printOn:
        print(f'Throttle (percent):            {u[0]}')
        print(f'Elevator (deg):                {u[1]}')
        print(f'Ailerons (deg):                {u[2]}')
        print(f'Rudder (deg):                  {u[3]}')
        print(f'Angle of Attack (deg):         {rtod*x[1]}')
        print(f'Sideslip Angle (deg):          {rtod*x[2]}')
        print(f'Pitch Angle (deg):             {rtod*x[4]}')
        print(f'Bank Angle (deg):              {rtod*x[3]}')

        amach, qbar = adc(x[0], x[11])
        print(f'Dynamic Pressure (psf):        {qbar}')
        print(f'Mach Number:                   {amach}')

    return x, u

