from numpy import sqrt

# Ref [1] page 634
# Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/adc.py
def adc(vt, alt):
    '''
    converts velocity (vt) and altitude (alt) to mach number (amach) and dynamic pressure (qbar)
    '''

    # vt = freestream air speed

    ro = 2.377e-3
    tfac = 1 - .703e-5 * alt

    if alt >= 35000: # in stratosphere
        t = 390
    else:
        t = 519 * tfac # 3 rankine per atmosphere (3 rankine per 1000 ft)

    # rho = freestream mass density
    rho = ro * tfac**4.14

    # a = speed of sound at the ambient conditions
    # speed of sound in a fluid is the sqrt of the quotient of the modulus of elasticity over the mass density
    a = sqrt(1.4 * 1716.3 * t)

    # amach = mach number
    amach = vt / a

    # qbar = dynamic pressure
    qbar = .5 * rho * vt * vt

    return amach, qbar

