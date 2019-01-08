import numpy as np


def propagate_spg4(i0, e0, omega0, phi0, ma0, n0):
    """Propagate using SPG4 method with the NORAD mean orbial elements
    
    Args:
        i0 (float): inclination at epoch
        e0 (float): eccentricity at epoch
        omega0 (float): argument of periapsis at epoch
        phi0 (float): longitude of ascending node at epoch
        M0 (float): mean anomaly at epoch
        n0 (float): mean motion at epoch

    Returns:
        float: r_ecef
        float: v_ecef
        float: time
    """
    aE = 6378.1363  # the equitorial radius of the earth [km]
    ke = 3.986004418e14  # [m^3/s^2] Earth's gravitational parameter
    k2 = 1/2*J2*aE**2  # factor from the second gravitational harmonic of the earth
    k4 = -3/8*J4*aE**4  # factor from the fourth gravitational harmonic of the earth
    A_30 = -J3*aE**3


    # common factors
    common_a = 3*np.cos(i0)**2-1
    common_b = 1-e0**2

    a1 = (ke/n0)**(2./3)
    delta1 = 1.5*k2/a1**2 * (common_a)/(common_b)**1.5
    a0 = a1*(1-1/3*delta1-delta1**2-134/81*delta1**3)
    delta0 = 1.5*k2/a0**2*(common_a)/(common_b)**1.5
    dd_n0 = n0/(1+delta0)
    dd_a0 = a/(1-delta0)

    if 98 <= rp <= 156:
        s = dd_a0*(1 - e0) - s + aE

    # Calculate Constants
    theta = np.cos(i0)
    epsilon = 1/(dd_a0-s)
    beta0 = common_b**0.5  # e0 = mean eccentricity at epoch

    eta = dd_a0*e0*epsilon

    common_c = (q0 - s)**4

    C2 = common_c * common_b**4*dd_n0*(1-eta**2)**(-7./2) * \
            (dd_a0*(1 + 1.5*eta**2 + 4*e0*eta + e0*eta**3) + \
            1.5*k2*epsilon/(1 - eta)*(-0.5 + 1.5*theta**2)*(8 + 24*eta**2 + 3*eta**4))

    B = 0.5*Cd*A/m  # ballistic coefficient of satellite
    C1 = B*C2

    C3 = (common_c * epsilon**5 * A_30 * dd_n0 * aE * np.sin(i0))/(k2*e0)

    C4 = 2*dd_n0*common_c*eta**4*dd_a0*beta0**2*(1-eta**2)**(-7./2)* \
            ((2*eta*(1+e0*eta) + 0.5*e0 + 0.5*eta**3) - 2*k2*epsilon/(dd_a0*(1-eta**2))* \
            (3*(1-3*theta**2)*(1+1.5*eta**2-2*e0*eta-0.5*e0*eta**3) + 3./4*(1-theta**2)*(2*eta**2 - e0*eta - e0*eta**3)*cos(2*omega0)))
    
    C5 = 2*common_c*epsilon**4*dd_a0*beta0**2*(1-eta**2)**(-7./2)*(1 + 11./4*eta*(eta+e0) + e0*eta**3)

    D2 = 4*dd_a0*epsilon*C1**2
    D3 = 4./3*dd_a0*epsilon**2*(17*dd_a0 + s)*C1**3
    D4 = 2./3*dd_a0*epsilon**3*(221*dd_a0 + 32*s)*C1**4


    # secular effects of atmospheric drag and gravitation
    Mdf = M0 + (1 + 3*k2*(-1 + 3*theta**2)/(2*dd_a0**2*beta0**3) + \
                3*k2*(13 - 78*theta**2 + 137*theta**4)/(16*dd_a0**4*beta0**7))*dd_n0*(t-t0)

    omgeadf = omega0 + (-3*k2*(1-5*theta**2)/(2*dd_a0**2*beta0**4) + \
                        3*k2**2*(7-114*theta**2 + 395*theta**4)/(16*dd_a0**4*beta0**8) + \
                        5*k4*(3 - 36*theta**2 + 49*theta**4)/(4*dd_a0**4*beta0**8))*dd_n0(t-t0)
    
    phidf = phi0 + (-3*k2*theta/(dd_a0**2*beta**4) + \
                    3*k2**2*(4*theta - 19*theta**3)/(2*dd_a0**4*beta0**8) + \
                    5*k4*theta*(3-6*theta**2)/(2*dd_a0**4*beta0**8))*dd_n0*(t-t0)
    
    domega = B*C3*np.cos(omega0)*(t-t0)

    dM = -2./3*common_c*B*epsilon**4*aE/(e0*eta)*((1+eta*np.cos(Mdf))**3 - (1+eta*np.cos(M0)**3)

    Mp = Mdf + domega *domega

    omega = omegadf - domega - dM
    
    phi = phidf - 21./2*dd_n0*k2*theta/(dd_a0**2*beta0**2)*C1*(t-t0)**2

    e = e0 - B*C4*(t-t0) - B*C5*(np.sin(Mp) - np.sin(M0))

    a = dd_a0*(1 - C1*(t-t0) - D2*(t - t0)**2 - D3*(t-t0)**3 - D4*(t-t0)**4)**2

    IL = Mp + omega + phi + dd_n0*(3./2*C1*(t-t0)**2 + (D2 + 2*C1**2)*(t-t0)**3 + \
            1./4*(3*D3 + 12*C1*D2 + 10*C1**3)*(t-t0)**4 + \
            1./5*(3*D4 + 12*C1*D3 + 6*D2**2 + 30*C1)
    
    beta = np.sqrt(1-e**2)

    n = ke/a**(3./2)

    # long-period periodic terms

    axN = e*np.coS(omega)
    IL_L = A_30*np.sin(i0)/(8*k2*a*beta**2)*(e*np.cos(omega)*(3 + 5*theta)/(1+theta))

    a_yNL = A_30*np.sin(i0)/(4*k2*a*beta**2)

    