from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import fire


def main(a, ecc, i=None, lon_ascend=None, arg_peri=None):
    """Calculate radius and velocity given orbital elements.

    Parameters
    ----------
    a : np.float
        semi-major axis
    ecc : np.float
        eccentricity
    i : np.float
        inclination
    lon_ascend : np.float
        longitude of the ascending node
    arg_peri : np.float
        argument of the periapsis
    """
    mu = 3.986012e5  # [km^3 s^-2]  gravitational parameter
    p = a*(1 - ecc**2)
    E = -mu/(2*a)

    rp = p/(1+ecc)  # radius at periapsis
    vp = np.sqrt(2*(E+mu/rp))

    print(vp, rp)


if __name__ == "__main__":
   fire.Fire(main)
