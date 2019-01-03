import numpy as np
import matplotlib.pyplot as plt


"""For now, lets just create an Earth ceneterd model for designing an Earth orbiting
constellation. Eventually, it will be necessary to support coordinate transformations between
ECI and ECEF at a minimum, but I'm going to start with ECI for now."""

class body():
    """Base class for orbital body."""
    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius


class satellite(body):
    def __init__(self, state):
        """
        parameters
        ----------
        state : dict, all the orbital elements of the satellite
            a : np.float, semi-major axis
            ecc : np.float, eccentricity
            i : np.float, inclination
            lon_ascend: np.float, longitude of the ascending node
            arg_peri: np.float, argument of the periapsis
            T: np.float, time of periapsis passage
        """
        np.state = state


def propagate(satellite):
    """J2 propagation. Assuming ECI coordinates for all functions."""


def main():
    m_e = 5.98e24
    r_e = 6378136.3 # meterss
    G = 6.670408e-11  # m^3 kg^-1 s^-2

    dt = 100 # s
    tf = 10*60
    time = np.linspace(0,tf,tf/dt+1)
    nsteps = len(time)

    ns = 9  # number of states (R3, position, velocity, acceleration)
    state = np.zeros((ns,1,nsteps))
    state[:,0,0] = np.array([])

    for idx in range(1,nsteps-1):
        mag_r = np.linalg.norm(state[0:3,0,idx-1])
        A = G*m_e/mag_r



if __name__ == "__main__":
    main()

