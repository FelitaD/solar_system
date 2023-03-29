"""
CompMod SCBF: This first part of the project is two basic functions that will form
part of the rest of the simulation code.

The first function, compute separations(particles), will compute the vectors be
-tween particles, for later use in calculating forces.
The second function, compute forces potential(particles, separations), will com
-pute the gravitational forces on the particles and the total gravitational potential
of the system.

Author: Charlie Fraser
Number: s1728085
Date: 16/2/23
"""
import numpy as np


def compute_separations(particles):
    """
    Creates a loop to calculate distance of the particles by their x, y, z components up to as many fed in.
    
    Assigns them to create n 3D arrays containing n null lines for a separation of a particle from itself.
    
    Parameters
    ----------
    particles: list
        Particles in the system with attributes of label, mass, position and velocity.
    
    Returns
    ----------
    separations: 3d array
        Vector position separation of each pair of particles.
    """
    n = len(particles)
    separations = np.zeros((n, n, 3))
    for i in range(n):
        for j in range(n):
            if j != i:
                separation = particles[i].position - particles[j].position
                separations[i, j] = separation
                separations[j, i] = -separation
    return separations


def compute_forces_potential(particles, separations):
    """
    Creates a loop to sum force acted upon particles by their x, y, z components up to as many fed in with a second loop to sum potential using said particles.
    
    Assigns them to n lines of a 3D array with the potential alongisde.
    
    Parameters
    ----------
    particles: list
        Particles in the system with attributes of label, mass, position and velocity.
    
    separations: 3d array
        Vector position separation of each pair of particles.
    
    Returns
    ----------
    forces: 2d array
        Force on each particle
    
    potential: float
        Total potential in the system
    """
    G = 8.887602593e-10  # Gravitational constant in units: AU^3 Mearth^-1 day^-2
    n = len(particles)
    forces = np.zeros((n, 3))
    for i in range(n):
        for j in range(n):
            if j > i:
                mod_separations = np.sqrt((separations[i, j, 0] ** 2) + (separations[i, j, 1] ** 2) + (
                            separations[i, j, 2] ** 2))  # For the denominator for universal gravitation equation
                force_normal = -G * particles[i].mass * np.sum(
                    particles[j].mass * separations[i, j] / (mod_separations ** 3))  # F = -Gmi∑jmjri-rj/|ri-rj|^3
                forces[i] += force_normal  # Using Broadcasting
                forces[j] -= force_normal
    for i in range(n):
        for j in range(n):
            if j > i:
                potential = np.sum(
                    -G * particles[i].mass * particles[j].mass / mod_separations)  # U = ∑i∑j>i-Gmimj/|ri-rj|
    return forces, potential
