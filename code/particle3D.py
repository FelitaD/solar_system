"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Charlie Fraser
Number: s1728085

"""
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in a3D space

    Attributes
    ----------
    label: name of the particle
    mass: mass of the particle
    position: position of the particle
    velocity: velocity of the particle

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
        """
        self.label = str(label)
        self.mass = float(mass)
        self.position = np.array(position)
        self.velocity = np.array(velocity)

    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        str
        """
        string = (self.label + ' ' + str(self.position[0]) + ' ' + str(self.position[1]) + ' ' + str(self.position[2]))
        return string

    def kinetic_energy(self):
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        ke: float
            1/2 m v**2
        """
        ke = 0.5 * self.mass * np.dot(self.velocity, self.velocity)
        return ke

    def momentum(self):
        """
        mom: float
            m v
        """
        mom = self.mass * self.velocity
        return mom

    def update_position_1st(self, dt):
        """
        A method update position 1st for a first-order update of the particle position
        for a timestep
        """
        self.position += dt * self.velocity
        ...

    def update_position_2nd(self, dt, force):
        """
        A method update position 2nd for a second-order update of the particle position
        for given timestep and force
        """
        self.position += dt * self.velocity + (dt ** 2) * force / 2 * self.mass
        ...

    def update_velocity(self, dt, force):
        """
        A method update velocity to update the velocity of a particle from a current
        time t to the next time t + δt, given a time-step dt and force f(t)
        """
        self.velocity += dt * force / self.mass
        ...


    @staticmethod
    def total_kinetic_energy(particles):
        """
        Computes the kinetic energy of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        total_kinetic_energy: float
        """
        total = 0
        for particle in particles:
            total += particle.kinetic_energy()
        return total

    @staticmethod
    def com_velocity(particles):
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity
        """
        # com = alice m * alice v + bob m * bob v + charlie m * charlie v / 
        top = 0
        bot = 0
        for i in range(len(particles)):
            top += particles[i].mass * particles[i].velocity
            bot += particles[i].mass
            com = top / bot
        return com

