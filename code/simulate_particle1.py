"""
Symplectic Euler and Velocity Verlet time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.

Author: FILL IN YOUR NAME HERE
Number: FILL IN YOUR MATRICULATION NUMBER HERE,
        THE ONES STARTING WITH AN S.

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from particle1D import Particle1D

def force_double_well(p1, a, b):
    """
    Return the force on a particle in a double well potential.

    The force is given by
        F(x) = -dV/dx = -4*a*x^3 + 2*b*x

    Parameters
    ----------
    p1: Particle1D
    a: float
    b: float

    Returns
    -------
    force: float
    """
    force = -4 * a * p1.position**3 + 2 * b * p1.position
    return force


def potential_double_well(p1, a, b):
    """
    Method to return potential energy 
    of particle in double-well potential
    V(x) = a*x^4 - b*x^2

    Parameters
    -----------
    p1: Particle1D
    a: float
    b: float

    Returns
    -------
    potential: float
    """
    potential = a * p1.position**4 - b * p1.position**2
    return potential


def main():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect three things:
    #    the name of this file
    #    euler or verlet
    #    the name of the output file the user wants to write to
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) != 3 :
        print("You left out the name of the output file when running.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <euler or verlet> <desired output file>")
        sys.exit(1)
    else:
        mode = sys.argv[1]
        outfile_name = sys.argv[2]

    # Open the output file for writing ("w")
    outfile = open(outfile_name, "w")

    # Choose our simulation parameters
    dt = 0.01
    numstep = 2000
    time = 0.0
    a = 0.1
    b = 1.0

    # Set up particle initial conditions:
    #  position x0 = 0.0
    #  velocity v0 = 1.0
    #  mass      m = 1.0
    p1 = Particle1D(0.0, 1.0, 1.0)


    # Get initial force
    force = force_double_well(p1, a, b)

    # Write out starting time, position, and energy values
    # to the output file.
    energy = p1.kinetic_energy() + potential_double_well(p1, a, b)
    outfile.write(f"{time}    {p1.position}    {energy}\n")

    # Initialise numpy arrays that we will plot later, to record
    # the trajectories of the particles.
    times = np.zeros(numstep)
    positions = np.zeros(numstep)
    energies = np.zeros(numstep)

    # Start the time integration loop
    for i in range(numstep):

        # Update the positions and velocities.
        # This will depend on whether we are doing an Euler or verlet integration
        if mode == "euler":
            # Update particle position
            p1.update_position_1st(dt)
            
            # Calculate force
            force = force_double_well(p1, a, b)

            # Update particle velocity 
            p1.update_velocity(dt, force)

        elif mode == "verlet":
            # Update particle position using previous force
            p1.update_position_2nd(dt, force)
            
            # Get the force value for the new positions
            force_new = force_double_well(p1, a, b)

            # Update particle velocity by averaging
            # current and new forces
            p1.update_velocity(dt, 0.5*(force+force_new))
            
            # Re-define force value for the next iteration
            force = force_new
        else:
            raise ValueError(f"Unknown mode {mode} - should be euler or verlet")

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + potential_double_well(p1, a, b)
        outfile.write(f"{time} {p1.position} {energy}\n")

        # Store the things we want to plot later in our arrays
        times[i] = time
        positions[i] = p1.position
        energies[i] = energy

    # Now the simulation has finished we can close our output file
    outfile.close()

    # Plot particle trajectory to screen. There are no units
    # here because it is an abstract simulation, but you should
    # include units in your plot labels!
    pyplot.figure()
    pyplot.title('Symplectic Euler: position vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Position')
    pyplot.plot(times, positions)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.figure()
    pyplot.title('Symplectic Euler: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(times, energies)
    pyplot.show()


# This strange but standard python idiom means that the main function will
# only be run if we run this file, not if we just import it from another
# python file. It is good practice to include it whenever your code can be
# run as a program.
if __name__ == "__main__":
    main()

