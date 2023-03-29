"""
Velocity Verlet time integration of
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
from particle3D import Particle3D


def force_morse(p1, p2, alpha, D_e, r_e):
    """
    Return the force on a particle in a double well potential.

    The force is given by
        F(x_1, x_2) = 2*alpha*D_e*[1 - e^(-alpha(|x_ij| - r_e))]e^(-alpha*|(x_ij| - x_e))x^_ij

        where x_ij = x_2 - x_1
        and x^_ij - x_ij / |x_ij|

    Parameters
    ----------
    p1: Particle1D
    a: float
    b: float

    Returns
    -------
    force: float
    """
    x_ij = p2.position - p1.position
    x_ij_norm = np.linalg.norm(x_ij)
    x_ij_hat = x_ij / x_ij_norm

    force = 2 * alpha * D_e * (1 - np.exp(-alpha * (x_ij_norm - r_e))) * np.exp(-alpha * (x_ij_norm - x_e)) * x_ij_hat
    return force


def potential_morse(p1, p2, alpha, D_e, r_e):
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
    potential = D_e * ((1 - np.exp(-alpha * (np.linalg.norm(p2.position - p1.poisition) - r_e))) ** 2 - 1)
    return potential


def peaks(y, x):
    """


    Parameters
    ----------
    y : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    peaks = find_peaks(y)
    peak_pos = x[peaks[0]]
    T = peak_pos[1] - peak_pos[0]
    T_seconds = T * 1.018050571e-14
    v = T_seconds ** -1
    v_bar = v / 3e8
    v_final = v_bar / 100
    return v_final


def compute_verlet():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect three things:
    #    the name of this file
    #    the name of the output file the user wants to write to
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) != 3:
        print("You left out the name of the output file when running.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <verlet> <desired output file>")
        sys.exit(1)
    else:
        infile_name = sys.argv[1]
        outfile_name = sys.argv[2]

    # Open the output file for writing ("w")
    outfile = open(outfile_name, "w")

    # Open the input file for reading ("r")
    infile = open(infile_name, "r")

    lines = infile.readlines()

    # Choose our simulation parameters
    dt = float(lines[1].split()[0])
    numstep = int(lines[1].split()[1])
    D_e = float(lines[3].split()[0])
    r_e = float(lines[3].split()[1])
    alpha = float(lines[3].split()[2])
    time = 0.0
    numstep = 60 * 60 * 24

    # Set up particle initial conditions:
    #  position x0 = 0.0
    #  velocity v0 = 1.0
    #  mass      m = 1.0
    p1 = Particle3D.read_line(lines[4])
    p2 = Particle3D.read_line(lines[5])

    # Get initial force
    force = force_morse(p1, p2, alpha, D_e, r_e)
    force_react = -force

    # Write out starting time, position, and energy values
    # to the output file.
    energy = Particle3D.total_kinetic_energy([p1, p2]) + potential_morse(p1, p2, alpha, D_e, r_e)
    outfile.write(f"{time}    {np.sqrt(np.sum((p2_position - p1_position) ** 2))}    {energy}\n")

    # Initialise numpy arrays that we will plot later, to record
    # the trajectories of the particles.
    times[i] = time
    separations[i] = np.zeros(numstep)
    energies[i] = energy

    # Start the time integration loop
    for i in range(numstep):
        # Update particle position using previous force
        p1.update_position_2nd(dt, force)
        p2.update_position_2nd(dt, force)

        # Get the force value for the new positions
        force_new_i = force_morse(p1, p2, alpha, D_e, r_e)
        force_new_j = -force_new_i

        # Update particle velocity by averaging
        # current and new forces
        p1.update_velocity(dt, 0.5 * (force_i + force_new_i))
        p2.update_velocity(dt, 0.5 * (force_j + force_new_j))

        # Re-define force value for the next iteration
        force_i = force_new_i
        force_j = force_new_j

        # Increase time
        time += dt

        # Output particle information
        energy = Particle3D.total_kinetic_energy([p1, p2]) + potential_morse(p1, p2, alpha, D_e, r_e)
        outfile.write(f"{time} {separation} {energy}\n")

        # Store the things we want to plot later in our arrays
        times[i] = time
        positions[i] = np.sqrt(np.sum((p2_position - p1_position) ** 2))
        energies[i] = energy

    # Now the simulation has finished we can close our output file
    outfile.close()

    # Plot particle trajectory to screen. There are no units
    # here because it is an abstract simulation, but you should
    # include units in your plot labels!
    pyplot.figure()
    pyplot.title('Velocity Verlet: separation vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Separation')
    pyplot.plot(times, separation)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.figure()
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(times, energies)
    pyplot.show()


# This strange but standard python idiom means that the main function will
# only be run if we run this file, not if we just import it from another
# python file. It is good practice to include it whenever your code can be
# run as a program.
if __name__ == "__main__":
    compute_verlet()
