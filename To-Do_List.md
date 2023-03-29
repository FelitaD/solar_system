# To-do list

##### 1 
implement a new read initial conditions function to set up the system 
to read particle positions from a file 
(You will need a function to open a file, read information from it, 
create a Particle3D object for each line of text in the file, 
and then return a list of all those objects. 
There are several ways to do this.
There are two examples of the kind of file the function should read included on 
Learn:solar system.txt and sun earth.txt. In each case their format is one body per line,
with each line taking the form:
Name mass x pos y pos z pos x vel y vel z vel
),
##### 2 
subtract the centre-of-mass velocity from a list of particles,
##### 3 
use the other two functions from the previous section to compute separations and forces,
##### 4 
update all the particle velocities at once,
##### 5 
update all the particle positions at once,
##### 6 
store the trajectory of the particles in a numpy array

****

`particle3D.py`

![img.png](data/img.png)

`simulate_particle.py`

![img_1.png](data/img_1.png)

`basic_functions_optimized.py`

![img_2.png](data/img_2.png)


