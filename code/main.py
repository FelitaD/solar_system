from particle3D import Particle3D
from simulate_particle import compute_verlet
from basic_functions_optimised import compute_separations, compute_forces_potential



def read_line(filename):
    """
    Creates a Particle3D instance given a line of text.

    The input line should be in the format:
    label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

    Parameters
    ----------
    filename: str
        Readable file handle in the above format

    Returns
    -------
    p: Particle3D
    """
    # Split the input line into individual items
    with open(filename, 'r') as f:
        data = f.read()

    items = data.split('\n')

    list_items = []
    for string_item in items:
        element = string_item.split(' ')
        list_items.append(element)

    particles = []
    # Extract the label, mass, position, and velocity from the input line
    for i in range(len(list_items) - 1):
        label = list_items[i][0]
        mass = float(list_items[i][1])
        position = list_items[i][2:5]
        velocity = list_items[i][5:]
        # Create a Particle3D object with the extracted items
        particles.append(Particle3D(label, mass, position, velocity))

    return particles


def main():
    sun, earth = read_line('data/sun_earth.txt')
    print(sun.label)
    print(sun.mass)
    print(sun.position)
    print(sun.velocity)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
