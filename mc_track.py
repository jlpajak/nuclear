import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D
from mc import *

if __name__ == "__main__":
    h1 = Nuclide(1, 1)
    o16 = Nuclide(16, 8)

    h1.add_cross_section('capture', 'h1_n_g.txt')
    h1.add_cross_section('elastic', 'h1_n_el.txt')
    o16.add_cross_section('capture', 'o16_n_g.txt')
    o16.add_cross_section('elastic', 'o16_n_el.txt')

    water = Material(1.0, "water")
    water.add_nuclide(2, h1)
    water.add_nuclide(1, o16) 

    neutron = Neutron(maxwell_energy(1.29), [0, 0, 0], None)
    
    position = []
    position.append(copy.deepcopy(neutron.position))
    while neutron.E > 25e-9:
        S_t =  water.macro_cross_section(neutron.E)
        r = water.free_path(1 / S_t)
        neutron.move(r)
        reaction, nuclide = water.random_event(neutron.E)
        position.append(copy.deepcopy(neutron.position))
        if reaction == 'elastic':
            neutron.scattering(nuclide)
        else:
            break
    position = numpy.array(position)
    fig = plt.figure(1, (6, 6))
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    ax1.plot(position[:, 0], position[:, 1], position[:, 2], 'o-', ms=3)
    plt.show()



