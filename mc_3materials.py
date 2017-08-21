import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D
from mc import *

if __name__ == "__main__":
    h1 = Nuclide(1, 1)
    d2 = Nuclide(2, 1)
    c12 = Nuclide(12, 6)
    o16 = Nuclide(16, 8)

    h1.add_cross_section('capture', 'h1_n_g.txt')
    h1.add_cross_section('elastic', 'h1_n_el.txt')
    d2.add_cross_section('capture', 'h2_n_g.txt')
    d2.add_cross_section('elastic', 'h2_n_el.txt')
    c12.add_cross_section('capture', 'c12_n_g.txt')
    c12.add_cross_section('elastic', 'c12_n_el.txt')
    o16.add_cross_section('capture', 'o16_n_g.txt')
    o16.add_cross_section('elastic', 'o16_n_el.txt')

    water = Material(1.0, "water")
    water.add_nuclide(2, h1)
    water.add_nuclide(1, o16)
    
    heavywater = Material(1.1, "d20")
    heavywater.add_nuclide(2, d2)
    heavywater.add_nuclide(1, o16) 

    graphite = Material(2.03, "graphite")
    graphite.add_nuclide(1, c12)
    
    iteracjawoda = 0.0
    zderzeniawoda = []

    while iteracjawoda < 1000:
        neutron = Neutron(maxwell_energy(1.29), [0, 0, 0], None)
        nzderzen = 0.0
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
                nzderzen += 1
            else:
                break
        iteracjawoda += 1
        zderzeniawoda.append(nzderzen)
        
    iteracjaciezkawoda = 0.0
    zderzeniaciezkawoda = []

    while iteracjaciezkawoda < 1000:
        neutron = Neutron(maxwell_energy(1.29), [0, 0, 0], None)
        nzderzen = 0.0
        position = []
        position.append(copy.deepcopy(neutron.position))
        while neutron.E > 25e-9:
            S_t =  water.macro_cross_section(neutron.E)
            r = heavywater.free_path(1 / S_t)
            neutron.move(r)
            reaction, nuclide = heavywater.random_event(neutron.E)
            position.append(copy.deepcopy(neutron.position))
            if reaction == 'elastic':
                neutron.scattering(nuclide)
                nzderzen += 1
            else:
                break
        iteracjaciezkawoda += 1
        zderzeniaciezkawoda.append(nzderzen)

    iteracjagrafit = 0.0
    zderzeniagrafit = []

    while iteracjagrafit < 1000:
        neutron = Neutron(maxwell_energy(1.29), [0, 0, 0], None)
        nzderzen = 0.0
        position = []
        position.append(copy.deepcopy(neutron.position))
        while neutron.E > 25e-9:
            S_t =  graphite.macro_cross_section(neutron.E)
            r = graphite.free_path(1 / S_t)
            neutron.move(r)
            reaction, nuclide = graphite.random_event(neutron.E)
            position.append(copy.deepcopy(neutron.position))
            if reaction == 'elastic':
                neutron.scattering(nuclide)
                nzderzen += 1
            else:
                break
        iteracjagrafit += 1
        zderzeniagrafit.append(nzderzen)

    plt.subplot(2, 2, 1)
    plt.hist(zderzeniawoda, bins=50, normed=False, range=[0.0, 50.0], facecolor='green', alpha=0.75)
    plt.xlabel('liczba zderzen')
    plt.ylabel('liczba neutronow')
    plt.subplot(2, 2, 2)
    plt.hist(zderzeniaciezkawoda, bins=60, normed=False, range=[0.0, 60.0], facecolor='green', alpha=0.75)
    plt.xlabel('liczba zderzen')
    plt.ylabel('liczba neutronow')
    plt.subplot(2, 2, 3)
    plt.hist(zderzeniagrafit, bins=160, normed=False, range=[0.0, 160.0], facecolor='green', alpha=0.75)
    plt.xlabel('liczba zderzen')
    plt.ylabel('liczba neutronow')
    plt.tight_layout()
    plt.show()
