import random
import numpy
import matplotlib.pyplot as plt
import math

def maxwell_energy(kT):
    while True:
        x = random.uniform(0,10)
        y = random.uniform(0,0.40)
        N = 2/(numpy.pi)**1/2 * (x/(kT)**3)**1/2 * math.exp(-x/kT);
        if y <=N:
           return x


class Neutron:

    MASS = 939.565 # MeV/c2


    def __init__(self, E, position, momentum):
        self.E = E
        self.position = position
        if momentum is None:
            theta, phi = self.isotropic_direction()
            self.momentum = [theta, phi]
        else:
            self.momentum = momentum

    def isotropic_direction (self):
        cos_teta = random.uniform(-1 , 1)
        fi = random.uniform(0, 2 * numpy.pi)
        teta = numpy.arccos(cos_teta)
        return teta, fi

    def scattering(self, nuclide):
        
        thetaCM, fiCM = Neutron.isotropic_direction(nuclide)
        
        self.E *= ((nuclide.A**2 + 1 + 2 * nuclide.A * numpy.cos(thetaCM))/
                ((nuclide.A + 1)**2))
        
        thetaCM = ( (1 + nuclide.A * numpy.cos(thetaCM)) / 
            ((1 + nuclide.A**2 + 2 * nuclide.A * numpy.cos(thetaCM))**(1/2)))
        
        fiCM = ((1 + nuclide.A * numpy.cos(fiCM)) / 
                ((1 + nuclide.A**2 + 2 * nuclide.A * 
                    numpy.cos(fiCM))**(1/2)))
        
        thetaCM = numpy.arccos(thetaCM)
        fiCM = numpy.arccos(fiCM)
        self.momentum = [thetaCM, fiCM]
        
        
        

    def move(self, r):
        self.position[0] += r * numpy.cos (self.momentum[0]) * numpy.cos (self.momentum[1])
        self.position[1] += r * numpy.cos (self.momentum[0]) * numpy.sin (self.momentum[1])
        self.position[2] += r * numpy.sin (self.momentum[0])

class Nuclide:

    def __init__(self, A, Z):
        self.A = A
        self.Z = Z
        self.cross_sections = {}


    def add_cross_section(self, reaction, file_name):
        data_file = open(file_name)

        data = []
        for line in data_file:
            if line.startswith('#'):
                continue
            Ei = float(line[0:13])
            Ci = float(line[13:25])
            data.append([Ei, Ci])

        data = numpy.array(data)
        self.cross_sections[reaction] = data


    def get_cross_section(self, reaction, E_neutron):
        data = self.cross_sections[reaction]
        E = data[:, 0]
        C = data[:, 1]

        j = 0
        k = 0
        for i, Ei in enumerate(E):
            if Ei > E_neutron:
                j = i - 1
                k = i
                break
        else:
            return C[-1]

        if k == 0:
            return C[0]
        else:
            a = (C[k] - C[j]) / (E[k] - E[j])
            b = C[j] - a * E[j]
            return a * E_neutron + b


class Material:
    def __init__(self, density, name):
        self.density = density
        self.name = name
        self.nuclides = []

    def add_nuclide(self, quantity, nuclide):
        self.nuclides.append([quantity, nuclide])

    def macro_cross_section(self, E_neutron):
        Na_times_barn = 0.602
        molar_mass = 0
        sum_C = 0
        for q, nuclide in self.nuclides:
            molar_mass += q * nuclide.A
            for reaction in nuclide.cross_sections.keys():
                sum_C += q * nuclide.get_cross_section(reaction, E_neutron)
        return self.density * Na_times_barn / molar_mass * sum_C

    def random_event(self, E_neutron):
        cs_total = 0
        process = []
        for nuclide in self.nuclides:
            for reaction in nuclide[1].cross_sections.keys():
                cs_total += (nuclide[0] * 
                        nuclide[1].get_cross_section(reaction, E_neutron))
                process.append([cs_total, reaction, nuclide[1]])

        r = random.random()
        for p in process:
            if p[0] / cs_total >= r:
                return p[1], p[2]
        else:
            return process[-1][1], process[1][2]


    def free_path(self, mean_free_path):
        r = random.random()    
        return - mean_free_path*numpy.log(1 - r)


if __name__ == '__main__':
    u235 = Nuclide(235, 92)
    u235.add_cross_section('fission', 'u235_n_f.txt')
    u235.add_cross_section('capture', 'u235_n_g.txt')
    u235.add_cross_section('elastic', 'u235_n_el.txt')
    print(u235.get_cross_section('fission', 1.0))


