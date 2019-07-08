import numpy as np
import ising
import lattice
import matplotlib.pyplot as plt

width = 500
height = 500
lat = lattice.Lattice_2D(np.array([[1., 0], [0, 1.]]))
ising_sim = ising.Ising_2D(lat, [width, height], True)
ising_sim.set_beta(1)
#0sing_sim.set_H(0.01)
ising_sim.set_Jij([1, 0.0, -0.2])

num_iterations = 100000000

for it in range(num_iterations):
    ising_sim.update()
    if(it % (num_iterations/100) == 0 ):
        print("Average site energy = " + repr(ising_sim.average_site_energy()))
        print("Magnetization = " + repr(ising_sim.magnetization()))
        plt.matshow(np.array(ising_sim.field).reshape((width, height)), cmap = "Blues_r")
        plt.savefig("Ising.png")
        plt.close()

