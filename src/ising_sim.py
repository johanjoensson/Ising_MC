import numpy as np
import ising
import lattice
import matplotlib.pyplot as plt

width = 500
height = 500
lat = lattice.Lattice_2D(np.array([[width, 0.], [0., height]]))
ising_sim = ising.Ising_2D(lat, [width, height], True)
ising_sim.set_beta(8)
#0sing_sim.set_H(0.01)
ising_sim.set_Jij([0.8, 0.0, -0.4])
ising_sim.add_spin_correlator(2, 3.)

num_iterations = width*height*100

corr = np.zeros((num_iterations, 3))
energy = np.zeros(num_iterations)
magnetization = np.zeros(num_iterations)

for it in range(num_iterations):
    ising_sim.update()
    r, _, si, sj = ising_sim.measure_spin_correlators()[0]
    corr[it] = np.array([r, si, sj])
    energy[it] = ising_sim.average_site_energy()
    magnetization[it] = ising_sim.magnetization()

#    if(it % (num_iterations/100) == 0 ):
#        print("Average site energy = " + repr(ising_sim.average_site_energy()))
#        print("Magnetization = " + repr(ising_sim.magnetization()))
#       plt.matshow(np.array(ising_sim.field).reshape((width, height)), cmap = "Blues_r")
#       plt.savefig("Ising.png")
#       plt.close()

plt.subplot(2, 1, 1)
plt.plot(energy)
plt.xlabel("iteration")
plt.ylabel("Average site energy")
plt.subplot(2, 1, 2)
plt.plot(magnetization)
plt.xlabel("iteration")
plt.ylabel("Magnetization")
plt.show()

