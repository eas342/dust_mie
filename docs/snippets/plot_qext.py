from dust_mie import calc_mie
import numpy as np
import matplotlib.pyplot as plt

wave = np.linspace(0.5,2.5,50)
qext, qsca, qback, g = calc_mie.get_mie_coeff(wave,r=1.0,material='SiO2')

plt.plot(wave,qext)
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Wavelength ($\mu$m)')
plt.savefig('extinct_func.png')