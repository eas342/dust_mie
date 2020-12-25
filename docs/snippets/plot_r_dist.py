from dust_mie import calc_mie
import matplotlib.pyplot as plt

median_r = 1.0
s = 0.5

r, dr = calc_mie.get_r_to_evaluate(r=median_r,s=s)

n = calc_mie.lognorm(r,s=s,med=median_r)

plt.plot(r,n)
plt.xlabel('Particle Radius ($\mu$m)')
plt.ylabel('Number')
plt.savefig('radius_distribution.png')