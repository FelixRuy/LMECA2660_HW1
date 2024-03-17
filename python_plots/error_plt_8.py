import numpy as np
import matplotlib.pyplot as plt
import re

def read_file(filename):
    time_array = []
    Ih_array = []
    Eh_array = []
    Rh_array = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.match(r'\(time, Ih, Eh, Rh\) : (\d+\.\d+), (\d+\.\d+), (\d+\.\d+), (\d+\.\d+)', line)
            if match:
                time_array.append(float(match.group(1)))
                Ih_array.append(float(match.group(2)))
                Eh_array.append(float(match.group(3)))
                Rh_array.append(float(match.group(4)))

    return time_array, Ih_array, Eh_array, Rh_array


tE2, IhE2, EhE2, RhE2 = read_file('diag_8_E2.txt')
tE4, IhE4, EhE4, RhE4 = read_file('diag_8_E4.txt')
tED, IhED, EhED, RhED = read_file('diag_8_ED.txt')
tI4, IhI4, EhI4, RhI4 = read_file('diag_8_I4.txt')

plt.figure(figsize=(16, 5))

# plot for Ih
plt.subplot(1, 3, 1)
plt.plot(tE2, IhE2, label='E2', linestyle="dashdot", color="blue")
plt.plot(tE4, IhE4, label='E4', color="gray")
plt.plot(tED, IhED, label='ED', linestyle="dotted", color="seaGreen")
plt.plot(tI4, IhI4, label='I4', linestyle="dashed", color="darkorchid")
plt.title(r'$I_h^n$')
plt.grid()
plt.xlabel(r'$\frac{ct}{L}$')
plt.legend()

# plot for Eh
plt.subplot(1, 3, 2)
plt.plot(tE2, EhE2, label='E2', linestyle="dashdot", color="blue")
plt.plot(tE4, EhE4, label='E4', color="gray")
plt.plot(tED, EhED, label='ED', linestyle="dotted", color="seaGreen")
plt.plot(tI4, EhI4, label='I4', linestyle="dashed", color="darkorchid")
plt.title(r'$E_h^n$')
plt.grid()
plt.xlabel(r'$\frac{ct}{L}$')
plt.legend()

# plot for Rh
plt.subplot(1, 3, 3)
plt.semilogy(tE2, RhE2, label='E2', linestyle="dashdot", color="blue")
plt.semilogy(tE4, RhE4, label='E4', color="gray")
plt.semilogy(tED, RhED, label='ED', linestyle="dotted", color="seaGreen")
plt.semilogy(tI4, RhI4, label='I4', linestyle="dashed", color="darkorchid")
plt.grid()
plt.legend()
plt.xlabel(r'$\frac{ct}{L}$')
plt.title(r'$R_h^n$')

#plt.savefig('../graphs/diag_plot_RESOL8.pdf', bbox_inches='tight')
plt.show()


