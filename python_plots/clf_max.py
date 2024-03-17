import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the functions
ED = lambda kh, CFL: -1/6 * (np.exp(-2j*kh) - 6*np.exp(-1j*kh) + 3 + 2*np.exp(1j*kh)) * CFL
RK4 = lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24

# Generate CFL and kh values
CFL_range = np.linspace(1, 2, 100)
kh = np.linspace(0, np.pi, 100)
CFL, kh = np.meshgrid(CFL_range, kh)

# Compute the values for the surface plot
z_values = np.abs(RK4(ED(kh, CFL)))

# Create the figure and 3D axis
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
surf = ax.plot_surface(CFL, kh, z_values, cmap='Greys_r', alpha=0.7)

# Highlight points where the value is more than 1
x_highlighted, y_highlighted = np.where(z_values > 1)

print(min(y_highlighted))
print(CFL_range[min(y_highlighted)])

ax.scatter(CFL[x_highlighted, y_highlighted], kh[x_highlighted, y_highlighted], z_values[x_highlighted, y_highlighted], color='r')
ax.view_init(25, 161)


# Add labels
ax.set_xlabel('CFL')
ax.set_ylabel('kh')
ax.set_zlabel(r'$|\rho|$')

plt.savefig("../graphs/maximum_clf_surf.pdf", bbox_inches='tight', pad_inches=0.3)
plt.show()
