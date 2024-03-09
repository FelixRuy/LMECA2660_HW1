import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import os


def read_vector_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    y = []
    x = []
    t = []
    for line in lines:
        parts = line.split(":")
        time = float(parts[0].split("=")[1].strip())
        values = [float(val.strip()) for val in parts[1].split(",")[:-1]]
        if time < 0:
            x = np.array(values)
        else:
            y.append(values)
            t.append(time)
    y = np.array(y)
    return x, y, np.array(t)

# Test the function
try:
    filename = "../data/comp_sol_phys.txt" 
    x_true, y_true, time_true = read_vector_file("../data/analytical_solution_phys.txt")
except FileNotFoundError:
    filename = "../data/comp_sol.txt"
    x_true, y_true, time_true = read_vector_file("../data/analytical_solution.txt")
x, y, time = read_vector_file(filename)
print(np.shape(x), np.shape(y), np.shape(time))
print(np.shape(x_true), np.shape(y_true), np.shape(time_true))

# Initialize the plot
fig, ax = plt.subplots()
line1, = ax.plot([], [], lw=2, marker='o', markersize=4, color='seagreen', label=r'$u^h(x)$')  # First line
line2, = ax.plot([], [], lw=2, color='gray', linestyle="-.", label='Analytical Solution')  # Second line

# Function to initialize the plot
def init():
    ax.set_xlim(-0.5, 0.5)  # Extend x-axis limit slightly to accommodate markers
    ax.set_ylim(-2, 2)
    ax.set_xlabel(r'$\frac{x}{L}$')
    ax.set_ylabel(r'$\frac{u}{U}$')
    ax.set_title(r'Convection over $\frac{ct}{L} \in [0, 1]$')
    ax.legend(loc='upper right')
    ax.grid(True)  # Adding grid
    return line1, line2

# Function to update the plot at each frame
def animate(i):
    line1.set_data(x, y[i])  # Update first line
    line2.set_data(x_true, y_true[i])  # Update second line
    return line1, line2

# Initialize the animation
ani = animation.FuncAnimation(fig, animate, frames=len(y), init_func=init, interval=50, blit=True)
#ani.save('convect.gif', writer='imagemagick', fps=10)
plt.show()