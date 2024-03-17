import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib.table import Table
import re

def read_file_error(filename):
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
    filename = "data/comp_sol_phys.txt" 
    x_true, y_true, time_true = read_vector_file("data/analytical_solution_phys.txt")
except FileNotFoundError:
    filename = "data/comp_sol.txt"
    x_true, y_true, time_true = read_vector_file("data/analytical_solution.txt")
x, y, time = read_vector_file(filename)

# Initialize the plot
fig, axs = plt.subplots(2, 2)
fig.set_size_inches(12, 10)

# PLOT 0 1
time1 = 0.25
time2 = 0.75

N = len(x)
U05 =  np.array([*y[time==time1][0], *y[time==time1][0]])
U05_t = np.array([*y_true[time_true==time1][0], *y_true[time_true==time1][0]])
U1 =  np.array([*y[time==time2][0], *y[time==time2][0]])
U1_t = np.array([*y_true[time_true==time2][0], *y_true[time_true==time2][0]])
X = np.linspace(-1/2, 3/2, 2*N)
axs[0, 1].plot(X, U1, color="seagreen", label=r"$\frac{ct}{L} = $ "+str(time2))
axs[0, 1].plot(X, U05, color="blue", label=r"$\frac{ct}{L} = $ "+str(time1))
axs[0, 1].plot(X, U1_t, color="gray", linestyle="dotted")#, label=r"$\frac{ct}{L} = 1$ (analytical)")
axs[0, 1].plot(X, U05_t, color="gray", linestyle="dashed")#, label=r"$\frac{ct}{L} = 0.5$ (analytical)")

axs[0, 1].legend(loc='upper left', ncol=2)
axs[0, 1].set_xlabel(r"$\frac{x}{L}$")
axs[0, 1].set_ylabel(r"$\frac{u}{U}$")
axs[0, 1].set_title(r'Results at fixed times')
axs[0, 1].grid()


# PLOT 1 1
# Table: Information
scheme = ""; ic = ""; grid = ""; res = ""; cfl = ""; time_exec = ""; mode = ""
with open("data/params.txt", 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if i == 0:
            cfl = line.split("=")[1].strip()
        elif i == 1:
            res = line.split("=")[1].strip()
        elif i == 2:
            scheme = line.split("=")[1].strip()
        elif i == 3:
            mode = line.split("=")[1].strip()
            if mode == "0":
                mode = "Uniform"
            else:
                mode = "Non-uniform"
        elif i == 4:
            ic = line.split("=")[1].strip()
            if(ic == "0"):
                ic = "Gaussian"
            else:
                ic = "Wave-packet"
        elif i == 5:
            time_exec = line.split("=")[1].strip()[:7]+' s'


info_data = [['Scheme selected', f'{scheme}'],
             ['Initial conditions', f'{ic}'],
             ['Grid type', f'{mode}'],
             [r'Resolution $\frac{h}{\sigma}$', f'{res}'],
             ['CFL', f'{cfl}'],
             ['Time of computation', f'{time_exec}'],]

table = axs[1, 1].table(cellText=info_data, loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)  # Adjust size of the table
axs[1, 1].axis('off')  # Turn off the axis for the table



# PLOT 1 0
t, Ih, Eh, Rh = read_file_error('data/diag.txt')

axs[1, 0].plot(t, Rh, label=scheme, linestyle="dashdot", color="blue")
axs[1, 0].set_xlabel(r'$\frac{ct}{L}$')
axs[1, 0].set_ylabel(r'$R_h^{0.5}$')
axs[1, 0].set_title(r'Global Error Squared over time')
axs[1, 0].grid()
axs[1, 0].legend()


# PLOT 0 0

line1, = axs[0,0].plot([], [], lw=2, marker='o', markersize=4, color='seagreen', label=r'$u^h(x)$')  # First line
line2, = axs[0,0].plot([], [], lw=2, color='gray', linestyle="-.", label='Analytical Solution')  # Second line

# Function to initialize the plot
def init():
    axs[0, 0].set_xlim(-0.5, 0.5)  # Extend x-axis limit slightly to accommodate markers
    axs[0, 0].set_ylim(-1, 1.5)   
    axs[0, 0].set_xlabel(r'$\frac{x}{L}$')
    axs[0, 0].set_ylabel(r'$\frac{u}{U}$')
    axs[0, 0].set_title(r'Convection over $\frac{ct}{L} \in [0, 1]$')
    axs[0, 0].legend(loc='upper right', ncol=2)
    axs[0, 0].grid(True)  # Adding grid
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