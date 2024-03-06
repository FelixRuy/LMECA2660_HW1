import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# Function to generate random data
def generate_data():
    x = np.linspace(0, 10, 100)
    y = np.random.rand(100)
    return x, y

# Initialize the plot
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)

# Function to initialize the plot
def init():
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    return line,

# Function to update the plot at each frame
def update(frame):
    x, y = generate_data()
    line.set_data(x, y)
    return line,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=50, init_func=init, blit=True)

plt.show()
