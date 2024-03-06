import numpy as np 
import matplotlib.pyplot as plt

a = 3/5
L = 1

map = lambda xi : xi - a*L/(2*np.pi) * np.sin(2*np.pi*xi/L)

xi = np.linspace(-L/2, L/2, 32)
x = map(xi)

plot = plt.plot(xi, np.zeros(32), 'o', markersize=2)
plt.plot(x, np.zeros(32), 'o', markersize=2)
plt.show()