import numpy as np
import matplotlib.pyplot as plt

CLF = 0.5 # CLF condition ?? 
# Eigenvalues (lamda * h/c * CLF) for partially decentered scheme (ED) 
l_ED = lambda kh : -1/6 * (np.exp(-2j*kh)-6*np.exp(-1j*kh)+3+2*np.exp(1j*kh)) * CLF 
kh = np.linspace(0, 2*np.pi, 5)
kh_cont = np.linspace(0, 2*np.pi, 100)
lamb = l_ED(kh_cont)
# Eigenvalues (lamda * h/c * CLF) for centered E2 scheme
l_E2 = lambda kh : (np.exp(1j*kh)-np.exp(-1j*kh)) * CLF
lamb_E2 = l_E2(kh)
# Eigenvalues (lamda * h/c * CLF) for centered E4 scheme
l_E4 = lambda kh : -1 * (2/3*(np.exp(1j*kh)-np.exp(-1j*kh)) - 1/12*(np.exp(2j*kh)-np.exp(-2j*kh))) * CLF
lamb_E4 = l_E4(kh)
# Eigenvalues (lamda * h/c * CLF) for implicit I4 scheme
l_I4 = lambda kh : -1j*((3/2*np.sin(kh)) / (1+1/2*np.cos(kh))) * CLF
lamb_I4 = l_I4(kh)
# Stability region for runge-kutta 4th order
l_RK4 = lambda z : 1 + z + z**2/2 + z**3/6 + z**4/24
x, y = np.meshgrid(np.arange(-3, 0.5, 0.01), np.arange(-3.5, 3.5, 0.01))
z = x + 1j * y
R = l_RK4(z)
zlevel = np.abs(R)

plt.plot(np.real(lamb), np.imag(lamb), color='gray', linestyle='dotted', linewidth=1, label=r'$ \frac{\lambda_k h}{c} * CLF$, ED')
plt.plot(np.real(lamb_E2), np.imag(lamb_E2), marker='x', color='black', linestyle='dotted', linewidth=1, label=r'$ \frac{\lambda_k h}{c} * CLF$, E2')
plt.plot(np.real(lamb_E4), np.imag(lamb_E4), marker='x', color='blue',linestyle='dashed', linewidth=1, label=r'$ \frac{\lambda_k h}{c} * CLF$, E4')
plt.plot(np.real(lamb_I4), np.imag(lamb_I4), marker='x', color='red', linestyle='dotted', linewidth=1, label=r'$ \frac{\lambda_k h}{c} * CLF$, I4')
plt.contourf(x, y, zlevel, levels=[0, 1], colors=[(0, 1, 0, 0.3)])
plt.grid()
plt.legend(loc='upper left')
plt.xlabel(r'$\Re\{\lambda \Delta t\}$')
plt.ylabel(r'$\Im\{\lambda \Delta t\}$')
plt.savefig(f'../graphs/CFL_{CLF}.pdf', bbox_inches='tight')
plt.show()




