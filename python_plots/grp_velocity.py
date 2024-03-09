import numpy as np
import matplotlib.pyplot as plt

U = 1 
L = 1 
sig = 1/16
h = sig/8
N = int(8/sig) 


kh = lambda j : 2*np.pi/L * j * h
c_E2 = lambda kh : np.cos(kh)
c_E4 = lambda kh : np.cos(kh)*(4/3 - 1/3*np.cos(kh))+1/3*np.sin(kh)**2
c_ED = lambda kh : 1/6*(-2*np.exp(-2*1j*kh) + np.exp(-1j*kh) + 2*np.exp(1j*kh))
c_I4 = lambda kh: 3/2 * (np.cos(kh) + 1/2 * (np.cos(kh)**2-np.sin(kh)**2)) / (1 + 1/2 * np.cos(kh))**2

j = np.linspace(-N/2, N/2, 100)

plt.figure(figsize=(7,5))
plt.plot(j, np.real(c_ED(kh(j))), color="seaGreen", linestyle="dotted", label=r'ED')
plt.plot(j, np.real(c_E2(kh(j))), color="blue", linestyle="dashdot", label=r'E2')
plt.plot(j, np.real(c_E4(kh(j))), color="gray", label=r'E4')
plt.plot(j, np.real(c_I4(kh(j))), color="darkorchid", linestyle="dashed", label=r'I4')
plt.xlabel(r'$j$') ; plt.ylabel(r'Amplitude Error : $\Re\{\frac{c^*_g}{c}\}$') ; plt.legend() ; plt.grid()
#plt.savefig(f'../graphs/grp_vel_re.pdf', bbox_inches='tight')
plt.show()


plt.figure(figsize=(7,5))
plt.plot(j, np.imag(c_ED(kh(j))), color="seaGreen", linestyle="dotted", label=r'ED')
plt.plot(j, np.imag(c_E2(kh(j))), color="blue", linestyle="dashdot", label=r'E2')
plt.plot(j, np.imag(c_E4(kh(j))), color="gray", label=r'E4')
plt.plot(j, np.imag(c_I4(kh(j))), color="darkorchid", linestyle="dashed", label=r'I4')
plt.xlabel(r'$j$') ; plt.ylabel(r'Phase Error : $\Im\{\frac{c^*_g}{c}\}$') ; plt.legend() ; plt.grid()
#plt.savefig(f'../graphs/group_vel_im.pdf', bbox_inches='tight')
plt.show()



