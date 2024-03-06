import numpy as np
import matplotlib.pyplot as plt



# Phase and amplification factor for partially decentered scheme (ED)

kh = np.linspace(0, np.pi, 100)
kmod_ED = lambda kh : 1/(6*1j) * (np.exp(-2*kh*1j) - 6*np.exp(-kh*1j) + 3 + 2*np.exp(kh*1j)) / np.pi
kmod_E2 = lambda kh : np.sin(kh) / np.pi
kmod_E4 = lambda kh : np.sin(kh)*(1+1/3*(1-np.cos(kh))) / np.pi
kmod_I4 = lambda kh : ( 22/14 * np.sin(kh) + 1/14 * np.sin(2*kh) ) / ( 1 + 10/14 * np.cos(kh) ) / np.pi

plt.figure(figsize=(6,5))
plt.plot(kh/np.pi, kh/np.pi,color="dimgray", label=r'$\Re\{\frac{kh}{\pi}\}$')
plt.plot(kh/np.pi, np.real(kmod_ED(kh)), color="seagreen", label=r'ED')
plt.plot(kh/np.pi, np.real(kmod_E2(kh)), color="blue", linestyle="dotted", label=r'E2')
plt.plot(kh/np.pi, np.real(kmod_E4(kh)), color="black", linestyle="dotted", label=r'E4')
plt.plot(kh/np.pi, np.real(kmod_I4(kh)), color="red", linestyle="dotted", label=r'I4')
plt.xlabel(r'$\frac{kh}{\pi}$') ; plt.ylabel(r'Amplitude error : $\Re\{k^*h\}$') ; plt.legend() ; plt.grid()
plt.savefig(f'../graphs/amp_error.pdf', bbox_inches='tight')
plt.show()


plt.figure(figsize=(6,5))
plt.plot(kh/np.pi, np.zeros(len(kh)),color="dimgray", label=r'$\Im\{\frac{kh}{\pi}\}$')
plt.plot(kh/np.pi, np.imag(kmod_ED(kh)),color="seagreen", label=r'ED')
plt.plot(kh/np.pi, np.imag(kmod_E2(kh)),color="blue", linestyle="dotted", label=r'E2')
plt.plot(kh/np.pi, np.imag(kmod_E4(kh)),color="black", linestyle="dotted", label=r'E4')
plt.plot(kh/np.pi, np.imag(kmod_I4(kh)),color="red", linestyle="dotted", label=r'I4')
plt.xlabel(r'$\frac{kh}{\pi}$') ; plt.ylabel(r'Phase Error : $\Im\{k^*h\}$') ; plt.legend() ; plt.grid()
plt.savefig(f'../graphs/pha_error.pdf', bbox_inches='tight')
plt.show()