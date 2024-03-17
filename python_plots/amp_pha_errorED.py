import numpy as np
import matplotlib.pyplot as plt



# Phase and amplification factor for partially decentered scheme (ED)

kh = np.linspace(0, np.pi, 100)

kmod_ED = lambda kh : 1/(6*1j) * (np.exp(-2*kh*1j) - 6*np.exp(-kh*1j) + 3 + 2*np.exp(kh*1j)) / np.pi
kmod_ED1 = lambda kh : (4/3*np.sin(kh)-1/3*np.cos(kh)*np.sin(kh))/np.pi + 1j*(-1/3 -1/3 * np.cos(kh)**2 + 2/3*np.cos(kh)) / np.pi
kmod_E2 = lambda kh : np.sin(kh) / np.pi
kmod_E4 = lambda kh : np.sin(kh)*(1+1/3*(1-np.cos(kh))) / np.pi
kmod_I4 = lambda kh : ( 3/2 * np.sin(kh) / (1+1/2*np.cos(kh)) ) / np.pi

plt.figure(figsize=(6,5))
plt.plot(kh/np.pi, kh/np.pi,color="black", label=r'$\Re\{\frac{kh}{\pi}\}$')
plt.plot(kh/np.pi, np.real(kmod_E2(kh)), color="blue", linestyle="dashdot", label=r'E2')
plt.plot(kh/np.pi, np.real(kmod_ED(kh)), color="seaGreen", linestyle="dotted", label=r'ED', marker="o", markersize=3, alpha=0.5)
plt.plot(kh/np.pi, np.real(kmod_ED1(kh)), color="red")
plt.plot(kh/np.pi, np.real(kmod_E4(kh)), color="gray", label=r'E4')
plt.plot(kh/np.pi, np.real(kmod_I4(kh)), color="darkorchid", linestyle="dashed", label=r'I4')
plt.xlabel(r'$\frac{kh}{\pi}$') ; plt.ylabel(r'Amplitude error : $\Re\{k^*h\}$') ; plt.legend() ; plt.grid()
#plt.savefig(f'../graphs/amp_error.pdf', bbox_inches='tight')
plt.show()


plt.figure(figsize=(6,5))
plt.plot(kh/np.pi, np.zeros(len(kh)),color="black", label=r'$\Im\{\frac{kh}{\pi}\}$')
plt.plot(kh/np.pi, np.imag(kmod_ED(kh)),color="seaGreen", linestyle="dotted", label=r'ED')
plt.plot(kh/np.pi, np.imag(kmod_E2(kh)),color="blue", linestyle="dashdot", label=r'E2')
plt.plot(kh/np.pi, np.imag(kmod_E4(kh)),color="gray", label=r'E4')
plt.plot(kh/np.pi, np.imag(kmod_ED1(kh)), color="red")
plt.plot(kh/np.pi, np.imag(kmod_I4(kh)),color="darkorchid", linestyle="dashed", label=r'I4')
plt.xlabel(r'$\frac{kh}{\pi}$') ; plt.ylabel(r'Phase Error : $\Im\{k^*h\}$') ; plt.legend() ; plt.grid()
#plt.savefig(f'../graphs/pha_error.pdf', bbox_inches='tight')
plt.show()