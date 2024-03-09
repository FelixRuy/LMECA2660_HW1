import numpy as np
import matplotlib.pyplot as plt 

# PARAMETERS 
U = 1 
L = 1 
sig = 1/16
h = sig/8
N = int(8/sig) 


############ GAUSSIAN ############
x = np.array([-L/2 + i*h for i in range(int(N))])
u_x0 = lambda x : U * np.exp(- (x*x) / (sig*sig))
y0 = u_x0(x)
yf = np.fft.fft(y0)/N
yf = np.fft.fftshift(yf)
u_analytic = lambda k : U * np.sqrt(np.pi) * sig * np.exp(-k*k*sig*sig/4)
j = np.linspace(-N/2, N/2, 100)
k_analytic = 2*np.pi/L * j
j_disc = np.arange(-N/2, N/2, 1)
k_discrerized = 2*np.pi/L * j_disc
plt.semilogy(j, u_analytic(k_analytic), color="gray", linestyle="dashed", label='analytic FT')
plt.semilogy(j_disc, np.abs(yf), 'o', color="seagreen", label='discrete FT', markersize=2)
plt.legend(loc='lower center')
plt.grid()
plt.xlabel('j')
plt.ylabel(r'$|\hat{u}(k, 0)|$')
#plt.savefig(f'../graphs/gauss_fft_sig{sig}.pdf', bbox_inches='tight')
plt.show() 



############ WAVE PACKET ############
kp = 2*np.pi*12 / L #specific mode
u_analytic = lambda k : 1/2 * U * np.sqrt(np.pi) * sig * (np.exp(-(k - kp)**2 * sig**2 / 4) + np.exp(-(k + kp)**2 * sig**2 / 4))
u_x0 = lambda x : U * np.cos( kp * x)*np.exp(- (x*x) / (sig*sig))
x = np.array([-L/2 + i*h for i in range(int(N))])
y0 = u_x0(x)
yf = np.fft.fft(y0)/N
yf = np.fft.fftshift(yf)
j = np.linspace(-N/2, N/2, 100)
k_analytic = 2*np.pi/L * j
j_disc = np.arange(-N/2, N/2, 1)
k_discrerized = 2*np.pi/L * j_disc
plt.semilogy(j_disc, np.abs(yf), 'o', color="seagreen", label='discrete FT', markersize=2)
plt.semilogy(j, u_analytic(k_analytic), color="gray", linestyle="dashed", label='analytic FT')
#plt.vlines(j_disc, 0, np.abs(yf), colors='seagreen', lw=0.7, alpha=0.5)
plt.grid()
plt.xlabel('j')
plt.ylabel(r'$|\hat{u}(k, 0)|$')
plt.legend()
#plt.savefig(f'../graphs/gauss_fft_sig{sig}_wave_log.pdf', bbox_inches='tight')
plt.show()



