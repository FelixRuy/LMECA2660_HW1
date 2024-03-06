import numpy as np
import matplotlib.pyplot as plt 

# PARAMETERS 
U = 1 
L = 1 
sig = 1/4
h = sig/8
N = int(8/sig) 
print(f"U={U}, L={L}, sig={sig}, h={h}, N={N}")

# GRID 
x = np.array([-L/2 + i*h for i in range(int(N))])

# INIT CONDITION 
u_x0 = lambda x : U * np.exp(- (x*x) / (sig*sig))

# EVAL OF N POINTS ON GRID 
y0 = u_x0(x)

# FFT 
yf = np.fft.fft(y0)/N
yf = np.fft.fftshift(yf)
print(yf)

u_analytic = lambda k : U * np.sqrt(np.pi) * sig * np.exp(-k*k*sig*sig/4)



j = np.linspace(-N/2, N/2, 100)
k_analytic = 2*np.pi/L * j
j_disc = np.arange(-N/2, N/2, 1)
k_discrerized = 2*np.pi/L * j_disc

plt.semilogy(j, u_analytic(k_analytic), label='analytic FT')
plt.semilogy(j_disc, np.abs(yf), 'o', color="red", label='discrete FT', markersize=2)
plt.legend(loc='lower center')
plt.grid()
plt.xlabel('j')
plt.ylabel(r'$|\hat{u}(k, 0)|$')
plt.savefig(f'../graphs/gauss_fft_sig{sig}.pdf', bbox_inches='tight')
plt.show()


