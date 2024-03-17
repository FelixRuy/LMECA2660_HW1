import matplotlib.pyplot as plt
import numpy as np


N = [32, 64, 128, 256, 512, 1024]
E_E2 = [0.508586, 0.096344, 0.007887, 0.000508, 0.0000318499, 0.0000019917578107]
E_E4 = [0.0750693758, 0.0010832376, 0.0000049695, 0.0000000200, 0.0000000000788115, 0.0000000000003084]
E_ED = [0.0842482943279914, 0.0069447411127837, 0.0001943319297496, 0.0000034119433961, 0.0000000542841410, 0.0000000008505934]
E_I4 = [0.0113350003994946, 0.0000488153468465, 0.0000001727436802, 0.0000000006573590, 0.0000000000025510, 0.0000000000000099]



plt.figure(figsize=(7, 5))

plt.loglog(N, E_E2, marker='^', linestyle="dashdot", color="blue", label='E2')
plt.loglog(N, E_E4, marker='^', color="gray", label='E4')
plt.loglog(N, E_ED, marker='^', linestyle="dotted", color="seaGreen", label='ED')
plt.loglog(N, E_I4, marker='^', linestyle="dashed", color="darkorchid", label='I4')


plt.plot(np.linspace(32, 1024, 100), 1e7*np.linspace(32, 1024, 100)**(-4), linestyle=(0, (1, 10)), color="blue", label=r'$\propto N^{-4}$')
plt.plot(np.linspace(32, 1024, 100), 1e8*np.linspace(32, 1024, 100)**(-6), linestyle=(0, (1, 10)), color="green", label=r'$\propto N^{-6}$')
plt.plot(np.linspace(32, 1024, 100), 1e9*np.linspace(32, 1024, 100)**(-8), linestyle=(0, (1, 10)), color="darkorchid", label=r'$\propto N^{-8}$')


plt.xlabel('N')
plt.ylabel(r'$R_h^{0.5}$')
plt.grid()
plt.legend()

plt.savefig("../graphs/error_order.pdf", bbox_inches='tight')
plt.show()