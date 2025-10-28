import numpy as np
import matplotlib.pyplot as plt

# Datos
T = np.array([400, 600, 800, 900, 1000, 1200])  # Temperatura [K]
k = np.array([1.16313, 1.15820, 1.15361, 1.15167, 1.15005, 1.14690])
dk = np.array([0.00052, 0.00051, 0.00047, 0.00049, 0.00046, 0.00050])

# Gráfica
plt.figure(figsize=(6,4))
plt.errorbar(T, k, yerr=dk, fmt='o-', capsize=4, lw=1.5)
plt.xlabel('Fuel temperature [K]', fontsize=11)
plt.ylabel(r'$k_{\mathrm{eff}}$', fontsize=11)
plt.grid(True, ls='--', alpha=0.6)
plt.tight_layout()
plt.show()
