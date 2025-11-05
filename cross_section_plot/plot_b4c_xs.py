# plot_b4c_mt102_107.py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_column(fname):
    with open(fname, "r", encoding="utf-8") as f:
        clean_lines = (line.replace(",", "").strip() for line in f if line.strip())
        return np.loadtxt(clean_lines)

# Archivos (mismo directorio que el script)
e107 = load_column("107_energy.txt")
xs107 = load_column("107_cross.txt")
e102 = load_column("102_energy.txt")
xs102 = load_column("102_cross.txt")

# Comprobaciones rápidas
assert e107.size == xs107.size, "107_energy y 107_cross deben tener misma longitud"
assert e102.size == xs102.size, "102_energy y 102_cross deben tener misma longitud"

plt.figure(figsize=(7, 4.5))
plt.loglog(e107, xs107, linewidth=1.8, label="Boron Carbide MT-107")
plt.loglog(e102, xs102, linewidth=1.8, label="Boron Carbide MT-102")

plt.xlabel("Energy (eV)")
plt.ylabel("Macroscopic cross section (cm$^{-1}$)")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("b4c_mt102_mt107.png", dpi=200)
plt.show()


