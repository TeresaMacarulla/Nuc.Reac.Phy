import openmc
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

sp = openmc.StatePoint("statepoint.100.h5")

# Recupera el tally por nombre (evita IDs frágiles)
tally = sp.get_tally(name="Energy")

# Localiza el EnergyFilter (puede no ser el [0] si añadiste otros filtros)
ef = next(f for f in tally.filters if isinstance(f, openmc.EnergyFilter))
edges = ef.bins                 # edges de energía (eV), longitud = N+1

y = tally.mean.ravel()          # valores del espectro, longitud = N
# Normaliza longitudes
if edges.ndim == 1:
    N = min(len(edges) - 1, len(y))
    edges = edges[:N+1]
    y     = y[:N]
    E_mid = np.sqrt(edges[:-1] * edges[1:])        # punto medio geométrico
elif edges.ndim == 2 and edges.shape[1] == 2:
    N = min(edges.shape[0], len(y))
    edges = edges[:N, :]
    y     = y[:N]
    E_mid = np.sqrt(edges[:, 0] * edges[:, 1])     # punto medio de cada par
else:
    raise RuntimeError(f"Formato inesperado de edges: shape={edges.shape}")

# Para log-log, descarta valores no positivos
mask = y > 0.0
E_plot = E_mid[mask]
y_plot = y[mask]

plt.figure()
print("E_plot:", E_plot.shape, " y_plot:", y_plot.shape)

plt.loglog(E_plot, y_plot)
ax = plt.gca()
print("series en el eje:", len(ax.lines))

plt.grid(True, which="both")
plt.xlabel("Energy [eV]")
plt.ylabel("Flux [n/cm$^2$/src]")   # unidades correctas
plt.title("Neutron energy spectrum in the fuel")
plt.tight_layout()
plt.savefig("flux_energy_PuO2_50%.pdf", dpi=300)
plt.show()
