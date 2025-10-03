import openmc
import numpy as np
import matplotlib.pyplot as plt

sp = openmc.StatePoint('statepoint.100.h5')
tally = sp.tallies[1]
flux = tally.get_slice(scores=['flux'])
prompt = tally.get_slice(scores=['prompt-nu-fission'])
flux.std_dev.shape = (100, 100)

flux.mean.shape = (100, 100)
prompt.std_dev.shape = (100, 100)
prompt.mean.shape = (100, 100)

fig = plt.subplot(121, title = "Neutron flux distribution")
im1=fig.imshow(flux.mean)
cbar1 = plt.colorbar(im1, ax=fig, fraction=0.046, pad=0.04)
#cbar1.set_label("Flux [n/cm-src]")   # normalizado por partícula fuente

fig2 = plt.subplot(122, title = "Fission site heatmap")
im2=fig2.imshow(prompt.mean)
cbar2 = plt.colorbar(im2, ax=fig2, fraction=0.046, pad=0.04)
#cbar2.set_label("Prompt ν-fission [per src]")

plt.tight_layout()
plt.savefig("flux_heatmap2.png", dpi=300)
plt.show()

