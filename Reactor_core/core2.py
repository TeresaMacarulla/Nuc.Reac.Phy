# conda install -c conda-forge openmc-plotter
# python -m pip install --upgrade "git+https://github.com/openmc-dev/plotter.git@develop"
# openmc-plotter
import openmc
import numpy as np
import matplotlib.pyplot as plt
from utilities import build_assembly_universe
import openmc.deplete as deplete
import math

# ============================================================
# Existing objects
# ============================================================

# assembly_universe
# water
# stainless_steel
# vessel_steel
assembly_uni, materials, light_water = build_assembly_universe(with_ctrl_3x3=True)
steel = next(mat for mat in materials if mat.name == "Steel")
steel2 = next(mat for mat in materials if mat.name == "Steel2")

# ============================================================
# Water universe
# ============================================================

water_cell = openmc.Cell(fill=light_water)
water_universe = openmc.Universe(name="Water universe", cells=[water_cell])

# ============================================================
# Core lattice
# ============================================================

pin_pitch = 1.4
pins_per_assembly = 17
assembly_pitch = pin_pitch * pins_per_assembly   # 23.8 cm
N = 15
core_array = np.full((N, N), water_universe, dtype=openmc.Universe)
center = (N - 1) / 2

core_array[5:10, 1] = assembly_uni
core_array[4:11, 2] = assembly_uni
core_array[3:12, 3] = assembly_uni
core_array[2:13, 4] = assembly_uni
core_array[1:14, 5:10] = assembly_uni
core_array[2:13, 10] = assembly_uni
core_array[3:12, 11] = assembly_uni
core_array[4:11, 12] = assembly_uni
core_array[5:10, 13] = assembly_uni

core_lattice = openmc.RectLattice(name="Core lattice")
core_lattice.pitch = (assembly_pitch, assembly_pitch)
core_lattice.lower_left = (-N * assembly_pitch / 2, -N * assembly_pitch / 2)
core_lattice.universes = core_array
core_lattice.outer = water_universe

# ============================================================
# Radial reactor structures
# ============================================================

barrel_inner = openmc.ZCylinder(r=180.0)
barrel_outer = openmc.ZCylinder(r=190.0)
vessel_inner = openmc.ZCylinder(r=190.0)
vessel_outer = openmc.ZCylinder(r=220.0, boundary_type="vacuum")

# ============================================================
# Axial boundaries
# ============================================================

core_height = 300.0
z_bottom = openmc.ZPlane(z0=-core_height / 2, boundary_type="vacuum")
z_top = openmc.ZPlane(z0=core_height / 2, boundary_type="vacuum")
axial = +z_bottom & -z_top

# ============================================================
# Cells
# ============================================================

core_cell = openmc.Cell(name="Core and moderator", fill=core_lattice, region=-barrel_inner & axial)
barrel_cell = openmc.Cell(name="Core barrel", fill=steel, region=+barrel_inner & -barrel_outer & axial)
vessel_cell = openmc.Cell(name="Pressure vessel", fill=steel2, region=+vessel_inner & -vessel_outer & axial)

# ============================================================
# Root universe / geometry / settings
# ============================================================

root = openmc.Universe(name="Reactor", cells=[core_cell, barrel_cell, vessel_cell])
geometry = openmc.Geometry(root)
geometry.export_to_xml()

core_radius= 220
hz = 150

settings = openmc.Settings()
settings.batches = 150
settings.inactive = 30
settings.particles = 10000 #30000
bounds = [-core_radius, -core_radius, -hz, core_radius, core_radius, hz]
src = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.IndependentSource(space=src)
settings.export_to_xml()

# ============================================================
# Neutron flux core
# ============================================================

mesh = openmc.RegularMesh()
mesh.dimension   = [200, 200, 1]
mesh.lower_left  = [-core_radius, -core_radius, -hz]
mesh.upper_right = [ core_radius,  core_radius,  hz]
mf = openmc.MeshFilter(mesh)

tallies = openmc.Tallies()

# Total flux
t_flux = openmc.Tally(name="core_flux_mesh")
t_flux.filters = [mf]
t_flux.scores  = ["flux"]
tallies.append(t_flux)
tallies.export_to_xml()

openmc.run()

# Open statepoint
sp = openmc.StatePoint("statepoint.150.h5")
# Get flux tally
tally = sp.get_tally(name="core_flux_mesh")
# Reshape according to mesh dimensions
flux = tally.get_reshaped_data(value="mean", expand_dims=True)

# Shape should essentially be:
# (200, 200, 1, 1, 1)
# x, y, z, nuclide, score
flux = flux.squeeze()      # -> (200, 200)

# The raw OpenMC flux tally is accumulated neutron track length per source particle. 
# We want volume-averaged flux with units
dx = core_radius*2 / 200
dy = core_radius*2 / 200
dz = hz / 1
voxel_volume = dx * dy * dz
flux = flux / voxel_volume

plt.figure(figsize=(8, 7))
plt.imshow(flux.T,  origin="lower", extent=[-core_radius, core_radius, -core_radius, core_radius])
plt.xlabel("x [cm]")
plt.ylabel("y [cm]")
plt.title("Neutron flux distribution")
plt.colorbar(label="Flux [1/cm²/source]")
plt.tight_layout()
plt.show()

# ============================================================
# Depletion
# ============================================================

"""# Designate the volume of the fuel 
for m in materials:
    if (m.name or "") == "Uranium dioxide":
        target = m

n_assemblies = 13*5 + 32*2
pins_per_assembly = 17 * 17
total_pins_all = n_assemblies * pins_per_assembly
total_ctrl_pins = n_assemblies * (3 * 3)
n_fuel_pins = total_pins_all - total_ctrl_pins
H_core_cm = 300
pin_area = math.pi * (0.4 ** 2)  # cm^2
V_fuel_cm3 = n_fuel_pins * pin_area * H_core_cm

target.volume = float(V_fuel_cm3)

# Depletion
reactor_model = openmc.Model()
reactor_model.materials = materials
reactor_model.geometry  = geometry
reactor_model.settings  = settings
reactor_model.tallies = tallies
chain = openmc.deplete.Chain.from_xml("./chain_endfb80_pwr.xml")
#Designating the depletion chain being used for our depletion calculation
operator = openmc.deplete.CoupledOperator(reactor_model, "./chain_endfb80_pwr.xml")
power = 1200e6 #Designating power settings in watts
time_steps = [30 * 24 * 60 * 60] * 6
#Number of timesteps and their length in this cae 6 steps of one month length in seconds.
integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power)
#Running the depletion calculation
integrator.integrate()

# Read the results
results = openmc.deplete.ResultsList.from_hdf5("./depletion_results.h5")
time, k = results.get_eigenvalue()
time /= (24 * 60 * 60)  # convert back to days from seconds

plt.errorbar(time, k[:, 0], yerr=k[:, 1])
plt.xlabel("Time [days]")
plt.ylabel("$k_{eff}\pm \sigma$");
plt.grid(True)
plt.show()

t, u5 = results.get_atoms("1", "U235")
t, pu239 = results.get_atoms("1", "Pu239")
t /= (24 * 60 * 60)  # convert back to days from seconds
plt.plot(t, u5, label="U235")
plt.plot(t, pu239, label="Pu239")"""