import openmc
import numpy as np
import matplotlib.pyplot as plt

# EXERCISE A---------------------------------------------------------------------------

#  UO2
urox = openmc.Material(name = "Uranium dioxide")
# This makes an instance of material that we may modify at will
urox.add_nuclide("U235", 0.02)
# This above adds uranium 235 and sets the relative amount of atoms with the float.
urox.add_nuclide("U238", 0.98)
urox.add_element("O", 2.0)
urox.set_density("g/cm3",10.97)
# Theoretical density of UO2 at room temperature, change when needed
urox.add_s_alpha_beta("c_O_in_UO2")
# Enhanced scattering cross sections thermal regions

# Light water (H2O)
light_water = openmc.Material(name="Light water")
light_water.add_nuclide("H1", 2.0)
light_water.add_element("O", 1.0)
light_water.set_density("g/cm3", 1.0)
light_water.add_s_alpha_beta("c_H_in_H2O")
# Enhanced scattering cross sections thermal regions

# Helium 
helium = openmc.Material(name="Helium")
helium.add_nuclide("He4", 1.0)
helium.set_density("g/cm3", 0.0001786)  

# Zircaloy (cladding)
zircaloy = openmc.Material(name="Zircaloy")
zircaloy.add_element("Zr", 1.0)
zircaloy.set_density("g/cm3", 6.55)  

materials = openmc.Materials([urox, helium, zircaloy, light_water])
materials.export_to_xml()

# EXERCISE B-----------------------------------------------------------------------------------------------

# Box of size 1.4 cm x 1.4 cm x 2 cm with “vacuum” boundary condition
box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type='vacuum')
zmin = openmc.ZPlane(-1.0, boundary_type='vacuum')
zmax = openmc.ZPlane(+1.0, boundary_type='vacuum')

# Concentric cylinders: fuel, cavity, coating
c_fuel = openmc.ZCylinder(r=0.40) # 4.0 mm
c_gap = openmc.ZCylinder(r=0.42) # 4.2 mm
c_clad = openmc.ZCylinder(r=0.45) # 4.5 mm

common_axial = (-box)&(+zmin)&(-zmax)

fuel_region = (-c_fuel)&common_axial # inside the fuel
gap_region = (+c_fuel)&(-c_gap)&common_axial # between fuel and gap
clad_region = (+c_gap)&(-c_clad)&common_axial # between gap and clad
water_region = (+c_clad)&common_axial # outside clad, inside the box

fuel_cell = openmc.Cell(name='fuel',fill=urox,region=fuel_region)
gap_cell = openmc.Cell(name='gap',fill=helium,region=gap_region)
clad_cell = openmc.Cell(name='clad',fill=zircaloy,region=clad_region)
mod_cell = openmc.Cell(name='water',fill=light_water,region=water_region)

# Universe and Geometry
root = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, mod_cell])
geom = openmc.Geometry(root)
geom.export_to_xml()

# Plots
colors = {urox: 'green', helium: 'brown', zircaloy: 'gray', light_water: 'blue'}

# Plot XY
plot_xy           = openmc.Plot()
plot_xy.filename  = 'pin_xy'
plot_xy.basis     = 'xy'
plot_xy.width     = (1.4, 1.4)
plot_xy.pixels    = (800, 800)
plot_xy.color_by  = 'material'
plot_xy.colors    = colors

# Plot XZ
plot_xz           = openmc.Plot()
plot_xz.filename  = 'pin_xz'
plot_xz.basis     = 'xz'
plot_xz.width     = (1.4, 2.0)    # ancho x, alto z
plot_xz.pixels    = (800, 800)
plot_xz.color_by  = 'material'
plot_xz.colors    = colors

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()

# Generar las imágenes
openmc.plot_geometry()

# EXERCISE C-----------------------------------------------------------------------------------------------

settings = openmc.Settings()
settings.batches = 100 # How many iterations of neutron generations
settings.inactive = 10 # How many generations used to determine sample spots
settings.particles = 10000 # Number of neutrons in each generation

materials = openmc.Materials([light_water, urox, helium, zircaloy])
materials.export_to_xml()

# The following section tells OpenMC where to sample particle sites.
bounds = [-0.7,-0.7,-1,0.7,0.7,1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()

if __name__ == "__main__":
 openmc.run(threads = 4) # Threads set for heplab primarily

# k value
sp = openmc.StatePoint('statepoint.100.h5')
print(sp.k_combined) # has returned 0.009488+/-0.000015