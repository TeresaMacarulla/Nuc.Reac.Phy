import openmc
import numpy as np
import matplotlib.pyplot as plt

# EXERCISE A---------------------------------------------------------------------------

#  UO2
urox = openmc.Material(name = "Uranium dioxide")
urox.add_nuclide("U235", 0.02)
urox.add_nuclide("U238", 0.98)
urox.add_element("O", 2.0)
urox.set_density("g/cm3",10.97)
# Enhanced scattering cross sections thermal regions
urox.add_s_alpha_beta("c_O_in_UO2")

# Light water (H2O)
light_water = openmc.Material(name="Light water")
light_water.add_nuclide("H1", 2.0)
light_water.add_element("O", 1.0)
light_water.set_density("g/cm3", 1.0)
# Enhanced scattering cross sections thermal regions
light_water.add_s_alpha_beta("c_H_in_H2O")

# Helium 
helium = openmc.Material(name="Helium")
helium.add_nuclide("He4", 1.0)
helium.set_density("g/cm3", 0.0001786)  

# Zircaloy (cladding)
zircaloy = openmc.Material(name="Zircaloy")
zircaloy.add_element("Zr", 1.0)
zircaloy.set_density("g/cm3", 6.52)  

materials = openmc.Materials([urox, helium, zircaloy, light_water])
materials.export_to_xml()

# EXERCISE B-----------------------------------------------------------------------------------------------

# Box of size 1.4 cm x 1.4 cm x 2 cm with “vacuum” or "reflective" boundary condition
box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type='reflective')
zmin = openmc.ZPlane(-1.0, boundary_type='reflective')
zmax = openmc.ZPlane(+1.0, boundary_type='reflective')

# Concentric cylinders: fuel, cavity, coating
c_fuel = openmc.ZCylinder(r=0.40, boundary_type='transmission') # 4.0 mm
c_gap = openmc.ZCylinder(r=0.42, boundary_type='transmission') # 4.2 mm
c_clad = openmc.ZCylinder(r=0.45, boundary_type='transmission') # 4.5 mm

box_outside_cylinder = (-box)&(+zmin)&(-zmax)#&(~c_fuel)
volume_c_fuel=(-c_fuel)&(+zmin)&(-zmax)
volume_c_gap=(-c_gap)&(+zmin)&(-zmax)
volume_c_clad=(-c_clad)&(+zmin)&(-zmax)

fuel_region = (-c_fuel)&box_outside_cylinder # inside the fuel
gap_region = (+c_fuel)&(-c_gap)&box_outside_cylinder # between fuel and gap
clad_region = (+c_gap)&(-c_clad)&box_outside_cylinder # between gap and clad
water_region = (+c_clad)&box_outside_cylinder # outside clad, inside the box

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
plot_xz.width     = (1.4, 2.0)    # width x, height z
plot_xz.pixels    = (800, 800)
plot_xz.color_by  = 'material'
plot_xz.colors    = colors

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()

# Plot generation
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

# k value obtained so far
# k-effective (Collision)     = 0.00954 +/- 0.00004
# k-effective (Track-length)  = 0.00950 +/- 0.00001
# k-effective (Absorption)    = 0.00964 +/- 0.00015
# Combined k-effective        = 0.00950 +/- 0.00001
# Leakage Fraction            = 0.99559 +/- 0.00006
# This happens because my conditions are 'vacuum' and my neutrons are scaping from the reactor

# EXERCISE D-----------------------------------------------------------------------------------------------

# Chage "vacuum" conditions to "reflective" conditions in lines 40, 41, 42.

# k-effective (Collision)     = 1.29182 +/- 0.00128
# k-effective (Track-length)  = 1.29260 +/- 0.00151
# k-effective (Absorption)    = 1.29092 +/- 0.00107
# Combined k-effective        = 1.29147 +/- 0.00086
# Leakage Fraction            = 0.00000 +/- 0.00000

# EXERCISE E-----------------------------------------------------------------------------------------------

mesh = openmc.RegularMesh() #Making a mesh instance
mesh.dimension = [100,100] #Setting the bins of the mesh to 100 x 100 bins
mesh.lower_left= [-0.7,-0.7] #Setting lower left coordinates
mesh.upper_right= [0.7,0.7] #Setting upper right coordinates

mesh_filter = openmc.MeshFilter(mesh) #Making a meshfilter
tallies_file = openmc.Tallies() #For our tallies to be included we have to put them into an input file
 
tally = openmc.Tally(name = "Tally") #Making a tally instance
tally.filters = [mesh_filter] #Designatingen the mesh for the tally
tally.scores = ["flux", 'prompt-nu-fission'] #Designating the objects to be tallied
tallies_file.append(tally) #Appending the tally to our tallies file

tallies_file.export_to_xml()

# EXERCISE F-----------------------------------------------------------------------------------------------

# Define an energy filter (example: from 1e-5 eV to 20 MeV, 100 bins log-spaced)
energy_bins = np.logspace(-5, 7, 200)   # Adjust the number of bins 
energy_filter = openmc.EnergyFilter(energy_bins)

tally2 = openmc.Tally(name = "Energy") #Making a new tally only looking at flux
tally2.filters = [energy_filter] #Designating the energy filter
tally2.scores = ["flux"]
tallies_file.append(tally2)

tallies_file.export_to_xml()
openmc.run()

# EXERCISE G and H-----------------------------------------------------------------------------------------------

# Change "Light water" for "Heavy water" or "Graphite". This exercise is in "heavy_water/Project2.py"

# EXERCISE I----------------------------------------------------------------------------------------------------

# Obtain the thermal utilization factor f vs. V from "Volume_variation.py"